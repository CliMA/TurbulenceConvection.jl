"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""

function handle_expr(expr::String; kwargs...)
    # see https://stackoverflow.com/a/57749395, takes a function that accepts kwargs
    # create tuple from args and create func
    
    if ~contains(expr, "->")
        expr  = "(" * join([string(x) for x in keys(kwargs)],',') *  ")" * " -> " * expr
    end

    expr  = eval(Meta.parse(expr))
    return Base.invokelatest(expr,values(NamedTuple(kwargs))...) # works but very slow to use invokelatest, maybe try https://github.com/SciML/RuntimeGeneratedFunctions.jl
end


function noneq_moisture_sources(param_set::APS, area::FT, ρ::FT, Δt::Real, ts, w, z; ts_LCL=nothing) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    if area > 0
        use_supersat = get(param_set.user_args, :use_supersat, false) # (:use_supersat in keys(param_set.user_args)) ? param_set.user_args.use_supersat : false # so we dont have to set everything we dont know is in user_args in the defaults...
        use_supersat = isa(use_supersat, Val) ? typeof(use_supersat).parameters[1] : use_supersat # extract the value from the Val type (we put it in there in parameter_set.jl to make it isbits)
        use_korolev_mazin = get(param_set.user_args, :use_korolev_mazin, false) # (:use_supersat in keys(param_set.user_args)) ? param_set.user_args.use_supersat : false # so we dont have to set everything we dont know is in user_args in the defaults...
        use_korolev_mazin = isa(use_korolev_mazin, Val) ? typeof(use_korolev_mazin).parameters[1] : use_korolev_mazin # extract the value from the Val type (we put it in there in parameter_set.jl to make it isbits)

        q = TD.PhasePartition(thermo_params, ts)
        T = TD.air_temperature(thermo_params, ts)
        p = TD.air_pressure(thermo_params, ts)
        q_vap = TD.vapor_specific_humidity(thermo_params, ts)

        if !(use_supersat == false) # it's either true or a specified value
            # use phase partition in case we wanna use the conv_q_vap fcn but maybe not best for supersat since it's not really a phase partition (all 3 are vapor amounts)
            D = FT((2.11 * 1e-5) * (T/273.15)^1.94 *  (p/101325))  # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
            
            if use_supersat == true || use_supersat == :Base
                # @info("using default supersaturation formulation")
                # nonmutable param_set so would have to edit tau usage everywhere -- let's just do in supersat
                tau_weights = get(param_set.user_aux,  :tau_weights, nothing) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
                if isnothing(tau_weights)
                    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
                    τ_ice = CMNe.τ_relax(microphys_params, ice_type)
                else
                    #use the weights to calculate tau
                    # @info("using τ from auxiliary weight function")
                    # τ_liq, τ_ice = 10 .^ tau_weights.liq, 10 .^ tau_weights.ice # log fcn
                    # tau_weights = eval(Meta.parse(tau_weights)) # test string version cause was crashing....
                    # τ_liq = handle_expr(string(tau_weights.liq.func_expr); liq_params) # works! (turn to string, parse to func and eval with the argument) (or use handle_expr above... :), # works but very slow to use invokelatest, maybe try https://github.com/SciML/RuntimeGeneratedFunctions.jl
                    # τ_ice = handle_expr(string(tau_weights.ice.func_expr); ice_params) # works! (turn to string, parse to func and eval with the argument) (or use handle_expr above... :)
                    liq_params, ice_params = tau_weights.liq.liq_params, tau_weights.ice.ice_params # gets used in eval below
                    τ_liq, τ_ice = 10 .^ liq_params.log10_tau_liq, 10 .^ ice_params.log10_tau_ice # log fcn hand implementation...
                end

                min_τ_liq = get(param_set.user_aux, :min_τ_liq, FT(1)) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
                min_τ_ice = get(param_set.user_aux, :min_τ_ice, FT(1)) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
                τ_liq = max(min_τ_liq, τ_liq) # make sure it's at least some value
                τ_ice = max(min_τ_ice, τ_ice) # make sure it's at least some value

            elseif use_supersat == :exponential_T_scaling_ice
                # TODO: Drop the constants from this and just subsume them into N (makes picking an initial c_1, c_2 harder though so maybe not? Also is much harder for more complicated N_INP expressions)
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice 
                # c_1 = get(param_set.user_aux, :T_scaling_c_1, 0.02 * exp(-0.6 * -273.15)) # Fletcher 1962 (values taken from Frostenberg 2022) (high power makes calibrating hard)
                c_1 = get(param_set.user_aux, :T_scaling_c_1, 0.02) # Fletcher 1962 (values taken from Frostenberg 2022)
                c_2 = get(param_set.user_aux, :T_scaling_c_2, -0.6) # Fletcher 1962 (values taken from Frostenberg 2022)
                N = c_1 * exp(c_2*(T-273.15))
                r_r = FT(20 * 1e-6) # 20 microns
                r_0 = FT(.2 * 1e-6) # .2 micron base aerosol
                # r_0 = FT(20. * 10^-6) # 20 micron
                r_r = (q.ice / (4/3 * π * N * ρ_i)) ^ (1/3) # (the mass-diameter relationship is poorly defined anyway for ice crystals) -- if q.ice is 0 this goes to 0... making it hard to generate ice
                r_r = max(r_r, r_0) # bound to be at least ~micron size...something like kohler crit radius
                τ_ice = 1 / (4 * π * D * N * r_r)
                τ_liq = CMNe.τ_relax(microphys_params, liq_type)

            elseif use_supersat == :powerlaw_T_scaling_ice
                # @info("using $use_supersat supersaturation formulation")
                error("NotImplmentedError: This supersat_type functionality has not been implemented yet")
            elseif use_supersat == :exponential_times_powerlaw_scaling_ice # Demott 2010
                # @info("using $use_supersat supersaturation formulation")
                error("NotImplmentedError: This supersat_type functionality has not been implemented yet")

            elseif use_supersat == :geometric_liq # scaling on q that impacts liquid
                r_0 = FT(20. * 10^-6) # 20 micron
                ρ_l = FT(1000.) # kg m^-3, CLIMAParameters default for cloud_liquid
                c_1 = get(param_set.user_aux, :geometric_liq_c_1, 1/(4/3 * π * ρ_l * r_0^2)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                c_2 = get(param_set.user_aux, :geometric_liq_c_2, 2/3.) # Halfway between 1/3 and 1
                τ_liq = 1 / (4 * π * D * c_1 * q.liq^(c_2)) # let  be Nr = c_1 * q^(c_2) 
                τ_ice = CMNe.τ_relax(microphys_params, liq_type)

            elseif use_supersat == :geometric_liq__geometric_ice # scaling on q that impacts liquid
                ρ_l = FT(1000.) # kg m^-3, CLIMAParameters default for cloud_liquid
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice 
                r_r = FT(20 * 1e-6) # 20 microns
                r_0 = FT(.2 * 1e-6) # .2 micron base aerosol
                N_l = FT(1e-5 / (4/3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q_liq.. (N = N_r in homogenous)
                N_i = FT(1e-7 / (4/3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q_ice... (N = N_r + N_0)

                # need a threshold bc otherwise, if there no condensate, τ will be infinity and we'll never get any since we don't have separate activation scheme...
                # q_liq = max(q.liq, FT(1e-12)) # replaced with a tunable parameter
                # q_ice = max(q.ice, FT(1e-12))

                q_liq = q.liq
                q_ice = q.ice

                c_1l = get(param_set.user_aux, :geometric_liq_c_1, FT(1/(4/3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                c_2l = get(param_set.user_aux, :geometric_liq_c_2, FT(2/3.)) # Halfway between 1/3 and 1
                c_3l = get(param_set.user_aux, :geometric_liq_c_3, FT(N_l * r_0) ) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                τ_liq = 1 / (4 * π * D * (c_1l * q_liq^(c_2l) + c_3l)) # let  be Nr = c_1 * q^(c_2) + c_3

                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice
                c_1i = get(param_set.user_aux, :geometric_ice_c_1, 1/(4/3 * π * ρ_i * r_r^2)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                c_2i = get(param_set.user_aux, :geometric_ice_c_2, 2/3.) # Halfway between 1/3 and 1
                c_3i = get(param_set.user_aux, :geometric_ice_c_3, FT(N_i * r_0) ) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

                τ_ice = 1 / (4 * π * D * (c_1i * q_ice^(c_2i) + c_3i)) # let  be Nr = c_1 * q^(c_2) + c_3

            elseif use_supersat == :geometric_liq__exponential_T_scaling_ice #
                r_0 = FT(20. * 10^-6) # 20 micron
                ρ_l = FT(1000.) # kg m^-3, CLIMAParameters default for cloud_liquid
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice
                c_1g = get(param_set.user_aux, :geometric_c_1, 1/(4/3 * π * ρ_l * r_0^2)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                c_2g = get(param_set.user_aux, :geometric_c_2, 2/3.) # Halfway between 1/3 and 1
                τ_liq = 1 / (4 * π * D * (c_1g * q.liq^(c_2g) + c_3g)) # let  be Nr = c_1 * q^(c_2)

                c_1 = get(param_set.user_aux, :T_scaling_c_1, 0.02) # Fletcher 1962 (values taken from Frostenberg 2022)
                c_2 = get(param_set.user_aux, :T_scaling_c_2, -0.6) # Fletcher 1962 (values taken from Frostenberg 2022)
                N = c_1 * exp(c_2*(T-273.15))
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice 
                r_0 = (q.ice / (4/3 * π * N * ρ_i)) ^ (1/3) # (the mass-diameter relationship is poorly defined anyway for ice crystals)
                τ_ice = 1 / (4 * π * D * N * r_0)

            elseif use_supersat == :geometric_liq__exponential_T_scaling_and_geometric_ice # scaling on q that impacts liquid and ice

                q_liq = q.liq
                q_ice = q.ice

                ρ_l = FT(1000.) # kg m^-3, CLIMAParameters default for cloud_liquid
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice

                r_r = 20 * 1e-6 # 20 microns
                r_0 = .2 * 1e-6 # .2 micron base aerosol
                N_l = FT(1e-5 / (4/3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q_liq.. (N = N_r in homogenous)
                N_i = FT(1e-7 / (4/3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q_ice... (N = N_r + N_0)

                c_1l = get(param_set.user_aux, :geometric_liq_c_1, FT(1/(4/3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                c_2l = get(param_set.user_aux, :geometric_liq_c_2, 2/3.) # Halfway between 1/3 and 1
                c_3l = get(param_set.user_aux, :geometric_liq_c_3, FT(N_l * r_0) ) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
                τ_liq = 1 / (4 * π * D * (c_1l * q_liq^(c_2l) + c_3l)) # let  be Nr = c_1 * q^(c_2) 

                c_1i = get(param_set.user_aux, :geometric_and_T_scaling_ice_c_1, FT((4*π*D) * ((4/3 * π * ρ_i)^(-1/3)  * (N_i)^(2/3) * (0.02)^(2/3) + (N_i * r_0))) ) # Yeahhh.... idk for this one lol... just combined them serially from the homogenous case where c_3 is -1/3, and used .02 as the prefactor
                c_2i = get(param_set.user_aux, :geometric_and_T_scaling_ice_c_2, FT(2/3.)) # Halfway between 1/3 and 1 -- should this be the same as c_2g? It's the same mixing... 
                c_3i = get(param_set.user_aux, :geometric_and_T_scaling_ice_c_3, FT((4*π*D) * r_0 * .02 ) ) # Fletcher 1962 (values taken from Frostenberg 2022) and used .02 as the prefactor
                c_4i = get(param_set.user_aux, :geometric_and_T_scaling_ice_c_4, FT(-0.6) ) # Fletcher 1962 (values taken from Frostenberg 2022)

                τ_ice = 1/ ((c_1i * q_ice^(c_2i) * exp(c_2i*c_4i*(T-273.15))  + c_3i) * exp(c_4i*(T-273.15)) )

                τ_liq = min(τ_liq, FT(1e10)) # ensure we don't just get a huge number all the time
                τ_ice = min(τ_ice, FT(1e10)) # ensure we don't just get a huge number all the time


            elseif use_supersat == :raymond_ice_test # should be closer to DeMott 2015 but...
                N_r_closure = get(param_set.user_args, :N_r_closure, :monodisperse)
                N_0 = get(param_set.user_aux,  :N_0, nothing) |> x ->  ( isnothing(x) ? get(param_set.user_args, :N_0, FT(100*10e6)) : pyinterp(z, x.z, x.values) ) # if defined in aux, interp at current z, otherwise default to user_args value, otherwise default to a value
                # @info("using $use_supersat supersaturation formulation")
                N_m = get(param_set.user_args, :N_m, FT(-0.2)) # log slope https://doi.org/10.1073/pnas.1514034112
                N_b = get(param_set.user_args, :N_b, FT(-5 - 273.15 * N_m)) # -5 - 273.15 * N_m
                N_INP = 10^(N_m*T + N_b) * 10^3 # per liter to per m^3
                # -- hopefully this gets around the problem of draining water vapor at initialization of clouds but also allows speedup as droplets grow (assuming fixed drop concenctration)
                # R    =  max(((q.liq + q.ice)/(4/3*π*ρ_l*N_0))^(1/3), FT(0.2*10^-6)) # bound to be at least ~micron size...something like kohler crit radius

                if N_r_closure == :inhomogeneous
                    N_0,R_liq = NR_inhomogeneous_mixing_liquid(thermo_params, N_0, TD.air_pressure(thermo_params,ts), q.liq, ts_LCL) # testing inhomogeneous mixing would have r fixed and then let N vary... set r based on adiabatic extrapolation from cloud base 
                    _,R_ice = NR_monodisperse(N_INP,q.ice)
                    # maybe look into using a mixed model where N/R are partly towards the inhomogenous value depending on the true entrainment/mixing params... see literature on this
                    # NOTE, ON DYCOMS ADIABATIC R W/ ORIGINAL N_0 WORKED BETTER... HMMM (though that's not using DYCOMS N)

                elseif N_r_closure == :monodisperse # uniform size for all droplets, liquid and Ice I guess
                    _,R_liq = NR_monodisperse(N_0  ,q.liq)
                    _,R_ice = NR_monodisperse(N_INP,q.ice)
                else
                    error("Unsupported size distribution closure (N_r_closure): $(N_r_closure)")
                end

                base = 1/(4*π*D) # as q goes up, R goes 
                τ_liq = base / (N_0 * R_liq)
                τ_ice = base / (N_INP * R_ice)
                # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster...
                # println("effective τ_liq = ",τ_liq_eff, " effective τ_ice = ",τ_ice_eff, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )


            elseif use_supersat == :linear_combination
                ρ_l = FT(1000.) # kg m^-3, CLIMAParameters default for cloud_liquid
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice
                r_r = 20 * 1e-6 # 20 microns
                r_0 = .2 * 1e-6 # .2 micron base aerosol
                N_l = FT(1e-5 / (4/3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q_liq.. (N = N_r in homogenous)
                N_i = FT(1e-7 / (4/3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q_ice... (N = N_r + N_0)
                #
                c_1l = get(param_set.user_aux, :linear_combination_liq_c_1, FT(N_l * r_0)  ) # I think at q=0, we need c_1 from linear = c_1 from geometric...
                c_2l = get(param_set.user_aux, :linear_combination_liq_c_2, FT(2/3) ) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
                c_3l = get(param_set.user_aux, :linear_combination_liq_c_3, FT(-1e-10)) # asssume nothing here? (keep 0 as upper bound?) 
                τ_liq = 1/( c_1l * exp(c_2l*(T-273.15) + c_3l*q.liq))
                #
                c_1i = get(param_set.user_aux, :linear_combination_ice_c_1, FT(N_i * r_0)  ) # I think at q=0, we need c_1 from linear = c_1 from geometric...
                c_2i = get(param_set.user_aux, :linear_combination_ice_c_2, FT(2/3) ) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
                c_3i = get(param_set.user_aux, :linear_combination_ice_c_3, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
                τ_ice = 1/( c_1i * exp(c_2i*(T-273.15) + c_3i*q.ice))

            elseif use_supersat == :linear_combination_with_w
                ρ_l = FT(1000.) # kg m^-3, CLIMAParameters default for cloud_liquid
                ρ_i = FT(916.7) # CLIMAParameters default for cloud_ice
                r_r = 20 * 1e-6 # 20 microns
                r_0 = .2 * 1e-6 # .2 micron base aerosol
                N_l = FT(1e-5 / (4/3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q_liq.. (N = N_r in homogenous)
                N_i = FT(1e-7 / (4/3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q_ice... (N = N_r + N_0)
                #
                c_1l = get(param_set.user_aux, :linear_combination_liq_c_1, FT(N_l * r_0)  ) # I think at q=0, we need c_1 from linear = c_1 from geometric...
                c_2l = get(param_set.user_aux, :linear_combination_liq_c_2, FT(2/3) ) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
                c_3l = get(param_set.user_aux, :linear_combination_liq_c_3, FT(0)) # asssume nothing here? (keep 0 as upper bound?) 
                c_4l = get(param_set.user_aux, :linear_combination_liq_c_4, FT(0)) # start at 0
                τ_liq = 1/( c_1l * exp(c_2l*(T-273.15) + c_3l*q.liq + c_4l*w))
                #
                c_1i = get(param_set.user_aux, :linear_combination_ice_c_1, FT(N_i * r_0)  ) # I think at q=0, we need c_1 from linear = c_1 from geometric...
                c_2i = get(param_set.user_aux, :linear_combination_ice_c_2, FT(2/3) ) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
                c_3i = get(param_set.user_aux, :linear_combination_ice_c_3, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
                c_4i = get(param_set.user_aux, :linear_combination_ice_c_4, FT(0)) # start at 0
                τ_ice = 1/( c_1i * exp(c_2i*(T-273.15) + c_3i*q.ice + c_4i*w))


            elseif use_supersat == :neural_network
                neural_network_params       = param_set.user_aux.neural_microphysics_relaxation_network # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
                neural_network_params       = Float32.(collect(neural_network_params)) # collect from ntuple, then convert to Float32 for NN since that's what it's supposed to be (save eltype in the jld2?)
                model_x_0_characteristic    = get(param_set.user_aux, :model_x_0_characteristic, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)

                TC = TurbulenceConvection
                if isdefined(TC, :neural_network) # do this so we don't have to read from disk and recreate the network evry time bc that's super slow
                    TC.neural_network = TC.re(neural_network_params) # set the parameters to the ones we just read in                       
                else
                    model_re_location           = get(param_set.user_aux, :model_re_location, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)
                    model_re_location           = isa(model_re_location, Val) ? string(typeof(model_re_location).parameters[1]) : model_re_location # extract the value from the Val type (we put it in there in parameter_set.jl to make it isbits)
                    re                          = model_destructure_re_from_file(model_re_location) # get the reconstruction function from the file ( this is probably slow, we could pass the string repr)
                    neural_network              = vec_to_NN(neural_network_params, re) # construct the NN from the parameters
                    TC.re                       = re # store it in the TC so we don't have to reconstruct it every time
                    TC.neural_network           = neural_network # store it in the TC so we don't have to reconstruct it every time
                end
                τ_liq, τ_ice = predict_τ(ρ, T, q, w, TC.neural_network; norm=model_x_0_characteristic) # pass in the NN and get the τs out

                τ_liq = min(τ_liq, FT(1e10)) # ensure we don't just get a huge number all the time
                τ_ice = min(τ_ice, FT(1e10)) # ensure we don't just get a huge number all the time
            else
                error("Unsupported supersaturation type: $(use_supersat)")
            end

            # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster... change to be dt related)

            min_τ_liq = get(param_set.user_aux, :min_τ_liq, FT(1)) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
            min_τ_ice = get(param_set.user_aux, :min_τ_ice, FT(1)) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)

            # @info("τ_liq", (τ_liq, min_τ_liq)) 
            # @info("τ_ice", (τ_ice, min_τ_ice)) 

            # if (0 < τ_liq) && (τ_liq < min_τ_liq) # fast source (convert to chaining operators)
            #     τ_liq = min_τ_liq
            # elseif (0 < τ_ice) && (τ_ice < min_τ_ice) # fast source
            #     τ_ice = min_τ_ice
            # else
            # end
            τ_liq = max(min_τ_liq, τ_liq) # make sure it's at least some value
            τ_ice = max(min_τ_ice, τ_ice) # make sure it's at least some value

            q_eq = TD.PhasePartition(q.tot, TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()), TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())) # all 3 are vapor amounts
            S_ql = (q_vap - q_eq.liq) / τ_liq # | microphys_params.τ_cond_evap | CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, TD.PhasePartition(FT(0),q_vap,FT(0)))  
            S_qi = (q_vap - q_eq.ice) / τ_ice # -(source to vapor) = source to condensate
        elseif use_korolev_mazin
            # need to get w into here somewhere...
            S_ql,S_qi = korolev_mazin_2007(param_set, area, ρ, Δt, ts, w)
            q_eq = TD.PhasePartition(q.tot, TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()), TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())) # all 3 are vapor amounts
            # println("effective τ_liq = ",(q_vap - q_eq.liq)/S_ql, " effective τ_ice = ",(q_vap - q_eq.ice)/S_qi, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )

            # korolev_mazin fix effective tau > 1 (or whatever timestep was realy but we run w/ 1 mostly rn)
            τ_liq_eff = (q_vap - q_eq.liq)/S_ql 
            τ_ice_eff = (q_vap - q_eq.ice)/S_qi
            # println("effective τ_liq = ",τ_liq_eff, " effective τ_ice = ",τ_ice_eff, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )
            if     (0 < τ_liq_eff) && (τ_liq_eff < 1) # fast source
                # @show("rate limiting cond: ", S_ql, (q_vap - q_eq.liq) / 1)
                S_ql = (q_vap - q_eq.liq) / 1
            elseif (0 < τ_ice_eff) && (τ_ice_eff < 1) # fast source
                # @show("rate limiting dep: ", S_qi, (q_vap - q_eq.ice) / 1)
                S_qi = (q_vap - q_eq.ice) / 1
            elseif (-1 < τ_liq_eff) && (τ_liq_eff < 0) # fast sink
                # @show("rate limiting evap: ", S_ql, (q_vap - q_eq.liq) / 1)
                S_ql = (q_vap - q_eq.liq) / -1
            elseif (-1 < τ_ice_eff) && (τ_ice_eff < 0) # fast sink
                # @show("rate limiting sub: ", S_qi, (q_vap - q_eq.ice) / 1)
                S_qi = (q_vap - q_eq.ice) / -1
            else
            end
            
        else # basic noneq (no supersat formulation so not likely to be right)
            # TODO - is that the state we want to be relaxing to?
            ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q.tot)
            q_eq = TD.PhasePartition(thermo_params, ts_eq)

            S_ql = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, q) 
            S_qi = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, ice_type, q_eq, q)
        end

        
        # melting
        if T >= thermo_params.T_freeze # could go insisde supersat bc base eq is already 0 above freezng
            S_qi = min(0,S_qi) # allow sublimation but not new ice formation ( or should we just melt and then let evaporation work? )
            
            # send existing ice to liquid
            S_qi -= q.ice/Δt # send any existing ice to liquid
            S_ql += q.ice/Δt # send any existing ice to liquid
        # homogenous freezing
        elseif (T < thermo_params.T_icenuc)
            if (S_ql > 0)
                S_qi += S_ql # any vapor coming to liquid goes to ice instead (smoother in total condensate than setting it to 0 suddenly?)
                S_ql = 0
            end
            if q.liq > 0
                S_ql -= q.liq/Δt # send any existing liquid to ice (maybe make this a rate later)
                S_qi += q.liq/Δt # send any existing liquid to ice (maybe make this a rate later)
            end
        end

        # LIMITER
        if !(use_supersat==false) # might need to do these first bc ql,qc tendencies maybe are applied individually and can still crash the code if one is too large...

            # Don't let depletion exceed condensate (esp after melting and homogenous freezing adjustments above) -- do this first before growth limiter since growth through vapor is potentially limited by condensate consumption given no direct liquid-ice pathway
            if S_ql < FT(0)
                S_ql = -min(-S_ql, q.liq / (2*Δt)) # limit to liquid amount (2*Δt for stability thoug hthat implies you can never fully deplete the condensate lol, just Δt doesn't work)
            end
            if S_qi < FT(0)
                S_qi = -min(-S_qi, q.ice / (2*Δt))
            end

            # Don't let growth exceed vapor (or existing vapor plus vapor gain from condensate)
            S = max(0, S_ql) + max(0,S_qi) # only add if positive sources to condensate
            Qv = q_vap / (2 * Δt) # should this be 2*Δt for stability like with condensate? we still have some crashes...

            # Too much vapor consumption (note in a supersaturation context, that's a little bit ridiculous to worry about since supersaturation is usually single digit percent at most and timesteps are short-ish)
            if S > ( Qv - min(0,S_ql) - min(0,S_qi)  ) # only add positive sources to vapor (i.e. subtract if negative S_q)
                if (S_qi > 0) && (S_ql > 0)
                        S_ql *= Qv/S
                        S_qi *= Qv/S
                elseif (S_qi > 0) && (S_ql < 0)
                    S_qi *= (Qv - S_ql)/S # source to ice not to exceed vapor plus addition from liquid... (S=S_qi here) # (are these stable?, theyre potentially big leaps)
                elseif (S_qi < 0) && (S_ql > 0)
                    S_ql *= (Qv - S_qi)/S # source to liquid not to exceed vapor plus addition from ice... (S=S_ql here) # (are these stable?, theyre potentially big leaps)
                end # otherwise we have the negative limiters below for sublimation, evaporation... if both are neg that's sufficient....
            end

            # # NOTE - TEST - DONT LOET EITHER S_QL OR S_QI APPROACH QVAP ON THEIR OWN ( in case one gets processed before the other -- but then why did it fail before?")
            # if S_ql >= FT(0)
            #     S_ql = min(S_ql, Qv)
            # end
            # if S_qi >= FT(0)
            #     S_qi = min(S_qi, Qv)
            # end

            

        else # let other limiters do their thing...
            # TODO - handle limiters elswhere (note, these break our allowance for liquid to ice compensation, no?)
            if S_ql >= FT(0)
                S_ql = min(S_ql, q_vap / Δt)
            else
                S_ql = -min(-S_ql, q.liq / Δt)
            end
            if S_qi >= FT(0)
                S_qi = min(S_qi, q_vap / Δt)
            else
                S_qi = -min(-S_qi, q.ice / Δt)
            end
        end

        ql_tendency += S_ql
        qi_tendency += S_qi
        # if area < 1e-2
            # @info( "current", Qv, S_ql, S_qi, q.liq / Δt, q.ice / Δt, T, area)
        # end
    end
    return NoneqMoistureSources{FT}(ql_tendency, qi_tendency)
end

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function precipitation_formation(
    param_set::APS,
    precip_model::AbstractPrecipitationModel,
    rain_formation_model::AbstractRainFormationModel,
    qr::FT,
    qs::FT,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts,
    precip_fraction,
) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)

    microphys_params = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    qt_tendency = FT(0)
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    qr_tendency = FT(0)
    qs_tendency = FT(0)
    θ_liq_ice_tendency = FT(0)

    α_acnv = TCP.microph_scaling_acnv(param_set)
    α_accr = TCP.microph_scaling_accr(param_set)

    if area > 0

        q = TD.PhasePartition(thermo_params, ts)

        Π_m = TD.exner(thermo_params, ts)
        c_pm = TD.cp_m(thermo_params, ts)
        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)
        I_i = TD.internal_energy_ice(thermo_params, ts)
        I = TD.internal_energy(thermo_params, ts)

        if precip_model isa Clima0M
            qsat = TD.q_vap_saturation(thermo_params, ts)
            λ = TD.liquid_fraction(thermo_params, ts)

            S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(microphys_params, q, qsat))

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1 - λ)
            qt_tendency += S_qt
            ql_tendency += S_qt * λ
            qi_tendency += S_qt * (1 - λ)
            θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))
        end

        if precip_model isa Clima1M
            T = TD.air_temperature(thermo_params, ts)
            T_fr = TCP.T_freeze(param_set)
            c_vl = TCP.cv_l(param_set)
            c_vm = TD.cv_m(thermo_params, ts)
            Rm = TD.gas_constant_air(thermo_params, ts)
            Lf = TD.latent_heat_fusion(thermo_params, ts)

            # TODO - limiters and positivity checks should be done elsewhere
            qr = max(qr, FT(0)) / precip_fraction
            qs = max(qs, FT(0)) / precip_fraction

            # Autoconversion of cloud ice to snow is done with a simplified rate.
            # The saturation adjustment scheme prevents using the
            # 1-moment snow autoconversion rate that assumes
            # that the supersaturation is present in the domain.
            if rain_formation_model isa Clima1M_default
                S_qt_rain = -min(q.liq / Δt, α_acnv * CM1.conv_q_liq_to_q_rai(microphys_params, q.liq))
            elseif rain_formation_model isa Clima2M
                S_qt_rain =
                    -min(
                        q.liq / Δt,
                        α_acnv * CM2.conv_q_liq_to_q_rai(
                            microphys_params,
                            rain_formation_model.type,
                            q.liq,
                            ρ,
                            N_d = rain_formation_model.prescribed_Nd,
                        ),
                    )
            else
                error("Unrecognized rain formation model")
            end
            S_qt_snow = -min(q.ice / Δt, α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice))
            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            ql_tendency += S_qt_rain
            qi_tendency += S_qt_snow

            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

            # accretion cloud water + rain
            if rain_formation_model isa Clima1M_default
                S_qr =
                    min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ)) *
                    precip_fraction
            elseif rain_formation_model isa Clima2M
                if rain_formation_model.type isa CMT.LD2004Type
                    S_qr =
                        min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ)) *
                        precip_fraction
                elseif rain_formation_model.type isa CMT.KK2000Type || rain_formation_model.type isa CMT.B1994Type
                    S_qr =
                        min(
                            q.liq / Δt,
                            α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr, ρ),
                        ) * precip_fraction
                elseif rain_formation_model.type isa CMT.TC1980Type
                    S_qr =
                        min(
                            q.liq / Δt,
                            α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr),
                        ) * precip_fraction
                else
                    error("Unrecognized 2-moment rain formation model type")
                end
            else
                error("Unrecognized rain formation model")
            end
            qr_tendency += S_qr
            qt_tendency -= S_qr
            ql_tendency -= S_qr
            θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

            # accretion cloud ice + snow
            S_qs =
                min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, snow_type, q.ice, qs, ρ)) *
                precip_fraction
            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

            # sink of cloud water via accretion cloud water + snow
            S_qt =
                -min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, snow_type, q.liq, qs, ρ)) *
                precip_fraction
            if T < T_fr # cloud droplets freeze to become snow)
                qs_tendency -= S_qt
                qt_tendency += S_qt
                ql_tendency += S_qt
                θ_liq_ice_tendency -= S_qt / Π_m / c_pm * Lf * (1 + Rm / c_vm)
            else # snow melts, both cloud water and snow become rain
                α::FT = c_vl / Lf * (T - T_fr)
                qt_tendency += S_qt
                ql_tendency += S_qt
                qs_tendency += S_qt * α
                qr_tendency -= S_qt * (1 + α)
                θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1 + Rm / c_vm) * α - L_v0)
            end

            # sink of cloud ice via accretion cloud ice - rain
            S_qt =
                -min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, rain_type, q.ice, qr, ρ)) *
                precip_fraction
            # sink of rain via accretion cloud ice - rain
            S_qr = -min(qr / Δt, α_accr * CM1.accretion_rain_sink(microphys_params, q.ice, qr, ρ)) * precip_fraction
            qt_tendency += S_qt
            qi_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)

            # accretion rain - snow
            if T < T_fr
                S_qs =
                    min(qr / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, snow_type, rain_type, qs, qr, ρ)) *
                    precip_fraction *
                    precip_fraction
            else
                S_qs =
                    -min(qs / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, rain_type, snow_type, qr, qs, ρ)) *
                    precip_fraction *
                    precip_fraction
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
        end
    end
    return PrecipFormation{FT}(θ_liq_ice_tendency, qt_tendency, ql_tendency, qi_tendency, qr_tendency, qs_tendency)
end
