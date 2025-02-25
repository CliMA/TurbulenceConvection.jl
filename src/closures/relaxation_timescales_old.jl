"""
Get timescales for relaxation of the supersaturation
"""
# ======================================================================================================================================== #
D_func(T::FT, p::FT) where {FT} = FT((2.11 * 1e-5) * (T / FT(273.15))^1.94 * (p / 101325))::FT  # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226


function get_τ(
    param_set::APS,
    microphys_params::ACMP,
    nonequilibrium_moisture_scheme::Symbol,
    q::TD.PhasePartition,
    T::FT,
    p::FT,
    ρ::FT,
    w::FT,
    z::FT,
) where {FT}


    # ===== Are these annotations faster? ===== #
    local τ_liq::FT
    local τ_ice::FT

    local ρ_l::FT
    local ρ_i::FT
    local D::FT
    local r_0::FT
    local r_i::FT
    local r_l::FT
    local r_is::FT
    local N::FT

    local q_i_0::FT
    local q_0::FT
    # ========================================= #


    if nonequilibrium_moisture_scheme === :Base
        tau_weights = TC.get_isbits_nt_nt(param_set.user_params, :tau_weights, nothing) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
        if isnothing(tau_weights)
            τ_liq = CMNe.τ_relax(microphys_params, liq_type)
            τ_ice = CMNe.τ_relax(microphys_params, ice_type)
        else
            liq_params, ice_params = tau_weights.liq.liq_params, tau_weights.ice.ice_params # gets used in eval below
            τ_liq, τ_ice = 10 .^ liq_params.log10_tau_liq, 10 .^ ice_params.log10_tau_ice # log fcn hand implementation...
        end

    elseif nonequilibrium_moisture_scheme === :exponential_T_scaling_ice
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        # TODO: Drop the constants from this and just subsume them into N (makes picking an initial c_1, c_2 harder though so maybe not? Also is much harder for more complicated N_INP expressions)
        ρ_i = CMP.ρ_cloud_ice(microphys_params)
        r_is = CMP.r_ice_snow(microphys_params)

        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)

        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        # r_i = FT((q.ice / ((4/3) * π * N_i * ρ_i))^(1/3)) # (the mass-diameter relationship is poorly defined anyway for ice crystals) -- if q.ice is 0 this goes to 0... making it hard to generate ice
        r_i = r_from_qN(q.ice, N_i, ρ_i)
        r_i = max(r_i, r_0) # bound to be at least ~micron size...something like kohler crit radius


        τ_ice = FT(1 / (4 * π * D * N_i * r_i))
        τ_liq = CMNe.τ_relax(microphys_params, liq_type)

    elseif nonequilibrium_moisture_scheme === :exponential_T_scaling_ice_raw
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        # this one is just the exponential part of the T scaling ice, no N or r_i complications
        q_0 = FT(1e-6)
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        # r_0 = FT((q_i_0 / ((4/3) * π * N_0 * ρ_i))^(1/3)) # (the mass-diameter relationship is poorly defined anyway for ice crystals) -- if q.ice is 0 this goes to 0... making it hard to generate ice
        # r_0 = r_from_qN(q_0, N_0, ρ_i)
        r_0 = FT(20.0 * 10^-6) # 20 micron
        # derived from typical values and assume q_i = 4/3 * π * r^3 * ρ_i * N_i and N_i = c_1 e^(c_2 T)

        c_1 = TC.get_isbits_nt_nt(
            param_set.user_params,
            :exponential_T_scaling_ice_raw_c_1,
            FT(4 * π * D * (q_i_0 / (4 / 3) * π * ρ_i)^(1 / 3) * (0.02)^(2 / 3)),
        ) # Fletcher 1962 (values taken from Frostenberg 2022)
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_ice_raw_c_2, FT(-0.6 * 2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022)
        τ_ice = FT(1 / (c_1 * exp(c_2 * (T - T_fr))))
        τ_liq = CMNe.τ_relax(microphys_params, liq_type)

    elseif (nonequilibrium_moisture_scheme === :T_scaling_ice_F23) || (nonequilibrium_moisture_scheme === :powerlaw_T_scaling_ice) # this should have a gentler curve though the cloud
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        if T >= T_fr
            τ_ice = FT(Inf)
        else
            ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
            r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
            N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)

            # r_i = FT((q.ice / ((4/3) * π * N_i * ρ_i))^(1/3)) # (the mass-diameter relationship is poorly defined anyway for ice crystals) -- if q.ice is 0 this goes to 0... making it hard to generate ice
            r_i = r_from_qN(q.ice, N_i, ρ_i)
            r_i = max(r_i, r_0) # bound to be at least ~micron size...something like kohler crit radius

            τ_ice = FT(1 / (4 * π * D * N_i * r_i))


        end
        τ_liq = CMNe.τ_relax(microphys_params, liq_type)

    elseif nonequilibrium_moisture_scheme === :exponential_times_powerlaw_scaling_ice # Demott 2010
        error("NotImplmentedError: This nonequilibrium_moisture_scheme functionality has not been implemented yet")


    elseif nonequilibrium_moisture_scheme === :geometric_liq__geometric_ice # scaling on q that impacts liquid
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice 
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)   # only for prior

        # c_1l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        # c_2l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3)) # Halfway between 1/3 and 1
        # c_3l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_l = r_from_qN(q.liq, N_l, ρ_l, r_min = r_0)
        τ_liq = FT(1 / (4 * π * D * N_l * r_l))
        # τ_liq = FT(1 / (4 * π * D * (c_1l * q.liq^(c_2l) + c_3l))) # let  be Nr = c_1 * q^(c_2) + c_3

        # ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        # c_1i = TC.get_isbits_nt_nt(param_set.user_params, :geometric_ice_c_1, FT(1 / (4 / 3 * π * ρ_i * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        # c_2i = TC.get_isbits_nt_nt(param_set.user_params, :geometric_ice_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
        # c_3i = TC.get_isbits_nt_nt(param_set.user_params, :geometric_ice_c_3, FT(N_i0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0)
        τ_ice = FT(1 / (4 * π * D * N_i * r_i))
        # τ_ice = FT(1 / (4 * π * D * (c_1i * q.ice^(c_2i) + c_3i))) # let  be Nr = c_1 * q^(c_2) + c_3

    elseif nonequilibrium_moisture_scheme === :geometric_liq__exponential_T_scaling_ice #
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        r_r = FT(20.0 * 10^-6) # 20 micron
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        c_1l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        c_2l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
        c_3l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_3, FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_i))) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0)
        τ_liq = FT(1 / (4 * π * D * N_l * r_l))
        # τ_liq = FT(1 / (4 * π * D * (c_1l * q.liq^(c_2l) + c_3l))) # let  be Nr = c_1 * q^(c_2)

        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0)
        τ_ice = FT(1 / (4 * π * D * N_i * r_i))

    elseif nonequilibrium_moisture_scheme === :geometric_liq__powerlaw_T_scaling_ice # scaling on q that impacts liquid and ice
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)

        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol

        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

        # liq 
        # c_1l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_0^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        # c_2l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
        # c_3l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0*r_0)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0)
        τ_liq = FT(1 / (4 * π * D * N_l * r_l))
        # τ_liq = FT(1 / (4 * π * D * (c_1l * q.liq^(c_2l) + c_3l))) # let  be Nr = c_1 * q^(c_2)

        # ice
        if T >= T_fr
            τ_ice = FT(Inf)
        else
            N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
            r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0)
            τ_ice = FT(1 / (4 * π * D * N_i * r_i))

        end




    elseif nonequilibrium_moisture_scheme === :geometric_liq__exponential_T_scaling_and_geometric_ice # scaling on q that impacts liquid and ice
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)

        ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice

        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

        # c_1l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        # c_2l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
        # c_3l = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0)
        τ_liq = FT(1 / (4 * π * D * N_l * r_l))
        # τ_liq = FT(1 / (4 * π * D * (c_1l * q.liq^(c_2l) + c_3l))) # let  be Nr = c_1 * q^(c_2) 

        c_1i = get(
            param_set.user_params,
            :exponential_T_scaling_and_geometric_ice_c_1,
            FT((4 * π * D) * ((4 / 3 * π * ρ_i)^(-1 / 3) * (N_i0)^(2 / 3) * (0.02)^(2 / 3) + (N_i0 * r_0))),
        ) # Yeahhh.... idk for this one lol... just combined them serially from the homogenous case where c_3 is -1/3, and used .02 as the prefactor
        # c_2i = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1 -- should this be the same as c_2g? It's the same mixing... 
        # c_3i = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_3, FT((4 * π * D) * r_0 * 0.02)) # Fletcher 1962 (values taken from Frostenberg 2022) and used .02 as the prefactor
        # c_4i = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_4, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)


        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0)
        τ_ice = 1 / (4 * π * D * N_i * r_i)
        # τ_ice = 1 / ((c_1i * q.ice^(c_2i) * exp(c_2i * c_4i * (T - T_fr)) + c_3i) * exp(c_4i * (T - T_fr)))
        # τ_ice = FT(1 / ((c_1i * q.ice^(c_2i) + c_3i) * exp(c_4i * (T - T_fr))))

        # not sure what N_i would be here

        # τ_liq = min(τ_liq, FT(1e10)) # ensure we don't just get a huge number all the time
        # τ_ice = min(τ_ice, FT(1e10)) # ensure we don't just get a huge number all the time

    elseif nonequilibrium_moisture_scheme === :linear_combination
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)

        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
        #
        # c_1l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_1, FT(N_l0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        # c_2l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_2, FT(-1e-10)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        # c_3l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_3, FT(2 / 3)) # asssume nothing here? (keep 0 as upper bound?) 

        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0)
        τ_liq = FT(1 / (4 * π * D * N_l * r_l))
        # τ_liq = FT(1 / ( exp(c_1l + c_2l * (T - T_fr) + c_3l * log10(q.liq))))
        #
        # c_1i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_1, FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        # c_2i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        # c_3i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_3, FT(2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0)
        τ_ice = FT(1 / (4 * π * D * N_i * r_i))
        # τ_ice = FT(1 / ( exp(c_1i + c_2i * (T - T_fr) + c_3i * log10(q.ice))))

        

    elseif nonequilibrium_moisture_scheme === :linear_combination_with_w
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
        w_0 = FT(1e-3) # 1 mm/s
        #
        # c_1l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_1, FT(N_l0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        # c_2l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        # c_3l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_3, FT(-10)) # asssume nothing here? (keep 0 as upper bound?) 
        # c_4l = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_4, FT(0)) # start at 0
        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0)
        τ_liq = FT(1 / (4 * π * D * N_l * r_l))
        # τ_liq = FT(1 / ( exp( c_1l + c_2l * (T - T_fr) + c_3l * log10(q.liq) + c_4l * w/w_0)))
        #
        # c_1i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_1, FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        # c_2i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        # c_3i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_3, FT(-10)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
        # c_4i = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_4, FT(0)) # start at 0
        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0)
        τ_ice = FT(1 / (4 * π * D * N_i * r_i))
        # τ_ice = FT(1 / (exp(c_1i + c_2i * (T - T_fr) + c_3i * log10(q.ice) + c_4i * w/w_0)))


    elseif nonequilibrium_moisture_scheme === :neural_network
        neural_network_params = param_set.user_params.neural_microphysics_relaxation_network # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
        neural_network_params = Float32.(collect(neural_network_params)) # collect from ntuple, then convert to Float32 for NN since that's what it's supposed to be (save eltype in the jld2?)
        model_x_0_characteristic = TC.get_isbits_nt_nt(param_set.user_params, :model_x_0_characteristic, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)

        TC = TurbulenceConvection
        if isdefined(TC, :neural_network) # do this so we don't have to read from disk and recreate the network evry time bc that's super slow
            TC.neural_network = TC.re(neural_network_params) # set the parameters to the ones we just read in                       
        else
            model_re_location = TC.get_isbits_nt_nt(param_set.user_params, :model_re_location, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)
            if !isnothing(model_re_location)
                model_re_location = string(model_re_location) # convert symbol to string...
            end

            re = model_destructure_re_from_file(model_re_location) # get the reconstruction function from the file ( this is probably slow, we could pass the string repr)
            neural_network = vec_to_NN(neural_network_params, re) # construct the NN from the parameters
            TC.re = re # store it in the TC so we don't have to reconstruct it every time
            TC.neural_network = neural_network # store it in the TC so we don't have to reconstruct it every time
        end
        τ_liq, τ_ice = predict_τ(ρ, T, q, w, TC.neural_network; norm = model_x_0_characteristic) # pass in the NN and get the τs out

        τ_liq = min(τ_liq, FT(1e10)) # ensure we don't just get a huge number all the time
        τ_ice = min(τ_ice, FT(1e10)) # ensure we don't just get a huge number all the time   


    elseif nonequilibrium_moisture_scheme === :raymond_ice_test # should be closer to DeMott 2015 but...
        T_fr = TCP.T_freeze(param_set)

        N_r_closure = TC.get_isbits_nt_nt(param_set.user_args, :N_r_closure, :monodisperse)
        N_0 =
            TC.get_isbits_nt_nt(param_set.user_params, :N_0, nothing) |>
            x -> (isnothing(x) ? TC.get_isbits_nt_nt(param_set.user_args, :N_0, FT(100 * 10e6)) : pyinterp(z, x.z, x.values)) # if defined in aux, interp at current z, otherwise default to user_args value, otherwise default to a value
        # @info("using $nonequilibrium_moisture_scheme supersaturation formulation")
        # N_m = TC.get_isbits_nt_nt(param_set.user_args, :N_m, FT(-0.2)) # log slope https://doi.org/10.1073/pnas.1514034112
        # N_b = TC.get_isbits_nt_nt(param_set.user_args, :N_b, FT(-5 - T_fr * N_m)) # -5 - T_fr * N_m
        # N_INP = FT(10^(N_m * T + N_b) * 10^3) # per liter to per m^3
        # -- hopefully this gets around the problem of draining water vapor at initialization of clouds but also allows speedup as droplets grow (assuming fixed drop concenctration)
        # R    =  max(((q.liq + q.ice)/(4/3*π*ρ_l*N_0))^(1/3), FT(0.2*10^-6)) # bound to be at least ~micron size...something like kohler crit radius

        N_l = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
        N_i = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)

        if N_r_closure === :inhomogeneous
            N_l, R_liq = NR_inhomogeneous_mixing_liquid(thermo_params, N_l, p, q.liq, ts_LCL) # testing inhomogeneous mixing would have r fixed and then let N vary... set r based on adiabatic extrapolation from cloud base 
            _, R_ice = NR_monodisperse(N_i, q.ice)
            # maybe look into using a mixed model where N/R are partly towards the inhomogenous value depending on the true entrainment/mixing params... see literature on this
            # NOTE, ON DYCOMS ADIABATIC R W/ ORIGINAL N_0 WORKED BETTER... HMMM (though that's not using DYCOMS N)

        elseif N_r_closure === :monodisperse # uniform size for all droplets, liquid and Ice I guess
            _, R_liq = NR_monodisperse(N_l, q.liq)
            _, R_ice = NR_monodisperse(N_i, q.ice)
        else
            error("Unsupported size distribution closure (N_r_closure): $(N_r_closure)")
        end

        base = FT(1 / (4 * π * D)) # as q goes up, R goes 
        τ_liq = base / (N_l * R_liq)
        τ_ice = base / (N_i * R_ice)
        # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster...
        # println("effective τ_liq = ",τ_liq_eff, " effective τ_ice = ",τ_ice_eff, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )

    else
        error("Unsupported supersaturation type: $(nonequilibrium_moisture_scheme)")
    end


    # Athough we have greatly improved stability with our timestepper, extremely fast timescales can still cause oscillations in S_ql, S_qi from morrison_milbrandt_2015_style, and ceven worse can still cause unrecoverable buoyancy shockwaves in the model... (not sure why, you can imagine massive heating etc, particularly during spinup...
    # The use of an absolute τ q limiter is thus nearly unavoidable...
    # This poses a trade-off between stabilizing the model (esp during spinup) and discouraging bad parameters in calibration...

    # while we should perhaps converge to saturation adjustment, the fact of the matter is we do not...

    # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster... change to be dt related) [used to have FT(1)]
    min_τ_liq::FT = TC.get_isbits_nt_nt(param_set.user_params, :min_τ_liq, eps(FT)) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
    min_τ_ice::FT = TC.get_isbits_nt_nt(param_set.user_params, :min_τ_ice, eps(FT)) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
    τ_liq = max(τ_liq, min_τ_liq)
    τ_ice = max(τ_ice, min_τ_ice)

    max_τ_liq::FT = TC.get_isbits_nt_nt(param_set.user_params, :max_τ_liq, Inf) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
    max_τ_ice::FT = TC.get_isbits_nt_nt(param_set.user_params, :max_τ_ice, Inf) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)

    τ_liq = min(τ_liq, max_τ_liq)
    τ_ice = min(τ_ice, max_τ_ice)

    return τ_liq, τ_ice
end



function get_N_i(param_set::APS, nonequilibrium_moisture_scheme::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT}

    thermo_params::TDPS = TCP.thermodynamics_params(param_set)

    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts)

    return get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
end

function get_N_l(param_set::APS, nonequilibrium_moisture_scheme::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT}

    thermo_params::TDPS = TCP.thermodynamics_params(param_set)

    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts)


    return get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
end


function get_Ns(param_set::APS, nonequilibrium_moisture_scheme::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT}

    # FT = eltype(param_set)

    thermo_params::TDPS = TCP.thermodynamics_params(param_set)

    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts)

    return get_Ns(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
end


function get_N_i(param_set::APS, nonequilibrium_moisture_scheme::Symbol, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}

    # ===== Are these annotations faster? ===== #
    local τ_liq::FT
    local τ_ice::FT

    local ρ_l::FT
    local ρ_i::FT
    local D::FT
    local r_0::FT
    local r_i::FT
    local r_is::FT

    local q_i_0::FT
    local q_0::FT

    # local N_i::Union{FT ,Nothing}
    local N_i::FT # testing for type stability

    microphys_params::ACMP = TCP.microphysics_params(param_set)
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    # ========================================= #



    if nonequilibrium_moisture_scheme === :Base
        # N_i = nothing
        N_i = FT(NaN) # testing NaN over nothing for type stability
    elseif nonequilibrium_moisture_scheme === :exponential_T_scaling_ice
        T_fr = TCP.T_freeze(param_set)
        # TODO: Drop the constants from this and just subsume them into N (makes picking an initial c_1, c_2 harder though so maybe not? Also is much harder for more complicated N_INP expressions)
        # ρ_i = CMP.ρ_cloud_ice(microphys_params)
        # q_is_0 = FT(ρ_i * 4/3 * π * CMP.r_ice_snow(microphys_params)^3) # r_is radius ice sphere --
        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_ice_c_1, FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_ice_c_2, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)
        N_i = c_1 * exp(c_2 * (T - T_fr))

        # if TC.get_isbits_nt_nt(param_set.user_params, :adjust_ice_N, false)
        #     N_i = max(N_i, q.ice / q_is_0) # make sure existing ice contributes to the size, particularly w/ sedimentation
        # end            

    elseif nonequilibrium_moisture_scheme === :exponential_T_scaling_ice_raw
        # N_i = nothing
        N_i = FT(NaN) # testing NaN over nothing for type stability

    elseif (nonequilibrium_moisture_scheme === :T_scaling_ice_F23) || (nonequilibrium_moisture_scheme === :powerlaw_T_scaling_ice) # this should have a gentler curve though the cloud
        # ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        # q_is_0 = FT(ρ_i * 4/3 * π * CMP.r_ice_snow(microphys_params)^3) # r_is radius ice sphere --
        T_fr = TCP.T_freeze(param_set)
        # c_3i = TC.get_isbits_nt_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_3, -c_2i) # F23 (values taken from Frostenberg 2022) -- is redundant w/ c_1i
        if T >= T_fr
            N_i = FT(0) # cant do powerlaw otherwise...
        else
            c_1 = TC.get_isbits_nt_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_1, FT(-9)) # F23 (values taken from Frostenberg 2022)
            c_2 = TC.get_isbits_nt_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_2, FT(9)) # F23 (values taken from Frostenberg 2022)
            N_i = (10^c_1) * (-(T - T_fr))^c_2
        end

        # if TC.get_isbits_nt_nt(param_set.user_params, :adjust_ice_N, false)
        #     N_i = max(N_i, q.ice / q_is_0) # make sure existing ice contributes to the size, particularly w/ sedimentation
        # end

    elseif nonequilibrium_moisture_scheme === :exponential_times_powerlaw_scaling_ice # Demott 2010
        error("NotImplementedError: This nonequilibrium_moisture_scheme functionality has not been implemented yet")

    elseif nonequilibrium_moisture_scheme === :geometric_liq__geometric_ice # scaling on q that impacts liquid
        # error("NotImplmentedError: This nonequilibrium_moisture_scheme functionality has not been implemented yet. We predict Nr but not N")

        # if monodisperse/N fixed, then the default of q_thresh is probably fine, we're not accounting for T
        # N_i = nothing
        # N_i = FT(NaN) # testing NaN over nothing for type stability

        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)   # only for prior
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice

        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :geometric_ice_c_1, FT(1 / (4 / 3 * π * ρ_i * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :geometric_ice_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
        c_3 = TC.get_isbits_nt_nt(param_set.user_params, :geometric_ice_c_3, FT(N_i0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

        N_i = c_1 * q.ice^(c_2) + c_3

    elseif nonequilibrium_moisture_scheme === :geometric_liq__exponential_T_scaling_ice #
        # D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        # q_is_0 = FT(ρ_i * 4/3 * π * CMP.r_ice_snow(microphys_params)^3) # r_is radius ice sphere --
        r_0 = FT(20.0 * 10^-6) # 20 micron
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice

        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_ice_c_1, FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_ice_c_2, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)
        N_i = FT(c_1 * exp(c_2 * (T - T_fr)))
        # ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice 

        # if TC.get_isbits_nt_nt(param_set.user_params, :adjust_ice_N, false)
        #     N_i = max(N_i, q.ice / q_is_0) # make sure existing ice contributes to the size, particularly w/ sedimentation
        # end

    elseif nonequilibrium_moisture_scheme === :geometric_liq__powerlaw_T_scaling_ice # scaling on q that impacts liquid and ice
        T_fr = TCP.T_freeze(param_set)
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        # q_is_0 = FT(ρ_i * 4/3 * π * CMP.r_ice_snow(microphys_params)^3) # r_is radius ice sphere --


        if T >= T_fr
            N_i = FT(0) # cant do powerlaw otherwise...
        else
            # ice
            c_1 = TC.get_isbits_nt_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_1, -9) # F23 (values taken from Frostenberg 2022)
            c_2 = TC.get_isbits_nt_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_2, FT(9)) # F23 (values taken from Frostenberg 2022)
            # c_3i = TC.get_isbits_nt_nt(param_set.user_params, :T_scaling_ice_c_3, -c_2i) # F23 (values taken from Frostenberg 2022) -- is redundant w/ c_1i
            N_i = FT((10^c_1) * (-(T - T_fr))^c_2)

            # if TC.get_isbits_nt_nt(param_set.user_params, :adjust_ice_N, false)
            #     N_i = FT(max(N_i, q.ice / q_is_0)) # make sure existing ice contributes to the size, particularly w/ sedimentation
            # end
        end


    elseif nonequilibrium_moisture_scheme === :geometric_liq__exponential_T_scaling_and_geometric_ice # scaling on q that impacts liquid and ice
        T_fr = TCP.T_freeze(param_set)
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

        c_1i = get(
            param_set.user_params,
            :exponential_T_scaling_and_geometric_ice_c_1,
            FT((4 * π * D) * ((4 / 3 * π * ρ_i)^(-1 / 3) * (N_i0)^(2 / 3) * (0.02)^(2 / 3) + (N_i0 * r_0))),
        ) # Yeahhh.... idk for this one lol... just combined them serially from the homogenous case where c_3 is -1/3, and used .02 as the prefactor
        c_2i = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1 -- should this be the same as c_2g? It's the same mixing... 
        c_3i = TC.get_isbits_nt_nt(
            param_set.user_params,
            :exponential_T_scaling_and_geometric_ice_c_3,
            FT((4 * π * D) * r_0 * 0.02),
        ) # Fletcher 1962 (values taken from Frostenberg 2022) and used .02 as the prefactor
        c_4i = TC.get_isbits_nt_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_4, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)

        # Not really clear what is N and what is r with this mixed mode... certainly there could be some threshold effect but it's not clear what it is...
        # I think the T-scaling part of N should be the same as in just exponential_T_scaling... the actual droplet dist is less important for thresh?
        # but not clear if this actually works...
        # N_i = c_3i/(4*π*D*r_0) * exp(c_4i * (T - T_fr))  # but then what is r? 
        N_i = (c_1i * q.ice^(c_2i) + c_3i) * exp(c_4i * (T - T_fr))

    elseif nonequilibrium_moisture_scheme === :linear_combination
        T_fr = TCP.T_freeze(param_set)
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        # N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

        # we could try something like N_i = c_1i * exp(c_2i * c_3i * (T - T_fr)) but no idea on the constants...
        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_1, FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        c_3 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_3, FT(2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
        # N_i = FT(NaN) # testing NaN over nothing for type stability
        N_i = exp(c_1 + c_2 * (T - T_fr) + c_3 * q.ice / 1e-7)

    elseif nonequilibrium_moisture_scheme === :linear_combination_with_w
        T_fr = TCP.T_freeze(param_set)
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        # N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
        w_0 = FT(1e-3) # 1 mm/s

        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_1, FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        c_3 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_3, FT(-10)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
        c_4 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_ice_c_4, FT(0)) # start at 0
        # N_i = FT(NaN) # testing NaN over nothing for type stability
        N_i = exp(c_1 + c_2 * (T - T_fr) + c_3 * q.ice / 1e-7 + c_4 * w / w_0)

    elseif nonequilibrium_moisture_scheme === :neural_network
        # N_i = nothing
        N_i = FT(NaN) # testing NaN over nothing for type stability

    elseif nonequilibrium_moisture_scheme === :raymond_ice_test # should be closer to DeMott 2015 but...
        T_fr = TCP.T_freeze(param_set)

        N_r_closure = TC.get_isbits_nt_nt(param_set.user_args, :N_r_closure, :monodisperse)
        N_0 =
            TC.get_isbits_nt_nt(param_set.user_params, :N_0, nothing) |>
            x -> (isnothing(x) ? TC.get_isbits_nt_nt(param_set.user_args, :N_0, FT(100 * 10e6)) : pyinterp(z, x.z, x.values)) # if defined in aux, interp at current z, otherwise default to user_args value, otherwise default to a value
        # @info("using $nonequilibrium_moisture_scheme supersaturation formulation")
        N_m = TC.get_isbits_nt_nt(param_set.user_args, :N_m, FT(-0.2)) # log slope https://doi.org/10.1073/pnas.1514034112
        N_b = TC.get_isbits_nt_nt(param_set.user_args, :N_b, FT(-5 - T_fr * N_m)) # -5 - T_fr * N_m
        N_i = FT(10^(N_m * T + N_b) * 10^3) # per liter to per m^3

    else
        error("Unsupported supersaturation type: $(nonequilibrium_moisture_scheme)")
    end

    # adjust N_i if desired
    if TC.get_isbits_nt_nt(param_set.user_params, :adjust_ice_N, false)
        N_i = adjust_ice_N(param_set, microphys_params, N_i, q, T, p; r_min=FT(0.2 * 1e-6))
    end

    return N_i
end


function get_N_l(param_set::APS, nonequilibrium_moisture_scheme::Symbol, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}

    # local N_l::Union{FT ,Nothing}
    local N_l::FT

    microphys_params::ACMP = TCP.microphysics_params(param_set)
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)

    if nonequilibrium_moisture_scheme === :raymond_ice_test # should be closer to DeMott 2015 but...
        N_r_closure = TC.get_isbits_nt_nt(param_set.user_args, :N_r_closure, :monodisperse)
        N_l =
            TC.get_isbits_nt_nt(param_set.user_params, :N_0, nothing) |>
            x -> (isnothing(x) ? TC.get_isbits_nt_nt(param_set.user_args, :N_0, FT(100 * 10e6)) : pyinterp(z, x.z, x.values)) # if defined in aux, interp at current z, otherwise default to user_args value, otherwise default to a value

        # -- hopefully this gets around the problem of draining water vapor at initialization of clouds but also allows speedup as droplets grow (assuming fixed drop concenctration)
        # R    =  max(((q.liq + q.ice)/(4/3*π*ρ_l*N_0))^(1/3), FT(0.2*10^-6)) # bound to be at least ~micron size...something like kohler crit radius

        if N_r_closure === :inhomogeneous
            N_l, R_liq = NR_inhomogeneous_mixing_liquid(thermo_params, N_l, p, q.liq, ts_LCL) # testing inhomogeneous mixing would have r fixed and then let N vary... set r based on adiabatic extrapolation from cloud base 
        # maybe look into using a mixed model where N/R are partly towards the inhomogenous value depending on the true entrainment/mixing params... see literature on this
        # NOTE, ON DYCOMS ADIABATIC R W/ ORIGINAL N_0 WORKED BETTER... HMMM (though that's not using DYCOMS N)

        elseif N_r_closure === :monodisperse # uniform size for all droplets, liquid and Ice I guess
        #
        else
            error("Unsupported size distribution closure (N_r_closure): $(N_r_closure)")
        end

    elseif occursin("geometric_liq", string(nonequilibrium_moisture_scheme))
        D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
        T_fr = TCP.T_freeze(param_set)
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)

        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
        c_3 = TC.get_isbits_nt_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

        N_l = c_1 * q.liq^(c_2) + c_3

    elseif occursin("linear_combination", string(nonequilibrium_moisture_scheme))
        T_fr = TCP.T_freeze(param_set)
        ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
        r_r = FT(20 * 1e-6) # 20 microns
        r_0 = FT(0.2 * 1e-6) # .2 micron base aerosol
        N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
        w_0 = FT(1e-3) # 1 mm/s

        c_1 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_1, FT(N_l0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
        c_2 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
        c_3 = TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_3, FT(-10)) # asssume nothing here? (keep 0 as upper bound?) 
        c_4 =
            nonequilibrium_moisture_scheme === :linear_combination_with_w ?
            TC.get_isbits_nt_nt(param_set.user_params, :linear_combination_liq_c_4, FT(0)) : FT(0) # start at 0

        N_l = exp(c_1 + c_2 * (T - T_fr) + c_3 * q.liq + c_4 * w / w_0)

    else
        # N_l = nothing
        N_l = FT(NaN) # testing NaN over nothing for type stability
    end
    return N_l

end

function get_Ns(param_set::APS, nonequilibrium_moisture_scheme::Symbol, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}

    N_l::FT = get_N_l(param_set, nonequilibrium_moisture_scheme, q, T, p, w)
    N_i::FT = get_N_i(param_set, nonequilibrium_moisture_scheme, q, T, p, w)

    return (; liq = N_l, ice = N_i)
end
