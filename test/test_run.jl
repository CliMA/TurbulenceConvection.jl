## Run test on cases and plot relative to LES

using Random
debug = false

FT = Float64
forcing_type = :obs_data

flight_numbers = 10
# flight_numbers = 9
# flight_numbers = [1, 9, 10, 12, 13]
# flight_numbers = 1


# convert flight_number to array if it's not already
if isa(flight_numbers, Int)
    flight_numbers = [flight_numbers]
end




if forcing_type === :obs_data
    forcing_str = "Obs"
elseif forcing_type === :ERA5_data
    forcing_str = "ERA5"
else
    error("forcing_type must be :obs_data or :ERA5_data")
end


reload_environment = false
if reload_environment || !isdefined(Main, :TurbulenceConvection) || !isdefined(Main, :main1d)
    using Pkg
    # Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/")) # this break precomilation of TC.jl for some reason...
    # Pkg.develop(path = expanduser("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl"))
    # Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
    # Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/"))
    # Pkg.develop(path = expanduser("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl"))


    # using Revise
    # Pkg.resolve()
    # Pkg.instantiate()

    using Revise
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
    using Revise
    # Pkg.resolve()
    # Pkg.instantiate()

    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests"))
    using Revise
    # Pkg.resolve()
    # Pkg.instantiate()
    using Revise
    using TurbulenceConvection
    tc = pkgdir(TurbulenceConvection)
    # include("/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/src/TurbulenceConvection.jl")
    # tc = "/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl"
    # includet(joinpath(tc, "driver", "main_revise.jl"))
    # includet(joinpath(tc, "driver", "generate_namelist.jl"))
    includet(joinpath(tc, "driver", "main.jl"))
    includet(joinpath(tc, "driver", "generate_namelist.jl"))

end

# ------------------ #
reload_param_set = false
if reload_param_set | !isdefined(Main, :param_set)
    using Pkg
    using Revise
    import CLIMAParameters as CP # use CLIMAParameters = "0.7, 0.8, 0.9, 0.10"
    # import ClimaParams as CPP # would using this trouble w/ TC.jl? it's a different uuid technically..., use ClimaParams = "0.10"
    import Thermodynamics as TD
    import Thermodynamics.Parameters as TDP
    # FT = Float64

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias") # CP 0.7 and below, Thermodynamics 0.11 and above
    aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
end
# ------------------ #

for flight_number in flight_numbers

    setup_environment = true # leave true
    if setup_environment
        case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * forcing_str * "_data"

        if !@isdefined(NameList) # For if we didn't use reload_enviornment (I'm not sure tbf if this also checks for modules lol, I guess we'll find out)
            @info("NameList not defined by module, defining now as simple Dict()") # This allows plotting even if you haven't gone through the environment loading by allowing you to get the simname/uuid etc in namelist 
            namelist = Dict{String, Union{Dict, NamedTuple}}(
                "time_stepping" => Dict(),
                "stats_io" => Dict(),
                "microphysics" => Dict(),
                "user_args" => Dict(),
                "user_params" => Dict(),
                "turbulence" => Dict("EDMF_PrognosticTKE" => Dict()),
                "thermodynamics" => Dict(),
                "output" => Dict(),
                "meta" => Dict("simname" => case_name, "uuid" => ""),
            )
        else
            @info("NameList defined by module, using NameList.default_namelist()")
            namelist = NameList.default_namelist(
                case_name;
                root = "./",  # controls where the namelist gets written if write=true
                # write = parsed_args["write_namelist"], # controls if the namelist gets written (I'm not sure if this is repetitive -- does it get written again after the namelist overwrite? I think this is the wrong namelist that kept getting written when I ran things earlier separate from the results...)
                write = false, # don't write the namelist because we only care to write the later one once we've done our overwrites since that's the one that actually runs, though idk where that one gets written tbh...
            )
            namelist["user_params"] = Dict()
        end



        # namelist["time_stepping"]["dt_min"] = 1.0
        # namelist["time_stepping"]["dt_max"] = 2.0
        namelist["time_stepping"]["dt_min"] = 20.0
        namelist["time_stepping"]["dt_max"] = 10.0
        namelist["time_stepping"]["adapt_dt"] = true
        namelist["time_stepping"]["use_tendency_timestep_limiter"] = true
        namelist["time_stepping"]["N_dt_max_edmf_violate_dt_min"] = 1e4
        # namelist["time_stepping"]["spinup_dt_factor"] = 0.2
        # namelist["time_stepping"]["spinup_half_t_max"] = 3600.0 * 1.0 # 1 hour
        namelist["time_stepping"]["t_max"] = 3600.0 * 7.0
        namelist["grid"]["dz_min"] = 5.0 * (2 * namelist["time_stepping"]["dt_min"] / 0.5)
        namelist["time_stepping"]["allow_spinup_adapt_dt"] = true
        # namelist["stats_io"]["frequency"] = 600.0
        namelist["stats_io"]["frequency"] = 60.0

        # namelist["time_stepping"]["dt_min"] = 0.1
        # namelist["time_stepping"]["dt_max"] = 1.0
        namelist["stats_io"]["frequency"] = 30.0
        # namelist["time_stepping"]["t_max"] = 3600.0 * 3.5

        # nonequilibrium_moisture_scheme = :Base
        nonequilibrium_moisture_scheme = :geometric_liq__powerlaw_T_scaling_ice
        # namelist["microphysics"]["τ_sub_dep"] = 100000.0
        # namelist["microphysics"]["τ_cond_evap"] = 1.0
        # namelist["user_args"] = (; nonequilibrium_moisture_scheme = nonequilibrium_moisture_scheme, τ_use = :morrison_milbrandt_2015_style, use_sedimentation = false, grid_mean_sedimentation = false)
        # namelist["user_args"] = (;nonequilibrium_moisture_scheme=nonequilibrium_moisture_scheme, τ_use=:morrison_milbrandt_2015_style_exponential_part_only) 
        # namelist["user_args"] = (;nonequilibrium_moisture_scheme=nonequilibrium_moisture_scheme, τ_use=:standard) 




        # namelist["user_params"]["initial_profile_updraft_area"] = 0.33 # i think either dict or named tuple is fine, we dispatch depending... # doesn't seem to work anyway...


        # namelist["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.6
        # ("turbulence", "EDMF_PrognosticTKE", "max_area", FT(.3)), # stability limiting...
        namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area_bc"] = "Prognostic" # unstable w/o  setting other params (didn't help lol, was 0 at sfc)
        # 


        # surface area  
        namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.45) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)
        namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = FT(1e-10) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = FT(0.3) # testing
        namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.3

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = 2 # testing lol
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = [0.25, 0.1] # testing lol, maybe they'll combine nicely


        # namelist["turbulence"]["EDMF_PrognosticTKE"]["limit_min_area"] = true # testing
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area_limiter_scale"] = 5 # testing strong detrainment...
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area_limiter_power"] = 2000 # testing strong detrainment...


        # # testing broad strong detrainment throughout (kills updraft)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 5 # testing strong detrainment...
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 4 # too 

        # # testing broad weak detrainment throughout (still killed very low)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 1 # testing strong detrainment...
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 3 # too 

        # # testing broad very weak detrainment throughout (leads to instability no matter what updraft area is seemingly)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = .1 # testing strong detrainment...
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 10 # too

        # # testing weak detrainment until just before limit (kills updraft very low)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 20 # testing strong detrainment...
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 30 # too 

        # testing mild detrainment until through limit then strong (runs slower bc updraft?) (big spike at top but at least makes it there -- runs super slow in adaptive mode)...
        namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 1 # testing strong detrainment...
        namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 10 # too 


        # careful, we're getting qt + precip right but not qt, probably bc we're raining/snowing out our cloud lol... (the snow timescale is fast)
        namelist["microphysics"]["τ_acnv_rai"] = 25000.0
        namelist["microphysics"]["τ_acnv_sno"] = 10000.0
        namelist["microphysics"]["q_liq_threshold"] = 1e-4
        namelist["microphysics"]["q_ice_threshold"] = 1e-5
    
        # ========================================================================================================================= #

        # nonequilibrium_moisture_scheme = :exponential_T_scaling_ice
        nonequilibrium_moisture_scheme = :Base
        # nonequilibrium_moisture_scheme = :geometric_liq__powerlaw_T_scaling_ice

        # nonequilibrium_moisture_scheme = :geometric_liq__exponential_T_scaling_and_geometric_ice

        # calibrated params (taken from a best calibration, but unstable for some reason lol (without sedimentation) --(probably some other params I didn't set)
        # also they're probably not ideal w/ sedimentation anyway bc they'd be too slow at the top for falling sediment...
        # namelist["user_params"]["T_scaling_ice_c_1"] = 0.0013737960245654272
        # namelist["user_params"]["T_scaling_ice_c_2"] = -2.6671518848963918
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 2.0261087615972144
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 0.015498615940217064
        # namelist["microphysics"]["q_liq_threshold"] = 5.096061819670483e-06
        # namelist["microphysics"]["q_ice_threshold"] = 1.2852512726660437e-07
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.030642886967713
        # namelist["microphysics"]["τ_acnv_rai"] = 58765.07428891631
        # namelist["microphysics"]["τ_acnv_sno"] =  53529.186891693455
        # namelist["microphysics"]["τ_cond_evap"] = 0.07565030997988167


        # namelist["user_params"]["powerlaw_T_scaling_ice_c_1"] = -2.933615846075211 
        # namelist["user_params"]["powerlaw_T_scaling_ice_c_2"] = 5.358678935403965
        namelist["user_params"]["powerlaw_T_scaling_ice_c_1"] = -9.
        namelist["user_params"]["powerlaw_T_scaling_ice_c_2"] = 9.
        namelist["user_params"]["ice_sedimentation_scaling_factor"] = 3.2942321652880064 
        namelist["user_params"]["rain_sedimentation_scaling_factor"] = 1.0813243749787707 
        namelist["user_params"]["snow_sedimentation_scaling_factor"] = 1.030624346659287 

        namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 35.52149505301267 
        namelist["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 0.20274259096885316 
        namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] =  0.03176056443055847 

        namelist["microphysics"]["τ_acnv_rai"] = 16.44967299681959 
        namelist["microphysics"]["τ_acnv_sno"] = 31.87217008160495 
        namelist["microphysics"]["q_liq_threshold"] = 4.273820398645035e-8 
        namelist["microphysics"]["q_ice_threshold"] = 1.6602807618941017e-6 
        namelist["microphysics"]["τ_cond_evap"]  = 5.684678316404151 
    
        

        # nonequilibrium_moisture_scheme = :powerlaw_T_scaling_ice
        namelist["user_params"]["T_scaling_ice_c_1"] = 2e-6 * 4
        namelist["user_params"]["T_scaling_ice_c_2"] = 7

        # sedimentation_integration_method = :upwinding # leads to spike at bottom but...
        sedimentation_integration_method = :right_biased
        namelist["user_args"] = (;
            nonequilibrium_moisture_scheme = nonequilibrium_moisture_scheme,
            # nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style,
            nonequilibrium_moisture_sources_limiter_type = :none,
            use_sedimentation = true, 
            grid_mean_sedimentation = false,
            sedimentation_integration_method = sedimentation_integration_method,
            use_heterogeneous_ice_nucleation = false,
            sedimentation_ice_number_concentration = nonequilibrium_moisture_scheme,
            sedimentation_liq_number_concentration = nonequilibrium_moisture_scheme, # idk if this is good or bad lol...
            liq_terminal_velocity_scheme = :Chen2022Vel,
            ice_terminal_velocity_scheme = :Chen2022Vel,
            rain_terminal_velocity_scheme = :Chen2022Vel,
            snow_terminal_velocity_scheme = :Chen2022Vel,
            )

            

        namelist["user_params"]["ice_sedimentation_Dmax"] = 62.5e-6 # 62.5 microns cutoff from CM

        namelist["user_params"]["ice_sedimentation_Dmax"] = Inf # 62.5 microns cutoff from CM
        
        namelist["user_params"]["adjust_ice_N"] = true
        namelist["microphysics"]["r_ice_snow"] = 62.5e-6


        # namelist["user_params"]["min_τ_liq"] = 1.0
        # namelist["user_params"]["min_τ_ice"] = 1.0

        # ========================================================================================================================= #
        # ========================================================================================================================= #
        reload_CalibrateEDMF = false
        if reload_CalibrateEDMF | !isdefined(Main, :CalibrateEDMF)
            Pkg.activate(expanduser("~/Research_Schneider/CliMA/CalibrateEDMF.jl/"))
            using CalibrateEDMF
            Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
        end
        path_to_Costa_SOTA = "/groups/esm/cchristo/cedmf_results/james_v1_runs/results_Inversion_p22_e300_i15_mb_LES_2024-03-15_10-08_Vxo_longer_long_run/Diagnostics.nc"
        CEDMF_dir = pkgdir(CalibrateEDMF.ModelTypes) # any submodule should work the same
        include(joinpath(CEDMF_dir, "tools", "DiagnosticsTools.jl")) # provides optimal_parameters()
        Costa_SOTA = optimal_parameters(path_to_Costa_SOTA, method = "last_nn_particle_mean")
        Costa_SOTA = Dict(zip(Costa_SOTA...)) # turn to dict

        NUM_NN_PARAMS = 12 # number of neural network parameters from Costa
        Costa_SOTA_linear_ent_params = [Costa_SOTA["linear_ent_params_{$(i)}"] for i in 1:NUM_NN_PARAMS]

        # Entrainment Params
        namelist["turbulence"]["EDMF_PrognosticTKE"]["linear_ent_params"] = Costa_SOTA_linear_ent_params
        # Diffusion Params
        namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = Costa_SOTA["tke_ed_coeff"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = Costa_SOTA["tke_diss_coeff"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = Costa_SOTA["static_stab_coeff"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] = Costa_SOTA["tke_surf_scale"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] = Costa_SOTA["Prandtl_number_0"]

        # Momentum Exchange parameters
        namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = Costa_SOTA["pressure_normalmode_adv_coeff"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = Costa_SOTA["pressure_normalmode_buoy_coeff1"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = Costa_SOTA["pressure_normalmode_drag_coeff"]

        # Area Limiters
        namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area_limiter_scale"] = Costa_SOTA["min_area_limiter_scale"]
        namelist["turbulence"]["EDMF_PrognosticTKE"]["min_area_limiter_power"] = Costa_SOTA["min_area_limiter_power"]
        # ========================================================================================================================= #
        # ========================================================================================================================= #




        # namelist["thermodynamics"]["potential_temperature_reference_pressure"] = 985. * 100 # try changing this...

        # namelist["thermodynamics"]["moisture_model"] = "equilibrium"
        namelist["thermodynamics"]["moisture_model"] = "nonequilibrium"
        namelist["thermodynamics"]["sgs"] = "mean"

        # ========================================================================================================================= #
        # ========================================================================================================================= #
        # ========================================================================================================================= #
        # search_locator__1
        using JLD2
        original_namelist = namelist
        # param_results_path = "/central/groups/esm/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/slurm/param_results/jl_bzlskt/param_results.jld2"
        # setup =  JLD2.load(param_results_path)["param_results"]
        # namelist = setup["namelist"]
        # namelist["user_params"] = namelist["user_aux"] # backport fix
        # dt = 0.5
        # dt = 5.0
        dt = 2.0
        # dt = 20.0
        # namelist["stats_io"]["frequency"] = dt
        # namelist["stats_io"]["frequency"] = dt * 10
        namelist["stats_io"]["frequency"] = 600.
        namelist["stats_io"]["frequency"] = 60.
        # namelist["stats_io"]["frequency"] = dt
        # namelist["user_params"]["ice_sedimentation_scaling_factor"] = 1.0

        namelist["user_params"]["liq_ice_collision_efficiency"] = 1.0 * 0.0
        namelist["user_params"]["liq_ice_collision_scaling_factor"] = 1.0 * 0.0

        namelist["meta"]["simname"] = original_namelist["meta"]["simname"]
        namelist["stats_io"]["calibrate_io"] = false
        # namelist["user_params"]["print_taus"] = false
        # namelist["user_params"]["print_sources"] = false

        # namelist["user_params"]["min_τ_liq"] = 100.
        # namelist["user_params"]["min_τ_ice"] = 200.

        namelist["user_params"]["mix_τ_liq"] = 1e20
        namelist["user_params"]["max_τ_ice"] = 1e20
        namelist["time_stepping"]["algorithm"] = "Euler";  namelist["time_stepping"]["dt_min"] = dt ; namelist["time_stepping"]["dt_max"] = 2*dt; # jacobian too slow
        # namelist["time_stepping"]["algorithm"] = "ImplicitEuler"; namelist["time_stepping"]["dt_min"] = 1e-10 ; namelist["time_stepping"]["dt_max"] = 100.;  namelist["time_stepping"]["abstol"] = 1e-3;  namelist["time_stepping"]["reltol"] = 1e-1 # jacobian too slow
        # namelist["time_stepping"]["algorithm"] = "KenCarp47"; namelist["time_stepping"]["dt_min"] = 1e-10 ; namelist["time_stepping"]["dt_max"] = 100.;  namelist["time_stepping"]["abstol"] = 1e-3;  namelist["time_stepping"]["reltol"] = 1e-1 # jacobian too slow
        # namelist["time_stepping"]["algorithm"] = "Heun" ; namelist["time_stepping"]["dt_min"] = 1e-10 ; namelist["time_stepping"]["dt_max"] = 100.;  namelist["time_stepping"]["abstol"] = 1e-3;  namelist["time_stepping"]["reltol"] = 1e-1
        # namelist["time_stepping"]["algorithm"] = "RK4"; namelist["time_stepping"]["dt_min"] = eps(FT) ; namelist["time_stepping"]["dt_max"] = 20.
        # namelist["time_stepping"]["algorithm"] = "Tsit5"; namelist["time_stepping"]["dt_min"] = 1e-10 ; namelist["time_stepping"]["dt_max"] = 20.;  namelist["time_stepping"]["abstol"] = 1e-3;  namelist["time_stepping"]["reltol"] = 1e-1

        # namelist["time_stepping"]["adaptive_depth_limit"] = -1
        # namelist["time_stepping"]["dt_min"] = dt
        # namelist["time_stepping"]["dt_max"] = 2*dt

        namelist["time_stepping"]["spinup_dt_factor"] = 0.25
        namelist["time_stepping"]["spinup_half_t_max"] = 3600.0 * 1.0 # 1 hour

        namelist["time_stepping"]["t_max"] = 3600.0 * 7.0

        namelist["user_params"]["ice_sedimentation_scaling_factor"] = 3.2942321652880064
        # namelist["user_params"]["ice_sedimentation_scaling_factor"] = 3.2942321652880064e3
        namelist["user_params"]["q_min"] = 0. #eps(FT)

        # namelist["grid"]["dz_min"] = 50.0
        namelist["grid"]["dz_min"] = 5.0 * (2 * namelist["time_stepping"]["dt_min"] / 0.5)
        # namelist["grid"]["dz_min"] = 5.0 * (2 * 0.5 / 0.5)


        namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.9) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)

        namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_type"] = "fractional" # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)

        namelist["time_stepping"]["use_tendency_timestep_limiter"] = false  # testing
        namelist["time_stepping"]["N_dt_max_edmf_violate_dt_min"] = 1e3
        namelist["time_stepping"]["allow_cfl_dt_max_violate_dt_min"] = false   
        namelist["time_stepping"]["use_fallback_during_spinup"] = true
        namelist["time_stepping"]["allow_spinup_adapt_dt"] = false

        # namelist["time_stepping"]["limit_tendencies_using_dt_min_factor"] = 1.0 # margin of safety 



        nonequilibrium_moisture_scheme = :neural_network
        calibration_setup = "tau_autoconv_noneq"
        global pkg_dir = pkgdir(CalibrateEDMF)
        global main_experiment_dir = joinpath(pkg_dir, "experiments", "SOCRATES")
        global experiment_dir = joinpath(pkg_dir, "experiments", "SOCRATES", "subexperiments", "SOCRATES_"*string(nonequilibrium_moisture_scheme))
        global calibration_setup = "tau_autoconv_noneq"
        global const no_constraint = CalibrateEDMF.KalmanProcessUtils.EnsembleKalmanProcesses.no_constraint
        include("/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/subexperiments/SOCRATES_neural_network/Calibrate_and_Run/tau_autoconv_noneq/calibration_parameters_and_namelist.jl")
        namelist["user_params"]["neural_microphysics_relaxation_network"] = calibration_parameters__experiment_setup["neural_microphysics_relaxation_network"]["prior_mean"]
        namelist["user_params"]["model_re_location"] = nn_path
        namelist["user_params"]["model_x_0_characteristic"] = FT.(nn_pretrained_x_0_characteristic)
        namelist["thermodynamics"]["moisture_model"] = "nonequilibrium"



        # nonequilibrium_moisture_scheme = :geometric_liq__powerlaw_T_scaling_ice
        # nonequilibrium_moisture_scheme = :Base
        # namelist["microphysics"]["τ_sub_dep"] =  1e6
        # namelist["microphysics"]["τ_cond_evap"] = 100

        # nonequilibrium_moisture_scheme = :exponential_T_scaling_ice
        # calibration_setup = "tau_autoconv_noneq"
        # namelist["thermodynamics"]["moisture_model"] = "nonequilibrium"
        # namelist["user_params"]["exponential_T_scaling_ice_c_1"] = 0.049
        # namelist["user_params"]["exponential_T_scaling_ice_c_2"] = -0.49


        nonequilibrium_moisture_scheme = :geometric_liq__exponential_T_scaling_and_geometric_ice
        

        
        namelist["user_args"] = (; 
            nonequilibrium_moisture_scheme = nonequilibrium_moisture_scheme,

            # nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style,
            # nonequilibrium_moisture_sources_limiter_type = :none,
            # nonequilibrium_moisture_sources_limiter_type = :basic,

            # nonequilibrium_moisture_sources_limiter_type = :none, # can go past supersaturation but tendency limited timesteps prevents
            # entr_detr_limiter_type = :none,
            # precipitation_tendency_limiter_type = :none,
            # default_tendency_limiter_type = :none,

            # nonequilibrium_moisture_sources_limiter_type = :standard_supersaturation,
            # nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style,
            nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style_exponential_part_only,
            entr_detr_limiter_type = :basic,
            precipitation_tendency_limiter_type = :basic,
            default_tendency_limiter_type = :basic,

            # fallback_nonequilibrium_moisture_sources_limiter_type = :basic, # can go past supersaturation but tendency limited timesteps should prevent...
            # fallback_nonequilibrium_moisture_sources_limiter_type = :standard_supersaturation, # has limiter issues...
            # fallback_nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style,
            fallback_nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style_exponential_part_only,
            fallback_entr_detr_limiter_type = :basic,
            fallback_precipitation_tendency_limiter_type = :basic,
            fallback_default_tendency_limiter_type = :basic,

            use_sedimentation = true, 
            grid_mean_sedimentation = false,
            sedimentation_integration_method = sedimentation_integration_method,
            use_heterogeneous_ice_nucleation = false,
            sedimentation_ice_number_concentration = nonequilibrium_moisture_scheme,
            sedimentation_liq_number_concentration = nonequilibrium_moisture_scheme, # idk if this is good or bad lol...
            liq_terminal_velocity_scheme = :Chen2022Vel,
            ice_terminal_velocity_scheme = :Chen2022Vel,
            rain_terminal_velocity_scheme = :Chen2022Vel,
            snow_terminal_velocity_scheme = :Chen2022Vel,
            )
            




        global debug = true # allow to run multiple simulaitons at once if we're just debugging and wan't to go faster, puts them in random subdirectories in target location with random uuid suffix
        # ========================================================================================================================= #
        # ========================================================================================================================= #
        # ========================================================================================================================= #

        

        if debug
            Random.seed!() # switch back to random entropy
            namelist["meta"]["uuid"] = randstring(4) # allow to run multiple simulaitons at once if we're just debugging and wan't to go faster
        else
            namelist["meta"]["uuid"] = "" # no uuid
        end
        uuid = string(namelist["meta"]["uuid"])
        simname = namelist["meta"]["simname"]
        @info("simname", simname)
        if debug
            namelist["output"]["output_root"] = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/test/debug")
        else
            namelist["output"]["output_root"] = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/test")
        end
        outpath = joinpath(namelist["output"]["output_root"], "Output.$simname.$uuid") # no uuid
        simulation_outpath = joinpath(outpath, "stats/Stats." * case_name * ".nc")
        @info("outpath", outpath)


        # output_relpath   =  "Output."*case_name*".out/stats/Stats."*case_name*".nc" # could try Glob.jl, see https://discourse.julialang.org/t/rdir-search-recursive-for-files-with-a-given-name-pattern/75605/20
        # simulation_outpath = joinpath(namelist["output"]["output_root"], output_relpath), 

        # @info(namelist)

    end

    
    run_simulation = true
    if run_simulation
        @info("Namelist: ", namelist)
        Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/")) # make sure were in this for revise maybe?
        main1d(namelist)
    end


    SOCRATES_truth_path = "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/Reference/Atlas_LES/"
    truth_file = joinpath(
        SOCRATES_truth_path,
        "RF" * string(flight_number, pad = 2) * "_" * forcing_str * "_SOCRATES_128x128_100m_10s_rad10_vg_M2005_aj.nc",
    )

    SOCRATES_input_path = "/home/jbenjami/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/Data/Atlas_LES_Profiles/Input_Data/"

    LES_forcing_str = forcing_str == "Obs" ? "obs" : "ERA5"
    LES_forcing_suffix = forcing_str == "Obs" ? ".nc" : "_mar18_2022.nc"
    LES_input_file =
        "RF" * string(flight_number, pad = 2) * "_" * LES_forcing_str * "-based_SAM_input" * LES_forcing_suffix



    make_plots = true
    # ============================================================================================================================================================= #
    # create figure
    simul_ind = 10000 # time index to plot
    # ============================================================================================================================================================= #

    if make_plots
        # Read output file from SCM run and plot


        using Dierckx
        function pyinterp(x, xp, fp; bc = "error")
            spl = Dierckx.Spline1D(xp, fp; k = 1, bc = bc)
            return spl(vec(x))
        end

        using Plots
        ENV["GKSwstype"]="nul"
        using NCDatasets
        # create directory if doesn't exist
        if !isdir(joinpath(outpath, "Figures"))
            mkdir(joinpath(outpath, "Figures"))
        end

        truth_data = NCDatasets.Dataset(truth_file, "r")
        simul_data = NCDatasets.Dataset(simulation_outpath, "r")
        z_simul = simul_data.group["profiles"]["zc"][:] # should be same for both but the truth one doesn't match the file atlas gave...
        z_simul_f = simul_data.group["profiles"]["zf"][:] # should be same for both but the truth one doesn't match the file atlas gave...

        z_truth = truth_data["z"][:] # should be same for both but the truth one doesn't match the file atlas gave...
        z = z_truth

        lt = length(simul_data.group["timeseries"]["t"])
        if lt < simul_ind
            @warn(
                "Simulation index $simul_ind is out of bounds for simulation data of time index length $lt, reverting to $lt"
            )
            simul_ind = lt
        end


        t_simul = simul_data.group["timeseries"]["t"][:]
        t_simul_ind = t_simul[simul_ind]

        t_truth = truth_data["time"][:]
        t_truth = (t_truth .- t_truth[1]) .* (24 * 3600) # make it relative to start time
        truth_ind = argmin(abs.(t_truth .- t_simul_ind))
        t_truth_ind = t_truth[truth_ind]


        LES_in_data = NCDatasets.Dataset(joinpath(SOCRATES_input_path, LES_input_file), "r") # technically there's some time editing we should do for this but we can put it off for obs


        # ============================================================================================================================================================= #
        # test box
        # import CLIMAParameters as CP
        # # # import Thermodynamics as TD
        # TC = TurbulenceConvection # set if loaded driver/main files
        # CM = TC.CM
        # TD = TC.TD
        # TDP = TC.TD.Parameters

        # FT = Float64
        # toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
        # aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
        # param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
        # thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
        # param_set = create_parameter_set(namelist, toml_dict)

        # aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
        # param_pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
        # # append :thermo_params => thermo_params to the param_pairs
        # param_pairs = vcat(param_pairs, :thermo_params => thermo_params)
        # TP = typeof(thermo_params)
        # microphys_params = CM.Parameters.CloudMicrophysicsParameters{FT, TP}(; param_pairs..., thermo_params)


        # Tg_simul_start = simul_data.group["timeseries"]["Tsurface"][1]
        # Tg_simul_start = simul_data.group["profiles"]["temperature_mean"][:,1][1]
        # pg = simul_data.group["reference"]["p_f"][:][1] # no timeseries
        # ρg = TD.air_density(thermo_params,Tg_simul_start,pg)
        # q_sfc = TD.q_vap_saturation_generic(thermo_params, Tg_simul_start, ρg, TD.Liquid()) # surface specific humidity
        # q_sfc = TD.q_vap_saturation_from_pressure(thermo_params,0. , pg, Tg_simul_start, TD.PhaseEquil) # surface specific humidity

        # ============================================================================================================================================================= #


        # figure out simulation time we're plotting and what time that is in the truth_data (start should be same...)

        ### temperature ###
        data_simul = simul_data.group["profiles"]["temperature_mean"][:, simul_ind]
        data_truth = truth_data["TABS"][:, truth_ind]

        ymax = 4000.0
        # ymax = nothing
        ymax = 4800
        ymin = -100.0
        ymaxs = Dict(1 => 4800, 9 => 4800, 10 => 4800, 12 => 3000, 13 => 3000)
        dpi = 600
        minorgrid = true

        # --------------------------------------------------------------------------------------------- #

        # temperature
        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :topright,
            dpi = dpi,
            marker = :circle,
            markersize = 0.5,
            markerstrokewidth = 0.2,
            xlabel = "T (K)",
            ylabel = "Height (m)",
            title = "Temperature (K)",
            ylim = (ymin, ymaxs[flight_number]),
            minorgrid = minorgrid,
        )
        p_s = plot!(
            data_simul,
            z_simul,
            label = "Simulation",
            marker = :circle,
            markersize = 0.5,
            markerstrokewidth = 0.2,
        ) # add simulation data to same plot
        # vertical line at 5.22 + 273.15 degrees
        # surface temp dot 
        plot!(
            simul_data.group["timeseries"]["Tsurface"][[simul_ind]],
            [0],
            label = "Surface Temp",
            seriestype = :scatter,
            markersize = 5,
            markercolor = :red,
        )
        plot!(
            truth_data["SST"][[truth_ind]],
            [0],
            label = "SST",
            seriestype = :scatter,
            markersize = 5,
            markercolor = :blue,
        )



        savefig(joinpath(outpath, "Figures", "temperature.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # temperature diff
        p_t = plot(
            data_simul .- pyinterp(z_simul, z_truth, data_truth, bc = "extrapolate"),
            z_simul,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "T diff (K)",
            ylabel = "Height (m)",
            title = "Temperature (K) diff",
            minorgrid = minorgrid,
        )
        savefig(joinpath(outpath, "Figures", "temperature_diff.png")) # save to file

        # --------------------------------------------------------------------------------------------- #
        # temperature on pressure grid
        data_simul = simul_data.group["profiles"]["temperature_mean"][:, simul_ind]
        data_truth = truth_data["TABS"][:, truth_ind]

        p_truth = truth_data["PRES"][:, truth_ind]
        p_truth_alt = truth_data["p"][:]

        p_t = plot(
            data_truth,
            p_truth,
            label = "Truth",
            legend = :topright,
            dpi = dpi,
            marker = :circle,
            markersize = 0.4,
            markerstrokewidth = 0.05,
            xlabel = "T (K)",
            ylabel = "Pressure (Pa)",
            title = "Temperature (K)",
            yflip = true,
            ylim = (500, 1020),
            minorgrid = minorgrid,
        )
        p_s = plot!(
            data_simul,
            simul_data.group["reference"]["p_c"] / 100.0,
            label = "Simulation",
            marker = :circle,
            markersize = 0.4,
            markerstrokewidth = 0.05,
        ) # add simulation data to same plot

        p_s = plot!(
            data_truth,
            p_truth_alt,
            label = "Truth alt",
            dpi = dpi,
            marker = :circle,
            markersize = 0.25,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :blue,
        )


        input_T_L = LES_in_data["T"][1, 1, :, 1][:]
        input_p_L = LES_in_data["lev"][:] / 100

        input_T_L = input_T_L[input_p_L .> minimum(p_truth)]
        input_p_L = input_p_L[input_p_L .> minimum(p_truth)]

        p_s = plot!(
            input_T_L,
            input_p_L,
            label = "Input",
            dpi = dpi,
            marker = :circle,
            markersize = 0.75,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :red,
        )

        savefig(joinpath(outpath, "Figures", "temperature_p.png")) # save to file


        # --------------------------------------------------------------------------------------------- #

        ### qt_mean ###
        data_simul = simul_data.group["profiles"]["qt_mean"][:, simul_ind]
        data_truth = truth_data["QT"][:, truth_ind] ./ 1000 # probably is mixing ratio
        data_truth = data_truth ./ (1 .+ data_truth) # convert to specific humidity

        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomleft,
            marker = :circle,
            markersize = 0.5,
            markerstrokewidth = 0.05,
            linewidth = 0.25,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qt (kg/kg)",
            ylabel = "Height (m)",
            title = "QT (kg/kg)",
            minorgrid = minorgrid,
        )



        p_s = plot!(
            data_simul,
            z_simul,
            label = "Simulation",
            marker = :square,
            markersize = 0.5,
            markerstrokewidth = 0.05,
            linewidth = 0.25,
            color = :orange,
        ) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "qt.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # qt diff 
        p_t = plot(
            data_simul .- pyinterp(z_simul, z_truth, data_truth, bc = "extrapolate"),
            z_simul,
            label = "diff",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qt diff (K)",
            ylabel = "Height (m)",
            title = "qt_diff",
            minorgrid = minorgrid,
        )
        savefig(joinpath(outpath, "Figures", "qt_diff.png")) # save to file

        # --------------------------------------------------------------------------------------------- #
        # qt on pressure grid
        data_simul = simul_data.group["profiles"]["qt_mean"][:, simul_ind]
        data_truth = truth_data["QT"][:, truth_ind] / 1000.0

        p_truth = truth_data["PRES"][:, truth_ind]
        p_truth_alt = truth_data["p"][:]

        p_t = plot(
            data_truth,
            p_truth,
            label = "Truth",
            legend = :topright,
            dpi = dpi,
            marker = :circle,
            markersize = 0.4,
            markerstrokewidth = 0.05,
            xlabel = "qt (kg/kg)",
            ylabel = "Pressure (Pa)",
            title = "Temperature (K)",
            yflip = true,
            ylim = (500, 1020),
            minorgrid = minorgrid,
        )
        p_s = plot!(
            data_simul,
            simul_data.group["reference"]["p_c"] / 100.0,
            label = "Simulation",
            marker = :circle,
            markersize = 0.4,
            markerstrokewidth = 0.05,
        ) # add simulation data to same plot

        p_s = plot!(
            data_truth,
            p_truth_alt,
            label = "Truth alt",
            dpi = dpi,
            marker = :circle,
            markersize = 0.25,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :blue,
        )


        input_q_L = LES_in_data["q"][1, 1, :, 1][:]
        input_p_L = LES_in_data["lev"][:] / 100

        input_q_L = input_q_L[input_p_L .> minimum(p_truth) - 50]
        input_p_L = input_p_L[input_p_L .> minimum(p_truth) - 50]

        p_s = plot!(
            input_q_L,
            input_p_L,
            label = "Input",
            dpi = dpi,
            marker = :circle,
            markersize = 0.75,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :red,
        )


        # --------------------------------------------------------------------------------------------- #


        ### ql_mean ### (need to figure out indexing here in particular..., interpolate in t)
        ql_simul_ind = simul_ind
        ql_truth_ind = truth_ind
        # ql_ind = 1
        data_simul = simul_data.group["profiles"]["ql_mean"][:, ql_simul_ind]
        data_truth = truth_data["QCL"][:, ql_truth_ind] ./ 1000

        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "ql (kg/kg)",
            ylabel = "Height (m)",
            title = "QL (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "ql.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        p_t = plot(
            data_simul .- pyinterp(z_simul, z_truth, data_truth, bc = "extrapolate"),
            z_simul,
            label = "diff",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "ql diff (K)",
            ylabel = "Height (m)",
            title = "ql_diff",
            minorgrid = minorgrid,
        )
        savefig(joinpath(outpath, "Figures", "ql_diff.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        ### qi_mean ### (need to figure out indexing here in particular..., interpolate in t)
        qi_simul_ind = simul_ind
        qi_truth_ind = truth_ind
        # ql_ind = 1
        data_simul = simul_data.group["profiles"]["qi_mean"][:, qi_simul_ind]
        data_truth = truth_data["QCI"][:, qi_truth_ind] ./ 1000

        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qi (kg/kg)",
            ylabel = "Height (m)",
            title = "QI (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "qi.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        p_t = plot(
            data_simul .- pyinterp(z_simul, z_truth, data_truth; bc = "extrapolate"),
            z_simul,
            label = "diff",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qi diff (K)",
            ylabel = "Height (m)",
            title = "qi_diff",
            minorgrid = minorgrid,
        )
        savefig(joinpath(outpath, "Figures", "qi_diff.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        ### qc_mean ### (need to figure out indexing here in particular..., interpolate in t)
        qc_simul_ind = simul_ind
        qc_truth_ind = truth_ind
        # ql_ind = 1
        data_simul =
            simul_data.group["profiles"]["ql_mean"][:, qc_simul_ind] +
            simul_data.group["profiles"]["qi_mean"][:, qc_simul_ind]
        data_truth = truth_data["QCL"][:, qc_truth_ind] ./ 1000 + truth_data["QCI"][:, qc_truth_ind] ./ 1000

        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qc (kg/kg)",
            ylabel = "Height (m)",
            title = "QC (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation") # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "qc.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        p_t = plot(
            data_simul .- pyinterp(z_simul, z_truth, data_truth; bc = "extrapolate"),
            z_simul,
            label = "diff",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qc diff (K)",
            ylabel = "Height (m)",
            title = "qc_diff",
            minorgrid = minorgrid,
        )
        savefig(joinpath(outpath, "Figures", "qc_diff.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # liquid potential Temp
        data_simul = simul_data.group["profiles"]["thetal_mean"][:, simul_ind]
        data_truth = truth_data["THETAL"][:, truth_ind]




        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi * 4,
            marker = :circle,
            markersize = 0.5,
            markerstrokewidth = 0.2,
            linewidth = 0.25,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "theta_l (K)",
            ylabel = "Height (m)",
            title = "theta_l (K)",
            minorgrid = minorgrid,
        )
        p_s = plot!(
            data_simul,
            z_simul,
            label = "Simulation",
            marker = :square,
            markersize = 0.5,
            markerstrokewidth = 0.2,
            linewidth = 0.25,
            color = :orange,
        ) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "theta_l.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        p_t = plot(
            data_simul .- pyinterp(z_simul, z_truth, data_truth; bc = "extrapolate"),
            z_simul,
            label = "diff",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "theta_l diff (K)",
            ylabel = "Height (m)",
            title = "thetal_diff",
            minorgrid = minorgrid,
        )
        savefig(joinpath(outpath, "Figures", "theta_l_diff.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # ql and qi combined on separate axes
        data_simul_ql = simul_data.group["profiles"]["ql_mean"][:, simul_ind]
        data_truth_ql = truth_data["QCL"][:, truth_ind] ./ 1000

        data_simul_qi = simul_data.group["profiles"]["qi_mean"][:, simul_ind]
        data_truth_qi = truth_data["QCI"][:, truth_ind] ./ 1000

        p_t = plot(
            data_truth_ql,
            z_truth,
            label = "Truth liq",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            color = :green,
            xlabel = "q (kg/kg)",
            ylabel = "Height (m)",
            title = "QL and QI (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul_ql, z_simul, label = "Simulation liq", color = :lime) # add simulation data to same plot

        # add ice data on same plot with twin x axis
        xaxis2 = Plots.twiny()
        p_s = plot!(
            xaxis2,
            data_truth_qi,
            z_truth,
            label = "Truth ice",
            ylim = (ymin, ymaxs[flight_number]),
            legend = :bottom,
            color = :blue,
            minorgrid = minorgrid,
        )
        p_s = plot!(xaxis2, data_simul_qi, z_simul, label = "Simulation ice", color = :cyan) # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "ql_qi.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # ql, qi updraft/env
        data_simul_ql_up = simul_data.group["profiles"]["updraft_ql"][:, simul_ind]
        data_simul_qi_up = simul_data.group["profiles"]["updraft_qi"][:, simul_ind]
        data_simul_ql_env = simul_data.group["profiles"]["env_ql"][:, simul_ind]
        data_simul_qi_env = simul_data.group["profiles"]["env_qi"][:, simul_ind]

        plot(
            data_simul_ql_up,
            z_simul,
            label = "Updraft liq",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            color = :green,
            xlabel = "q (kg/kg)",
            ylabel = "Height (m)",
            title = "QL and QI (kg/kg)",
            minorgrid = minorgrid,
        )
        plot!(data_simul_ql_env, z_simul, label = "Environment liq", color = :green, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_qi_up, z_simul, label = "Updraft ice", color = :blue) # add simulation data to same plot
        plot!(data_simul_qi_env, z_simul, label = "Environment ice", color = :blue, linestyle = :dash) # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "ql_qi_updraft_env.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # updraft and environmental temperature
        data_simul_up = simul_data.group["profiles"]["updraft_temperature"][:, simul_ind]
        data_simul_env = simul_data.group["profiles"]["env_temperature"][:, simul_ind]
        data_simul = simul_data.group["profiles"]["temperature_mean"][:, simul_ind]

        data_truth = truth_data["TABS"][:, truth_ind]

        p_t = plot(
            data_simul_up,
            z_simul,
            label = "Updraft",
            legend = :topright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            color = :red,
            xlabel = "T (K)",
            ylabel = "Height (m)",
            title = "Updraft and Environmental Temperature (K)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot
        p_s = plot!(data_simul, z_simul, label = "Simulation", color = :brown, xlim = (250, 280)) # add simulation data to same plot

        plot!(
            simul_data.group["timeseries"]["Tsurface"][[simul_ind]],
            [0],
            label = "Surface Temp",
            seriestype = :scatter,
            markersize = 5,
            markercolor = :red,
        )
        plot!(
            truth_data["SST"][[truth_ind]],
            [0],
            label = "SST",
            seriestype = :scatter,
            markersize = 2,
            markercolor = :blue,
        )
        plot!(data_truth, z_truth, label = "Truth mean", color = :green) # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "updraft_env_temp.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # updraft and environmental temperature start and now
        data_simul_up = simul_data.group["profiles"]["updraft_temperature"][:, simul_ind]
        data_simul_env = simul_data.group["profiles"]["env_temperature"][:, simul_ind]
        data_truth = truth_data["TABS"][:, truth_ind]

        data_simul_up_start = simul_data.group["profiles"]["updraft_temperature"][:, 1]
        data_simul_env_start = simul_data.group["profiles"]["env_temperature"][:, 1]
        data_truth_start = truth_data["TABS"][:, 1]



        p_t = plot(
            data_simul_up,
            z_simul,
            label = "Updraft",
            legend = :topright,
            ylim = (ymin, ymaxs[flight_number]),
            color = :red,
            dpi = dpi,
            xlabel = "T (K)",
            ylabel = "Height (m)",
            title = "Updraft and Environmental Temperature (K)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot

        plot!(
            simul_data.group["timeseries"]["Tsurface"][[simul_ind]],
            [0],
            label = "Surface Temp",
            seriestype = :scatter,
            markersize = 5,
            markercolor = :red,
        )
        plot!(
            truth_data["SST"][[truth_ind]],
            [0],
            label = "SST",
            seriestype = :scatter,
            markersize = 2,
            markercolor = :blue,
        )
        plot!(data_truth, z_truth, label = "Truth mean", color = :green) # add simulation data to same plot

        plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :pink, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_env_start, z_simul, label = "Environment start", color = :cyan, linestyle = :dash) # add simulation data to same plot
        plot!(data_truth_start, z_truth, label = "Truth start", color = :lime, linestyle = :dash, xlim = (250, 280)) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "updraft_env_temp_start_now.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # qt start vs now
        data_simul = simul_data.group["profiles"]["qt_mean"][:, simul_ind]
        data_simul_plus_precip =
            simul_data.group["profiles"]["qt_mean"][:, simul_ind] +
            simul_data.group["profiles"]["qs_mean"][:, simul_ind] +
            simul_data.group["profiles"]["qr_mean"][:, simul_ind]


        data_simul_up = simul_data.group["profiles"]["updraft_qt"][:, simul_ind]
        data_simul_env = simul_data.group["profiles"]["env_qt"][:, simul_ind]
        data_truth = truth_data["QT"][:, truth_ind] ./ 1000 # + 
        # truth_data["QS"][:, truth_ind] ./ 1000 +
        # truth_data["QR"][:, truth_ind] ./ 1000

        data_simul_start = simul_data.group["profiles"]["qt_mean"][:, 1]
        data_simul_up_start = simul_data.group["profiles"]["updraft_qt"][:, 1]
        data_simul_env_start = simul_data.group["profiles"]["env_qt"][:, 1]
        data_truth_start = truth_data["QT"][:, 1] ./ 1000


        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :topright,
            ylim = (ymin, ymaxs[flight_number]),
            dpi = dpi,
            marker = :circle,
            markersize = 1.5,
            xlabel = "qt (kg/kg)",
            ylabel = "Height (m)",
            title = "QT (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 1.5) # add simulation data to same plot
        plot!(data_simul_up, z_simul, label = "Updraft", color = :red) # add simulation data to same plot
        plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot
        plot!(data_simul_plus_precip, z_simul, label = "simul+Precip", color = :black) # add simulation data to same plot

        plot!(data_truth_start, z_truth, label = "Truth start", color = :green, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_start, z_simul, label = "Simulation start", color = :pink, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :red, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_env_start, z_simul, label = "Environment start", color = :blue, linestyle = :dash) # add simulation data to same plot


        data_simul_minus_cloud =
            simul_data.group["profiles"]["qt_mean"][:, simul_ind] -
            simul_data.group["profiles"]["qi_mean"][:, simul_ind] #- simul_data.group["profiles"]["qi_mean"][:, simul_ind]
        plot!(data_simul_minus_cloud, z_simul, label = "simul - cloud", color = :brown, linestyle = :dash) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "qt_start_now.png")) # save to file

        # --------------------------------------------------------------------------------------------- #

        # qt start vs now vs pressure
        # data_simul = simul_data.group["profiles"]["qt_mean"][:,simul_ind] 
        # data_simul_plus_precip = simul_data.group["profiles"]["qt_mean"][:,simul_ind]  + simul_data.group["profiles"]["qs_mean"][:,simul_ind]  + simul_data.group["profiles"]["qr_mean"][:,simul_ind]

        # data_simul_up = simul_data.group["profiles"]["updraft_qt"][:,simul_ind]
        # data_simul_env = simul_data.group["profiles"]["env_qt"][:,simul_ind]
        data_truth = truth_data["QT"][:, truth_ind] ./ 1000

        # data_simul_start = simul_data.group["profiles"]["qt_mean"][:,1]
        # data_simul_up_start = simul_data.group["profiles"]["updraft_qt"][:,1]
        # data_simul_env_start = simul_data.group["profiles"]["env_qt"][:,1]
        data_truth_start = truth_data["QT"][:, 1] ./ 1000

        p_truth_start = truth_data["PRES"][:, 1]

        p_t = plot(
            data_truth_start,
            p_truth_start,
            label = "Truth",
            legend = :bottomleft,
            dpi = dpi,
            marker = :circle,
            markersize = 1.5,
            yflip = true,
            xlabel = "qt (kg/kg)",
            ylabel = "Height (m)",
            title = "QT (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_truth, p_truth_start, label = "Truth now", marker = :circle, markersize = 1.5)



        savefig(joinpath(outpath, "Figures", "qt_start_now_vs_p.png")) # save to file

        # --------------------------------------------------------------------------------------------- #
        # ql start vs now
        data_simul = simul_data.group["profiles"]["ql_mean"][:, simul_ind]
        data_simul_plus_precip = simul_data.group["profiles"]["ql_mean"][:, simul_ind]
        # simul_data.group["profies"]["qr_mean"][:, simul_ind]

        data_simul_up = simul_data.group["profiles"]["updraft_ql"][:, simul_ind]
        data_simul_env = simul_data.group["profiles"]["env_ql"][:, simul_ind]
        data_truth = truth_data["QCL"][:, truth_ind] ./ 1000
        # truth_data["QS"][:, truth_ind] ./ 1000 +

        data_simul_start = simul_data.group["profiles"]["ql_mean"][:, 1]
        data_simul_up_start = simul_data.group["profiles"]["updraft_ql"][:, 1]
        data_simul_env_start = simul_data.group["profiles"]["env_ql"][:, 1]
        data_truth_start = truth_data["QCL"][:, 1] ./ 1000


        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            ylim = (ymin, ymaxs[flight_number]),
            dpi = dpi,
            marker = :circle,
            markersize = 1.5,
            xlabel = "qt (kg/kg)",
            ylabel = "Height (m)",
            title = "QL (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 1.5) # add simulation data to same plot
        plot!(data_simul_up, z_simul, label = "Updraft", color = :red) # add simulation data to same plot
        plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot
        plot!(data_simul_plus_precip, z_simul, label = "simul+Precip", color = :black) # add simulation data to same plot

        plot!(data_truth_start, z_truth, label = "Truth start", color = :green, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_start, z_simul, label = "Simulation start", color = :pink, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :red, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_env_start, z_simul, label = "Environment start", color = :blue, linestyle = :dash) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "ql_start_now.png")) # save to file


        # --------------------------------------------------------------------------------------------- #
        # liquid potential Temp start vs now
        data_simul = simul_data.group["profiles"]["thetal_mean"][:, simul_ind]

        data_simul_up = simul_data.group["profiles"]["updraft_thetal"][:, simul_ind]
        data_simul_env = simul_data.group["profiles"]["env_thetal"][:, simul_ind]
        data_truth = truth_data["THETAL"][:, truth_ind]

        data_simul_start = simul_data.group["profiles"]["thetal_mean"][:, 1]
        data_simul_up_start = simul_data.group["profiles"]["updraft_thetal"][:, 1]
        data_simul_env_start = simul_data.group["profiles"]["env_thetal"][:, 1]
        data_truth_start = truth_data["THETAL"][:, 1]


        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            ylim = (ymin, ymaxs[flight_number]),
            dpi = dpi,
            marker = :circle,
            markersize = 1.5,
            xlabel = "T (K)",
            ylabel = "Height (m)",
            title = "THETAL (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 1.5) # add simulation data to same plot
        plot!(data_simul_up, z_simul, label = "Updraft", color = :red) # add simulation data to same plot
        plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot

        plot!(data_truth_start, z_truth, label = "Truth start", color = :green, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_start, z_simul, label = "Simulation start", color = :pink, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :red, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_env_start, z_simul, label = "Environment start", color = :blue, linestyle = :dash) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "thetal_start_now.png")) # save to file

        # --------------------------------------------------------------------------------------------- #
        # temp start vs now 
        data_simul = simul_data.group["profiles"]["temperature_mean"][:, simul_ind]

        data_simul_up = simul_data.group["profiles"]["updraft_temperature"][:, simul_ind]
        data_simul_env = simul_data.group["profiles"]["env_temperature"][:, simul_ind]
        data_truth = truth_data["TABS"][:, truth_ind]

        data_simul_start = simul_data.group["profiles"]["temperature_mean"][:, 1]
        data_simul_up_start = simul_data.group["profiles"]["updraft_temperature"][:, 1]
        data_simul_env_start = simul_data.group["profiles"]["env_temperature"][:, 1]
        data_truth_start = truth_data["TABS"][:, 1]

        p_t = plot(
            data_truth,
            z_truth,
            label = "Truth",
            legend = :topright,
            ylim = (ymin, ymaxs[flight_number]),
            dpi = dpi,
            marker = :circle,
            markersize = 1.5,
            xlabel = "T (k)",
            ylabel = "Height (m)",
            title = "Temperature (K)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul, z_simul, label = "Simulation", marker = :circle, markersize = 1.5) # add simulation data to same plot
        plot!(data_simul_up, z_simul, label = "Updraft", color = :red) # add simulation data to same plot
        plot!(data_simul_env, z_simul, label = "Environment", color = :blue) # add simulation data to same plot

        plot!(data_truth_start, z_truth, label = "Truth start", color = :green, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_start, z_simul, label = "Simulation start", color = :pink, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_up_start, z_simul, label = "Updraft start", color = :red, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_env_start, z_simul, label = "Environment start", color = :blue, linestyle = :dash) # add simulation data to same plot

        savefig(joinpath(outpath, "Figures", "temp_start_now.png")) # save to file


        # --------------------------------------------------------------------------------------------- #

        # updraft area
        data_simul = simul_data.group["profiles"]["updraft_area"][:, simul_ind]

        plot(
            data_simul,
            z_simul,
            label = "Updraft Area",
            legend = :topright,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "Area",
            ylabel = "Height (m)",
            title = "Updraft Area",
        )
        savefig(joinpath(outpath, "Figures", "updraft_area.png")) # save to file

        # --------------------------------------------------------------------------------------------- #


        # rain snow (as a qt check)
        data_simul_rain = simul_data.group["profiles"]["qr_mean"][:, simul_ind]
        data_simul_snow = simul_data.group["profiles"]["qs_mean"][:, simul_ind]
        data_truth_rain = truth_data["QR"][:, truth_ind] ./ 1000
        data_truth_snow = truth_data["QS"][:, truth_ind] ./ 1000

        p_t = plot(
            data_truth_rain,
            z_truth,
            label = "Truth Rain",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qr (kg/kg)",
            ylabel = "Height (m)",
            title = "QR (kg/kg)",
            color = :green,
            linestyle = :dash,
            minorgrid = minorgrid,
        )
        plot!(data_simul_rain, z_simul, label = "Simulation Rain", color = :lime)
        plot!(data_truth_snow, z_truth, label = "Truth Snow", color = :blue, linestyle = :dash) # add simulation data to same plot
        plot!(data_simul_snow, z_simul, label = "Simulation Snow", color = :cyan) # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "qr_qs.png")) # save to file


        # reference states p

        p_truth = truth_data["PRES"][:, truth_ind]
        p_truth_ref = truth_data["p"][:]

        p_ref_simul = simul_data.group["reference"]["p_c"] / 100.0
        p_ref_simul_f = simul_data.group["reference"]["p_f"] / 100.0


        p = plot(
            p_truth,
            z_truth,
            label = "Truth",
            legend = :topright,
            dpi = dpi * 4,
            marker = :circle,
            markersize = 0.5,
            markerstrokewidth = 0.05,
            xlabel = "Pressure (Pa)",
            ylabel = "Height (m)",
            title = "Pressure (Pa)",
            minorgrid = minorgrid,
        )
        plot!(
            p_truth_ref,
            z_truth,
            label = "Truth ref",
            dpi = dpi,
            marker = :circle,
            markersize = 0.25,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :blue,
        ) # add simulation data to same plot
        plot!(
            p_ref_simul,
            z_simul,
            label = "Simulation ref",
            dpi = dpi,
            marker = :circle,
            markersize = 0.25,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :red,
        ) # add simulation data to same plot
        plot!(
            p_ref_simul_f,
            z_simul_f,
            label = "Simulation ref f",
            dpi = dpi,
            marker = :circle,
            markersize = 0.25,
            markerstrokewidth = 0.05,
            linewidth = 0.1,
            color = :green,
        ) # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "p.png")) # save to file

        # --------------------------------------------------------------------------------------------- #
        # ql + qr
        data_simul_lr = simul_data.group["profiles"]["ql_mean"][:, simul_ind] + simul_data.group["profiles"]["qr_mean"][:, simul_ind]
        data_truth_lr = truth_data["QCL"][:, truth_ind] ./ 1000 + truth_data["QR"][:, truth_ind] ./ 1000

        p_t = plot(
            data_truth_lr,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "ql + qr (kg/kg)",
            ylabel = "Height (m)",
            title = "QL + QR (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul_lr, z_simul, label = "Simulation") # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "ql+qr.png")) # save to file


        # --------------------------------------------------------------------------------------------- #
        # qi + qs
        data_simul_is = simul_data.group["profiles"]["qi_mean"][:, simul_ind] + simul_data.group["profiles"]["qs_mean"][:, simul_ind]
        data_truth_is = truth_data["QCI"][:, truth_ind] ./ 1000 + truth_data["QS"][:, truth_ind] ./ 1000

        p_t = plot(
            data_truth_is,
            z_truth,
            label = "Truth",
            legend = :bottomright,
            dpi = dpi,
            ylim = (ymin, ymaxs[flight_number]),
            xlabel = "qi + qs (kg/kg)",
            ylabel = "Height (m)",
            title = "QI + QS (kg/kg)",
            minorgrid = minorgrid,
        )
        p_s = plot!(data_simul_is, z_simul, label = "Simulation") # add simulation data to same plot
        savefig(joinpath(outpath, "Figures", "qi+qs.png")) # save to file

    end
end
