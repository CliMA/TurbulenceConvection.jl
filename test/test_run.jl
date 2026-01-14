## Run test on cases and plot relative to LES

using JLD2
using Random
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1) # set to 1 thread for reproducibility and to stop slow bug

function edit_nt_key(nt::NamedTuple, key::Symbol, value)
    return merge(Base.structdiff(nt, NamedTuple{(key,)}(nt)), NamedTuple{(key,)}((value,)))
end

debug = false

FT = Float64
forcing_type = :obs_data

# flight_numbers = 1
flight_numbers = 9
# flight_numbers = 10
# flight_numbers = [1, 9, 10, 12, 13]

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

    using Revise
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
    using Revise
    using JLD2
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests"))
    using Revise
    using TurbulenceConvection
    using NCDatasets
    tc = pkgdir(TurbulenceConvection)
    Revise.includet(joinpath(tc, "driver", "main.jl"))
    Revise.includet(joinpath(tc, "driver", "generate_namelist.jl"))
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
reload_CalibrateEDMF = false
if reload_CalibrateEDMF | !isdefined(Main, :CalibrateEDMF)
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/CalibrateEDMF.jl/"))
    using CalibrateEDMF
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
    CEDMF_dir = pkgdir(CalibrateEDMF) # any submodule should work the same
end
# ------------------ #

namelist_stash = []
use_default_namelist = false # if true, use the default namelist, otherwise use the postprocessing namelist
for flight_number in flight_numbers
    setup_environment = true # leave true
    if setup_environment
        case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * forcing_str * "_data"

        # ========================================================================================================================= #
        # -- load from postprocessing namelist -- #

        if use_default_namelist

            if !@isdefined(NameList) # For if we didn't use reload_enviornment (I'm not sure tbf if this also checks for modules lol, I guess we'll find out)
                @info("NameList not defined by module, defining now as simple Dict()") # This allows plotting even if you haven't gone through the environment loading by allowing you to get the simname/uuid etc in namelist 
                namelist = Dict{String, Union{Dict, NamedTuple}}(
                    "time_stepping" => Dict(),
                    "stats_io" => Dict(),
                    "microphysics" => Dict(),
                    "user_args" => Dict(),
                    "user_params" => Dict(),
                    "relaxation_timescale_params" => Dict(),
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
                namelist["relaxation_timescale_params"] = Dict()
            end

            if namelist["user_args"] isa NamedTuple
                @warn ("Converting user_args from NamedTuple to Dict")
                namelist["user_args"] = Dict(pairs(namelist["user_args"])) # convert to dict
            else
                @info ("user_args is already a Dict")
            end

            nonequilibrium_moisture_scheme = :Base

        else
            # ================================================================================================= #
            # -- run from postprocessing namelist -- #

            # nonequilibrium_moisture_scheme = :Base
            # dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0"
            # method = "best_particle_final"
            # flight_number = 12
            # case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
            # namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "pow_icenuc_autoconv_eq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
            # namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
            # namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
            # default_namelist = default_namelist = NameList.default_namelist(case_name)
            # NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type

            # nonequilibrium_moisture_scheme = :geometric_liq__exponential_T_scaling_and_geometric_ice
            # nonequilibrium_moisture_scheme = :geometric_liq__exponential_T_scaling_ice
            nonequilibrium_moisture_scheme = :exponential_T_scaling_ice
            dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0"
            method = "best_particle_final"
            flight_number = 9
            case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
            namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
            namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
            namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
            default_namelist = default_namelist = NameList.default_namelist(case_name)
            NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type
            print(namelist_path)

            # nonequilibrium_moisture_scheme = :neural_network
            # # dt_string = "adapt_dt__dt_min_2.0__dt_max_4.0"
            # dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0" # this is the one we want to run
            # method = "best_particle_final"
            # case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
            # namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
            # namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
            # namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
            # default_namelist = default_namelist = NameList.default_namelist(case_name)
            # NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type


            # nonequilibrium_moisture_scheme = :linear_combination
            # # dt_string = "adapt_dt__dt_min_2.0__dt_max_4.0"
            # dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0" # this is the one we want to run
            # method = "best_particle_final"
            # case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
            # namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
            # namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
            # namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
            # default_namelist = default_namelist = NameList.default_namelist(case_name)
            # NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type

            # nonequilibrium_moisture_scheme = :neural_network # random init [[ sets below]]
            # # dt_string = "adapt_dt__dt_min_2.0__dt_max_4.0"
            # dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0" # this is the one we want to run
            # method = "best_particle_final"
            # case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
            # namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
            # namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
            # namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
            # default_namelist = default_namelist = NameList.default_namelist(case_name)
            # NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type
            # nonequilibrium_moisture_scheme = :neural_network_random_init
            # namelist["user_args"]["nonequilibrium_moisture_scheme"] = nonequilibrium_moisture_scheme    


            # nonequilibrium_moisture_scheme = :neural_network_pca_noise
            # dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0" # this is the one we want to run
            # method = "best_particle_final"
            # case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
            # namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
            # namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
            # namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
            # default_namelist = default_namelist = NameList.default_namelist(case_name)
            # NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type
            # nonequilibrium_moisture_scheme = :neural_network_random_init
            # namelist["user_args"]["nonequilibrium_moisture_scheme"] = nonequilibrium_moisture_scheme    
            # # 
            # include(joinpath(CEDMF_dir, "experiments", "SOCRATES", "subexperiments", "SOCRATES_neural_network", "Calibrate_and_Run", "NN_pca.jl")) # this will load the NN PCA functions and add them to the valid experiment setups
            # experiment_dir = joinpath(CEDMF_dir, "experiments", "SOCRATES", "subexperiments", "SOCRATES_neural_network_pca_noise") # this is the path to the experiment directory, we will use it in the namelist
            # neural_network_dir = joinpath(dirname(experiment_dir), "SOCRATES_neural_network") # experiment_dir has no trailing slash so this should work
            # pretrained_NN_checkpoints_path = joinpath(neural_network_dir, "Calibrate_and_Run", "pretrained_NN_checkpoints.jld2")
            # _, model_params_matrix, _ = JLD2.load(pretrained_NN_checkpoints_path, "model_checkpoint_states", "model_checkpoint_params", "model_checkpoint_losses") # load the previous states
            # #
            # nn_path = joinpath(neural_network_dir, "Calibrate_and_Run", "pretrained_NN.jld2") # this is the path to the pretrained NN file, we will use it in the namelist
            # nn_pretrained_params, nn_pretrained_repr, nn_pretrained_x_0_characteristic = JLD2.load(nn_path, "params", "re", "x_0_characteristic")
            # #
            # nn_path = joinpath(neural_network_dir, "Calibrate_and_Run", "pretrained_NN_simple_chain.jld2") # needed to be this for simple_chain loading later.
            # #
            # i_checkpoint = 500 # let's trial starting from scratch
            # model_params = model_params_matrix[:, i_checkpoint+1]
            # #
            # n_pca = 12 # number of PCA components to use, this is the number of principal components we will use as noise
            # (; input_train, truth_train) = get_NN_training_data()
            # by_output = true; reference_to_truth=false
            # # pca_components, pca_explained_variance, pca_mean_vec = get_NN_pca_by_output(model_params, nn_pretrained_repr, input_train, truth_train; number_components=n_pca, filter_threshold=CalibrateEDMF.FTNN(0), scale_gradients=false, reference_to_truth=reference_to_truth)
            # pca_components, pca_explained_variance, pca_mean_vec = get_output_specific_directions_by_output(model_params, nn_pretrained_repr, input_train, truth_train; number_components=n_pca, filter_threshold=CalibrateEDMF.FTNN(0), scale_gradients=false, reference_to_truth=reference_to_truth)
            # #
            # namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_number_components"] = n_pca # set the number of PCA components to use
            # namelist["relaxation_timescale_params"]["reference_to_truth"] = reference_to_truth # set the explained variance of the PCA components
            # namelist["relaxation_timescale_params"]["pca_by_output"] = by_output # set whether the PCA is by output or not
            # namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network"] = model_params
            # namelist["relaxation_timescale_params"]["specific_toggle"] = true # maximize toggling by direction by solving the inverse problem.
            # namelist["relaxation_timescale_params"]["fine_tuning_only"] = true
            # namelist["relaxation_timescale_params"]["model_x_0_characteristic"] = nn_pretrained_x_0_characteristic
            # namelist["relaxation_timescale_params"]["model_re_location"] = nn_path # set the path to the pretrained NN file
            # #
            # if !by_output
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_weights"] = randn(n_pca+1) # initialize the PCA weights to random values
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_mean"] = pca_mean_vec # set the PCA mean to the mean of the training data
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_components"] = pca_components # set the PCA components to the components we just calculated 
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_explained_variance"] = pca_explained_variance # set the PCA explained variance to the explained variance we just calculated
            # else
            #     n_outputs = 4
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_weights"] = randn(n_pca+n_outputs) # initialize the PCA weights to random values
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_mean"] = [FT.(pca_mean_vec[i_o][:]) for i_o in 1:n_outputs] # set the PCA mean to the mean of the training data
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_components"] = [FT.(pca_components[i_o]) for i_o in 1:n_outputs] # set the PCA components to the components we just calculated
            #     namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_explained_variance"] = [FT.(pca_explained_variance[i_o]) for i_o in 1:n_outputs] # set the PCA explained variance to the explained variance we just calculated
            # end

            # ================================================================================================= #
            # -- take parameters from a Diagnostics.nc (here you still need to set manually) -- #

            # nonequilibrium_moisture_scheme = :geometric_liq__powerlaw_T_scaling_ice
            # path_to_my_SOTA = "/resnick/groups/esm/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/subexperiments/SOCRATES_geometric_liq__powerlaw_T_scaling_ice/Calibrate_and_Run/tau_autoconv_noneq/adapt_dt__dt_min_5.0__dt_max_10.0/iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean/calibrate/output/Atlas_LES/RFAll_obs/Diagnostics.nc"
            # my_SOTA = optimal_parameters(path_to_my_SOTA, method = "best_particle_final")
            # my_SOTA = Dict(zip(my_SOTA...)) # turn to dict
            # namelist["relaxation_timescale_params"]["powerlaw_T_scaling_ice_c_1"] = my_SOTA["powerlaw_T_scaling_ice_c_1"]
            # namelist["relaxation_timescale_params"]["powerlaw_T_scaling_ice_c_2"] = my_SOTA["powerlaw_T_scaling_ice_c_2"]

            # namelist["relaxation_timescale_params"]["geometric_liq_c_1"] = my_SOTA["geometric_liq_c_1"]
            # namelist["relaxation_timescale_params"]["geometric_liq_c_2"] = my_SOTA["geometric_liq_c_2"]
            # namelist["relaxation_timescale_params"]["geometric_liq_c_3"] = my_SOTA["geometric_liq_c_3"]
            # namelist["relaxation_timescale_params"]["geometric_liq_c_4"] = my_SOTA["geometric_liq_c_4"] 
            
            
            # # ================================================================================================= #
            # # --- run this jld2 --- #

            # param_results_path = "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_tCb99B_simulation_crashed/" # Output.SOCRATES_RF10_obs_data.jl_tCb99B/stats/Stats.SOCRATES_RF10_obs_data.nc"
            # param_results_path = "/resnick/scratch/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_djcmkh_simulation_crashed/" # Output.SOCRATES_RF12_obs_data.jl_djcmkh/stats/Stats.SOCRATES_RF12_obs_data.nc"
            # param_results_path = "/resnick/scratch/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_rEeUQl_simulation_crashed/" # Output.SOCRATES_RF10_obs_data.jl_rEeUQl/stats/Stats.SOCRATES_RF10_obs_data.nc"
            # param_results_path = "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_6vcECq_simulation_crashed/" #Output.SOCRATES_RF13_obs_data.jl_6vcECq/stats/Stats.SOCRATES_RF13_obs_data.nc"
            # param_results_path = "/resnick/scratch/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_xiIubQ_simulation_crashed" #/param_results.jld2"
            # param_results_path = "/resnick/scratch/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_Hty4CZ_simulation_crashed"
            # param_results_path = "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/julia_parallel/slurm/param_results/jl_SNkTvT_simulation_crashed/" #Output.SOCRATES_RF10_obs_data.jl_SNkTvT/stats/Stats.SOCRATES_RF10_obs_data.nc"

            # setup =  JLD2.load(param_results_path * "/param_results.jld2")["param_results"]
            # namelist = setup["namelist"]
            # nonequilibrium_moisture_scheme = namelist["user_args"]["nonequilibrium_moisture_scheme"]
            # dt = namelist["time_stepping"]["dt_min"]

            # change key nonequilibrium_moisture_sources_limiter_type in namelist["user_args"] (immutable named tuple so need to reconstruct)
            # namelist["user_args"]["nonequilibrium_moisture_sources_limiter_type"] = "morrison_milbrandt_2015_style"
            # namelist["user_args"]["nonequilibrium_moisture_sources_limiter_type"] = "standard_supersaturation"
            # namelist["user_args"]["fallback_nonequilibrium_moisture_sources_limiter_type"] = namelist["user_args"]["nonequilibrium_moisture_sources_limiter_type"]
            # # ============================================================================================================================ #

            apply_edits = false
        end 

        # [[ -- fixes -- ]]

        if apply_edits
            if isnothing(get(namelist["grid"], "conservative_interp_kwargs", nothing)) || isnothing(get(namelist["grid"]["conservative_interp_kwargs"], "in", nothing))
                namelist["grid"]["conservative_interp_kwargs"] = Dict("in" => Dict{Symbol, Union{Bool, String}}(), "out" => Dict{Symbol, Union{Bool, String}}())
            end
            namelist["grid"]["conservative_interp_kwargs"]["in"] = Dict{Symbol, Union{Bool, String}}(
                :conservative_interp => true, # leave off by default
                :preserve_monotonicity => false,
                :enforce_positivity => true,
                :nnls_alg => "pivot",
                :enforce_conservation => false,
                :integrate_method => "invert",
                # :integrate_method => :integrate,
            )
        end

        original_namelist = namelist
        namelist["meta"]["simname"] = original_namelist["meta"]["simname"]


        # ========================================================================================================================= #
        # -- Costa_SOTA parameters -- #
        if apply_edits
            # path_to_Costa_SOTA = "/groups/esm/cchristo/cedmf_results/james_v1_runs/results_Inversion_p22_e300_i15_mb_LES_2024-03-15_10-08_Vxo_longer_long_run/Diagnostics.nc"
            path_to_Costa_SOTA = joinpath(CEDMF_dir, "experiments", "SOCRATES", "Reference", "Costa_SOTA_Diagnostics.nc") #  "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/Reference/Costa_SOTA_Diagnostics.nc"
            Revise.includet(joinpath(CEDMF_dir, "tools", "DiagnosticsTools.jl")) # provides optimal_parameters()
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

            # args from costa....
            namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_type"] = "total_rate"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] = "w_height"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["detr_dim_scale"] = "mf_grad"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.0
            namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "None"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["ml_entrainment"] = "Linear"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["entr_pi_subset"] = [1, 2, 3, 4, 6]
            namelist["turbulence"]["EDMF_PrognosticTKE"]["pi_norm_consts"] = [100.0, 2.0, 1.0, 1.0, 1.0, 1.0]
            namelist["turbulence"]["EDMF_PrognosticTKE"]["linear_ent_biases"] = true

            # ========================================================================================================================= #
            # ========================================================================================================================= #
            # ========================================================================================================================= #

            # namelist["thermodynamics"]["moisture_model"] = "equilibrium"
            # namelist["thermodynamics"]["moisture_model"] = "nonequilibrium"
            # namelist["thermodynamics"]["sgs"] = "mean"


            user_args = get!(namelist, "user_args", Dict()) # make sure user_args exists, otherwise it will error out later

            merge!(user_args, Dict(
                "nonequilibrium_moisture_scheme" => nonequilibrium_moisture_scheme,
                # "nonequilibrium_moisture_sources_limiter_type" = :none, # can go past supersaturation but tendency limited timesteps prevents
                # "entr_detr_limiter_type" = :none,
                # "precipitation_tendency_limiter_type" = :none,
                # "default_tendency_limiter_type" = :none,
                #
                # nonequilibrium_moisture_sources_limiter_type = :none,
                # nonequilibrium_moisture_sources_limiter_type = :basic,
                # nonequilibrium_moisture_sources_limiter_type = :standard_supersaturation,
                # nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style,
                "nonequilibrium_moisture_sources_limiter_type"  => :morrison_milbrandt_2015_style_exponential_part_only,
                #
                # fallback_nonequilibrium_moisture_sources_limiter_type = :basic, # can go past supersaturation but tendency limited timesteps should prevent...
                # fallback_nonequilibrium_moisture_sources_limiter_type = :standard_supersaturation, # has limiter issues...
                # fallback_nonequilibrium_moisture_sources_limiter_type = :morrison_milbrandt_2015_style,
                "fallback_nonequilibrium_moisture_sources_limiter_type"  => :morrison_milbrandt_2015_style_exponential_part_only,
                "fallback_to_standard_supersaturation_limiter"  => true,
                # precipitation_tendency_limiter_type = :basic,
                # default_tendency_limiter_type = :basic.
                # fallback_precipitation_tendency_limiter_type = :basic,
                # fallback_default_tendency_limiter_type = :basic,
                #
                "precipitation_tendency_limiter_type"  => :truncated_basic,
                "fallback_precipitation_tendency_limiter_type"  => :truncated_basic,
                "default_tendency_limiter_type"  => :truncated_basic,
                "fallback_default_tendency_limiter_type"  => :truncated_basic,
                "truncated_basic_limiter_factor"  => FT(0.9),
                #
                # "entr_detr_limiter_type"  => :basic, # kills updrafts
                # "fallback_entr_detr_limiter_type"  => :basic,
                "entr_detr_limiter_type"  => :none,
                "fallback_entr_detr_limiter_type"  => :none,
                #
                "use_sedimentation"  => true, 
                "grid_mean_sedimentation"  => false,
                "sedimentation_integration_method"  => :upwinding, # could be :right_biased I suppose.
                "use_heterogeneous_ice_nucleation"  => true, # i think we do need this at correct values to make sure ice doesn't lag around forever nd gets the jumpstart it needs...
                # "sedimentation_ice_number_concentration" => nonequilibrium_moisture_scheme,
                # "sedimentation_liq_number_concentration" => nonequilibrium_moisture_scheme, # idk if this is good or bad lol...
                "liq_terminal_velocity_scheme"  => :Chen2022Vel,
                "ice_terminal_velocity_scheme"  => :Chen2022Vel,
                "rain_terminal_velocity_scheme"  => :Chen2022Vel,
                "snow_terminal_velocity_scheme"  => :Chen2022Vel,

                "adjust_liq_N" => true, # adjust liquid number concentration
                "adjust_ice_N" => true,
                )
            )

        end
            
        # ---------------------------------------------------------------- #

        # namelist["time_stepping"]["dt_min"] = 50.
        # namelist["time_stepping"]["dt_max"] = 100.
        # namelist["time_stepping"]["spinup_half_t_max"] = FT(3600) * 0

        # namelist["stats_io"]["frequency"] = .1
        # namelist["stats_io"]["frequency"] = 30.0
        # namelist["stats_io"]["frequency"] = namelist["time_stepping"]["dt_min"] # 1 timestep
        namelist["stats_io"]["frequency"] = 60.
        namelist["stats_io"]["calibrate_io"] = false

        # namelist["stats_io"]["frequency"] = 600.
        # namelist["stats_io"]["calibrate_io"] = true
        # namelist["thermodynamics"]["moisture_model"] = "equilibrium"

        # nonequilibrium_moisture_scheme = :neural_network_pca_noise
        # namelist["user_args"]["nonequilibrium_moisture_scheme"] = nonequilibrium_moisture_scheme # set the nonequilibrium moisture scheme to neural_network_pca_noise



        # namelist["user_args"]["adjust_ice_N"] = true # adjust ice number concentration

        # namelist["user_params"]["initial_profile_updraft_area"] = FT(0)

        # namelist["user_args"]["tendency_resolver_setup"] = :none # i think either dict or named tuple is fine, we dispatch depending... # doesn't seem to work anyway...
        # namelist["user_args"]["tendency_resolver_setup"] = :normal # i think either dict or named tuple is fine, we dispatch depending... # doesn't seem to work anyway...
        # namelist["user_args"]["tendency_resolver_setup"] = :full_tendencies # i think either dict or named tuple is fine, we dispatch depending... # doesn't seem to work anyway...

        # namelist["user_args"]["entr_detr_limiter_type"] = "none" # no entrainment/detrainment limiter
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["stalled_updraft_handler"] = "kill"
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["stalled_updraft_handler"] = "none" # kill stalled updrafts
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["stalled_updraft_handler"] = "detrain_downdrafts" # kill stalled updrafts
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["base_detrainment_rate_inv_s"] = FT(1/(60.0*60.0)) # 10 minutes
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["base_entrainment_rate_inv_s"] = FT(.0) # strong... hopefully doesn't trigger max area limiters
        # namelist["user_params"]["ice_sedimentation_Dmax"] = FT(Inf)
        # namelist["user_args"]["use_heterogeneous_ice_nucleation"] = true

        # # namelist["user_args"]["truncated_basic_limiter_factor"] =  FT(0.9) # make things take more than 1 timestep...
        # namelist["user_args"]["entr_detr_limiter_type"] = "none"

        # path_to_Costa_SOTA = "/groups/esm/cchristo/cedmf_results/james_v1_runs/results_Inversion_p22_e300_i15_mb_LES_2024-03-15_10-08_Vxo_longer_long_run/Diagnostics.nc"
        path_to_Costa_SOTA = joinpath(CEDMF_dir, "experiments", "SOCRATES", "Reference", "Costa_SOTA_Diagnostics.nc") #  "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/Reference/Costa_SOTA_Diagnostics.nc"
        Pkg.activate(CEDMF_dir)
        Revise.includet(joinpath(CEDMF_dir, "tools", "DiagnosticsTools.jl")) # provides optimal_parameters()
        Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl")) # go back to climajl environment
        Costa_SOTA = optimal_parameters(path_to_Costa_SOTA, method = "last_nn_particle_mean")
        Costa_SOTA = Dict(zip(Costa_SOTA...)) # turn to dict

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = Costa_SOTA["pressure_normalmode_adv_coeff"]
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = .5
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = Costa_SOTA["pressure_normalmode_drag_coeff"] / 5


        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = Costa_SOTA["tke_ed_coeff"] / 2
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = Costa_SOTA["tke_diss_coeff"] / 2
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = Costa_SOTA["static_stab_coeff"] / 2
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] = Costa_SOTA["tke_surf_scale"] * 2

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.33) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)

        # namelist["meta"]["simname"] = "SOCRATES_RF01_obs_data"
        # namelist["meta"]["flight_number"] = 01

        namelist["meta"]["simname"] = "SOCRATES_RF09_obs_data"
        namelist["meta"]["flight_number"] = 09

        # namelist["user_params"]["χm_liq"] = FT(1)
        # namelist["microphysics"]["χm_ice"] = FT(1) # inv(namelist["user_params"]["mean_r_factor_ice"])^3
        # namelist["user_params"]["massflux_N_i_boost_factor"] = FT(2.0)
        # namelist["user_params"]["sedimentation_N_i_boost_factor"] = FT(0.2)
        # namelist["user_params"]["apply_massflux_N_i_boost"] = true
        # namelist["user_params"]["apply_sedimentation_N_i_boost"] = false
        # namelist["user_args"]["use_ice_mult"] = false
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["base_detrainment_rate_inv_s"] = FT(1/(0.1*3600.0)) # 6 hours
        # namelist["relaxation_timescale_params"]["τ_sub_dep_scaling_factor"] = FT(500)
        # namelist["user_params"]["ice_dep_acnv_scaling_factor"] = FT(1.0)
        # namelist["user_params"]["ice_dep_acnv_scaling_factor_above"] = FT(1.0)
        namelist["thermodynamics"]["sgs"] = "mean"
        # namelist["thermodynamics"]["sgs"] = "mean_w_quadrature_adjusted_noneq_moisture_sources"
        # namelist["thermodynamics"]["sgs"] = "quadrature"
        # namelist["thermodynamics"]["quadrature_order"] = 12

        # namelist["microphysics"]["τ_cond_evap"] = FT(1.)
        # namelist["microphysics"]["χv_ice"] = FT(0.0001)
        # namelist["user_params"]["ice_sedimentation_scaling_factor"] = FT(0.3)


        # namelist["user_params"]["initial_profile_updraft_area"] = FT(0.0001)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.5) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area_bc"] = "Fixed"

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_partition_model"] = "standard"
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["apply_second_order_flux_correction"] = false

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["area_partition_model"] = "core_cloak"
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["cloak_area_factor"] = 0.01
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["cloak_dn_area_factor"] = 10.
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["cloak_mix_factor"] = 0.05
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["cloak_dn_area_factor"] = 1.0
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["confine_all_downdraft_to_cloak"] = false # true leads to stronger transport across z = 22, combined w/ second order correction seems to lead to better qt mixing?
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["apply_second_order_flux_correction"] = true # true w/ sgs crashed with LBF RBF logic [[ I can't remember why i thought this was necessary, it only really matters for area, no? Having it on leads to huge temp problems at the inversion, though we also do not get as much qt downward transport...]]
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["second_order_correction_limit_factor"] = Inf # Especially above 0.8, this can clash with quadrature... but it is a huge help in actually allowing for downwards qt advection w/o the need for prognostic cloaks.
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["apply_cloak_to_condensate_formation"] = false # true w/ sgs crashed with LBF RBF logic [[ I can't remember why i thought this was necessary, it only really matters for area, no? Having it on leads to huge temp problems at the inversion, though we also do not get as much qt downward transport...]]
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["fraction_of_area_above_max_area_allowed_in_cloak"] = 0.0
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["l_max"] = 1e4 # default was 1e6



        namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_buoyancy_coeff"] = FT(.25)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_advection_coeff"] = FT(0.1)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_dissipation_coeff"] = FT(0.05) # smaller than generation
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_self_dissipation_coeff"] = FT(0.002) # extra dissipation when tke is convectively generated [ this one is the problem ]
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_transport_tke_by_advection"] = true
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_transport_conserved_by_advection"] = true
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_transport_condensed_by_advection"] = true

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_convective_max_scaling_factor"] = FT(0.66)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_ed_scaling_factor"] = (1/(namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"])) / 2 # so that ed is same as default at max convective tke
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_ed_scaling_factor"] = FT(20)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_ed_scaling_factor"] = FT(1)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_conv_entr_detr_rate_inv_s"] = FT(0)

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_model_type"] = "convective_tke_production_and_graft_only"
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_model_type"] = "convective_tke"

        # namelist["user_params"]["S_ice_min_activation"] = FT(0.1) # default 0.01, higher means less ice nucleation


        delete!.(Ref(namelist["user_params"]), ["convective_tke_buoyancy_coeff", "convective_tke_advection_coeff", "convective_tke_dissipation_coeff", "convective_tke_self_dissipation_coeff", "entr_detr_rate_inv_s", "convective_tke_ed_scaling_factor", "tke_conv_entr_inv_s"])

        # namelist["user_params"]["apply_nudging_to_updraft"] = true # if this is true, the updraft will feel the nudging just as much as the environment, otherwise the env would feel all of it... I think we need this. The LES updraft for RF09 doesn't go that high. without this i think we risk overshoots when say a warm updraft impinges, the environment takes the nudging cooling in response, and the updraft continues to accelerate.
        # namelist["user_params"]["condensate_qt_SD"] = -.26

        # namelist["user_args"]["nonequilibrium_moisture_sources_limiter_type"] = "none"
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["stalled_updraft_handler"] = "kill"

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["base_detrainment_rate_inv_s"] = FT(1/(0.1*3600.0)) # 360 seconds


        # namelist["user_params"]["massflux_N_i_boost_max_ratio"] = 1.0
        # namelist["user_params"]["massflux_N_i_boost_progress_fraction"] = 1.0

        # I think through compute_en_tendencies!(), entr/detr mediate tke, which in turn mediates mixing length and tke contribution to diffusive fluxes.... I don't understand why z should matter for tke.
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] = "abs_massflux" # trying to use this intead of w_height... idk. we basically get no entrainment except at the sfc and then huge detraiment at cloud top. might explain our tiny areas. This is also bad because our diffusive flux is tied to entrainment it seems? idk...
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["detr_dim_scale"] = "mf_grad"
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["detr_dim_scale"] = "b_sqrt_tke"
        namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.1
        # namelist["user_params"]["initial_profile_updraft_area"] = FT(0.9*namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"]) * 0.0001
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.6) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_conv_entr_detr_rate_inv_s"] = FT(1/(0.0100*3600.0)) # 6 minutes
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_conv_entr_detr_rate_inv_s"] = FT(1e-2)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_conv_entr_detr_rate_inv_s"] = FT(10)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_conv_entr_detr_rate_inv_s"] = FT(0)

        # parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_ed_coeff")
        namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.0002
        # namelist["user_params"]["stable_updraft_area_reduction_factor"] = FT(5.0)
        # namelist["microphysics"]["microph_scaling_dep_sub"] = FT(2)
        # namelist["microphysics"]["microph_scaling_melt"] = FT(.1)
        # namelist["user_params"]["q_min"] = eps(FT)

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KQl"] = FT(1)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KQi"] = FT(1)

        # namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KTKEql"] = FT(.1)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KTKEqi"] = FT(.1)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KTKEqt"] = FT(.1)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KTKEh"]  = FT(.1)
        namelist["turbulence"]["EDMF_PrognosticTKE"]["c_KTKEqs"]  = FT(.5)
        # namelist["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] = 0.001

        # namelist["user_args"]["snow_terminal_velocity_scheme"] = "Blk1MVel"

        # namelist["microphysics"]["χv_sno"] = FT(0.5) # inv(namelist["user_params"]["mean_r_factor_ice"])^3 ("microphysics", "χv_sno"
        # namelist["user_params"]["zrough"] = FT(5e-5)
        # namelist["user_params"]["τ_acnv_sno_threshold"] = FT(200.0)

        if apply_edits
            dt = namelist["time_stepping"]["dt_min"]
            namelist["time_stepping"]["dt_min"] = dt
            namelist["time_stepping"]["dt_max"] = 2*dt

            namelist["time_stepping"]["algorithm"] = "Euler"; # jacobian too slow
            namelist["time_stepping"]["dt_limit_tendencies_factor"] = 2.0 # for safety (though now it'll probably go up to 20.0 sec)
            namelist["time_stepping"]["use_tendency_timestep_limiter"] = false  # testing
            namelist["time_stepping"]["N_dt_max_edmf_violate_dt_min"] = 1e3
            namelist["time_stepping"]["allow_cfl_dt_max_violate_dt_min"] = false   
            namelist["time_stepping"]["use_fallback_during_spinup"] = true
            namelist["time_stepping"]["allow_spinup_adapt_dt"] = true
            # namelist["time_stepping"]["limit_tendencies_using_dt_min_factor"] = 1.0 # margin of safety [[ depreacated ]]
            namelist["time_stepping"]["spinup_dt_factor"] = 1.00
            namelist["time_stepping"]["spinup_half_t_max"] = 3600.0 * 1.0 # 1 hour

            namelist["grid"]["dz_min"] = 5.0 * (2 * namelist["time_stepping"]["dt_min"] / 0.5)

            # namelist["user_params"]["print_taus"] = false
            # namelist["user_params"]["print_sources"] = false


            namelist["user_params"]["q_min"] = 0. #eps(FT)
            namelist["user_params"]["min_τ_liq"] = eps(FT)
            namelist["user_params"]["max_τ_ice"] = 1/eps(FT)

            namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = FT(0.9) # stability (maybe we need to use the limiter instead tho to not get flat cloud tops?)

            # ------------------------------------------------------------ #

            namelist["user_params"]["reweight_processes_for_grid"] = false # whatever is going on in reweight make things much worse... (but it's still happening)
            namelist["user_params"]["reweight_extrema_only"] = false
            # namelist["turbulence"]["EDMF_PrognosticTKE"]["stalled_updraft_handler"] = "kill"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["stalled_updraft_handler"] = "mix_to_grid_mean"
            namelist["turbulence"]["EDMF_PrognosticTKE"]["base_detrainment_rate_inv_s"] = FT(1/(6.0*3600.0)) # 6 hours
            # namelist["user_params"]["tendency_resolver_setup"] = "none" # i think either dict or named tuple is fine, we dispatch depending... # doesn't seem to work anyway...

            # namelist["user_args"] =  CalibrateEDMF.HelperFuncs.set_property(namelist["user_args"], :truncated_basic_limiter_factor,  FT(0.9)) # make things take more than 1 timestep...
            # namelist["user_args"] =  CalibrateEDMF.HelperFuncs.set_property(namelist["user_args"], :fallback_to_standard_supersaturation_limiter, true) # make things take more than 1 timestep...
            # namelist["user_args"] =  CalibrateEDMF.HelperFuncs.set_property(namelist["user_args"], :nonequilibrium_moisture_sources_limiter_type, :standard_supersaturation) # make things take more than 1 timestep...

            if namelist["user_params"]["reweight_processes_for_grid"] 
                @warn ("reweight_processes_for_grid is true, this is a test")
            end
        end

        global debug = true # allow to run multiple simulaitons at once if we're just debugging and wan't to go faster, puts them in random subdirectories in target location with random uuid suffix\
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
        case_name = simname # in case the jld2 changed it
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

        # ============================================================================================================================ #
        # -- print namelist -- #
        # search_locator__1
        TC.full_print(namelist); println("\n\n")
        # ============================================================================================================================ #

    end


 push!(namelist_stash, namelist)

# ========================================================================================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================================================================================== #

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


# ========================================================================================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================================================================================== #


    make_plots = false
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

        # toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
        # aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
        # param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
        # thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
        # param_set = create_parameter_set(namelist, toml_dict, FT)

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
