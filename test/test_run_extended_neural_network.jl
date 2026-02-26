## Standalone run script for SOCRATES using :neural_network_extended
## Keeps test_run.jl untouched.
using Pkg
Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
using JLD2
using Random
using LinearAlgebra

using Pkg
using Revise

Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))

Pkg.activate(expanduser("~/Research_Schneider/CliMA/CalibrateEDMF.jl/"))
using CalibrateEDMF

Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests"))
using Revise
using TurbulenceConvection

LinearAlgebra.BLAS.set_num_threads(1)

FT = Float64
forcing_type = :obs_data
flight_number = 9

tc = pkgdir(TurbulenceConvection)
CEDMF_dir = pkgdir(CalibrateEDMF)

Revise.includet(joinpath(tc, "driver", "main.jl"))
Revise.includet(joinpath(tc, "driver", "generate_namelist.jl"))

forcing_str = forcing_type === :obs_data ? "Obs" : "ERA5"

# Load base neural network namelist
nonequilibrium_moisture_scheme = :neural_network_extended
dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0"
method = "best_particle_final"
case_name = "SOCRATES_RF" * string(flight_number, pad=2) * "_" * lowercase(forcing_str) * "_data"
namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments", "SOCRATES_" * string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")

namelist = JSON.parsefile(namelist_path, dicttype=Dict, inttype=Int64, null=FT(NaN))
namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"])
default_namelist = NameList.default_namelist(case_name)
NameList.convert_namelist_types_to_default!(namelist, default_namelist)

# Override with extended neural network
namelist["user_args"]["nonequilibrium_moisture_scheme"] = :neural_network_extended


conservative_interp_kwargs_in = Dict{String,Union{Bool,String,FT}}( # Dict for easier namelist parsing later since that's all a dict.
    "conservative_interp" => true, # leave off by default (testing leaving off to see if our extrema fix works better w/o it. pointwise is ok if you've a point in each layer...)
    # "conservative_interp" => false, # leave off by default
    "preserve_monotonicity" => false,
    "enforce_positivity" => true, # gets auto-turned of as needed in SSCF
    "nnls_alg" => "pivot", # is faster, use as default
    "nnls_tol" => FT(1e-12),
    # "enforce_conservation" => false,
    # "enforce_conservation" => true, # testing if this is a layer killer ( update, it is )
    # "integrate_method" => "invert",
    "integrate_method" => "integrate", # this kills layers at low res unless you set reweight_processes_for_grid = true...
)

conservative_interp_kwargs_out = Dict{String,Union{Bool,String,FT}}( # Dict for easier namelist parsing later since that's all a dict.
    "conservative_interp" => true, # leave off by default (testing leaving off to see if our extrema fix works better w/o it. pointwise is ok if you've a point in each layer...)
    # "conservative_interp" => false, # leave off by default
    "preserve_monotonicity" => false, # I just don't think we need this, een for downsampling, at k=1...
    "enforce_positivity" => false, # this will get autoset as needed. and while CEDMF has fallbacks, we can't apply them so trivially everywhere...
    "nnls_alg" => "pivot", # is faster, use as default
    "nnls_tol" => FT(1e-12),
    # "enforce_conservation" => true, # let's just eschew this...
    # "integrate_method" => "invert",
    "integrate_method" => "invert",
)

namelist["grid"]["conservative_interp_kwargs"] = Dict(
    "in" => conservative_interp_kwargs_in,
    "out" => conservative_interp_kwargs_out
)

dz_min = namelist["grid"]["dz_min"]
if dz_min <= FT(100)
    namelist["grid"]["conservative_interp_kwargs"]["in"]["conservative_interp"] = false
    namelist["grid"]["conservative_interp_kwargs"]["out"]["conservative_interp"] = false
end


apply_edits = true
if apply_edits
    #  namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_model_type"] = "convective_tke"
     namelist["turbulence"]["EDMF_PrognosticTKE"]["convective_tke_model_type"] = "no_convective_tke"
end

dt = namelist["time_stepping"]["dt_min"]
namelist["time_stepping"]["dt_max"] = 2 * dt
namelist["time_stepping"]["algorithm"] = "Euler"
namelist["time_stepping"]["dt_limit_tendencies_factor"] = 2.0
namelist["time_stepping"]["use_tendency_timestep_limiter"] = false
namelist["time_stepping"]["N_dt_max_edmf_violate_dt_min"] = 1e3
namelist["time_stepping"]["allow_cfl_dt_max_violate_dt_min"] = false
namelist["time_stepping"]["use_fallback_during_spinup"] = true
namelist["time_stepping"]["allow_spinup_adapt_dt"] = true
namelist["time_stepping"]["spinup_dt_factor"] = 1.00
namelist["time_stepping"]["spinup_half_t_max"] = 3600.0 * 1.0

namelist["grid"]["dz_min"] = 5.0 * (2 * namelist["time_stepping"]["dt_min"] / 0.5)

calibrate_and_run_dir = joinpath(
    expanduser("~/Research_Schneider/CliMA/CalibrateEDMF.jl"),
    "experiments",
    "SOCRATES",
    "subexperiments",
    "SOCRATES_neural_network",
    "Calibrate_and_Run",
)

nn_params_file = joinpath(calibrate_and_run_dir, "pretrained_NN_extended.jld2")
nn_simple_chain_file = joinpath(calibrate_and_run_dir, "pretrained_NN_simple_chain_extended.jld2")

isfile(nn_params_file) || error("Missing extended NN parameter file: $nn_params_file")
isfile(nn_simple_chain_file) || error("Missing extended NN simple-chain file: $nn_simple_chain_file")
Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
using SimpleChains
Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests"))
nn_params, nn_x0 = JLD2.load(nn_params_file, "params", "x_0_characteristic")

# @warn "nn_params = $nn_params, nn_x0 = $nn_x0"

length(nn_x0) == 9 || error("Expected 9 extended input normalization values, got $(length(nn_x0)).")

namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network"] = nn_params
namelist["relaxation_timescale_params"]["model_x_0_characteristic"] = FT.(nn_x0)
namelist["relaxation_timescale_params"]["model_re_location"] = nn_simple_chain_file

# Optional quick-run defaults; tune as needed.
namelist["stats_io"]["frequency"] = 60.0
# namelist["stats_io"]["frequency"] = 600.0
# namelist["stats_io"]["frequency"] = namelist["time_stepping"]["dt_min"] # 1 timestep

namelist["stats_io"]["calibrate_io"] = false

# Put outputs in test/debug with a short uuid.
Random.seed!()
namelist["meta"]["uuid"] = randstring(4)
namelist["output"]["output_root"] = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/test/debug")

@info "Running case" case_name scheme = nonequilibrium_moisture_scheme
@info "Using NN files" nn_params_file nn_simple_chain_file

main1d(namelist)
