include(joinpath(@__DIR__, "common.jl"))
import Profile

# case_name = "Bomex"
# namelist = NameList.default_namelist(case_name)
# namelist["time_stepping"]["t_max"] = 15 * namelist["time_stepping"]["dt_max"]
# namelist["stats_io"]["frequency"] = 0 # io at every step




# nonequilibrium_moisture_scheme = :geometric_liq__exponential_T_scaling_and_geometric_ice
nonequilibrium_moisture_scheme = :geometric_liq__exponential_T_scaling_ice
# nonequilibrium_moisture_scheme = :exponential_T_scaling_ice
dt_string = "adapt_dt__dt_min_5.0__dt_max_10.0"
method = "best_particle_final"
flight_number = 9
forcing_str = "Obs"
CEDMF_dir = expanduser("~/Research_Schneider/CliMA/CalibrateEDMF.jl")
FT = Float64
case_name = "SOCRATES_RF" * string(flight_number, pad = 2) * "_" * lowercase(forcing_str) * "_data" # can't recall why it's lower here lol
namelist_path = joinpath(CEDMF_dir, "experiments", "SOCRATES_postprocess_runs_storage", "subexperiments","SOCRATES_"*string(nonequilibrium_moisture_scheme), "Calibrate_and_Run", "tau_autoconv_noneq", dt_string, "iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean", "postprocessing", "output", "Atlas_LES", "RFAll_obs", method, "data", "Output.$case_name.1_1", "namelist_SOCRATES.in")
namelist = JSON.parsefile(namelist_path, dicttype = Dict, inttype = Int64, null = FT(NaN))
namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"]) # this gets messed up for some reason...
default_namelist = default_namelist = NameList.default_namelist(case_name)
NameList.convert_namelist_types_to_default!(namelist, default_namelist) # coerce remaining type
print(namelist_path)
global debug = true # allow to run multiple simulaitons at once if we're just debugging and wan't to go faster, puts them in random subdirectories in target location with random uuid suffix\
# ========================================================================================================================= #
# ========================================================================================================================= #
# ========================================================================================================================= #

if debug
    Random.seed!() # switch back to random entropy
    namelist["meta"]["uuid"] = Random.randstring(4) # allow to run multiple simulaitons at once if we're just debugging and wan't to go faster
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
TC.full_print(namelist); println("\n\n")


namelist["meta"]["uuid"] = "01_flame_io"
sim = Simulation1d(namelist)
initialize(sim)
open_files(sim) # force compilation
close_files(sim) # force compilation
# open_files(sim)
# (prob, alg, kwargs) = solve_args(sim)
# integrator = ODE.init(prob, alg; kwargs...)

# ODE.step!(integrator) # force compilation
# # callbacks not called after first `step!`
# # call, so need to call `step!` again:
# ODE.step!(integrator) # force compilation

Profile.clear_malloc_data()
Profile.init(10_000_000, 0.001)  # n = buffer size, delay = sampling interval in seconds
prof = Profile.@profile begin
    # for _ in 1:10
    #     ODE.step!(integrator)
    # end

    main1d(namelist)
end
close_files(sim)
∑tendencies!
import PProf
PProf.pprof()
# http://localhost:57599/ui/flamegraph?tf
# import ProfileView
# ProfileView.view()
