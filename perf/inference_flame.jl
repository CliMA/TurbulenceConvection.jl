include(joinpath(@__DIR__, "common.jl"))
import SnoopCompileCore

# case_name = "Bomex"
# println("Running $case_name...")


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
default_namelist = NameList.default_namelist(case_name)
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

namelist["meta"]["uuid"] = "01_invalidations"
# sim = Simulation1d(namelist)
# initialize(sim)
sim = init_sim(namelist; skip_io = true, single_timestep = true, prefix = "" )




# sim = init_sim(case_name)
sim.skip_io || open_files(sim) # #removeVarsHack
(prob, alg, kwargs) = solve_args(sim)

tinf = SnoopCompileCore.@snoop_inference begin
    sol = SciMLBase.solve(prob, alg; kwargs...)
    # integrator, ds_tc_filenames, return_code = main(namelist)
end;

sim.skip_io || close_files(sim) # #removeVarsHack

# import ProfileView
import SnoopCompile # need SnoopCompile to iterate over InferenceTimingNode's
import FlameGraphs
fg = FlameGraphs.flamegraph(tinf)

using Serialization
# thisdir = @__DIR__
thisdir = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/perf/")
folder = joinpath(thisdir, "flame_output")
mkpath(folder)

flamegraph_jls_file = joinpath(folder, "flamegraph.jls")
open(flamegraph_jls_file, "w") do io
    serialize(io, fg)
    println("File written successfully")
end
# read back:
read_back = false
if read_back
    begin
        using Serialization
        using ProfileView
        thisdir = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/perf/")
        folder = joinpath(thisdir, "flame_output")
        flamegraph_jls_file = joinpath(folder, "flamegraph.jls")
        fg = open(flamegraph_jls_file, "r") do io
            deserialize(io)
        end
        # ProfileView.view(fg) # looks good, even without initial compiled run
        PProf.pprof(fg; web=true, webhost="127.0.0.1", webport=57600)
    end
end


# It would have been nice to auto-generate these flame graphs
# as a part of CI, but they're really large and slow to load / navigate.
# ProfileView works much better.
import ProfileSVG
folder = joinpath(thisdir, "flame_output")
mkpath(folder)
ProfileSVG.save(joinpath(folder, "flame.svg"), fg; maxframes = 40000, maxdepth = 100)


# == Helping us Fix Inference Triggers == #
itrigs = SnoopCompile.inference_triggers(tinf)
SnoopCompile.suggest(itrigs[end])

SnoopCompile.accumulate_by_source(SnoopCompile.flatten(tinf))