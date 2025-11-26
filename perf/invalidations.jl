# From: https://timholy.github.io/SnoopCompile.jl/stable/snoopr/
using SnoopCompileCore
invalidations = SnoopCompileCore.@snoop_invalidations begin
    import Pkg
    import TurbulenceConvection

    const tc_dir = pkgdir(TurbulenceConvection)
    include(joinpath(tc_dir, "driver", "main.jl"))
    include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
    import .NameList

    # case_name = "Bomex"
    # println("Running $case_name...")
    # namelist = NameList.default_namelist(case_name)
    # namelist["meta"]["uuid"] = "01_invalidations"

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


    namelist["meta"]["uuid"] = "01_invalidations"

    integrator, ds_tc_filenames, return_code = main(namelist)
end;

# import ReportMetrics
# ReportMetrics.report_invalidations(;
#     job_name = "invalidations",
#     invalidations,
#     process_filename = x -> last(split(x, "packages/")),
# )

# open("invalidations.txt", "w") do io
#     SnoopCompile.report_invalidations(
#         io; 
#         invalidations=invalidations, 
#         n_rows=20, 
#         process_filename=true
#     )
# end


# Define a function for process_filename (identity, or shorten path)
process_filename_fn = x -> last(split(x, "/"))  # prints just the filename

# # File output
# open("invalidations.txt", "w") do io
#     SnoopCompile.report_invalidations(
#         io;
#         invalidations=invalidations,
#         n_rows=20,
#         process_filename=process_filename_fn,
#     )
# end

# Print to REPL
using SnoopCompile: SnoopCompile
SnoopCompile.report_invalidations(
    Base.stdout;
    invalidations=invalidations,
    n_rows=20,
    process_filename=process_filename_fn,
)