include(joinpath(dirname(@__DIR__), "integration_tests", "cli_options.jl"))
if !@isdefined parsed_args
    (s, parsed_args) = parse_commandline()
end
job_id = parsed_args["job_id"]

include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = job_id # TODO: make these orthogonal
sim = init_sim(case_name)
(prob, alg, kwargs) = solve_args(sim)
integrator = ODE.init(prob, alg; kwargs...)

function do_work!(integrator)
    Logging.with_logger(Logging.NullLogger()) do
        for _ in 1:1000
            ODE.step!(integrator)
        end
    end
end

do_work!(integrator) # force compilation
Profile.clear_malloc_data()
prof = Profile.@profile begin
    do_work!(integrator)
end

if !isempty(get(ENV, "CI_PERF_CPUPROFILE", ""))

    import ChromeProfileFormat
    output_path = job_id
    mkpath(output_path)
    cpufile = "flame.cpuprofile"
    ChromeProfileFormat.save_cpuprofile(joinpath(output_path, cpufile))

    if !isempty(get(ENV, "BUILDKITE", ""))
        import URIs

        print_link_url(url) = print("\033]1339;url='$(url)'\a\n")

        profiler_url(uri) = URIs.URI("https://profiler.firefox.com/from-url/$(URIs.escapeuri(uri))")

        # copy the file to the clima-ci bucket
        buildkite_pipeline = ENV["BUILDKITE_PIPELINE_SLUG"]
        buildkite_buildnum = ENV["BUILDKITE_BUILD_NUMBER"]
        buildkite_step = ENV["BUILDKITE_STEP_KEY"]

        profile_uri = "$buildkite_pipeline/build/$buildkite_buildnum/$buildkite_step/$cpufile"
        gs_profile_uri = "gs://clima-ci/$profile_uri"
        dl_profile_uri = "https://storage.googleapis.com/clima-ci/$profile_uri"

        # sync to bucket
        Base.run(`gsutil cp $(joinpath(output_path, cpufile)) $gs_profile_uri`)

        # print link
        println("+++ Profiler link for '$profile_uri': ")
        print_link_url(profiler_url(dl_profile_uri))
    end
else
    import PProf
    PProf.pprof()
    # http://localhost:57599/ui/flamegraph?tf
end
