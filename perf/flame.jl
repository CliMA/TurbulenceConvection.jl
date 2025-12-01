include(joinpath(dirname(@__DIR__), "integration_tests", "cli_options.jl"))
if !@isdefined parsed_args
    (s, parsed_args) = parse_commandline()
end
job_id = parsed_args["job_id"]

include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = parsed_args["case"]
sim = init_sim(case_name)
(prob, alg, kwargs) = solve_args(sim)
integrator = SciMLBase.init(prob, alg; kwargs...)

function do_work!(integrator)
    Logging.with_logger(Logging.NullLogger()) do
        for _ in 1:1000
            SciMLBase.step!(integrator)
        end
    end
end

do_work!(integrator) # force compilation
Profile.clear_malloc_data()
prof = Profile.@profile begin
    do_work!(integrator)
end

import ProfileCanvas

if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = job_id
    mkdir(output_dir)
    ProfileCanvas.html_file(joinpath(output_dir, "flame.html"))
else
    ProfileCanvas.view(Profile.fetch())
end
