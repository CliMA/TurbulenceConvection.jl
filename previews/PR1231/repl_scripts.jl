const tc_dir = joinpath(@__DIR__, "..", "..")
include(joinpath(tc_dir, "integration_tests", "cli_options.jl"))
using PrettyTables
(s, parsed_args) = parse_commandline();

buildkite_commands = readlines(joinpath(tc_dir, ".buildkite", "pipeline.yml"));
filter!(x -> occursin("driver.jl", x), buildkite_commands)

@assert length(buildkite_commands) > 0

buildkite_flags = Dict()
for bkcs in buildkite_commands
    println("### Buildkite `$bkcs`")
    print_repl_script(bkcs)
    println()
end
