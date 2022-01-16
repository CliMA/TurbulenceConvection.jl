if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
using Test
import SnoopCompileCore
import TurbulenceConvection
tc_dir_glob = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir_glob, "perf", "common.jl"))

case_name = "Bomex"
println("Running $case_name...")
sim = init_sim(case_name)
sim.skip_io || open_files(sim.Stats) # #removeVarsHack
(prob, alg, kwargs) = solve_args(sim)

tinf = SnoopCompileCore.@snoopi_deep begin
    sol = ODE.solve(prob, alg; kwargs...)
end

sim.skip_io || close_files(sim.Stats) # #removeVarsHack

import SnoopCompile # need SnoopCompile to iterate over InferenceTimingNode's

itrigs = SnoopCompile.inference_triggers(tinf)
@show length(itrigs)
mtrigs = SnoopCompile.accumulate_by_source(Method, itrigs)
pmtrigs = SnoopCompile.parcel(mtrigs)
# filtered_mods = (
#     :Base,
#     :OrdinaryDiffEq,
#     :Broadcast,
#     :Printf,
#     :ProgressMeter,
#     :TerminalLoggers,
#     )
# filter!(x->!any(nameof(x.first) == y for y in filtered_mods), pmtrigs)
filter!(x -> (x.first == TurbulenceConvection), pmtrigs)

@testset "Zero inference triggers" begin
    @test isempty(pmtrigs)
end

# If we can't help but have inference triggers:

# tc_trigs = first(pmtrigs).second
# @show length(tc_trigs)
# for tc_trig in tc_trigs
#     println("-------------")
#     summary(tc_trig)
# end

# @testset "Number of inference triggers" begin
#     # These tests are not strictly needed,
#     # but they will hopefully prevent performance
#     # regressions
#     @test length(tc_trigs) ≤ 2
#     @test length(itrigs) ≤ 36
# end
