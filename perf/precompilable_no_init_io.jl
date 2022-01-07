if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
using Test
import TurbulenceConvection
tc_dir_glob = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir_glob, "perf", "common.jl"))

case_name = "Bomex"
println("Running $case_name...")

sim = init_sim(case_name; single_timestep = false, prefix = "pc_no_init_io1")
(prob, alg, kwargs) = solve_args(sim)
t_precompile = @elapsed ODE.solve(prob, alg; kwargs...)

sim = init_sim(case_name; single_timestep = false, prefix = "pc_no_init_io2")
(prob, alg, kwargs) = solve_args(sim)
t_precompiled = @elapsed ODE.solve(prob, alg; kwargs...)

@info "Precompiling run: $(t_precompile)"
@info "Precompiled  run: $(t_precompiled)"
@info "precompiled/precompiling: $(t_precompiled/t_precompile))"

@testset "Test runtime" begin
    @test t_precompiled / t_precompile < 0.08
end
