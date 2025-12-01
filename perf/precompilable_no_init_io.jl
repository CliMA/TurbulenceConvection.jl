include(joinpath(@__DIR__, "common.jl"))
using Test

case_name = "Bomex"
println("Running $case_name...")

sim = init_sim(case_name; single_timestep = false, prefix = "pc_no_init_io1")
(prob, alg, kwargs) = solve_args(sim)
integrator = SciMLBase.init(prob, alg; kwargs...)
t_precompile = @elapsed SciMLBase.solve!(integrator)

sim = init_sim(case_name; single_timestep = false, prefix = "pc_no_init_io2")
(prob, alg, kwargs) = solve_args(sim)
integrator = SciMLBase.init(prob, alg; kwargs...)
t_precompiled = @elapsed SciMLBase.solve!(integrator)

@info "Precompiling run: $(t_precompile)"
@info "Precompiled  run: $(t_precompiled)"
@info "precompiled/precompiling: $(t_precompiled/t_precompile))"

@testset "Test runtime" begin
    @test t_precompiled / t_precompile < 0.6
end
