include(joinpath(@__DIR__, "common.jl"))
import BenchmarkTools
# Use a NullLogger to avoid highly excessive progress bar printing
Logging.with_logger(Logging.NullLogger()) do

    sim = init_sim("Bomex"; prefix = "bm_step_Bomex", single_timestep = false)
    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)
    ODE.step!(integrator) # compile first
    b = BenchmarkTools.@benchmarkable ODE.step!($integrator)
    tfinal = integrator.sol.prob.tspan[2]
    trial = BenchmarkTools.run(b, samples = floor(Int, tfinal / integrator.dt))
    show(stdout, MIME("text/plain"), trial)
    println()

    sim = init_sim("TRMM_LBA"; prefix = "bm_step_TRMM_LBA", single_timestep = false)
    (prob, alg, kwargs) = solve_args(sim)
    ODE.step!(integrator) # compile first
    b = BenchmarkTools.@benchmarkable ODE.step!($integrator)
    tfinal = integrator.sol.prob.tspan[2]
    trial = BenchmarkTools.run(b, samples = floor(Int, tfinal / integrator.dt))
    show(stdout, MIME("text/plain"), trial)
    println()

    sim = init_sim("Rico"; prefix = "bm_step_Rico", single_timestep = false)
    (prob, alg, kwargs) = solve_args(sim)
    ODE.step!(integrator) # compile first
    b = BenchmarkTools.@benchmarkable ODE.step!($integrator)
    tfinal = integrator.sol.prob.tspan[2]
    trial = BenchmarkTools.run(b, samples = floor(Int, tfinal / integrator.dt))
    show(stdout, MIME("text/plain"), trial)
    println()
end
