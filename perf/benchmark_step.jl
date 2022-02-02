include(joinpath(@__DIR__, "common.jl"))
import BenchmarkTools
# Use a NullLogger to avoid highly excessive progress bar printing
Logging.with_logger(Logging.NullLogger()) do

    sim = init_sim("Bomex")
    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)
    trial = BenchmarkTools.@benchmark ODE.step!($integrator)
    show(stdout, MIME("text/plain"), trial)
    println()

    sim = init_sim("TRMM_LBA")
    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)
    trial = BenchmarkTools.@benchmark ODE.step!($integrator)
    show(stdout, MIME("text/plain"), trial)
    println()

    sim = init_sim("Rico")
    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)
    trial = BenchmarkTools.@benchmark ODE.step!($integrator)
    show(stdout, MIME("text/plain"), trial)
    println()
end
