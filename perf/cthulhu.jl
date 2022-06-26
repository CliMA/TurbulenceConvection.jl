include(joinpath(@__DIR__, "common.jl"))
import SnoopCompileCore
import Cthulhu

case_name = "Bomex"
println("Running $case_name...")
sim = init_sim(case_name)
(prob, alg, kwargs) = solve_args(sim)

sol = Cthulhu.@descend ODE.solve(prob, alg; kwargs...)
