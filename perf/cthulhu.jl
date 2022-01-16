if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import SnoopCompileCore
import Cthulhu
import TurbulenceConvection
tc_dir_glob = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir_glob, "perf", "common.jl"))

case_name = "Bomex"
println("Running $case_name...")
sim = init_sim(case_name)
sim.skip_io || open_files(sim.Stats) # #removeVarsHack
(prob, alg, kwargs) = solve_args(sim)

sol = Cthulhu.@descend ODE.solve(prob, alg; kwargs...)

sim.skip_io || close_files(sim.Stats) # #removeVarsHack
