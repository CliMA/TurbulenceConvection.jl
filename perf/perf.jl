include("common.jl")
import BenchmarkTools

sim = init_sim("Bomex")
tendencies = copy(sim.state.prog)
update_n(sim, tendencies, 1) # compile first
BenchmarkTools.@benchmark update_n($sim, $tendencies, 1)
