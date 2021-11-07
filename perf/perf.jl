include("common.jl")
import BenchmarkTools

sim = init_sim("Bomex")
update_n(sim, 1) # compile first
BenchmarkTools.@benchmark update_n($sim, 1)
