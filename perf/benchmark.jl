if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
include("common.jl")
import BenchmarkTools

sim = init_sim("Bomex")
tendencies = copy(sim.state.prog)
update_n(sim, tendencies, 1) # compile first
trial = BenchmarkTools.@benchmark update_n($sim, $tendencies, 1)
show(stdout, MIME("text/plain"), trial)
println()

sim = init_sim("TRMM_LBA")
tendencies = copy(sim.state.prog)
update_n(sim, tendencies, 1) # compile first
trial = BenchmarkTools.@benchmark update_n($sim, $tendencies, 1)
show(stdout, MIME("text/plain"), trial)
println()

sim = init_sim("Rico")
tendencies = copy(sim.state.prog)
update_n(sim, tendencies, 1) # compile first
trial = BenchmarkTools.@benchmark update_n($sim, $tendencies, 1)
show(stdout, MIME("text/plain"), trial)
println()
