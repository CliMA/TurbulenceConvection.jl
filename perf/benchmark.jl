include(joinpath(@__DIR__, "common.jl"))
import BenchmarkTools

sim = init_sim("Bomex")
(; tendencies, prog, params, TS) = unpack_params(sim)
∑tendencies!(tendencies, prog, params, TS.t) # compile first
trial = BenchmarkTools.@benchmark ∑tendencies!($tendencies, $prog, $params, $(TS.t))
show(stdout, MIME("text/plain"), trial)
println()

sim = init_sim("TRMM_LBA")
(; tendencies, prog, params, TS) = unpack_params(sim)
∑tendencies!(tendencies, prog, params, TS.t) # compile first
trial = BenchmarkTools.@benchmark ∑tendencies!($tendencies, $prog, $params, $(TS.t))
show(stdout, MIME("text/plain"), trial)
println()

sim = init_sim("Rico")
(; tendencies, prog, params, TS) = unpack_params(sim)
∑tendencies!(tendencies, prog, params, TS.t) # compile first
trial = BenchmarkTools.@benchmark ∑tendencies!($tendencies, $prog, $params, $(TS.t))
show(stdout, MIME("text/plain"), trial)
println()
