include(joinpath(@__DIR__, "common.jl"))
import BenchmarkTools

sim = init_sim("Bomex"; prefix = "bm_Bomex", single_timestep = false)
(; tendencies, prog, params, TS) = unpack_params(sim)
∑tendencies!(tendencies, prog, params, TS.t) # compile first
b = BenchmarkTools.@benchmarkable ∑tendencies!($tendencies, $prog, $params, $(TS.t))
trial = BenchmarkTools.run(b, samples = floor(Int, TS.t_max / TS.dt))
show(stdout, MIME("text/plain"), trial)
println()

sim = init_sim("TRMM_LBA"; prefix = "bm_TRMM_LBA", single_timestep = false)
(; tendencies, prog, params, TS) = unpack_params(sim)
∑tendencies!(tendencies, prog, params, TS.t) # compile first
b = BenchmarkTools.@benchmarkable ∑tendencies!($tendencies, $prog, $params, $(TS.t))
trial = BenchmarkTools.run(b, samples = floor(Int, TS.t_max / TS.dt))
show(stdout, MIME("text/plain"), trial)
println()

sim = init_sim("Rico"; prefix = "bm_Rico", single_timestep = false)
(; tendencies, prog, params, TS) = unpack_params(sim)
∑tendencies!(tendencies, prog, params, TS.t) # compile first
b = BenchmarkTools.@benchmarkable ∑tendencies!($tendencies, $prog, $params, $(TS.t))
trial = BenchmarkTools.run(b, samples = floor(Int, TS.t_max / TS.dt))
show(stdout, MIME("text/plain"), trial)
println()
