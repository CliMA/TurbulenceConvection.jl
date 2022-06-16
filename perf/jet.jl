include(joinpath(@__DIR__, "common.jl"))
import JET

# Don't error on type-unstable branching in Thermodynamics
import Thermodynamics as TD
TD.error_on_non_convergence() = false

case_name = "Bomex"
sim = init_sim(case_name)
(; tendencies, prog, params, TS) = unpack_params(sim)

∑tendencies!(tendencies, prog, params, TS.t) # compile first

JET.@test_opt broken = true ∑tendencies!(tendencies, prog, params, TS.t)
