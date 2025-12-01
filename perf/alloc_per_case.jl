# Launch with `julia --project --track-allocation=user`
include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = ENV["ALLOCATION_CASE_NAME"]
@info "Recording allocations for `$case_name`"
sim = init_sim(case_name; prefix = "alloc_$case_name")
(prob, alg, kwargs) = solve_args(sim)
integrator = SciMLBase.init(prob, alg; kwargs...)

SciMLBase.step!(integrator) # compile first
Profile.clear_malloc_data()
SciMLBase.step!(integrator)
