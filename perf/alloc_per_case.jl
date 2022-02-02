# Launch with `julia --project --track-allocation=user`
include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = ENV["ALLOCATION_CASE_NAME"]
@info "Recording allocations for $case_name"
sim = init_sim(case_name)
tendencies = copy(sim.state.prog)
update_n(sim, tendencies, 1) # compile first
Profile.clear_malloc_data()
update_n(sim, tendencies, 1)
