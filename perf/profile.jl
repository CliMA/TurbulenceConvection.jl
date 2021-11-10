include("common.jl")
import ProfileView
import Profile

sim = init_sim("Bomex")
tendencies = copy(sim.state.prog)

Profile.@profile update_n(sim, tendencies, 100)
Profile.print()
# Profile.print(; format = :flat, sortedby = :count)
ProfileView.@profview update_n(sim, tendencies, 1000) # compile first
ProfileView.@profview update_n(sim, tendencies, 1000)
Profile.clear_malloc_data()
