include("common.jl")
import ProfileView
import Profile

Profile.@profile update_n(sim, 100)
Profile.print()
# Profile.print(; format = :flat, sortedby = :count)
ProfileView.@profview update_n(sim, 1000) # compile first
ProfileView.@profview update_n(sim, 1000)
Profile.clear_malloc_data()
