include(joinpath(@__DIR__, "common.jl"))
import ProfileView
import Profile

sim = init_sim("Bomex")
(; tendencies, prog, params, TS) = unpack_params(sim)

ProfileView.@profview begin # compile first
    for i in 1:100
        ∑tendencies!(tendencies, prog, params, TS.t)
    end
end

Profile.clear_malloc_data()
ProfileView.@profview begin
    for i in 1:1000
        ∑tendencies!(tendencies, prog, params, TS.t)
    end
end
