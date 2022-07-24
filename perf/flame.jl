include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = "Bomex"
sim = init_sim(case_name)
(prob, alg, kwargs) = solve_args(sim)
integrator = ODE.init(prob, alg; kwargs...)

function do_work!(integrator)
    for _ in 1:1000
        ODE.step!(integrator)
    end
end

do_work!(integrator) # force compilation
Profile.clear_malloc_data()
prof = Profile.@profile begin
    do_work!(integrator)
end

import PProf
PProf.pprof()
# http://localhost:57599/ui/flamegraph?tf
# import ProfileView
# ProfileView.view()
