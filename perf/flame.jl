import Pkg
Pkg.develop(path = ".")
import Profile

include("common.jl")
case_name = "Bomex"
sim = init_sim(case_name)
(prob, alg, kwargs) = solve_args(sim)
integrator = ODE.init(prob, alg; kwargs...)

ODE.step!(integrator) # force compilation
prof = Profile.@profile begin
    for _ in 1:100_000
        ODE.step!(integrator)
    end
end

import PProf
PProf.pprof()
# http://localhost:57599/ui/flamegraph?tf
# import ProfileView
# ProfileView.view()
