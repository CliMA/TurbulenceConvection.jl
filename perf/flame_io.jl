include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = "Bomex"
namelist = NameList.default_namelist(case_name)
namelist["time_stepping"]["t_max"] = 15 * namelist["time_stepping"]["dt_max"]
namelist["stats_io"]["frequency"] = 0 # io at every step
namelist["meta"]["uuid"] = "01_flame_io"
sim = Simulation1d(namelist)
initialize(sim)
open_files(sim.Stats) # force compilation
close_files(sim.Stats) # force compilation
open_files(sim.Stats)
(prob, alg, kwargs) = solve_args(sim)
integrator = ODE.init(prob, alg; kwargs...)

ODE.step!(integrator) # force compilation
# callbacks not called after first `step!`
# call, so need to call `step!` again:
ODE.step!(integrator) # force compilation

Profile.clear_malloc_data()
prof = Profile.@profile begin
    for _ in 1:10
        ODE.step!(integrator)
    end
end
close_files(sim.Stats)
∑tendencies!
import PProf
PProf.pprof()
# http://localhost:57599/ui/flamegraph?tf
# import ProfileView
# ProfileView.view()
