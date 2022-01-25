import Pkg
Pkg.develop(path = ".")
import Profile

include("common.jl")
case_name = "Bomex"
sim = init_sim(case_name)

prog = sim.state.prog
tendencies = copy(prog)
params = (; edmf = sim.edmf, grid = sim.grid, gm = sim.gm, case = sim.case, TS = sim.TS, aux = sim.state.aux)

∑tendencies!(tendencies, prog, params, sim.TS.t) # force compilation
prof = Profile.@profile begin
    for _ in 1:100_000
        ∑tendencies!(tendencies, prog, params, sim.TS.t)
    end
end

import PProf
PProf.pprof()
# http://localhost:57599/ui/flamegraph?tf
# import ProfileView
# ProfileView.view()
