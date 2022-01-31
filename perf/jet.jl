include(joinpath(@__DIR__, "common.jl"))
import JET

case_name = "Bomex"
sim = init_sim(case_name)
prog = sim.state.prog
tendencies = copy(prog)
params = (; edmf = sim.edmf, grid = sim.grid, gm = sim.gm, case = sim.case, TS = sim.TS, aux = sim.state.aux)
JET.@test_opt ∑tendencies!(tendencies, prog, params, sim.TS.t)
