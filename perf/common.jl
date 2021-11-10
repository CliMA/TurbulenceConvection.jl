import TurbulenceConvection

const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "integration_tests", "utils", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "Cases.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "parameter_set.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "main.jl"))
import .NameList

TurbulenceConvection.initialize_io(sim::Simulation1d) = nothing
TurbulenceConvection.io(sim::Simulation1d) = nothing

update_n(sim, N::Int) = update_n(sim, Val(N))

function update_n(sim, ::Val{N}) where {N}
    grid = sim.grid
    TS = sim.TS
    prog = sim.state.prog
    aux = sim.state.aux
    tendencies = sim.state.tendencies
    params = (; edmf = sim.Turb, grid = grid, gm = sim.GMV, case = sim.Case, TS = TS, aux = aux)
    for i in 1:N
        TC.step!(tendencies, prog, params, TS.t)
    end
    return nothing
end

function init_sim(case_name)
    @info "Initializing $case_name for single timestep, with no IO."
    @info "call update_n(sim, n) to run update n-times"
    namelist = NameList.default_namelist(case_name)
    namelist["time_stepping"]["t_max"] = namelist["time_stepping"]["dt"]
    namelist["stats_io"]["frequency"] = 10.0e10
    namelist["stats_io"]["skip"] = true
    namelist["meta"]["uuid"] = "01"
    sim = Simulation1d(namelist)

    Cases.initialize_profiles(sim.Case, sim.grid, sim.GMV, sim.state)
    TC.satadjust(sim.GMV, sim.grid, sim.state)

    Cases.initialize_surface(sim.Case, sim.grid, sim.state, sim.param_set)
    Cases.initialize_forcing(sim.Case, sim.grid, sim.state, sim.GMV, sim.param_set)
    Cases.initialize_radiation(sim.Case, sim.grid, sim.state, sim.GMV, sim.param_set)

    initialize_edmf(sim.Turb, sim.grid, sim.state, sim.Case, sim.GMV, sim.TS)
    return sim
end
