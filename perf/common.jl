if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
const TC = TurbulenceConvection

const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "driver", "Cases.jl"))
include(joinpath(tc_dir, "driver", "parameter_set.jl"))
include(joinpath(tc_dir, "driver", "dycore.jl"))
include(joinpath(tc_dir, "driver", "main.jl"))
import .NameList

update_n(sim, tendencies, N::Int) = update_n(sim, tendencies, Val(N))

function update_n(sim, tendencies, ::Val{N}) where {N}
    grid = sim.grid
    TS = sim.TS
    prog = sim.state.prog
    aux = sim.state.aux
    params = (; edmf = sim.edmf, grid = grid, gm = sim.gm, case = sim.case, TS = TS, aux = aux)
    for i in 1:N
        ∑tendencies!(tendencies, prog, params, TS.t)
    end
    return nothing
end

function init_sim(case_name; skip_io = true, single_timestep = true, prefix = "")
    if single_timestep && skip_io
        @info "Initializing $case_name for single timestep, with no IO."
    elseif single_timestep
        @info "Initializing $case_name for single timestep."
    elseif skip_io
        @info "Initializing $case_name with no IO."
    else
        @info "Initializing $case_name with IO."
    end
    @info "call update_n(sim, tendencies, n) to run update n-times"
    namelist = NameList.default_namelist(case_name)
    if single_timestep
        namelist["time_stepping"]["t_max"] = namelist["time_stepping"]["dt_max"]
    end
    if skip_io
        namelist["stats_io"]["frequency"] = 10.0e10
        namelist["stats_io"]["skip"] = true
    end
    namelist["meta"]["uuid"] = "$(prefix)01"
    sim = Simulation1d(namelist)
    initialize(sim)
    return sim
end
