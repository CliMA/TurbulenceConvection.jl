import TurbulenceConvection
const TC = TurbulenceConvection

const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "driver", "main.jl"))
import .NameList

function unpack_params(sim)
    tendencies = copy(sim.prog)
    TS = sim.TS
    prog = sim.prog
    calibrate_io = sim.calibrate_io
    aux = sim.aux
    params = (;
        calibrate_io,
        edmf = sim.edmf,
        precip_model = sim.precip_model,
        param_set = sim.param_set,
        case = sim.case,
        radiation = sim.radiation,
        forcing = sim.forcing,
        surf_params = sim.surf_params,
        TS = TS,
        aux = aux,
    )
    return (; tendencies, prog, params, TS)
end

function init_sim(case_name; skip_io = true, single_timestep = true, prefix = "")
    if single_timestep && skip_io
        @info "Initializing `$case_name` for single timestep, with no IO."
    elseif single_timestep
        @info "Initializing `$case_name` for single timestep."
    elseif skip_io
        @info "Initializing `$case_name` with no IO."
    else
        @info "Initializing `$case_name` with IO."
    end
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
