import TurbulenceConvection
const TC = TurbulenceConvection

const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "integration_tests", "cli_options.jl"))
include(joinpath(tc_dir, "integration_tests", "overwrite_namelist.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "driver", "main.jl"))
import .NameList

function unpack_params(sim)
    (; aux, calibrate_io, TS, prog, edmf, precip_model, param_set, case, radiation, forcing, surf_params, aux) = sim
    tendencies = similar(prog)
    parent(tendencies) .= 0
    params = (; calibrate_io, edmf, precip_model, param_set, case, radiation, forcing, surf_params, TS, aux)
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
    # if single_timestep
    #     namelist["time_stepping"]["t_max"] = namelist["time_stepping"]["dt_max"]
    # end
    # if skip_io
    #     namelist["stats_io"]["frequency"] = 10.0e10
    #     namelist["stats_io"]["skip"] = true
    # end
    # namelist["meta"]["uuid"] = "$(prefix)01"
    # sim = Simulation1d(namelist)
    # initialize(sim)
    # return sim
    return init_sim(namelist; skip_io = skip_io, single_timestep = single_timestep, prefix = prefix)
end

function init_sim(namelist::Dict; skip_io = true, single_timestep = true, prefix = "")
    case_name = namelist["meta"]["simname"]
    if single_timestep && skip_io
        @info "Initializing `$case_name` for single timestep, with no IO."
    elseif single_timestep
        @info "Initializing `$case_name` for single timestep."
    elseif skip_io
        @info "Initializing `$case_name` with no IO."
    else
        @info "Initializing `$case_name` with IO."
    end
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
