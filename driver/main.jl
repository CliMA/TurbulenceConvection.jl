import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

import UnPack
import JSON
import ArgParse
import TurbulenceConvection

import CloudMicrophysics
const CM = CloudMicrophysics
const CM0 = CloudMicrophysics.Microphysics_0M
const CM1 = CloudMicrophysics.Microphysics_1M

import ClimaCore
const CC = ClimaCore
import SciMLBase

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq
import StaticArrays: SVector

const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))

include(joinpath(tc_dir, "driver", "initial_conditions.jl"))
include(joinpath(tc_dir, "diagnostics", "compute_diagnostics.jl"))
include(joinpath(tc_dir, "driver", "parameter_set.jl"))
include(joinpath(tc_dir, "driver", "Cases.jl"))
include(joinpath(tc_dir, "driver", "dycore.jl"))
import .Cases

struct Simulation1d{IONT, G, S, GM, C, EDMF, D, TIMESTEPPING, STATS, PS}
    io_nt::IONT
    grid::G
    state::S
    gm::GM
    Case::C
    edmf::EDMF
    diagnostics::D
    TS::TIMESTEPPING
    Stats::STATS
    param_set::PS
    skip_io::Bool
    adapt_dt::Bool
    cfl_limit::Float64
    dt_min::Float64
end

function Simulation1d(namelist)
    TC = TurbulenceConvection
    param_set = create_parameter_set(namelist)

    FT = Float64
    skip_io = namelist["stats_io"]["skip"]
    adapt_dt = namelist["time_stepping"]["adapt_dt"]
    cfl_limit = namelist["time_stepping"]["cfl_limit"]
    dt_min = namelist["time_stepping"]["dt_min"]

    grid = TC.Grid(FT(namelist["grid"]["dz"]), namelist["grid"]["nz"])
    Stats = skip_io ? nothing : TC.NetCDFIO_Stats(namelist, grid)
    case = Cases.get_case(namelist)
    ref_params = Cases.reference_params(case, grid, param_set, namelist)

    gm = TC.GridMeanVariables(namelist, grid, param_set)
    Sur = TC.SurfaceBase(Cases.get_surface_type(case); namelist, ref_params)
    Fo = TC.ForcingBase{Cases.get_forcing_type(case)}()
    Rad = TC.RadiationBase{Cases.get_radiation_type(case)}()

    Case = Cases.CasesBase(case, namelist, grid, param_set, Sur, Fo, Rad)
    edmf = TC.EDMF_PrognosticTKE(namelist, grid, param_set)
    TS = TC.TimeStepping(namelist)

    N_up = TC.n_updrafts(edmf)

    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)

    cent_prog_fields() = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars(FT, N_up))
    face_prog_fields() = TC.FieldFromNamedTuple(fspace, face_prognostic_vars(FT, N_up))
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT, N_up))
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT, N_up))
    diagnostic_cent_fields = TC.FieldFromNamedTuple(cspace, cent_diagnostic_vars(FT, N_up))
    diagnostic_face_fields = TC.FieldFromNamedTuple(fspace, face_diagnostic_vars(FT, N_up))

    prog = CC.Fields.FieldVector(cent = cent_prog_fields(), face = face_prog_fields())
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    diagnostics = CC.Fields.FieldVector(cent = diagnostic_cent_fields, face = diagnostic_face_fields)

    # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
    state = TC.State(prog, aux, nothing)

    TC.compute_ref_state!(state, grid, param_set; ref_params...)

    io_nt = (;
        ref_state = TC.io_dictionary_ref_state(),
        aux = TC.io_dictionary_aux(),
        diagnostics = io_dictionary_diagnostics(),
    )

    return Simulation1d(
        io_nt,
        grid,
        state,
        gm,
        Case,
        edmf,
        diagnostics,
        TS,
        Stats,
        param_set,
        skip_io,
        adapt_dt,
        cfl_limit,
        dt_min,
    )
end

function TurbulenceConvection.initialize(sim::Simulation1d, namelist)
    TC = TurbulenceConvection
    state = sim.state
    FT = eltype(sim.grid)
    t = FT(0)
    Cases.initialize_profiles(sim.Case, sim.grid, sim.gm, state)
    satadjust(sim.gm, sim.grid, sim.state)

    Cases.initialize_surface(sim.Case, sim.grid, state, sim.param_set)
    Cases.initialize_forcing(sim.Case, sim.grid, state, sim.gm, sim.param_set)
    Cases.initialize_radiation(sim.Case, sim.grid, state, sim.gm, sim.param_set)

    initialize_edmf(sim.edmf, sim.grid, state, sim.Case, sim.gm, t)

    sim.skip_io && return nothing
    TC.initialize_io(sim.io_nt.ref_state, sim.Stats)
    TC.io(sim.io_nt.ref_state, sim.Stats, state) # since the reference prog is static

    TC.initialize_io(sim.io_nt.aux, sim.Stats)
    TC.initialize_io(sim.io_nt.diagnostics, sim.Stats)

    # TODO: depricate
    TC.initialize_io(sim.gm, sim.Stats)
    TC.initialize_io(sim.Case, sim.Stats)
    TC.initialize_io(sim.edmf, sim.Stats)

    TC.open_files(sim.Stats)
    TC.write_simulation_time(sim.Stats, t)

    TC.io(sim.io_nt.aux, sim.Stats, state)
    TC.io(sim.io_nt.diagnostics, sim.Stats, sim.diagnostics)

    # TODO: depricate
    TC.io(sim.Case, sim.grid, state, sim.Stats)
    TC.io(sim.edmf, sim.grid, state, sim.Stats)
    TC.close_files(sim.Stats)

    return
end

include("callbacks.jl")

function run(sim::Simulation1d; time_run = true)
    TC = TurbulenceConvection
    iter = 0
    grid = sim.grid
    state = sim.state
    prog = state.prog
    aux = state.aux
    TS = sim.TS
    diagnostics = sim.diagnostics
    sim.skip_io || TC.open_files(sim.Stats) # #removeVarsHack

    t_span = (0.0, sim.TS.t_max)
    params = (;
        edmf = sim.edmf,
        grid = grid,
        gm = sim.gm,
        aux = aux,
        io_nt = sim.io_nt,
        case = sim.Case,
        diagnostics = diagnostics,
        TS = sim.TS,
        Stats = sim.Stats,
        skip_io = sim.skip_io,
        adapt_dt = sim.adapt_dt,
        cfl_limit = sim.cfl_limit,
        dt_min = sim.dt_min,
    )

    callback_io = ODE.DiscreteCallback(condition_io, affect_io!; save_positions = (false, false))
    callback_cfl = ODE.DiscreteCallback(condition_every_iter, monitor_cfl!; save_positions = (false, false))
    callback_cfl = sim.edmf.Precip.precipitation_model == "clima_1m" ? (callback_cfl,) : ()
    callback_dtmax = ODE.DiscreteCallback(condition_every_iter, dt_max!; save_positions = (false, false))
    callback_filters = ODE.DiscreteCallback(condition_every_iter, affect_filter!; save_positions = (false, false))
    callback_adapt_dt = ODE.DiscreteCallback(condition_every_iter, adaptive_dt!; save_positions = (false, false))
    callback_adapt_dt = sim.adapt_dt ? (callback_adapt_dt,) : ()

    callbacks = ODE.CallbackSet(callback_adapt_dt..., callback_dtmax, callback_cfl..., callback_filters, callback_io)

    prob = ODE.ODEProblem(TC.∑tendencies!, state.prog, t_span, params; dt = sim.TS.dt)

    # TODO: LES_driven_SCM is currently unstable w.r.t. higher order moments (HOM).
    # So, we tell OrdinaryDiffEq.jl to not perform NaNs check on the solution
    # so that it doesn't abort early (as the HOM prognostic variables are 1-way coupled)
    unstable_check_kwarg(::Cases.LES_driven_SCM) = (; unstable_check = (dt, u, p, t) -> false)
    unstable_check_kwarg(case) = ()

    kwargs = (;
        progress_steps = 100,
        save_start = false,
        saveat = last(t_span),
        callback = callbacks,
        progress = true,
        unstable_check_kwarg(sim.Case.case)...,
        progress_message = (dt, u, p, t) -> t,
    )
    alg = ODE.Euler()

    if time_run
        sol = @timev ODE.solve(prob, alg; kwargs...)
    else
        sol = ODE.solve(prob, alg; kwargs...)
    end

    sim.skip_io || TC.close_files(sim.Stats) # #removeVarsHack
    if first(sol.t) == sim.TS.t_max
        return :success
    else
        return :simulation_aborted
    end
end

main(namelist; kwargs...) = @timev main1d(namelist; kwargs...)

nc_results_file(stats::TC.NetCDFIO_Stats) = stats.path_plus_file
nc_results_file(::Nothing) = @info "The simulation was run without IO, so no nc files were exported"

function main1d(namelist; time_run = true)
    # TODO: generalize convesion of arrays from namelist to `SVector`s.
    _p = namelist["turbulence"]["EDMF_PrognosticTKE"]["general_ent_params"]
    namelist["turbulence"]["EDMF_PrognosticTKE"]["general_ent_params"] = SVector{length(_p)}(_p)
    sim = Simulation1d(namelist)
    TurbulenceConvection.initialize(sim, namelist)
    if time_run
        return_code = @timev run(sim; time_run)
    else
        return_code = run(sim)
    end
    return_code == :success && println("The simulation has completed.")
    return nc_results_file(sim.Stats), return_code
end


function parse_commandline()
    s = ArgParse.ArgParseSettings(; description = "Run case input")

    ArgParse.@add_arg_table! s begin
        "case_name"
        help = "The case name"
        arg_type = String
        required = true
    end

    return ArgParse.parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__

    args = parse_commandline()
    case_name = args["case_name"]

    namelist = open("namelist_" * "$case_name.in", "r") do io
        JSON.parse(io; dicttype = Dict, inttype = Int64)
    end
    main(namelist)
end
