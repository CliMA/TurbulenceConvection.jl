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

const tc_dir = pkgdir(TurbulenceConvection)

struct DiffusivityModel{FT}
    diffusivity::FT
    function DiffusivityModel(namelist)
        diffusivity = namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"]
        return new{typeof(diffusivity)}(diffusivity)
    end
end
abstract type AbstractPrecipitationModel end
struct NoPrecipitation <: AbstractPrecipitationModel end
struct CutoffPrecipitation <: AbstractPrecipitationModel end
struct Clima1M <: AbstractPrecipitationModel end


include(joinpath(tc_dir, "driver", "NetCDFIO.jl"))
include(joinpath(tc_dir, "driver", "initial_conditions.jl"))
include(joinpath(tc_dir, "driver", "compute_diagnostics.jl"))
include(joinpath(tc_dir, "driver", "parameter_set.jl"))
include(joinpath(tc_dir, "driver", "Cases.jl"))
include(joinpath(tc_dir, "driver", "dycore.jl"))
include(joinpath(tc_dir, "driver", "TimeStepping.jl"))
include(joinpath(tc_dir, "driver", "Surface.jl"))
import .Cases

struct Simulation1d{IONT, G, S, C, TCModel, PM, D, TIMESTEPPING, STATS, PS}
    io_nt::IONT
    grid::G
    state::S
    case::C
    turb_conv::TCModel
    precip_model::PM
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
    turb_conv_model = namelist["turbulence"]["turbulence_convection_model"]

    Δz = FT(namelist["grid"]["dz"])
    nz = namelist["grid"]["nz"]
    z₀, z₁ = FT(0), FT(nz * Δz)
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(z₀),
        CC.Geometry.ZPoint{FT}(z₁),
        boundary_tags = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain, nelems = nz)
    grid = TC.Grid(mesh)

    Stats = skip_io ? nothing : NetCDFIO_Stats(namelist, grid)
    case_type = Cases.get_case(namelist)
    surf_ref_state = Cases.surface_ref_state(case_type, param_set, namelist)

    Fo = TC.ForcingBase(case_type, param_set; Cases.forcing_kwargs(case_type, namelist)...)
    Rad = TC.RadiationBase(case_type)
    TS = TimeStepping(namelist)

    # Create the class for precipitation

    precip_name = TC.parse_namelist(
        namelist,
        "microphysics",
        "precipitation_model";
        default = "None",
        valid_options = ["None", "cutoff", "clima_1m"],
    )

    precip_model = if precip_name == "None"
        TC.NoPrecipitation()
    elseif precip_name == "cutoff"
        TC.CutoffPrecipitation()
    elseif precip_name == "clima_1m"
        TC.Clima1M()
    else
        error("Invalid precip_name $(precip_name)")
    end

    edmf = TC.EDMFModel(namelist, precip_model)
    isbits(edmf) || error("Something non-isbits was added to edmf and needs to be fixed.")
    N_up = TC.n_updrafts(edmf)

    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)

    if turb_conv_model == "eddy_diffusivity"
        turb_conv = DiffusivityModel(namelist)
        # N_up = nothing
        cent_prog_fields = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars_gm(FT))
        face_prog_fields = TC.FieldFromNamedTuple(fspace, face_prognostic_vars_gm(FT))
        aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT))
        aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT))
        diagnostic_cent_fields = TC.FieldFromNamedTuple(cspace, cent_diagnostic_vars_gm(FT))
        diagnostic_face_fields = TC.FieldFromNamedTuple(fspace, face_diagnostic_vars_gm(FT))
    elseif turb_conv_model == "EDMF"
        turb_conv = TC.EDMFModel(namelist)
        N_up = TC.n_updrafts(turb_conv)
        cent_prog_fields = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars(FT, N_up))
        face_prog_fields = TC.FieldFromNamedTuple(fspace, face_prognostic_vars(FT, N_up))
        aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT, N_up))
        aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT, N_up))
        diagnostic_cent_fields = TC.FieldFromNamedTuple(cspace, cent_diagnostic_vars(FT, N_up))
        diagnostic_face_fields = TC.FieldFromNamedTuple(fspace, face_diagnostic_vars(FT, N_up))
    end
    isbits(turb_conv) || error("Something non-isbits was added to edmf and needs to be fixed.")

    prog = CC.Fields.FieldVector(cent = cent_prog_fields, face = face_prog_fields)
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    diagnostics = CC.Fields.FieldVector(cent = diagnostic_cent_fields, face = diagnostic_face_fields)

    # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
    state = TC.State(prog, aux, nothing)
    compute_ref_state!(state, grid, param_set; ts_g = surf_ref_state)

    Ri_bulk_crit = namelist["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"]
    spk = Cases.surface_param_kwargs(case_type, namelist)
    surf_params = Cases.surface_params(case_type, grid, surf_ref_state, param_set; Ri_bulk_crit = Ri_bulk_crit, spk...)
    inversion_type = Cases.inversion_type(case_type)
    case = Cases.CasesBase(case_type; inversion_type, surf_params, Fo, Rad, spk...)

    io_nt = (;
        ref_state = TC.io_dictionary_ref_state(),
        aux = TC.io_dictionary_aux(),
        diagnostics = io_dictionary_diagnostics(),
    )

    return Simulation1d(
        io_nt,
        grid,
        state,
        case,
        turb_conv,
        precip_model,
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

function initialize(sim::Simulation1d)
    TC = TurbulenceConvection
    state = sim.state
    FT = eltype(sim.grid)
    t = FT(0)
    Cases.initialize_profiles(sim.case, sim.grid, sim.param_set, state)
    satadjust(sim.param_set, sim.grid, sim.state)

    Cases.initialize_forcing(sim.case, sim.grid, state, sim.param_set)
    Cases.initialize_radiation(sim.case, sim.grid, state, sim.param_set)

    initialize_turb_conv(sim.turb_conv, sim.grid, state, sim.case, sim.param_set, t)

    sim.skip_io && return nothing
    initialize_io(sim.io_nt.ref_state, sim.Stats)
    io(sim.io_nt.ref_state, sim.Stats, state) # since the reference prog is static

    initialize_io(sim.io_nt.aux, sim.Stats)
    initialize_io(sim.io_nt.diagnostics, sim.Stats)

    # TODO: deprecate
    initialize_io(sim.Stats)
    initialize_io(sim.edmf, sim.Stats)

    open_files(sim.Stats)
    write_simulation_time(sim.Stats, t)

    io(sim.io_nt.aux, sim.Stats, state)
    io(sim.io_nt.diagnostics, sim.Stats, sim.diagnostics)

    # TODO: deprecate
    surf = get_surface(sim.case.surf_params, sim.grid, state, t, sim.param_set)
    io(surf, sim.case.surf_params, sim.grid, state, sim.Stats, t)
    close_files(sim.Stats)

    return
end

include("callbacks.jl")

function solve_args(sim::Simulation1d)
    TC = TurbulenceConvection
    grid = sim.grid
    state = sim.state
    prog = state.prog
    aux = state.aux
    TS = sim.TS
    diagnostics = sim.diagnostics

    t_span = (0.0, sim.TS.t_max)
    params = (;
        turb_conv = sim.turb_conv,
        precip_model = sim.precip_model,
        grid = grid,
        param_set = sim.param_set,
        aux = aux,
        io_nt = sim.io_nt,
        case = sim.case,
        diagnostics = diagnostics,
        TS = sim.TS,
        Stats = sim.Stats,
        skip_io = sim.skip_io,
        cfl_limit = sim.cfl_limit,
        dt_min = sim.dt_min,
    )

    callback_io = ODE.DiscreteCallback(condition_io, affect_io!; save_positions = (false, false))
    callback_io = sim.skip_io ? () : (callback_io,)
    callback_cfl = ODE.DiscreteCallback(condition_every_iter, monitor_cfl!; save_positions = (false, false))
    callback_cfl = sim.precip_model isa TC.Clima1M ? (callback_cfl,) : ()
    dt_max! = sim.turb_conv isa TC.EDMFModel ? edmf_dt_max! : diffusivity_dt_max!
    callback_dtmax = ODE.DiscreteCallback(condition_every_iter, dt_max!; save_positions = (false, false))
    callback_filters = ODE.DiscreteCallback(condition_every_iter, affect_filter!; save_positions = (false, false))
    callback_adapt_dt = ODE.DiscreteCallback(condition_every_iter, adaptive_dt!; save_positions = (false, false))
    callback_adapt_dt = sim.adapt_dt ? (callback_adapt_dt,) : ()

    callbacks = ODE.CallbackSet(callback_adapt_dt..., callback_dtmax, callback_cfl..., callback_filters, callback_io...)

    prob = ODE.ODEProblem(∑tendencies!, state.prog, t_span, params; dt = sim.TS.dt)

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
        unstable_check_kwarg(sim.case.case)...,
        progress_message = (dt, u, p, t) -> t,
    )
    alg = ODE.Euler()
    return (prob, alg, kwargs)
end

function run(sim::Simulation1d; time_run = true)
    TC = TurbulenceConvection
    sim.skip_io || open_files(sim.Stats) # #removeVarsHack
    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)
    if time_run
        sol = @timev ODE.solve!(integrator)
    else
        sol = ODE.solve!(integrator)
    end

    sim.skip_io || close_files(sim.Stats) # #removeVarsHack
    if first(sol.t) == sim.TS.t_max
        return :success
    else
        return :simulation_aborted
    end
end

main(namelist; kwargs...) = @timev main1d(namelist; kwargs...)

nc_results_file(stats::NetCDFIO_Stats) = stats.path_plus_file
nc_results_file(::Nothing) = @info "The simulation was run without IO, so no nc files were exported"

function main1d(namelist; time_run = true)
    # TODO: generalize conversion of arrays from namelist to `SVector`s.
    for param_name in ["general_ent_params", "general_stochastic_ent_params", "fno_ent_params"]
        _p = namelist["turbulence"]["EDMF_PrognosticTKE"][param_name]
        namelist["turbulence"]["EDMF_PrognosticTKE"][param_name] = SVector{length(_p)}(_p)
    end
    sim = Simulation1d(namelist)
    initialize(sim)
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
