import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

import UnPack
import JSON
import ArgParse
import TurbulenceConvection

import CloudMicrophysics
const CM = CloudMicrophysics
const CMNe = CloudMicrophysics.MicrophysicsNonEq
const CM0 = CloudMicrophysics.Microphysics0M
const CM1 = CloudMicrophysics.Microphysics1M

import ClimaCore
const CC = ClimaCore
import SciMLBase

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq
import StochasticDiffEq
const SDE = StochasticDiffEq
import StaticArrays: SVector

const tc_dir = pkgdir(TurbulenceConvection)

include(joinpath(tc_dir, "driver", "NetCDFIO.jl"))
include(joinpath(tc_dir, "driver", "initial_conditions.jl"))
include(joinpath(tc_dir, "driver", "compute_diagnostics.jl"))
include(joinpath(tc_dir, "driver", "parameter_set.jl"))
include(joinpath(tc_dir, "driver", "Cases.jl"))
include(joinpath(tc_dir, "driver", "dycore.jl"))
include(joinpath(tc_dir, "driver", "TimeStepping.jl"))
include(joinpath(tc_dir, "driver", "Surface.jl"))
import .Cases

import DiffEqNoiseProcess
import Random
function DiffEqNoiseProcess.wiener_randn!(rng::Random.AbstractRNG, rand_vec::CC.Fields.FieldVector)
    # TODO: fix this hack. ClimaCore's axes(::FieldVector) cannot
    #       return a space (since it has multiple spaces)
    parent(rand_vec.cent) .= Random.randn.()
    parent(rand_vec.face) .= Random.randn.()
end

struct Simulation1d{IONT, P, A, C, EDMF, PM, D, TIMESTEPPING, STATS, PS, SRS}
    io_nt::IONT
    prog::P
    aux::A
    case::C
    edmf::EDMF
    precip_model::PM
    diagnostics::D
    TS::TIMESTEPPING
    Stats::STATS
    param_set::PS
    skip_io::Bool
    calibrate_io::Bool
    adapt_dt::Bool
    cfl_limit::Float64
    dt_min::Float64
    truncate_stack_trace::Bool
    surf_ref_state::SRS
end

function open_files(sim::Simulation1d)
    for inds in TC.iterate_columns(sim.prog.cent)
        open_files(sim.Stats[inds...])
    end
end
function close_files(sim::Simulation1d)
    for inds in TC.iterate_columns(sim.prog.cent)
        close_files(sim.Stats[inds...])
    end
end

include("common_spaces.jl")

function Simulation1d(namelist)
    TC = TurbulenceConvection
    param_set = create_parameter_set(namelist)

    FT = Float64
    skip_io = namelist["stats_io"]["skip"]
    calibrate_io = namelist["stats_io"]["calibrate_io"]
    adapt_dt = namelist["time_stepping"]["adapt_dt"]
    cfl_limit = namelist["time_stepping"]["cfl_limit"]
    dt_min = namelist["time_stepping"]["dt_min"]
    truncate_stack_trace = namelist["logging"]["truncate_stack_trace"]

    cspace, fspace, svpc_space = get_spaces(namelist, param_set, FT)

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
        TC.Clima0M()
    elseif precip_name == "clima_1m"
        TC.Clima1M()
    else
        error("Invalid precip_name $(precip_name)")
    end

    edmf = TC.EDMFModel(namelist, precip_model)
    isbits(edmf) || error("Something non-isbits was added to edmf and needs to be fixed.")
    N_up = TC.n_updrafts(edmf)

    cent_prog_fields = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars, FT, edmf)
    face_prog_fields = TC.FieldFromNamedTuple(fspace, face_prognostic_vars, FT, edmf)
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT, edmf))
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT, edmf))
    diagnostic_cent_fields = TC.FieldFromNamedTuple(cspace, cent_diagnostic_vars(FT, edmf))
    diagnostic_face_fields = TC.FieldFromNamedTuple(fspace, face_diagnostic_vars(FT, edmf))
    diagnostics_single_value_per_col =
        TC.FieldFromNamedTuple(svpc_space, single_value_per_col_diagnostic_vars(FT, edmf))

    prog = CC.Fields.FieldVector(cent = cent_prog_fields, face = face_prog_fields)
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    diagnostics = CC.Fields.FieldVector(;
        cent = diagnostic_cent_fields,
        face = diagnostic_face_fields,
        svpc = diagnostics_single_value_per_col,
    )

    # TODO: clean up
    frequency = namelist["stats_io"]["frequency"]
    nc_filename, outpath = nc_fileinfo(namelist)
    nc_filename_suffix = if namelist["config"] == "column"
        (fn, inds) -> fn
    elseif namelist["config"] == "sphere"
        (fn, inds) -> first(split(fn, ".nc")) * "_col=$(inds).nc"
    end

    Stats = if skip_io
        nothing
    else
        map(TC.iterate_columns(prog.cent)) do inds
            col_state = TC.column_prog_aux(prog, aux, inds...)
            grid = TC.Grid(col_state)
            ncfn = nc_filename_suffix(nc_filename, inds)
            NetCDFIO_Stats(ncfn, frequency, grid)
        end
    end
    casename = namelist["meta"]["casename"]
    open(joinpath(outpath, "namelist_$casename.in"), "w") do io
        JSON.print(io, namelist, 4)
    end

    case_type = Cases.get_case(namelist)
    surf_ref_state = Cases.surface_ref_state(case_type, param_set, namelist)

    Fo = TC.ForcingBase(case_type, param_set; Cases.forcing_kwargs(case_type, namelist)...)
    Rad = TC.RadiationBase(case_type)
    TS = TimeStepping(namelist)

    Ri_bulk_crit = namelist["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"]
    spk = Cases.surface_param_kwargs(case_type, namelist)
    surf_params = Cases.surface_params(case_type, surf_ref_state, param_set; Ri_bulk_crit = Ri_bulk_crit, spk...)
    inversion_type = Cases.inversion_type(case_type)
    case = Cases.CasesBase(case_type; inversion_type, surf_params, Fo, Rad, spk...)

    calibrate_io = namelist["stats_io"]["calibrate_io"]
    aux_dict = calibrate_io ? TC.io_dictionary_aux_calibrate() : TC.io_dictionary_aux()
    diagnostics_dict = calibrate_io ? Dict() : io_dictionary_diagnostics()

    io_nt = (; aux = aux_dict, diagnostics = diagnostics_dict)

    return Simulation1d(
        io_nt,
        prog,
        aux,
        case,
        edmf,
        precip_model,
        diagnostics,
        TS,
        Stats,
        param_set,
        skip_io,
        calibrate_io,
        adapt_dt,
        cfl_limit,
        dt_min,
        truncate_stack_trace,
        surf_ref_state,
    )
end

function initialize(sim::Simulation1d)
    TC = TurbulenceConvection

    (; prog, aux, edmf, case, param_set, skip_io, Stats, io_nt, diagnostics, truncate_stack_trace, surf_ref_state) = sim
    FT = eltype(edmf)
    t = FT(0)

    ts_gm = ["Tsurface", "shf", "lhf", "ustar", "wstar", "lwp_mean", "iwp_mean"]
    ts_edmf = [
        "cloud_base_mean",
        "cloud_top_mean",
        "cloud_cover_mean",
        "env_cloud_base",
        "env_cloud_top",
        "env_cloud_cover",
        "env_lwp",
        "env_iwp",
        "updraft_cloud_cover",
        "updraft_cloud_base",
        "updraft_cloud_top",
        "updraft_lwp",
        "updraft_iwp",
        "rwp_mean",
        "swp_mean",
        "cutoff_precipitation_rate",
        "Hd",
    ]
    ts_list = vcat(ts_gm, ts_edmf)

    # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
    for inds in TC.iterate_columns(prog.cent)
        state = TC.column_prog_aux(prog, aux, inds...)
        diagnostics_col = TC.column_diagnostics(diagnostics, inds...)
        grid = TC.Grid(state)
        compute_ref_state!(state, grid, param_set; ts_g = surf_ref_state)
        if !skip_io
            stats = Stats[inds...]
            NC.Dataset(stats.nc_filename, "a") do ds
                group = "reference"
                add_write_field(ds, "ρ_f", vec(TC.face_aux_grid_mean(state).ρ), group, ("zf",))
                add_write_field(ds, "ρ_c", vec(TC.center_prog_grid_mean(state).ρ), group, ("zc",))
                add_write_field(ds, "p_f", vec(TC.face_aux_grid_mean(state).p), group, ("zf",))
                add_write_field(ds, "p_c", vec(TC.center_aux_grid_mean(state).p), group, ("zc",))
            end
        end
        Cases.initialize_profiles(case, grid, param_set, state)
        set_thermo_state_pθq!(state, grid, edmf.moisture_model, param_set)
        set_grid_mean_from_thermo_state!(param_set, state, grid)
        assign_thermo_aux!(state, grid, edmf.moisture_model, param_set)
        Cases.initialize_forcing(case, grid, state, param_set)
        Cases.initialize_radiation(case, grid, state, param_set)
        initialize_edmf(edmf, grid, state, case, param_set, t)
        skip_io && return nothing
        stats = Stats[inds...]
        initialize_io(stats.nc_filename, io_nt.aux, io_nt.diagnostics)
        initialize_io(stats.nc_filename, ts_list)
        surf = get_surface(case.surf_params, grid, state, t, param_set)
    end

    return
end

include("callbacks.jl")

function solve_args(sim::Simulation1d)
    TC = TurbulenceConvection
    calibrate_io = sim.calibrate_io
    prog = sim.prog
    aux = sim.aux
    TS = sim.TS
    diagnostics = sim.diagnostics

    t_span = (0.0, sim.TS.t_max)
    params = (;
        calibrate_io,
        edmf = sim.edmf,
        precip_model = sim.precip_model,
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
    callback_dtmax = ODE.DiscreteCallback(condition_every_iter, dt_max!; save_positions = (false, false))
    callback_filters = ODE.DiscreteCallback(condition_every_iter, affect_filter!; save_positions = (false, false))
    callback_adapt_dt = ODE.DiscreteCallback(condition_every_iter, adaptive_dt!; save_positions = (false, false))
    callback_adapt_dt = sim.adapt_dt ? (callback_adapt_dt,) : ()
    if sim.edmf.entr_closure isa TC.PrognosticNoisyRelaxationProcess && sim.adapt_dt
        @warn("The prognostic noisy relaxation process currently uses a Euler-Maruyama time stepping method,
              which does not support adaptive time stepping. Adaptive time stepping disabled.")
        callback_adapt_dt = ()
    end

    callbacks = ODE.CallbackSet(callback_adapt_dt..., callback_dtmax, callback_cfl..., callback_filters, callback_io...)

    if sim.edmf.entr_closure isa TC.PrognosticNoisyRelaxationProcess
        prob = SDE.SDEProblem(∑tendencies!, ∑stoch_tendencies!, prog, t_span, params; dt = sim.TS.dt)
        alg = SDE.EM()
    else
        prob = ODE.ODEProblem(∑tendencies!, prog, t_span, params; dt = sim.TS.dt)
        alg = ODE.Euler()
    end

    # TODO: LES_driven_SCM is currently unstable w.r.t. higher order moments (HOM).
    # So, we tell OrdinaryDiffEq.jl to not perform NaNs check on the solution
    # so that it doesn't abort early (as the HOM prognostic variables are 1-way coupled)
    unstable_check_kwarg(::Cases.LES_driven_SCM) = (; unstable_check = (dt, u, p, t) -> false)
    unstable_check_kwarg(case) = ()

    kwargs = (;
        progress_steps = 100,
        save_start = true,
        saveat = last(t_span),
        callback = callbacks,
        progress = true,
        unstable_check_kwarg(sim.case.case)...,
        progress_message = (dt, u, p, t) -> t,
    )
    return (prob, alg, kwargs)
end

function run(sim::Simulation1d; time_run = true)
    TC = TurbulenceConvection
    sim.skip_io || open_files(sim)
    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)

    local sol
    try
        if time_run
            sol = @timev ODE.solve!(integrator)
        else
            sol = ODE.solve!(integrator)
        end
    catch e
        if sim.truncate_stack_trace
            @error "TurbulenceConvection simulation crashed. $(e)"
        else
            @error "TurbulenceConvection simulation crashed. Stacktrace for failed simulation" exception =
                (e, catch_backtrace())
        end
        # "Stacktrace for failed simulation" exception = (e, catch_backtrace())
        return (integrator, :simulation_crashed)
    finally
        sim.skip_io || close_files(sim)
    end

    if first(sol.t) == sim.TS.t_max
        return (integrator, :success)
    else
        return (integrator, :simulation_aborted)
    end
end

main(namelist; kwargs...) = @timev main1d(namelist; kwargs...)

nc_results_file(stats) = map(x -> x.nc_filename, stats)
nc_results_file(stats::NetCDFIO_Stats) = stats.nc_filename
function nc_results_file(::Nothing)
    @info "The simulation was run without IO, so no nc files were exported"
    return ""
end

to_svec(x::AbstractArray) = SA.SVector{length(x)}(x)
to_svec(x::Tuple) = SA.SVector{length(x)}(x)

function main1d(namelist; time_run = true)
    edmf_turb_dict = namelist["turbulence"]["EDMF_PrognosticTKE"]
    for key in keys(edmf_turb_dict)
        entry = edmf_turb_dict[key]
        if entry isa AbstractArray || entry isa Tuple
            edmf_turb_dict[key] = to_svec(entry)
        end
    end
    sim = Simulation1d(namelist)
    initialize(sim)
    if time_run
        (Y, return_code) = @timev run(sim; time_run)
    else
        (Y, return_code) = run(sim; time_run)
    end
    return_code == :success && @info "The simulation has completed."
    return Y, nc_results_file(sim.Stats), return_code
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
