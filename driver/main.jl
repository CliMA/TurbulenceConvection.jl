import Logging
import ForwardDiff
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

import UnPack
import JSON
import ArgParse
import TurbulenceConvection

import CloudMicrophysics as CM
import CloudMicrophysics.MicrophysicsNonEq as CMNe
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.CommonTypes as CMT

import ClimaCore as CC
import SciMLBase

import OrdinaryDiffEq as ODE
import StochasticDiffEq as SDE
import StaticArrays: SVector
import StaticArrays as SA  # moved out of driver/Surface.jl due to https://github.com/CliMA/SurfaceFluxes.jl/pull/128/commits 

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

struct Simulation1d{IONT, P, A, C, F, R, SM, EDMF, PM, RFM, D, TIMESTEPPING, STATS, PS, SRS, AUXDK}
    io_nt::IONT
    prog::P
    aux::A
    case::C
    forcing::F
    radiation::R
    surf_params::SM
    edmf::EDMF
    precip_model::PM
    rain_formation_model::RFM
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
    aux_data_kwarg::AUXDK
end

function open_files(sim::Simulation1d)
    CC.Fields.bycolumn(axes(sim.prog.cent)) do colidx
        open_files(sim.Stats[colidx])
    end
end
function close_files(sim::Simulation1d)
    CC.Fields.bycolumn(axes(sim.prog.cent)) do colidx
        close_files(sim.Stats[colidx])
    end
end

include("common_spaces.jl")

function Simulation1d(namelist)
    TC = TurbulenceConvection

    FT = namelist["float_type"] == "Float32" ? Float32 : Float64
    # This is used for testing Duals
    FTD = namelist["test_duals"] ? typeof(ForwardDiff.Dual{Nothing}(FT(1), FT(0))) : FT

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    # TODO: the namelist should override whatever
    # we need to override for calibration

    param_set = create_parameter_set(namelist, toml_dict, FTD) # can this be const? or do we edit it somehow later? I think it's supposed to be tuples anyway

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
        TC.Clima1M(param_set, namelist)
    else
        error("Invalid precip_name $(precip_name)")
    end

    if precip_name == "clima_1m"
        rain_formation_name =
            TC.parse_namelist(namelist, "microphysics", "rain_formation_scheme"; default = "clima_1m_default")
        prescribed_Nd = TC.parse_namelist(namelist, "microphysics", "prescribed_Nd"; default = 1e8)

        rain_formation_model = if rain_formation_name == "clima_1m_default"
            TC.Clima1M_default()
        elseif rain_formation_name == "KK2000"
            TC.Clima2M(prescribed_Nd, CMT.KK2000Type())
        elseif rain_formation_name == "B1994"
            TC.Clima2M(prescribed_Nd, CMT.B1994Type())
        elseif rain_formation_name == "TC1980"
            TC.Clima2M(prescribed_Nd, CMT.TC1980Type())
        elseif rain_formation_name == "LD2004"
            TC.Clima2M(prescribed_Nd, CMT.LD2004Type())
        else
            error("Invalid rain_formation_name $(rain_formation_name)")
        end
    else
        rain_formation_model = TC.NoRainFormation()
    end

    edmf = TC.EDMFModel(FTD, namelist, precip_model, rain_formation_model, param_set)
    # RF contains a very large number of parameters,
    # using SVectors / Tuples for this is very expensive
    # for the compiler, so we'll just accept Arrays for now.
    @info "edmf = \n$(summary(edmf))"
    if !isbits(edmf) && !(edmf.ml_entr_closure isa TC.RFEntr)
        @show edmf
        error("Something non-isbits was added to edmf and needs to be fixed.")
    end
    # N_up = TC.n_updrafts(edmf) # not used?

    cent_prog_fields = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars, FT, edmf; calibrate_io=calibrate_io)
    face_prog_fields = TC.FieldFromNamedTuple(fspace, face_prognostic_vars, FT, edmf; calibrate_io=calibrate_io)
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars, FT, edmf; calibrate_io=calibrate_io)
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars, FT, edmf; calibrate_io=calibrate_io)
    diagnostic_cent_fields = TC.FieldFromNamedTuple(cspace, cent_diagnostic_vars, FT, edmf; calibrate_io=calibrate_io)
    diagnostic_face_fields = TC.FieldFromNamedTuple(fspace, face_diagnostic_vars, FT, edmf; calibrate_io=calibrate_io)
    diagnostics_single_value_per_col =
        TC.FieldFromNamedTuple(svpc_space, single_value_per_col_diagnostic_vars(FT, edmf, Val{calibrate_io}()))

    prog = CC.Fields.FieldVector(cent = cent_prog_fields, face = face_prog_fields)
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    diagnostics = CC.Fields.FieldVector(;
        cent = diagnostic_cent_fields,
        face = diagnostic_face_fields,
        svpc = diagnostics_single_value_per_col,
    )

    if namelist["test_duals"]
        prog = ForwardDiff.Dual{Nothing}.(prog, 1.0)
        aux = ForwardDiff.Dual{Nothing}.(aux, 1.0)
        diagnostics = ForwardDiff.Dual{Nothing}.(diagnostics, 1.0)
    end

    # TODO: clean up
    frequency = namelist["stats_io"]["frequency"]
    nc_filename, outpath = nc_fileinfo(namelist)
    nc_filename_suffix = if namelist["config"] == "column"
        (fn, colidx) -> fn
    elseif namelist["config"] == "sphere"
        (fn, colidx) -> first(split(fn, ".nc")) * "_col=$(colidx).nc"
    end

    Stats = if skip_io
        nothing
    else
        colidx_type = TC.column_idx_type(axes(prog.cent))
        stats = Dict{colidx_type, NetCDFIO_Stats{FT}}()
        CC.Fields.bycolumn(axes(prog.cent)) do colidx
            col_state = TC.column_prog_aux(prog, aux, colidx, calibrate_io)
            grid = TC.Grid(col_state)
            ncfn = nc_filename_suffix(nc_filename, colidx)
            stats[colidx] = NetCDFIO_Stats(ncfn, frequency, grid)
        end
        stats
    end
    casename = namelist["meta"]["casename"]
    open(joinpath(outpath, "namelist_$casename.in"), "w") do io
        JSON.print(io, namelist, 4)
    end

    case = Cases.get_case(namelist)
    surf_ref_state = Cases.surface_ref_state(case, param_set, namelist) # this doesn't accept aux_data_kwarg which would be useful for socrates... didn't want to change the call for all cases though...

    forcing = Cases.ForcingBase(case, FT; Cases.forcing_kwargs(case, namelist)...)

    radiation = Cases.RadiationBase(case, FT)
    TS = TimeStepping(FT, namelist)

    @info "TS" TS 

    Ri_bulk_crit::FTD = namelist["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"]


    if case isa Cases.SOCRATES
        # Helper to get the first column so we can get the grid
        first_col = let f = nothing; 
            try
                CC.Fields.bycolumn(axes(prog.cent)) do col
                    f = col
                    throw(:stop)
                end
            catch e
                e === :stop || rethrow()
            end
            f
        end

        sample_state = TC.column_prog_aux(prog, aux, first_col, calibrate_io)
        sample_grid = TC.Grid(sample_state) # Create a saple grid so we can only calculate SSCF forcing functions once
        # @warn "first_col = $first_col; sample_state = $sample_state; sample_grid = $sample_grid;"
        aux_data_kwarg = Cases.aux_data_kwarg(case, namelist, param_set, sample_grid)
    else
        aux_data_kwarg = Cases.aux_data_kwarg(case, namelist)
    end
    surf_params = Cases.surface_params(case, surf_ref_state, param_set; Ri_bulk_crit = Ri_bulk_crit, aux_data_kwarg...)

    calibrate_io = namelist["stats_io"]["calibrate_io"]
    aux_dict = calibrate_io ? TC.io_dictionary_aux_calibrate() : TC.io_dictionary_aux(edmf)
    diagnostics_dict = calibrate_io ? Dict() : io_dictionary_diagnostics()

    io_nt = (; aux = aux_dict, diagnostics = diagnostics_dict)

    return Simulation1d(
        io_nt,
        prog,
        aux,
        case,
        forcing,
        radiation,
        surf_params,
        edmf,
        precip_model,
        rain_formation_model,
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
        aux_data_kwarg,
    )
end

function initialize(sim::Simulation1d)
    TC = TurbulenceConvection

    (; prog, aux, edmf, case, forcing, radiation, surf_params, param_set, surf_ref_state) = sim
    (; aux_data_kwarg, skip_io, Stats, io_nt, calibrate_io, diagnostics, truncate_stack_trace) = sim

    @info "calibrate_io = $calibrate_io"

    ts_gm = ["Tsurface", "qtsurface", "shf", "lhf", "ustar", "wstar", "lwp_mean", "iwp_mean"]
    ts_edmf = calibrate_io ?  ["rwp_mean", "swp_mean"] : [
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
        "integ_total_flux_qt",
        "integ_total_flux_s",
    ]

    ts_list = vcat(ts_gm, ts_edmf)


    # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
    CC.Fields.bycolumn(axes(prog.cent)) do colidx
        state = TC.column_prog_aux(prog, aux, colidx, calibrate_io)
        diagnostics_col = TC.column_diagnostics(diagnostics, colidx)
        grid = TC.Grid(state)
        FT = TC.float_type(state)
        t = FT(0)
        compute_ref_state!(state, param_set; ts_g = surf_ref_state)
        Cases.overwrite_ref_state_from_file!(case, state, param_set)  # if we have an external reference state we may want to do this instead, if isn't socrates should do nothing
        if !skip_io
            stats = Stats[colidx]
            NC.Dataset(stats.nc_filename, "a") do ds
                group = "reference"
                add_write_field(ds, "ρ_f", vec(TC.face_aux_grid_mean(state).ρ), group, ("zf",))
                add_write_field(ds, "ρ_c", vec(TC.center_prog_grid_mean(state).ρ), group, ("zc",))
                add_write_field(ds, "p_f", vec(TC.face_aux_grid_mean(state).p), group, ("zf",))
                add_write_field(ds, "p_c", vec(TC.center_aux_grid_mean(state).p), group, ("zc",))
            end
        end
        Cases.initialize_profiles(case, param_set, state; aux_data_kwarg...)
        set_thermo_state_from_aux!(state, edmf.moisture_model, param_set)
        set_grid_mean_from_thermo_state!(param_set, state)
        assign_thermo_aux!(state, param_set)
        Cases.initialize_forcing(case, forcing, state, param_set; aux_data_kwarg...)
        Cases.initialize_radiation(case, radiation, state, param_set; aux_data_kwarg...)
        initialize_edmf(edmf, state, surf_params, param_set, t, case)
        if !skip_io
            stats = Stats[colidx]
            initialize_io(stats.nc_filename, eltype(grid), io_nt.aux, io_nt.diagnostics)
            initialize_io(stats.nc_filename, eltype(grid), ts_list) # technically w/ calibrate_io this should be pared down to match what's listed in affect_io!()
        end
    end

    (prob, alg, kwargs) = solve_args(sim)
    integrator = ODE.init(prob, alg; kwargs...)
    skip_io && return (integrator, :success)

    open_files(sim)
    try
        affect_io!(integrator)
    catch e
        if truncate_stack_trace
            @warn "IO during initialization failed."
        else
            @warn "IO during initialization failed." exception = (e, catch_backtrace())
        end
        return (integrator, :simulation_crashed)
    finally
        close_files(sim)
    end

    return (integrator, :success)
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
        rain_formation_model = sim.rain_formation_model,
        param_set = sim.param_set,
        aux = aux,
        io_nt = sim.io_nt,
        case = sim.case,
        forcing = sim.forcing,
        surf_params = sim.surf_params,
        radiation = sim.radiation,
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
    callback_adapt_dt = (sim.adapt_dt || sim.TS.spinup_half_t_max > 0 || sim.TS.use_tendency_timestep_limiter) ? (callback_adapt_dt,) : ()
    # callback_spinup_dt = ODE.DiscreteCallback(condition_every_iter, spinup_dt!; save_positions = (false, false)) 
    # callback_spinup_dt = (sim.TS.spinup_half_t_max > 0) ? (callback_spinup_dt,) : () # we condition on spinup_half_t_max > 0 | put in tuple so can splat if it doesnt exist
    if sim.edmf.entr_closure isa TC.PrognosticNoisyRelaxationProcess && sim.adapt_dt
        @warn("The prognostic noisy relaxation process currently uses a Euler-Maruyama time stepping method,
              which does not support adaptive time stepping. Adaptive time stepping disabled.")
        callback_adapt_dt = ()
    end

    # callback_∑tendencies! = sim.TS.use_tendency_timestep_limiter ? (ODE.DiscreteCallback(condition_every_iter, call_∑tendencies!; save_positions = (false, false)),) : ()
    callback_∑tendencies! = sim.TS.use_tendency_timestep_limiter ? (ODE.DiscreteCallback(condition_every_iter, call_∑tendencies!; save_positions = (false, false)),) : (ODE.DiscreteCallback(condition_every_iter, call_∑tendencies!; save_positions = (false, false)),)
    callback_call_update_aux_caller! = sim.TS.use_tendency_timestep_limiter ? () : (ODE.DiscreteCallback(condition_every_iter, call_update_aux_caller!; save_positions = (false, false)),)

    callback_reset_dt = sim.adapt_dt ? ODE.DiscreteCallback(condition_every_iter, reset_dt!; save_positions = (false, false)) : ()
    

    # Maybe we need to edit TS.dt here at the start? idk when the callback is called but seems to be every 100 steps... perhaps?
    # The ode callback will adjust sim.TS.dt but i'm not sure letting dycore do it is the right move...

    # callbacks = ODE.CallbackSet(callback_dtmax, callback_adapt_dt..., callback_cfl..., callback_filters, callback_io...)
    local callbacks::ODE.CallbackSet

    if sim.edmf.entr_closure isa TC.PrognosticNoisyRelaxationProcess
        prob = SDE.SDEProblem(∑tendencies!, ∑stoch_tendencies!, prog, t_span, params; dt = sim.TS.dt)
        alg = SDE.EM()
    else
        if sim.TS.explicit_solver
            if sim.TS.fixed_step_solver
                if sim.TS.use_tendency_timestep_limiter
                    @info "Using ∑tendencies_null! for ODEProblem, callback will calculate tendencies"
                    # Before, we would update_aux then calculate tendencies in ∑tendencies!, then take the step, then calculate a dt in callbacks, then filter then do io in callbacks. 
                    # To match, we should filter then do io in callbacks, then update aux then calculate tendencies in callback_∑tendencies!, then calculate the needed timestep in dt_max/adapt_dt, then take the step
                    # if we left the original order, we'd calculate our tendencies then do all our filtering, then take the step...
                    callbacks = ODE.CallbackSet(callback_filters, callback_io..., callback_reset_dt, callback_∑tendencies!..., callback_dtmax, callback_adapt_dt..., callback_cfl..., )  
                    prob = ODE.ODEProblem{true, SciMLBase.FullSpecialize}(∑tendencies_null!, prog, t_span, params; dt = sim.TS.dt) # we will calculate tendencies in the callback so that we can use them for the most up to date du so we don't need another call to ∑tendencies!, so we use ∑tendencies_null!
                    # The upside is we get the dt we want, the downside is we calculate the tendencies in the dt_max! callback so when the step is taken, we return to the io callback without having dones any filtering/etc
                else
                    #= 
                        We do not have a separate fcn to filter tendencies. 
                        Right now, however, sum_tendencies calls update_aux(), which is relevant for calculating cfl, as well as many subtendency timestep limits.
                        When we call update_aux() however, those functions want the upcoming Δt, before the step... The problem is things are so intertwined in update_aux() that it's hard to separate... for example q_liq and q_ice get set in saturation_adjustment(), but precipitation rates are also calculated here.

                        1) Filter Prognostic
                        ---> Ideally here we'd want a way to back out the new aux tracer values... without calculating tendencies.... [[ we added this as callback_call_update_aux_caller!() ]]
                        2) Calculate dt limits [cfl, adapt dt, dtmax] [ this will not be tendency dependent]
                            - Monitor CFL sets cfl_dt_max which we store in sim.TS.cfl_dt_max ⟹ used in dt_max! and adaptive_dt!
                            - The CFL error check will be thrown if Δt cannot be set fast enough to meet cfl condition for example. If not using adaptive_dt!, 
                        2) ∑tendencies! [which calls update_aux() etc]
                        4) step

                    =#
                    # callbacks = ODE.CallbackSet(callback_dtmax, callback_adapt_dt..., callback_cfl..., callback_filters, callback_io...) # original
                    # prob = ODE.ODEProblem(∑tendencies!, prog, t_span, params; dt = sim.TS.dt)
                    callbacks = ODE.CallbackSet( callback_filters, callback_call_update_aux_caller!..., callback_dtmax, callback_adapt_dt..., callback_cfl..., callback_∑tendencies!..., callback_io...,)
                    # prob = ODE.ODEProblem(∑tendencies_null!, prog, t_span, params; dt = sim.TS.dt) # we will calculate tendencies in the step so we don't need another call to ∑tendencies! so we use ∑tendencies_null!
                    prob = ODE.ODEProblem{true, SciMLBase.FullSpecialize}(∑tendencies_null!, prog, t_span, params; dt = sim.TS.dt) # full specialize to try to avoid inference issues
                end
                alg = if sim.TS.algorithm isa Val{:Euler}
                    ODE.Euler()
                else
                    # alg_string = String(typeof(alg).parameters[1])
                    alg_string = String(TC.Parameters.unwrap_val(alg))
                    error("Unsupported explicit fixed-step algorithm $(alg_string)()")
                end

            else

                callbacks = ODE.CallbackSet(callback_cfl..., callback_filters, callback_io...) # drop the dt stuff since the adaptive solver will do it by itself
                prob = ODE.ODEProblem{true, SciMLBase.FullSpecialize}(∑tendencies_robust!, prog, t_span, params; dt = sim.TS.dt) # use robust version that won't crash when out of domain.

                alg = if sim.TS.algorithm isa Val{:Heun}
                    ODE.Heun()
                elseif sim.TS.algorithm isa Val{:RK4}
                    alg = ODE.RK4()
                elseif sim.TS.algorithm isa Val{:Tsit5}
                    alg = ODE.Tsit5()
                else
                    # alg_string = String(typeof(a).parameters[1])
                    alg_string = String(TC.Parameters.unwrap_val(alg))
                    error("Unsupported explicit adaptive algorithm $(alg_string)()")
                end
            end

        else  # implicit solver
            if sim.TS.fixed_step_solver
                error("No implicit fixed-step solver currently supported")
            else
                callbacks = ODE.CallbackSet(callback_cfl..., callback_filters, callback_io...) # drop the dt stuff since the implicit solver will do it by itself. leave cfl... hopefully the model stability also prevents cfl violations.... (we still also have buoyancy limits)
                prob = ODE.ODEProblem{true, SciMLBase.FullSpecialize}(∑tendencies!, prog, t_span, params; dt = sim.TS.dt)
                if sim.TS.algorithm isa Val{:ImplicitEuler} # adaptive solver
                    alg = ODE.ImplicitEuler(; autodiff = false) # this is far too slow to actually use...
                elseif sim.TS.algorithm isa Val{:KenCarp47}
                    alg = ODE.KenCarp47(; autodiff = false, linsolve = ODE.KrylovJL_GMRES()) # see https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/#Using-Jacobian-Free-Newton-Krylov
                else
                    # alg_string = String(typeof(a).parameters[1])
                    alg_string = String(TC.Parameters.unwrap_val(alg))
                    error("Unsupported implicit algorithm $(alg_string)()")
                end
            end

        end
    end
    
    kwargs = (;
        progress_steps = 100,
        save_start = false,
        saveat = last(t_span),
        callback = callbacks,
        progress = true,
        progress_message = (dt, u, p, t) -> t,
        # my additions
        isoutofdomain = sim.TS.use_isoutofdomain_limiter ? isoutofdomain : ODE.DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN, # use our isoutofdomain if we set it, otherwise use their default which just returns false
        unstable_check = my_unstable_check,
        (isnan(sim.TS.abstol) ? (;) : (;abstol = sim.TS.abstol))...,
        (isnan(sim.TS.reltol) ? (;) : (;reltol = sim.TS.reltol))...,
        (sim.TS.fixed_step_solver ? (;) : (;))..., # use their default for dtmin (set nothing, don't set for fixed-step solver)
        ((!sim.TS.fixed_step_solver &&  (sim.TS.adaptive_depth_limit ≥ 0)) ? (;maxiters = sim.TS.adaptive_depth_limit) : (;) )..., # use our adaptive_depth_limit if we set it, otherwise use their default
        ((!sim.TS.fixed_step_solver && !isnan(sim.TS.dt_max) && !isinf(sim.TS.dt_max)) ? (;dtmax = sim.TS.dt_max) : (;))..., # use our dt_max if we set it, otherwise use their default
        force_dtmin = true, # Allow solver to continue at dtmin failure in adapt
    )
    return (prob, alg, kwargs)
end


"""
    default is any(isnan, u)
"""
function my_unstable_check(dt, u, p, t ) # see https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/
    unstable = false

    # if any(isnan, u)
    if any(!isfinite, u) # this is more strict than isnan, but also catches Inf
        @warn "Unstable: NaN detected in solution, typeof(u) = $(typeof(u)); u.cent = $(u.cent); u.face = $(u.face)"
        println(u)
        flush(stdout); flush(stderr)

        # prog = integrator.u
        # (; edmf, param_set, aux, case, surf_params) = integrator.p

        summary(stdout, u)
        summary(stdout, p.aux)

        flush(stdout); flush(stderr)

        # if isdefined(Main, :Infiltrator)
        #     Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
        # end

        unstable = true
        # unstable = false # test just ignoring
    end
    return unstable
end

function sim_run(sim::Simulation1d, integrator; time_run = true) # named sim_run to not clash with Base.run()
    TC = TurbulenceConvection
    sim.skip_io || open_files(sim)

    local sol
    try
        if time_run
            # sol = @timev ODE.solve!(integrator)
        sol = begin
            t = @timed ODE.solve!(integrator)
            @warn "ODE.solve!() completed in $(t.time) seconds, memory used: $(round(t.bytes / 1024^3, digits=2)) GiB"
            t.value
        end

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

# main(namelist; kwargs...) = @timev main1d(namelist; kwargs...)
main(namelist; kwargs...) = begin
    t = @timed main1d(namelist; kwargs...)
    @warn "main1d() completed in $(t.time) seconds, memory used: $(round(t.bytes / 1024^3, digits=2)) GiB"
    t.value
end


nc_results_file(stats) = map(x -> stats[x].nc_filename, collect(keys(stats)))
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
    if time_run
        # (integrator, return_init) = @timev initialize(sim)
        t = @timed initialize(sim)
        @warn "initialize() completed in $(t.time) seconds, memory used: $(round(t.bytes / 1024^3, digits=2)) GiB"
        (integrator, return_init) = t.value
    else
        (integrator, return_init) = initialize(sim)
    end    
    return_init === :success && @info "The initialization has completed."

    if time_run
        (Y, return_code) = begin
            # (Y, return_code) = @timev sim_run(sim, integrator; time_run)
            t = @timed sim_run(sim, integrator; time_run)
            @warn "sim_run() completed in $(t.time) seconds, memory used: $(round(t.bytes / 1024^3, digits=2)) GiB"
            t.value
        end
    else
        (Y, return_code) = sim_run(sim, integrator; time_run)
    end
    return_code === :success && @info "The simulation has completed."
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
