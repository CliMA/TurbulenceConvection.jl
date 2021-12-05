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

const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))

include("initial_conditions.jl")
include(joinpath(tc_dir, "diagnostics", "compute_diagnostics.jl"))
include("parameter_set.jl")
include("Cases.jl")
import .Cases

struct Simulation1d
    io_nt::NamedTuple
    grid
    state
    GMV
    Case
    Turb
    diagnostics
    TS
    Stats
    param_set
    skip_io::Bool
    adapt_dt::Bool
    cfl_limit::Float64
    dt_min::Float64
end

# TODO: not quite sure, yet, where this all should live.

#####
##### Auxiliary fields
#####
# Face & Center
aux_vars_ref_state(FT) = (; ref_state = (ρ0 = FT(0), α0 = FT(0), p0 = FT(0)))

# Center only
cent_aux_vars_gm(FT) = (;
    tke = FT(0),
    Hvar = FT(0),
    QTvar = FT(0),
    HQTcov = FT(0),
    q_liq = FT(0),
    q_ice = FT(0),
    RH = FT(0),
    s = FT(0),
    T = FT(0),
    buoy = FT(0),
    cloud_fraction = FT(0),
    H_third_m = FT(0),
    W_third_m = FT(0),
    QT_third_m = FT(0),
    # From RadiationBase
    dTdt_rad = FT(0), # horizontal advection temperature tendency
    dqtdt_rad = FT(0), # horizontal advection moisture tendency
    # From ForcingBase
    subsidence = FT(0), #Large-scale subsidence
    dTdt = FT(0), #Large-scale temperature tendency
    dqtdt = FT(0), #Large-scale moisture tendency
    dTdt_hadv = FT(0), #Horizontal advection of temperature
    H_nudge = FT(0), #Reference H profile for relaxation tendency
    dTdt_fluc = FT(0), #Vertical turbulent advection of temperature
    dqtdt_hadv = FT(0), #Horizontal advection of moisture
    qt_nudge = FT(0), #Reference qt profile for relaxation tendency
    dqtdt_fluc = FT(0), #Vertical turbulent advection of moisture
    u_nudge = FT(0), #Reference u profile for relaxation tendency
    v_nudge = FT(0), #Reference v profile for relaxation tendency
    ug = FT(0), #Geostrophic u velocity
    vg = FT(0), #Geostrophic v velocity
    ∇θ_liq_ice_gm = FT(0),
    ∇q_tot_gm = FT(0),
)
cent_aux_vars_en_2m(FT) = (;
    dissipation = FT(0),
    shear = FT(0),
    entr_gain = FT(0),
    detr_loss = FT(0),
    press = FT(0),
    buoy = FT(0),
    interdomain = FT(0),
    rain_src = FT(0),
)
cent_aux_vars_up(FT) = (;
    q_liq = FT(0),
    q_ice = FT(0),
    T = FT(0),
    RH = FT(0),
    s = FT(0),
    buoy = FT(0),
    area = FT(0),
    q_tot = FT(0),
    θ_liq_ice = FT(0),
    θ_liq_ice_tendency_precip_formation = FT(0),
    qt_tendency_precip_formation = FT(0),
    entr_sc = FT(0),
    detr_sc = FT(0),
    frac_turb_entr = FT(0),
    entr_turb_dyn = FT(0),
    detr_turb_dyn = FT(0),
    asp_ratio = FT(0),
    horiz_K_eddy = FT(0),
)
cent_aux_vars_edmf(FT, n_up) = (;
    turbconv = (;
        bulk = (;
            area = FT(0),
            θ_liq_ice = FT(0),
            RH = FT(0),
            buoy = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            T = FT(0),
            cloud_fraction = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            qt_tendency_precip_formation = FT(0),
            Ri = FT(0),
            D_env = FT(0),
        ),
        up = ntuple(i -> cent_aux_vars_up(FT), n_up),
        en = (;
            w = FT(0),
            area = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            θ_liq_ice = FT(0),
            θ_virt = FT(0),
            θ_dry = FT(0),
            RH = FT(0),
            s = FT(0),
            T = FT(0),
            buoy = FT(0),
            cloud_fraction = FT(0),
            tke = FT(0),
            Hvar = FT(0),
            QTvar = FT(0),
            HQTcov = FT(0),
            qt_tendency_precip_formation = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            unsat = (; q_tot = FT(0), θ_dry = FT(0), θ_virt = FT(0)),
            sat = (; T = FT(0), q_vap = FT(0), q_tot = FT(0), θ_dry = FT(0), θ_liq_ice = FT(0)),
            Hvar_rain_dt = FT(0),
            QTvar_rain_dt = FT(0),
            HQTcov_rain_dt = FT(0),
        ),
        θ_liq_ice_tendency_precip_sinks = FT(0),
        qt_tendency_precip_sinks = FT(0),
        qr_tendency_evap = FT(0),
        qs_tendency_melt = FT(0),
        qs_tendency_dep_sub = FT(0),
        qr_tendency_advection = FT(0),
        qs_tendency_advection = FT(0),
        en_2m = (;
            tke = cent_aux_vars_en_2m(FT),
            Hvar = cent_aux_vars_en_2m(FT),
            QTvar = cent_aux_vars_en_2m(FT),
            HQTcov = cent_aux_vars_en_2m(FT),
        ),
        term_vel_rain = FT(0),
        term_vel_snow = FT(0),
        KM = FT(0),
        KH = FT(0),
        mixing_length = FT(0),
        m_entr_detr = FT(0),
        ∇m_entr_detr = FT(0),
        θ_virt = FT(0),
        ∇θ_virt = FT(0),
        ∇Ri_bulk = FT(0),
        massflux_tendency_h = FT(0),
        massflux_tendency_qt = FT(0),
        diffusive_tendency_h = FT(0),
        diffusive_tendency_qt = FT(0),
        prandtl_nvec = FT(0),
        # Added by Ignacio : Length scheme in use (mls), and smooth min effect (ml_ratio)
        # Variable Prandtl number initialized as neutral value.
        mls = FT(0),
        b_exch = FT(0),
        ml_ratio = FT(0),
        w_up_c = FT(0),
        w_en_c = FT(0),
        Shear² = FT(0),
        ∂θv∂z = FT(0),
        ∂qt∂z = FT(0),
        ∂θl∂z = FT(0),
        ∂θv∂z_unsat = FT(0),
        ∂qt∂z_sat = FT(0),
        ∂θl∂z_sat = FT(0),
        l_entdet = FT(0),
        ϕ_gm = FT(0), # temporary for grid-mean variables
        ϕ_gm_cov = FT(0), # temporary for grid-mean covariance variables
        ϕ_en_cov = FT(0), # temporary for environmental covariance variables
        ϕ_up_cubed = FT(0), # temporary for cubed updraft variables in grid mean 3rd moment functions
    ),
)
cent_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., cent_aux_vars_gm(FT)..., cent_aux_vars_edmf(FT, n_up)...)

# Face only
face_aux_vars_gm(FT) = (; massflux_s = FT(0), diffusive_flux_s = FT(0), total_flux_s = FT(0), f_rad = FT(0))
face_aux_vars_up(FT) = (;
    w = FT(0),
    nh_pressure = FT(0),
    nh_pressure_b = FT(0),
    nh_pressure_adv = FT(0),
    nh_pressure_drag = FT(0),
    massflux = FT(0),
)

face_aux_vars_edmf(FT, n_up) = (;
    turbconv = (;
        bulk = (; w = FT(0)),
        ρ_ae_KM = FT(0),
        ρ_ae_KH = FT(0),
        ρ_ae_K = FT(0),
        en = (; w = FT(0)),
        up = ntuple(i -> face_aux_vars_up(FT), n_up),
        massflux_h = FT(0),
        massflux_qt = FT(0),
        ϕ_temporary = FT(0),
        diffusive_flux_h = FT(0),
        diffusive_flux_qt = FT(0),
        diffusive_flux_u = FT(0),
        diffusive_flux_v = FT(0),
    ),
)
face_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., face_aux_vars_gm(FT)..., face_aux_vars_edmf(FT, n_up)...)

#####
##### Diagnostic fields
#####

# Center only
cent_diagnostic_vars_gm(FT) = ()
cent_diagnostic_vars_edmf(FT, n_up) = (;
    turbconv = (;
        asp_ratio = FT(0),
        entr_sc = FT(0),
        detr_sc = FT(0),
        massflux = FT(0),
        frac_turb_entr = FT(0),
        horiz_K_eddy = FT(0),
    ),
)
cent_diagnostic_vars(FT, n_up) = (; cent_diagnostic_vars_gm(FT)..., cent_diagnostic_vars_edmf(FT, n_up)...)

# Face only
face_diagnostic_vars_gm(FT) = ()
face_diagnostic_vars_edmf(FT, n_up) =
    (; turbconv = (; nh_pressure = FT(0), nh_pressure_adv = FT(0), nh_pressure_drag = FT(0), nh_pressure_b = FT(0)))
face_diagnostic_vars(FT, n_up) = (; face_diagnostic_vars_gm(FT)..., face_diagnostic_vars_edmf(FT, n_up)...)

#####
##### Prognostic fields
#####

# Center only
cent_prognostic_vars(FT, n_up) = (; cent_prognostic_vars_gm(FT)..., cent_prognostic_vars_edmf(FT, n_up)...)
cent_prognostic_vars_gm(FT) = (; u = FT(0), v = FT(0), θ_liq_ice = FT(0), q_tot = FT(0))
cent_prognostic_vars_up(FT) = (; ρarea = FT(0), ρaθ_liq_ice = FT(0), ρaq_tot = FT(0))
cent_prognostic_vars_en(FT) = (; ρatke = FT(0), ρaHvar = FT(0), ρaQTvar = FT(0), ρaHQTcov = FT(0))
cent_prognostic_vars_edmf(FT, n_up) = (;
    turbconv = (;
        en = cent_prognostic_vars_en(FT),
        up = ntuple(i -> cent_prognostic_vars_up(FT), n_up),
        pr = (; q_rai = FT(0), q_sno = FT(0)),
    ),
)
# cent_prognostic_vars_edmf(FT, n_up) = (;) # could also use this for empty model

# Face only
face_prognostic_vars(FT, n_up) = (; w = FT(0), face_prognostic_vars_edmf(FT, n_up)...)
face_prognostic_vars_up(FT) = (; ρaw = FT(0))
face_prognostic_vars_edmf(FT, n_up) = (; turbconv = (; up = ntuple(i -> face_prognostic_vars_up(FT), n_up)))
# face_prognostic_vars_edmf(FT, n_up) = (;) # could also use this for empty model

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

    GMV = TC.GridMeanVariables(namelist, grid, param_set)
    Sur = TC.SurfaceBase(Cases.get_surface_type(case); namelist, ref_params)
    Fo = TC.ForcingBase{Cases.get_forcing_type(case)}()
    Rad = TC.RadiationBase{Cases.get_radiation_type(case)}()

    Case = Cases.CasesBase(case, namelist, grid, param_set, Sur, Fo, Rad)
    Turb = TC.EDMF_PrognosticTKE(namelist, grid, param_set)
    TS = TC.TimeStepping(namelist)

    n_updrafts = Turb.n_updrafts

    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)

    cent_prog_fields() = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars(FT, n_updrafts))
    face_prog_fields() = TC.FieldFromNamedTuple(fspace, face_prognostic_vars(FT, n_updrafts))
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT, n_updrafts))
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT, n_updrafts))
    diagnostic_cent_fields = TC.FieldFromNamedTuple(cspace, cent_diagnostic_vars(FT, n_updrafts))
    diagnostic_face_fields = TC.FieldFromNamedTuple(fspace, face_diagnostic_vars(FT, n_updrafts))

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
        GMV,
        Case,
        Turb,
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
    Cases.initialize_profiles(sim.Case, sim.grid, sim.GMV, state)
    TC.satadjust(sim.GMV, sim.grid, sim.state)

    Cases.initialize_surface(sim.Case, sim.grid, state, sim.param_set)
    Cases.initialize_forcing(sim.Case, sim.grid, state, sim.GMV, sim.param_set)
    Cases.initialize_radiation(sim.Case, sim.grid, state, sim.GMV, sim.param_set)

    initialize_edmf(sim.Turb, sim.grid, state, sim.Case, sim.GMV, sim.TS)

    sim.skip_io && return nothing
    TC.initialize_io(sim.io_nt.ref_state, sim.Stats)
    TC.io(sim.io_nt.ref_state, sim.Stats, state) # since the reference prog is static

    TC.initialize_io(sim.io_nt.aux, sim.Stats)
    TC.initialize_io(sim.io_nt.diagnostics, sim.Stats)

    # TODO: depricate
    TC.initialize_io(sim.GMV, sim.Stats)
    TC.initialize_io(sim.Case, sim.Stats)
    TC.initialize_io(sim.Turb, sim.Stats)

    TC.open_files(sim.Stats)
    TC.write_simulation_time(sim.Stats, sim.TS.t)

    TC.io(sim.io_nt.aux, sim.Stats, state)
    TC.io(sim.io_nt.diagnostics, sim.Stats, sim.diagnostics)

    # TODO: depricate
    TC.io(sim.GMV, sim.grid, state, sim.Stats)
    TC.io(sim.Case, sim.grid, state, sim.Stats)
    TC.io(sim.Turb, sim.grid, state, sim.Stats, sim.TS, sim.param_set)
    TC.close_files(sim.Stats)

    return
end

function affect_io!(integrator)
    UnPack.@unpack edmf, aux, grid, io_nt, diagnostics, case, gm, TS, Stats, skip_io = integrator.p
    skip_io && return nothing

    state = TC.State(integrator.u, aux, integrator.du)

    param_set = TC.parameter_set(gm)
    # TODO: is this the best location to call diagnostics?
    compute_diagnostics!(edmf, gm, grid, state, diagnostics, case, TS)

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack
    # TurbulenceConvection.io(sim) # #removeVarsHack
    TC.write_simulation_time(Stats, TS.t) # #removeVarsHack

    TC.io(io_nt.aux, Stats, state)
    TC.io(io_nt.diagnostics, Stats, diagnostics)

    TC.io(gm, grid, state, Stats) # #removeVarsHack
    TC.io(case, grid, state, Stats) # #removeVarsHack
    TC.io(edmf, grid, state, Stats, TS, param_set) # #removeVarsHack

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end

function affect_filter!(integrator)
    UnPack.@unpack edmf, grid, gm, aux, case, TS = integrator.p
    state = TC.State(integrator.u, aux, integrator.du)
    TC.affect_filter!(edmf, grid, state, gm, case, TS)

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end

function adaptive_dt!(integrator)
    UnPack.@unpack edmf, TS, dt_min = integrator.p
    TS.dt = min(TS.dt_max, max(edmf.dt_max, dt_min))
    SciMLBase.set_proposed_dt!(integrator, TS.dt)
    ODE.u_modified!(integrator, false)
end

function dt_max!(integrator)
    UnPack.@unpack gm, grid, edmf, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, integrator.du)
    Δz = grid.Δz
    CFL_limit = TS.cfl_limit
    N_up = TC.n_updrafts(edmf)

    dt_max = TS.dt_max # initialize dt_max

    aux_tc = TC.center_aux_turbconv(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_en_f = TC.face_aux_environment(state)
    KM = aux_tc.KM
    KH = aux_tc.KH

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_face_indices(grid)
        TC.is_surface_face(grid, k) && continue
        @inbounds for i in 1:N_up
            dt_max = min(dt_max, CFL_limit * Δz / (abs(aux_up_f[i].w[k]) + eps(Float32)))
        end
        dt_max = min(dt_max, CFL_limit * Δz / (abs(aux_en_f.w[k]) + eps(Float32)))
    end
    @inbounds for k in TC.real_center_indices(grid)
        vel_max = max(term_vel_rain[k], term_vel_snow[k])
        # Check terminal rain/snow velocity CFL
        dt_max = min(dt_max, CFL_limit * Δz / (vel_max + eps(Float32)))
        # Check diffusion CFL (i.e., Fourier number)
        dt_max = min(dt_max, CFL_limit * Δz^2 / (max(KH[k], KM[k]) + eps(Float32)))
    end
    edmf.dt_max = dt_max

    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(integrator)
    UnPack.@unpack gm, grid, edmf, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, integrator.du)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = TS.cfl_limit

    aux_tc = TC.center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_center_indices(grid)
        # check stability criterion
        CFL_out_rain = Δt / Δz * term_vel_rain[k]
        CFL_out_snow = Δt / Δz * term_vel_snow[k]
        if TC.is_toa_center(grid, k)
            CFL_in_rain = 0.0
            CFL_in_snow = 0.0
        else
            CFL_in_rain = Δt / Δz * term_vel_rain[k + 1]
            CFL_in_snow = Δt / Δz * term_vel_snow[k + 1]
        end
        if max(CFL_in_rain, CFL_in_snow, CFL_out_rain, CFL_out_snow) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end
    end

    ODE.u_modified!(integrator, false)
end


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
        edmf = sim.Turb,
        grid = grid,
        gm = sim.GMV,
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

    function condition_io(u, t, integrator)
        UnPack.@unpack TS, Stats = integrator.p
        TS.dt_io += TS.dt
        io_flag = false
        if TS.dt_io > Stats.frequency
            TS.dt_io = 0
            io_flag = true
        end
        return io_flag || t ≈ 0 || t ≈ TS.t_max
    end
    condition_every_iter(u, t, integrator) = true

    callback_io = ODE.DiscreteCallback(condition_io, affect_io!; save_positions = (false, false))
    callback_cfl = ODE.DiscreteCallback(condition_every_iter, monitor_cfl!; save_positions = (false, false))
    callback_cfl = sim.Turb.Precip.precipitation_model == "clima_1m" ? (callback_cfl,) : ()
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
