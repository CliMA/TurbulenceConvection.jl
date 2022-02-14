import UnPack

import TurbulenceConvection
const TC = TurbulenceConvection

import Thermodynamics
const TD = Thermodynamics

import ClimaCore
const CC = ClimaCore

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

import CLIMAParameters
const CPP = CLIMAParameters.Planet

#####
##### Fields
#####

##### Auxiliary fields

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
    ∇u_gm = FT(0),
    ∇v_gm = FT(0),
)
cent_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., cent_aux_vars_gm(FT)..., TC.cent_aux_vars_edmf(FT, n_up)...)

# Face only
face_aux_vars_gm(FT) = (;
    massflux_s = FT(0),
    diffusive_flux_s = FT(0),
    total_flux_s = FT(0),
    f_rad = FT(0),
    sgs_flux_θ_liq_ice = FT(0),
    sgs_flux_q_tot = FT(0),
    sgs_flux_u = FT(0),
    sgs_flux_v = FT(0),
)
face_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., face_aux_vars_gm(FT)..., TC.face_aux_vars_edmf(FT, n_up)...)

##### Diagnostic fields

# Center only
cent_diagnostic_vars_gm(FT) = ()
cent_diagnostic_vars(FT, n_up) = (; cent_diagnostic_vars_gm(FT)..., TC.cent_diagnostic_vars_edmf(FT, n_up)...)

# Face only
face_diagnostic_vars_gm(FT) = ()
face_diagnostic_vars(FT, n_up) = (; face_diagnostic_vars_gm(FT)..., TC.face_diagnostic_vars_edmf(FT, n_up)...)

##### Prognostic fields

# Center only
cent_prognostic_vars(FT, n_up) = (; cent_prognostic_vars_gm(FT)..., TC.cent_prognostic_vars_edmf(FT, n_up)...)
cent_prognostic_vars_gm(FT) = (; u = FT(0), v = FT(0), θ_liq_ice = FT(0), q_tot = FT(0))

# Face only
face_prognostic_vars(FT, n_up) = (; w = FT(0), TC.face_prognostic_vars_edmf(FT, n_up)...)
# TC.face_prognostic_vars_edmf(FT, n_up) = (;) # could also use this for empty model


#####
##### Methods
#####

####
#### Reference state
####

"""
    compute_ref_state!(
        state,
        grid::Grid,
        param_set::PS;
        Pg::FT,
        Tg::FT,
        qtg::FT,
    ) where {PS, FT}

TODO: add better docs once the API converges

The reference profiles, given
 - `grid` the grid
 - `param_set` the parameter set
"""
function compute_ref_state!(state, grid::TC.Grid, param_set::PS; Pg::FT, Tg::FT, qtg::FT) where {PS, FT}

    p0_c = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
    α0_c = TC.center_ref_state(state).α0
    p0_f = TC.face_ref_state(state).p0
    ρ0_f = TC.face_ref_state(state).ρ0
    α0_f = TC.face_ref_state(state).α0

    q_pt_g = TD.PhasePartition(qtg)
    ts_g = TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
    θ_liq_ice_g = TD.liquid_ice_pottemp(ts_g)

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    logp = log(Pg)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure
    function rhs(logp, u, z)
        p_ = exp(logp)
        ts = TD.PhaseEquil_pθq(param_set, p_, θ_liq_ice_g, qtg)
        R_m = TD.gas_constant_air(ts)
        T = TD.air_temperature(ts)
        return -FT(CPP.grav(param_set)) / (T * R_m)
    end

    # Perform the integration
    z_span = (grid.zmin, grid.zmax)
    @show z_span
    prob = ODE.ODEProblem(rhs, logp, z_span)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)

    parent(p0_f) .= sol.(vec(grid.zf))
    parent(p0_c) .= sol.(vec(grid.zc))

    p0_f .= exp.(p0_f)
    p0_c .= exp.(p0_c)

    # Compute reference state thermodynamic profiles
    @inbounds for k in TC.real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], θ_liq_ice_g, qtg)
        α0_c[k] = TD.specific_volume(ts)
    end

    @inbounds for k in TC.real_face_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_f[k], θ_liq_ice_g, qtg)
        α0_f[k] = TD.specific_volume(ts)
    end

    ρ0_f .= 1 ./ α0_f
    ρ0_c .= 1 ./ α0_c
    return nothing
end

function satadjust(gm::TC.GridMeanVariables, grid, state)
    p0_c = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    param_set = TC.parameter_set(gm)
    @inbounds for k in TC.real_center_indices(grid)
        θ_liq_ice = prog_gm.θ_liq_ice[k]
        q_tot = prog_gm.q_tot[k]
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], θ_liq_ice, q_tot)
        aux_gm.q_liq[k] = TD.liquid_specific_humidity(ts)
        aux_gm.q_ice[k] = TD.ice_specific_humidity(ts)
        aux_gm.T[k] = TD.air_temperature(ts)
        ρ = TD.air_density(ts)
        aux_gm.buoy[k] = TC.buoyancy_c(param_set, ρ0_c[k], ρ)
        aux_gm.RH[k] = TD.relative_humidity(ts)
    end
    return
end


# Compute the sum of tendencies for the scheme
function ∑tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    UnPack.@unpack turb_conv, grid, gm, case, aux, TS = params

    state = TC.State(prog, aux, tendencies)

    Δt = TS.dt
    param_set = TC.parameter_set(gm)
    surf = get_surface(case.surf_params, grid, state, gm, t, param_set)
    force = case.Fo
    radiation = case.Rad
    en_thermo = turb_conv.en_thermo
    precip_model = turb_conv.precip_model

    TC.affect_filter!(turb_conv, grid, state, gm, surf, case.casename, t)

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.
    Cases.update_forcing(case, grid, state, gm, t, param_set)
    Cases.update_radiation(case.Rad, grid, state, gm, param_set)

    TC.update_aux!(turb_conv, gm, grid, state, case, surf, param_set, t, Δt)

    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0

    # causes division error in dry bubble first time step
    TC.compute_precipitation_formation_tendencies(grid, state, turb_conv, precip_model, Δt, param_set)

    TC.microphysics(en_thermo, grid, state, precip_model, Δt, param_set)
    TC.compute_precipitation_sink_tendencies(precip_model, grid, state, gm, Δt)
    TC.compute_precipitation_advection_tendencies(precip_model, turb_conv, grid, state, gm)

    # compute tendencies
    compute_gm_tendencies!(turb_conv, grid, state, surf, radiation, force, gm)
    compute_turbconv_tendencies!(turb_conv, param_set, grid, state, gm, surf)

    return nothing
end

function compute_turbconv_tendencies!(turb_conv::TC.EDMFModel, param_set, grid, state, gm, surf)
    TC.compute_up_tendencies!(turb_conv, grid, state, gm, surf)
    TC.compute_en_tendencies!(turb_conv, grid, state, param_set, Val(:tke), Val(:ρatke))
    TC.compute_en_tendencies!(turb_conv, grid, state, param_set, Val(:Hvar), Val(:ρaHvar))
    TC.compute_en_tendencies!(turb_conv, grid, state, param_set, Val(:QTvar), Val(:ρaQTvar))
    TC.compute_en_tendencies!(turb_conv, grid, state, param_set, Val(:HQTcov), Val(:ρaHQTcov))
    return nothing
end

function compute_turbconv_tendencies!(turb_conv::TC.DiffusivityModel, param_set, grid, state, gm, surf)
    return nothing
end


function compute_sgs_tendencies!(turb_conv::TC.EDMFModel, param_set, grid, state, gm, surf)
    TC.compute_sgs_flux!(turb_conv, grid, state, surf, gm)

    tendencies_gm = TC.center_tendencies_grid_mean(state)
    kf_surf = TC.kf_surface(grid)
    aux_gm_f = TC.face_aux_grid_mean(state)
    α0_c = TC.center_ref_state(state).α0
    wvec = CC.Geometry.WVector

    sgs_flux_θ_liq_ice = aux_gm_f.sgs_flux_θ_liq_ice
    sgs_flux_q_tot = aux_gm_f.sgs_flux_q_tot
    sgs_flux_u = aux_gm_f.sgs_flux_u
    sgs_flux_v = aux_gm_f.sgs_flux_v
    # apply surface BC as SGS flux at lowest level
    sgs_flux_θ_liq_ice[kf_surf] = surf.ρθ_liq_ice_flux
    sgs_flux_q_tot[kf_surf] = surf.ρq_tot_flux
    sgs_flux_u[kf_surf] = surf.ρu_flux
    sgs_flux_v[kf_surf] = surf.ρv_flux

    tends_θ_liq_ice = tendencies_gm.θ_liq_ice
    tends_q_tot = tendencies_gm.q_tot
    tends_u = tendencies_gm.u
    tends_v = tendencies_gm.v

    ∇θ_liq_ice_sgs = CCO.DivergenceF2C()
    ∇q_tot_sgs = CCO.DivergenceF2C()
    ∇u_sgs = CCO.DivergenceF2C()
    ∇v_sgs = CCO.DivergenceF2C()

    @. tends_θ_liq_ice += -α0_c * ∇θ_liq_ice_sgs(wvec(sgs_flux_θ_liq_ice))
    @. tends_q_tot += -α0_c * ∇q_tot_sgs(wvec(sgs_flux_q_tot))
    @. tends_u += -α0_c * ∇u_sgs(wvec(sgs_flux_u))
    @. tends_v += -α0_c * ∇v_sgs(wvec(sgs_flux_v))
    return nothing
end

function compute_sgs_tendencies!(turb_conv::TC.DiffusivityModel, param_set, grid, state, gm, surf)
    FT = eltype(grid)
    zf = grid.zf
    kf_surf = TC.kf_surface(grid)
    ρ0_f = TC.face_ref_state(state).ρ0
    wvec = CC.Geometry.WVector
    prog_gm = TC.center_prog_grid_mean(state)
    tendencies_gm = TC.center_tendencies_grid_mean(state)
    # ∇θ_liq_ice_gm = TC.center_aux_grid_mean(state).∇θ_liq_ice_gm
    # ∇q_tot_gm = TC.center_aux_grid_mean(state).∇q_tot_gm
    # ∇u_gm = TC.center_aux_grid_mean(state).∇u_gm
    # ∇v_gm = TC.center_aux_grid_mean(state).∇v_gm

    ν = turb_conv.diff_coeff
    grad_θ = CCO.GradientC2F(; bottom = CCO.SetGradient(wvec(surf.ρθ_liq_ice_flux/ρ0_f[kf_surf])), top = CCO.SetGradient(wvec(FT(0))))
    grad_q = CCO.GradientC2F(; bottom = CCO.SetGradient(wvec(surf.ρq_tot_flux/ρ0_f[kf_surf])), top = CCO.SetGradient(wvec(FT(0))))
    grad_u = CCO.GradientC2F(; bottom = CCO.SetGradient(wvec(surf.ρu_flux/ρ0_f[kf_surf])), top = CCO.SetGradient(wvec(FT(0))))
    grad_v = CCO.GradientC2F(; bottom = CCO.SetGradient(wvec(surf.ρv_flux/ρ0_f[kf_surf])), top = CCO.SetGradient(wvec(FT(0))))
    ∇c = CCO.DivergenceF2C()

    @. tendencies_gm.θ_liq_ice += ν*sqrt(∇u_gm^2+∇v_gm^2)*∇c(grad_θ(prog_gm.θ_liq_ice))
    @. tendencies_gm.q_tot += ν*sqrt(∇u_gm^2+∇v_gm^2)*∇c(grad_q(prog_gm.q_tot))
    @. tendencies_gm.u += ν*sqrt(∇u_gm^2+∇v_gm^2)*∇c(grad_θ(prog_gm.θ_liq_ice))
    @. tendencies_gm.v += ν*sqrt(∇u_gm^2+∇v_gm^2)*∇c(grad_q(prog_gm.q_tot))

    return nothing
end


function compute_gm_tendencies!(
    turb_conv::ATCM,
    grid::TC.Grid,
    state::TC.State,
    surf::TC.SurfaceBase,
    radiation::TC.RadiationBase,
    force::TC.ForcingBase,
    gm::TC.GridMeanVariables,
) where {ATCM}
    tendencies_gm = TC.center_tendencies_grid_mean(state)
    kc_toa = TC.kc_top_of_atmos(grid)
    FT = eltype(grid)
    param_set = TC.parameter_set(gm)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    ∇θ_liq_ice_gm = TC.center_aux_grid_mean(state).∇θ_liq_ice_gm
    ∇q_tot_gm = TC.center_aux_grid_mean(state).∇q_tot_gm
    aux_en = TC.center_aux_environment(state)
    aux_bulk = TC.center_aux_bulk(state)
    p0_c = TC.center_ref_state(state).p0
    aux_tc = TC.center_aux_turbconv(state)

    θ_liq_ice_gm_toa = prog_gm.θ_liq_ice[kc_toa]
    q_tot_gm_toa = prog_gm.q_tot[kc_toa]
    RBθ = CCO.RightBiasedC2F(; top = CCO.SetValue(θ_liq_ice_gm_toa))
    RBq = CCO.RightBiasedC2F(; top = CCO.SetValue(q_tot_gm_toa))
    wvec = CC.Geometry.WVector
    ∇c = CCO.DivergenceF2C()
    @. ∇θ_liq_ice_gm = ∇c(wvec(RBθ(prog_gm.θ_liq_ice)))
    @. ∇q_tot_gm = ∇c(wvec(RBq(prog_gm.q_tot)))

    @inbounds for k in TC.real_center_indices(grid)
        # Apply large-scale horizontal advection tendencies
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(ts)

        if force.apply_coriolis
            tendencies_gm.u[k] -= force.coriolis_param * (aux_gm.vg[k] - prog_gm.v[k])
            tendencies_gm.v[k] += force.coriolis_param * (aux_gm.ug[k] - prog_gm.u[k])
        end
        if TC.rad_type(radiation) <: Union{TC.RadiationDYCOMS_RF01, TC.RadiationLES}
            tendencies_gm.θ_liq_ice[k] += aux_gm.dTdt_rad[k] / Π
        end

        if TC.force_type(force) <: TC.ForcingDYCOMS_RF01
            tendencies_gm.q_tot[k] += aux_gm.dqtdt[k]
            # Apply large-scale subsidence tendencies
            tendencies_gm.θ_liq_ice[k] -= ∇θ_liq_ice_gm[k] * aux_gm.subsidence[k]
            tendencies_gm.q_tot[k] -= ∇q_tot_gm[k] * aux_gm.subsidence[k]
        end

        if TC.force_type(force) <: TC.ForcingStandard
            if force.apply_subsidence
                tendencies_gm.θ_liq_ice[k] -= ∇θ_liq_ice_gm[k] * aux_gm.subsidence[k]
                tendencies_gm.q_tot[k] -= ∇q_tot_gm[k] * aux_gm.subsidence[k]
            end
            tendencies_gm.θ_liq_ice[k] += aux_gm.dTdt[k] / Π
            tendencies_gm.q_tot[k] += aux_gm.dqtdt[k]
        end

        if TC.force_type(force) <: TC.ForcingLES
            H_horz_adv = aux_gm.dTdt_hadv[k] / Π
            H_fluc = aux_gm.dTdt_fluc[k] / Π

            gm_U_nudge_k = (aux_gm.u_nudge[k] - prog_gm.u[k]) / force.nudge_tau
            gm_V_nudge_k = (aux_gm.v_nudge[k] - prog_gm.v[k]) / force.nudge_tau

            Γᵣ = TC.compute_les_Γᵣ(grid.zc[k])
            if Γᵣ != 0
                tau_k = 1 / Γᵣ
                gm_H_nudge_k = (aux_gm.H_nudge[k] - prog_gm.θ_liq_ice[k]) / tau_k
                gm_q_tot_nudge_k = (aux_gm.qt_nudge[k] - prog_gm.q_tot[k]) / tau_k
            else
                gm_H_nudge_k = 0.0
                gm_q_tot_nudge_k = 0.0
            end

            if force.apply_subsidence
                # Apply large-scale subsidence tendencies
                gm_H_subsidence_k = -∇θ_liq_ice_gm[k] * aux_gm.subsidence[k]
                gm_QT_subsidence_k = -∇q_tot_gm[k] * aux_gm.subsidence[k]
            else
                gm_H_subsidence_k = 0.0
                gm_QT_subsidence_k = 0.0
            end

            tendencies_gm.θ_liq_ice[k] += H_horz_adv + gm_H_nudge_k + H_fluc + gm_H_subsidence_k
            tendencies_gm.q_tot[k] +=
                aux_gm.dqtdt_hadv[k] + gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k] + gm_QT_subsidence_k

            tendencies_gm.u[k] += gm_U_nudge_k
            tendencies_gm.v[k] += gm_V_nudge_k
        end
        tendencies_gm.q_tot[k] +=
            aux_bulk.qt_tendency_precip_formation[k] +
            aux_en.qt_tendency_precip_formation[k] +
            aux_tc.qt_tendency_precip_sinks[k]
        tendencies_gm.θ_liq_ice[k] +=
            aux_bulk.θ_liq_ice_tendency_precip_formation[k] +
            aux_en.θ_liq_ice_tendency_precip_formation[k] +
            aux_tc.θ_liq_ice_tendency_precip_sinks[k]
    end
    compute_sgs_tendencies!(turb_conv, param_set, grid, state, gm, surf)
    return nothing
end
