import UnPack

import TurbulenceConvection
const TC = TurbulenceConvection
const TCP = TC.TurbulenceConvectionParameters

import Thermodynamics
const TD = Thermodynamics

import ClimaCore
const CC = ClimaCore
const CCG = CC.Geometry

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

import CLIMAParameters
const APS = CLIMAParameters.AbstractEarthParameterSet

#####
##### Fields
#####

##### Auxiliary fields

# Center only
cent_aux_vars_gm_moisture(FT, ::TC.NonEquilibriumMoisture) = (;
    ∇q_liq_gm = FT(0),
    ∇q_ice_gm = FT(0),
    dqldt_rad = FT(0),
    dqidt_rad = FT(0),
    dqldt = FT(0),
    dqidt = FT(0),
    dqldt_hadv = FT(0),
    dqidt_hadv = FT(0),
    ql_nudge = FT(0),
    qi_nudge = FT(0),
    dqldt_fluc = FT(0),
    dqidt_fluc = FT(0),
)
cent_aux_vars_gm_moisture(FT, ::TC.EquilibriumMoisture) = NamedTuple()
cent_aux_vars_gm(FT, edmf) = (;
    ts = TC.thermo_state(FT, edmf.moisture_model),
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
    dTdt_hadv = FT(0), #Horizontal advection of temperature
    dqtdt_hadv = FT(0), #Horizontal advection of moisture
    T_nudge = FT(0), #Reference T profile for relaxation tendency
    qt_nudge = FT(0), #Reference qt profile for relaxation tendency
    dTdt_fluc = FT(0), #Vertical turbulent advection of temperature
    dqtdt_fluc = FT(0), #Vertical turbulent advection of moisture
    u_nudge = FT(0), #Reference u profile for relaxation tendency
    v_nudge = FT(0), #Reference v profile for relaxation tendency
    ug = FT(0), #Geostrophic u velocity
    vg = FT(0), #Geostrophic v velocity
    ∇MSE_gm = FT(0),
    ∇q_tot_gm = FT(0),
    cent_aux_vars_gm_moisture(FT, edmf.moisture_model)...,
    θ_virt = FT(0),
    Ri = FT(0),
    θ_liq_ice = FT(0),
    q_tot = FT(0),
    p = FT(0),
    e_kin = FT(0),
    h_tot = FT(0),
)
cent_aux_vars(FT, edmf) = (; cent_aux_vars_gm(FT, edmf)..., TC.cent_aux_vars_edmf(FT, edmf)...)

# Face only
face_aux_vars_gm_moisture(FT, ::TC.NonEquilibriumMoisture) = (; sgs_flux_q_liq = FT(0), sgs_flux_q_ice = FT(0))
face_aux_vars_gm_moisture(FT, ::TC.EquilibriumMoisture) = NamedTuple()
face_aux_vars_gm(FT, edmf) = (;
    massflux_s = FT(0),
    diffusive_flux_s = FT(0),
    total_flux_s = FT(0),
    f_rad = FT(0),
    sgs_flux_h_tot = FT(0),
    sgs_flux_q_tot = FT(0),
    face_aux_vars_gm_moisture(FT, edmf.moisture_model)...,
    sgs_flux_u = FT(0),
    sgs_flux_v = FT(0),
    p = FT(0),
    ρ = FT(0),
)
face_aux_vars(FT, edmf) = (; face_aux_vars_gm(FT, edmf)..., TC.face_aux_vars_edmf(FT, edmf)...)

##### Diagnostic fields

# Center only
cent_diagnostic_vars_gm(FT) = NamedTuple()
cent_diagnostic_vars(FT, edmf) = (; cent_diagnostic_vars_gm(FT)..., TC.cent_diagnostic_vars_edmf(FT, edmf)...)

# Face only
face_diagnostic_vars_gm(FT) = NamedTuple()
face_diagnostic_vars(FT, edmf) = (; face_diagnostic_vars_gm(FT)..., TC.face_diagnostic_vars_edmf(FT, edmf)...)

# Single value per column diagnostic variables
single_value_per_col_diagnostic_vars_gm(FT) = (;
    Tsurface = FT(0),
    shf = FT(0),
    lhf = FT(0),
    ustar = FT(0),
    wstar = FT(0),
    lwp_mean = FT(0),
    iwp_mean = FT(0),
    rwp_mean = FT(0),
    swp_mean = FT(0),
    cutoff_precipitation_rate = FT(0),
    cloud_base_mean = FT(0),
    cloud_top_mean = FT(0),
    cloud_cover_mean = FT(0),
)
single_value_per_col_diagnostic_vars(FT, edmf) =
    (; single_value_per_col_diagnostic_vars_gm(FT)..., TC.single_value_per_col_diagnostic_vars_edmf(FT, edmf)...)

##### Prognostic fields

# Center only
cent_prognostic_vars(::Type{FT}, local_geometry, edmf) where {FT} =
    (; cent_prognostic_vars_gm(FT, local_geometry, edmf)..., TC.cent_prognostic_vars_edmf(FT, edmf)...)
cent_prognostic_vars_gm_moisture(::Type{FT}, ::TC.NonEquilibriumMoisture) where {FT} = (; q_liq = FT(0), q_ice = FT(0))
cent_prognostic_vars_gm_moisture(::Type{FT}, ::TC.EquilibriumMoisture) where {FT} = NamedTuple()
cent_prognostic_vars_gm(::Type{FT}, local_geometry, edmf) where {FT} = (;
    ρ = FT(0),
    uₕ = CCG.Covariant12Vector(CCG.UVVector(FT(0), FT(0)), local_geometry),
    ρe_tot = FT(0),
    ρq_tot = FT(0),
    cent_prognostic_vars_gm_moisture(FT, edmf.moisture_model)...,
)

# Face only
face_prognostic_vars(::Type{FT}, local_geometry, edmf) where {FT} =
    (; w = CCG.Covariant3Vector(FT(0)), TC.face_prognostic_vars_edmf(FT, local_geometry, edmf)...)

# TC.face_prognostic_vars_edmf(FT, edmf) = (;) # could also use this for empty model


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
        ts_g,
    ) where {PS}

TODO: add better docs once the API converges

The reference profiles, given
 - `grid` the grid
 - `param_set` the parameter set
 - `ts_g` the surface reference state (a thermodynamic state)
"""
function compute_ref_state!(state, grid::TC.Grid, param_set::PS; ts_g) where {PS}
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    p_f = aux_gm_f.p
    ρ_f = aux_gm_f.ρ
    compute_ref_state!(p_c, ρ_c, p_f, ρ_f, grid, param_set; ts_g)
end

function compute_ref_state!(
    p_c::CC.Fields.Field,
    ρ_c::CC.Fields.Field,
    p_f::CC.Fields.Field,
    ρ_f::CC.Fields.Field,
    grid::TC.Grid,
    param_set::PS;
    ts_g,
) where {PS}
    FT = eltype(grid)
    qtg = TD.total_specific_humidity(param_set, ts_g)
    θ_liq_ice_g = TD.liquid_ice_pottemp(param_set, ts_g)
    Pg = TD.air_pressure(param_set, ts_g)

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    logp = log(Pg)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure
    function rhs(logp, u, z)
        p_ = exp(logp)
        ts = TD.PhaseEquil_pθq(param_set, p_, θ_liq_ice_g, qtg)
        R_m = TD.gas_constant_air(param_set, ts)
        T = TD.air_temperature(param_set, ts)
        return -FT(TCP.grav(param_set)) / (T * R_m)
    end

    # Perform the integration
    z_span = (grid.zmin, grid.zmax)
    prob = ODE.ODEProblem(rhs, logp, z_span)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    parent(p_f) .= sol.(vec(grid.zf.z))
    parent(p_c) .= sol.(vec(grid.zc.z))

    p_f .= exp.(p_f)
    p_c .= exp.(p_c)

    # Compute reference state thermodynamic profiles
    @inbounds for k in TC.real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p_c[k], θ_liq_ice_g, qtg)
        ρ_c[k] = TD.air_density(param_set, ts)
    end

    @inbounds for k in TC.real_face_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p_f[k], θ_liq_ice_g, qtg)
        ρ_f[k] = TD.air_density(param_set, ts)
    end
    return nothing
end


function set_thermo_state_peq!(state, grid, moisture_model, param_set)
    Ic = CCO.InterpolateF2C()
    FT = eltype(grid)
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm_u = TC.grid_mean_u(state)
    prog_gm_v = TC.grid_mean_v(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    w_c = copy(prog_gm.ρe_tot)
    @. w_c = Ic(FT(0) + prog_gm_f.w)
    @inbounds for k in TC.real_center_indices(grid)
        thermo_args = if moisture_model isa TC.EquilibriumMoisture
            ()
        elseif moisture_model isa TC.NonEquilibriumMoisture
            (prog_gm.q_liq[k], prog_gm.q_ice[k])
        else
            error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
        end
        aux_gm.e_kin[k] = TC.kinetic_energy(prog_gm_u[k], prog_gm_v[k], w_c[k])
        e_pot = TC.geopotential(param_set, grid.zc.z[k])
        e_int = prog_gm.ρe_tot[k] / ρ_c[k] - aux_gm.e_kin[k] - e_pot
        ts_gm[k] = TC.thermo_state_peq(param_set, p_c[k], e_int, aux_gm.q_tot[k], thermo_args...)
        aux_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp(param_set, ts_gm[k])
        aux_gm.q_tot[k] = prog_gm.ρq_tot[k] / ρ_c[k]
        aux_gm.h_tot[k] = TC.total_enthalpy(param_set, prog_gm.ρe_tot[k] / ρ_c[k], ts_gm[k])
    end
    return nothing
end

function set_thermo_state_pθq!(state, grid, moisture_model, param_set)
    Ic = CCO.InterpolateF2C()
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    p_c = aux_gm.p
    @inbounds for k in TC.real_center_indices(grid)
        thermo_args = if moisture_model isa TC.EquilibriumMoisture
            ()
        elseif moisture_model isa TC.NonEquilibriumMoisture
            (prog_gm.q_liq[k], prog_gm.q_ice[k])
        else
            error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
        end
        ts_gm[k] = TC.thermo_state_pθq(param_set, p_c[k], aux_gm.θ_liq_ice[k], aux_gm.q_tot[k], thermo_args...)
    end
    return nothing
end

function set_grid_mean_from_thermo_state!(param_set, state, grid)
    Ic = CCO.InterpolateF2C()
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm_u = TC.grid_mean_u(state)
    prog_gm_v = TC.grid_mean_v(state)
    ρ_c = prog_gm.ρ
    FT = eltype(grid)
    w_c = copy(prog_gm.ρe_tot)
    @. w_c = Ic(FT(0) + prog_gm_f.w)
    @inbounds for k in TC.real_center_indices(grid)
        e_kin = TC.kinetic_energy(prog_gm_u[k], prog_gm_v[k], w_c[k])
        e_pot = TC.geopotential(param_set, grid.zc.z[k])
        prog_gm.ρe_tot[k] = ρ_c[k] * TD.total_energy(param_set, ts_gm[k], e_kin, e_pot)
        prog_gm.ρq_tot[k] = ρ_c[k] * aux_gm.q_tot[k]
    end
    return nothing
end

function assign_thermo_aux!(state, grid, moisture_model, param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    ρ_c = prog_gm.ρ
    @inbounds for k in TC.real_center_indices(grid)
        ts = ts_gm[k]
        aux_gm.q_liq[k] = TD.liquid_specific_humidity(param_set, ts)
        aux_gm.q_ice[k] = TD.ice_specific_humidity(param_set, ts)
        aux_gm.T[k] = TD.air_temperature(param_set, ts)
        ρ = TD.air_density(param_set, ts)
        aux_gm.buoy[k] = TC.buoyancy_c(param_set, ρ_c[k], ρ)
        aux_gm.RH[k] = TD.relative_humidity(param_set, ts)
    end
    return
end

function ∑stoch_tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    UnPack.@unpack edmf, param_set, case, aux = params

    # set all tendencies to zero
    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0

    for inds in TC.iterate_columns(prog.cent)

        state = TC.column_state(prog, aux, tendencies, inds...)
        grid = TC.Grid(state)
        surf = get_surface(case.surf_params, grid, state, t, param_set)

        # compute updraft stochastic tendencies
        TC.compute_up_stoch_tendencies!(edmf, grid, state, param_set, surf)
    end
end

# Compute the sum of tendencies for the scheme
function ∑tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    UnPack.@unpack edmf, precip_model, param_set, case, aux, TS = params

    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0

    for inds in TC.iterate_columns(prog.cent)
        state = TC.column_state(prog, aux, tendencies, inds...)
        grid = TC.Grid(state)

        set_thermo_state_peq!(state, grid, edmf.moisture_model, param_set)

        aux_gm = TC.center_aux_grid_mean(state)

        @. aux_gm.θ_virt = TD.virtual_pottemp(param_set, aux_gm.ts)

        Δt = TS.dt
        surf = get_surface(case.surf_params, grid, state, t, param_set)
        force = case.Fo
        radiation = case.Rad

        TC.affect_filter!(edmf, grid, state, param_set, surf, case.casename, t)

        # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
        # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
        # treat them as auxiliary variables for now, until we disentangle the tendency computations.
        Cases.update_forcing(case, grid, state, t, param_set)
        Cases.update_radiation(case.Rad, grid, state, t, param_set)

        TC.update_aux!(edmf, grid, state, surf, param_set, t, Δt)

        en_thermo = edmf.en_thermo

        # compute tendencies
        # causes division error in dry bubble first time step
        TC.compute_precipitation_formation_tendencies(grid, state, edmf, precip_model, Δt, param_set)
        TC.microphysics(en_thermo, grid, state, edmf, precip_model, Δt, param_set)
        TC.compute_precipitation_sink_tendencies(precip_model, edmf, grid, state, param_set, Δt)
        TC.compute_precipitation_advection_tendencies(precip_model, edmf, grid, state, param_set)

        TC.compute_turbconv_tendencies!(edmf, grid, state, param_set, surf, Δt)
        compute_gm_tendencies!(edmf, grid, state, surf, radiation, force, param_set)
    end

    return nothing
end

function compute_gm_tendencies!(
    edmf::TC.EDMFModel,
    grid::TC.Grid,
    state::TC.State,
    surf::TC.SurfaceBase,
    radiation::TC.RadiationBase,
    force::TC.ForcingBase,
    param_set::APS,
)
    Ic = CCO.InterpolateF2C()
    R_d = TCP.R_d(param_set)
    T_0 = TCP.T_0(param_set)
    Lv_0 = TCP.LH_v0(param_set)
    tendencies_gm = TC.center_tendencies_grid_mean(state)
    kc_toa = TC.kc_top_of_atmos(grid)
    kf_surf = TC.kf_surface(grid)
    FT = eltype(grid)
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    ∇MSE_gm = TC.center_aux_grid_mean(state).∇MSE_gm
    ∇q_tot_gm = TC.center_aux_grid_mean(state).∇q_tot_gm
    aux_en = TC.center_aux_environment(state)
    aux_en_f = TC.face_aux_environment(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_bulk = TC.center_aux_bulk(state)
    prog_gm_u = TC.grid_mean_u(state)
    prog_gm_v = TC.grid_mean_v(state)
    ρ_f = aux_gm_f.ρ
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_tc = TC.center_aux_turbconv(state)
    ts_gm = TC.center_aux_grid_mean(state).ts

    MSE_gm_toa = aux_gm.h_tot[kc_toa] - aux_gm.e_kin[kc_toa]
    q_tot_gm_toa = prog_gm.ρq_tot[kc_toa] / ρ_c[kc_toa]
    RBe = CCO.RightBiasedC2F(; top = CCO.SetValue(MSE_gm_toa))
    RBq = CCO.RightBiasedC2F(; top = CCO.SetValue(q_tot_gm_toa))
    wvec = CC.Geometry.WVector
    ∇c = CCO.DivergenceF2C()
    @. ∇MSE_gm = ∇c(wvec(RBe(aux_gm.h_tot - aux_gm.e_kin)))
    @. ∇q_tot_gm = ∇c(wvec(RBq(prog_gm.ρq_tot / ρ_c)))

    if edmf.moisture_model isa TC.NonEquilibriumMoisture
        ∇q_liq_gm = TC.center_aux_grid_mean(state).∇q_liq_gm
        ∇q_ice_gm = TC.center_aux_grid_mean(state).∇q_ice_gm
        q_liq_gm_toa = prog_gm.q_liq[kc_toa]
        q_ice_gm_toa = prog_gm.q_ice[kc_toa]
        RBq_liq = CCO.RightBiasedC2F(; top = CCO.SetValue(q_liq_gm_toa))
        RBq_ice = CCO.RightBiasedC2F(; top = CCO.SetValue(q_ice_gm_toa))
        @. ∇q_liq_gm = ∇c(wvec(RBq(prog_gm.q_liq)))
        @. ∇q_ice_gm = ∇c(wvec(RBq(prog_gm.q_ice)))
    end

    # Apply forcing and radiation
    prog_gm_u = TC.grid_mean_u(state)
    prog_gm_v = TC.grid_mean_v(state)
    tendencies_gm_u = TC.tendencies_grid_mean_u(state)
    tendencies_gm_v = TC.tendencies_grid_mean_v(state)
    @inbounds for k in TC.real_center_indices(grid)
        cp_m = TD.cp_m(param_set, ts_gm[k])
        cp_v = TCP.cp_v(param_set)
        cv_v = TCP.cv_v(param_set)
        R_v = TCP.R_v(param_set)
        cv_m = TD.cv_m(param_set, ts_gm[k])
        h_v = cv_v * (aux_gm.T[k] - T_0) + Lv_0 - R_v * T_0

        # Coriolis
        if force.apply_coriolis
            tendencies_gm_u[k] -= force.coriolis_param * (aux_gm.vg[k] - prog_gm_v[k])
            tendencies_gm_v[k] += force.coriolis_param * (aux_gm.ug[k] - prog_gm_u[k])
        end
        # LS Subsidence
        tendencies_gm.ρe_tot[k] -= ρ_c[k] * aux_gm.subsidence[k] * ∇MSE_gm[k]
        tendencies_gm.ρq_tot[k] -= ρ_c[k] * aux_gm.subsidence[k] * ∇q_tot_gm[k]
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tendencies_gm.q_liq[k] -= ∇q_liq_gm[k] * aux_gm.subsidence[k]
            tendencies_gm.q_ice[k] -= ∇q_ice_gm[k] * aux_gm.subsidence[k]
        end
        # Radiation
        if TC.rad_type(radiation) <: Union{TC.RadiationDYCOMS_RF01, TC.RadiationLES, TC.RadiationTRMM_LBA}
            tendencies_gm.ρe_tot[k] += ρ_c[k] * cv_m * aux_gm.dTdt_rad[k]
        end
        # LS advection
        tendencies_gm.ρq_tot[k] += ρ_c[k] * aux_gm.dqtdt_hadv[k]
        if !(TC.force_type(force) <: TC.ForcingDYCOMS_RF01)
            tendencies_gm.ρe_tot[k] += ρ_c[k] * (cp_m * aux_gm.dTdt_hadv[k] + h_v * aux_gm.dqtdt_hadv[k])
        end
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tendencies_gm.q_liq[k] += aux_gm.dqldt[k]
            tendencies_gm.q_ice[k] += aux_gm.dqidt[k]
        end

        # LES specific forcings
        if TC.force_type(force) <: TC.ForcingLES
            T_fluc = aux_gm.dTdt_fluc[k]

            gm_U_nudge_k = (aux_gm.u_nudge[k] - prog_gm_u[k]) / force.wind_nudge_τᵣ
            gm_V_nudge_k = (aux_gm.v_nudge[k] - prog_gm_v[k]) / force.wind_nudge_τᵣ

            Γᵣ = TC.compute_les_Γᵣ(grid.zc[k].z, force.scalar_nudge_τᵣ, force.scalar_nudge_zᵢ, force.scalar_nudge_zᵣ)
            gm_T_nudge_k = Γᵣ * (aux_gm.T_nudge[k] - aux_gm.T[k])
            gm_q_tot_nudge_k = Γᵣ * (aux_gm.qt_nudge[k] - aux_gm.q_tot[k])
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                gm_q_liq_nudge_k = Γᵣ * (aux_gm.ql_nudge[k] - prog_gm.q_liq[k])
                gm_q_ice_nudge_k = Γᵣ * (aux_gm.qi_nudge[k] - prog_gm.q_ice[k])
            end

            tendencies_gm.ρq_tot[k] += ρ_c[k] * (gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k])
            tendencies_gm.ρe_tot[k] +=
                ρ_c[k] * (cv_m * (gm_T_nudge_k + T_fluc) + h_v * (gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k]))
            tendencies_gm_u[k] += gm_U_nudge_k
            tendencies_gm_v[k] += gm_V_nudge_k
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                tendencies_gm.q_liq[k] += aux_gm.dqldt_hadv[k] + gm_q_liq_nudge_k + aux_gm.dqldt_fluc[k]
                tendencies_gm.q_ice[k] += aux_gm.dqidt_hadv[k] + gm_q_ice_nudge_k + aux_gm.dqidt_fluc[k]
            end
        end

        # Apply precipitation tendencies
        tendencies_gm.ρq_tot[k] +=
            ρ_c[k] * (
                aux_bulk.qt_tendency_precip_formation[k] +
                aux_en.qt_tendency_precip_formation[k] +
                aux_tc.qt_tendency_precip_sinks[k]
            )

        tendencies_gm.ρe_tot[k] +=
            ρ_c[k] * (
                aux_bulk.e_tot_tendency_precip_formation[k] +
                aux_en.e_tot_tendency_precip_formation[k] +
                aux_tc.e_tot_tendency_precip_sinks[k]
            )

        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tendencies_gm.q_liq[k] += aux_bulk.ql_tendency_precip_formation[k] + aux_en.ql_tendency_precip_formation[k]
            tendencies_gm.q_ice[k] += aux_bulk.qi_tendency_precip_formation[k] + aux_en.qi_tendency_precip_formation[k]
        end
    end

    TC.compute_sgs_flux!(edmf, grid, state, surf, param_set)
    sgs_flux_h_tot = aux_gm_f.sgs_flux_h_tot
    sgs_flux_q_tot = aux_gm_f.sgs_flux_q_tot
    sgs_flux_u = aux_gm_f.sgs_flux_u
    sgs_flux_v = aux_gm_f.sgs_flux_v
    # apply surface BC as SGS flux at lowest level
    sgs_flux_h_tot[kf_surf] = surf.ρe_tot_flux
    sgs_flux_q_tot[kf_surf] = surf.ρq_tot_flux
    sgs_flux_u[kf_surf] = surf.ρu_flux
    sgs_flux_v[kf_surf] = surf.ρv_flux

    tends_ρe_tot = tendencies_gm.ρe_tot
    tends_ρq_tot = tendencies_gm.ρq_tot
    tends_u = TC.tendencies_grid_mean_u(state)
    tends_v = TC.tendencies_grid_mean_v(state)

    ∇sgs = CCO.DivergenceF2C()
    @. tends_ρe_tot += -∇sgs(wvec(sgs_flux_h_tot))
    @. tends_ρq_tot += -∇sgs(wvec(sgs_flux_q_tot))
    @. tends_u += -∇sgs(wvec(sgs_flux_u)) / ρ_c
    @. tends_v += -∇sgs(wvec(sgs_flux_v)) / ρ_c

    if edmf.moisture_model isa TC.NonEquilibriumMoisture
        sgs_flux_q_liq = aux_gm_f.sgs_flux_q_liq
        sgs_flux_q_ice = aux_gm_f.sgs_flux_q_ice
        sgs_flux_q_liq[kf_surf] = surf.ρq_liq_flux
        sgs_flux_q_ice[kf_surf] = surf.ρq_ice_flux
        tends_q_liq = tendencies_gm.q_liq
        tends_q_ice = tendencies_gm.q_ice
        @. tends_q_liq += -∇sgs(wvec(sgs_flux_q_liq)) / ρ_c
        @. tends_q_ice += -∇sgs(wvec(sgs_flux_q_ice)) / ρ_c
    end

    return nothing
end
