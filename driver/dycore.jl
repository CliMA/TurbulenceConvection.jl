import UnPack
import ClimaCore
const CC = ClimaCore

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
)
cent_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., cent_aux_vars_gm(FT)..., TC.cent_aux_vars_edmf(FT, n_up)...)

# Face only
face_aux_vars_gm(FT) = (; massflux_s = FT(0), diffusive_flux_s = FT(0), total_flux_s = FT(0), f_rad = FT(0))
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

function satadjust(gm::TC.GridMeanVariables, grid, state)
    p0_c = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    param_set = TC.parameter_set(gm)
    @inbounds for k in TC.real_center_indices(grid)
        θ_liq_ice = prog_gm.θ_liq_ice[k]
        q_tot = prog_gm.q_tot[k]
        ts = TC.thermo_state_pθq(param_set, p0_c[k], θ_liq_ice, q_tot)
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
    UnPack.@unpack edmf, grid, gm, case, aux, TS = params

    state = TC.State(prog, aux, tendencies)

    Δt = TS.dt
    param_set = TC.parameter_set(gm)
    surf = get_surface(case.surf_params, grid, state, gm, t, param_set)
    force = case.Fo
    radiation = case.Rad
    en_thermo = edmf.en_thermo
    precip_model = edmf.precip_model

    TC.affect_filter!(edmf, grid, state, gm, surf, case.casename, t)

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.

    TC.update_aux!(edmf, gm, grid, state, case, surf, param_set, t, Δt)

    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0

    # causes division error in dry bubble first time step
    TC.compute_precipitation_formation_tendencies(grid, state, edmf, precip_model, Δt, param_set)

    TC.microphysics(en_thermo, grid, state, precip_model, Δt, param_set)
    TC.compute_precipitation_sink_tendencies(precip_model, grid, state, gm, Δt)
    TC.compute_precipitation_advection_tendencies(precip_model, edmf, grid, state, gm)

    # compute tendencies
    TC.compute_gm_tendencies!(edmf, grid, state, surf, radiation, force, gm)
    TC.compute_up_tendencies!(edmf, grid, state, gm, surf)

    TC.compute_en_tendencies!(edmf, grid, state, param_set, Val(:tke), Val(:ρatke))
    TC.compute_en_tendencies!(edmf, grid, state, param_set, Val(:Hvar), Val(:ρaHvar))
    TC.compute_en_tendencies!(edmf, grid, state, param_set, Val(:QTvar), Val(:ρaQTvar))
    TC.compute_en_tendencies!(edmf, grid, state, param_set, Val(:HQTcov), Val(:ρaHQTcov))

    return nothing
end
