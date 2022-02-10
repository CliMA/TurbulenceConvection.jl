import StaticArrays
const SA = StaticArrays

import SurfaceFluxes
const SF = SurfaceFluxes
const UF = SF.UniversalFunctions

function get_surface(
    surf_params::TC.FixedSurfaceFlux,
    grid::TC.Grid,
    state::TC.State,
    gm::TC.GridMeanVariables,
    t::Real,
    param_set::CP.AbstractEarthParameterSet,
)
    FT = eltype(grid)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    p0_f_surf = TC.face_ref_state(state).p0[kf_surf]
    p0_c_surf = TC.center_ref_state(state).p0[kc_surf]
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    Tsurface = TC.surface_temperature(surf_params, t)
    qsurface = TC.surface_q_tot(surf_params, t)
    shf = TC.sensible_heat_flux(surf_params, t)
    lhf = TC.latent_heat_flux(surf_params, t)
    Ri_bulk_crit = surf_params.Ri_bulk_crit
    zrough = surf_params.zrough

    ts_sfc = TD.PhaseEquil_pTq(param_set, p0_f_surf, Tsurface, qsurface)
    ts_in = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    universal_func = UF.Businger()
    scheme = SF.FVScheme()

    bflux = SF.compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)
    zi = TC.get_inversion(grid, state, param_set, Ri_bulk_crit)
    convective_vel = TC.get_wstar(bflux, zi) # yair here zi in TRMM should be adjusted

    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    kwargs = (;
        state_in = vals_int,
        state_sfc = vals_sfc,
        shf = shf,
        lhf = lhf,
        z0m = zrough,
        z0b = zrough,
        gustiness = convective_vel,
    )
    sc = if TC.fixed_ustar(surf_params)
        SF.FluxesAndFrictionVelocity{FT}(; kwargs..., ustar = surf_params.ustar)
    else
        SF.Fluxes{FT}(; kwargs...)
    end
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    return TC.SurfaceBase{FT}(;
        shf = shf,
        lhf = lhf,
        ustar = TC.fixed_ustar(surf_params) ? surf_params.ustar : result.ustar,
        bflux = bflux,
        obukhov_length = result.L_MO,
        cm = result.Cd,
        ch = result.Ch,
        ρu_flux = surf_params.zero_uv_fluxes ? FT(0) : result.ρτxz,
        ρv_flux = surf_params.zero_uv_fluxes ? FT(0) : result.ρτyz,
        ρθ_liq_ice_flux = shf / TD.cp_m(ts_in),
        ρq_tot_flux = lhf / TD.latent_heat_vapor(ts_in),
        wstar = convective_vel,
    )
end

function get_surface(
    surf_params::TC.FixedSurfaceCoeffs,
    grid::TC.Grid,
    state::TC.State,
    gm::TC.GridMeanVariables,
    t::Real,
    param_set::CP.AbstractEarthParameterSet,
)
    FT = eltype(grid)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    p0_f_surf = TC.face_ref_state(state).p0[kf_surf]
    p0_c_surf = TC.center_ref_state(state).p0[kc_surf]
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    Tsurface = TC.surface_temperature(surf_params, t)
    qsurface = TC.surface_q_tot(surf_params, t)
    zrough = surf_params.zrough
    cm = surf_params.cm
    ch = surf_params.ch
    Ri_bulk_crit = surf_params.Ri_bulk_crit

    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_sfc = TD.PhaseEquil_pθq(param_set, p0_f_surf, Tsurface, qsurface)
    ts_in = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.Coefficients{FT}(state_in = vals_int, state_sfc = vals_sfc, Cd = cm, Ch = ch, z0m = zrough, z0b = zrough)
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    lhf = result.lhf
    shf = result.shf

    zi = TC.get_inversion(grid, state, param_set, Ri_bulk_crit)
    convective_vel = TC.get_wstar(result.buoy_flux, zi)
    return TC.SurfaceBase{FT}(;
        cm = result.Cd,
        ch = result.Ch,
        obukhov_length = result.L_MO,
        lhf = lhf,
        shf = shf,
        ustar = result.ustar,
        ρu_flux = result.ρτxz,
        ρv_flux = result.ρτyz,
        ρθ_liq_ice_flux = shf / TD.cp_m(ts_in),
        ρq_tot_flux = lhf / TD.latent_heat_vapor(ts_in),
        bflux = result.buoy_flux,
        wstar = convective_vel,
    )
end

function get_surface(
    surf_params::TC.MoninObukhovSurface,
    grid::TC.Grid,
    state::TC.State,
    gm::TC.GridMeanVariables,
    t::Real,
    param_set::CP.AbstractEarthParameterSet,
)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = TC.face_ref_state(state).p0[kf_surf]
    p0_c_surf = TC.center_ref_state(state).p0[kc_surf]
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    Tsurface = TC.surface_temperature(surf_params, t)
    qsurface = TC.surface_q_tot(surf_params, t)
    shf = TC.sensible_heat_flux(surf_params, t)
    lhf = TC.latent_heat_flux(surf_params, t)
    zrough = surf_params.zrough
    Ri_bulk_crit = surf_params.Ri_bulk_crit

    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    ts_sfc = TD.PhaseEquil_pTq(param_set, p0_f_surf, Tsurface, qsurface)
    ts_in = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, qsurface)

    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(state_in = vals_int, state_sfc = vals_sfc, z0m = zrough, z0b = zrough)
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    lhf = result.lhf
    shf = result.shf
    zi = TC.get_inversion(grid, state, param_set, Ri_bulk_crit)
    convective_vel = TC.get_wstar(result.buoy_flux, zi)
    return TC.SurfaceBase{FT}(;
        cm = result.Cd,
        ch = result.Ch,
        obukhov_length = result.L_MO,
        lhf = lhf,
        shf = shf,
        ustar = result.ustar,
        ρu_flux = result.ρτxz,
        ρv_flux = result.ρτyz,
        ρθ_liq_ice_flux = shf / TD.cp_m(ts_in),
        ρq_tot_flux = lhf / TD.latent_heat_vapor(ts_in),
        bflux = result.buoy_flux,
        wstar = convective_vel,
    )
end

function get_surface(
    surf_params::TC.SullivanPattonSurface,
    grid::TC.Grid,
    state::TC.State,
    gm::TC.GridMeanVariables,
    t::Real,
    param_set::CP.AbstractEarthParameterSet,
)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = TC.face_ref_state(state).p0[kf_surf]
    p0_c_surf = TC.center_ref_state(state).p0[kc_surf]
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    Tsurface = surface_temperature(surf_params, t)
    zrough = surf_params.zrough
    Ri_bulk_crit = surf_params.Ri_bulk_crit

    phase_part = TD.PhasePartition(q_tot_gm_surf, 0.0, 0.0)
    pvg = TD.saturation_vapor_pressure(param_set, TD.PhaseEquil, Tsurface)
    qsurface = TD.q_vap_saturation_from_density(param_set, Tsurface, ρ0_f_surf, pvg)
    θ_star = TD.liquid_ice_pottemp_given_pressure(param_set, Tsurface, p0_f_surf, phase_part)

    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    ts_sfc = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_star, qsurface)
    ts_in = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, qsurface)

    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(state_in = vals_int, state_sfc = vals_sfc, z0m = zrough, z0b = zrough)
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    lhf = result.lhf
    shf = result.shf
    zi = TC.get_inversion(grid, state, param_set, Ri_bulk_crit)
    convective_vel = TC.get_wstar(result.buoy_flux, zi)
    return TC.SurfaceBase{FT}(;
        cm = result.Cd,
        ch = result.Ch,
        obukhov_length = result.L_MO,
        lhf = lhf * 0.0, # TODO: why are we zero-ing this out?
        shf = shf,
        ustar = result.ustar,
        ρu_flux = result.ρτxz,
        ρv_flux = result.ρτyz,
        ρθ_liq_ice_flux = shf / TD.cp_m(ts_in),
        ρq_tot_flux = lhf / TD.latent_heat_vapor(ts_in),
        bflux = result.buoy_flux,
        wstar = convective_vel,
    )
end
