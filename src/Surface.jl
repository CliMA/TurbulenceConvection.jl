
function update(surf::SurfaceBase{SurfaceFixedFlux}, grid, state, gm::GridMeanVariables, param_set)
    FT = eltype(grid)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    ts_sfc = TD.PhaseEquil_pTq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    universal_func = UF.Businger()
    scheme = SF.FVScheme()

    surf.bflux = SF.compute_buoyancy_flux(param_set, surf.shf, surf.lhf, ts_in, ts_sfc, scheme)
    zi = get_inversion(grid, state, param_set, surf.Ri_bulk_crit)
    convective_vel = get_wstar(surf.bflux, zi) # yair here zi in TRMM should be adjusted

    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    if surf.ustar_fixed
        sc = SF.FluxesAndFrictionVelocity{FT}(;
            state_in = vals_int,
            state_sfc = vals_sfc,
            shf = surf.shf,
            lhf = surf.lhf,
            ustar = surf.ustar,
            z0m = surf.zrough,
            z0b = surf.zrough,
            gustiness = convective_vel,
        )
        result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    else
        sc = SF.Fluxes{FT}(
            state_in = vals_int,
            state_sfc = vals_sfc,
            shf = surf.shf,
            lhf = surf.lhf,
            z0m = surf.zrough,
            z0b = surf.zrough,
            gustiness = convective_vel,
        )
        result = SF.surface_conditions(param_set, sc, universal_func, scheme)
        surf.ustar = result.ustar
    end
    surf.obukhov_length = result.L_MO
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux = surf.shf / TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf / TD.latent_heat_vapor(ts_in)
    return
end

function update(surf::SurfaceBase{SurfaceFixedCoeffs}, grid, state, gm::GridMeanVariables, param_set)
    FT = eltype(grid)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_sfc = thermo_state_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.Coefficients{FT}(
        state_in = vals_int,
        state_sfc = vals_sfc,
        Cd = surf.cm,
        Ch = surf.ch,
        z0m = surf.zrough,
        z0b = surf.zrough,
    )
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux = surf.shf / TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf / TD.latent_heat_vapor(ts_in)
    surf.bflux = result.buoy_flux
    return
end

function update(surf::SurfaceBase{SurfaceMoninObukhov}, grid, state, gm::GridMeanVariables, param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    ts_sfc = thermo_state_pTq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)

    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(state_in = vals_int, state_sfc = vals_sfc, z0m = surf.zrough, z0b = surf.zrough)
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux = surf.shf / TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf / TD.latent_heat_vapor(ts_in)
    surf.bflux = result.buoy_flux
    return
end

function update(surf::SurfaceBase{SurfaceMoninObukhovDry}, grid, state, gm::GridMeanVariables, param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    ts_sfc = TD.PhaseEquil_pTq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)
    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(state_in = vals_int, state_sfc = vals_sfc, z0m = surf.zrough, z0b = surf.zrough)
    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf * 0.0
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux = surf.shf / TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf / TD.latent_heat_vapor(ts_in)
    surf.bflux = result.buoy_flux
    return
end

function update(surf::SurfaceBase{SurfaceSullivanPatton}, grid, state, gm::GridMeanVariables, param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    phase_part = TD.PhasePartition(q_tot_gm_surf, 0.0, 0.0)
    pvg = TD.saturation_vapor_pressure(param_set, TD.PhaseEquil, surf.Tsurface)
    surf.qsurface = TD.q_vap_saturation_from_density(param_set, surf.Tsurface, ρ0_f_surf, pvg)
    θ_star = TD.liquid_ice_pottemp_given_pressure(param_set, surf.Tsurface, p0_f_surf, phase_part)

    universal_func = UF.Businger()
    scheme = SF.FVScheme()
    ts_sfc = thermo_state_pθq(param_set, p0_f_surf, θ_star, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)

    u_sfc = SA.SVector{2, FT}(0, 0)
    u_in = SA.SVector{2, FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.InteriorValues(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(state_in = vals_int, state_sfc = vals_sfc, z0m = surf.zrough, z0b = surf.zrough)
    result = SF.surface_conditions(param_set, sc, universal_func, scheme)
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf * 0.0
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux = surf.shf / TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf / TD.latent_heat_vapor(ts_in)
    surf.bflux = result.buoy_flux
    return
end
