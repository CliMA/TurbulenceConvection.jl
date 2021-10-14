#####
##### BaseCase methods
#####

function free_convection_windspeed(surf::SurfaceBase, grid, state, gm::GridMeanVariables, param_set, ::BaseCase)
    θ_ρ = center_field(grid)
    p0_c = center_ref_state(state).p0

    # Need to get θ_ρ
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], gm.H.values[k], gm.QT.values[k])
        θ_ρ[k] = TD.virtual_pottemp(ts)
    end
    zi = get_inversion(param_set, θ_ρ, gm.U.values, gm.V.values, grid, surf.Ri_bulk_crit)
    wstar = get_wstar(surf.bflux, zi) # yair here zi in TRMM should be adjusted
    surf.windspeed = sqrt(surf.windspeed * surf.windspeed + (1.2 * wstar) * (1.2 * wstar))
    return
end

#####
##### Default SurfaceBase behavior
#####

free_convection_windspeed(surf::SurfaceBase, grid, state, gm::GridMeanVariables, param_set) =
    free_convection_windspeed(surf, grid, state, gm, param_set, BaseCase())

#####
##### SurfaceNone
#####

function update(surf::SurfaceBase{SurfaceNone}, grid, state, gm::GridMeanVariables, param_set)
    # JH: assigning small fluxed so that simulation won"t crash when computing mixing length
    kc_surf = kc_surface(grid)
    surf.windspeed = 0.0001
    surf.zrough = 1e-4
    surf.bflux = 1e-4
    surf.ustar = compute_friction_velocity(param_set, surf.windspeed, surf.bflux, surf.zrough, grid.zc[kc_surf])
    return
end
free_convection_windspeed(surf::SurfaceBase{SurfaceNone}, grid, state, gm::GridMeanVariables, param_set) = nothing

function update(surf::SurfaceBase{SurfaceFixedFlux}, grid, state, gm::GridMeanVariables, param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    zb = grid.zc[kc_surf]
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    u_gm_surf = gm.U.values[kc_surf]
    v_gm_surf = gm.V.values[kc_surf]
    q_tot_gm_surf = gm.QT.values[kc_surf]
    θ_liq_ice_gm_surf = gm.H.values[kc_surf]
    T_gm_surf = gm.T.values[kc_surf]

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    rho_tflux = surf.shf / TD.cp_m(ts)

    surf.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    surf.rho_qtflux = surf.lhf / TD.latent_heat_vapor(param_set, surf.Tsurface)

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    Π = TD.exner(ts)
    surf.rho_hflux = rho_tflux / Π
    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts)

    if !surf.ustar_fixed
        # Correction to windspeed for free convective cases (Beljaars, QJRMS (1994), 121, pp. 255-270)
        # Value 1.2 is empirical, but should be O(1)
        if surf.windspeed < 0.1  # Limit here is heuristic
            if surf.bflux > 0.0
                free_convection_windspeed(surf, grid, state, gm, param_set)
            else
                print("WARNING: Low windspeed + stable conditions, need to check ustar computation")
                print("surf.bflux ==>", surf.bflux)
                print("surf.shf ==>", surf.shf)
                print("surf.lhf ==>", surf.lhf)
                print("u_gm_surf ==>", u_gm_surf)
                print("v_gm_surf ==>", v_gm_surf)
                print("q_tot_gm_surf ==>", q_tot_gm_surf)
                print("α0_f_surf ==>", α0_f_surf)
                error("Terminating execution.")
            end
        end

        surf.ustar = compute_friction_velocity(param_set, surf.windspeed, surf.bflux, surf.zrough, zb)
    end

    von_karman_const = CPSGS.von_karman_const(param_set)
    surf.obukhov_length = obukhov_length(param_set, surf.ustar, surf.bflux)
    surf.rho_uflux = -ρ0_f_surf * surf.ustar^2 / surf.windspeed * u_gm_surf
    surf.rho_vflux = -ρ0_f_surf * surf.ustar^2 / surf.windspeed * v_gm_surf
    return
end

function update(surf::SurfaceBase{SurfaceFixedCoeffs}, grid, state, gm::GridMeanVariables, param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    u_gm_surf = gm.U.values[kc_surf]
    v_gm_surf = gm.V.values[kc_surf]
    T_gm_surf = gm.T.values[kc_surf]
    q_tot_gm_surf = gm.QT.values[kc_surf]
    θ_liq_ice_gm_surf = gm.H.values[kc_surf]

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    surf.windspeed = max(sqrt(u_gm_surf^2 + v_gm_surf^2), 0.01)
    cp_ = TD.cp_m(ts)
    lv = TD.latent_heat_vapor(ts)

    surf.rho_qtflux = -surf.cq * surf.windspeed * (q_tot_gm_surf - surf.qsurface) * ρ0_f_surf
    surf.lhf = lv * surf.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    Π = TD.exner(ts)
    surf.rho_hflux = -surf.ch * surf.windspeed * (θ_liq_ice_gm_surf - surf.Tsurface / Π) * ρ0_f_surf
    surf.shf = cp_ * surf.rho_hflux

    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts)

    surf.ustar = friction_velocity_given_windspeed(surf.cm, surf.windspeed)
    surf.obukhov_length = obukhov_length(param_set, surf.ustar, surf.bflux)

    surf.rho_uflux = -ρ0_f_surf * surf.ustar^2 / surf.windspeed * u_gm_surf
    surf.rho_vflux = -ρ0_f_surf * surf.ustar^2 / surf.windspeed * v_gm_surf
    return
end

function update(surf::SurfaceBase{SurfaceMoninObukhov}, grid, state, gm::GridMeanVariables, param_set)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    zb = grid.zc[kc_surf]
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    u_gm_surf = gm.U.values[kc_surf]
    v_gm_surf = gm.V.values[kc_surf]
    q_tot_gm_surf = gm.QT.values[kc_surf]
    θ_liq_ice_gm_surf = gm.H.values[kc_surf]
    T_gm_surf = gm.T.values[kc_surf]
    Pg = surf.ref_params.Pg

    pvg = TD.saturation_vapor_pressure(param_set, surf.Tsurface, TD.Liquid())
    surf.qsurface = TD.q_vap_saturation_from_density(param_set, surf.Tsurface, ρ0_f_surf, pvg)
    lv = TD.latent_heat_vapor(param_set, T_gm_surf)

    phase_part = TD.PhasePartition(surf.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, surf.Tsurface, Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, Pg, h_star, surf.qsurface)
    θ_ρ_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)
    θ_ρ_b = TD.virtual_pottemp(ts_b)

    surf.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    Nb2 = g / θ_ρ_g * (θ_ρ_b - θ_ρ_g) / zb
    Ri = Nb2 * zb * zb / (surf.windspeed * surf.windspeed)

    surf.cm, surf.ch, surf.obukhov_length = exchange_coefficients_byun(param_set, Ri, zb, surf.zrough)
    surf.rho_uflux = -surf.cm * surf.windspeed * u_gm_surf * ρ0_f_surf
    surf.rho_vflux = -surf.cm * surf.windspeed * v_gm_surf * ρ0_f_surf

    surf.rho_hflux = -surf.ch * surf.windspeed * (θ_liq_ice_gm_surf - h_star) * ρ0_f_surf
    surf.rho_qtflux = -surf.ch * surf.windspeed * (q_tot_gm_surf - surf.qsurface) * ρ0_f_surf
    surf.lhf = lv * surf.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    surf.shf = TD.cp_m(ts) * surf.rho_hflux

    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts)
    surf.ustar = friction_velocity_given_windspeed(surf.cm, surf.windspeed)

    return
end

function update(surf::SurfaceBase{SurfaceMoninObukhovDry}, grid, state, gm::GridMeanVariables, param_set)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)

    zb = grid.zc[kc_surf]
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    u_gm_surf = gm.U.values[kc_surf]
    v_gm_surf = gm.V.values[kc_surf]
    q_tot_gm_surf = gm.QT.values[kc_surf]
    θ_liq_ice_gm_surf = gm.H.values[kc_surf]
    T_gm_surf = gm.T.values[kc_surf]
    Pg = surf.ref_params.Pg

    pvg = TD.saturation_vapor_pressure(param_set, surf.Tsurface, TD.Liquid())
    surf.qsurface = TD.q_vap_saturation_from_density(param_set, surf.Tsurface, ρ0_f_surf, pvg)
    lv = TD.latent_heat_vapor(param_set, T_gm_surf)

    phase_part = TD.PhasePartition(surf.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, surf.Tsurface, Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, Pg, h_star, surf.qsurface)
    θ_ρ_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)
    θ_ρ_b = TD.virtual_pottemp(ts_b)

    surf.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    Nb2 = g / θ_ρ_g * (θ_ρ_b - θ_ρ_g) / zb
    Ri = Nb2 * zb * zb / (surf.windspeed * surf.windspeed)

    surf.cm, surf.ch, surf.obukhov_length = exchange_coefficients_byun(param_set, Ri, zb, surf.zrough)
    surf.rho_uflux = -surf.cm * surf.windspeed * u_gm_surf * ρ0_f_surf
    surf.rho_vflux = -surf.cm * surf.windspeed * v_gm_surf * ρ0_f_surf

    surf.rho_hflux = -surf.ch * surf.windspeed * (θ_liq_ice_gm_surf - h_star) * ρ0_f_surf
    surf.rho_qtflux = -surf.ch * surf.windspeed * (q_tot_gm_surf - surf.qsurface) * ρ0_f_surf
    surf.lhf = lv * 0.0

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    surf.shf = TD.cp_m(ts) * surf.rho_hflux

    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts)
    surf.ustar = friction_velocity_given_windspeed(surf.cm, surf.windspeed)

    return
end

function update(surf::SurfaceBase{SurfaceSullivanPatton}, grid, state, gm::GridMeanVariables, param_set)
    g = CPP.grav(param_set)
    R_d = CPP.R_d(param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)

    zb = grid.zc[kc_surf]

    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_c_surf = center_ref_state(state).α0[kc_surf]
    u_gm_surf = gm.U.values[kc_surf]
    v_gm_surf = gm.V.values[kc_surf]
    q_tot_gm_surf = gm.QT.values[kc_surf]
    θ_liq_ice_gm_surf = gm.H.values[kc_surf]
    T_gm_surf = gm.T.values[kc_surf]
    Pg = surf.ref_params.Pg

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    lv = TD.latent_heat_vapor(param_set, T_gm_surf)
    T0 = p0_c_surf * α0_c_surf / R_d # TODO: can we use a thermo state here?

    θ_flux = 0.24
    phase_part = TD.PhasePartition(q_tot_gm_surf, 0.0, 0.0)
    Π = TD.exner_given_pressure(param_set, p0_f_surf, phase_part)
    surf.bflux = g * θ_flux * Π / T0
    pvg = TD.saturation_vapor_pressure(param_set, surf.Tsurface, TD.Liquid())
    surf.qsurface = TD.q_vap_saturation_from_density(param_set, surf.Tsurface, ρ0_f_surf, pvg)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, surf.Tsurface, Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, Pg, h_star, surf.qsurface)
    θ_ρ_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)
    θ_ρ_b = TD.virtual_pottemp(ts_b)

    surf.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    Nb2 = g / θ_ρ_g * (θ_ρ_b - θ_ρ_g) / zb
    Ri = Nb2 * zb * zb / (surf.windspeed * surf.windspeed)

    surf.cm, surf.ch, surf.obukhov_length = exchange_coefficients_byun(param_set, Ri, zb, surf.zrough)
    surf.rho_uflux = -surf.cm * surf.windspeed * u_gm_surf * ρ0_f_surf
    surf.rho_vflux = -surf.cm * surf.windspeed * v_gm_surf * ρ0_f_surf

    surf.rho_hflux = -surf.ch * surf.windspeed * (θ_liq_ice_gm_surf - h_star) * ρ0_f_surf
    surf.rho_qtflux = -surf.ch * surf.windspeed * (q_tot_gm_surf - surf.qsurface) * ρ0_f_surf
    surf.lhf = lv * surf.rho_qtflux
    surf.shf = TD.cp_m(ts) * surf.rho_hflux

    surf.ustar = friction_velocity_given_windspeed(surf.cm, surf.windspeed)
    return
end
