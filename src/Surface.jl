#####
##### BaseCase methods
#####

function free_convection_windspeed(self::SurfaceBase, GMV::GridMeanVariables, ::BaseCase)
    θ_ρ = center_field(self.grid)
    param_set = parameter_set(GMV)

    # Need to get θ_ρ
    @inbounds for k in real_center_indices(self.grid)
        ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[k], GMV.H.values[k], GMV.QT.values[k])
        θ_ρ[k] = TD.virtual_pottemp(ts)
    end
    zi = get_inversion(param_set, θ_ρ, GMV.U.values, GMV.V.values, self.grid, self.Ri_bulk_crit)
    wstar = get_wstar(self.bflux, zi) # yair here zi in TRMM should be adjusted
    self.windspeed = sqrt(self.windspeed * self.windspeed + (1.2 * wstar) * (1.2 * wstar))
    return
end

#####
##### Default SurfaceBase behavior
#####

free_convection_windspeed(self::SurfaceBase, GMV::GridMeanVariables) = free_convection_windspeed(self, GMV, BaseCase())

initialize(self::SurfaceBase) = nothing

#####
##### SurfaceNone
#####

function update(self::SurfaceBase{SurfaceNone}, GMV::GridMeanVariables)
    # JH: assigning small fluxed so that simulation won"t crash when computing mixing length
    param_set = parameter_set(GMV)
    kc_surf = kc_surface(self.grid)
    self.windspeed = 0.0001
    self.zrough = 1e-4
    self.bflux = 1e-4
    self.ustar = compute_ustar(param_set, self.windspeed, self.bflux, self.zrough, self.grid.zc[kc_surf])
    return
end
free_convection_windspeed(self::SurfaceBase{SurfaceNone}, GMV::GridMeanVariables) = nothing

function update(self::SurfaceBase{SurfaceFixedFlux}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)
    zb = self.grid.zc[kc_surf]
    p0_f_surf = self.ref_state.p0[kf_surf]
    ρ0_f_surf = self.ref_state.rho0[kf_surf]
    α0_f_surf = self.ref_state.alpha0[kf_surf]
    u_gm_surf = GMV.U.values[kc_surf]
    v_gm_surf = GMV.V.values[kc_surf]
    q_tot_gm_surf = GMV.QT.values[kc_surf]
    θ_liq_ice_gm_surf = GMV.H.values[kc_surf]
    T_gm_surf = GMV.T.values[kc_surf]

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, self.Tsurface, self.qsurface)
    rho_tflux = self.shf / TD.cp_m(ts)

    self.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    self.rho_qtflux = self.lhf / TD.latent_heat_vapor(param_set, self.Tsurface)

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    Π = TD.exner(ts)
    self.rho_hflux = rho_tflux / Π
    self.bflux = buoyancy_flux(param_set, self.shf, self.lhf, T_gm_surf, q_tot_gm_surf, α0_f_surf, ts)

    if !self.ustar_fixed
        # Correction to windspeed for free convective cases (Beljaars, QJRMS (1994), 121, pp. 255-270)
        # Value 1.2 is empirical, but should be O(1)
        if self.windspeed < 0.1  # Limit here is heuristic
            if self.bflux > 0.0
                free_convection_windspeed(self, GMV)
            else
                print("WARNING: Low windspeed + stable conditions, need to check ustar computation")
                print("self.bflux ==>", self.bflux)
                print("self.shf ==>", self.shf)
                print("self.lhf ==>", self.lhf)
                print("u_gm_surf ==>", u_gm_surf)
                print("v_gm_surf ==>", v_gm_surf)
                print("q_tot_gm_surf ==>", q_tot_gm_surf)
                print("α0_f_surf ==>", α0_f_surf)
                error("Terminating execution.")
            end
        end

        self.ustar = compute_ustar(param_set, self.windspeed, self.bflux, self.zrough, zb)
    end

    von_karman_const = CPSGS.von_karman_const(param_set)
    self.obukhov_length = -self.ustar * self.ustar^2 / self.bflux / von_karman_const
    self.rho_uflux = -ρ0_f_surf * self.ustar^2 / self.windspeed * u_gm_surf
    self.rho_vflux = -ρ0_f_surf * self.ustar^2 / self.windspeed * v_gm_surf
    return
end

# Cases such as Rico which provide values of transfer coefficients
function initialize(self::SurfaceBase{SurfaceFixedCoeffs})
    param_set = parameter_set(self.ref_state)
    kf_surf = kf_surface(self.grid)
    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], self.Tsurface, self.qsurface)
    pvg = TD.saturation_vapor_pressure(ts, TD.Liquid())
    self.qsurface = TD.q_vap_saturation(ts)
    return
end

function update(self::SurfaceBase{SurfaceFixedCoeffs}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)
    p0_f_surf = self.ref_state.p0[kf_surf]
    ρ0_f_surf = self.ref_state.rho0[kf_surf]
    α0_f_surf = self.ref_state.alpha0[kf_surf]
    u_gm_surf = GMV.U.values[kc_surf]
    v_gm_surf = GMV.V.values[kc_surf]
    T_gm_surf = GMV.T.values[kc_surf]
    q_tot_gm_surf = GMV.QT.values[kc_surf]
    θ_liq_ice_gm_surf = GMV.H.values[kc_surf]

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, self.Tsurface, self.qsurface)
    windspeed = max(sqrt(u_gm_surf^2 + v_gm_surf^2), 0.01)
    cp_ = TD.cp_m(ts)
    lv = TD.latent_heat_vapor(ts)

    self.rho_qtflux = -self.cq * windspeed * (q_tot_gm_surf - self.qsurface) * ρ0_f_surf
    self.lhf = lv * self.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    Π = TD.exner(ts)
    self.rho_hflux = -self.ch * windspeed * (θ_liq_ice_gm_surf - self.Tsurface / Π) * ρ0_f_surf
    self.shf = cp_ * self.rho_hflux

    self.bflux = buoyancy_flux(param_set, self.shf, self.lhf, T_gm_surf, q_tot_gm_surf, α0_f_surf, ts)

    self.ustar = sqrt(self.cm) * windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar^3 / self.bflux / von_karman_const
    end

    self.rho_uflux = -ρ0_f_surf * self.ustar^2 / windspeed * u_gm_surf
    self.rho_vflux = -ρ0_f_surf * self.ustar^2 / windspeed * v_gm_surf
    return
end

function update(self::SurfaceBase{SurfaceMoninObukhov}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)
    zb = self.grid.zc[kc_surf]
    p0_f_surf = self.ref_state.p0[kf_surf]
    ρ0_f_surf = self.ref_state.rho0[kf_surf]
    α0_f_surf = self.ref_state.alpha0[kf_surf]
    u_gm_surf = GMV.U.values[kc_surf]
    v_gm_surf = GMV.V.values[kc_surf]
    q_tot_gm_surf = GMV.QT.values[kc_surf]
    θ_liq_ice_gm_surf = GMV.H.values[kc_surf]
    T_gm_surf = GMV.T.values[kc_surf]

    pvg = TD.saturation_vapor_pressure(param_set, self.Tsurface, TD.Liquid())
    self.qsurface = TD.q_vap_saturation_from_pressure(param_set, self.Tsurface, ρ0_f_surf, pvg)
    lv = TD.latent_heat_vapor(param_set, T_gm_surf)

    phase_part = TD.PhasePartition(self.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, self.Tsurface, self.ref_state.Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, self.ref_state.Pg, h_star, self.qsurface)
    θ_ρ_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[kc_surf], θ_liq_ice_gm_surf, self.qsurface)
    θ_ρ_b = TD.virtual_pottemp(ts_b)

    self.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    Nb2 = g / θ_ρ_g * (θ_ρ_b - θ_ρ_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(param_set, Ri, zb, self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * u_gm_surf * ρ0_f_surf
    self.rho_vflux = -self.cm * self.windspeed * v_gm_surf * ρ0_f_surf

    self.rho_hflux = -self.ch * self.windspeed * (θ_liq_ice_gm_surf - h_star) * ρ0_f_surf
    self.rho_qtflux = -self.ch * self.windspeed * (q_tot_gm_surf - self.qsurface) * ρ0_f_surf
    self.lhf = lv * self.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, self.Tsurface, self.qsurface)
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.bflux = buoyancy_flux(param_set, self.shf, self.lhf, T_gm_surf, q_tot_gm_surf, α0_f_surf, ts)
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar^3 / self.bflux / von_karman_const
    end

    return
end

function update(self::SurfaceBase{SurfaceMoninObukhovDry}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)

    zb = self.grid.zc[kc_surf]
    p0_f_surf = self.ref_state.p0[kf_surf]
    p0_c_surf = self.ref_state.p0_half[kf_surf]
    ρ0_f_surf = self.ref_state.rho0[kf_surf]
    α0_f_surf = self.ref_state.alpha0[kf_surf]
    u_gm_surf = GMV.U.values[kc_surf]
    v_gm_surf = GMV.V.values[kc_surf]
    q_tot_gm_surf = GMV.QT.values[kc_surf]
    θ_liq_ice_gm_surf = GMV.H.values[kc_surf]
    T_gm_surf = GMV.T.values[kc_surf]

    pvg = TD.saturation_vapor_pressure(param_set, self.Tsurface, TD.Liquid())
    self.qsurface = TD.q_vap_saturation_from_pressure(param_set, self.Tsurface, ρ0_f_surf, pvg)
    lv = TD.latent_heat_vapor(param_set, T_gm_surf)

    phase_part = TD.PhasePartition(self.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, self.Tsurface, self.ref_state.Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, self.ref_state.Pg, h_star, self.qsurface)
    θ_ρ_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, self.qsurface)
    θ_ρ_b = TD.virtual_pottemp(ts_b)

    self.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    Nb2 = g / θ_ρ_g * (θ_ρ_b - θ_ρ_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(param_set, Ri, zb, self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * u_gm_surf * ρ0_f_surf
    self.rho_vflux = -self.cm * self.windspeed * v_gm_surf * ρ0_f_surf

    self.rho_hflux = -self.ch * self.windspeed * (θ_liq_ice_gm_surf - h_star) * ρ0_f_surf
    self.rho_qtflux = -self.ch * self.windspeed * (q_tot_gm_surf - self.qsurface) * ρ0_f_surf
    self.lhf = lv * 0.0

    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], self.Tsurface, self.qsurface)
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.bflux = buoyancy_flux(param_set, self.shf, self.lhf, T_gm_surf, 0.0, α0_f_surf, ts)
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar^3 / self.bflux / von_karman_const
    end

    return
end

function update(self::SurfaceBase{SurfaceSullivanPatton}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    R_d = CPP.R_d(param_set)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)

    zb = self.grid.zc[kc_surf]
    p0_f_surf = self.ref_state.p0[kf_surf]
    p0_c_surf = self.ref_state.p0_half[kf_surf]
    ρ0_f_surf = self.ref_state.rho0[kf_surf]
    α0_c_surf = self.ref_state.alpha0_half[kc_surf]
    u_gm_surf = GMV.U.values[kc_surf]
    v_gm_surf = GMV.V.values[kc_surf]
    q_tot_gm_surf = GMV.QT.values[kc_surf]
    θ_liq_ice_gm_surf = GMV.H.values[kc_surf]
    T_gm_surf = GMV.T.values[kc_surf]

    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, self.Tsurface, self.qsurface)
    lv = TD.latent_heat_vapor(param_set, T_gm_surf)
    T0 = p0_c_surf * α0_c_surf / R_d # TODO: can we use a thermo state here?

    θ_flux = 0.24
    phase_part = TD.PhasePartition(q_tot_gm_surf, 0.0, 0.0)
    Π = TD.exner_given_pressure(param_set, p0_f_surf, phase_part)
    self.bflux = g * θ_flux * Π / T0
    pvg = TD.saturation_vapor_pressure(param_set, self.Tsurface, TD.Liquid())
    self.qsurface = TD.q_vap_saturation_from_pressure(param_set, self.Tsurface, ρ0_f_surf, pvg)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, self.Tsurface, self.ref_state.Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, self.ref_state.Pg, h_star, self.qsurface)
    θ_ρ_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, self.qsurface)
    θ_ρ_b = TD.virtual_pottemp(ts_b)

    self.windspeed = sqrt(u_gm_surf^2 + v_gm_surf^2)
    Nb2 = g / θ_ρ_g * (θ_ρ_b - θ_ρ_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(param_set, Ri, zb, self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * u_gm_surf * ρ0_f_surf
    self.rho_vflux = -self.cm * self.windspeed * v_gm_surf * ρ0_f_surf

    self.rho_hflux = -self.ch * self.windspeed * (θ_liq_ice_gm_surf - h_star) * ρ0_f_surf
    self.rho_qtflux = -self.ch * self.windspeed * (q_tot_gm_surf - self.qsurface) * ρ0_f_surf
    self.lhf = lv * self.rho_qtflux
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar^3 / self.bflux / von_karman_const
    end
    return
end
