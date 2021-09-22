#####
##### BaseCase methods
#####

function free_convection_windspeed(self::SurfaceBase, GMV::GridMeanVariables, ::BaseCase)
    theta_rho = center_field(self.grid)
    param_set = parameter_set(GMV)

    # Need to get theta_rho
    @inbounds for k in real_center_indices(self.grid)
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[k], GMV.H.values[k], GMV.QT.values[k])
        theta_rho[k] = TD.virtual_pottemp(ts)
    end
    zi = get_inversion(param_set, theta_rho, GMV.U.values, GMV.V.values, self.grid, self.Ri_bulk_crit)
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
    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], self.Tsurface, self.qsurface)
    rho_tflux = self.shf / TD.cp_m(ts)

    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    self.rho_qtflux = self.lhf / TD.latent_heat_vapor(param_set, self.Tsurface)

    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], GMV.H.values[kc_surf], GMV.QT.values[kc_surf])
    Π = TD.exner(ts)
    self.rho_hflux = rho_tflux / Π
    self.bflux = buoyancy_flux(
        param_set,
        self.shf,
        self.lhf,
        GMV.T.values[kc_surf],
        GMV.QT.values[kc_surf],
        self.ref_state.alpha0[kf_surf],
        ts,
    )

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
                print("GMV.U.values[kc_surf] ==>", GMV.U.values[kc_surf])
                print("GMV.v.values[kc_surf] ==>", GMV.V.values[kc_surf])
                print("GMV.QT.values[kc_surf] ==>", GMV.QT.values[kc_surf])
                print("self.ref_state.alpha0[kf_surf] ==>", self.ref_state.alpha0[kf_surf])
            end
        end

        self.ustar = compute_ustar(param_set, self.windspeed, self.bflux, self.zrough, self.grid.zc[kc_surf])
    end

    von_karman_const = CPSGS.von_karman_const(param_set)
    self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / von_karman_const
    self.rho_uflux = -self.ref_state.rho0[kf_surf] * self.ustar * self.ustar / self.windspeed * GMV.U.values[kc_surf]
    self.rho_vflux = -self.ref_state.rho0[kf_surf] * self.ustar * self.ustar / self.windspeed * GMV.V.values[kc_surf]
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
    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], self.Tsurface, self.qsurface)
    windspeed =
        max(sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf]), 0.01)
    cp_ = TD.cp_m(ts)
    lv = TD.latent_heat_vapor(ts)

    self.rho_qtflux = -self.cq * windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.ref_state.rho0[kf_surf]
    self.lhf = lv * self.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], GMV.H.values[kc_surf], GMV.QT.values[kc_surf])
    Π = TD.exner(ts)
    self.rho_hflux = -self.ch * windspeed * (GMV.H.values[kc_surf] - self.Tsurface / Π) * self.ref_state.rho0[kf_surf]
    self.shf = cp_ * self.rho_hflux

    self.bflux = buoyancy_flux(
        param_set,
        self.shf,
        self.lhf,
        GMV.T.values[kc_surf],
        GMV.QT.values[kc_surf],
        self.ref_state.alpha0[kf_surf],
        ts,
    )

    self.ustar = sqrt(self.cm) * windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / von_karman_const
    end

    self.rho_uflux = -self.ref_state.rho0[kf_surf] * self.ustar * self.ustar / windspeed * GMV.U.values[kc_surf]
    self.rho_vflux = -self.ref_state.rho0[kf_surf] * self.ustar * self.ustar / windspeed * GMV.V.values[kc_surf]
    return
end

function update(self::SurfaceBase{SurfaceMoninObukhov}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)
    pvg = TD.saturation_vapor_pressure(param_set, self.Tsurface, TD.Liquid())
    self.qsurface = TD.q_vap_saturation_from_pressure(param_set, self.Tsurface, self.ref_state.rho0[kf_surf], pvg)
    zb = self.grid.zc[kc_surf]
    lv = TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])

    phase_part = TD.PhasePartition(self.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, self.Tsurface, self.ref_state.Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, self.ref_state.Pg, h_star, self.qsurface)
    theta_rho_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[kc_surf], GMV.H.values[kc_surf], self.qsurface)
    theta_rho_b = TD.virtual_pottemp(ts_b)

    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(param_set, Ri, zb, self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[kc_surf]) * self.ref_state.rho0[kf_surf]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[kc_surf]) * self.ref_state.rho0[kf_surf]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[kc_surf] - h_star) * self.ref_state.rho0[kf_surf]
    self.rho_qtflux =
        -self.ch * self.windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.ref_state.rho0[kf_surf]
    self.lhf = lv * self.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], self.Tsurface, self.qsurface)
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.bflux = buoyancy_flux(
        param_set,
        self.shf,
        self.lhf,
        GMV.T.values[kc_surf],
        GMV.QT.values[kc_surf],
        self.ref_state.alpha0[kf_surf],
        ts,
    )
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / von_karman_const
    end

    return
end

function update(self::SurfaceBase{SurfaceMoninObukhovDry}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(self.grid)
    kf_surf = kf_surface(self.grid)
    pvg = TD.saturation_vapor_pressure(param_set, self.Tsurface, TD.Liquid())
    self.qsurface = TD.q_vap_saturation_from_pressure(param_set, self.Tsurface, self.ref_state.rho0[kf_surf], pvg)
    zb = self.grid.zc[kc_surf]
    lv = TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])

    phase_part = TD.PhasePartition(self.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, self.Tsurface, self.ref_state.Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, self.ref_state.Pg, h_star, self.qsurface)
    theta_rho_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[kc_surf], GMV.H.values[kc_surf], self.qsurface)
    theta_rho_b = TD.virtual_pottemp(ts_b)

    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(param_set, Ri, zb, self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[kc_surf]) * self.ref_state.rho0[kf_surf]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[kc_surf]) * self.ref_state.rho0[kf_surf]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[kc_surf] - h_star) * self.ref_state.rho0[kf_surf]
    self.rho_qtflux =
        -self.ch * self.windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.ref_state.rho0[kf_surf]
    self.lhf = lv * 0.0

    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kf_surf], self.Tsurface, self.qsurface)
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.bflux =
        buoyancy_flux(param_set, self.shf, self.lhf, GMV.T.values[kc_surf], 0.0, self.ref_state.alpha0[kf_surf], ts)
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / von_karman_const
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
    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0[kc_surf], self.Tsurface, self.qsurface)
    lv = TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])
    T0 = self.ref_state.p0_half[kc_surf] * self.ref_state.alpha0_half[kc_surf] / R_d

    theta_flux = 0.24
    phase_part = TD.PhasePartition(GMV.QT.values[kc_surf], 0.0, 0.0)
    Π = TD.exner_given_pressure(param_set, self.ref_state.p0[kf_surf], phase_part)
    self.bflux = g * theta_flux * Π / T0
    pvg = TD.saturation_vapor_pressure(param_set, self.Tsurface, TD.Liquid())
    self.qsurface = TD.q_vap_saturation_from_pressure(param_set, self.Tsurface, self.ref_state.rho0[kf_surf], pvg)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, self.Tsurface, self.ref_state.Pg, phase_part)

    ts_g = TD.PhaseEquil_pθq(param_set, self.ref_state.Pg, h_star, self.qsurface)
    theta_rho_g = TD.virtual_pottemp(ts_g)

    ts_b = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[kc_surf], GMV.H.values[kc_surf], self.qsurface)
    theta_rho_b = TD.virtual_pottemp(ts_b)


    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(param_set, Ri, zb, self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[kc_surf]) * self.ref_state.rho0[kf_surf]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[kc_surf]) * self.ref_state.rho0[kf_surf]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[kc_surf] - h_star) * self.ref_state.rho0[kf_surf]
    self.rho_qtflux =
        -self.ch * self.windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.ref_state.rho0[kf_surf]
    self.lhf = lv * self.rho_qtflux
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / von_karman_const
    end
    return
end
