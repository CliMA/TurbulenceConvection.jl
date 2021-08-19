initialize(self::SurfaceBase, ::BaseCase) = nothing
update(self::SurfaceBase, GMV::GridMeanVariables, ::BaseCase) = nothing

function free_convection_windspeed(self::SurfaceBase, GMV::GridMeanVariables, ::BaseCase)
    theta_rho = center_field(self.Gr)
    param_set = parameter_set(GMV)

    # Need to get theta_rho
    @inbounds for k in center_indicies(self.Gr)
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        theta_rho[k] = theta_rho_c(self.Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], qv)
    end
    zi = get_inversion(param_set, theta_rho, GMV.U.values, GMV.V.values, self.Gr, self.Ri_bulk_crit)
    wstar = get_wstar(self.bflux, zi) # yair here zi in TRMM should be adjusted
    self.windspeed = sqrt(self.windspeed * self.windspeed + (1.2 * wstar) * (1.2 * wstar))
    return
end

initialize(self::SurfaceBase{SurfaceNone}) = nothing

function update(self::SurfaceBase{SurfaceNone}, GMV::GridMeanVariables)
    # JH: assigning small fluxed so that simulation won"t crash when computing mixing length
    kc_surf = kc_surface(self.Gr)
    self.windspeed = 0.0001
    self.zrough = 1e-4
    self.bflux = 1e-4
    self.ustar = compute_ustar(self.windspeed, self.bflux, self.zrough, self.Gr.z_half[kc_surf])
    return
end
free_convection_windspeed(self::SurfaceBase{SurfaceNone}, GMV::GridMeanVariables) = nothing

initialize(self::SurfaceBase{SurfaceFixedFlux}) = nothing

function update(self::SurfaceBase{SurfaceFixedFlux}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    kc_surf = kc_surface(self.Gr)
    kf_surf = kf_surface(self.Gr)
    ts = TD.PhaseEquil_pθq(param_set, self.Ref.p0[kf_surf], self.Tsurface, self.qsurface)
    rho_tflux = self.shf / TD.cp_m(ts)

    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    self.rho_qtflux = self.lhf / TD.latent_heat_vapor(param_set, self.Tsurface)

    self.rho_hflux = rho_tflux / exner_c(self.Ref.Pg)
    self.bflux = buoyancy_flux(
        param_set,
        self.shf,
        self.lhf,
        GMV.T.values[kc_surf],
        GMV.QT.values[kc_surf],
        self.Ref.alpha0[kf_surf],
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
                print("self.Ref.alpha0[kf_surf] ==>", self.Ref.alpha0[kf_surf])
            end
        end

        self.ustar = compute_ustar(self.windspeed, self.bflux, self.zrough, self.Gr.z_half[kc_surf])
    end

    self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    self.rho_uflux = -self.Ref.rho0[kf_surf] * self.ustar * self.ustar / self.windspeed * GMV.U.values[kc_surf]
    self.rho_vflux = -self.Ref.rho0[kf_surf] * self.ustar * self.ustar / self.windspeed * GMV.V.values[kc_surf]
    return
end

function free_convection_windspeed(self::SurfaceBase{SurfaceFixedFlux}, GMV::GridMeanVariables)
    free_convection_windspeed(self, GMV, BaseCase())
    return
end

# Cases such as Rico which provide values of transfer coefficients
function initialize(self::SurfaceBase{SurfaceFixedCoeffs})
    pvg = pv_star(self.Tsurface)
    pdg = self.Ref.Pg - pvg
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    self.s_surface = (1.0 - self.qsurface) * sd_c(pdg, self.Tsurface) + self.qsurface * sv_c(pvg, self.Tsurface)
    return
end

function update(self::SurfaceBase{SurfaceFixedCoeffs}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    kc_surf = kc_surface(self.Gr)
    kf_surf = kf_surface(self.Gr)
    ts = TD.PhaseEquil_pθq(param_set, self.Ref.p0[kf_surf], self.Tsurface, self.qsurface)
    windspeed =
        max(sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf]), 0.01)
    cp_ = TD.cp_m(ts)
    lv = TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])

    self.rho_qtflux = -self.cq * windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.Ref.rho0[kf_surf]
    self.lhf = lv * self.rho_qtflux
    self.rho_hflux =
        -self.ch * windspeed * (GMV.H.values[kc_surf] - self.Tsurface / exner_c(self.Ref.Pg)) * self.Ref.rho0[kf_surf]
    self.shf = cp_ * self.rho_hflux

    self.bflux = buoyancy_flux(
        param_set,
        self.shf,
        self.lhf,
        GMV.T.values[kc_surf],
        GMV.QT.values[kc_surf],
        self.Ref.alpha0[kf_surf],
        ts,
    )

    self.ustar = sqrt(self.cm) * windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    end

    self.rho_uflux = -self.Ref.rho0[kf_surf] * self.ustar * self.ustar / windspeed * GMV.U.values[kc_surf]
    self.rho_vflux = -self.Ref.rho0[kf_surf] * self.ustar * self.ustar / windspeed * GMV.V.values[kc_surf]
    return
end

function free_convection_windspeed(self::SurfaceBase{SurfaceFixedCoeffs}, GMV::GridMeanVariables)
    free_convection_windspeed(self, GMV, BaseCase())
    return
end

initialize(self::SurfaceBase{SurfaceMoninObukhov}) = nothing

function update(self::SurfaceBase{SurfaceMoninObukhov}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(self.Gr)
    kf_surf = kf_surface(self.Gr)
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    zb = self.Gr.z_half[kc_surf]
    theta_rho_g = theta_rho_c(self.Ref.Pg, self.Tsurface, self.qsurface, self.qsurface)
    theta_rho_b = theta_rho_c(self.Ref.p0_half[kc_surf], GMV.T.values[kc_surf], self.qsurface, self.qsurface)
    TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])

    h_star = t_to_thetali_c(param_set, self.Ref.Pg, self.Tsurface, self.qsurface, 0.0, 0.0)

    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(Ri, self.Gr.z_half[kc_surf], self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[kc_surf]) * self.Ref.rho0[kf_surf]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[kc_surf]) * self.Ref.rho0[kf_surf]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[kc_surf] - h_star) * self.Ref.rho0[kf_surf]
    self.rho_qtflux = -self.ch * self.windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.Ref.rho0[kf_surf]
    self.lhf = lv * self.rho_qtflux

    ts = TD.PhaseEquil_pθq(param_set, self.Ref.p0[kf_surf], self.Tsurface, self.qsurface)
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.bflux = buoyancy_flux(
        param_set,
        self.shf,
        self.lhf,
        GMV.T.values[kc_surf],
        GMV.QT.values[kc_surf],
        self.Ref.alpha0[kf_surf],
        ts,
    )
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    end

    return
end

function free_convection_windspeed(self::SurfaceBase{SurfaceMoninObukhov}, GMV::GridMeanVariables)
    free_convection_windspeed(self, GMV, BaseCase())
    return
end

initialize(self::SurfaceBase{SurfaceMoninObukhovDry}) = nothing

function update(self::SurfaceBase{SurfaceMoninObukhovDry}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    kc_surf = kc_surface(self.Gr)
    kf_surf = kf_surface(self.Gr)
    zb = self.Gr.z_half[kc_surf]
    theta_rho_g = theta_rho_c(self.Ref.Pg, self.Tsurface, self.qsurface, self.qsurface)
    theta_rho_b = theta_rho_c(self.Ref.p0_half[kc_surf], GMV.T.values[kc_surf], self.qsurface, self.qsurface)
    lv = TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])

    h_star = t_to_thetali_c(param_set, self.Ref.Pg, self.Tsurface, self.qsurface, 0.0, 0.0)

    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(Ri, self.Gr.z_half[kc_surf], self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[kc_surf]) * self.Ref.rho0[kf_surf]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[kc_surf]) * self.Ref.rho0[kf_surf]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[kc_surf] - h_star) * self.Ref.rho0[kf_surf]
    self.rho_qtflux = -self.ch * self.windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.Ref.rho0[kf_surf]
    self.lhf = lv * 0.0

    ts = TD.PhaseEquil_pθq(param_set, self.Ref.p0[kf_surf], self.Tsurface, self.qsurface)
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.bflux = buoyancy_flux(param_set, self.shf, self.lhf, GMV.T.values[kc_surf], 0.0, self.Ref.alpha0[kf_surf], ts)
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    end

    return
end

free_convection_windspeed(self::SurfaceBase{SurfaceMoninObukhovDry}, GMV::GridMeanVariables) =
    free_convection_windspeed(self, GMV, BaseCase())

initialize(self::SurfaceBase{SurfaceSullivanPatton}) = nothing

function update(self::SurfaceBase{SurfaceSullivanPatton}, GMV::GridMeanVariables)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(self.Gr)
    kf_surf = kf_surface(self.Gr)
    zb = self.Gr.z_half[kc_surf]
    theta_rho_g = theta_rho_c(self.Ref.Pg, self.Tsurface, self.qsurface, self.qsurface)
    theta_rho_b = theta_rho_c(self.Ref.p0_half[kc_surf], GMV.T.values[kc_surf], self.qsurface, self.qsurface)
    ts = TD.PhaseEquil_pθq(param_set, self.Ref.p0[kc_surf], self.Tsurface, self.qsurface)
    lv = TD.latent_heat_vapor(param_set, GMV.T.values[kc_surf])
    T0 = self.Ref.p0_half[kc_surf] * self.Ref.alpha0_half[kc_surf] / Rd

    theta_flux = 0.24
    self.bflux = g * theta_flux * exner_c(self.Ref.p0_half[kc_surf]) / T0
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    h_star = t_to_thetali_c(param_set, self.Ref.Pg, self.Tsurface, self.qsurface, 0.0, 0.0)


    self.windspeed = sqrt(GMV.U.values[kc_surf] * GMV.U.values[kc_surf] + GMV.V.values[kc_surf] * GMV.V.values[kc_surf])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(Ri, self.Gr.z_half[kc_surf], self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[kc_surf]) * self.Ref.rho0[kf_surf]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[kc_surf]) * self.Ref.rho0[kf_surf]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[kc_surf] - h_star) * self.Ref.rho0[kf_surf]
    self.rho_qtflux = -self.ch * self.windspeed * (GMV.QT.values[kc_surf] - self.qsurface) * self.Ref.rho0[kf_surf]
    self.lhf = lv * self.rho_qtflux
    self.shf = TD.cp_m(ts) * self.rho_hflux

    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if abs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    end
    return
end

free_convection_windspeed(self::SurfaceBase{SurfaceSullivanPatton}, GMV::GridMeanVariables) =
    free_convection_windspeed(self, GMV, BaseCase())
