
initialize(self::SurfaceBase, ::BaseCase) = nothing
update(self::SurfaceBase, GMV::GridMeanVariables, ::BaseCase) = nothing

function free_convection_windspeed(self::SurfaceBase, GMV::GridMeanVariables, ::BaseCase)
    theta_rho = pyzeros(self.Gr.nzg)

    # Need to get theta_rho
    @inbounds for k in center_indicies(self.Gr)
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        theta_rho[k] = theta_rho_c(self.Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], qv)
    end
    zi = get_inversion(theta_rho, GMV.U.values, GMV.V.values, self.Gr, self.Ri_bulk_crit)
    wstar = get_wstar(self.bflux, zi) # yair here zi in TRMM should be adjusted
    self.windspeed = sqrt(self.windspeed * self.windspeed + (1.2 * wstar) * (1.2 * wstar))
    return
end

initialize(self::SurfaceBase{SurfaceNone}) = nothing

function update(self::SurfaceBase{SurfaceNone}, GMV::GridMeanVariables)
    # JH: assigning small fluxed so that simulation won"t crash when computing mixing length
    gw = self.Gr.gw
    self.windspeed = 0.0001
    self.zrough = 1e-4
    self.bflux = 1e-4
    self.ustar = compute_ustar(self.windspeed, self.bflux, self.zrough, self.Gr.z_half[gw])
    return
end
free_convection_windspeed(self::SurfaceBase{SurfaceNone}, GMV::GridMeanVariables) = nothing

initialize(self::SurfaceBase{SurfaceFixedFlux}) = nothing

function update(self::SurfaceBase{SurfaceFixedFlux}, GMV::GridMeanVariables)
    gw = self.Gr.gw
    rho_tflux = self.shf / (cpm_c(self.qsurface))

    self.windspeed = sqrt(GMV.U.values[gw] * GMV.U.values[gw] + GMV.V.values[gw] * GMV.V.values[gw])
    self.rho_qtflux = self.lhf / (latent_heat(self.Tsurface))

    self.rho_hflux = rho_tflux / exner_c(self.Ref.Pg)
    self.bflux = buoyancy_flux(self.shf, self.lhf, GMV.T.values[gw], GMV.QT.values[gw], self.Ref.alpha0[gw - 1])

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
                print("GMV.U.values[gw] ==>", GMV.U.values[gw])
                print("GMV.v.values[gw] ==>", GMV.V.values[gw])
                print("GMV.QT.values[gw] ==>", GMV.QT.values[gw])
                print("self.Ref.alpha0[gw-1] ==>", self.Ref.alpha0[gw - 1])
            end
        end

        self.ustar = compute_ustar(self.windspeed, self.bflux, self.zrough, self.Gr.z_half[gw])
    end

    self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    self.rho_uflux = -self.Ref.rho0[gw - 1] * self.ustar * self.ustar / self.windspeed * GMV.U.values[gw]
    self.rho_vflux = -self.Ref.rho0[gw - 1] * self.ustar * self.ustar / self.windspeed * GMV.V.values[gw]
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
    gw = self.Gr.gw
    windspeed = max(sqrt(GMV.U.values[gw] * GMV.U.values[gw] + GMV.V.values[gw] * GMV.V.values[gw]), 0.01)
    cp_ = cpm_c(GMV.QT.values[gw])
    lv = latent_heat(GMV.T.values[gw])

    self.rho_qtflux = -self.cq * windspeed * (GMV.QT.values[gw] - self.qsurface) * self.Ref.rho0[gw - 1]
    self.lhf = lv * self.rho_qtflux
    self.rho_hflux =
        -self.ch * windspeed * (GMV.H.values[gw] - self.Tsurface / exner_c(self.Ref.Pg)) * self.Ref.rho0[gw - 1]
    self.shf = cp_ * self.rho_hflux

    self.bflux = buoyancy_flux(self.shf, self.lhf, GMV.T.values[gw], GMV.QT.values[gw], self.Ref.alpha0[gw - 1])

    self.ustar = sqrt(self.cm) * windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if fabs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    end

    self.rho_uflux = -self.Ref.rho0[gw - 1] * self.ustar * self.ustar / windspeed * GMV.U.values[gw]
    self.rho_vflux = -self.Ref.rho0[gw - 1] * self.ustar * self.ustar / windspeed * GMV.V.values[gw]
    return
end

function free_convection_windspeed(self::SurfaceBase{SurfaceFixedCoeffs}, GMV::GridMeanVariables)
    free_convection_windspeed(self, GMV, BaseCase())
    return
end

initialize(self::SurfaceBase{SurfaceMoninObukhov}) = nothing

function update(self::SurfaceBase{SurfaceMoninObukhov}, GMV::GridMeanVariables)
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    gw = self.Gr.gw
    zb = self.Gr.z_half[gw]
    theta_rho_g = theta_rho_c(self.Ref.Pg, self.Tsurface, self.qsurface, self.qsurface)
    theta_rho_b = theta_rho_c(self.Ref.p0_half[gw], GMV.T.values[gw], self.qsurface, self.qsurface)
    lv = latent_heat(GMV.T.values[gw])

    h_star = t_to_thetali_c(self.Ref.Pg, self.Tsurface, self.qsurface, 0.0, 0.0)

    self.windspeed = sqrt(GMV.U.values[gw] * GMV.U.values[gw] + GMV.V.values[gw] * GMV.V.values[gw])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(Ri, self.Gr.z_half[gw], self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[gw]) * self.Ref.rho0[gw - 1]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[gw]) * self.Ref.rho0[gw - 1]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[gw] - h_star) * self.Ref.rho0[gw - 1]
    self.rho_qtflux = -self.ch * self.windspeed * (GMV.QT.values[gw] - self.qsurface) * self.Ref.rho0[gw - 1]
    self.lhf = lv * self.rho_qtflux

    self.shf = cpm_c(GMV.QT.values[gw]) * self.rho_hflux

    self.bflux = buoyancy_flux(self.shf, self.lhf, GMV.T.values[gw], GMV.QT.values[gw], self.Ref.alpha0[gw - 1])
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if fabs(self.bflux) < 1e-10
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
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    gw = self.Gr.gw
    zb = self.Gr.z_half[gw]
    theta_rho_g = theta_rho_c(self.Ref.Pg, self.Tsurface, self.qsurface, self.qsurface)
    theta_rho_b = theta_rho_c(self.Ref.p0_half[gw], GMV.T.values[gw], self.qsurface, self.qsurface)
    lv = latent_heat(GMV.T.values[gw])

    h_star = t_to_thetali_c(self.Ref.Pg, self.Tsurface, self.qsurface, 0.0, 0.0)

    self.windspeed = sqrt(GMV.U.values[gw] * GMV.U.values[gw] + GMV.V.values[gw] * GMV.V.values[gw])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(Ri, self.Gr.z_half[gw], self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[gw]) * self.Ref.rho0[gw - 1]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[gw]) * self.Ref.rho0[gw - 1]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[gw] - h_star) * self.Ref.rho0[gw - 1]
    self.rho_qtflux = -self.ch * self.windspeed * (GMV.QT.values[gw] - self.qsurface) * self.Ref.rho0[gw - 1]
    self.lhf = lv * 0.0

    self.shf = cpm_c(0.0) * self.rho_hflux

    self.bflux = buoyancy_flux(self.shf, self.lhf, GMV.T.values[gw], 0.0, self.Ref.alpha0[gw - 1])
    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if fabs(self.bflux) < 1e-10
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
    gw = self.Gr.gw
    zb = self.Gr.z_half[gw]
    theta_rho_g = theta_rho_c(self.Ref.Pg, self.Tsurface, self.qsurface, self.qsurface)
    theta_rho_b = theta_rho_c(self.Ref.p0_half[gw], GMV.T.values[gw], self.qsurface, self.qsurface)
    lv = latent_heat(GMV.T.values[gw])
    g = 9.81
    T0 = self.Ref.p0_half[gw] * self.Ref.alpha0_half[gw] / Rd

    theta_flux = 0.24
    self.bflux = g * theta_flux * exner_c(self.Ref.p0_half[gw]) / T0
    self.qsurface = qv_star_t(self.Ref.Pg, self.Tsurface)
    h_star = t_to_thetali_c(self.Ref.Pg, self.Tsurface, self.qsurface, 0.0, 0.0)


    self.windspeed = sqrt(GMV.U.values[gw] * GMV.U.values[gw] + GMV.V.values[gw] * GMV.V.values[gw])
    Nb2 = g / theta_rho_g * (theta_rho_b - theta_rho_g) / zb
    Ri = Nb2 * zb * zb / (self.windspeed * self.windspeed)

    #TODO: make sure pass by reference: &self.cm, &self.ch, &self.obukhov_length
    self.cm, self.ch, self.obukhov_length =
        exchange_coefficients_byun(Ri, self.Gr.z_half[gw], self.zrough, self.cm, self.ch, self.obukhov_length)
    self.rho_uflux = -self.cm * self.windspeed * (GMV.U.values[gw]) * self.Ref.rho0[gw - 1]
    self.rho_vflux = -self.cm * self.windspeed * (GMV.V.values[gw]) * self.Ref.rho0[gw - 1]

    self.rho_hflux = -self.ch * self.windspeed * (GMV.H.values[gw] - h_star) * self.Ref.rho0[gw - 1]
    self.rho_qtflux = -self.ch * self.windspeed * (GMV.QT.values[gw] - self.qsurface) * self.Ref.rho0[gw - 1]
    self.lhf = lv * self.rho_qtflux

    self.shf = cpm_c(GMV.QT.values[gw]) * self.rho_hflux

    self.ustar = sqrt(self.cm) * self.windspeed
    # CK--testing this--EDMF scheme checks greater or less than zero,
    if fabs(self.bflux) < 1e-10
        self.obukhov_length = 0.0
    else
        self.obukhov_length = -self.ustar * self.ustar * self.ustar / self.bflux / vkb
    end
    return
end

free_convection_windspeed(self::SurfaceBase{SurfaceSullivanPatton}, GMV::GridMeanVariables) =
    free_convection_windspeed(self, GMV, BaseCase())
