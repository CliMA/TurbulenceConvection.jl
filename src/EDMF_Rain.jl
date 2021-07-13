
function set_bcs(self::RainVariable, Gr::Grid)
    @inbounds for k in xrange(Gr.gw)
        self.values[Gr.nzg-Gr.gw+k] = self.values[Gr.nzg-Gr.gw - 1 - k]
        self.values[Gr.gw-1-k]      = self.values[Gr.gw + k]
    end
    return
end

function RainVariables(namelist, Gr::Grid)
    nzg = Gr.nzg

    QR           = RainVariable(nzg, "qr_mean",       "kg/kg")
    # temporary variables for diagnostics to know where the rain is coming from
    Upd_QR       = RainVariable(nzg, "upd_qr",        "kg/kg")
    Env_QR       = RainVariable(nzg, "env_qr",        "kg/kg")
    # in the future we could test prognostic equations for stratiform and updraft rain
    RainArea     = RainVariable(nzg, "rain_area",     "rain_area_fraction [-]" )
    Upd_RainArea = RainVariable(nzg, "upd_rain_area", "updraft_rain_area_fraction [-]" )
    Env_RainArea = RainVariable(nzg, "env_rain_area", "environment_rain_area_fraction [-]" )

    mean_rwp = 0.
    upd_rwp = 0.
    env_rwp = 0.

    rain_model = try
        string(namelist["microphysics"]["rain_model"])
    catch
        println("EDMF_Rain: defaulting to no rain")
        "None"
    end

    if !(rain_model in ["None", "cutoff", "clima_1m"])
        error("rain model not recognized")
    end

    return RainVariables(;rain_model,
        mean_rwp,
        env_rwp,
        upd_rwp,
        Gr,
        QR,
        RainArea,
        Upd_QR,
        Upd_RainArea,
        Env_QR,
        Env_RainArea)
end

function initialize_io(self::RainVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "qr_mean")
    add_profile(Stats, "updraft_qr")
    add_profile(Stats, "env_qr")
    add_profile(Stats, "rain_area")
    add_profile(Stats, "updraft_rain_area")
    add_profile(Stats, "env_rain_area")
    add_ts(Stats, "rwp_mean")
    add_ts(Stats, "updraft_rwp")
    add_ts(Stats, "env_rwp")
    add_ts(Stats, "cutoff_rain_rate")
    return
end

function io(
    self::RainVariables,
    Stats::NetCDFIO_Stats,
    Ref::ReferenceState,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics,
    TS::TimeStepping)

    cinterior = self.Gr.cinterior
    finterior = self.Gr.finterior

    write_profile(Stats, "qr_mean",           self.QR.values[cinterior])
    write_profile(Stats, "updraft_qr",        self.Upd_QR.values[cinterior])
    write_profile(Stats, "env_qr",            self.Env_QR.values[cinterior])
    write_profile(Stats, "rain_area",         self.RainArea.values[cinterior])
    write_profile(Stats, "updraft_rain_area", self.Upd_RainArea.values[cinterior])
    write_profile(Stats, "env_rain_area",     self.Env_RainArea.values[cinterior])

    rain_diagnostics(self, Ref, UpdThermo, EnvThermo, TS)
    write_ts(Stats, "rwp_mean", self.mean_rwp)
    write_ts(Stats, "updraft_rwp", self.upd_rwp)
    write_ts(Stats, "env_rwp", self.env_rwp)
    write_ts(Stats, "cutoff_rain_rate", self.cutoff_rain_rate)
    #TODO - change to rain rate that depends on rain model choice
    return
end

function rain_diagnostics(
    self::RainVariables,
    Ref::ReferenceState,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics,
    TS::TimeStepping
    )
    self.upd_rwp  = 0.
    self.env_rwp  = 0.
    self.mean_rwp = 0.
    self.cutoff_rain_rate = 0.

  @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        self.upd_rwp  += Ref.rho0_half[k] * self.Upd_QR.values[k] * self.Upd_RainArea.values[k] * self.Gr.dz
        self.env_rwp  += Ref.rho0_half[k] * self.Env_QR.values[k] * self.Env_RainArea.values[k] * self.Gr.dz
        self.mean_rwp += Ref.rho0_half[k] * self.QR.values[k]     * self.RainArea.values[k]     * self.Gr.dz

        # rain rate from cutoff microphysics scheme defined as a total amount of removed water
        # per timestep per EDMF surface area [mm/h]
        if(self.rain_model == "cutoff")
            self.cutoff_rain_rate -= (EnvThermo.prec_source_qt[k] + UpdThermo.prec_source_qt_tot[k]) *
                                     Ref.rho0_half[k] * self.Gr.dz / TS.dt / rho_cloud_liq *
                                     3.6 * 1e6
        end
    end
    return
end

function sum_subdomains_rain(
    self::RainVariables,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics)
  @inbounds for k in xrange(self.Gr.nzg)

        self.QR.values[k]       -= (EnvThermo.prec_source_qt[k] + UpdThermo.prec_source_qt_tot[k])
        self.Upd_QR.values[k]   -= UpdThermo.prec_source_qt_tot[k]
        self.Env_QR.values[k]   -= EnvThermo.prec_source_qt[k]

        # TODO Assuming that updraft and environment rain area fractions are either 1 or 0.
        if self.QR.values[k] > 0.
            self.RainArea.values[k] = 1
        end
        if self.Upd_QR.values[k] > 0.
            self.Upd_RainArea.values[k] = 1
        end
        if self.Env_QR.values[k] > 0.
            self.Env_RainArea.values[k] = 1
        end
    end
    return
end

function solve_rain_fall(
    self::RainPhysics,
    GMV::GridMeanVariables,
    TS::TimeStepping,
    QR::RainVariable,
    RainArea::RainVariable
)
    gw  = self.Gr.gw
    nzg = self.Gr.nzg
    dz = self.Gr.dz
    dt_model = TS.dt
    CFL_limit = 0.5

    term_vel     = pyzeros(nzg)
    term_vel_new = pyzeros(nzg)

    t_elapsed = 0.

    # helper to calculate the rain velocity
    # TODO: assuming GMV.W = 0
    # TODO: verify translation
  @inbounds for k in revxrange(nzg - gw - 1, gw, -1)
        term_vel[k] = terminal_velocity(
                          QR.values[k],
                          self.Ref.rho0_half[k]
                       )
    end

    # calculate the allowed timestep (CFL_limit >= v dt / dz)
    # TODO: report bug: dt_rain has no value in else-case
    if !(maximum(term_vel) ≈ 0)
        dt_rain = min(dt_model, CFL_limit * self.Gr.dz / maximum(term_vel))
    else
        dt_rain = dt_model
    end

    # rain falling through the domain
    while t_elapsed < dt_model
        # TODO: verify translation
      @inbounds for k in revxrange(nzg - gw - 1, gw, -1)

            CFL_out = dt_rain / dz * term_vel[k]

            if k == (nzg - gw - 1)
                CFL_in = 0.
            else
                CFL_in = dt_rain / dz * term_vel[k+1]
            end

            rho_frac  = self.Ref.rho0_half[k+1] / self.Ref.rho0_half[k]
            area_frac = 1. # RainArea.values[k] / RainArea.new[k]

            QR.new[k] = (QR.values[k]   * (1 - CFL_out) +
                         QR.values[k+1] * CFL_in * rho_frac) * area_frac
            if QR.new[k] != 0.
                RainArea.new[k] = 1.
            end

            term_vel_new[k] = terminal_velocity(
                                  QR.new[k],
                                  self.Ref.rho0_half[k]
                              )
        end

        t_elapsed += dt_rain

        QR.values .= QR.new
        RainArea.values .= RainArea.new

        term_vel .= term_vel_new

        if maximum(abs.(term_vel)) > eps(Float64)
            dt_rain = min(dt_model - t_elapsed,
                     CFL_limit * self.Gr.dz / maximum(term_vel)
                    )
        else
            dt_rain = dt_model - t_elapsed
        end
    end

    return
end

function solve_rain_evap(
    self::RainPhysics,
    GMV::GridMeanVariables,
    TS::TimeStepping,
    QR::RainVariable,
    RainArea::RainVariable
)
    gw  = self.Gr.gw
    nzg = self.Gr.nzg

    dz = self.Gr.dz
    dt_model = TS.dt

    flag_evaporate_all = false

  @inbounds for k in xrange(gw, nzg - gw)

        flag_evaporate_all = false

        tmp_evap = max(0, conv_q_rai_to_q_vap(QR.values[k],
                                              GMV.QT.values[k],
                                              GMV.QL.values[k],
                                              GMV.T.values[k],
                                              self.Ref.p0_half[k],
                                              self.Ref.rho0[k]
                                            ) * dt_model)

        if tmp_evap > QR.values[k]
            flag_evaporate_all = true
            tmp_evap = QR.values[k]
        end

        self.rain_evap_source_qt[k] = tmp_evap * RainArea.values[k]

        self.rain_evap_source_h[k]  = rain_source_to_thetal(
            self.Ref.p0[k],
            GMV.T.values[k],
            - tmp_evap
        ) * RainArea.values[k]

        if flag_evaporate_all
            QR.values[k] = 0.
            RainArea.values[k] = 0.
        else
            # TODO: assuming that rain evaporation doesn"t change
            # rain area fraction
            # (should be changed for prognostic rain area fractions)
            QR.values[k] -= tmp_evap
        end
    end
    return
end
