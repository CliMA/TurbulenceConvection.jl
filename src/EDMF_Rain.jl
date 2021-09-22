
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
    ref_state::ReferenceState,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics,
    TS::TimeStepping,
)

    write_profile(Stats, "qr_mean", self.QR.values)
    write_profile(Stats, "updraft_qr", self.Upd_QR.values)
    write_profile(Stats, "env_qr", self.Env_QR.values)
    write_profile(Stats, "rain_area", self.RainArea.values)
    write_profile(Stats, "updraft_rain_area", self.Upd_RainArea.values)
    write_profile(Stats, "env_rain_area", self.Env_RainArea.values)

    rain_diagnostics(self, ref_state, UpdThermo, EnvThermo)
    write_ts(Stats, "rwp_mean", self.mean_rwp)
    write_ts(Stats, "updraft_rwp", self.upd_rwp)
    write_ts(Stats, "env_rwp", self.env_rwp)
    write_ts(Stats, "cutoff_rain_rate", self.cutoff_rain_rate)
    #TODO - change to rain rate that depends on rain model choice
    return
end

function rain_diagnostics(
    self::RainVariables,
    ref_state::ReferenceState,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics,
)
    self.upd_rwp = 0.0
    self.env_rwp = 0.0
    self.mean_rwp = 0.0
    self.cutoff_rain_rate = 0.0
    grid = self.grid

    @inbounds for k in real_center_indices(grid)
        self.upd_rwp += ref_state.rho0_half[k] * self.Upd_QR.values[k] * self.Upd_RainArea.values[k] * grid.Δz
        self.env_rwp += ref_state.rho0_half[k] * self.Env_QR.values[k] * self.Env_RainArea.values[k] * grid.Δz
        self.mean_rwp += ref_state.rho0_half[k] * self.QR.values[k] * self.RainArea.values[k] * grid.Δz

        # rain rate from cutoff microphysics scheme defined as a total amount of removed water
        # per timestep per EDMF surface area [mm/h]
        if (self.rain_model == "cutoff")
            self.cutoff_rain_rate -=
                (EnvThermo.prec_source_qt[k] + UpdThermo.prec_source_qt_tot[k]) * ref_state.rho0_half[k] * grid.Δz /
                rho_cloud_liq *
                3.6 *
                1e6
        end
    end
    return
end

function sum_subdomains_rain(
    self::RainVariables,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics,
    TS::TimeStepping,
)
    @inbounds for k in real_center_indices(self.grid)
        self.QR.values[k] -= (EnvThermo.prec_source_qt[k] + UpdThermo.prec_source_qt_tot[k]) * TS.dt
        self.Upd_QR.values[k] -= UpdThermo.prec_source_qt_tot[k] * TS.dt
        self.Env_QR.values[k] -= EnvThermo.prec_source_qt[k] * TS.dt

        # TODO Assuming that updraft and environment rain area fractions are either 1 or 0.
        if self.QR.values[k] > 0.0
            self.RainArea.values[k] = 1
        end
        if self.Upd_QR.values[k] > 0.0
            self.Upd_RainArea.values[k] = 1
        end
        if self.Env_QR.values[k] > 0.0
            self.Env_RainArea.values[k] = 1
        end
    end
    return
end

function solve_rain_fall(
    self::RainPhysics,
    Rain::RainVariables,
    GMV::GridMeanVariables,
    TS::TimeStepping,
    QR::RainVariable,
    RainArea::RainVariable,
)
    param_set = parameter_set(GMV)
    grid = get_grid(GMV)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = 0.5

    term_vel = center_field(grid)
    term_vel_new = center_field(grid)
    ρ_0_c_field = self.ref_state.rho0_half

    # helper to calculate the rain velocity
    # TODO: assuming GMV.W = 0
    # TODO: verify translation
    @inbounds for k in real_center_indices(grid)
        term_vel[k] = CM1.terminal_velocity(param_set, rai_type, self.ref_state.rho0_half[k], QR.values[k])
    end

    # rain falling through the domain
    @inbounds for k in reverse(real_center_indices(grid))
        CFL_out = Δt / Δz * term_vel[k]

        if is_toa_center(grid, k)
            CFL_in = 0.0
        else
            CFL_in = Δt / Δz * term_vel[k + 1]
        end

        if max(CFL_in, CFL_out) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end
        ρ_0_c = self.ref_state.rho0_half[k]

        ρ_0_cut = ccut_downwind(ρ_0_c_field, grid, k)
        QR_cut = ccut_downwind(QR.values, grid, k)
        w_cut = ccut_downwind(term_vel, grid, k)
        ρQRw_cut = ρ_0_cut .* QR_cut .* w_cut
        ∇ρQRw = c∇_downwind(ρQRw_cut, grid, k; bottom = FreeBoundary(), top = SetValue(0))

        # TODO: incorporate area fraction into this equation (right now assume a = 1)
        QR.new[k] = (QR.values[k] * ρ_0_c + Δt * ∇ρQRw) / ρ_0_c
        QR.new[k] = max(QR.new[k], 0)
        if !(QR.new[k] ≈ 0)
            RainArea.new[k] = 1.0
        end

        term_vel_new[k] = CM1.terminal_velocity(param_set, rai_type, ρ_0_c, QR.new[k])
    end

    QR.values .= QR.new
    RainArea.values .= RainArea.new

    term_vel .= term_vel_new
    return
end

function solve_rain_evap(
    self::RainPhysics,
    Rain::RainVariables,
    GMV::GridMeanVariables,
    TS::TimeStepping,
    QR::RainVariable,
    RainArea::RainVariable,
)
    param_set = parameter_set(GMV)
    Δt = TS.dt
    flag_evaporate_all = false

    @inbounds for k in real_center_indices(self.grid)
        flag_evaporate_all = false

        q = TD.PhasePartition(GMV.QT.values[k], GMV.QL.values[k], 0.0)

        tmp_evap_rate =
            -CM1.evaporation_sublimation(
                param_set,
                rai_type,
                q,
                QR.values[k],
                self.ref_state.rho0_half[k],
                GMV.T.values[k],
            )

        if tmp_evap_rate * Δt > QR.values[k]
            flag_evaporate_all = true
            tmp_evap_rate = QR.values[k] / Δt
        end

        self.rain_evap_source_qt[k] = tmp_evap_rate * RainArea.values[k]

        # TODO add ice
        rain_source_to_thetal(
            param_set,
            self.ref_state.p0_half[k],
            GMV.T.values[k],
            GMV.QT.values[k],
            GMV.QL.values[k],
            -tmp_evap_rate,
        ) * RainArea.values[k]

        if flag_evaporate_all
            QR.values[k] = 0.0
            RainArea.values[k] = 0.0
        else
            # TODO: assuming that rain evaporation doesn"t change
            # rain area fraction
            # (should be changed for prognostic rain area fractions)
            # TODO  - new should be removed
            QR.values[k] -= tmp_evap_rate * Δt
        end
    end
    return
end
