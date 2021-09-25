
function initialize_io(self::RainVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "qr_mean")
    add_ts(Stats, "rwp_mean")
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

    rain_diagnostics(self, ref_state, UpdThermo, EnvThermo)
    write_ts(Stats, "rwp_mean", self.mean_rwp)

    #TODO - change to rain rate that depends on rain model choice
    write_ts(Stats, "cutoff_rain_rate", self.cutoff_rain_rate)
    return
end

function rain_diagnostics(
    self::RainVariables,
    ref_state::ReferenceState,
    UpdThermo::UpdraftThermodynamics,
    EnvThermo::EnvironmentThermodynamics,
)
    self.mean_rwp = 0.0
    self.cutoff_rain_rate = 0.0
    grid = self.grid

    @inbounds for k in real_center_indices(grid)
        self.mean_rwp += ref_state.rho0_half[k] * self.QR.values[k] * grid.Δz

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
    end
    return
end

function solve_rain_fall(self::RainPhysics, GMV::GridMeanVariables, TS::TimeStepping, QR::RainVariable)
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
        term_vel[k] = CM1.terminal_velocity(param_set, rain_type, self.ref_state.rho0_half[k], QR.values[k])
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

        term_vel_new[k] = CM1.terminal_velocity(param_set, rain_type, ρ_0_c, QR.new[k])
    end

    QR.values .= QR.new

    term_vel .= term_vel_new
    return
end

function solve_rain_evap(self::RainPhysics, GMV::GridMeanVariables, TS::TimeStepping, QR::RainVariable)
    param_set = parameter_set(GMV)
    Δt = TS.dt
    flag_evaporate_all = false

    @inbounds for k in real_center_indices(self.grid)
        flag_evaporate_all = false

        q = TD.PhasePartition(GMV.QT.values[k], GMV.QL.values[k], 0.0)

        tmp_evap_rate =
            -CM1.evaporation_sublimation(
                param_set,
                rain_type,
                q,
                QR.values[k],
                self.ref_state.rho0_half[k],
                GMV.T.values[k],
            )

        if tmp_evap_rate * Δt > QR.values[k]
            flag_evaporate_all = true
            tmp_evap_rate = QR.values[k] / Δt
        end

        self.rain_evap_source_qt[k] = tmp_evap_rate

        # TODO add ice
        rain_source_to_thetal(
            param_set,
            self.ref_state.p0_half[k],
            GMV.T.values[k],
            GMV.QT.values[k],
            GMV.QL.values[k],
            0.0,
            -tmp_evap_rate,
        )

        if flag_evaporate_all
            QR.values[k] = 0.0
        else
            QR.values[k] -= tmp_evap_rate * Δt
        end

    end
    return
end
