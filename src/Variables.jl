function update(self::GridMeanVariables, TS::TimeStepping)
    grid = self.grid
    @inbounds for k in real_center_indices(grid)
        self.U.values[k] += self.U.tendencies[k] * TS.dt
        self.V.values[k] += self.V.tendencies[k] * TS.dt
        self.H.values[k] += self.H.tendencies[k] * TS.dt
        self.QT.values[k] += self.QT.tendencies[k] * TS.dt
    end
    return
end

function initialize_io(self::GridMeanVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "u_mean")
    add_profile(Stats, "v_mean")
    add_profile(Stats, "qt_mean")
    add_profile(Stats, "RH_mean")
    add_profile(Stats, "thetal_mean")
    add_profile(Stats, "temperature_mean")
    add_profile(Stats, "buoyancy_mean")
    add_profile(Stats, "ql_mean")
    add_profile(Stats, "tke_mean")
    add_profile(Stats, "Hvar_mean")
    add_profile(Stats, "QTvar_mean")
    add_profile(Stats, "HQTcov_mean")

    add_profile(Stats, "W_third_m")
    add_profile(Stats, "H_third_m")
    add_profile(Stats, "QT_third_m")

    add_profile(Stats, "cloud_fraction_mean")

    add_ts(Stats, "lwp_mean")
    add_ts(Stats, "cloud_base_mean")
    add_ts(Stats, "cloud_top_mean")
    add_ts(Stats, "cloud_cover_mean")
    return
end

function io(self::GridMeanVariables, Stats::NetCDFIO_Stats)
    write_profile(Stats, "u_mean", self.U.values)
    write_profile(Stats, "v_mean", self.V.values)
    write_profile(Stats, "qt_mean", self.QT.values)
    write_profile(Stats, "ql_mean", self.QL.values)
    write_profile(Stats, "temperature_mean", self.T.values)
    write_profile(Stats, "RH_mean", self.RH.values)
    write_profile(Stats, "buoyancy_mean", self.B.values)
    write_profile(Stats, "thetal_mean", self.H.values)
    write_profile(Stats, "tke_mean", self.TKE.values)
    write_profile(Stats, "W_third_m", self.W_third_m.values)
    write_profile(Stats, "Hvar_mean", self.Hvar.values)
    write_profile(Stats, "QTvar_mean", self.QTvar.values)
    write_profile(Stats, "HQTcov_mean", self.HQTcov.values)

    write_profile(Stats, "H_third_m", self.H_third_m.values)
    write_profile(Stats, "QT_third_m", self.QT_third_m.values)

    write_profile(Stats, "cloud_fraction_mean", self.cloud_fraction.values)
    write_ts(Stats, "cloud_cover_mean", self.cloud_cover)

    mean_cloud_diagnostics(self)
    write_ts(Stats, "lwp_mean", self.lwp)
    write_ts(Stats, "cloud_base_mean", self.cloud_base)
    write_ts(Stats, "cloud_top_mean", self.cloud_top)
    return
end

function mean_cloud_diagnostics(self)
    self.lwp = 0.0
    kc_toa = kc_top_of_atmos(self.grid)
    self.cloud_base = self.grid.zc[kc_toa]
    self.cloud_top = 0.0

    @inbounds for k in real_center_indices(self.grid)
        self.lwp += self.ref_state.rho0_half[k] * self.QL.values[k] * self.grid.Δz

        if self.QL.values[k] > 1e-8
            self.cloud_base = min(self.cloud_base, self.grid.zc[k])
            self.cloud_top = max(self.cloud_top, self.grid.zc[k])
        end
    end
    return
end

function satadjust(self::GridMeanVariables)
    sa = eos_struct()
    param_set = parameter_set(self)
    @inbounds for k in real_center_indices(self.grid)
        h = self.H.values[k]
        qt = self.QT.values[k]
        p0 = self.ref_state.p0_half[k]
        sa = eos(param_set, p0, qt, h)
        self.QL.values[k] = sa.ql
        self.T.values[k] = sa.T
        qv = qt - sa.ql
        rho = rho_c(p0, sa.T, qt, qv)
        self.B.values[k] = buoyancy_c(param_set, self.ref_state.rho0_half[k], rho)
        ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[k], h, qt)
        self.RH.values[k] = TD.relative_humidity(ts)
    end
    return
end
