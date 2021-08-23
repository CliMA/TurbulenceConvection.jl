function zero_tendencies(self::VariablePrognostic, Gr::Grid)
    @inbounds for k in center_indicies(Gr)
        self.tendencies[k] = 0.0
    end
    return
end

function zero_tendencies(self::GridMeanVariables)
    zero_tendencies(self.U, self.Gr)
    zero_tendencies(self.V, self.Gr)
    zero_tendencies(self.QT, self.Gr)
    zero_tendencies(self.H, self.Gr)
    return
end

function update(self::GridMeanVariables, TS::TimeStepping)
    grid = self.Gr
    @inbounds for k in center_indicies(grid)
        self.U.values[k] += self.U.tendencies[k] * TS.dt
        self.V.values[k] += self.V.tendencies[k] * TS.dt
        self.H.values[k] += self.H.tendencies[k] * TS.dt
        self.QT.values[k] += self.QT.tendencies[k] * TS.dt
    end

    set_bcs(self.U, grid)
    set_bcs(self.V, grid)
    set_bcs(self.H, grid)
    set_bcs(self.QT, grid)
    set_bcs(self.TKE, grid)

    set_bcs(self.QTvar, grid)
    set_bcs(self.Hvar, grid)
    set_bcs(self.HQTcov, grid)

    zero_tendencies(self)
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
    cinterior = self.Gr.cinterior
    finterior = self.Gr.finterior

    write_profile(Stats, "u_mean", self.U.values[cinterior])
    write_profile(Stats, "v_mean", self.V.values[cinterior])
    write_profile(Stats, "qt_mean", self.QT.values[cinterior])
    write_profile(Stats, "ql_mean", self.QL.values[cinterior])
    write_profile(Stats, "temperature_mean", self.T.values[cinterior])
    write_profile(Stats, "RH_mean", self.RH.values[cinterior])
    write_profile(Stats, "buoyancy_mean", self.B.values[cinterior])
    write_profile(Stats, "thetal_mean", self.H.values[cinterior])
    write_profile(Stats, "tke_mean", self.TKE.values[cinterior])
    write_profile(Stats, "W_third_m", self.W_third_m.values[cinterior])
    write_profile(Stats, "Hvar_mean", self.Hvar.values[cinterior])
    write_profile(Stats, "QTvar_mean", self.QTvar.values[cinterior])
    write_profile(Stats, "HQTcov_mean", self.HQTcov.values[cinterior])

    write_profile(Stats, "H_third_m", self.H_third_m.values[cinterior])
    write_profile(Stats, "QT_third_m", self.QT_third_m.values[cinterior])

    write_profile(Stats, "cloud_fraction_mean", self.cloud_fraction.values[cinterior])
    write_ts(Stats, "cloud_cover_mean", self.cloud_cover)

    mean_cloud_diagnostics(self)
    write_ts(Stats, "lwp_mean", self.lwp)
    write_ts(Stats, "cloud_base_mean", self.cloud_base)
    write_ts(Stats, "cloud_top_mean", self.cloud_top)
    return
end

function mean_cloud_diagnostics(self)
    self.lwp = 0.0
    kc_toa = kc_top_of_atmos(self.Gr)
    self.cloud_base = self.Gr.z_half[kc_toa]
    self.cloud_top = 0.0

    @inbounds for k in real_center_indicies(self.Gr)
        self.lwp += self.Ref.rho0_half[k] * self.QL.values[k] * self.Gr.dz

        if self.QL.values[k] > 1e-8
            self.cloud_base = min(self.cloud_base, self.Gr.z_half[k])
            self.cloud_top = max(self.cloud_top, self.Gr.z_half[k])
        end
    end
    return
end

function satadjust(self::GridMeanVariables)

    param_set = parameter_set(self)

    @inbounds for k in center_indicies(self.Gr)
        p0 = self.Ref.p0_half[k]

        h = self.H.values[k]
        qt = self.QT.values[k]
        ts = TD.PhaseEquil_pθq(param_set, p0, h, qt)

        self.QL.values[k] = TD.liquid_specific_humidity(ts)
        self.T.values[k] = TD.air_temperature(ts)
        qv = TD.vapor_specific_humidity(ts)
        rho = TD.air_density(ts)

        self.B.values[k] = buoyancy_c(param_set, self.Ref.rho0_half[k], rho)
        self.RH.values[k] = TD.relative_humidity(ts)
    end
    return
end
