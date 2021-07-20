function zero_tendencies(self::VariablePrognostic, Gr::Grid)
    @inbounds for k in xrange(Gr.nzg)
        self.tendencies[k] = 0.0
    end
    return
end

function set_bcs(self::VariablePrognostic, Gr::Grid)
    start_low = Gr.gw - 1
    start_high = Gr.nzg - Gr.gw - 1

    if self.bc == "sym"
        @inbounds for k in xrange(Gr.gw)
            self.values[start_high + k + 1] = self.values[start_high - k]
            self.values[start_low - k] = self.values[start_low + 1 + k]

            self.mf_update[start_high + k + 1] = self.mf_update[start_high - k]
            self.mf_update[start_low - k] = self.mf_update[start_low + 1 + k]

            self.new[start_high + k + 1] = self.new[start_high - k]
            self.new[start_low - k] = self.new[start_low + 1 + k]
        end
    else
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0

        self.mf_update[start_high] = 0.0
        self.mf_update[start_low] = 0.0

        self.new[start_high] = 0.0
        self.new[start_low] = 0.0

        @inbounds for k in xrange(1, Gr.gw)
            self.values[start_high + k] = -self.values[start_high - k]
            self.values[start_low - k] = -self.values[start_low + k]

            self.mf_update[start_high + k] = -self.mf_update[start_high - k]
            self.mf_update[start_low - k] = -self.mf_update[start_low + k]

            self.new[start_high + k] = -self.new[start_high - k]
            self.new[start_low - k] = -self.new[start_low + k]
        end
    end

    return
end

function set_bcs(self::VariableDiagnostic, Gr::Grid)
    start_low = Gr.gw - 1
    start_high = Gr.nzg - Gr.gw

    if self.bc == "sym"
        @inbounds for k in xrange(Gr.gw)
            self.values[start_high + k] = self.values[start_high - 1]
            self.values[start_low - k] = self.values[start_low + 1]
        end

    else
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0
        @inbounds for k in xrange(1, Gr.gw)
            self.values[start_high + k] = 0.0  #-self.values[start_high - k ]
            self.values[start_low - k] = 0.0 #-self.values[start_low + k ]
        end
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
    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg - self.Gr.gw)
        self.U.values[k] += self.U.tendencies[k] * TS.dt
        self.V.values[k] += self.V.tendencies[k] * TS.dt
        self.H.values[k] += self.H.tendencies[k] * TS.dt
        self.QT.values[k] += self.QT.tendencies[k] * TS.dt
    end

    set_bcs(self.U, self.Gr)
    set_bcs(self.V, self.Gr)
    set_bcs(self.H, self.Gr)
    set_bcs(self.QT, self.Gr)

    if self.calc_tke
        set_bcs(self.TKE, self.Gr)
    end

    if self.calc_scalar_var
        set_bcs(self.QTvar, self.Gr)
        set_bcs(self.Hvar, self.Gr)
        set_bcs(self.HQTcov, self.Gr)
    end

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
    if self.calc_tke
        add_profile(Stats, "tke_mean")
    end
    if self.calc_scalar_var
        add_profile(Stats, "Hvar_mean")
        add_profile(Stats, "QTvar_mean")
        add_profile(Stats, "HQTcov_mean")

        add_profile(Stats, "W_third_m")
        add_profile(Stats, "H_third_m")
        add_profile(Stats, "QT_third_m")
    end

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
    if self.calc_tke
        write_profile(Stats, "tke_mean", self.TKE.values[cinterior])
        write_profile(Stats, "W_third_m", self.W_third_m.values[cinterior])
    end
    if self.calc_scalar_var
        write_profile(Stats, "Hvar_mean", self.Hvar.values[cinterior])
        write_profile(Stats, "QTvar_mean", self.QTvar.values[cinterior])
        write_profile(Stats, "HQTcov_mean", self.HQTcov.values[cinterior])

        write_profile(Stats, "H_third_m", self.H_third_m.values[cinterior])
        write_profile(Stats, "QT_third_m", self.QT_third_m.values[cinterior])
    end

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
    self.cloud_base = self.Gr.z_half[self.Gr.nzg - self.Gr.gw - 1]
    self.cloud_top = 0.0

    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg - self.Gr.gw)
        self.lwp += self.Ref.rho0_half[k] * self.QL.values[k] * self.Gr.dz

        if self.QL.values[k] > 1e-8
            self.cloud_base = fmin(self.cloud_base, self.Gr.z_half[k])
            self.cloud_top = fmax(self.cloud_top, self.Gr.z_half[k])
        end
    end
    return
end

function satadjust(self::GridMeanVariables)
    sa = eos_struct()
    @inbounds for k in xrange(self.Gr.nzg)
        h = self.H.values[k]
        qt = self.QT.values[k]
        p0 = self.Ref.p0_half[k]
        sa = eos(self.t_to_prog_fp, self.prog_to_t_fp, p0, qt, h)
        self.QL.values[k] = sa.ql
        self.T.values[k] = sa.T
        qv = qt - sa.ql
        self.THL.values[k] = t_to_thetali_c(p0, sa.T, qt, sa.ql, 0.0)
        rho = rho_c(p0, sa.T, qt, qv)
        self.B.values[k] = buoyancy_c(self.Ref.rho0_half[k], rho)
        self.RH.values[k] = relative_humidity_c(self.Ref.p0_half[k], qt, qt - qv, 0.0, self.T.values[k])
    end
    return
end
