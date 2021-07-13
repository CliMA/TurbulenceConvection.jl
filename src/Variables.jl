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
            self.values[start_high+k+1] = self.values[start_high-k]
            self.values[start_low-k] = self.values[start_low+1+k]

            self.mf_update[start_high+k+1] = self.mf_update[start_high-k]
            self.mf_update[start_low-k] = self.mf_update[start_low+1+k]

            self.new[start_high+k+1] = self.new[start_high-k]
            self.new[start_low-k] = self.new[start_low+1+k]
        end
    else
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0

        self.mf_update[start_high] = 0.0
        self.mf_update[start_low] = 0.0

        self.new[start_high] = 0.0
        self.new[start_low] = 0.0

      @inbounds for k in xrange(1,Gr.gw)
            self.values[start_high+ k] = -self.values[start_high - k ]
            self.values[start_low- k] = -self.values[start_low + k  ]

            self.mf_update[start_high+ k] = -self.mf_update[start_high - k ]
            self.mf_update[start_low- k] = -self.mf_update[start_low + k  ]

            self.new[start_high+ k] = -self.new[start_high - k ]
            self.new[start_low- k] = -self.new[start_low + k  ]
        end
    end

    return
end

function set_bcs(self::VariableDiagnostic,Gr::Grid)
    start_low = Gr.gw - 1
    start_high = Gr.nzg - Gr.gw

    if self.bc == "sym"
      @inbounds for k in xrange(Gr.gw)
            self.values[start_high + k] = self.values[start_high  - 1]
            self.values[start_low - k] = self.values[start_low + 1]
        end

    else
        self.values[start_high] = 0.0
        self.values[start_low] = 0.0
      @inbounds for k in xrange(1,Gr.gw)
            self.values[start_high+ k] = 0.0  #-self.values[start_high - k ]
            self.values[start_low- k] = 0.0 #-self.values[start_low + k ]
        end
    end

    return
end

function GridMeanVariables(namelist, Gr::Grid, Ref::ReferenceState)
    lwp = 0.
    cloud_base   = 0.
    cloud_top    = 0.
    cloud_cover  = 0.

    U = VariablePrognostic(Gr.nzg, "half", "velocity", "sym","u", "m/s" )
    V = VariablePrognostic(Gr.nzg, "half", "velocity","sym", "v", "m/s" )
    # Just leave this zero for now!
    W = VariablePrognostic(Gr.nzg, "full", "velocity","asym", "v", "m/s" )

    # Create thermodynamic variables
    QT = VariablePrognostic(Gr.nzg, "half", "scalar","sym", "qt", "kg/kg")
    RH = VariablePrognostic(Gr.nzg, "half", "scalar","sym", "RH", "%")

    if namelist["thermodynamics"]["thermal_variable"] == "entropy"
        H = VariablePrognostic(Gr.nzg, "half", "scalar", "sym","s", "J/kg/K" )
        t_to_prog_fp = t_to_entropy_c
        prog_to_t_fp = eos_first_guess_entropy
    elseif namelist["thermodynamics"]["thermal_variable"] == "thetal"
        H = VariablePrognostic(Gr.nzg, "half", "scalar", "sym","thetal", "K")
        t_to_prog_fp = t_to_thetali_c
        prog_to_t_fp = eos_first_guess_thetal
    else
        error("Did not recognize thermal variable " + namelist["thermodynamics"]["thermal_variable"])
    end

    # Diagnostic Variables--same class as the prognostic variables, but we append to diagnostics list
    QL  = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "ql",              "kg/kg")
    T   = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "temperature",     "K")
    B   = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "buoyancy",        "m^2/s^3")
    THL = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "thetal",          "K")

    cloud_fraction  = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "cloud fraction", "-")

    # TKE   TODO   repeated from EDMF_Environment.pyx logic
    if  namelist["turbulence"]["scheme"] == "EDMF_PrognosticTKE"
        calc_tke = true
    else
        calc_tke = false
    end
    try
        calc_tke = namelist["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"]
    catch
    end

    calc_scalar_var = try
        namelist["turbulence"]["EDMF_PrognosticTKE"]["calc_scalar_var"]
    catch
        false
    end

    EnvThermo_scheme = try
        string(namelist["thermodynamics"]["sgs"])
    catch
        "mean"
    end

    #Now add the 2nd moment variables
    if calc_tke
        TKE = VariableDiagnostic(Gr.nzg, "half", "scalar","sym", "tke","m^2/s^2" )
        W_third_m = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "W_third_m", "m^3/s^3")
    end

    if calc_scalar_var
        QTvar = VariableDiagnostic(Gr.nzg, "half", "scalar","sym", "qt_var","kg^2/kg^2" )
        QT_third_m = VariableDiagnostic(Gr.nzg, "half", "scalar","sym", "qt_third_m","kg^3/kg^3" )
        if namelist["thermodynamics"]["thermal_variable"] == "entropy"
            Hvar = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "s_var", "(J/kg/K)^2")
            H_third_m = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "s__third_m", "-")
            HQTcov = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym" ,"s_qt_covar", "(J/kg/K)(kg/kg)" )
        elseif namelist["thermodynamics"]["thermal_variable"] == "thetal"
            Hvar = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym" ,"thetal_var", "K^2")
            H_third_m = VariableDiagnostic(Gr.nzg, "half", "scalar", "sym", "thetal_third_m", "-")
            HQTcov = VariableDiagnostic(Gr.nzg, "half", "scalar","sym" ,"thetal_qt_covar", "K(kg/kg)" )
        end
    end

    return GridMeanVariables(;Gr,
        Ref,
        lwp,
        cloud_base,
        cloud_top,
        cloud_cover,
        U,
        V,
        W,
        QT,
        RH,
        H,
        t_to_prog_fp,
        prog_to_t_fp,
        QL,
        T,
        B,
        THL,
        cloud_fraction,
        calc_tke,
        calc_scalar_var,
        EnvThermo_scheme,
        TKE,
        W_third_m,
        QTvar,
        QT_third_m,
        Hvar,
        H_third_m,
        HQTcov)

end

function zero_tendencies(self::GridMeanVariables)
    zero_tendencies(self.U, self.Gr)
    zero_tendencies(self.V, self.Gr)
    zero_tendencies(self.QT, self.Gr)
    zero_tendencies(self.H, self.Gr)
    return
end

function update(self::GridMeanVariables, TS::TimeStepping)
  @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        self.U.values[k]  +=  self.U.tendencies[k] * TS.dt
        self.V.values[k]  +=  self.V.tendencies[k] * TS.dt
        self.H.values[k]  +=  self.H.tendencies[k] * TS.dt
        self.QT.values[k] +=  self.QT.tendencies[k] * TS.dt
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
    if self.H.name == "s"
        add_profile(Stats, "s_mean")
        add_profile(Stats, "thetal_mean")
    elseif self.H.name == "thetal"
        add_profile(Stats, "thetal_mean")
    end

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
    write_profile(Stats, "v_mean",self.V.values[cinterior])
    write_profile(Stats, "qt_mean",self.QT.values[cinterior])
    write_profile(Stats, "ql_mean",self.QL.values[cinterior])
    write_profile(Stats, "temperature_mean",self.T.values[cinterior])
    write_profile(Stats, "RH_mean",self.RH.values[cinterior])
    write_profile(Stats, "buoyancy_mean",self.B.values[cinterior])
    if self.H.name == "s"
        write_profile(Stats, "s_mean",self.H.values[cinterior])
        write_profile(Stats, "thetal_mean",self.THL.values[cinterior])
    elseif self.H.name == "thetal"
        write_profile(Stats, "thetal_mean",self.H.values[cinterior])
    end
    if self.calc_tke
        write_profile(Stats, "tke_mean",self.TKE.values[cinterior])
        write_profile(Stats, "W_third_m",self.W_third_m.values[cinterior])
    end
    if self.calc_scalar_var
        write_profile(Stats, "Hvar_mean",self.Hvar.values[cinterior])
        write_profile(Stats, "QTvar_mean",self.QTvar.values[cinterior])
        write_profile(Stats, "HQTcov_mean",self.HQTcov.values[cinterior])

        write_profile(Stats, "H_third_m",self.H_third_m.values[cinterior])
        write_profile(Stats, "QT_third_m",self.QT_third_m.values[cinterior])
    end

    write_profile(Stats, "cloud_fraction_mean",self.cloud_fraction.values[cinterior])
    write_ts(Stats, "cloud_cover_mean", self.cloud_cover)

    mean_cloud_diagnostics(self)
    write_ts(Stats, "lwp_mean", self.lwp)
    write_ts(Stats, "cloud_base_mean",  self.cloud_base)
    write_ts(Stats, "cloud_top_mean",   self.cloud_top)
    return
end

function mean_cloud_diagnostics(self)
    self.lwp = 0.
    self.cloud_base   = self.Gr.z_half[self.Gr.nzg - self.Gr.gw - 1]
    self.cloud_top    = 0.

  @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        self.lwp += self.Ref.rho0_half[k] * self.QL.values[k] * self.Gr.dz

        if self.QL.values[k] > 1e-8
            self.cloud_base  = fmin(self.cloud_base,  self.Gr.z_half[k])
            self.cloud_top   = fmax(self.cloud_top,   self.Gr.z_half[k])
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
        sa = eos(self.t_to_prog_fp,self.prog_to_t_fp, p0, qt, h )
        self.QL.values[k] = sa.ql
        self.T.values[k] = sa.T
        qv = qt - sa.ql
        self.THL.values[k] = t_to_thetali_c(p0, sa.T, qt, sa.ql,0.0)
        rho = rho_c(p0, sa.T, qt, qv)
        self.B.values[k] = buoyancy_c(self.Ref.rho0_half[k], rho)
        self.RH.values[k] = relative_humidity_c(self.Ref.p0_half[k], qt, qt-qv, 0.0, self.T.values[k])
    end
    return
end
