function initialize_io(self::EnvironmentVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "env_w")
    add_profile(Stats, "env_qt")
    add_profile(Stats, "env_ql")
    add_profile(Stats, "env_area")
    add_profile(Stats, "env_temperature")
    add_profile(Stats, "env_RH")
    add_profile(Stats, "env_thetal")
    add_profile(Stats, "env_tke")
    add_profile(Stats, "env_Hvar")
    add_profile(Stats, "env_QTvar")
    add_profile(Stats, "env_HQTcov")

    add_profile(Stats, "env_cloud_fraction")

    add_ts(Stats, "env_cloud_base")
    add_ts(Stats, "env_cloud_top")
    add_ts(Stats, "env_cloud_cover")
    add_ts(Stats, "env_lwp")

    return
end

function io(self::EnvironmentVariables, Stats::NetCDFIO_Stats, Ref::ReferenceState)
    cinterior = self.Gr.cinterior
    finterior = self.Gr.finterior
    write_profile(Stats, "env_w", self.W.values[finterior])
    write_profile(Stats, "env_qt", self.QT.values[cinterior])
    write_profile(Stats, "env_ql", self.QL.values[cinterior])
    write_profile(Stats, "env_area", self.Area.values[cinterior])
    write_profile(Stats, "env_temperature", self.T.values[cinterior])
    write_profile(Stats, "env_RH", self.RH.values[cinterior])
    write_profile(Stats, "env_thetal", self.H.values[cinterior])
    write_profile(Stats, "env_tke", self.TKE.values[cinterior])
    write_profile(Stats, "env_Hvar", self.Hvar.values[cinterior])
    write_profile(Stats, "env_QTvar", self.QTvar.values[cinterior])
    write_profile(Stats, "env_HQTcov", self.HQTcov.values[cinterior])

    write_profile(Stats, "env_cloud_fraction", self.cloud_fraction.values[cinterior])

    env_cloud_diagnostics(self, Ref)
    # Assuming amximum overlap in environmental clouds
    write_ts(Stats, "env_cloud_cover", self.cloud_cover)
    write_ts(Stats, "env_cloud_base", self.cloud_base)
    write_ts(Stats, "env_cloud_top", self.cloud_top)
    write_ts(Stats, "env_lwp", self.lwp)
    return
end

function env_cloud_diagnostics(self::EnvironmentVariables, Ref::ReferenceState)
    self.cloud_top = 0.0
    self.cloud_base = zc_toa(self.Gr)
    self.cloud_cover = 0.0
    self.lwp = 0.0

    @inbounds for k in real_center_indicies(self.Gr)
        self.lwp += Ref.rho0_half[k] * self.QL.values[k] * self.Area.values[k] * self.Gr.dz

        if self.QL.values[k] > 1e-8 && self.Area.values[k] > 1e-3
            self.cloud_base = fmin(self.cloud_base, self.Gr.z_half[k])
            self.cloud_top = fmax(self.cloud_top, self.Gr.z_half[k])
            self.cloud_cover = fmax(self.cloud_cover, self.Area.values[k] * self.cloud_fraction.values[k])
        end
    end
    return
end

function update_EnvVar(self::EnvironmentThermodynamics, k, EnvVar::EnvironmentVariables, T, H, qt, ql, rho)
    EnvVar.T.values[k] = T
    EnvVar.H.values[k] = H
    EnvVar.QT.values[k] = qt
    EnvVar.QL.values[k] = ql
    EnvVar.B.values[k] = buoyancy_c(self.Ref.rho0_half[k], rho)
    EnvVar.RH.values[k] = relative_humidity_c(self.Ref.p0_half[k], qt, ql, 0.0, T)
    return
end

function update_EnvRain_sources(self::EnvironmentThermodynamics, k, EnvVar::EnvironmentVariables, qr_src, thl_rain_src)

    self.prec_source_qt[k] = -qr_src * EnvVar.Area.values[k]
    self.prec_source_h[k] = thl_rain_src * EnvVar.Area.values[k]

    return
end

function update_cloud_dry(self::EnvironmentThermodynamics, k, EnvVar::EnvironmentVariables, T, th, qt, ql, qv)
    if ql > 0.0
        EnvVar.cloud_fraction.values[k] = 1.0
        self.th_cloudy[k] = th
        self.t_cloudy[k] = T
        self.qt_cloudy[k] = qt
        self.qv_cloudy[k] = qv
    else
        EnvVar.cloud_fraction.values[k] = 0.0
        self.th_dry[k] = th
        self.qt_dry[k] = qt
    end
    return
end

function saturation_adjustment(self::EnvironmentThermodynamics, EnvVar::EnvironmentVariables)

    sa = eos_struct()
    mph = mph_struct()

    @inbounds for k in real_center_indicies(self.Gr)
        sa = eos(self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k])

        EnvVar.T.values[k] = sa.T
        EnvVar.QL.values[k] = sa.ql
        rho = rho_c(
            self.Ref.p0_half[k],
            EnvVar.T.values[k],
            EnvVar.QT.values[k],
            EnvVar.QT.values[k] - EnvVar.QL.values[k],
        )
        EnvVar.B.values[k] = buoyancy_c(self.Ref.rho0_half[k], rho)

        update_cloud_dry(
            self,
            k,
            EnvVar,
            EnvVar.T.values[k],
            EnvVar.H.values[k],
            EnvVar.QT.values[k],
            EnvVar.QL.values[k],
            EnvVar.QT.values[k] - EnvVar.QL.values[k],
        )
    end
    return
end


function sgs_mean(self::EnvironmentThermodynamics, EnvVar::EnvironmentVariables, Rain::RainVariables, dt)

    sa = eos_struct()
    mph = mph_struct()

    if EnvVar.H.name != "thetal"
        error("EDMF_Environment: rain source terms are defined for thetal as model variable")
    end

    @inbounds for k in real_center_indicies(self.Gr)
        # condensation
        sa = eos(self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k])
        # autoconversion and accretion
        mph = microphysics_rain_src(
            Rain.rain_model,
            Rain.max_supersaturation,
            Rain.C_drag,
            Rain.MP_n_0,
            Rain.q_liq_threshold,
            Rain.tau_acnv,
            Rain.E_col,
            EnvVar.QT.values[k],
            sa.ql,
            Rain.Env_QR.values[k],
            EnvVar.Area.values[k],
            sa.T,
            self.Ref.p0_half[k],
            self.Ref.rho0_half[k],
            dt,
        )
        update_EnvVar(self, k, EnvVar, sa.T, mph.thl, mph.qt, mph.ql, mph.rho)
        update_cloud_dry(self, k, EnvVar, sa.T, mph.th, mph.qt, mph.ql, mph.qv)
        update_EnvRain_sources(self, k, EnvVar, mph.qr_src, mph.thl_rain_src)
    end
    return
end

function sgs_quadrature(self::EnvironmentThermodynamics, EnvVar::EnvironmentVariables, Rain::RainVariables, dt)
    # TODO: double check this python-> julia translation
    # a, w = np.polynomial.hermite.hermgauss(self.quadrature_order)
    a, w = FastGaussQuadrature.gausshermite(self.quadrature_order)

    #TODO - remember you output source terms multipierd by dt (bec. of instanteneous autoconcv)
    #TODO - add tendencies for GMV H, QT and QR due to rain
    #TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    abscissas = a
    weights = w
    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)
    env_len = 10
    src_len = 6

    sqpi_inv = 1.0 / sqrt(π)
    sqrt2 = sqrt(2.0)
    sa = eos_struct()
    mph = mph_struct()

    epsilon = 10e-14 #np.finfo(np.float).eps
    if EnvVar.H.name != "thetal"
        error("EDMF_Environment: rain source terms are only defined for thetal as model variable")
    end

    # initialize the quadrature points and their labels
    inner_env = zeros(env_len)
    outer_env = zeros(env_len)
    inner_src = zeros(src_len)
    outer_src = zeros(src_len)
    i_ql, i_T, i_thl, i_rho, i_cf, i_qt_cld, i_qt_dry, i_T_cld, i_T_dry, i_rf = xrange(env_len)
    i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH = xrange(src_len)

    @inbounds for k in real_center_indicies(self.Gr)
        if (
            EnvVar.QTvar.values[k] > epsilon &&
            EnvVar.Hvar.values[k] > epsilon &&
            fabs(EnvVar.HQTcov.values[k]) > epsilon &&
            EnvVar.QT.values[k] > epsilon &&
            sqrt(EnvVar.QTvar.values[k]) < EnvVar.QT.values[k]
        )

            if self.quadrature_type == "log-normal"
                # Lognormal parameters (mu, sd) from mean and variance
                sd_q = sqrt(log(EnvVar.QTvar.values[k] / EnvVar.QT.values[k] / EnvVar.QT.values[k] + 1.0))
                sd_h = sqrt(log(EnvVar.Hvar.values[k] / EnvVar.H.values[k] / EnvVar.H.values[k] + 1.0))
                # Enforce Schwarz"s inequality
                corr = fmax(
                    fmin(EnvVar.HQTcov.values[k] / sqrt(EnvVar.Hvar.values[k] * EnvVar.QTvar.values[k]), 1.0),
                    -1.0,
                )
                sd2_hq = log(
                    corr * sqrt(EnvVar.Hvar.values[k] * EnvVar.QTvar.values[k]) / EnvVar.H.values[k] /
                    EnvVar.QT.values[k] + 1.0,
                )
                sd_cond_h_q = sqrt(fmax(sd_h * sd_h - sd2_hq * sd2_hq / sd_q / sd_q, 0.0))
                mu_q = log(
                    EnvVar.QT.values[k] * EnvVar.QT.values[k] /
                    sqrt(EnvVar.QT.values[k] * EnvVar.QT.values[k] + EnvVar.QTvar.values[k]),
                )
                mu_h = log(
                    EnvVar.H.values[k] * EnvVar.H.values[k] /
                    sqrt(EnvVar.H.values[k] * EnvVar.H.values[k] + EnvVar.Hvar.values[k]),
                )
            else
                sd_q = sqrt(EnvVar.QTvar.values[k])
                sd_h = sqrt(EnvVar.Hvar.values[k])
                corr = fmax(fmin(EnvVar.HQTcov.values[k] / fmax(sd_h * sd_q, 1e-13), 1.0), -1.0)

                # limit sd_q to prevent negative qt_hat
                sd_q_lim = (1e-10 - EnvVar.QT.values[k]) / (sqrt2 * abscissas[1])
                # walking backwards to assure your q_t will not be smaller than 1e-10
                # TODO - check
                # TODO - change 1e-13 and 1e-10 to some epislon
                sd_q = fmin(sd_q, sd_q_lim)
                qt_var = sd_q * sd_q
                sigma_h_star = sqrt(fmax(1.0 - corr * corr, 0.0)) * sd_h
            end

            # zero outer quadrature points
            for idx in xrange(env_len)
                outer_env[idx] = 0.0
            end
            for idx in xrange(src_len)
                outer_src[idx] = 0.0
            end

            for m_q in xrange(self.quadrature_order)
                if self.quadrature_type == "log-normal"
                    qt_hat = exp(mu_q + sqrt2 * sd_q * abscissas[m_q])
                    mu_h_star = mu_h + sd2_hq / sd_q / sd_q * (log(qt_hat) - mu_q)
                else
                    qt_hat = EnvVar.QT.values[k] + sqrt2 * sd_q * abscissas[m_q]
                    mu_h_star = EnvVar.H.values[k] + sqrt2 * corr * sd_h * abscissas[m_q]
                end

                # zero inner quadrature points
                for idx in xrange(env_len)
                    inner_env[idx] = 0.0
                end
                for idx in xrange(src_len)
                    inner_src[idx] = 0.0
                end

                for m_h in xrange(self.quadrature_order)
                    if self.quadrature_type == "log-normal"
                        h_hat = exp(mu_h_star + sqrt2 * sd_cond_h_q * abscissas[m_h])
                    else
                        h_hat = sqrt2 * sigma_h_star * abscissas[m_h] + mu_h_star
                    end

                    # condensation
                    sa = eos(self.Ref.p0_half[k], qt_hat, h_hat)
                    # autoconversion and accretion
                    mph = microphysics_rain_src(
                        Rain.rain_model,
                        Rain.max_supersaturation,
                        Rain.C_drag,
                        Rain.MP_n_0,
                        Rain.q_liq_threshold,
                        Rain.tau_acnv,
                        Rain.E_col,
                        qt_hat,
                        sa.ql,
                        Rain.Env_QR.values[k],
                        EnvVar.Area.values[k],
                        sa.T,
                        self.Ref.p0_half[k],
                        self.Ref.rho0_half[k],
                        dt,
                    )

                    # environmental variables
                    inner_env[i_ql] += mph.ql * weights[m_h] * sqpi_inv
                    inner_env[i_T] += sa.T * weights[m_h] * sqpi_inv
                    inner_env[i_thl] += mph.thl * weights[m_h] * sqpi_inv
                    inner_env[i_rho] += mph.rho * weights[m_h] * sqpi_inv
                    # rain area fraction
                    if mph.qr_src > 0.0
                        inner_env[i_rf] += weights[m_h] * sqpi_inv
                    end
                    # cloudy/dry categories for buoyancy in TKE
                    if mph.ql > 0.0
                        inner_env[i_cf] += weights[m_h] * sqpi_inv
                        inner_env[i_qt_cld] += mph.qt * weights[m_h] * sqpi_inv
                        inner_env[i_T_cld] += sa.T * weights[m_h] * sqpi_inv
                    else
                        inner_env[i_qt_dry] += mph.qt * weights[m_h] * sqpi_inv
                        inner_env[i_T_dry] += sa.T * weights[m_h] * sqpi_inv
                    end
                    # products for variance and covariance source terms
                    inner_src[i_Sqt] += -mph.qr_src * weights[m_h] * sqpi_inv
                    inner_src[i_SH] += mph.thl_rain_src * weights[m_h] * sqpi_inv
                    inner_src[i_Sqt_H] += -mph.qr_src * mph.thl * weights[m_h] * sqpi_inv
                    inner_src[i_Sqt_qt] += -mph.qr_src * mph.qt * weights[m_h] * sqpi_inv
                    inner_src[i_SH_H] += mph.thl_rain_src * mph.thl * weights[m_h] * sqpi_inv
                    inner_src[i_SH_qt] += mph.thl_rain_src * mph.qt * weights[m_h] * sqpi_inv
                end

                for idx in xrange(env_len)
                    outer_env[idx] += inner_env[idx] * weights[m_q] * sqpi_inv
                end
                for idx in xrange(src_len)
                    outer_src[idx] += inner_src[idx] * weights[m_q] * sqpi_inv
                end
            end

            # update environmental variables
            update_EnvVar(
                self,
                k,
                EnvVar,
                outer_env[i_T],
                outer_env[i_thl],
                outer_env[i_qt_cld] + outer_env[i_qt_dry],
                outer_env[i_ql],
                outer_env[i_rho],
            )
            update_EnvRain_sources(self, k, EnvVar, -outer_src[i_Sqt], outer_src[i_SH])

            # update cloudy/dry variables for buoyancy in TKE
            EnvVar.cloud_fraction.values[k] = outer_env[i_cf]
            self.qt_dry[k] = outer_env[i_qt_dry]
            self.th_dry[k] = theta_c(self.Ref.p0_half[k], outer_env[i_T_dry])
            self.t_cloudy[k] = outer_env[i_T_cld]
            self.qv_cloudy[k] = outer_env[i_qt_cld] - outer_env[i_ql]
            self.qt_cloudy[k] = outer_env[i_qt_cld]
            self.th_cloudy[k] = theta_c(self.Ref.p0_half[k], outer_env[i_T_cld])

            # update var/covar rain sources
            self.Hvar_rain_dt[k] = outer_src[i_SH_H] - outer_src[i_SH] * EnvVar.H.values[k]
            self.QTvar_rain_dt[k] = outer_src[i_Sqt_qt] - outer_src[i_Sqt] * EnvVar.QT.values[k]
            self.HQTcov_rain_dt[k] =
                outer_src[i_SH_qt] - outer_src[i_SH] * EnvVar.QT.values[k] + outer_src[i_Sqt_H] -
                outer_src[i_Sqt] * EnvVar.H.values[k]

        else
            # if variance and covariance are zero do the same as in SA_mean
            sa = eos(self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k])
            mph = microphysics_rain_src(
                Rain.rain_model,
                Rain.max_supersaturation,
                Rain.C_drag,
                Rain.MP_n_0,
                Rain.q_liq_threshold,
                Rain.tau_acnv,
                Rain.E_col,
                EnvVar.QT.values[k],
                sa.ql,
                Rain.Env_QR.values[k],
                EnvVar.Area.values[k],
                sa.T,
                self.Ref.p0_half[k],
                self.Ref.rho0_half[k],
                dt,
            )
            update_EnvVar(self, k, EnvVar, sa.T, mph.thl, mph.qt, mph.ql, mph.rho)
            update_EnvRain_sources(self, k, EnvVar, mph.qr_src, mph.thl_rain_src)
            update_cloud_dry(self, k, EnvVar, sa.T, mph.th, mph.qt, mph.ql, mph.qv)

            self.Hvar_rain_dt[k] = 0.0
            self.QTvar_rain_dt[k] = 0.0
            self.HQTcov_rain_dt[k] = 0.0
        end
    end

    return
end

function microphysics(self::EnvironmentThermodynamics, EnvVar::EnvironmentVariables, Rain::RainVariables, dt)

    if EnvVar.EnvThermo_scheme == "mean"
        sgs_mean(self, EnvVar, Rain, dt)

    elseif EnvVar.EnvThermo_scheme == "quadrature"
        sgs_quadrature(self, EnvVar, Rain, dt)

    else
        error("EDMF_Environment: Unrecognized EnvThermo_scheme. Possible options: mean, quadrature")
    end

    return
end
