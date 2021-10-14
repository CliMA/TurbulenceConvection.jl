function initialize_io(en::EnvironmentVariables, Stats::NetCDFIO_Stats)
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

function io(en::EnvironmentVariables, grid, state, Stats::NetCDFIO_Stats)
    write_profile(Stats, "env_w", en.W.values)
    write_profile(Stats, "env_qt", en.QT.values)
    write_profile(Stats, "env_ql", en.QL.values)
    write_profile(Stats, "env_area", en.Area.values)
    write_profile(Stats, "env_temperature", en.T.values)
    write_profile(Stats, "env_RH", en.RH.values)
    write_profile(Stats, "env_thetal", en.H.values)
    write_profile(Stats, "env_tke", en.TKE.values)
    write_profile(Stats, "env_Hvar", en.Hvar.values)
    write_profile(Stats, "env_QTvar", en.QTvar.values)
    write_profile(Stats, "env_HQTcov", en.HQTcov.values)

    write_profile(Stats, "env_cloud_fraction", en.cloud_fraction.values)

    env_cloud_diagnostics(en, grid, state)
    # Assuming amximum overlap in environmental clouds
    write_ts(Stats, "env_cloud_cover", en.cloud_cover)
    write_ts(Stats, "env_cloud_base", en.cloud_base)
    write_ts(Stats, "env_cloud_top", en.cloud_top)
    write_ts(Stats, "env_lwp", en.lwp)
    return
end

function env_cloud_diagnostics(en::EnvironmentVariables, grid, state)
    en.cloud_top = 0.0
    en.cloud_base = zc_toa(grid)
    en.cloud_cover = 0.0
    en.lwp = 0.0
    ρ0_c = center_ref_state(state).ρ0

    @inbounds for k in real_center_indices(grid)
        en.lwp += ρ0_c[k] * en.QL.values[k] * en.Area.values[k] * grid.Δz

        if en.QL.values[k] > 1e-8 && en.Area.values[k] > 1e-3
            en.cloud_base = min(en.cloud_base, grid.zc[k])
            en.cloud_top = max(en.cloud_top, grid.zc[k])
            en.cloud_cover = max(en.cloud_cover, en.Area.values[k] * en.cloud_fraction.values[k])
        end
    end
    return
end

function update_env_precip_tendencies(en_thermo::EnvironmentThermodynamics, k, en, qt_tendency, θ_liq_ice_tendency)

    en_thermo.qt_tendency_rain_formation[k] = qt_tendency * en.Area.values[k]
    en_thermo.θ_liq_ice_tendency_rain_formation[k] = θ_liq_ice_tendency * en.Area.values[k]

    return
end

function update_cloud_dry(en_thermo::EnvironmentThermodynamics, k, en::EnvironmentVariables, ts)
    q_liq = TD.liquid_specific_humidity(ts)
    if q_liq > 0.0
        en.cloud_fraction.values[k] = 1.0
        en_thermo.th_cloudy[k] = TD.liquid_ice_pottemp(ts)
        en_thermo.t_cloudy[k] = TD.air_temperature(ts)
        en_thermo.qt_cloudy[k] = TD.total_specific_humidity(ts)
        en_thermo.qv_cloudy[k] = TD.vapor_specific_humidity(ts)
    else
        en.cloud_fraction.values[k] = 0.0
        en_thermo.th_dry[k] = TD.liquid_ice_pottemp(ts)
        en_thermo.qt_dry[k] = TD.total_specific_humidity(ts)
    end
    return
end

function sgs_mean(en_thermo::EnvironmentThermodynamics, grid, state, en, rain, dt, param_set)

    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0

    @inbounds for k in real_center_indices(grid)
        # condensation
        q_tot_en = en.QT.values[k]
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], en.H.values[k], q_tot_en)
        # autoconversion and accretion
        mph = precipitation_formation(param_set, rain.rain_model, rain.QR.values[k], en.Area.values[k], ρ0_c[k], dt, ts)
        update_cloud_dry(en_thermo, k, en, ts)
        update_env_precip_tendencies(en_thermo, k, en, mph.qt_tendency, mph.θ_liq_ice_tendency)
    end
    return
end

function sgs_quadrature(en_thermo::EnvironmentThermodynamics, grid, state, en, rain, dt, param_set)
    # TODO: double check this python-> julia translation
    # a, w = np.polynomial.hermite.hermgauss(en_thermo.quadrature_order)
    a, w = FastGaussQuadrature.gausshermite(en_thermo.quadrature_order)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0

    #TODO - remember you output source terms multipierd by dt (bec. of instanteneous autoconcv)
    #TODO - add tendencies for GMV H, QT and QR due to rain
    #TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    abscissas = a
    weights = w
    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)
    env_len = 8
    src_len = 6

    sqpi_inv = 1.0 / sqrt(π)
    sqrt2 = sqrt(2.0)

    epsilon = 10e-14 #np.finfo(np.float).eps

    # initialize the quadrature points and their labels
    inner_env = zeros(env_len)
    outer_env = zeros(env_len)
    inner_src = zeros(src_len)
    outer_src = zeros(src_len)
    i_ql, i_T, i_cf, i_qt_cld, i_qt_dry, i_T_cld, i_T_dry, i_rf = 1:env_len
    i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH = 1:src_len

    @inbounds for k in real_center_indices(grid)
        if (
            en.QTvar.values[k] > epsilon &&
            en.Hvar.values[k] > epsilon &&
            abs(en.HQTcov.values[k]) > epsilon &&
            en.QT.values[k] > epsilon &&
            sqrt(en.QTvar.values[k]) < en.QT.values[k]
        )

            if en_thermo.quadrature_type == "log-normal"
                # Lognormal parameters (mu, sd) from mean and variance
                sd_q = sqrt(log(en.QTvar.values[k] / en.QT.values[k] / en.QT.values[k] + 1.0))
                sd_h = sqrt(log(en.Hvar.values[k] / en.H.values[k] / en.H.values[k] + 1.0))
                # Enforce Schwarz"s inequality
                corr = max(min(en.HQTcov.values[k] / sqrt(en.Hvar.values[k] * en.QTvar.values[k]), 1.0), -1.0)
                sd2_hq =
                    log(corr * sqrt(en.Hvar.values[k] * en.QTvar.values[k]) / en.H.values[k] / en.QT.values[k] + 1.0)
                sd_cond_h_q = sqrt(max(sd_h * sd_h - sd2_hq * sd2_hq / sd_q / sd_q, 0.0))
                mu_q = log(
                    en.QT.values[k] * en.QT.values[k] / sqrt(en.QT.values[k] * en.QT.values[k] + en.QTvar.values[k]),
                )
                mu_h = log(en.H.values[k] * en.H.values[k] / sqrt(en.H.values[k] * en.H.values[k] + en.Hvar.values[k]))
            else
                sd_q = sqrt(en.QTvar.values[k])
                sd_h = sqrt(en.Hvar.values[k])
                corr = max(min(en.HQTcov.values[k] / max(sd_h * sd_q, 1e-13), 1.0), -1.0)

                # limit sd_q to prevent negative qt_hat
                sd_q_lim = (1e-10 - en.QT.values[k]) / (sqrt2 * abscissas[1])
                # walking backwards to assure your q_t will not be smaller than 1e-10
                # TODO - check
                # TODO - change 1e-13 and 1e-10 to some epislon
                sd_q = min(sd_q, sd_q_lim)
                qt_var = sd_q * sd_q
                sigma_h_star = sqrt(max(1.0 - corr * corr, 0.0)) * sd_h
            end

            # zero outer quadrature points
            for idx in 1:env_len
                outer_env[idx] = 0.0
            end
            for idx in 1:src_len
                outer_src[idx] = 0.0
            end

            for m_q in 1:(en_thermo.quadrature_order)
                if en_thermo.quadrature_type == "log-normal"
                    qt_hat = exp(mu_q + sqrt2 * sd_q * abscissas[m_q])
                    mu_h_star = mu_h + sd2_hq / sd_q / sd_q * (log(qt_hat) - mu_q)
                else
                    qt_hat = en.QT.values[k] + sqrt2 * sd_q * abscissas[m_q]
                    mu_h_star = en.H.values[k] + sqrt2 * corr * sd_h * abscissas[m_q]
                end

                # zero inner quadrature points
                for idx in 1:env_len
                    inner_env[idx] = 0.0
                end
                for idx in 1:src_len
                    inner_src[idx] = 0.0
                end

                for m_h in 1:(en_thermo.quadrature_order)
                    if en_thermo.quadrature_type == "log-normal"
                        h_hat = exp(mu_h_star + sqrt2 * sd_cond_h_q * abscissas[m_h])
                    else
                        h_hat = sqrt2 * sigma_h_star * abscissas[m_h] + mu_h_star
                    end

                    # condensation
                    ts = TD.PhaseEquil_pθq(param_set, p0_c[k], h_hat, qt_hat)
                    q_liq_en = TD.liquid_specific_humidity(ts)
                    T = TD.air_temperature(ts)
                    # autoconversion and accretion
                    mph = precipitation_formation(
                        param_set,
                        rain.rain_model,
                        rain.QR.values[k],
                        en.Area.values[k],
                        ρ0_c[k],
                        dt,
                        ts,
                    )

                    # environmental variables
                    inner_env[i_ql] += q_liq_en * weights[m_h] * sqpi_inv
                    inner_env[i_T] += T * weights[m_h] * sqpi_inv
                    # rain area fraction
                    if mph.qr_tendency > 0.0
                        inner_env[i_rf] += weights[m_h] * sqpi_inv
                    end
                    # cloudy/dry categories for buoyancy in TKE
                    if q_liq_en > 0.0
                        inner_env[i_cf] += weights[m_h] * sqpi_inv
                        inner_env[i_qt_cld] += qt_hat * weights[m_h] * sqpi_inv
                        inner_env[i_T_cld] += T * weights[m_h] * sqpi_inv
                    else
                        inner_env[i_qt_dry] += qt_hat * weights[m_h] * sqpi_inv
                        inner_env[i_T_dry] += T * weights[m_h] * sqpi_inv
                    end
                    # products for variance and covariance source terms
                    inner_src[i_Sqt] += mph.qt_tendency * weights[m_h] * sqpi_inv
                    inner_src[i_SH] += mph.θ_liq_ice_tendency * weights[m_h] * sqpi_inv
                    inner_src[i_Sqt_H] += mph.qt_tendency * h_hat * weights[m_h] * sqpi_inv
                    inner_src[i_Sqt_qt] += mph.qt_tendency * qt_hat * weights[m_h] * sqpi_inv
                    inner_src[i_SH_H] += mph.θ_liq_ice_tendency * h_hat * weights[m_h] * sqpi_inv
                    inner_src[i_SH_qt] += mph.θ_liq_ice_tendency * qt_hat * weights[m_h] * sqpi_inv
                end

                for idx in 1:env_len
                    outer_env[idx] += inner_env[idx] * weights[m_q] * sqpi_inv
                end
                for idx in 1:src_len
                    outer_src[idx] += inner_src[idx] * weights[m_q] * sqpi_inv
                end
            end

            # update environmental variables
            update_env_precip_tendencies(en_thermo, k, en, outer_src[i_Sqt], outer_src[i_SH])

            # update cloudy/dry variables for buoyancy in TKE
            en.cloud_fraction.values[k] = outer_env[i_cf]
            en_thermo.qt_dry[k] = outer_env[i_qt_dry]
            # TD.jl cannot compute θ_dry when T=0
            en_thermo.th_dry[k] = if outer_env[i_T_dry] > 0
                ts_dry = TD.PhaseEquil_pTq(param_set, p0_c[k], outer_env[i_T_dry], en_thermo.qt_dry[k])
                TD.dry_pottemp(ts_dry)
            else
                0
            end

            en_thermo.t_cloudy[k] = outer_env[i_T_cld]
            en_thermo.qv_cloudy[k] = outer_env[i_qt_cld] - outer_env[i_ql]
            en_thermo.qt_cloudy[k] = outer_env[i_qt_cld]
            ts_cld = TD.PhaseEquil_pTq(param_set, p0_c[k], en_thermo.t_cloudy[k], en_thermo.qt_cloudy[k])
            en_thermo.th_cloudy[k] = TD.dry_pottemp(ts_cld)

            # update var/covar rain sources
            en_thermo.Hvar_rain_dt[k] = outer_src[i_SH_H] - outer_src[i_SH] * en.H.values[k]
            en_thermo.QTvar_rain_dt[k] = outer_src[i_Sqt_qt] - outer_src[i_Sqt] * en.QT.values[k]
            en_thermo.HQTcov_rain_dt[k] =
                outer_src[i_SH_qt] - outer_src[i_SH] * en.QT.values[k] + outer_src[i_Sqt_H] -
                outer_src[i_Sqt] * en.H.values[k]

        else
            # if variance and covariance are zero do the same as in SA_mean
            ts = TD.PhaseEquil_pθq(param_set, p0_c[k], en.H.values[k], en.QT.values[k])
            mph = precipitation_formation(
                param_set,
                rain.rain_model,
                rain.QR.values[k],
                en.Area.values[k],
                ρ0_c[k],
                dt,
                ts,
            )
            update_env_precip_tendencies(en_thermo, k, en, mph.qt_tendency, mph.θ_liq_ice_tendency)
            update_cloud_dry(en_thermo, k, en, ts)

            en_thermo.Hvar_rain_dt[k] = 0.0
            en_thermo.QTvar_rain_dt[k] = 0.0
            en_thermo.HQTcov_rain_dt[k] = 0.0
        end
    end

    return
end

function microphysics(en_thermo::EnvironmentThermodynamics, grid, state, en, rain, dt, param_set)

    if en.EnvThermo_scheme == "mean"
        sgs_mean(en_thermo, grid, state, en, rain, dt, param_set)

    elseif en.EnvThermo_scheme == "quadrature"
        sgs_quadrature(en_thermo, grid, state, en, rain, dt, param_set)

    else
        error("EDMF_Environment: Unrecognized EnvThermo_scheme. Possible options: mean, quadrature")
    end

    return
end
