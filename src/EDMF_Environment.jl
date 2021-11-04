function update_env_precip_tendencies(state, k, qt_tendency, θ_liq_ice_tendency, qr_tendency, qs_tendency)

    aux_en = center_aux_environment(state)
    tendencies_pr = center_tendencies_precipitation(state)

    aux_en.qt_tendency_precip_formation[k] = qt_tendency * aux_en.area[k]
    aux_en.θ_liq_ice_tendency_precip_formation[k] = θ_liq_ice_tendency * aux_en.area[k]

    tendencies_pr.q_rai[k] += qr_tendency * aux_en.area[k]
    tendencies_pr.q_sno[k] += qs_tendency * aux_en.area[k]
    return
end

function update_sat_unsat(en_thermo::EnvironmentThermodynamics, state, k, ts)
    q_liq = TD.liquid_specific_humidity(ts)
    q_ice = TD.ice_specific_humidity(ts)
    aux_en = center_aux_environment(state)
    if TD.has_condensate(q_liq + q_ice)
        aux_en.cloud_fraction[k] = 1.0
        en_thermo.θ_sat[k] = TD.dry_pottemp(ts)
        en_thermo.θ_liq_ice_sat[k] = TD.liquid_ice_pottemp(ts)
        en_thermo.t_sat[k] = TD.air_temperature(ts)
        en_thermo.qt_sat[k] = TD.total_specific_humidity(ts)
        en_thermo.qv_sat[k] = TD.vapor_specific_humidity(ts)
    else
        aux_en.cloud_fraction[k] = 0.0
        en_thermo.θ_unsat[k] = TD.dry_pottemp(ts)
        en_thermo.θv_unsat[k] = TD.virtual_pottemp(ts)
        en_thermo.qt_unsat[k] = TD.total_specific_humidity(ts)
    end
    return
end

function sgs_mean(en_thermo::EnvironmentThermodynamics, grid, state, en, precip, dt, param_set)

    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)

    @inbounds for k in real_center_indices(grid)
        # condensation
        q_tot_en = aux_en.q_tot[k]
        ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], q_tot_en)
        # autoconversion and accretion
        mph = precipitation_formation(
            param_set,
            precip.precipitation_model,
            prog_pr.q_rai[k],
            prog_pr.q_sno[k],
            aux_en.area[k],
            ρ0_c[k],
            dt,
            ts,
        )
        update_sat_unsat(en_thermo, state, k, ts)
        update_env_precip_tendencies(
            state,
            k,
            mph.qt_tendency,
            mph.θ_liq_ice_tendency,
            mph.qr_tendency,
            mph.qs_tendency,
        )
    end
    return
end

function sgs_quadrature(en_thermo::EnvironmentThermodynamics, grid, state, en, precip, dt, param_set)
    # TODO: double check this python-> julia translation
    # a, w = np.polynomial.hermite.hermgauss(en_thermo.quadrature_order)
    a, w = FastGaussQuadrature.gausshermite(en_thermo.quadrature_order)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)

    #TODO - remember you output source terms multipierd by dt (bec. of instanteneous autoconcv)
    #TODO - add tendencies for gm H, QT and QR due to rain
    #TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    abscissas = a
    weights = w
    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)
    env_len = 8
    src_len = 8

    sqpi_inv = 1.0 / sqrt(π)
    sqrt2 = sqrt(2.0)

    epsilon = 10e-14 #np.finfo(np.float).eps

    # initialize the quadrature points and their labels
    inner_env = zeros(env_len)
    outer_env = zeros(env_len)
    inner_src = zeros(src_len)
    outer_src = zeros(src_len)
    i_ql, i_qi, i_T, i_cf, i_qt_sat, i_qt_unsat, i_T_sat, i_T_unsat = 1:env_len
    i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH, i_Sqr, i_Sqs = 1:src_len

    @inbounds for k in real_center_indices(grid)
        if (
            aux_en.QTvar[k] > epsilon &&
            aux_en.Hvar[k] > epsilon &&
            abs(aux_en.HQTcov[k]) > epsilon &&
            aux_en.q_tot[k] > epsilon &&
            sqrt(aux_en.QTvar[k]) < aux_en.q_tot[k]
        )

            if en_thermo.quadrature_type == "log-normal"
                # Lognormal parameters (mu, sd) from mean and variance
                sd_q = sqrt(log(aux_en.QTvar[k] / aux_en.q_tot[k] / aux_en.q_tot[k] + 1.0))
                sd_h = sqrt(log(aux_en.Hvar[k] / aux_en.θ_liq_ice[k] / aux_en.θ_liq_ice[k] + 1.0))
                # Enforce Schwarz"s inequality
                corr = max(min(aux_en.HQTcov[k] / sqrt(aux_en.Hvar[k] * aux_en.QTvar[k]), 1.0), -1.0)
                sd2_hq =
                    log(corr * sqrt(aux_en.Hvar[k] * aux_en.QTvar[k]) / aux_en.θ_liq_ice[k] / aux_en.q_tot[k] + 1.0)
                sd_cond_h_q = sqrt(max(sd_h * sd_h - sd2_hq * sd2_hq / sd_q / sd_q, 0.0))
                mu_q =
                    log(aux_en.q_tot[k] * aux_en.q_tot[k] / sqrt(aux_en.q_tot[k] * aux_en.q_tot[k] + aux_en.QTvar[k]))
                mu_h = log(
                    aux_en.θ_liq_ice[k] * aux_en.θ_liq_ice[k] /
                    sqrt(aux_en.θ_liq_ice[k] * aux_en.θ_liq_ice[k] + aux_en.Hvar[k]),
                )
            else
                sd_q = sqrt(aux_en.QTvar[k])
                sd_h = sqrt(aux_en.Hvar[k])
                corr = max(min(aux_en.HQTcov[k] / max(sd_h * sd_q, 1e-13), 1.0), -1.0)

                # limit sd_q to prevent negative qt_hat
                sd_q_lim = (1e-10 - aux_en.q_tot[k]) / (sqrt2 * abscissas[1])
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
                    qt_hat = aux_en.q_tot[k] + sqrt2 * sd_q * abscissas[m_q]
                    mu_h_star = aux_en.θ_liq_ice[k] + sqrt2 * corr * sd_h * abscissas[m_q]
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
                    ts = thermo_state_pθq(param_set, p0_c[k], h_hat, qt_hat)
                    q_liq_en = TD.liquid_specific_humidity(ts)
                    q_ice_en = TD.ice_specific_humidity(ts)
                    T = TD.air_temperature(ts)
                    # autoconversion and accretion
                    mph = precipitation_formation(
                        param_set,
                        precip.precipitation_model,
                        prog_pr.q_rai[k],
                        prog_pr.q_sno[k],
                        aux_en.area[k],
                        ρ0_c[k],
                        dt,
                        ts,
                    )

                    # environmental variables
                    inner_env[i_ql] += q_liq_en * weights[m_h] * sqpi_inv
                    inner_env[i_qi] += q_ice_en * weights[m_h] * sqpi_inv
                    inner_env[i_T] += T * weights[m_h] * sqpi_inv
                    # cloudy/dry categories for buoyancy in TKE
                    if TD.has_condensate(q_liq_en + q_ice_en)
                        inner_env[i_cf] += weights[m_h] * sqpi_inv
                        inner_env[i_qt_sat] += qt_hat * weights[m_h] * sqpi_inv
                        inner_env[i_T_sat] += T * weights[m_h] * sqpi_inv
                    else
                        inner_env[i_qt_unsat] += qt_hat * weights[m_h] * sqpi_inv
                        inner_env[i_T_unsat] += T * weights[m_h] * sqpi_inv
                    end
                    # products for variance and covariance source terms
                    inner_src[i_Sqt] += mph.qt_tendency * weights[m_h] * sqpi_inv
                    inner_src[i_Sqr] += mph.qr_tendency * weights[m_h] * sqpi_inv
                    inner_src[i_Sqs] += mph.qs_tendency * weights[m_h] * sqpi_inv
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
            update_env_precip_tendencies(
                state,
                k,
                outer_src[i_Sqt],
                outer_src[i_SH],
                outer_src[i_Sqr],
                outer_src[i_Sqs],
            )

            # update cloudy/dry variables for buoyancy in TKE
            aux_en.cloud_fraction[k] = outer_env[i_cf]
            en_thermo.qt_unsat[k] = outer_env[i_qt_unsat]
            # TD.jl cannot compute θ_unsat when T=0
            en_thermo.θ_unsat[k] = if outer_env[i_T_unsat] > 0
                ts_unsat = TD.PhaseEquil_pTq(param_set, p0_c[k], outer_env[i_T_unsat], en_thermo.qt_unsat[k])
                TD.dry_pottemp(ts_unsat)
            else
                0
            end

            en_thermo.t_sat[k] = outer_env[i_T_sat]
            en_thermo.qv_sat[k] = outer_env[i_qt_sat] - outer_env[i_ql]
            en_thermo.qt_sat[k] = outer_env[i_qt_sat]
            ts_sat = TD.PhaseEquil_pTq(param_set, p0_c[k], en_thermo.t_sat[k], en_thermo.qt_sat[k])
            en_thermo.θ_sat[k] = TD.dry_pottemp(ts_sat)
            en_thermo.θ_liq_ice_sat[k] = TD.liquid_ice_pottemp(ts_sat)

            # update var/covar rain sources
            en_thermo.Hvar_rain_dt[k] = outer_src[i_SH_H] - outer_src[i_SH] * aux_en.θ_liq_ice[k]
            en_thermo.QTvar_rain_dt[k] = outer_src[i_Sqt_qt] - outer_src[i_Sqt] * aux_en.q_tot[k]
            en_thermo.HQTcov_rain_dt[k] =
                outer_src[i_SH_qt] - outer_src[i_SH] * aux_en.q_tot[k] + outer_src[i_Sqt_H] -
                outer_src[i_Sqt] * aux_en.θ_liq_ice[k]

        else
            # if variance and covariance are zero do the same as in SA_mean
            ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])
            mph = precipitation_formation(
                param_set,
                precip.precipitation_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_en.area[k],
                ρ0_c[k],
                dt,
                ts,
            )
            update_env_precip_tendencies(
                state,
                k,
                mph.qt_tendency,
                mph.θ_liq_ice_tendency,
                mph.qr_tendency,
                mph.qs_tendency,
            )
            update_sat_unsat(en_thermo, state, k, ts)

            en_thermo.Hvar_rain_dt[k] = 0.0
            en_thermo.QTvar_rain_dt[k] = 0.0
            en_thermo.HQTcov_rain_dt[k] = 0.0
        end
    end

    return
end

function microphysics(en_thermo::EnvironmentThermodynamics, grid, state, en, precip, dt, param_set)

    if en.EnvThermo_scheme == "mean"
        sgs_mean(en_thermo, grid, state, en, precip, dt, param_set)

    elseif en.EnvThermo_scheme == "quadrature"
        sgs_quadrature(en_thermo, grid, state, en, precip, dt, param_set)

    else
        error("EDMF_Environment: Unrecognized EnvThermo_scheme. Possible options: mean, quadrature")
    end

    return
end
