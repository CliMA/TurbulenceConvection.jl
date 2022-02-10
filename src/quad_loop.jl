const i_ql, i_qi, i_T, i_cf, i_qt_sat, i_qt_unsat, i_T_sat, i_T_unsat = 1:8
const i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH, i_Sqr, i_Sqs = 1:8

function quad_loop(en_thermo::SGSQuadrature, precip_model, vars, param_set, Δt::Real)
    a = en_thermo.a
    w = en_thermo.w
    quadrature_type = en_thermo.quadrature_type
    quad_order = quadrature_order(en_thermo)

    UnPack.@unpack QTvar_en, q_tot_en, Hvar_en, θ_liq_ice_en, HQTcov_en, area_en, q_rai, q_sno, area_en, ρ0_c, p0_c =
        vars

    env_len = 8
    src_len = 8

    inner_env = zeros(env_len)
    outer_env = zeros(env_len)
    inner_src = zeros(src_len)
    outer_src = zeros(src_len)

    abscissas = a
    weights = w

    sqpi_inv = 1 / sqrt(π)
    sqrt2 = sqrt(2)

    if quadrature_type isa LogNormalQuad
        # Lognormal parameters (mu, sd) from mean and variance
        sd_q = sqrt(log(QTvar_en / q_tot_en / q_tot_en + 1))
        sd_h = sqrt(log(Hvar_en / θ_liq_ice_en / θ_liq_ice_en + 1))
        # Enforce Schwarz"s inequality
        corr = max(min(HQTcov_en / sqrt(Hvar_en * QTvar_en), 1), -1)
        sd2_hq = log(corr * sqrt(Hvar_en * QTvar_en) / θ_liq_ice_en / q_tot_en + 1)
        sd_cond_h_q = sqrt(max(sd_h * sd_h - sd2_hq * sd2_hq / sd_q / sd_q, 0))
        mu_q = log(q_tot_en * q_tot_en / sqrt(q_tot_en * q_tot_en + QTvar_en))
        mu_h = log(θ_liq_ice_en * θ_liq_ice_en / sqrt(θ_liq_ice_en * θ_liq_ice_en + Hvar_en))
    else
        sd_q = sqrt(QTvar_en)
        sd_h = sqrt(Hvar_en)
        corr = max(min(HQTcov_en / max(sd_h * sd_q, 1e-13), 1), -1)

        # limit sd_q to prevent negative qt_hat
        sd_q_lim = (1e-10 - q_tot_en) / (sqrt2 * abscissas[1])
        # walking backwards to assure your q_t will not be smaller than 1e-10
        # TODO - check
        # TODO - change 1e-13 and 1e-10 to some epislon
        sd_q = min(sd_q, sd_q_lim)
        qt_var = sd_q * sd_q
        σ_h_star = sqrt(max(1 - corr * corr, 0)) * sd_h
    end

    # zero outer quadrature points
    @inbounds for idx in 1:env_len
        outer_env[idx] = 0.0
    end
    @inbounds for idx in 1:src_len
        outer_src[idx] = 0.0
    end

    @inbounds for m_q in 1:quad_order
        if quadrature_type isa LogNormalQuad
            qt_hat = exp(mu_q + sqrt2 * sd_q * abscissas[m_q])
            μ_h_star = mu_h + sd2_hq / sd_q / sd_q * (log(qt_hat) - mu_q)
        else
            qt_hat = q_tot_en + sqrt2 * sd_q * abscissas[m_q]
            μ_h_star = θ_liq_ice_en + sqrt2 * corr * sd_h * abscissas[m_q]
        end

        # zero inner quadrature points
        for idx in 1:env_len
            inner_env[idx] = 0.0
        end
        for idx in 1:src_len
            inner_src[idx] = 0.0
        end

        for m_h in 1:quad_order
            if quadrature_type isa LogNormalQuad
                h_hat = exp(μ_h_star + sqrt2 * sd_cond_h_q * abscissas[m_h])
            else
                h_hat = sqrt2 * σ_h_star * abscissas[m_h] + μ_h_star
            end

            # condensation
            ts = thermo_state_pθq(param_set, p0_c, h_hat, qt_hat)
            q_liq_en = TD.liquid_specific_humidity(ts)
            q_ice_en = TD.ice_specific_humidity(ts)
            T = TD.air_temperature(ts)
            # autoconversion and accretion
            mph = precipitation_formation(param_set, precip_model, q_rai, q_sno, area_en, ρ0_c, Δt, ts)

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
    return outer_env, outer_src
end
