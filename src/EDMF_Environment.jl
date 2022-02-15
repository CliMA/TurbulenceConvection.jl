function microphysics(
    ::SGSMean,
    grid::Grid,
    state::State,
    precip_model::AbstractPrecipitationModel,
    Δt::Real,
    param_set::APS,
)

    tendencies_pr = center_tendencies_precipitation(state)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    aux_en_sat = aux_en.sat
    aux_en_unsat = aux_en.unsat

    @inbounds for k in real_center_indices(grid)
        # condensation
        q_tot_en = aux_en.q_tot[k]
        ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], q_tot_en)
        # autoconversion and accretion
        mph = precipitation_formation(
            param_set,
            precip_model,
            prog_pr.q_rai[k],
            prog_pr.q_sno[k],
            aux_en.area[k],
            ρ0_c[k],
            Δt,
            ts,
        )

        # update_sat_unsat
        if TD.has_condensate(ts)
            aux_en.cloud_fraction[k] = 1
            aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts)
            aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts)
            aux_en_sat.T[k] = TD.air_temperature(ts)
            aux_en_sat.q_tot[k] = TD.total_specific_humidity(ts)
            aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(ts)
        else
            aux_en.cloud_fraction[k] = 0
            aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts)
            aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(ts)
            aux_en_unsat.q_tot[k] = TD.total_specific_humidity(ts)
        end

        # update_env_precip_tendencies
        # TODO: move qt_tendency_precip_formation and θ_liq_ice_tendency_precip_formation
        # to diagnostics
        aux_en.qt_tendency_precip_formation[k] = mph.qt_tendency * aux_en.area[k]
        aux_en.θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_en.area[k]
        tendencies_pr.q_rai[k] += mph.qr_tendency * aux_en.area[k]
        tendencies_pr.q_sno[k] += mph.qs_tendency * aux_en.area[k]
    end
    return nothing
end

function quad_loop(en_thermo::SGSQuadrature, precip_model, vars, param_set, Δt::Real)
    i_ql, i_qi, i_T, i_cf, i_qt_sat, i_qt_unsat, i_T_sat, i_T_unsat = 1:8
    i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH, i_Sqr, i_Sqs = 1:8
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
    outer_src_nt = (;
        SH_qt = outer_src[i_SH_qt],
        Sqt_H = outer_src[i_Sqt_H],
        SH_H = outer_src[i_SH_H],
        Sqt_qt = outer_src[i_Sqt_qt],
        Sqt = outer_src[i_Sqt],
        SH = outer_src[i_SH],
        Sqr = outer_src[i_Sqr],
        Sqs = outer_src[i_Sqs],
    )
    outer_env_nt = (;
        ql = outer_env[i_ql],
        qi = outer_env[i_qi],
        T = outer_env[i_T],
        cf = outer_env[i_cf],
        qt_sat = outer_env[i_qt_sat],
        qt_unsat = outer_env[i_qt_unsat],
        T_sat = outer_env[i_T_sat],
        T_unsat = outer_env[i_T_unsat],
    )
    return outer_env_nt, outer_src_nt
end

function microphysics(en_thermo::SGSQuadrature, grid::Grid, state::State, precip_model, Δt::Real, param_set::APS)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    aux_en_unsat = aux_en.unsat
    aux_en_sat = aux_en.sat
    tendencies_pr = center_tendencies_precipitation(state)

    #TODO - remember you output source terms multiplied by Δt (bec. of instantaneous autoconv)
    #TODO - add tendencies for gm H, QT and QR due to rain
    #TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)

    epsilon = 10e-14 # eps(float)

    # initialize the quadrature points and their labels

    @inbounds for k in real_center_indices(grid)
        if (
            aux_en.QTvar[k] > epsilon &&
            aux_en.Hvar[k] > epsilon &&
            abs(aux_en.HQTcov[k]) > epsilon &&
            aux_en.q_tot[k] > epsilon &&
            sqrt(aux_en.QTvar[k]) < aux_en.q_tot[k]
        )
            vars = (;
                QTvar_en = aux_en.QTvar[k],
                q_tot_en = aux_en.q_tot[k],
                Hvar_en = aux_en.Hvar[k],
                θ_liq_ice_en = aux_en.θ_liq_ice[k],
                HQTcov_en = aux_en.HQTcov[k],
                area_en = aux_en.area[k],
                q_rai = prog_pr.q_rai[k],
                q_sno = prog_pr.q_sno[k],
                ρ0_c = ρ0_c[k],
                p0_c = p0_c[k],
            )
            outer_env, outer_src = quad_loop(en_thermo, precip_model, vars, param_set, Δt)

            # update environmental variables

            # update_env_precip_tendencies
            qt_tendency = outer_src.Sqt
            θ_liq_ice_tendency = outer_src.SH
            qr_tendency = outer_src.Sqr
            qs_tendency = outer_src.Sqs
            # TODO: move qt_tendency_precip_formation and θ_liq_ice_tendency_precip_formation
            # to diagnostics
            aux_en.qt_tendency_precip_formation[k] = qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = θ_liq_ice_tendency * aux_en.area[k]

            tendencies_pr.q_rai[k] += qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += qs_tendency * aux_en.area[k]

            # update cloudy/dry variables for buoyancy in TKE
            aux_en.cloud_fraction[k] = outer_env.cf
            if aux_en.cloud_fraction[k] < 1
                aux_en_unsat.q_tot[k] = outer_env.qt_unsat / (1 - aux_en.cloud_fraction[k])
                T_unsat = outer_env.T_unsat / (1 - aux_en.cloud_fraction[k])
                ts_unsat = TD.PhaseEquil_pTq(param_set, p0_c[k], T_unsat, aux_en_unsat.q_tot[k])
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts_unsat)
            else
                aux_en_unsat.q_tot[k] = 0
                aux_en_unsat.θ_dry[k] = 0
            end

            if aux_en.cloud_fraction[k] > 0
                aux_en_sat.T[k] = outer_env.T_sat / aux_en.cloud_fraction[k]
                aux_en_sat.q_tot[k] = outer_env.qt_sat / aux_en.cloud_fraction[k]
                aux_en_sat.q_vap[k] = (outer_env.qt_sat - outer_env.ql - outer_env.qi) / aux_en.cloud_fraction[k]
                ts_sat = TD.PhaseEquil_pTq(param_set, p0_c[k], aux_en_sat.T[k], aux_en_sat.q_tot[k])
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts_sat)
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts_sat)
            else
                aux_en_sat.T[k] = 0
                aux_en_sat.q_vap[k] = 0
                aux_en_sat.q_tot[k] = 0
                aux_en_sat.θ_dry[k] = 0
                aux_en_sat.θ_liq_ice[k] = 0
            end

            # update var/covar rain sources
            aux_en.Hvar_rain_dt[k] = outer_src.SH_H - outer_src.SH * aux_en.θ_liq_ice[k]
            aux_en.QTvar_rain_dt[k] = outer_src.Sqt_qt - outer_src.Sqt * aux_en.q_tot[k]
            aux_en.HQTcov_rain_dt[k] =
                outer_src.SH_qt - outer_src.SH * aux_en.q_tot[k] + outer_src.Sqt_H - outer_src.Sqt * aux_en.θ_liq_ice[k]

        else
            # if variance and covariance are zero do the same as in SA_mean
            ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])
            mph = precipitation_formation(
                param_set,
                precip_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_en.area[k],
                ρ0_c[k],
                Δt,
                ts,
            )

            # update_env_precip_tendencies
            # TODO: move qt_tendency_precip_formation and θ_liq_ice_tendency_precip_formation
            # to diagnostics
            aux_en.qt_tendency_precip_formation[k] = mph.qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_en.area[k]
            tendencies_pr.q_rai[k] += mph.qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += mph.qs_tendency * aux_en.area[k]

            # update_sat_unsat
            if TD.has_condensate(ts)
                aux_en.cloud_fraction[k] = 1
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts)
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts)
                aux_en_sat.T[k] = TD.air_temperature(ts)
                aux_en_sat.q_tot[k] = TD.total_specific_humidity(ts)
                aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(ts)
            else
                aux_en.cloud_fraction[k] = 0
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts)
                aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(ts)
                aux_en_unsat.q_tot[k] = TD.total_specific_humidity(ts)
            end

            aux_en.Hvar_rain_dt[k] = 0
            aux_en.QTvar_rain_dt[k] = 0
            aux_en.HQTcov_rain_dt[k] = 0
        end
    end

    return nothing
end
