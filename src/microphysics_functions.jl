
"""
Source term for thetal because of qr transitioning between the working fluid and rain
(simple version to avoid exponents)
"""
function rain_source_to_thetal(p0, T, qr)
    return latent_heat(T) * qr / exner_c(p0) / cpd
end

"""
Source term for thetal because of qr transitioning between the working fluid and rain
(more detailed version, but still ignoring dqt/dqr)
"""
function rain_source_to_thetal_detailed(p0, T, qt, ql, qr)
    L = latent_heat(T)

    old_source = L * qr / exner_c(p0) / cpd

    new_source = old_source / (1 - qt) * exp(-L * ql / T / cpd / (1 - qt))

    return new_source
end

# instantly convert all cloud water exceeding a threshold to rain water
# the threshold is specified as axcess saturation
# rain water is immediately removed from the domain
function acnv_instant(ql, qt, T, p0)

    psat = pv_star(T)
    qsat = qv_star_c(p0, qt, psat)

    return fmax(0.0, ql - max_supersaturation * qsat)
end
# CLIMA microphysics rates
function terminal_velocity_single_drop_coeff(rho)

    return sqrt(8 / 3 / C_drag * (rho_cloud_liq / rho - 1))
end

function terminal_velocity(q_rai, rho)

    v_c = terminal_velocity_single_drop_coeff(rho)
    gamma_9_2 = 11.631728396567448

    term_vel = 0.0

    if q_rai > 0
        lambda_param = (8.0 * π * rho_cloud_liq * MP_n_0 / rho / q_rai)^0.25
        term_vel = gamma_9_2 * v_c / 6.0 * sqrt(g / lambda_param)
    end

    return term_vel
end

function conv_q_vap_to_q_liq(q_sat_liq, q_liq)

    return (q_sat_liq - q_liq) / tau_cond_evap
end

function conv_q_liq_to_q_rai_acnv(q_liq)

    return fmax(0.0, q_liq - q_liq_threshold) / tau_acnv
end

function conv_q_liq_to_q_rai_accr(q_liq, q_rai, rho)

    v_c = terminal_velocity_single_drop_coeff(rho)
    gamma_7_2 = 3.3233509704478426

    accr_coeff = gamma_7_2 * 8^(-7 / 8) * π^(1 / 8) * v_c * E_col * (rho / rho_cloud_liq)^(7 / 8)

    return accr_coeff * MP_n_0^(1 / 8) * sqrt(g) * q_liq * q_rai^(7 / 8)
end

function conv_q_rai_to_q_vap(q_rai, q_tot, q_liq, T, p, rho)

    L = latent_heat(T)
    gamma_11_4 = 1.6083594219855457
    v_c = terminal_velocity_single_drop_coeff(rho)
    N_Sc = nu_air / D_vapor

    av_param = sqrt(2.0 * π) * a_vent * sqrt(rho / rho_cloud_liq)
    bv_param =
        2^(7 / 16) * gamma_11_4 * π^(5 / 16) * b_vent * N_Sc^(1 / 3) * sqrt(v_c) * (rho / rho_cloud_liq)^(11 / 16)

    p_vs = pv_star(T)
    qv_sat = qv_star_c(p, q_tot, p_vs)
    q_v = q_tot - q_liq
    S = q_v / qv_sat - 1

    G_param = 1 / (L / K_therm / T * (L / Rv / T - 1.0) + Rv * T / D_vapor / p_vs)

    F_param = av_param * sqrt(q_rai) + bv_param * g^0.25 / MP_n_0^(3 / 16) / sqrt(nu_air) * q_rai^(11 / 16)

    return S * F_param * G_param * sqrt(MP_n_0) / rho
end

"""
compute the autoconversion and accretion rate
return
  new values: qt, ql, qv, thl, th, alpha
  rates: qr_src, thl_rain_src
"""
function microphysics_rain_src(rain_model, qt, ql, qr, area, T, p0, rho, dt)

    # TODO assumes no ice
    _ret = mph_struct(0, 0, 0, 0, 0, 0, 0, 0, 0)
    _ret.qv = qt - ql
    _ret.thl = t_to_thetali_c(p0, T, qt, ql, 0.0)
    _ret.th = theta_c(p0, T)
    _ret.rho = rho_c(p0, T, qt, _ret.qv)

    #TODO - temporary way to handle different autoconversion rates
    # cython doesn"t allow for string comparison without gil
    tmp_clima_acnv_flag = false
    tmp_cutoff_acnv_flag = false
    tmp_no_acnv_flag = false
    if rain_model == "clima_1m"
        tmp_clima_acnv_flag = true
    elseif rain_model == "cutoff"
        tmp_cutoff_acnv_flag = true
    elseif rain_model == "None"
        tmp_no_acnv_flag = true
    else
        error("rain model not recognized")
    end

    if area > 0.0
        if tmp_clima_acnv_flag
            _ret.qr_src = fmin(ql, (conv_q_liq_to_q_rai_acnv(ql) + conv_q_liq_to_q_rai_accr(ql, qr, rho)) * dt)
        end

        if tmp_cutoff_acnv_flag
            _ret.qr_src = fmin(ql, acnv_instant(ql, qt, T, p0))
        end

        if tmp_no_acnv_flag
            _ret.qr_src = 0.0
        end

        _ret.thl_rain_src = rain_source_to_thetal(p0, T, _ret.qr_src)

    else
        _ret.qr_src = 0.0
        _ret.thl_rain_src = 0.0
    end

    _ret.qt = qt - _ret.qr_src
    _ret.ql = ql - _ret.qr_src

    _ret.thl += _ret.thl_rain_src

    return _ret
end

"""
Source terams for rain and rain area
assuming constant rain area fraction of 1
"""
function rain_area(source_area, source_qr, current_area, current_qr)
    _ret = rain_struct()

    if source_qr <= 0.0
        _ret.qr = current_qr
        _ret.ar = current_area
    else
        _ret.qr = current_qr + source_area * source_qr
        _ret.ar = 1.0
    end

    # sketch of what to do for prognostic rain area fraction

    #function a_big, q_big, a_sml, q_sml
    #function a_const = 0.2
    #function eps     = 1e-5

    #if source_qr ==  0.
    #    _ret.qr = current_qr
    #    _ret.ar = current_area
    #else
    #    if current_area != 0.
    #        if current_area >= source_area
    #            a_big = current_area
    #            q_big = current_qr
    #            a_sml = source_area
    #            q_sml = source_qr
    #        else
    #            a_sml = current_area
    #            q_sml = current_qr
    #            a_big = source_area
    #            q_big = source_qr

    #        _ret.qr = q_big + a_sml / a_big * q_sml
    #        _ret.ar = a_big

    #    else
    #        _ret.qr = source_qr
    #        _ret.ar = source_area

    return _ret
end
