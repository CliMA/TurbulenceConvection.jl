"""
Computes the tendency to θ_liq_ice due to water moving between
the working fluid and rain
"""
function θ_liq_ice_helper_rain(param_set::APS, ts, qr_tendency::FT) where {FT}
    Lv0 = CPP.LH_v0(param_set)
    exner_moist = TD.exner(ts)
    cp_moist = TD.cp_m(ts)
    return Lv0 / exner_moist / cp_moist * qr_tendency
end

"""
Computes the tendency to θ_liq_ice due to water moving between
the working fluid and snow
"""
function θ_liq_ice_helper_snow(param_set::APS, ts, qs_tendency::FT) where {FT}
    Ls0 = CPP.LH_s0(param_set)
    exner_moist = TD.exner(ts)
    cp_moist = TD.cp_m(ts)
    return Ls0 / exner_moist / cp_moist * qs_tendency
end

"""
Computes the tendency to θ_liq_ice due to changes in T from snow melting into rain
(heating sources outside of the working fluid)
"""
function θ_liq_ice_helper_snow_melt(param_set::APS, ts, qs_tendency::FT) where {FT}
    exner_moist = TD.exner(ts)
    cp_moist = TD.cp_m(ts)
    Lf0 = CPP.LH_f0(param_set)
    T_tendency = qs_tendency * Lf0 / cp_moist
    return T_tendency / exner_moist
end

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
"""
function precipitation_formation(
    param_set::APS,
    precipitation_model,
    qr::FT,
    qs::FT,
    area::FT,
    ρ0::FT,
    dt::FT,
    ts,
) where {FT}

    qr_tendency::FT = 0.0
    qs_tendency::FT = 0.0
    qt_tendency::FT = 0.0
    θ_liq_ice_tendency::FT = 0.0
    tmp::FT = 0.0

    if area > 0.0

        q = TD.PhasePartition(ts)

        if precipitation_model == "cutoff"
            qsat::FT = TD.q_vap_saturation(ts)

            # remove cloud condensate greater than a threshold
            tmp = min((q.liq + q.ice) / dt, -CM0.remove_precipitation(param_set, q, qsat))
            qt_tendnecy = -tmp
            θ_liq_ice_tendency = θ_liq_ice_helper_rain(param_set, ts, tmp)
        end

        if precipitation_model == "clima_1m"
            T::FT = TD.air_temperature(ts) # TODO -does it work for quadratures?
            _T_freeze::FT = CPP.T_freeze(param_set)
            _L_f::FT = TD.latent_heat_fusion(ts) # TODO -should I use Lf_0  to be consistent with the theta_liq_ice definition?
            _cv_l::FT = CPP.cv_l(param_set)

            # autoconversion cloud water to rain
            tmp = min(q.liq / dt, CM1.conv_q_liq_to_q_rai(param_set, q.liq))
            qr_tendency += tmp
            qt_tendency -= tmp
            θ_liq_ice_tendency += θ_liq_ice_helper_rain(param_set, ts, tmp)

            # autoconversion cloud ice to snow
            tmp = min(q.ice / dt, CM1.conv_q_ice_to_q_sno(param_set, q, ρ0, T))
            qs_tendency += tmp
            qt_tendency -= tmp
            θ_liq_ice_tendency += θ_liq_ice_helper_snow(param_set, ts, tmp)

            # accretion cloud water + rain
            tmp = min(q.liq / dt, CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ0))
            qr_tendency += tmp
            qt_tendency -= tmp
            θ_liq_ice_tendency += θ_liq_ice_helper_rain(param_set, ts, tmp)

            # accretion cloud ice + snow
            tmp = min(q.ice / dt, CM1.accretion(param_set, ice_type, snow_type, q.ice, qs, ρ0))
            qs_tendency += tmp
            qt_tendency -= tmp
            θ_liq_ice_tendency += θ_liq_ice_helper_snow(param_set, ts, tmp)

            # sink of cloud water via accretion cloud water + snow
            tmp = min(q.liq / dt, CM1.accretion(param_set, liq_type, snow_type, q.liq, qs, ρ0))
            if T < _T_freeze # cloud droplets freeze to become snow)
                qs_tendency += tmp
                qt_tendency -= tmp
                θ_liq_ice_tendency += θ_liq_ice_helper_snow(param_set, ts, tmp)
            else # snow melts, both cloud water and snow become rain
                α::FT = _cv_l / _L_f * (T - _T_freeze)
                qt_tendency -= tmp
                qs_tendency -= tmp * α
                qr_tendency += tmp * (1 + α)
                # TODO - check
                θ_liq_ice_tendency += θ_liq_ice_helper_rain(param_set, ts, tmp * (1 + α))
                θ_liq_ice_tendency -= θ_liq_ice_helper_snow(param_set, ts, α * tmp)
            end

            # sink of cloud ice via accretion cloud ice - rain
            tmp1::FT = min(q.ice / dt, CM1.accretion(param_set, ice_type, rain_type, q.ice, qr, ρ0))
            # sink of rain via accretion cloud ice - rain
            tmp2::FT = min(qr / dt, CM1.accretion_rain_sink(param_set, q.ice, qr, ρ0))
            qt_tendency -= tmp1
            qr_tendency -= tmp2
            qs_tendency += tmp1 + tmp2
            θ_liq_ice_tendency += θ_liq_ice_helper_snow(param_set, ts, tmp)
            # TODO
            #S_e += tmp2 * _L_f

            # accretion rain - snow
            if T < _T_freeze
                tmp = min(qr / dt, CM1.accretion_snow_rain(param_set, snow_type, rain_type, qs, qr, ρ0))
                qs_tendency += tmp
                qr_tendency -= tmp
                # TODO
                #S_e += tmp * _L_f
            else
                tmp = min(qs / dt, CM1.accretion_snow_rain(param_set, rain_type, snow_type, qr, qs, ρ0))
                qs_tendency -= tmp
                qr_tendency += tmp
                # TODO
                #S_e -= tmp * _L_f
            end
        end
    end

    return PrecipFormation(θ_liq_ice_tendency, qt_tendency, qr_tendency, qs_tendency)
end
