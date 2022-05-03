"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function noneq_moisture_sources(param_set::APS, area::FT, ρ0::FT, Δt::Real, ts) where {FT}

    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    if area > 0

        q = TD.PhasePartition(param_set, ts)
        T = TD.air_temperature(param_set, ts)
        q_vap = TD.vapor_specific_humidity(param_set, ts)

        # TODO - is that the state we want to be relaxing to?
        ts_eq = TD.PhaseEquil_ρTq(param_set, ρ0, T, q.tot)
        q_eq = TD.PhasePartition(param_set, ts_eq)

        S_ql = CM1.conv_q_vap_to_q_liq_ice(param_set, liq_type, q_eq, q)
        S_qi = CM1.conv_q_vap_to_q_liq_ice(param_set, ice_type, q_eq, q)

        # TODO - handle limiters elswhere
        if S_ql >= FT(0)
            S_ql = min(S_ql, q_vap / Δt)
        else
            S_ql = -min(-S_ql, q.liq / Δt)
        end
        if S_qi >= FT(0)
            S_qi = min(S_qi, q_vap / Δt)
        else
            S_qi = -min(-S_qi, q.ice / Δt)
        end

        ql_tendency += S_ql
        qi_tendency += S_qi
    end
    return NoneqMoistureSources{FT}(ql_tendency, qi_tendency)
end

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function precipitation_formation(
    param_set::APS,
    precip_model::AbstractPrecipitationModel,
    qr::FT,
    qs::FT,
    area::FT,
    ρ0::FT,
    z::FT,
    Δt::Real,
    ts,
) where {FT}

    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    qt_tendency = FT(0)
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    qr_tendency = FT(0)
    qs_tendency = FT(0)
    e_tot_tendency = FT(0)

    if area > 0

        q = TD.PhasePartition(param_set, ts)

        I_l = TD.internal_energy_liquid(ts)
        I_i = TD.internal_energy_ice(ts)
        #Φ = gravitational_potential(atmos.orientation, aux) #TODO how to use it here?
        g = CPP.grav(param_set)
        Φ = g * z

        if precip_model isa CutoffPrecipitation
            qsat = TD.q_vap_saturation(param_set, ts)
            λ = TD.liquid_fraction(param_set, ts)

            S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(param_set, q, qsat))

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1 - λ)
            qt_tendency += S_qt
            ql_tendency += S_qt * λ
            qi_tendency += S_qt * (1 - λ)
            e_tot_tendency -= (λ * I_l + (1 - λ) * I_i + Φ) * S_qt
        end

        if precip_model isa Clima1M
            T = TD.air_temperature(param_set, ts)
            T_fr = CPP.T_freeze(param_set)
            I_d = TD.internal_energy_dry(ts)
            I_v = TD.internal_energy_vapor(ts)
            Lf = TD.latent_heat_fusion(param_set, ts)
            c_vl = CPP.cv_l(param_set)

            # Rain autoconversion
            S_qt_rain = -min(q.liq / Δt, CM1.conv_q_liq_to_q_rai(param_set, q.liq))
            qr_tendency -= S_qt_rain
            ql_tendency += S_qt_rain
            qt_tendency += S_qt_rain
            e_tot_tendency += S_qt_rain * (I_l + Φ)

            # Autoconversion of cloud ice to snow is done with a simplified rate.
            # The saturation adjustment scheme prevents using the
            # 1-moment snow autoconversion rate that assumes
            # that the supersaturation is present in the domain.
            S_qt_snow = -min(q.ice / Δt, CM1.conv_q_ice_to_q_sno_no_supersat(param_set, q.ice))
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_snow
            qi_tendency += S_qt_snow
            e_tot_tendency += S_qt_snow * (I_i + Φ)

            # accretion cloud water + rain
            S_qr = min(q.liq / Δt, CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ0))
            qr_tendency += S_qr
            qt_tendency -= S_qr
            ql_tendency -= S_qr
            e_tot_tendency -= S_qr * (I_l + Φ)

            # accretion cloud ice + snow
            S_qs = min(q.ice / Δt, CM1.accretion(param_set, ice_type, snow_type, q.ice, qs, ρ0))
            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            e_tot_tendency -= S_qs * (I_i + Φ)

            # sink of cloud water via accretion cloud water + snow
            S_qt = -min(q.liq / Δt, CM1.accretion(param_set, liq_type, snow_type, q.liq, qs, ρ0))
            if T < T_fr # cloud droplets freeze to become snow)
                qs_tendency -= S_qt
                qt_tendency += S_qt
                ql_tendency += S_qt
                e_tot_tendency += S_qt * (I_i + Φ)
            else # snow melts, both cloud water and snow become rain
                α::FT = c_vl / Lf * (T - T_fr)
                qt_tendency += S_qt
                ql_tendency += S_qt
                qs_tendency += S_qt * α
                qr_tendency -= S_qt * (1 + α)
                e_tot_tendency += S_qt * ((1 + α) * I_l - α * I_i + Φ)
            end

            # sink of cloud ice via accretion cloud ice - rain
            S_qt = -min(q.ice / Δt, CM1.accretion(param_set, ice_type, rain_type, q.ice, qr, ρ0))
            # sink of rain via accretion cloud ice - rain
            S_qr = -min(qr / Δt, CM1.accretion_rain_sink(param_set, q.ice, qr, ρ0))
            qt_tendency += S_qt
            qi_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            e_tot_tendency += S_qt * (I_i + Φ)
            e_tot_tendency -= S_qr * Lf

            # accretion rain - snow
            if T < T_fr
                S_qs = min(qr / Δt, CM1.accretion_snow_rain(param_set, snow_type, rain_type, qs, qr, ρ0))
            else
                S_qs = -min(qs / Δt, CM1.accretion_snow_rain(param_set, rain_type, snow_type, qr, qs, ρ0))
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            e_tot_tendency += S_qs * Lf
        end
    end
    return PrecipFormation{FT}(e_tot_tendency, qt_tendency, ql_tendency, qi_tendency, qr_tendency, qs_tendency)
end
