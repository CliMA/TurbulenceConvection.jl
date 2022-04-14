"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function noneq_moisture_sources(
    param_set::APS,
    area::FT,
    ρ0::FT,
    Δt::Real,
    ts,
) where {FT}

    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    if area > 0

        q = TD.PhasePartition(ts)
        T = TD.air_temperature(ts)

        q_vap = TD.vapour_specific_humidity(ts)

        # TODO - is that the state we want to be relaxing to?
        ts_eq = TD.PhaseEquil_ρTq(param_set, ρ0, T, q.tot)
        q_eq = PhasePartition(ts_eq)

        S_ql = CM1.conv_q_vap_to_q_liq_ice(param_set, liquid_type, q_eq, q)
        S_qi = CM1.conv_q_vap_to_q_liq_ice(param_set, ice_type, q_eq, q)

        # TODO - handle limiters elswhere
        if S_ql >= FT(0)
            S_ql = min(S_ql, qv / Δt)
        else
            S_ql = -min(-S_ql, ql / Δt)
        end
        if S_qi >= FT(0)
            S_qi = min(S_qi, qv / Δt)
        else
            S_qi = -min(-S_qi, qi / Δt)
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
    θ_liq_ice_tendency = FT(0)

    if area > 0

        q = TD.PhasePartition(ts)

        Π_m = TD.exner(ts)
        c_pm = TD.cp_m(ts)
        L_v0 = CPP.LH_v0(param_set)
        L_s0 = CPP.LH_s0(param_set)

        if precip_model isa CutoffPrecipitation
            qsat = TD.q_vap_saturation(ts)
            λ = TD.liquid_fraction(ts)

            S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(param_set, q, qsat))

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1 - λ)
            qt_tendency += S_qt
            ql_tendency -= S_qt * λ # are these right? I'm assuming this here ignores cross terms....
            qi_tendency -= S_qt * (1 - λ)

            θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))
        end

        if precip_model isa Clima1M
            T = TD.air_temperature(ts)
            T_fr = CPP.T_freeze(param_set)
            c_vl = CPP.cv_l(param_set)
            c_vm = TD.cv_m(ts)
            Rm = TD.gas_constant_air(ts)
            Lf = TD.latent_heat_fusion(ts)

            # Autoconversion of cloud ice to snow is done with a simplified rate.
            # The saturation adjustment scheme prevents using the
            # 1-moment snow autoconversion rate that assumes
            # that the supersaturation is present in the domain.
            S_qt_rain = -min(q.liq / Δt, CM1.conv_q_liq_to_q_rai(param_set, q.liq))
            S_qt_snow = -min(q.ice / Δt, CM1.conv_q_ice_to_q_sno_no_supersat(param_set, q.ice))
            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            ql_tendency -= S_qt_rain  # same here, are these right? it's just auto conv so i presume
            qi_tendency -= S_qt_snow
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

            # accretion cloud water + rain
            S_qr = min(q.liq / Δt, CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ0))
            qr_tendency += S_qr
            qt_tendency -= S_qr
            ql_tendency -= S_qr
            θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

            # accretion cloud ice + snow
            S_qs = min(q.ice / Δt, CM1.accretion(param_set, ice_type, snow_type, q.ice, qs, ρ0))
            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

            # sink of cloud water via accretion cloud water + snow
            S_qt = -min(q.liq / Δt, CM1.accretion(param_set, liq_type, snow_type, q.liq, qs, ρ0))
            if T < T_fr # cloud droplets freeze to become snow)
                qs_tendency -= S_qt
                qt_tendency += S_qt
                ql_tendency += S_qt
                θ_liq_ice_tendency -= S_qt / Π_m / c_pm * Lf * (1 + Rm / c_vm)
            else # snow melts, both cloud water and snow become rain (# is this melting then accretion?, if so then the factor is some sort of scavenging factor?)
                α::FT = c_vl / Lf * (T - T_fr)
                qt_tendency += S_qt
                ql_tendency += S_qt # i think then the loss in liquid is just this?
                qs_tendency += S_qt * α
                qr_tendency -= S_qt * (1 + α)
                θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1 + Rm / c_vm) * α - L_v0)
            end

            # sink of cloud ice via accretion cloud ice - rain
            S_qt = -min(q.ice / Δt, CM1.accretion(param_set, ice_type, rain_type, q.ice, qr, ρ0))
            # sink of rain via accretion cloud ice - rain
            S_qr = -min(qr / Δt, CM1.accretion_rain_sink(param_set, q.ice, qr, ρ0))
            qt_tendency += S_qt
            qi_tendency += s_qt # added for sink of cloud ice via accretion cloud ice - rain
            qr_tendency += S_qr
            qi_tendency += S_qr # added for  # sink of rain via accretion cloud ice - rain (is necessary?)
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)

            # accretion rain - snow
            if T < T_fr
                S_qs = min(qr / Δt, CM1.accretion_snow_rain(param_set, snow_type, rain_type, qs, qr, ρ0))
            else
                S_qs = -min(qs / Δt, CM1.accretion_snow_rain(param_set, rain_type, snow_type, qr, qs, ρ0))
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
        end
    end
    return PrecipFormation{FT}(θ_liq_ice_tendency, qt_tendency, ql_tendency, qi_tendency qr_tendency, qs_tendency) # to confirm, theres no way to have some sort of vectorized q_pt_tendency?
end
