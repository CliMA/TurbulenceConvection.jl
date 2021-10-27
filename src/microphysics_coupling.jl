"""
Computes the tendency to θ_liq_ice due to qt moving between
the working fluid and precipitation
"""
function θ_liq_ice_helper(ts, qt_tendency::FT) where {FT}
    L = TD.latent_heat_vapor(ts)
    Π = TD.exner(ts)
    return -L * qt_tendency / Π / FT(CPP.cp_d(ts.param_set))
end

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
"""
function precipitation_formation(param_set::APS, precipitation_model, qr, qs, area, ρ0, dt, ts)

    qr_tendency = 0.0

    if area > 0.0

        q = TD.PhasePartition(ts)

        if precipitation_model == "clima_1m"
            qr_tendency = min(
                q.liq / dt,
                (
                    CM1.conv_q_liq_to_q_rai(param_set, q.liq) +
                    CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ0)
                ),
            )
        end
        if precipitation_model == "cutoff"
            qsat = TD.q_vap_saturation(ts)
            qr_tendency = min((q.liq + q.ice) / dt, -CM0.remove_precipitation(param_set, q, qsat))
        end
    end

    qt_tendency = -qr_tendency
    θ_liq_ice_tendency = θ_liq_ice_helper(ts, qt_tendency)
    qs_tendency = 0.0

    return PrecipFormation(θ_liq_ice_tendency, qt_tendency, qr_tendency, qs_tendency)
end
