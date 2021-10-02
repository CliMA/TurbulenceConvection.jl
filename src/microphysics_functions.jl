"""
It computes the tendency term to thetal because of qr transitioning between
the working fluid and rain
"""
function rain_source_to_thetal(ts, dqr_dt::FT) where {FT}
    L = TD.latent_heat_vapor(ts)
    Π = TD.exner(ts)
    return L * dqr_dt / Π / FT(CPP.cp_d(ts.param_set))
end

"""
It computes the tendencies to qr and thetal due to autoconversion and accretion

return: dqr_dt_rain, dthl_dt_rain
"""
function microphysics_rain_src(param_set::APS, rain_model, qr, area, ρ0, dt, ts)

    #TODO - temporary way to handle different autoconversion rates
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

            q = TD.PhasePartition(ts)

            qr_src = min(
                q.liq / dt,
                (
                    CM1.conv_q_liq_to_q_rai(param_set, q.liq) +
                    CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ0)
                ),
            )
        end

        if tmp_cutoff_acnv_flag

            qsat = TD.q_vap_saturation(ts)

            q = TD.PhasePartition(ts)

            qr_src = min(q.liq / dt, -CM0.remove_precipitation(param_set, q, qsat))
        end

        if tmp_no_acnv_flag
            qr_src = 0.0
        end

        # TODO add ice here
        thl_rain_src = rain_source_to_thetal(ts, qr_src)

    else
        qr_src = 0.0
        thl_rain_src = 0.0
    end
    return mph_struct(thl_rain_src, qr_src)
end

"""
Source terms for rain and rain area
assuming constant rain area fraction of 1
"""
function rain_area(source_area, source_qr, current_area, current_qr)
    if source_qr <= 0.0
        qr = current_qr
        ar = current_area
    else
        qr = current_qr + source_area * source_qr
        ar = 1.0
    end
    return rain_struct(qr, ar)
end
