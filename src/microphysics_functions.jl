# TODO - figure out which rain source to thetal we want to keep
"""
Source term for thetal because of qr transitioning between the working fluid and rain
(simple version to avoid exponents)
"""
function rain_source_to_thetal(param_set, p0::FT, T::FT, qt::FT, ql::FT, qr::FT) where {FT}
    qi = FT(0.0)
    L = TD.latent_heat_vapor(param_set, T)
    q = TD.PhasePartition(qt, ql, qi)
    Π = TD.exner_given_pressure(param_set, p0, q)
    cp_m = TD.cp_m(param_set, q)

    return L * qr / Π / cp_m
end
#"""
#Source term for thetal because of qr transitioning between the working fluid and rain
#(more detailed version, but still ignoring dqt/dqr)
#"""
#function rain_source_to_thetal_detailed(param_set, p0, T, qt, ql, qi, qr)
#    L = TD.latent_heat_vapor(param_set, T)
#    ts = TD.PhaseEquil_pTq(param_set, p0, T, qt)
#    Π = TD.exner(ts)
#    old_source = L * qr / Π / cpd
#
#    new_source = old_source / (1 - qt) * exp(-L * ql / T / cpd / (1 - qt))
#
#    return new_source
#end

"""
compute the autoconversion and accretion rate
return
  rates: qr_src, thl_rain_src
"""
function microphysics_rain_src(param_set::APS, rain_model, qt, ql, qr, area, T, p0, ρ0, dt)
    # TODO - change to pass in thermodynamics state?
    # TODO - would PhasePartition(ts) done for updrafts and environment
    # by themsleves return the correct ql and qi?
    # TODO - add ice
    qi = 0.0
    _ret = mph_struct(0, 0)

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

            q = TD.PhasePartition(qt, ql, qi)

            _ret.qr_src = min(
                q.liq / dt,
                (
                    CM1.conv_q_liq_to_q_rai(param_set, q.liq) +
                    CM1.accretion(param_set, liq_type, rai_type, q.liq, qr, ρ0)
                ),
            )
        end

        if tmp_cutoff_acnv_flag
            #qsat = TD.q_vap_saturation_generic(param_set, T, ρ0, TD.Liquid())
            q = TD.PhasePartition(qt, ql, qi)

            #_ret.qr_src = min(q.liq / dt, CM0.remove_precipitation(param_set, q, qsat))
            _ret.qr_src = min(q.liq / dt, -CM0.remove_precipitation(param_set, q))
        end

        if tmp_no_acnv_flag
            _ret.qr_src = 0.0
        end

        # TODO add ice here
        _ret.thl_rain_src = rain_source_to_thetal(param_set, p0, T, qt, ql, _ret.qr_src)

    else
        _ret.qr_src = 0.0
        _ret.thl_rain_src = 0.0
    end
    return _ret
end

"""
Source terms for rain and rain area
assuming constant rain area fraction of 1
"""
function rain_area(source_area, source_qr, current_area, current_qr)
    # TODO - get rid of it

    _ret = rain_struct()

    if source_qr <= 0.0
        _ret.qr = current_qr
        _ret.ar = current_area
    else
        _ret.qr = current_qr + source_area * source_qr
        _ret.ar = 1.0
    end
    return _ret
end
