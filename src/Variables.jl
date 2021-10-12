function initialize_io(gm::GridMeanVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "u_mean")
    add_profile(Stats, "v_mean")
    add_profile(Stats, "qt_mean")
    add_profile(Stats, "RH_mean")
    add_profile(Stats, "thetal_mean")
    add_profile(Stats, "temperature_mean")
    add_profile(Stats, "buoyancy_mean")
    add_profile(Stats, "ql_mean")
    add_profile(Stats, "tke_mean")
    add_profile(Stats, "Hvar_mean")
    add_profile(Stats, "QTvar_mean")
    add_profile(Stats, "HQTcov_mean")

    add_profile(Stats, "W_third_m")
    add_profile(Stats, "H_third_m")
    add_profile(Stats, "QT_third_m")

    add_profile(Stats, "cloud_fraction_mean")

    add_ts(Stats, "lwp_mean")
    add_ts(Stats, "cloud_base_mean")
    add_ts(Stats, "cloud_top_mean")
    add_ts(Stats, "cloud_cover_mean")
    return
end

function io(gm::GridMeanVariables, Stats::NetCDFIO_Stats)
    write_profile(Stats, "u_mean", gm.U.values)
    write_profile(Stats, "v_mean", gm.V.values)
    write_profile(Stats, "qt_mean", gm.QT.values)
    write_profile(Stats, "ql_mean", gm.QL.values)
    write_profile(Stats, "temperature_mean", gm.T.values)
    write_profile(Stats, "RH_mean", gm.RH.values)
    write_profile(Stats, "buoyancy_mean", gm.B.values)
    write_profile(Stats, "thetal_mean", gm.H.values)
    write_profile(Stats, "tke_mean", gm.TKE.values)
    write_profile(Stats, "W_third_m", gm.W_third_m.values)
    write_profile(Stats, "Hvar_mean", gm.Hvar.values)
    write_profile(Stats, "QTvar_mean", gm.QTvar.values)
    write_profile(Stats, "HQTcov_mean", gm.HQTcov.values)

    write_profile(Stats, "H_third_m", gm.H_third_m.values)
    write_profile(Stats, "QT_third_m", gm.QT_third_m.values)

    write_profile(Stats, "cloud_fraction_mean", gm.cloud_fraction.values)
    write_ts(Stats, "cloud_cover_mean", gm.cloud_cover)

    write_ts(Stats, "lwp_mean", gm.lwp)
    write_ts(Stats, "cloud_base_mean", gm.cloud_base)
    write_ts(Stats, "cloud_top_mean", gm.cloud_top)
    return
end

function satadjust(gm::GridMeanVariables, grid, state, param_set)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0

    @inbounds for k in real_center_indices(grid)
        h = gm.H.values[k]
        qt = gm.QT.values[k]
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], h, qt)
        gm.QL.values[k] = TD.liquid_specific_humidity(ts)
        gm.T.values[k] = TD.air_temperature(ts)
        ρ = TD.air_density(ts)
        gm.B.values[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
        gm.RH.values[k] = TD.relative_humidity(ts)
    end
    return
end
