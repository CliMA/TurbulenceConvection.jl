function initialize_io(up::UpdraftVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "updraft_cloud_fraction")

    add_ts(Stats, "updraft_cloud_cover")
    add_ts(Stats, "updraft_cloud_base")
    add_ts(Stats, "updraft_cloud_top")
    add_ts(Stats, "updraft_lwp")
    add_ts(Stats, "updraft_iwp")
    return
end

function io(up::UpdraftVariables, grid, state, Stats::NetCDFIO_Stats)
    upd_cloud_diagnostics(up, grid, state)
    write_profile(Stats, "updraft_cloud_fraction", up.cloud_fraction)
    # Note definition of cloud cover : each updraft is associated with a cloud cover equal to the maximum
    # area fraction of the updraft where ql > 0. Each updraft is assumed to have maximum overlap with respect to
    # itup (i.e. no consideration of tilting due to shear) while the updraft classes are assumed to have no overlap
    # at all. Thus total updraft cover is the sum of each updraft"s cover
    write_ts(Stats, "updraft_cloud_cover", sum(up.cloud_cover))
    write_ts(Stats, "updraft_cloud_base", minimum(abs.(up.cloud_base)))
    write_ts(Stats, "updraft_cloud_top", maximum(abs.(up.cloud_top)))
    write_ts(Stats, "updraft_lwp", up.lwp)
    write_ts(Stats, "updraft_iwp", up.iwp)
    return
end

function upd_cloud_diagnostics(up::UpdraftVariables, grid, state)
    up.lwp = 0.0
    up.iwp = 0.0

    aux_up = center_aux_updrafts(state)
    aux_up = center_aux_updrafts(state)
    ρ0_c = center_ref_state(state).ρ0
    @inbounds for i in 1:(up.n_updrafts)
        up.cloud_base[i] = zc_toa(grid)
        up.cloud_top[i] = 0.0
        up.updraft_top[i] = 0.0
        up.cloud_cover[i] = 0.0

        @inbounds for k in real_center_indices(grid)
            if aux_up[i].area[k] > 1e-3
                up.updraft_top[i] = max(up.updraft_top[i], grid.zc[k])
                up.lwp += ρ0_c[k] * aux_up[i].q_liq[k] * aux_up[i].area[k] * grid.Δz
                up.iwp += ρ0_c[k] * aux_up[i].q_ice[k] * aux_up[i].area[k] * grid.Δz

                if TD.has_condensate(aux_up[i].q_liq[k] + aux_up[i].q_ice[k])
                    up.cloud_base[i] = min(up.cloud_base[i], grid.zc[k])
                    up.cloud_top[i] = max(up.cloud_top[i], grid.zc[k])
                    up.cloud_cover[i] = max(up.cloud_cover[i], aux_up[i].area[k])
                end
            end
        end
    end
    return
end

"""
Computes tendencies to qt and θ_liq_ice due to precipitation formation
"""
function compute_precipitation_formation_tendencies(
    up_thermo::UpdraftThermodynamics,
    grid,
    state,
    up::UpdraftVariables,
    precip::PrecipVariables,
    dt,
    param_set,
)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_up = center_aux_updrafts(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    @inbounds for i in 1:(up.n_updrafts)
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_tot_up = aux_up[i].q_tot[k]
            ts_up = TD.PhaseEquil_pTq(param_set, p0_c[k], T_up, q_tot_up)

            # autoconversion and accretion
            mph = precipitation_formation(
                param_set,
                precip.precipitation_model,
                prog_pr.qr[k],
                aux_up[i].area[k],
                ρ0_c[k],
                dt,
                ts_up,
            )
            up_thermo.qt_tendency_rain_formation[i, k] = mph.qt_tendency * aux_up[i].area[k]
            up_thermo.θ_liq_ice_tendency_rain_formation[i, k] = mph.θ_liq_ice_tendency * aux_up[i].area[k]
        end
    end
    # TODO - to be deleted once we sum all tendencies elsewhere
    up_thermo.θ_liq_ice_tendency_rain_formation_tot .= up_sum(up_thermo.θ_liq_ice_tendency_rain_formation)
    up_thermo.qt_tendency_rain_formation_tot .= up_sum(up_thermo.qt_tendency_rain_formation)
    parent(tendencies_pr.qr) .+= -up_thermo.qt_tendency_rain_formation_tot
    return
end
