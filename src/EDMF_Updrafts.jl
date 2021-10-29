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
