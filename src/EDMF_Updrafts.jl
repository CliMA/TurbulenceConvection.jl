# TODO - move to EDMF_Precipitation
# TODO - do all three precip tendencies computations in one loop over k
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
    aux_tc = center_aux_turbconv(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    n_up = up.n_updrafts

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(n_up)
            T_up = aux_up[i].T[k]
            q_tot_up = aux_up[i].q_tot[k]
            ts_up = TD.PhaseEquil_pTq(param_set, p0_c[k], T_up, q_tot_up)

            # autoconversion and accretion
            mph = precipitation_formation(
                param_set,
                precip.precipitation_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_up[i].area[k],
                ρ0_c[k],
                dt,
                ts_up,
            )
            aux_up[i].θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_up[i].area[k]
            aux_up[i].qt_tendency_precip_formation[k] = mph.qt_tendency * aux_up[i].area[k]

            tendencies_pr.q_rai[k] += mph.qr_tendency * aux_up[i].area[k]
            tendencies_pr.q_sno[k] += mph.qs_tendency * aux_up[i].area[k]
        end
        aux_tc.bulk.θ_liq_ice_tendency_precip_formation_tot[k] =
            sum(i -> aux_up[i].θ_liq_ice_tendency_precip_formation[k], 1:n_up)
        aux_tc.bulk.qt_tendency_precip_formation_tot[k] = sum(i -> aux_up[i].qt_tendency_precip_formation[k], 1:n_up)
    end
    return
end
