"""
Computes tendencies to qt and θ_liq_ice due to precipitation formation
"""
function compute_precipitation_formation_tendencies(
    grid,
    state,
    up::UpdraftVariables,
    precip_model::AbstractPrecipitationModel,
    dt,
    param_set,
)
    N_up = n_updrafts(up)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    @inbounds for i in 1:N_up
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_tot_up = aux_up[i].q_tot[k]
            ts_up = TD.PhaseEquil_pTq(param_set, p0_c[k], T_up, q_tot_up)

            # autoconversion and accretion
            mph = precipitation_formation(
                param_set,
                precip_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_up[i].area[k],
                ρ0_c[k],
                dt,
                ts_up,
            )
            aux_up[i].qt_tendency_precip_formation[k] = mph.qt_tendency * aux_up[i].area[k]
            aux_up[i].θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_up[i].area[k]

            tendencies_pr.q_rai[k] += mph.qr_tendency * aux_up[i].area[k]
            tendencies_pr.q_sno[k] += mph.qs_tendency * aux_up[i].area[k]
        end
    end
    # TODO - to be deleted once we sum all tendencies elsewhere
    @inbounds for k in real_center_indices(grid)
        aux_bulk.θ_liq_ice_tendency_precip_formation[k] = 0
        aux_bulk.qt_tendency_precip_formation[k] = 0
        @inbounds for i in 1:N_up
            aux_bulk.θ_liq_ice_tendency_precip_formation[k] += aux_up[i].θ_liq_ice_tendency_precip_formation[k]
            aux_bulk.qt_tendency_precip_formation[k] += aux_up[i].qt_tendency_precip_formation[k]
        end
    end
    return nothing
end
