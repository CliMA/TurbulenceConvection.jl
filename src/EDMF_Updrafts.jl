"""
Computes tendencies to q_liq and q_ice due to
condensation, evaporation, deposition and sublimation
"""
function compute_nonequilibrium_moisture_tendencies!(
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    Δt::Real,
    param_set::EDMFPS,
)
    N_up = n_updrafts(edmf)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)

    @inbounds for i in 1:N_up
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            ts_up = TD.PhaseNonEquil_pTq(param_set.TPS, p0_c[k], T_up, q_up)

            # condensation/evaporation, deposition/sublimation
            mph = noneq_moisture_sources(param_set.PPS, aux_up[i].area[k], ρ0_c[k], Δt, ts_up)
            aux_up[i].ql_tendency_noneq[k] = mph.ql_tendency * aux_up[i].area[k]
            aux_up[i].qi_tendency_noneq[k] = mph.qi_tendency * aux_up[i].area[k]
        end
    end
    @inbounds for k in real_center_indices(grid)
        aux_bulk.ql_tendency_noneq[k] = 0
        aux_bulk.qi_tendency_noneq[k] = 0
        @inbounds for i in 1:N_up
            aux_bulk.ql_tendency_noneq[k] += aux_up[i].ql_tendency_noneq[k]
            aux_bulk.qi_tendency_noneq[k] += aux_up[i].qi_tendency_noneq[k]
        end
    end
    return nothing
end

"""
Computes tendencies to qt and θ_liq_ice due to precipitation formation
"""
function compute_precipitation_formation_tendencies(
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    precip_model::AbstractPrecipitationModel,
    Δt::Real,
    param_set::APS,
)
    N_up = n_updrafts(edmf)
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
            if edmf.moisture_model isa EquilibriumMoisture
                ts_up = TD.PhaseEquil_pTq(param_set, p0_c[k], T_up, q_tot_up)
            elseif edmf.moisture_model isa NonEquilibriumMoisture
                q_liq_up = aux_up[i].q_liq[k]
                q_ice_up = aux_up[i].q_ice[k]
                q = TD.PhasePartition(q_tot_up, q_liq_up, q_ice_up)
                ts_up = TD.PhaseNonEquil_pTq(param_set, p0_c[k], T_up, q)
            else
                error("Something went wrong in EDMF_Updrafts. The expected moisture model is Equilibrium or NonEquilibrium")
            end

            # autoconversion and accretion
            mph = precipitation_formation(
                param_set,
                precip_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_up[i].area[k],
                ρ0_c[k],
                Δt,
                ts_up,
            )
            aux_up[i].qt_tendency_precip_formation[k] = mph.qt_tendency * aux_up[i].area[k]
            aux_up[i].θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_up[i].area[k]
            if edmf.moisture_model isa NonEquilibriumMoisture
                aux_up[i].ql_tendency_precip_formation[k] = mph.ql_tendency * aux_up[i].area[k]
                aux_up[i].qi_tendency_precip_formation[k] = mph.qi_tendency * aux_up[i].area[k]
            end
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
        if edmf.moisture_model isa NonEquilibriumMoisture
            aux_bulk.ql_tendency_precip_formation[k] = 0
            aux_bulk.qi_tendency_precip_formation[k] = 0
            @inbounds for i in 1:N_up
                aux_bulk.ql_tendency_precip_formation[k] += aux_up[i].ql_tendency_precip_formation[k]
                aux_bulk.qi_tendency_precip_formation[k] += aux_up[i].qi_tendency_precip_formation[k]
            end
        end
    end
    return nothing
end
