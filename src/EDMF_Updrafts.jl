"""
Computes tendencies to q_liq and q_ice due to
condensation, evaporation, deposition and sublimation
"""
function compute_nonequilibrium_moisture_tendencies!(
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    Δt::Real,
    param_set::APS,
)
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    aux_gm = center_aux_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    prog_gm = center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    
    @inbounds for i in 1:N_up
        ts_LCL = cloud_base(aux_up[i], grid, TD.PhaseNonEquil_pTq.(thermo_params, p_c, aux_up[i].T, TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice)), :up)[:cloud_base_ts] # cloud base, only keep the thermodynamic state part
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            ts_up = TD.PhaseNonEquil_pTq(thermo_params, p_c[k], T_up, q_up)

            # condensation/evaporation, deposition/sublimation

            # if keep this w getting, you'd want this call to be moved outside the loop
            aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
            w = CCO.InterpolateF2C(aux_up_f[i].w) # how to access?
            k_int = getfield(k,Base.propertynames(k)[1]) # get the value k.i out
            w = Base.getindex(CC.Fields.field_values(getfield(w,Base.propertynames(w)[1])), k_int) # w first field is w.bcs = getfield(w,Base.propertynames(w)[1]) , then get the index from the field value
            zc = FT(grid.zc[k].z)

            mph = noneq_moisture_sources(param_set, aux_up[i].area[k], ρ_c[k], Δt, ts_up, w, zc; ts_LCL=ts_LCL)
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
    rain_formation_model::AbstractRainFormationModel,
    Δt::Real,
    param_set::APS,
)
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    precip_fraction = compute_precip_fraction(edmf, state)

    @inbounds for i in 1:N_up
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_tot_up = aux_up[i].q_tot[k]
            if edmf.moisture_model isa EquilibriumMoisture
                ts_up = TD.PhaseEquil_pTq(thermo_params, p_c[k], T_up, q_tot_up)
            elseif edmf.moisture_model isa NonEquilibriumMoisture
                q_liq_up = aux_up[i].q_liq[k]
                q_ice_up = aux_up[i].q_ice[k]
                q = TD.PhasePartition(q_tot_up, q_liq_up, q_ice_up)
                ts_up = TD.PhaseNonEquil_pTq(thermo_params, p_c[k], T_up, q)
            else
                error(
                    "Something went wrong in EDMF_Updrafts. The expected moisture model is Equilibrium or NonEquilibrium",
                )
            end

            # autoconversion and accretion
            mph = precipitation_formation(
                param_set,
                precip_model,
                rain_formation_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_up[i].area[k],
                ρ_c[k],
                Δt,
                ts_up,
                precip_fraction,
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