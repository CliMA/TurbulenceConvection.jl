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
    use_fallback_tendency_limiters::Bool,
)
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    aux_gm = center_aux_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_bulk = center_aux_bulk(state)
    prog_gm = center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    F2Cw::CCO.InterpolateF2C = CCO.InterpolateF2C(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0))) # shouldnt need bcs for interior interp right?
    nonequilibrium_moisture_scheme = edmf.moisture_model.scheme # my new addition
    # moisture_sources_limiter = edmf.tendency_limiters.moisture_sources_limiter
    moisture_sources_limiter = get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources), use_fallback_tendency_limiters)

    ε = eps(FT) # see EDMF_Environment.jl for explanation

    @inbounds for i in 1:N_up
        # ts_LCL = cloud_base(aux_up[i], grid, TD.PhaseNonEquil_pTq.(thermo_params, p_c, aux_up[i].T, TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice)), :up)[:cloud_base_ts] # cloud base, only keep the thermodynamic state part # deprecated for now
        w::CC.Fields.Field = F2Cw.(aux_up_f[i].w) # fill correct type  w/ interpolated w values

        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            ts_up = TD.PhaseNonEquil_pTq(thermo_params, p_c[k], T_up, q_up)

            # condensation/evaporation, deposition/sublimation

            zc = FT(grid.zc[k].z)

            q_vap_sat_liq = aux_up[i].q_vap_sat_liq[k]
            q_vap_sat_ice = aux_up[i].q_vap_sat_ice[k]
            mph = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, aux_up[i].area[k], ρ_c[k], Δt + ε, ts_up, w[k], zc, q_vap_sat_liq, q_vap_sat_ice) # ts_LCL = nothing)
            aux_up[i].ql_tendency_noneq[k] = mph.ql_tendency * aux_up[i].area[k]
            aux_up[i].qi_tendency_noneq[k] = mph.qi_tendency * aux_up[i].area[k]

            aux_up[i].ql_tendency_cond_evap[k] = aux_up[i].ql_tendency_noneq[k] # for storage
            aux_up[i].qi_tendency_sub_dep[k] = aux_up[i].qi_tendency_noneq[k] # for storage
        end
    end
    @inbounds for k in real_center_indices(grid)
        aux_bulk.ql_tendency_noneq[k] = FT(0)
        aux_bulk.qi_tendency_noneq[k] = FT(0)

        aux_bulk.ql_tendency_cond_evap[k] = FT(0) # for storage
        aux_bulk.qi_tendency_sub_dep[k] = FT(0) # for storage
        
        @inbounds for i in 1:N_up
            aux_bulk.ql_tendency_noneq[k] += aux_up[i].ql_tendency_noneq[k]
            aux_bulk.qi_tendency_noneq[k] += aux_up[i].qi_tendency_noneq[k]

            aux_bulk.ql_tendency_cond_evap[k] += aux_up[i].ql_tendency_cond_evap[k] # for storage
            aux_bulk.qi_tendency_sub_dep[k] += aux_up[i].qi_tendency_sub_dep[k] # for storage
        end
    end
    return nothing
end

"""
Computes tendencies to q_liq and q_ice due to...
"""
function compute_other_microphysics_tendencies!(grid::Grid, state::State, edmf::EDMFModel, Δt::Real, param_set::APS, use_fallback_tendency_limiters::Bool)
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    aux_gm = center_aux_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_bulk = center_aux_bulk(state)
    prog_gm = center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    F2Cw::CCO.InterpolateF2C = CCO.InterpolateF2C(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    local w::CC.Fields.Field

    # mstl = edmf.tendency_limiters.moisture_sources_limiter
    mstl = get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources), use_fallback_tendency_limiters)


    @inbounds for i in 1:N_up
    #     ts_LCL = cloud_base(
    #         aux_up[i],
    #         grid,
    #         TD.PhaseNonEquil_pTq.(
    #             thermo_params,
    #             p_c,
    #             aux_up[i].T,
    #             TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice),
    #         ),
    #         :up,
    #     )[:cloud_base_ts] # cloud base, only keep the thermodynamic state part # deprecate for now
        w = F2Cw.(aux_up_f[i].w) # fill correct type  w/ interpolated w values
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            ts_up = TD.PhaseNonEquil_pTq(thermo_params, p_c[k], T_up, q_up)

            zc = FT(grid.zc[k].z)

            S_ql_from_before::FT = iszero(aux_up[i].area[k]) ? 0 : aux_up[i].ql_tendency_noneq[k] / aux_up[i].area[k]
            S_qi_from_before::FT = iszero(aux_up[i].area[k]) ? 0 : aux_up[i].qi_tendency_noneq[k] / aux_up[i].area[k] # is it possible we lost some info here if area won't be 0 after? i dont think so cause you multiply byy 0 in the end anyway
            mph_other = other_microphysics_processes(
                param_set,
                edmf.moisture_model.heterogeneous_ice_nucleation,
                mstl,
                aux_up[i].area[k],
                ρ_c[k],
                Δt,
                ts_up,
                w[k],
                zc,
                S_ql_from_before,
                S_qi_from_before;
                # ts_LCL = nothing,
            ) # new way w/ no depdnedence on sources from anywhere else...

            aux_up[i].ql_tendency_noneq[k] += mph_other.ql_tendency * aux_up[i].area[k] # add to existing tendency (it's called 2nd to compute_nonequilibrium_moisture_tendencies!)
            aux_up[i].qi_tendency_noneq[k] += mph_other.qi_tendency * aux_up[i].area[k] # add to existing tendency (it's called 2nd to compute_nonequilibrium_moisture_tendencies!)
            aux_bulk.ql_tendency_noneq[k] += mph_other.ql_tendency * aux_up[i].area[k] # add here bc we're not gonna zero out again (we zero out in compute_nonequilibrium_moisture_tendencies!() )
            aux_bulk.qi_tendency_noneq[k] += mph_other.qi_tendency * aux_up[i].area[k] # add here bc we're not gonna zero out again (we zero out in compute_nonequilibrium_moisture_tendencies!() )


            aux_up[i].qi_tendency_hom_frz[k] = mph_other.qi_tendency_homogeneous_freezing * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_hom_frz[k] += mph_other.qi_tendency_homogeneous_freezing * aux_up[i].area[k] # for storage

            aux_up[i].qi_tendency_het_frz[k] = mph_other.qi_tendency_heterogeneous_freezing * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_het_frz[k] += mph_other.qi_tendency_heterogeneous_freezing * aux_up[i].area[k] # for storage

            aux_up[i].qi_tendency_het_nuc[k] = mph_other.qi_tendency_heterogeneous_icenuc * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_het_nuc[k] += mph_other.qi_tendency_heterogeneous_icenuc * aux_up[i].area[k] # for storage

            aux_up[i].qi_tendency_mlt[k] = mph_other.qi_tendency_melting * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_mlt[k] += mph_other.qi_tendency_melting * aux_up[i].area[k] # for storage



        end
    end

    return nothing
end



"""
Sedimentation for cloud condensate handled here...
"""
function compute_cloud_condensate_sedimentation_tendencies!(
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
    aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_bulk = center_aux_bulk(state)
    prog_gm = center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    # local liq_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type}
    # local ice_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type}
    # local liq_Dmax::FT
    # local ice_Dmax::FT
    # local scaling_factor::FT
    # local c_1::FT
    # local c_2::FT

    # F2Cw::CCO.InterpolateF2C = CCO.InterpolateF2C(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    # local w::CC.Fields.Field
    @inbounds for i in 1:N_up
        # ======================================================================== #
        if edmf.cloud_sedimentation_model isa CloudSedimentationModel &&  !edmf.cloud_sedimentation_model.grid_mean        
            # sedimentation (should this maybe be a grid mean tendency?)
            ts_up =
                TD.PhaseNonEquil_pTq.(
                    thermo_params,
                    p_c,
                    aux_up[i].T,
                    TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice),
                )

            # w = F2Cw.(aux_up_f[i].w) # fill correct type  w/ interpolated w values

            mph, mph_other = calculate_sedimentation_sources(
                param_set,
                ice_type,
                ρ_c,
                aux_up[i].q_ice,
                aux_up[i].term_vel_ice,
                aux_up[i].area,
                grid,
                CMP.ρ_cloud_ice(TCP.microphysics_params(param_set));
                differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
                velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme, # defined in update_aux
                grid_mean = false,
            ) # should this be a grid mean tendency?


            L_v0 = TCP.LH_v0(param_set)
            L_s0 = TCP.LH_s0(param_set)
            @inbounds for k in real_center_indices(grid)
                ql_tendency_sedimentation = FT(0)
                qi_tendency_sedimentation = mph[k].q_tendency
                qt_tendency_sedimentation = ql_tendency_sedimentation + qi_tendency_sedimentation
                aux_up[i].ql_tendency_sedimentation[k] += ql_tendency_sedimentation
                aux_up[i].qi_tendency_sedimentation[k] += qi_tendency_sedimentation
                aux_up[i].qt_tendency_sedimentation[k] += qt_tendency_sedimentation # = not += , cause this doesnt seem to get reset every iteration? (updated in update_aux)

                aux_bulk.ql_tendency_sedimentation[k] += ql_tendency_sedimentation
                aux_bulk.qi_tendency_sedimentation[k] += qi_tendency_sedimentation
                aux_bulk.qt_tendency_sedimentation[k] += qt_tendency_sedimentation

                Π_m = TD.exner(thermo_params, ts_up[k])
                c_pm = TD.cp_m(thermo_params, ts_up[k])
                θ_liq_ice_tendency_sedimentation =
                    1 / Π_m / c_pm * (L_v0 * ql_tendency_sedimentation + L_s0 * qi_tendency_sedimentation)
                aux_up[i].θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # adapted from microphysics_coupling.jl | precipitation_formation()
                aux_bulk.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # = not += caue these don't seem to get reset every iteration?

                # sedimentation loss into environment
                if aux_bulk[k].area < 1 # we have an environment
                    aux_en = center_aux_environment(state)
                    ql_tendency_sedimentation_other = FT(0)
                    qi_tendency_sedimentation_other = mph_other[k].q_tendency
                    qt_tendency_sedimentation_other = ql_tendency_sedimentation_other + qi_tendency_sedimentation_other
                    θ_liq_ice_tendency_sedimentation_other =
                        1 / Π_m / c_pm *
                        (L_v0 * ql_tendency_sedimentation_other + L_s0 * qi_tendency_sedimentation_other)
                    aux_en.qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                    aux_en.ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                    aux_en.qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                    aux_en.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other
                end
            end
        end
        # ======================================================================== #
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
    cloud_sedimentation_model::AbstractCloudSedimentationModel,
    rain_formation_model::AbstractRainFormationModel,
    snow_formation_model::AbstractSnowFormationModel,
    Δt::Real,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
)
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_tc = center_aux_turbconv(state)
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
                cloud_sedimentation_model,
                rain_formation_model,
                snow_formation_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_up[i].q_liq[k],
                aux_up[i].q_ice[k],
                aux_up[i].N_i[k],
                aux_up[i].term_vel_ice[k],
                aux_tc.term_vel_rain[k],
                aux_tc.term_vel_snow[k],
                aux_up[i].area[k],
                ρ_c[k],
                Δt,
                ts_up,
                precip_fraction,
                # edmf.tendency_limiters.precipitation_tendency_limiter,
                get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
            )
            aux_up[i].qt_tendency_precip_formation[k] = mph.qt_tendency * aux_up[i].area[k]
            aux_up[i].θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_up[i].area[k]
            if edmf.moisture_model isa NonEquilibriumMoisture
                aux_up[i].ql_tendency_precip_formation[k] = mph.ql_tendency * aux_up[i].area[k]
                aux_up[i].qi_tendency_precip_formation[k] = mph.qi_tendency * aux_up[i].area[k]
            end
            tendencies_pr.q_rai[k] += mph.qr_tendency * aux_up[i].area[k]
            tendencies_pr.q_sno[k] += mph.qs_tendency * aux_up[i].area[k]

            # store autoconversion and accretion for diagnostics (doens't mean much for Eq since liq/ice get set by sat adjust and T...)
            aux_up[i].ql_tendency_acnv[k] = mph.ql_tendency_acnv * aux_up[i].area[k]
            aux_up[i].qi_tendency_acnv[k] = mph.qi_tendency_acnv * aux_up[i].area[k]
            aux_up[i].ql_tendency_accr_liq_rai[k] = mph.ql_tendency_accr_liq_rai * aux_up[i].area[k]
            aux_up[i].ql_tendency_accr_liq_ice[k] = mph.ql_tendency_accr_liq_ice * aux_up[i].area[k]
            aux_up[i].ql_tendency_accr_liq_sno[k] = mph.ql_tendency_accr_liq_sno * aux_up[i].area[k]
            aux_up[i].qi_tendency_accr_ice_liq[k] = mph.qi_tendency_accr_ice_liq * aux_up[i].area[k]
            aux_up[i].qi_tendency_accr_ice_rai[k] = mph.qi_tendency_accr_ice_rai * aux_up[i].area[k]
            aux_up[i].qi_tendency_accr_ice_sno[k] = mph.qi_tendency_accr_ice_sno * aux_up[i].area[k]

        end
    end
    # TODO - to be deleted once we sum all tendencies elsewhere
    @inbounds for k in real_center_indices(grid)
        aux_bulk.θ_liq_ice_tendency_precip_formation[k] = 0
        aux_bulk.qt_tendency_precip_formation[k] = 0
        @inbounds for i in 1:N_up
            aux_bulk.θ_liq_ice_tendency_precip_formation[k] += aux_up[i].θ_liq_ice_tendency_precip_formation[k]
            aux_bulk.qt_tendency_precip_formation[k] += aux_up[i].qt_tendency_precip_formation[k]

            aux_bulk.ql_tendency_acnv[k] += aux_up[i].ql_tendency_acnv[k] # storage
            aux_bulk.qi_tendency_acnv[k] += aux_up[i].qi_tendency_acnv[k] # storage
            aux_bulk.ql_tendency_accr_liq_rai[k] += aux_up[i].ql_tendency_accr_liq_rai[k] # storage
            aux_bulk.ql_tendency_accr_liq_ice[k] += aux_up[i].ql_tendency_accr_liq_ice[k] # storage
            aux_bulk.ql_tendency_accr_liq_sno[k] += aux_up[i].ql_tendency_accr_liq_sno[k] # storage
            aux_bulk.qi_tendency_accr_ice_liq[k] += aux_up[i].qi_tendency_accr_ice_liq[k] # storage
            aux_bulk.qi_tendency_accr_ice_rai[k] += aux_up[i].qi_tendency_accr_ice_rai[k] # storage
            aux_bulk.qi_tendency_accr_ice_sno[k] += aux_up[i].qi_tendency_accr_ice_sno[k] # storage

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
