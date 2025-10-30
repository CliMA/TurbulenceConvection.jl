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

    nonequilibrium_moisture_scheme = edmf.moisture_model.scheme # my new addition
    # moisture_sources_limiter = edmf.tendency_limiters.moisture_sources_limiter
    moisture_sources_limiter = get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources), use_fallback_tendency_limiters)

    ε = eps(FT) # see EDMF_Environment.jl for explanation

    reweight_processes_for_grid::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_processes_for_grid, false)
    reweight_extrema_only::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false)


    @inbounds for i in 1:N_up
        w::CC.Fields.Field = Ic.(aux_up_f[i].w) # fill correct type  w/ interpolated w values

        @inbounds for k in real_center_indices(grid)
            # condensation/evaporation, deposition/sublimation

            mph = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, aux_up[i].area[k], ρ_c[k], aux_up[i].p[k], aux_up[i].T[k], Δt + ε, aux_up[i].ts[k], w[k],  aux_up[i].q_vap_sat_liq[k], aux_up[i].q_vap_sat_ice[k], aux_up[i].dqvdt[k], aux_up[i].dTdt[k], aux_up[i].τ_liq[k], aux_up[i].τ_ice[k]) # ts_LCL = nothing) # we currently don't store dqvdt for the updraft...


            if reweight_processes_for_grid
                mph = reweight_noneq_moisture_sources_for_grid(k, grid, param_set, thermo_params, aux_up[i], aux_up_f[i], mph, nonequilibrium_moisture_scheme, moisture_sources_limiter, Δt, ρ_c, p_c, w, aux_up[i].dqvdt, aux_up[i].dTdt; reweight_extrema_only = reweight_extrema_only)
            end


            aux_up[i].ql_tendency_noneq[k] = mph.ql_tendency * aux_up[i].area[k]
            aux_up[i].qi_tendency_noneq[k] = mph.qi_tendency * aux_up[i].area[k]

            aux_up[i].ql_tendency_cond_evap[k] = aux_up[i].ql_tendency_noneq[k] # for storage
            aux_up[i].qi_tendency_sub_dep[k] = aux_up[i].qi_tendency_noneq[k] # for storage
        end
    end
    @inbounds for k in real_center_indices(grid)
        aux_bulk.ql_tendency_noneq[k] = FT(0)
        aux_bulk.qi_tendency_noneq[k] = FT(0)

        aux_bulk.ql_tendency_cond_evap[k] = FT(0) # for storage [ this gets zeroed in update_aux tho]
        aux_bulk.qi_tendency_sub_dep[k] = FT(0) # for storage [ this gets zeroed in update_aux tho]
        
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

    local w::CC.Fields.Field
    nonequilibrium_moisture_scheme = edmf.moisture_model.scheme # my new addition

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
        w = Ic.(aux_up_f[i].w) # fill correct type  w/ interpolated w values
        @inbounds for k in real_center_indices(grid)
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])


            S_ql_from_before::FT = iszero(aux_up[i].area[k]) ? FT(0) : aux_up[i].ql_tendency_noneq[k] / aux_up[i].area[k]
            S_qi_from_before::FT = iszero(aux_up[i].area[k]) ? FT(0) : aux_up[i].qi_tendency_noneq[k] / aux_up[i].area[k] # is it possible we lost some info here if area won't be 0 after? i dont think so cause you multiply byy 0 in the end anyway
            mph_neq_other = other_microphysics_processes(
                param_set,
                edmf.moisture_model,
                # edmf.moisture_model.heterogeneous_ice_nucleation,
                nonequilibrium_moisture_scheme,
                mstl,
                aux_up[i].area[k],
                ρ_c[k],
                aux_up[i].p[k],
                aux_up[i].T[k],
                Δt,
                aux_up[i].ts[k],
                w[k],
                S_ql_from_before,
                S_qi_from_before;
                N_l = aux_up[i].N_l[k],
                N_i = aux_up[i].N_i[k],
            ) # new way w/ no depdnedence on sources from anywhere else...

            aux_up[i].ql_tendency_noneq[k] += mph_neq_other.ql_tendency * aux_up[i].area[k] # add to existing tendency (it's called 2nd to compute_nonequilibrium_moisture_tendencies!)
            aux_up[i].qi_tendency_noneq[k] += mph_neq_other.qi_tendency * aux_up[i].area[k] # add to existing tendency (it's called 2nd to compute_nonequilibrium_moisture_tendencies!)
            aux_bulk.ql_tendency_noneq[k] += mph_neq_other.ql_tendency * aux_up[i].area[k] # add here bc we're not gonna zero out again (we zero out in compute_nonequilibrium_moisture_tendencies!() )
            aux_bulk.qi_tendency_noneq[k] += mph_neq_other.qi_tendency * aux_up[i].area[k] # add here bc we're not gonna zero out again (we zero out in compute_nonequilibrium_moisture_tendencies!() )


            aux_up[i].qi_tendency_hom_frz[k] = mph_neq_other.qi_tendency_homogeneous_freezing * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_hom_frz[k] += mph_neq_other.qi_tendency_homogeneous_freezing * aux_up[i].area[k] # for storage

            aux_up[i].qi_tendency_het_frz[k] = mph_neq_other.qi_tendency_heterogeneous_freezing * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_het_frz[k] += mph_neq_other.qi_tendency_heterogeneous_freezing * aux_up[i].area[k] # for storage

            aux_up[i].qi_tendency_het_nuc[k] = mph_neq_other.qi_tendency_heterogeneous_icenuc * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_het_nuc[k] += mph_neq_other.qi_tendency_heterogeneous_icenuc * aux_up[i].area[k] # for storage

            aux_up[i].qi_tendency_melt[k] = mph_neq_other.qi_tendency_melting * aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_melt[k] += mph_neq_other.qi_tendency_melting * aux_up[i].area[k] # for storage




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
    prog_up = center_prog_updrafts(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ


    @inbounds for i in 1:N_up
        # ======================================================================== #
        if edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean        
            # sedimentation (should this maybe be a grid mean tendency?)
            ts_up = aux_up[i].ts # this is the thermodynamic state of the updraft, which is used to compute the sedimentation sources


            wvec = CC.Geometry.WVector
            # ∇c = CCO.DivergenceF2C()
            LBF = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))


            mph_sed_ice, mph_sed_ice_other = calculate_sedimentation_sources(
                param_set,
                # ice_type,
                ρ_c,
                aux_up[i].q_ice,
                aux_up[i].term_vel_ice,
                aux_up_f[i].w, # could just be 0 bc we don't use it w/ use_relative_w = false, bc advection is handled elsewhere and we argue that diffusion is maybe ok since we don't model the size distribution
                aux_up[i].area,
                grid,
                # CMP.ρ_cloud_ice(TCP.microphysics_params(param_set));
                differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
                # velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme, # defined in update_aux
                grid_mean = false,
                use_relative_w = false, # we don't use w for sedimentation in the updrafts, so we can just set it to 0
            ) # should this be a grid mean tendency?

            # mph_sed_ice_advection_only, mph_sed_ice_other_advection_only = calculate_sedimentation_sources(
            #     param_set,
            #     # ice_type,
            #     ρ_c,
            #     aux_up[i].q_ice,
            #     aux_up[i].term_vel_ice .* 0,
            #     aux_up_f[i].w,
            #     aux_up[i].area,
            #     grid,
            #     # CMP.ρ_cloud_ice(TCP.microphysics_params(param_set));
            #     differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
            #     # velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme, # defined in update_aux
            #     grid_mean = false,
            #     use_relative_w = true, # we don't use w for sedimentation in the updrafts, so we can just set it to 0
            # ) # should this be a grid mean tendency?
            # mph_sed_ice_advection_only = @.  -∇c(wvec(LBF(Ic(aux_up_f[i].w) * prog_up[i].ρaq_ice))) / ρ_c # this is what EDMF_functions() will add
            # @. mph_sed_ice -= mph_sed_ice_advection_only
            # @. mph_sed_ice_other -= mph_sed_ice_other_advection_only


            mph_sed_liq, mph_sed_liq_other = calculate_sedimentation_sources(
                param_set,
                # liq_type,
                ρ_c,
                aux_up[i].q_liq,
                aux_up[i].term_vel_liq,
                aux_up_f[i].w, # could just be 0 bc we don't use it w/ use_relative_w = false, bc advection is handled elsewhere and we argue that diffusion is maybe ok since we don't model the size distribution
                aux_up[i].area,
                grid,
                # CMP.ρ_cloud_liq(TCP.microphysics_params(param_set));
                differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
                # velo_scheme = edmf.cloud_sedimentation_model.liq_terminal_velocity_scheme, # defined in update_aux
                grid_mean = false,
                use_relative_w = false,
            ) # should this be a grid mean tendency?

            # mph_sed_liq_advection_only, mph_sed_liq_other_advection_only = calculate_sedimentation_sources(
            #     param_set,
            #     # liq_type,
            #     ρ_c,
            #     aux_up[i].q_liq,
            #     aux_up[i].term_vel_liq .* 0,
            #     aux_up_f[i].w,
            #     aux_up[i].area,
            #     grid,
            #     # CMP.ρ_cloud_liq(TCP.microphysics_params(param_set));
            #     differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
            #     # velo_scheme = edmf.cloud_sedimentation_model.liq_terminal_velocity_scheme, # defined in update_aux
            #     grid_mean = false,
            #     use_relative_w = true, # we don't use w for sedimentation in the updrafts, so we can just set it to 0
            # ) # should this be a grid mean tendency?
            # mph_sed_liq_advection_only = @. -∇c(wvec(LBF(Ic(aux_up_f[i].w) * prog_up[i].ρaq_liq))) / ρ_c
            # @. mph_sed_liq -= mph_sed_liq_advection_only
            # @. mph_sed_liq_other -= mph_sed_liq_other_advection_only





            L_v0 = TCP.LH_v0(param_set)
            L_s0 = TCP.LH_s0(param_set)
            @inbounds for k in real_center_indices(grid)
                ql_tendency_sedimentation = mph_sed_liq[k].q_tendency
                qi_tendency_sedimentation = mph_sed_ice[k].q_tendency
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
                    ql_tendency_sedimentation_other = mph_sed_liq_other[k].q_tendency
                    qi_tendency_sedimentation_other = mph_sed_ice_other[k].q_tendency
                    qt_tendency_sedimentation_other = ql_tendency_sedimentation_other + qi_tendency_sedimentation_other
                    θ_liq_ice_tendency_sedimentation_other =
                        1 / Π_m / c_pm *
                        (L_v0 * ql_tendency_sedimentation_other + L_s0 * qi_tendency_sedimentation_other)
                    aux_en.qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                    aux_en.ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                    aux_en.qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                    aux_en.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other


                    # store contributions to other
                    aux_en.ql_tendency_sedimentation_other[k] += ql_tendency_sedimentation_other # for storage
                    aux_en.qi_tendency_sedimentation_other[k] += qi_tendency_sedimentation_other # for storage
                    aux_en.qt_tendency_sedimentation_other[k] += qt_tendency_sedimentation_other # for storage
                    aux_en.θ_liq_ice_tendency_sedimentation_other[k] += θ_liq_ice_tendency_sedimentation_other # for storage

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
    moisture_model::AbstractMoistureModel,
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
    aux_en = center_aux_environment(state)
    aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_en_f = face_aux_environment(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_tc = center_aux_turbconv(state)
    aux_bulk = center_aux_bulk(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    precip_fraction = compute_precip_fraction(edmf, state)

    @inbounds for i in 1:N_up
        if moisture_model isa NonEquilibriumMoisture
            w::CC.Fields.Field = Ic.(aux_up_f[i].w) # fill correct type w/ interpolated w values
        end
        @inbounds for k in real_center_indices(grid)

            q = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            S_i = TD.supersaturation(thermo_params, q, ρ_c[k], aux_up[i].T[k], TD.Ice())
            N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, aux_up[i].T[k], ρ_c[k], w[k])

            if moisture_model isa NonEquilibriumMoisture
                # autoconversion and accretion
                mph = precipitation_formation(
                    param_set,
                    moisture_model,
                    precip_model,
                    cloud_sedimentation_model,
                    rain_formation_model,
                    snow_formation_model,
                    edmf.area_partition_model,
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
                    aux_up[i].p[k],
                    aux_up[i].T[k],
                    Δt,
                    aux_up[i].ts[k],
                    (aux_up[i].area[k] > FT(0)) ? aux_up[i].ql_tendency_cond_evap[k] / aux_up[i].area[k] : FT(0), # this is the tendency from compute_nonequilibrium_moisture_tendencies!() above
                    (aux_up[i].area[k] > FT(0)) ? aux_up[i].qi_tendency_sub_dep[k] / aux_up[i].area[k] : FT(0),
                    (aux_up[i].area[k] > FT(0)) ? aux_up[i].qi_tendency_sedimentation[k] / aux_up[i].area[k] : FT(0),
                    precip_fraction,
                    # edmf.tendency_limiters.precipitation_tendency_limiter,
                    get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                    FT(0), # no dissipation in updrafts
                    aux_up[i].N_i[k], # no boost in updraft bc no downdraft in updraft
                    S_i,
                    aux_up[i].τ_ice[k],
                    aux_up[i].dN_i_dz[k],
                    aux_up[i].dqidz[k],
                    N_INP,
                    FT(0), # massflux doesnt count in updraft [[ we pass domain as an arg now anyway]]
                    Up,
                )

            else
                mph = precipitation_formation(
                    param_set,
                    moisture_model,
                    precip_model,
                    cloud_sedimentation_model,
                    rain_formation_model,
                    snow_formation_model,
                    edmf.area_partition_model,
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
                    aux_up[i].p[k],
                    aux_up[i].T[k],
                    Δt,
                    aux_up[i].ts[k],
                    precip_fraction,
                    # edmf.tendency_limiters.precipitation_tendency_limiter,
                    get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                    FT(0); # no dissipation in updrafts
                    N_i_no_boost = aux_up[i].N_i[k], # no boost in updraft bc no downdraft in updraft
                    qi_tendency_sed =  (aux_up[i].area[k] > FT(0)) ? aux_up[i].qi_tendency_sedimentation[k] / aux_up[i].area[k] : FT(0),
                    S_i = S_i,
                    τ_sub_dep = aux_up[i].τ_ice[k],
                    dN_i_dz = aux_up[i].dN_i_dz[k],
                    dqidz = aux_up[i].dqidz[k],
                    N_INP = N_INP,
                    massflux = FT(0), # massflux doesnt count in updraft
                    domain = Up
                )
            end



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
            aux_up[i].qi_tendency_acnv_dep[k] = mph.qi_tendency_acnv_dep * aux_up[i].area[k]
            aux_up[i].qi_tendency_acnv_dep_is[k] = mph.qi_tendency_acnv_dep_is * aux_up[i].area[k]
            aux_up[i].qi_tendency_acnv_dep_above[k] = mph.qi_tendency_acnv_dep_above * aux_up[i].area[k]
            aux_up[i].qi_tendency_acnv_agg_mix[k] = mph.qi_tendency_acnv_agg_mix * aux_up[i].area[k] # this should have already ben set in compute_nonequilibrium_moisture_tendencies() but doesn't hurt
            aux_up[i].qi_tendency_acnv_thresh[k] = mph.qi_tendency_acnv_thresh * aux_up[i].area[k] # this should have already ben set in compute_nonequilibrium_moisture_tendencies() but doesn't hurt
            aux_up[i].ql_tendency_accr_liq_rai[k] = mph.ql_tendency_accr_liq_rai * aux_up[i].area[k]
            aux_up[i].ql_tendency_accr_liq_ice[k] = mph.ql_tendency_accr_liq_ice * aux_up[i].area[k]
            aux_up[i].ql_tendency_accr_liq_sno[k] = mph.ql_tendency_accr_liq_sno * aux_up[i].area[k]
            aux_up[i].qi_tendency_accr_ice_liq[k] = mph.qi_tendency_accr_ice_liq * aux_up[i].area[k]
            aux_up[i].qi_tendency_accr_ice_rai[k] = mph.qi_tendency_accr_ice_rai * aux_up[i].area[k]
            aux_up[i].qi_tendency_accr_ice_sno[k] = mph.qi_tendency_accr_ice_sno * aux_up[i].area[k]
            #
            aux_tc.qs_tendency_accr_rai_sno[k] += mph.qs_tendency_accr_rai_sno * aux_up[i].area[k] # we calculate in aux/en for the temperature dependence but store a combined output
        end
    end
    # TODO - to be deleted once we sum all tendencies elsewhere
    @inbounds for k in real_center_indices(grid)
        aux_bulk.θ_liq_ice_tendency_precip_formation[k] = FT(0)
        aux_bulk.qt_tendency_precip_formation[k] = FT(0)
        @inbounds for i in 1:N_up
            aux_bulk.θ_liq_ice_tendency_precip_formation[k] += aux_up[i].θ_liq_ice_tendency_precip_formation[k]
            aux_bulk.qt_tendency_precip_formation[k] += aux_up[i].qt_tendency_precip_formation[k]

            aux_bulk.ql_tendency_acnv[k] += aux_up[i].ql_tendency_acnv[k] # storage
            aux_bulk.qi_tendency_acnv[k] += aux_up[i].qi_tendency_acnv[k] # storage
            aux_bulk.qi_tendency_acnv_dep[k] += aux_up[i].qi_tendency_acnv_dep[k] # storage
            aux_bulk.qi_tendency_acnv_dep_is[k] += aux_up[i].qi_tendency_acnv_dep_is[k] # storage
            aux_bulk.qi_tendency_acnv_dep_above[k] += aux_up[i].qi_tendency_acnv_dep_above[k] # storage
            aux_bulk.qi_tendency_acnv_agg_mix[k] += aux_up[i].qi_tendency_acnv_agg_mix[k] # storage
            aux_bulk.qi_tendency_acnv_thresh[k] += aux_up[i].qi_tendency_acnv_thresh[k] # storage
            aux_bulk.ql_tendency_accr_liq_rai[k] += aux_up[i].ql_tendency_accr_liq_rai[k] # storage
            aux_bulk.ql_tendency_accr_liq_ice[k] += aux_up[i].ql_tendency_accr_liq_ice[k] # storage
            aux_bulk.ql_tendency_accr_liq_sno[k] += aux_up[i].ql_tendency_accr_liq_sno[k] # storage
            aux_bulk.qi_tendency_accr_ice_liq[k] += aux_up[i].qi_tendency_accr_ice_liq[k] # storage
            aux_bulk.qi_tendency_accr_ice_rai[k] += aux_up[i].qi_tendency_accr_ice_rai[k] # storage
            aux_bulk.qi_tendency_accr_ice_sno[k] += aux_up[i].qi_tendency_accr_ice_sno[k] # storage

        end
        if edmf.moisture_model isa NonEquilibriumMoisture
            aux_bulk.ql_tendency_precip_formation[k] = FT(0)
            aux_bulk.qi_tendency_precip_formation[k] = FT(0)

            @inbounds for i in 1:N_up
                aux_bulk.ql_tendency_precip_formation[k] += aux_up[i].ql_tendency_precip_formation[k]
                aux_bulk.qi_tendency_precip_formation[k] += aux_up[i].qi_tendency_precip_formation[k]
            end
        end
    end
    return nothing
end
