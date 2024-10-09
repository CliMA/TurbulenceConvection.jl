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
        ts_LCL = cloud_base(
            aux_up[i],
            grid,
            TD.PhaseNonEquil_pTq.(
                thermo_params,
                p_c,
                aux_up[i].T,
                TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice),
            ),
            :up,
        )[:cloud_base_ts] # cloud base, only keep the thermodynamic state part
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            ts_up = TD.PhaseNonEquil_pTq(thermo_params, p_c[k], T_up, q_up)

            # condensation/evaporation, deposition/sublimation

            # if keep this w getting, you'd want this call to be moved outside the loop
            aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
            w = CCO.InterpolateF2C(aux_up_f[i].w) # how to access?
            k_int = getfield(k, Base.propertynames(k)[1]) # get the value k.i out
            w = Base.getindex(CC.Fields.field_values(getfield(w, Base.propertynames(w)[1])), k_int) # w first field is w.bcs = getfield(w,Base.propertynames(w)[1]) , then get the index from the field value (would toscalar(x::CCG.Covariant3Vector) = x.u₃ have worked?  i guess not cause it's an interp object      )
            zc = FT(grid.zc[k].z)

            mph = noneq_moisture_sources(param_set, aux_up[i].area[k], ρ_c[k], Δt, ts_up, w, zc; ts_LCL = ts_LCL)
            aux_up[i].ql_tendency_noneq[k] = mph.ql_tendency * aux_up[i].area[k]
            aux_up[i].qi_tendency_noneq[k] = mph.qi_tendency * aux_up[i].area[k]

            aux_up[i].ql_tendency_cond_evap[k] = aux_up[i].ql_tendency_noneq[k] # for storage
            aux_up[i].qi_tendency_sub_dep[k] = aux_up[i].qi_tendency_noneq[k] # for storage
        end
    end
    @inbounds for k in real_center_indices(grid)
        aux_bulk.ql_tendency_noneq[k] = 0
        aux_bulk.qi_tendency_noneq[k] = 0

        aux_bulk.ql_tendency_cond_evap[k] = 0 # for storage
        aux_bulk.qi_tendency_sub_dep[k] = 0 # for storage
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
function compute_other_microphysics_tendencies!(
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
        ts_LCL = cloud_base(
            aux_up[i],
            grid,
            TD.PhaseNonEquil_pTq.(
                thermo_params,
                p_c,
                aux_up[i].T,
                TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice),
            ),
            :up,
        )[:cloud_base_ts] # cloud base, only keep the thermodynamic state part
        @inbounds for k in real_center_indices(grid)
            T_up = aux_up[i].T[k]
            q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            ts_up = TD.PhaseNonEquil_pTq(thermo_params, p_c[k], T_up, q_up)

            # if keep this w getting, you'd want this call to be moved outside the loop
            aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
            w = CCO.InterpolateF2C(aux_up_f[i].w) # how to access?
            k_int = getfield(k, Base.propertynames(k)[1]) # get the value k.i out
            w = Base.getindex(CC.Fields.field_values(getfield(w, Base.propertynames(w)[1])), k_int) # w first field is w.bcs = getfield(w,Base.propertynames(w)[1]) , then get the index from the field value (would toscalar(x::CCG.Covariant3Vector) = x.u₃ have worked?  i guess not cause it's an interp object      )
            zc = FT(grid.zc[k].z)

            S_ql_from_before::FT = iszero(aux_up[i].area[k]) ? 0 : aux_up[i].ql_tendency_noneq[k] / aux_up[i].area[k]
            S_qi_from_before::FT = iszero(aux_up[i].area[k]) ? 0 : aux_up[i].qi_tendency_noneq[k] / aux_up[i].area[k] # is it possible we lost some info here if area won't be 0 after? i dont think so cause you multiply byy 0 in the end anyway
            mph_other = other_microphysics_processes(param_set, aux_up[i].area[k], ρ_c[k], Δt, ts_up, w, zc, S_ql_from_before, S_qi_from_before; ts_LCL = ts_LCL) # new way w/ no depdnedence on sources from anywhere else...

            aux_up[i].ql_tendency_noneq[k] += mph_other.ql_tendency * aux_up[i].area[k] # add to existing tendency (it's called 2nd to compute_nonequilibrium_moisture_tendencies!)
            aux_up[i].qi_tendency_noneq[k] += mph_other.qi_tendency * aux_up[i].area[k] # add to existing tendency (it's called 2nd to compute_nonequilibrium_moisture_tendencies!)
            aux_bulk.ql_tendency_noneq[k] += mph_other.ql_tendency * aux_up[i].area[k] # add here bc we're not gonna zero out again (we zero out in compute_nonequilibrium_moisture_tendencies!() )
            aux_bulk.qi_tendency_noneq[k] += mph_other.qi_tendency * aux_up[i].area[k] # add here bc we're not gonna zero out again (we zero out in compute_nonequilibrium_moisture_tendencies!() )

            aux_up[i].qi_tendency_het_nuc[k] = mph_other.qi_tendency_heterogeneous_icenuc *  aux_up[i].area[k] # for storage
            aux_bulk.qi_tendency_het_nuc[k] = mph_other.qi_tendency_heterogeneous_icenuc *  aux_up[i].area[k] # for storage
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
    aux_bulk = center_aux_bulk(state)
    prog_gm = center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    local liq_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type}
    local ice_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} 
    local liq_Dmax::FT
    local ice_Dmax::FT
    local scaling_factor::FT
    local c_1::FT
    local c_2::FT

    @inbounds for i in 1:N_up
        # ======================================================================== #
        if get_isbits_nt(param_set.user_args, :use_sedimentation, false) && !get_isbits_nt(param_set.user_args, :grid_mean_sedimentation, false) # drop this eventually?
            # sedimentation (should this maybe be a grid mean tendency?)
            ts_up = TD.PhaseNonEquil_pTq.(
                thermo_params,
                p_c,
                aux_up[i].T,
                TD.PhasePartition.(aux_up[i].q_tot, aux_up[i].q_liq, aux_up[i].q_ice),
            )

            aux_up_f = face_aux_updrafts(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
            w = CCO.InterpolateF2C(aux_up_f[i].w) # how to access?
            w = [ Base.getindex(CC.Fields.field_values(getfield(w, Base.propertynames(w)[1])), getfield(k, Base.propertynames(k)[1]) )  for k in real_center_indices(grid)] # w first field is w.bcs = getfield(w,Base.propertynames(w)[1]) , then get the index from the field value

            # sedimentation_liq_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_liq_number_concentration, nothing)
            # sedimentation_ice_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_ice_number_concentration, nothing)

            sedimentation_liq_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_liq_number_concentration, FT(NaN)) # testing NaN over nothing for type stability
            sedimentation_ice_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_ice_number_concentration, FT(NaN)) # testing NaN over nothing for type stability
    

            # get liquid number concentration
            if isa(sedimentation_liq_number_concentration, Number)
                N_l = FT(sedimentation_liq_number_concentration)
            elseif isa(sedimentation_liq_number_concentration, Symbol)
                # N_l = get_N_l.(param_set, sedimentation_liq_number_concentration, ts_up, w)
                N_l = get_N_l.(param_set, sedimentation_liq_number_concentration, [ts_up[k] for k in real_center_indices(grid)], w)
                # N_l = FT(NaN)
            # elseif isnothing(sedimentation_liq_number_concentration)
                # N_l = sedimentation_liq_number_concentration
            else
                error("Unsupported liquid number concentration")
            end

            # get ice number concentration
            if isa(sedimentation_ice_number_concentration, Number)
                N_i = FT(sedimentation_ice_number_concentration)
            elseif isa(sedimentation_ice_number_concentration, Symbol)
                # N_i = get_N_i.(param_set, sedimentation_ice_number_concentration, ts_up, w)
                N_i = get_N_i.(param_set, sedimentation_ice_number_concentration, [ts_up[k] for k in real_center_indices(grid)], w)
                # N_i = FT(NaN)
            # elseif isnothing(sedimentation_ice_number_concentration)
                # N_i = sedimentation_ice_number_concentration
            else
                error("Unsupported ice number concentration") 
            end



            liq_velo_scheme = get_termvel_type(get_isbits_nt(param_set.user_args, :liq_velo_scheme, :Blk1MVel)) # we dont have this anyway
            ice_velo_scheme = get_termvel_type(get_isbits_nt(param_set.user_args, :ice_velo_scheme, :Chen2022Vel)) # Blk1MVel was too fast I believe


            sedimentation_integration_method::Symbol = get_isbits_nt(param_set.user_args, :sedimentation_integration_method, :upwinding)
            liq_Dmax = get_isbits_nt(param_set.user_aux, :liq_sedimentation_Dmax, FT(Inf))
            ice_Dmax = get_isbits_nt(param_set.user_aux, :ice_sedimentation_Dmax, FT(62.5e-6)) # should this default to r_ice_snow?
            liq_sedimentation_scaling_factor = get_isbits_nt(param_set.user_aux, :liq_sedimentation_scaling_factor, FT(1.0))
            ice_sedimentation_scaling_factor = get_isbits_nt(param_set.user_aux, :ice_sedimentation_scaling_factor, FT(1.0))
            mph, mph_other = calculate_sedimentation_sources(param_set, grid, ρ_c, ts_up;
                w = w,
                area = aux_up[i].area,
                grid_mean = false,
                integration_method = sedimentation_integration_method, 
                liq_velo_scheme = liq_velo_scheme, # defined in update_aux
                ice_velo_scheme = ice_velo_scheme, # defined in update_aux
                liq_Dmax = liq_Dmax, 
                ice_Dmax = ice_Dmax,
                liq_scaling_factor = liq_sedimentation_scaling_factor,
                ice_scaling_factor = ice_sedimentation_scaling_factor,
                Nl = N_l,
                Ni = N_i,
                )

            L_v0 = TCP.LH_v0(param_set)
            L_s0 = TCP.LH_s0(param_set)
            @inbounds for k in real_center_indices(grid)
                ql_tendency_sedimentation = mph[k].ql_tendency 
                qi_tendency_sedimentation = mph[k].qi_tendency
                qt_tendency_sedimentation = ql_tendency_sedimentation + qi_tendency_sedimentation
                aux_up[i].ql_tendency_sedimentation[k] += ql_tendency_sedimentation
                aux_up[i].qi_tendency_sedimentation[k] += qi_tendency_sedimentation
                aux_up[i].qt_tendency_sedimentation[k] += qt_tendency_sedimentation # = not += , cause this doesnt seem to get reset every iteration? (updated in update_aux)

                aux_bulk.ql_tendency_sedimentation[k] += ql_tendency_sedimentation
                aux_bulk.qi_tendency_sedimentation[k] += qi_tendency_sedimentation
                aux_bulk.qt_tendency_sedimentation[k] += qt_tendency_sedimentation

                Π_m = TD.exner(thermo_params, ts_up[k])
                c_pm = TD.cp_m(thermo_params, ts_up[k])
                θ_liq_ice_tendency_sedimentation = 1 / Π_m / c_pm * ( L_v0 * ql_tendency_sedimentation + L_s0 * qi_tendency_sedimentation )
                aux_up[i].θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # adapted from microphysics_coupling.jl | precipitation_formation()
                aux_bulk.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # = not += caue these don't seem to get reset every iteration?

                # sedimentation loss into environment
                if aux_bulk[k].area < 1 # we have an environment
                    aux_en = center_aux_environment(state)
                    ql_tendency_sedimentation_other = mph_other[k].ql_tendency 
                    qi_tendency_sedimentation_other = mph_other[k].qi_tendency 
                    qt_tendency_sedimentation_other = ql_tendency_sedimentation_other + qi_tendency_sedimentation_other
                    θ_liq_ice_tendency_sedimentation_other = 1 / Π_m / c_pm * ( L_v0 * ql_tendency_sedimentation_other + L_s0 * qi_tendency_sedimentation_other )
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

             # store autoconversion and accretion for diagnostics (doens't mean much for Eq since liq/ice get set by sat adjust and T...)
             aux_up[i].ql_tendency_acnv[k] = mph.ql_tendency_acnv *  aux_up[i].area[k]
             aux_up[i].qi_tendency_acnv[k] = mph.qi_tendency_acnv *  aux_up[i].area[k]
             aux_up[i].ql_tendency_accr_liq_rai[k] = mph.ql_tendency_accr_liq_rai *  aux_up[i].area[k]
             aux_up[i].ql_tendency_accr_liq_ice[k] = mph.ql_tendency_accr_liq_ice *  aux_up[i].area[k]
             aux_up[i].ql_tendency_accr_liq_sno[k] = mph.ql_tendency_accr_liq_sno *  aux_up[i].area[k]
             aux_up[i].qi_tendency_accr_ice_liq[k] = mph.qi_tendency_accr_ice_liq *  aux_up[i].area[k]
             aux_up[i].qi_tendency_accr_ice_rai[k] = mph.qi_tendency_accr_ice_rai *  aux_up[i].area[k]
             aux_up[i].qi_tendency_accr_ice_sno[k] = mph.qi_tendency_accr_ice_sno *  aux_up[i].area[k]
           
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
