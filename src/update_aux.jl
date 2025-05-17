function update_aux!(edmf::EDMFModel, grid::Grid, state::State, surf::SurfaceBase, param_set::APS, t::Real, Δt::Real, cfl_limit::Real, use_fallback_tendency_limiters::Bool)
    #####
    ##### Unpack common variables
    #####
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params = TCP.microphysics_params(param_set)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    c_m = mixing_length_params(edmf).c_m
    Le = mixing_length_params(edmf).Le
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    KQ = center_aux_turbconv(state).KQ
    obukhov_length = surf.obukhov_length
    FT = float_type(state)
    prog_gm = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_tc = center_aux_turbconv(state)
    aux_bulk = center_aux_bulk(state)
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    ρ_f = aux_gm_f.ρ
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_en_unsat = aux_en.unsat
    aux_en_sat = aux_en.sat
    wvec = CC.Geometry.WVector
    max_area = edmf.max_area
    ts_gm = center_aux_grid_mean(state).ts
    ts_env = center_aux_environment(state).ts

    massflux = aux_tc_f.massflux
    ∂M∂z = aux_tc.∂M∂z
    ∂lnM∂z = aux_tc.∂lnM∂z
    ∂w∂z = aux_tc.∂w∂z
    ∂M∂z_ρa = aux_tc.∂M∂z_ρa
    massflux_c = aux_tc.massflux
    parent(massflux) .= 0
    parent(massflux_c) .= 0
    parent(∂M∂z) .= 0
    parent(∂lnM∂z) .= 0
    parent(∂w∂z) .= 0
    parent(∂M∂z_ρa) .= 0


    prog_gm_uₕ = grid_mean_uₕ(state)
    Ic = CCO.InterpolateF2C()
    #####
    ##### center variables
    #####

    @inbounds for k in real_center_indices(grid)
        #####
        ##### Set primitive variables
        #####
        e_pot = geopotential(param_set, grid.zc.z[k])
        @inbounds for i in 1:N_up
            if prog_up[i].ρarea[k] / ρ_c[k] > edmf.minimum_area # this condition is problematic for limiters, should never rely on aux_up or aux_bulk for limiting...

                #=
                    note this goes wrong if the tracer value is 0, θ_liq_ice = 0 is definitely broken, and q_tot = 0 almost certainly is too...
                    we fall back to gm there (we fixed the filter on prog so this shouldn't come up as much now but...)
                    This doesn't imapct env as much there's no prognostic anything for env, so if the implied tracer in env is 0, prog up is unrealistically large but most of prog_up's tendencies feed into gm so env is buffered.

                    This situation would typically arise due to decoupled entr/detr between area and tracers, and is most visible at small areas where it is easy for the tracer to be depleted.
                =#
                # aux_up[i].θ_liq_ice[k] = prog_up[i].ρaθ_liq_ice[k] / prog_up[i].ρarea[k]
                # aux_up[i].q_tot[k] = prog_up[i].ρaq_tot[k] / prog_up[i].ρarea[k]
                aux_up[i].θ_liq_ice[k] = iszero(prog_up[i].ρaθ_liq_ice[k]) ? aux_gm.θ_liq_ice[k] : prog_up[i].ρaθ_liq_ice[k] / prog_up[i].ρarea[k] # if tracer is 0 but area isn't, fall back to gm
                aux_up[i].q_tot[k] = iszero(prog_up[i].ρaq_tot[k]) ? aux_gm.q_tot[k] : prog_up[i].ρaq_tot[k] / prog_up[i].ρarea[k] # if tracer is 0 but area isn't, fall back to gm
                ##
                aux_up[i].area[k] = prog_up[i].ρarea[k] / ρ_c[k]
            else
                aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
                aux_up[i].q_tot[k] = aux_gm.q_tot[k]
                aux_up[i].area[k] = 0
                # aux_up[i].area[k] = edmf.minimum_area # my test area [decided not to, area can just be 0 in aux bc we should be using prog_up for real calcs...]
            end
            thermo_args = ()
            if edmf.moisture_model isa NonEquilibriumMoisture
                if prog_up[i].ρarea[k] / ρ_c[k] > edmf.minimum_area
                    aux_up[i].q_liq[k] = prog_up[i].ρaq_liq[k] / prog_up[i].ρarea[k]
                    aux_up[i].q_ice[k] = prog_up[i].ρaq_ice[k] / prog_up[i].ρarea[k]
                else
                    aux_up[i].q_liq[k] = prog_gm.q_liq[k] # what's the point of these? If we do this should we not just set ρarea to 0
                    aux_up[i].q_ice[k] = prog_gm.q_ice[k] # what's the point of  these?
                    # if prog_up[i].ρarea[k] / ρ_c[k] > 0
                    #     aux_up[i].q_liq[k] = prog_up[i].ρaq_liq[k] / prog_up[i].ρarea[k]
                    #     aux_up[i].q_ice[k] = prog_up[i].ρaq_ice[k] / prog_up[i].ρarea[k]
                    # else
                    #     aux_up[i].q_liq[k] = prog_gm.q_liq[k]
                    #     aux_up[i].q_ice[k] = prog_gm.q_ice[k]
                    # end
                end
                thermo_args = (aux_up[i].q_liq[k], aux_up[i].q_ice[k])
            end
            # ts_up_i = thermo_state_pθq(param_set, p_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k], thermo_args...)
        end

        #####
        ##### compute bulk (aux)
        #####
        aux_bulk.q_tot[k] = 0
        aux_bulk.θ_liq_ice[k] = 0
        aux_bulk.area[k] = sum(i -> aux_up[i].area[k], 1:N_up)
        if aux_bulk.area[k] > 0
            @inbounds for i in 1:N_up
                a_k = aux_up[i].area[k]
                a_bulk_k = aux_bulk.area[k]
                aux_bulk.q_tot[k] += a_k * aux_up[i].q_tot[k] / a_bulk_k
                aux_bulk.θ_liq_ice[k] += a_k * aux_up[i].θ_liq_ice[k] / a_bulk_k
            end
        else
            aux_bulk.q_tot[k] = aux_gm.q_tot[k]
            aux_bulk.θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
        end
        if edmf.moisture_model isa NonEquilibriumMoisture
            aux_bulk.q_liq[k] = 0
            aux_bulk.q_ice[k] = 0
            if aux_bulk.area[k] > 0
                @inbounds for i in 1:N_up
                    a_k = aux_up[i].area[k]
                    a_bulk_k = aux_bulk.area[k]
                    aux_bulk.q_liq[k] += a_k * aux_up[i].q_liq[k] / a_bulk_k
                    aux_bulk.q_ice[k] += a_k * aux_up[i].q_ice[k] / a_bulk_k
                end
            else
                aux_bulk.q_liq[k] = prog_gm.q_liq[k]
                aux_bulk.q_ice[k] = prog_gm.q_ice[k]
            end
        end
        aux_en.area[k] = 1 - aux_bulk.area[k]
        aux_en.tke[k] = prog_en.ρatke[k] / (ρ_c[k] * aux_en.area[k])
        if edmf.thermo_covariance_model isa PrognosticThermoCovariances
            aux_en.Hvar[k] = prog_en.ρaHvar[k] / (ρ_c[k] * aux_en.area[k])
            aux_en.QTvar[k] = prog_en.ρaQTvar[k] / (ρ_c[k] * aux_en.area[k])
            aux_en.HQTcov[k] = prog_en.ρaHQTcov[k] / (ρ_c[k] * aux_en.area[k])
        end

        #####
        ##### decompose_environment
        #####
        a_bulk_c = aux_bulk.area[k]
        val1 = 1 / (1 - a_bulk_c)
        val2 = a_bulk_c * val1
        aux_en.q_tot[k] = max(val1 * aux_gm.q_tot[k] - val2 * aux_bulk.q_tot[k], 0) #Yair - this is here to prevent negative QT

        # if (val1 * aux_gm.q_tot[k] - val2 * aux_bulk.q_tot[k]) < 0
        #     env_bad_area = (val1 * aux_gm.q_tot[k] - val2 * aux_bulk.q_tot[k])
        #     @warn "This is gonna crash, environmental q_tot was gonna be $env_bad_area"
        # end

        aux_en.θ_liq_ice[k] = val1 * aux_gm.θ_liq_ice[k] - val2 * aux_bulk.θ_liq_ice[k]
        if edmf.moisture_model isa NonEquilibriumMoisture # otherwise it's diagnosed in thermo_state_pθq() calculation...
            aux_en.q_liq[k] = max(val1 * prog_gm.q_liq[k] - val2 * aux_bulk.q_liq[k], 0)
            aux_en.q_ice[k] = max(val1 * prog_gm.q_ice[k] - val2 * aux_bulk.q_ice[k], 0)
        end

        #=
            NOTE: We could apply a fix here to ensure that if aux_en area is greater than 0, θ_liq_ice and q_tot are not 0...
            We did the same for aux_up above and prog_up in updraft filters...
                Note though that applying the fix does change the gm, and means the prognostic updraft value is bad... 
                It was one thing to try to fix small bad excess tracer detrainment at small area, but tracer going to 0 prolly means the updraft temp blew up or something...
                It's good to consider the imapct of giving that a pass -- still, it's prolly not gonna help calibrations so the improved stability could help CalibrateEDMF.jl
        =#
        if (aux_en.area[k] > 0)
            if iszero(aux_en.θ_liq_ice[k]) # 0 temp is bad
                aux_en.θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
            end
            if iszero(aux_en.q_tot[k]) # debatably bad but qt=0 is almost certainly not good
                aux_en.q_tot[k] = aux_gm.q_tot[k]
            end
        end

        #####
        ##### condensation, etc (done via saturation_adjustment or non-equilibrium) and buoyancy
        #####
        thermo_args = if edmf.moisture_model isa EquilibriumMoisture
            ()
        elseif edmf.moisture_model isa NonEquilibriumMoisture
            (aux_en.q_liq[k], aux_en.q_ice[k])
        else
            error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
        end

        ts_env[k] = thermo_state_pθq(param_set, p_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k], thermo_args...)
        ts_en = ts_env[k]
        aux_en.T[k] = TD.air_temperature(thermo_params, ts_en)

        # @assert aux_en.T[k] > 0 "Backed out environment has negative temperature, this is bad. Status is z = $(grid.zc.z[k]), T = $(aux_en.T[k]), p = $(p_c[k]), θ = $(aux_en.θ_liq_ice[k]), q = $(aux_en.q_tot[k]),  ts = $ts_en"

        aux_en.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts_en)
        aux_en.θ_dry[k] = TD.dry_pottemp(thermo_params, ts_en)
        rho = TD.air_density(thermo_params, ts_en)
        aux_en.buoy[k] = buoyancy_c(param_set, ρ_c[k], rho)
        aux_en.RH[k] = TD.relative_humidity(thermo_params, ts_en)
    end


    #####
    #####  update prognostic bulk (these mirror the prognostic variables exactly, no minimum_area or other clippings -- filtering already happened so should be good)
    ### - these are important because while sometimes tendencies are calculated using aux, they are always applied to prog, so for limiting we need to know the exact prognostic values
    #####
    prog_bulk = center_prog_bulk(state)
    prog_bulk.ρarea .= FT(0)
    prog_bulk.ρaθ_liq_ice .= FT(0)
    prog_bulk.ρaq_tot .= FT(0)
    if edmf.moisture_model isa NonEquilibriumMoisture
        prog_bulk.ρaq_liq .= FT(0)
        prog_bulk.ρaq_ice .= FT(0)
    end
    prog_bulk_f = face_prog_bulk(state)
    prog_bulk_f.ρaw .= FT(0)

    @inbounds for i in 1:N_up
        @. prog_bulk.ρarea += prog_up[i].ρarea
        @. prog_bulk.ρaθ_liq_ice += prog_up[i].ρaθ_liq_ice
        @. prog_bulk.ρaq_tot += prog_up[i].ρaq_tot
        if edmf.moisture_model isa NonEquilibriumMoisture
            @. prog_bulk.ρaq_liq += prog_up[i].ρaq_liq
            @. prog_bulk.ρaq_ice += prog_up[i].ρaq_ice
        end
        @. prog_bulk_f.ρaw += prog_up_f[i].ρaw
    end

    # track prog_en exactly (not to be confused with aux_en which is limited by clippings on aux_up etc)
    prog_en_gm = center_prog_environment_up_gm_version(state)
    @. prog_en_gm.ρarea = ρ_c - prog_bulk.ρarea
    @. prog_en_gm.ρaθ_liq_ice = prog_gm.ρθ_liq_ice - prog_bulk.ρaθ_liq_ice
    @. prog_en_gm.ρaq_tot = prog_gm.ρq_tot - prog_bulk.ρaq_tot
    if edmf.moisture_model isa NonEquilibriumMoisture
        @. prog_en_gm.ρaq_liq = ρ_c * prog_gm.q_liq - prog_bulk.ρaq_liq
        @. prog_en_gm.ρaq_ice = ρ_c * prog_gm.q_ice - prog_bulk.ρaq_ice
    end
    prog_en_gm_f = face_prog_environment_up_gm_version(state)
    @. prog_en_gm_f.ρaw = -prog_bulk_f.ρaw

    # ======================================================== # zero out noneq tendencies before calling microphysics again (does this interfere w/ writing to disk? -- needed bc env/updraft are cross writing to each other... #
    # aux_en.ql_tendency_noneq .= FT(0) # don't zero this out  bc it seeemd to break the output writing (some order of read, calculate gm, write, zero out problem probably...)
    # aux_en.qi_tendency_noneq .= FT(0) # don't zero this out  bc it seeemd to break the output writing (some order of read, calculate gm, write, zero out problem probably...)
    #
    aux_en.ql_tendency_cond_evap .= FT(0)
    aux_en.qi_tendency_sub_dep .= FT(0)
    #
    aux_en.ql_tendency_sedimentation .= FT(0)
    aux_en.qi_tendency_sedimentation .= FT(0)
    aux_en.qt_tendency_sedimentation .= FT(0)
    aux_en.θ_liq_ice_tendency_sedimentation .= FT(0)
    #
    aux_en.ql_tendency_acnv .= FT(0)
    aux_en.qi_tendency_acnv .= FT(0)
    aux_en.qi_tendency_acnv_dep .= FT(0)
    aux_en.qi_tendency_acnv_agg .= FT(0)
    #
    aux_en.ql_tendency_accr_liq_rai .= FT(0)
    aux_en.ql_tendency_accr_liq_ice .= FT(0)
    aux_en.ql_tendency_accr_liq_sno .= FT(0)
    #
    aux_en.qi_tendency_accr_ice_liq .= FT(0)
    aux_en.qi_tendency_accr_ice_rai .= FT(0)
    aux_en.qi_tendency_accr_ice_sno .= FT(0)
    #
    aux_en.qi_tendency_hom_frz .= FT(0)
    aux_en.qi_tendency_het_frz .= FT(0)
    aux_en.qi_tendency_het_nuc .= FT(0)
    aux_en.qi_tendency_mlt .= FT(0)
    #
    # aux_en.qi_tendency_vert_adv .= FT(0) # gm only
    # aux_en.qi_tendency_ls_vert_adv .= FT(0) # gm only
    # aux_en.qi_tendency_sgs .= FT(0) # gm only




    @inbounds for i in 1:N_up
        # aux_up[i].ql_tendency_noneq .= FT(0)
        # aux_up[i].qi_tendency_noneq .= FT(0)
        #
        aux_up[i].ql_tendency_cond_evap .= FT(0)
        aux_up[i].qi_tendency_sub_dep .= FT(0)
        #
        aux_up[i].ql_tendency_sedimentation .= FT(0)
        aux_up[i].qi_tendency_sedimentation .= FT(0)
        aux_up[i].qt_tendency_sedimentation .= FT(0)
        aux_up[i].θ_liq_ice_tendency_sedimentation .= FT(0)
        #
        #
        aux_up[i].ql_tendency_acnv .= FT(0)
        aux_up[i].qi_tendency_acnv .= FT(0)
        aux_up[i].qi_tendency_acnv_dep .= FT(0)
        aux_up[i].qi_tendency_acnv_agg .= FT(0)
        #
        aux_up[i].ql_tendency_accr_liq_rai .= FT(0)
        aux_up[i].ql_tendency_accr_liq_ice .= FT(0)
        aux_up[i].ql_tendency_accr_liq_sno .= FT(0)
        #
        aux_up[i].qi_tendency_accr_ice_liq .= FT(0)
        aux_up[i].qi_tendency_accr_ice_rai .= FT(0)
        aux_up[i].qi_tendency_accr_ice_sno .= FT(0)
        #
        aux_up[i].qi_tendency_hom_frz .= FT(0)
        aux_up[i].qi_tendency_het_frz .= FT(0)
        aux_up[i].qi_tendency_het_nuc .= FT(0)
        aux_up[i].qi_tendency_mlt .= FT(0)
        #
        # aux_up[i].qi_tendency_vert_adv .= FT(0) # gm only
        # aux_up[i].qi_tendency_ls_vert_adv .= FT(0) # gm only
        # aux_up[i].qi_tendency_sgs .= FT(0) # gm only
    end
    # aux_bulk.ql_tendency_noneq .= FT(0)
    # aux_bulk.qi_tendency_noneq .= FT(0)
    #
    aux_bulk.ql_tendency_cond_evap .= FT(0)
    aux_bulk.qi_tendency_sub_dep .= FT(0)
    #
    aux_bulk.ql_tendency_sedimentation .= FT(0)
    aux_bulk.qi_tendency_sedimentation .= FT(0)
    aux_bulk.qt_tendency_sedimentation .= FT(0)
    aux_bulk.θ_liq_ice_tendency_sedimentation .= FT(0)
    #
    aux_bulk.ql_tendency_acnv .= FT(0)
    aux_bulk.qi_tendency_acnv .= FT(0)
    aux_bulk.qi_tendency_acnv_dep .= FT(0)
    aux_bulk.qi_tendency_acnv_agg .= FT(0)
    #
    aux_bulk.ql_tendency_accr_liq_rai .= FT(0)
    aux_bulk.ql_tendency_accr_liq_ice .= FT(0)
    aux_bulk.ql_tendency_accr_liq_sno .= FT(0)
    #
    aux_bulk.qi_tendency_accr_ice_liq .= FT(0)
    aux_bulk.qi_tendency_accr_ice_rai .= FT(0)
    aux_bulk.qi_tendency_accr_ice_sno .= FT(0)
    #
    aux_bulk.qi_tendency_hom_frz .= FT(0)
    aux_bulk.qi_tendency_het_frz .= FT(0)
    aux_bulk.qi_tendency_het_nuc .= FT(0)
    aux_bulk.qi_tendency_mlt .= FT(0)
    #
    # aux_bulk.qi_tendency_vert_adv .= FT(0) # gm only
    # aux_bulk.qi_tendency_ls_vert_adv .= FT(0) # gm only
    # aux_bulk.qi_tendency_sgs .= FT(0) # gm only


    # zero out our bulk tendency trackers... these are used in limit_up_tendencies!() to limit tendencies for model stability...
    tendencies_bulk = center_tendencies_bulk(state)
    tendencies_bulk.ρarea .= FT(0)
    tendencies_bulk.ρaθ_liq_ice .= FT(0)
    tendencies_bulk.ρaq_tot .= FT(0)
    if edmf.moisture_model isa NonEquilibriumMoisture
        tendencies_bulk.ρaq_liq .= FT(0)
        tendencies_bulk.ρaq_ice .= FT(0)
    end
    tendencies_bulk_f = face_tendencies_bulk(state)
    tendencies_bulk_f.ρaw .= FT(0)


    # these get overwritten technically so don't actually need zeroing out
    # tendencies_bulk_adjustments = center_tendencies_bulk_adjustments(state)
    # tendencies_bulk_adjustments.ρarea .= FT(0)
    # tendencies_bulk_adjustments.ρaθ_liq_ice .= FT(0)
    # tendencies_bulk_adjustments.ρaq_tot .= FT(0)
    # tendencies_bulk_adjustments.ρaq_liq .= FT(0)
    # tendencies_bulk_adjustments.ρaq_ice .= FT(0)
    # tendencies_bulk_adjustments_f = face_tendencies_bulk_adjustments(state)
    # tendencies_bulk_adjustments_f.ρaw .= FT(0)

    # ======================================================== # 

    # calculate and  q_vap_sat_liq, q_vap_sat_ice (maybe move to assing_thermo_aux?)
    if edmf.moisture_model isa NonEquilibriumMoisture
        # en has been updated, call before microphysics(). Use aux_en bc that's what microphysics uses
        @. aux_en.q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, aux_en.T, ρ_c, TD.Liquid())
        @. aux_en.q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, aux_en.T, ρ_c, TD.Ice())
    end

    # ======================================================== #

    microphysics(edmf.en_thermo, grid, state, edmf, edmf.moisture_model, edmf.precip_model, edmf.cloud_sedimentation_model, edmf.rain_formation_model, edmf.snow_formation_model, Δt, param_set, use_fallback_tendency_limiters) # set env tendencies for microphysics
    
    reweight_processes_for_grid::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_processes_for_grid, false)
    reweight_extrema_only::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false)

    @inbounds for k in real_center_indices(grid)
        a_bulk_c = aux_bulk.area[k]
        @inbounds for i in 1:N_up
            if aux_up[i].area[k] < edmf.minimum_area && k > kc_surf && aux_up[i].area[k - 1] > 0.0 # buoyancy hack for left biased buoyancy calc... should ensure paw, buoyancy never goes below 0.
                qt = aux_up[i].q_tot[k - 1]
                h = aux_up[i].θ_liq_ice[k - 1]
                if edmf.moisture_model isa EquilibriumMoisture
                    ts_up = thermo_state_pθq(param_set, p_c[k], h, qt)
                elseif edmf.moisture_model isa NonEquilibriumMoisture
                    ql = aux_up[i].q_liq[k - 1]
                    qi = aux_up[i].q_ice[k - 1]
                    ts_up = thermo_state_pθq(param_set, p_c[k], h, qt, ql, qi)
                else
                    error("Something went wrong. emdf.moisture_model options are equilibrium or nonequilibrium")
                end
            else
                if edmf.moisture_model isa EquilibriumMoisture
                    ts_up = thermo_state_pθq(param_set, p_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k])
                elseif edmf.moisture_model isa NonEquilibriumMoisture
                    ts_up = thermo_state_pθq(
                        param_set,
                        p_c[k],
                        aux_up[i].θ_liq_ice[k],
                        aux_up[i].q_tot[k],
                        aux_up[i].q_liq[k],
                        aux_up[i].q_ice[k],
                    )
                else
                    error("Something went wrong. emdf.moisture_model options are equilibrium or nonequilibrium")
                end
            end
            aux_up[i].q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts_up)
            aux_up[i].q_ice[k] = TD.ice_specific_humidity(thermo_params, ts_up)


            if reweight_processes_for_grid && (edmf.moisture_model isa EquilibriumMoisture) # really this needs qt to be updated at all k so you arent overreaching right?
                q = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k]) # already did the calc, let's reuse
                q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_up[i], q, p_c; reweight_extrema_only = reweight_extrema_only)
                aux_up[i].q_liq[k] = q_liq
                aux_up[i].q_ice[k] = q_ice
            end

            aux_up[i].T[k] = TD.air_temperature(thermo_params, ts_up)
            ρ = TD.air_density(thermo_params, ts_up)
            aux_up[i].buoy[k] = buoyancy_c(param_set, ρ_c[k], ρ)
            aux_up[i].RH[k] = TD.relative_humidity(thermo_params, ts_up)
        end
        aux_gm.buoy[k] = (1.0 - aux_bulk.area[k]) * aux_en.buoy[k]
        @inbounds for i in 1:N_up
            aux_gm.buoy[k] += aux_up[i].area[k] * aux_up[i].buoy[k]
        end
        @inbounds for i in 1:N_up
            aux_up[i].buoy[k] -= aux_gm.buoy[k]
        end
        aux_en.buoy[k] -= aux_gm.buoy[k]


        #####
        ##### compute bulk thermodynamics
        #####
        aux_bulk.q_liq[k] = 0
        aux_bulk.q_ice[k] = 0
        aux_bulk.T[k] = 0
        aux_bulk.RH[k] = 0
        aux_bulk.buoy[k] = 0
        if a_bulk_c > 0
            @inbounds for i in 1:N_up
                aux_bulk.q_liq[k] += aux_up[i].area[k] * aux_up[i].q_liq[k] / a_bulk_c
                aux_bulk.q_ice[k] += aux_up[i].area[k] * aux_up[i].q_ice[k] / a_bulk_c
                aux_bulk.T[k] += aux_up[i].area[k] * aux_up[i].T[k] / a_bulk_c
                aux_bulk.RH[k] += aux_up[i].area[k] * aux_up[i].RH[k] / a_bulk_c
                aux_bulk.buoy[k] += aux_up[i].area[k] * aux_up[i].buoy[k] / a_bulk_c
            end
        else
            aux_bulk.RH[k] = aux_en.RH[k]
            aux_bulk.T[k] = aux_en.T[k]
        end

        #####
        ##### update_GMV_diagnostics
        #####
        aux_gm.q_liq[k] = (aux_bulk.area[k] * aux_bulk.q_liq[k] + (1 - aux_bulk.area[k]) * aux_en.q_liq[k])
        aux_gm.q_ice[k] = (aux_bulk.area[k] * aux_bulk.q_ice[k] + (1 - aux_bulk.area[k]) * aux_en.q_ice[k])
        aux_gm.T[k] = (aux_bulk.area[k] * aux_bulk.T[k] + (1 - aux_bulk.area[k]) * aux_en.T[k])
        aux_gm.buoy[k] = (aux_bulk.area[k] * aux_bulk.buoy[k] + (1 - aux_bulk.area[k]) * aux_en.buoy[k])

        has_condensate = TD.has_condensate(aux_bulk.q_liq[k] + aux_bulk.q_ice[k])
        aux_bulk.cloud_fraction[k] = if has_condensate && a_bulk_c > 1e-3
            1
        else
            0
        end
    end
    #####
    ##### face variables: diagnose primitive, diagnose env and compute bulk
    #####
    # TODO: figure out why `ifelse` is allocating
    # clip updraft w below minimum area threshold
    @inbounds for i in 1:N_up
        a_min = edmf.minimum_area
        # a_up = aux_up[i].area
        # @. aux_up_f[i].w = ifelse(ᶠinterp_a(a_up) > a_min, max(prog_up_f[i].ρaw / (ρ_f * ᶠinterp_a(a_up)), 0), FT(0)) # should we have some cfl type limiter on this?

        @. aux_up_f[i].area = ᶠinterp_a(aux_up[i].area)  # aux_up_f[i].area is set in callbacks but condition_io() is called irregularly so maybe don't depend on it...
        a_up = aux_up_f[i].area
        @. aux_up_f[i].w = ifelse(a_up > a_min, max(prog_up_f[i].ρaw / (ρ_f * a_up), 0), FT(0)) # should we have some cfl type limiter on this?
        w_max = cfl_limit * min(grid.Δz) / Δt

        for k in real_face_indices(grid)
            if aux_up_f[i].w[k] > w_max
                z = grid.zf.z[k]
                # @warn "Updraft $i at z=$z has w = $(aux_up_f[i].w[k]) > CFL w_max $w_max, reducing" # should the reducing be done in limiters? It can't be bc limiters is called before this but this is needed to calculate tendencies...
                aux_up_f[i].w[k] = w_max
                adj = ρ_f[k] * a_up[k] * (w_max - aux_up_f[i].w[k])
                prog_up_f[i].ρaw[k] += adj
                prog_bulk_f.ρaw[k] += adj
                prog_en_gm_f.ρaw[k] -= adj
            end
        end
    end
    @inbounds for i in 1:N_up
        aux_up_f[i].w[kf_surf] = w_surface_bc(surf)
    end

    parent(aux_tc_f.bulk.w) .= 0
    @inbounds for i in 1:N_up
        a_up = aux_up[i].area
        @. aux_tc_f.bulk.w +=
            ifelse(ᶠinterp_a(aux_bulk.area) > 0, ᶠinterp_a(a_up) * aux_up_f[i].w / ᶠinterp_a(aux_bulk.area), FT(0))
    end
    # Assuming w_gm = 0!
    @. aux_en_f.w = -1 * ᶠinterp_a(aux_bulk.area) / (1 - ᶠinterp_a(aux_bulk.area)) * aux_tc_f.bulk.w

    #####
    #####  diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, grid, state, Val(:Hvar), Val(:θ_liq_ice), Val(:θ_liq_ice))
    get_GMV_CoVar(edmf, grid, state, Val(:QTvar), Val(:q_tot), Val(:q_tot))
    get_GMV_CoVar(edmf, grid, state, Val(:HQTcov), Val(:θ_liq_ice), Val(:q_tot))



    # ======================================================== # 
    # calculate and q_vap_sat_liq, q_vap_sat_ice (maybe move to assing_thermo_aux?)
    # en was already done, now do up/bulk before caling compute_nonequilibrium_moisture_tendencies!() now that aux up/bulk have ben updated.
    # use aux bc that's what compute_nonequilibrium_moisture_tendencies!() uses
    if edmf.moisture_model isa NonEquilibriumMoisture
        @inbounds for i in 1:N_up
            @. aux_up[i].q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, aux_up[i].T, ρ_c, TD.Liquid())
            @. aux_up[i].q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, aux_up[i].T, ρ_c, TD.Ice())
        end
        # This is nonlinear so really there's not a good answer here... just never use this variable I suppose...
        # @. aux_bulk.q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, aux_bulk.T, ρ_c, TD.Liquid())
        # @. aux_bulk.q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, aux_bulk.T, ρ_c, TD.Ice())
    end
    # ======================================================== #

    #####
    ##### compute_updraft_closures
    #####
    #TODO - AJ add the non-equilibrium tendency computation here
    if edmf.moisture_model isa NonEquilibriumMoisture
        compute_nonequilibrium_moisture_tendencies!(grid, state, edmf, Δt, param_set, use_fallback_tendency_limiters)
        compute_other_microphysics_tendencies!(grid, state, edmf, Δt, param_set, use_fallback_tendency_limiters) # i think only in noneq case is ok... idk... these are tendencies that would be overwritten in equilibrium case, and ql_tendency_noneq/qi_tendency_noneq don't exist in equilibrium case
    end
    compute_cloud_condensate_sedimentation_tendencies!(grid, state, edmf, Δt, param_set) # not sure on the merits of doing it here vs at the end w/ precipitation tendencies... leaving here bc originally i had it in compute_nonequilibrium_moisture_tendencies!()

    # update massflux quantities
    w_gm = prog_gm_f.w
    ∇c = CCO.DivergenceF2C()
    @inbounds for i in 1:N_up
        massflux_face_i = aux_up_f[i].massflux
        parent(massflux_face_i) .= 0
        aux_up_i = aux_up[i]
        a_up = aux_up[i].area
        w_up_i = aux_up_f[i].w
        @. aux_up_f[i].massflux =
            ifelse(ᶠinterp_a(aux_bulk.area) > 0, ρ_f * ᶠinterp_a(a_up) * (w_up_i - toscalar(w_gm)), FT(0))
        @. massflux_c += Ic(aux_up_f[i].massflux)
    end

    @inbounds for i in 1:N_up
        @. massflux += aux_up_f[i].massflux
    end
    @. ∂M∂z = ∇c(wvec(massflux))
    @. ∂M∂z_ρa = ifelse(aux_bulk.area > 0, ∂M∂z / (ρ_c * aux_bulk.area), 0)

    @inbounds for k in real_center_indices(grid)
        aux_tc.∂lnM∂z[k] = ∂M∂z[k] / (massflux_c[k] + eps(FT))
    end

    @. ∂w∂z = ∇c(wvec(aux_tc_f.bulk.w))

    compute_turb_entr!(state, grid, edmf)
    compute_phys_entr_detr!(state, grid, edmf, param_set, surf, Δt, edmf.entr_closure)
    compute_ml_entr_detr!(state, grid, edmf, param_set, surf, Δt, edmf.ml_entr_closure)
    compute_nh_pressure!(state, grid, edmf, surf)

    #####
    ##### compute_eddy_diffusivities_tke
    #####

    # Subdomain exchange term
    ∇c = CCO.DivergenceF2C()
    Ic = CCO.InterpolateF2C()
    b_exch = center_aux_turbconv(state).b_exch
    parent(b_exch) .= 0
    a_en = aux_en.area
    w_en = aux_en_f.w
    tke_en = aux_en.tke
    @inbounds for i in 1:N_up
        a_up = aux_up[i].area
        w_up = aux_up_f[i].w
        δ_dyn = aux_up[i].detr_sc
        δ_ml = aux_up[i].detr_ml
        ε_turb = aux_up[i].frac_turb_entr
        @. b_exch +=
            a_up * Ic(w_up) * (δ_dyn + δ_ml) / a_en * (1 / 2 * (Ic(w_up) - Ic(w_en))^2 - tke_en) -
            a_up * Ic(w_up) * (Ic(w_up) - Ic(w_en)) * ε_turb * Ic(w_en) / a_en
    end

    Shear² = center_aux_turbconv(state).Shear²
    ∂qt∂z = center_aux_turbconv(state).∂qt∂z
    ∂θl∂z = center_aux_turbconv(state).∂θl∂z
    ∂θv∂z = center_aux_turbconv(state).∂θv∂z
    ∂qt∂z_sat = center_aux_turbconv(state).∂qt∂z_sat
    ∂θl∂z_sat = center_aux_turbconv(state).∂θl∂z_sat
    ∂θv∂z_unsat = center_aux_turbconv(state).∂θv∂z_unsat

    ∇0_bcs = (; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    If0 = CCO.InterpolateC2F(; ∇0_bcs...)
    C123 = CCG.Covariant123Vector

    uₕ_gm = grid_mean_uₕ(state)
    w_en = aux_en_f.w
    # compute shear

    # TODO: Will need to be changed with topography
    local_geometry = CC.Fields.local_geometry_field(axes(ρ_c))
    k̂ = center_aux_turbconv(state).k̂
    @. k̂ = CCG.Contravariant3Vector(CCG.WVector(FT(1)), local_geometry)
    Ifuₕ = uₕ_bcs()
    ∇uvw = CCO.GradientF2C()
    uvw = face_aux_turbconv(state).uvw
    Shear² = center_aux_turbconv(state).Shear²
    @. uvw = C123(Ifuₕ(uₕ_gm)) + C123(wvec(w_en))
    @. Shear² = LA.norm_sqr(adjoint(∇uvw(uvw)) * k̂)

    q_tot_en = aux_en.q_tot
    θ_liq_ice_en = aux_en.θ_liq_ice
    θ_virt_en = aux_en.θ_virt
    @. ∂qt∂z = ∇c(wvec(If0(q_tot_en)))
    @. ∂θl∂z = ∇c(wvec(If0(θ_liq_ice_en)))
    @. ∂θv∂z = ∇c(wvec(If0(θ_virt_en)))

    # Second order approximation: Use dry and cloudy environmental fields.
    mix_len_params = mixing_length_params(edmf)

    @. ∂qt∂z_sat = ∇c(wvec(If0(aux_en_sat.q_tot)))
    @. ∂θl∂z_sat = ∇c(wvec(If0(aux_en_sat.θ_liq_ice)))
    @. ∂θv∂z_unsat = ∇c(wvec(If0(aux_en_unsat.θ_virt)))

    @inbounds for k in real_center_indices(grid)

        # buoyancy_gradients
        if edmf.bg_closure == BuoyGradMean()
            # First order approximation: Use environmental mean fields.
            ts_en = ts_env[k]
            bg_kwargs = (;
                t_sat = aux_en.T[k],
                qv_sat = TD.vapor_specific_humidity(thermo_params, ts_en),
                qt_sat = aux_en.q_tot[k],
                θ_sat = aux_en.θ_dry[k],
                θ_liq_ice_sat = aux_en.θ_liq_ice[k],
                ∂θv∂z_unsat = ∂θv∂z[k],
                ∂qt∂z_sat = ∂qt∂z[k],
                ∂θl∂z_sat = ∂θl∂z[k],
                p = p_c[k],
                en_cld_frac = aux_en.cloud_fraction[k],
                ρ = ρ_c[k],
            )
            bg_model = EnvBuoyGrad(edmf.bg_closure; bg_kwargs...)

        elseif edmf.bg_closure == BuoyGradQuadratures()

            bg_kwargs = (;
                t_sat = aux_en_sat.T[k],
                qv_sat = aux_en_sat.q_vap[k],
                qt_sat = aux_en_sat.q_tot[k],
                θ_sat = aux_en_sat.θ_dry[k],
                θ_liq_ice_sat = aux_en_sat.θ_liq_ice[k],
                ∂θv∂z_unsat = ∂θv∂z_unsat[k],
                ∂qt∂z_sat = ∂qt∂z_sat[k],
                ∂θl∂z_sat = ∂θl∂z_sat[k],
                p = p_c[k],
                en_cld_frac = aux_en.cloud_fraction[k],
                ρ = ρ_c[k],
            )
            bg_model = EnvBuoyGrad(edmf.bg_closure; bg_kwargs...)
        else
            error("Something went wrong. The buoyancy gradient model is not specified")
        end
        bg = buoyancy_gradients(param_set, bg_model)

        # Limiting stratification scale (Deardorff, 1976)
        # compute ∇Ri and Pr
        ∇_Ri = gradient_Richardson_number(mix_len_params, bg.∂b∂z, Shear²[k], FT(eps(FT)))
        aux_tc.prandtl_nvec[k] = turbulent_Prandtl_number(mix_len_params, obukhov_length, ∇_Ri)

        ml_model = MinDisspLen{FT}(;
            z = FT(grid.zc[k].z),
            obukhov_length = obukhov_length,
            tke_surf = aux_en.tke[kc_surf],
            ustar = surf.ustar,
            Pr = aux_tc.prandtl_nvec[k],
            p = p_c[k],
            ∇b = bg,
            Shear² = Shear²[k],
            tke = aux_en.tke[k],
            b_exch = b_exch[k],
        )

        ml = mixing_length(mix_len_params, param_set, ml_model)
        aux_tc.mls[k] = ml.min_len_ind
        aux_tc.mixing_length[k] = ml.mixing_length
        aux_tc.ml_ratio[k] = ml.ml_ratio

        KM[k] = c_m * ml.mixing_length * sqrt(max(aux_en.tke[k], 0))
        KH[k] = KM[k] / aux_tc.prandtl_nvec[k]
        KQ[k] = KH[k] / Le

        aux_en_2m.tke.buoy[k] = -aux_en.area[k] * ρ_c[k] * KH[k] * bg.∂b∂z
    end

    #####
    ##### compute covariances tendencies
    #####
    tke_press = aux_en_2m.tke.press
    w_en = aux_en_f.w
    parent(tke_press) .= 0
    @inbounds for i in 1:N_up
        w_up = aux_up_f[i].w
        nh_press = aux_up_f[i].nh_pressure
        @. tke_press += (Ic(w_en) - Ic(w_up)) * Ic(nh_press)
    end

    compute_covariance_entr(edmf, grid, state, Val(:tke), Val(:w), Val(:w))
    compute_covariance_entr(edmf, grid, state, Val(:Hvar), Val(:θ_liq_ice), Val(:θ_liq_ice))
    compute_covariance_entr(edmf, grid, state, Val(:QTvar), Val(:q_tot), Val(:q_tot))
    compute_covariance_entr(edmf, grid, state, Val(:HQTcov), Val(:θ_liq_ice), Val(:q_tot))
    compute_covariance_shear(edmf, grid, state, Val(:tke), Val(:w), Val(:w))
    compute_covariance_shear(edmf, grid, state, Val(:Hvar), Val(:θ_liq_ice), Val(:θ_liq_ice))
    compute_covariance_shear(edmf, grid, state, Val(:QTvar), Val(:q_tot), Val(:q_tot))
    compute_covariance_shear(edmf, grid, state, Val(:HQTcov), Val(:θ_liq_ice), Val(:q_tot))
    compute_covariance_dissipation(edmf, grid, state, Val(:tke), param_set)
    compute_covariance_dissipation(edmf, grid, state, Val(:Hvar), param_set)
    compute_covariance_dissipation(edmf, grid, state, Val(:QTvar), param_set)
    compute_covariance_dissipation(edmf, grid, state, Val(:HQTcov), param_set)

    # TODO defined again in compute_covariance_shear and compute_covaraince
    @inbounds for k in real_center_indices(grid)
        aux_en_2m.tke.rain_src[k] = 0
        aux_en_2m.Hvar.rain_src[k] = ρ_c[k] * aux_en.area[k] * 2 * aux_en.Hvar_rain_dt[k]
        aux_en_2m.QTvar.rain_src[k] = ρ_c[k] * aux_en.area[k] * 2 * aux_en.QTvar_rain_dt[k]
        aux_en_2m.HQTcov.rain_src[k] = ρ_c[k] * aux_en.area[k] * aux_en.HQTcov_rain_dt[k]
    end

    get_GMV_CoVar(edmf, grid, state, Val(:tke), Val(:w), Val(:w))

    compute_diffusive_fluxes(edmf, grid, state, surf, param_set)

    # TODO: use dispatch
    if edmf.precip_model isa Clima1M
        # helper to calculate the rain velocity
        # TODO: assuming w_gm = 0
        # TODO: verify translation
        term_vel_rain = aux_tc.term_vel_rain
        term_vel_snow = aux_tc.term_vel_snow
        prog_pr = center_prog_precipitation(state)

        #precip_fraction = compute_precip_fraction(edmf, state)

        @inbounds for k in real_center_indices(grid)
            # term_vel_rain[k] = CM1.terminal_velocity(microphys_params, rain_type, rain_velo_scheme, ρ_c[k], prog_pr.q_rai[k])
            # term_vel_snow[k] = CM1.terminal_velocity(microphys_params, snow_type, snow_velo_scheme, ρ_c[k], prog_pr.q_sno[k])
            # use my own that's nan-safe
            term_vel_rain[k] =
                my_terminal_velocity(microphys_params, rain_type, edmf.precip_model.rain_terminal_velocity_scheme, ρ_c[k], prog_pr.q_rai[k]) .* edmf.precip_model.rain_sedimentation_scaling_factor
            term_vel_snow[k] =
                my_terminal_velocity(microphys_params, snow_type, edmf.precip_model.snow_terminal_velocity_scheme, ρ_c[k], prog_pr.q_sno[k]) .* edmf.precip_model.snow_sedimentation_scaling_factor
        end
    end

    if edmf.cloud_sedimentation_model isa CloudSedimentationModel
        ρ_i = CMP.ρ_cloud_ice(microphys_params)
        F2Cw::CCO.InterpolateF2C = CCO.InterpolateF2C(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
        if edmf.cloud_sedimentation_model.grid_mean
            error("Not impelemented yet")
            # TODO: Impelement this. term_vel_ice is stored in aux_gm here I believe.
        else # TODO : Implement the same below for N_l
            w::CC.Fields.Field = F2Cw.(aux_en_f.w)
            term_vel_ice = aux_en.term_vel_ice
            N_i = aux_en.N_i
            @inbounds for k in real_center_indices(grid)
                N_i[k] = get_N_i(param_set, edmf.cloud_sedimentation_model.sedimentation_ice_number_concentration, ts_env[k], w[k])

                term_vel_ice[k] = calculate_sedimentation_velocity(
                    microphys_params,
                    q_effective_nan_N_safe(aux_en.q_ice[k], N_i[k], ρ_i, param_set.user_params.particle_min_radius), # if N is NaN this just returns q, which should be fine? the Chen2022Type makes some assumption about what N is... which is probably trustworthy? assume no NaNs in the vector...
                    ρ_c[k], # air density
                    ice_type,
                    N_i[k]; # broadcasting is using _first and not going through properly for some reason...
                    velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme,
                    Dmax = edmf.cloud_sedimentation_model.ice_Dmax,
                ) .* edmf.cloud_sedimentation_model.ice_sedimentation_scaling_factor 
            end

            @. aux_bulk.N_i = FT(0) # reset N_i
            @. aux_bulk.term_vel_ice = FT(0) # reset term_vel_ice for bulk so we can sum
            @inbounds for i in 1:N_up
                w = F2Cw.(aux_up_f[i].w)
                term_vel_ice = aux_up[i].term_vel_ice
                N_i = aux_up[i].N_i
                @inbounds for k in real_center_indices(grid)
                    T_up = aux_up[i].T[k]
                    q_up = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k])
                    ts_up = TD.PhaseNonEquil_ρTq(thermo_params, ρ_c[k], T_up, q_up)
                    N_i[k] = get_N_i(param_set, edmf.cloud_sedimentation_model.sedimentation_ice_number_concentration, ts_up, w[k])
                    
                    term_vel_ice[k] = calculate_sedimentation_velocity(
                        microphys_params,
                        q_effective_nan_N_safe(aux_up[i].q_ice[k], N_i[k], ρ_i, param_set.user_params.particle_min_radius), # if N is NaN this just returns q, which should be fine? the Chen2022Type makes some assumption about what N is... which is probably trustworthy? assume no NaNs in the vector...
                        ρ_c[k], # air density
                        ice_type,
                        N_i[k]; # broadcasting is using _first and not going through properly for some reason...
                        velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme,
                        Dmax = edmf.cloud_sedimentation_model.ice_Dmax,
                    ) .* edmf.cloud_sedimentation_model.ice_sedimentation_scaling_factor

                    
                    if aux_bulk.area[k] > FT(0)
                        if !isnan(N_i[k])
                            aux_bulk.N_i[k] += N_i[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add N weighted by updraft fraction
                        else
                            # is 0 right here? or something else?
                        end
                        aux_bulk.term_vel_ice[k] += term_vel_ice[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add sedimentation velocity weighted by updraft fraction
                    else
                        if !isnan(N_i[k])
                            aux_bulk.N_i[k] += N_i[k] * (1. / N_up)
                        else
                            # is 0 right here? or something else?
                        end
                        aux_bulk.term_vel_ice[k] += term_vel_ice[k] * (1. / N_up)
                    end

                end
            end
            
            # Save N_i in grid mean being sensitive to what could be NaNs in one but not the other
            @inbounds for k in real_center_indices(grid)
                if aux_bulk.area[k] > FT(0)
                    if aux_en.area[k] > FT(0)
                        aux_gm.N_i[k] = (aux_bulk.N_i[k] * aux_bulk.area[k]) + (aux_en.N_i[k] * aux_en.area[k])
                        aux_gm.term_vel_ice[k] = (aux_bulk.term_vel_ice[k] * aux_bulk.area[k] * aux_bulk.q_ice[k]  + aux_en.term_vel_ice[k] * aux_en.area[k] * aux_en.q_ice[k]) / aux_gm.q_ice[k] # mass weighted, area sums and cancels out... as does density
                    else
                        aux_gm.N_i[k] = aux_bulk.N_i[k]
                        aux_gm.term_vel_ice[k] = aux_bulk.term_vel_ice[k]
                    end
                else
                    aux_gm.N_i[k] = aux_en.N_i[k]
                    aux_gm.term_vel_ice[k] = aux_en.term_vel_ice[k]
                end
            end
        end
        # should we save N grid mean? Could be a useful thing to look at ...
    end

    ### Diagnostic thermodynamiccovariances
    if edmf.thermo_covariance_model isa DiagnosticThermoCovariances
        flux1 = surf.shf / TD.cp_m(thermo_params, ts_gm[kc_surf])
        flux2 = surf.ρq_tot_flux
        zLL::FT = grid.zc[kc_surf].z
        ustar = surf.ustar
        oblength = surf.obukhov_length
        prog_gm = center_prog_grid_mean(state)
        ρLL = prog_gm.ρ[kc_surf]
        update_diagnostic_covariances!(edmf, grid, state, param_set, Val(:Hvar))
        update_diagnostic_covariances!(edmf, grid, state, param_set, Val(:QTvar))
        update_diagnostic_covariances!(edmf, grid, state, param_set, Val(:HQTcov))
        @inbounds for k in real_center_indices(grid)
            aux_en.Hvar[k] = max(aux_en.Hvar[k], 0)
            aux_en.QTvar[k] = max(aux_en.QTvar[k], 0)
            aux_en.HQTcov[k] = max(aux_en.HQTcov[k], -sqrt(aux_en.Hvar[k] * aux_en.QTvar[k]))
            aux_en.HQTcov[k] = min(aux_en.HQTcov[k], sqrt(aux_en.Hvar[k] * aux_en.QTvar[k]))
        end
        ae_surf = 1 - aux_bulk.area[kc_surf]
        aux_en.Hvar[kc_surf] = ae_surf * get_surface_variance(flux1 / ρLL, flux1 / ρLL, ustar, zLL, oblength)
        aux_en.QTvar[kc_surf] = ae_surf * get_surface_variance(flux2 / ρLL, flux2 / ρLL, ustar, zLL, oblength)
        aux_en.HQTcov[kc_surf] = ae_surf * get_surface_variance(flux1 / ρLL, flux2 / ρLL, ustar, zLL, oblength)
    end

    compute_precipitation_formation_tendencies(
        grid,
        state,
        edmf,
        edmf.moisture_model,
        edmf.precip_model,
        edmf.cloud_sedimentation_model,
        edmf.rain_formation_model,
        edmf.snow_formation_model,
        Δt,
        param_set,
        use_fallback_tendency_limiters,
    )


    # limit env tendencies (up tendencies done after compute_up_tendencies!())
    # if !isa(edmf.tendency_limiters.up_en_tendency_limiter, NoTendencyLimiter)
        # limit_env_tendencies!(grid, state, edmf, Δt, param_set)
    # end


    return nothing
end
