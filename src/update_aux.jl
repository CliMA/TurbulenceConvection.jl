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

    # Bringing earlier for dTdz calc
    ∇0_bcs = (; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    If0 = CCO.InterpolateC2F(; ∇0_bcs...)
    RB = CCO.RightBiasedC2F(; top = CCO.Extrapolate()) # top = CCO.SetValue(T_toa)) # right biased
    ∇c = CCO.DivergenceF2C()


    w_gm = prog_gm_f.w
    w_en = aux_en_f.w



    prog_gm_uₕ = grid_mean_uₕ(state)
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
            elseif edmf.moisture_model isa EquilibriumMoisture
                thermo_args = ()
            else
                error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
            end
            
            aux_up[i].ts[k] = thermo_state_pθq(param_set, p_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k], thermo_args...)
            aux_up[i].T[k] = TD.air_temperature(thermo_params, aux_up[i].ts[k])
            aux_up[i].p[k] = TD.air_pressure(thermo_params, aux_up[i].ts[k])

            if edmf.moisture_model isa EquilibriumMoisture
                aux_up[i].q_liq[k] = TD.liquid_specific_humidity(thermo_params, aux_up[i].ts[k])
                aux_up[i].q_ice[k] = TD.ice_specific_humidity(thermo_params, aux_up[i].ts[k])
            end
        end
    end
    @inbounds for i in 1:N_up
        @. aux_up[i].dTdz = ∇c(wvec(If0(aux_up[i].T))) # compute dTdz for each updraft [[ we dont use it for boost though ]]
        # @. aux_up[i].dTdz = ∇c(wvec(RB(aux_up[i].T))) # right biased version, maybe better for downwards fluxes for boosting
    end

    # reweight in equilibrium... (env call is in microphysics())
    reweight_processes_for_grid::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_processes_for_grid, false)
    reweight_extrema_only::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false)
    if reweight_processes_for_grid && (edmf.moisture_model isa EquilibriumMoisture) # really this needs qt to be updated at all k so you arent overreaching right?
        @inbounds for i in 1:N_up
            @inbounds for k in real_center_indices(grid)
                q = TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k]) # already did the calc, let's reuse
                q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_up[i], q, p_c; reweight_extrema_only = reweight_extrema_only)
                aux_up[i].q_liq[k] = q_liq
                aux_up[i].q_ice[k] = q_ice
            end
        end
    end

    # compute bulk and decompose environment
    @inbounds for k in real_center_indices(grid)
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



        # aux_en.θ_liq_ice[k] = val1 * aux_gm.θ_liq_ice[k] - val2 * aux_bulk.θ_liq_ice[k]
        aux_en.θ_liq_ice[k] = max(val1 * aux_gm.θ_liq_ice[k] - val2 * aux_bulk.θ_liq_ice[k], 0) #Jordan - this is here to prevent negative θ_liq_ice -- why would this happen? some tendency mismatch...

        # if aux_en.θ_liq_ice[k] < 0
        #     error("aux_en.θ_liq_ice[$k] = $(aux_en.θ_liq_ice[k]); aux_gm.θ_liq_ice[$k] = $(aux_gm.θ_liq_ice[k]); aux_bulk.θ_liq_ice[$k] = $(aux_bulk.θ_liq_ice[k]); val1 = $val1; val2 = $val2")
        # end

        if edmf.moisture_model isa NonEquilibriumMoisture # otherwise it's diagnosed in thermo_state_pθq() calculation...
            # aux_en.q_liq[k] = max(val1 * prog_gm.q_liq[k] - val2 * aux_bulk.q_liq[k], 0) # in equilibrium, these get diagnosed in microphysics()
            # aux_en.q_ice[k] = max(val1 * prog_gm.q_ice[k] - val2 * aux_bulk.q_ice[k], 0) # in equilibrium, these get diagnosed in microphysics()
            aux_en.q_liq[k] = clamp(val1 * aux_gm.q_liq[k] - val2 * aux_bulk.q_liq[k], FT(0), FT(0.5) - eps(FT))
            aux_en.q_ice[k] = clamp(val1 * aux_gm.q_ice[k] - val2 * aux_bulk.q_ice[k], FT(0), FT(0.5) - eps(FT))

            aux_en.q_tot[k] = max(aux_en.q_tot[k], aux_en.q_liq[k] + aux_en.q_ice[k] + eps(FT)) # ensure that total specific humidity is at least the sum of liquid and ice

            if isnan(aux_en.q_tot[k])
                error("aux_en.q_tot[$k] is NaN; aux_gm.q_tot[$k] = $(aux_gm.q_tot[k]); aux_bulk.q_tot[$k] = $(aux_bulk.q_tot[k]); val1 = $val1; val2 = $val2")
            end

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
        aux_en.p[k] = TD.air_pressure(thermo_params, ts_en)


        ρ = TD.air_density(thermo_params, ts_en)
        aux_en.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts_en)
        aux_en.θ_dry[k] = TD.dry_pottemp(thermo_params, ts_en)
        aux_en.buoy[k] = buoyancy_c(param_set, ρ_c[k], ρ)
        aux_en.RH[k] = TD.relative_humidity(thermo_params, ts_en)
        aux_en.RH_liq[k] = relative_humidity_over_liquid(thermo_params, ts_en)
        aux_en.RH_ice[k] = relative_humidity_over_ice(thermo_params, ts_en)
    end

    @. aux_en.dTdz = ∇c(wvec(If0(aux_en.T))) # compute dTdz for env
    # @. aux_en.dTdz = ∇c(wvec(RB(aux_en.T))) # right biased version, maybe better for downwards fluxes for boosting





    #####
    #####  update prognostic bulk (these mirror the prognostic variables exactly, no minimum_area or other clippings -- filtering already happened so should be good)
    ### - these are important because while sometimes tendencies are calculated using aux, they are always applied to prog, so for limiting we need to know the exact prognostic values
    #####


    # ==================================================================================== #

    @inbounds for k in real_center_indices(grid)
        a_bulk_c = aux_bulk.area[k]
        @inbounds for i in 1:N_up
            if aux_up[i].area[k] < edmf.minimum_area && k > kc_surf && aux_up[i].area[k - 1] > 0.0 # buoyancy hack for left biased buoyancy calc... should ensure paw, buoyancy never goes below 0.
                qt = aux_up[i].q_tot[k - 1]
                h = aux_up[i].θ_liq_ice[k - 1]
                if edmf.moisture_model isa EquilibriumMoisture
                    ts_up = thermo_state_pθq(param_set, p_c[k], h, qt) # see https://github.com/CliMA/TurbulenceConvection.jl/pull/1056
                elseif edmf.moisture_model isa NonEquilibriumMoisture
                    ql = aux_up[i].q_liq[k - 1]
                    qi = aux_up[i].q_ice[k - 1]
                    ts_up = thermo_state_pθq(param_set, p_c[k], h, qt, ql, qi) # see https://github.com/CliMA/TurbulenceConvection.jl/pull/1056
                else
                    error("Something went wrong. emdf.moisture_model options are equilibrium or nonequilibrium")
                end
                # aux_up[i].q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts_up) # holdover from buoyancy hack days... but it was only done for liq,ice not qt... not sure why.. so i'm turning it off
                # aux_up[i].q_ice[k] = TD.ice_specific_humidity(thermo_params, ts_up) # holdover from buoyancy hack days... but it was only done for liq,ice not qt... not sure why.. so i'm turning it off
            else
                ts_up = aux_up[i].ts[k] # this is already set above, so we can just use it
            end

            ρ = TD.air_density(thermo_params, ts_up)
            aux_up[i].buoy[k] = buoyancy_c(param_set, ρ_c[k], ρ)
            aux_up[i].RH[k] = TD.relative_humidity(thermo_params, ts_up)
            aux_up[i].RH_liq[k] = relative_humidity_over_liquid(thermo_params, ts_up)
            aux_up[i].RH_ice[k] = relative_humidity_over_ice(thermo_params, ts_up)
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
        aux_bulk.RH_liq[k] = 0
        aux_bulk.RH_ice[k] = 0
        aux_bulk.buoy[k] = 0
        if a_bulk_c > 0
            @inbounds for i in 1:N_up
                aux_bulk.q_liq[k] += aux_up[i].area[k] * aux_up[i].q_liq[k] / a_bulk_c
                aux_bulk.q_ice[k] += aux_up[i].area[k] * aux_up[i].q_ice[k] / a_bulk_c
                aux_bulk.T[k] += aux_up[i].area[k] * aux_up[i].T[k] / a_bulk_c
                aux_bulk.RH[k] += aux_up[i].area[k] * aux_up[i].RH[k] / a_bulk_c
                aux_bulk.RH_liq[k] += aux_up[i].area[k] * aux_up[i].RH_liq[k] / a_bulk_c
                aux_bulk.RH_ice[k] += aux_up[i].area[k] * aux_up[i].RH_ice[k] / a_bulk_c
                aux_bulk.buoy[k] += aux_up[i].area[k] * aux_up[i].buoy[k] / a_bulk_c
            end
        else
            aux_bulk.RH[k] = aux_en.RH[k]
            aux_bulk.RH_liq[k] = aux_en.RH_liq[k]
            aux_bulk.RH_ice[k] = aux_en.RH_ice[k]
            aux_bulk.T[k] = aux_en.T[k]
        end

        #####
        ##### update_GMV_diagnostics (can't call in Eq until after microphysics())
        #####
        if edmf.moisture_model isa NonEquilibriumMoisture
            @inbounds for k in real_center_indices(grid)
                a_bulk_c = aux_bulk.area[k]
                #####
                ##### update_GMV_diagnostics (can't call in Eq until after microphysics())
                #####
                aux_gm.q_liq[k] = (aux_bulk.area[k] * aux_bulk.q_liq[k] + (1 - aux_bulk.area[k]) * aux_en.q_liq[k])
                aux_gm.q_ice[k] = (aux_bulk.area[k] * aux_bulk.q_ice[k] + (1 - aux_bulk.area[k]) * aux_en.q_ice[k])
                aux_gm.T[k] = (aux_bulk.area[k] * aux_bulk.T[k] + (1 - aux_bulk.area[k]) * aux_en.T[k])
                aux_gm.buoy[k] = (aux_bulk.area[k] * aux_bulk.buoy[k] + (1 - aux_bulk.area[k]) * aux_en.buoy[k])

                if iszero(aux_gm.T[k]) # print aux_gm and aux_bulk and aux_en values for q_tot, q_liq, q_ice, T, buoy and area in the error message...
                    error("aux_gm.T[$k] is 0. Other values are aux_gm.q_tot[k] = $(aux_gm.q_tot[k]), aux_gm.q_liq[k] = $(aux_gm.q_liq[k]); aux_gm.q_ice[k] = $(aux_gm.q_ice[k]); aux_gm.buoy[k] = $(aux_gm.buoy[k]); aux_bulk.area[k] = $(aux_bulk.area[k]); aux_bulk.T[k] = $(aux_bulk.T[k]); aux_en.area[k] = $(aux_en.area[k]); aux_en.T[k] = $(aux_en.T[k]); aux_en.q_tot[k] = $(aux_en.q_tot[k]); aux_en.q_liq[k] = $(aux_en.q_liq[k]); aux_en.q_ice[k] = $(aux_en.q_ice[k])")
                end

                has_condensate = TD.has_condensate(aux_bulk.q_liq[k] + aux_bulk.q_ice[k])
                aux_bulk.cloud_fraction[k] = if has_condensate && a_bulk_c > 1e-3
                    1
                else
                    0
                end
            end
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

        @inbounds for k in real_face_indices(grid)
            if aux_up_f[i].w[k] > w_max
                z = grid.zf.z[k]
                # @debug "Updraft $i at z=$z has w = $(aux_up_f[i].w[k]) > CFL w_max $w_max, reducing" # should the reducing be done in limiters? It can't be bc limiters is called before this but this is needed to calculate tendencies...
                adj = ρ_f[k] * a_up[k] * (w_max - aux_up_f[i].w[k])
                aux_up_f[i].w[k] = w_max
                prog_up_f[i].ρaw[k] += adj
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
    # @. aux_en.w = Ic(aux_en_f.w) # face to center interp [[ gonna deprecate this for space savings ]]

    #####
    #####  diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, grid, state, Val(:Hvar), Val(:θ_liq_ice), Val(:θ_liq_ice))
    get_GMV_CoVar(edmf, grid, state, Val(:QTvar), Val(:q_tot), Val(:q_tot))
    get_GMV_CoVar(edmf, grid, state, Val(:HQTcov), Val(:θ_liq_ice), Val(:q_tot))



    # === CLOAK ========================================================================================================================= #
    if edmf.area_partition_model isa CoreCloakAreaPartitionModel
        f_cm = edmf.area_partition_model.cloak_mix_factor
        f_aamaaic = edmf.area_partition_model.fraction_of_area_above_max_area_allowed_in_cloak

        a_cloak_up = aux_en.a_cloak_up
        a_cloak_dn = aux_en.a_cloak_dn
        a_en_remaining = aux_en.a_en_remaining
        #
        w_cloak_up  = aux_en_f.w_cloak_up
        w_cloak_dn  = aux_en_f.w_cloak_dn
        #
        q_tot_cloak_up = aux_en.q_tot_cloak_up
        q_tot_cloak_dn = aux_en.q_tot_cloak_dn
        #
        h_tot_cloak_up = aux_en.θ_liq_ice_cloak_up
        h_tot_cloak_dn = aux_en.θ_liq_ice_cloak_dn
        #
        if edmf.moisture_model isa NonEquilibriumMoisture # For EquilibriumMoisture, we can't update these immediately. We'll set them in the call to microphysics() based on q_tot, h_tot, etc... The grid-mean ql, qi aren't updated till later after microphysics() is called (see defereed update_GMV_diagnostics later) so we couldn't use them here even if we wanted to
            q_liq_cloak_up = aux_en.q_liq_cloak_up
            q_liq_cloak_dn = aux_en.q_liq_cloak_dn
            #
            q_ice_cloak_up = aux_en.q_ice_cloak_up
            q_ice_cloak_dn = aux_en.q_ice_cloak_dn
        end

        ## Calculate areas
        @inbounds for k in real_center_indices(grid)
            max_combined_up_area = max(aux_bulk.area[k], edmf.max_area + f_aamaaic * (one(FT) - edmf.max_area)) # allow some area above max area to go into cloak
            a_cloak_up[k] = min(aux_bulk.area[k] * edmf.area_partition_model.cloak_area_factor, max(max_combined_up_area - aux_bulk.area[k], 0)) # limit updraft cloak area to ensure combined_up_area <= max_combined_up_area
            a_cloak_dn[k] = min(aux_bulk.area[k] + a_cloak_up[k], one(FT) - (aux_bulk.area[k] + a_cloak_up[k])) # limit downdraft cloak area to ensure total area <= 1
            combined_area = aux_bulk.area[k] + a_cloak_up[k] + a_cloak_dn[k]
            a_en_remaining[k] =  max(one(FT) - combined_area, zero(FT)) # remaining env area not in cloaks
        end

        # calculate ws
        w_up = aux_tc_f.bulk.w # should be the same as aux_bulk_f.w [[ this isn't set till line 600 rn...]]
        @. w_cloak_up = f_cm * w_up # + (1 - f_cm) * w_gm # use w_gm to ensure positivity, w_gm should be 0...

        if edmf.area_partition_model.confine_all_downdraft_to_cloak
            @. w_cloak_dn = ifelse( ᶠinterp_a(a_cloak_dn) > edmf.minimum_area, (toscalar(w_gm) - (ᶠinterp_a(a_cloak_up) * w_cloak_up) - (ᶠinterp_a(aux_bulk.area) * w_up)) / ᶠinterp_a(a_cloak_dn), w_en) # for w_dn, leave nothing for env (w_en -> 0 ). w_en can be (is) negative, so do not clamp. If the cloak area is too small, just set to w_en env value, 
        else # w_en should remain unchanged, so we need to included a_en_remaining
            @. w_cloak_dn = ifelse( ᶠinterp_a(a_cloak_dn) > edmf.minimum_area, (toscalar(w_gm) - (ᶠinterp_a(a_cloak_up) * w_cloak_up) - (ᶠinterp_a(aux_bulk.area) * w_up) - (ᶠinterp_a(a_en_remaining) * w_en)) / ᶠinterp_a(a_cloak_dn), w_en) # if the cloak area is too small, just set to env value
        end


        # Calculate tracers
        @inbounds for k in real_center_indices(grid) # for loop to remove allocations for max/min checks (doesn't work for w because those use face interp)
            q_tot_cloak_max = (a_cloak_up[k] > edmf.minimum_area) ? max((aux_gm.q_tot[k] - (aux_bulk.area[k] * aux_bulk.q_tot[k]) - (a_en_remaining[k] * aux_en.q_tot[k])) / a_cloak_up[k], FT(0)) : aux_en.q_tot[k]
            h_tot_cloak_max = (a_cloak_up[k] > edmf.minimum_area) ? max((aux_gm.θ_liq_ice[k] - (aux_bulk.area[k] * aux_bulk.θ_liq_ice[k]) - (a_en_remaining[k] * aux_en.θ_liq_ice[k])) / a_cloak_up[k], FT(0)) : aux_en.θ_liq_ice[k] # limit to ensure h_cloak_dn >= 0

            q_tot_cloak_up[k] = clamp(f_cm * aux_bulk.q_tot[k] + (1 - f_cm) * aux_en.q_tot[k], 0, q_tot_cloak_max)
            h_tot_cloak_up[k] = clamp(f_cm * aux_bulk.θ_liq_ice[k] + (1 - f_cm) * aux_en.θ_liq_ice[k], 0, h_tot_cloak_max)

            # For the downdraft, we just have to close the budget. we need qi_mean to remain unchanged, and w_mean to remain unchanged
            q_tot_cloak_dn[k] = (a_cloak_dn[k] > edmf.minimum_area) ? max((aux_gm.q_tot[k] - (a_cloak_up[k] * q_tot_cloak_up[k]) - (aux_bulk.area[k] * aux_bulk.q_tot[k]) - (a_en_remaining[k] * aux_en.q_tot[k])) / a_cloak_dn[k], FT(0)) : aux_en.q_tot[k] # if the cloak area is too small, just set to env
            h_tot_cloak_dn[k] = (a_cloak_dn[k] > edmf.minimum_area) ? max((aux_gm.θ_liq_ice[k] - (a_cloak_up[k] * h_tot_cloak_up[k]) - (aux_bulk.area[k] * aux_bulk.θ_liq_ice[k]) - (a_en_remaining[k] * aux_en.θ_liq_ice[k])) / a_cloak_dn[k], FT(0)) : aux_en.θ_liq_ice[k]

            if edmf.moisture_model isa NonEquilibriumMoisture
                q_liq_up = aux_bulk.q_liq
                q_ice_up = aux_bulk.q_ice
                q_liq_cloak_max = (a_cloak_up[k] > edmf.minimum_area) ? max((aux_gm.q_liq[k] - (aux_bulk.area[k] * q_liq_up[k]) - (a_en_remaining[k] * aux_en.q_liq[k])) / a_cloak_up[k], FT(0)) : aux_en.q_liq[k]
                q_ice_cloak_max = (a_cloak_up[k] > edmf.minimum_area) ? max((aux_gm.q_ice[k] - (aux_bulk.area[k] * q_ice_up[k]) - (a_en_remaining[k] * aux_en.q_ice[k])) / a_cloak_up[k], FT(0)) : aux_en.q_ice[k] # limit to ensure q_cloak_dn >= 0

                q_liq_cloak_up[k] = clamp(f_cm * q_liq_up[k] + (1 - f_cm) * aux_en.q_liq[k], FT(0), q_liq_cloak_max)
                q_ice_cloak_up[k] = clamp(f_cm * q_ice_up[k] + (1 - f_cm) * aux_en.q_ice[k], FT(0), q_ice_cloak_max)
                # For the downdraft, we just have to close the budget. we need qi_mean to remain unchanged, and w_mean to remain unchanged
                q_liq_cloak_dn[k] = (a_cloak_dn[k] > edmf.minimum_area) ? max((aux_gm.q_liq[k] - (a_cloak_up[k] * q_liq_cloak_up[k]) - (aux_bulk.area[k] * q_liq_up[k]) - (a_en_remaining[k] * aux_en.q_liq[k])) / a_cloak_dn[k], FT(0)) : aux_en.q_liq[k] # if the cloak area is too small, just set to env
                q_ice_cloak_dn[k] = (a_cloak_dn[k] > edmf.minimum_area) ? max((aux_gm.q_ice[k] - (a_cloak_up[k] * q_ice_cloak_up[k]) - (aux_bulk.area[k] * q_ice_up[k]) - (a_en_remaining[k] * aux_en.q_ice[k])) / a_cloak_dn[k], FT(0)) : aux_en.q_ice[k]
            end


            # Derived variables

            # Thermo State
            aux_en.ts_cloak_up[k] =  (edmf.moisture_model isa EquilibriumMoisture) ? thermo_state_pθq(param_set, p_c[k], h_tot_cloak_up[k], q_tot_cloak_up[k]) : thermo_state_pθq(param_set, p_c[k], h_tot_cloak_up[k], q_tot_cloak_up[k], q_liq_cloak_up[k], q_ice_cloak_up[k])
            aux_en.ts_cloak_dn[k] =  (edmf.moisture_model isa EquilibriumMoisture) ? thermo_state_pθq(param_set, p_c[k], h_tot_cloak_dn[k], q_tot_cloak_dn[k]) : thermo_state_pθq(param_set, p_c[k], h_tot_cloak_dn[k], q_tot_cloak_dn[k], q_liq_cloak_dn[k], q_ice_cloak_dn[k])
            # Temperature
            aux_en.T_cloak_up[k] = TD.air_temperature(thermo_params, aux_en.ts_cloak_up[k])
            aux_en.T_cloak_dn[k] = TD.air_temperature(thermo_params, aux_en.ts_cloak_dn[k])
        end


    end
    # ===================================================================================================================================== #




    ### == ZERO OUT MIROPHYSICS TENDENCIES ========================================================================================== #
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
    aux_en.ql_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    aux_en.qi_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    aux_en.qt_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    aux_en.θ_liq_ice_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    #
    aux_en.ql_tendency_acnv .= FT(0)
    aux_en.qi_tendency_acnv .= FT(0)
    aux_en.qi_tendency_acnv_dep .= FT(0)
    aux_en.qi_tendency_acnv_dep_is .= FT(0)
    aux_en.qi_tendency_acnv_dep_above .= FT(0)
    aux_en.qi_tendency_acnv_agg .= FT(0)
    aux_en.qi_tendency_acnv_agg_other .= FT(0) # this is the aggregation contribution to precip formation from the other updrafts
    aux_en.qi_tendency_acnv_agg_mix .= FT(0) # this is the aggregation tendency that goes into snow (i.e. not the one that goes into rain)
    aux_en.qi_tendency_acnv_thresh .= FT(0) # this is the aggregation tendency that goes into snow (i.e. not the one that goes into rain)
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
    aux_en.qi_tendency_melt .= FT(0)
    #
    # aux_en.dqvdt .= FT(0) # zero out [ don't do this bc we need it for later calculations... in microphysics!() for example... it'll just have to be from the last timestep...]




    @inbounds for i in 1:N_up
        # aux_up[i].ql_tendency_noneq .= FT(0) # don't zero this out  bc it seeemd to break the output writing (some order of read, calculate gm, write, zero out problem probably...)
        # aux_up[i].qi_tendency_noneq .= FT(0) # don't zero this out  bc it seeemd to break the output writing (some order of read, calculate gm, write, zero out problem probably...)
        #
        aux_up[i].ql_tendency_cond_evap .= FT(0)
        aux_up[i].qi_tendency_sub_dep .= FT(0)
        #
        aux_up[i].ql_tendency_sedimentation .= FT(0)
        aux_up[i].qi_tendency_sedimentation .= FT(0)
        aux_up[i].qt_tendency_sedimentation .= FT(0)
        aux_up[i].θ_liq_ice_tendency_sedimentation .= FT(0)
        #
        aux_up[i].ql_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
        aux_up[i].qi_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
        aux_up[i].qt_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
        aux_up[i].θ_liq_ice_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
        #
        #
        aux_up[i].ql_tendency_acnv .= FT(0)
        aux_up[i].qi_tendency_acnv .= FT(0)
        aux_up[i].qi_tendency_acnv_dep .= FT(0)
        aux_up[i].qi_tendency_acnv_dep_is .= FT(0)
        aux_up[i].qi_tendency_acnv_dep_above .= FT(0)
        aux_up[i].qi_tendency_acnv_agg .= FT(0)
        aux_up[i].qi_tendency_acnv_agg_other .= FT(0) # this is the aggregation contribution to precip formation from the other updrafts
        aux_up[i].qi_tendency_acnv_agg_mix .= FT(0) # this is the aggregation tendency that goes into snow (i.e. not the one that goes into rain)
        aux_up[i].qi_tendency_acnv_thresh .= FT(0) # this is the aggregation tendency that goes into snow (i.e. not the one that goes into rain)
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
        aux_up[i].qi_tendency_melt .= FT(0)
        #
        # aux_up[i].dqvdt .= FT(0) # zero out [ don't do this bc we need it for later calculations... in microphysics!() for example... it'll just have to be from the last timestep...]
    end
    #
    aux_bulk.ql_tendency_cond_evap .= FT(0)
    aux_bulk.qi_tendency_sub_dep .= FT(0)
    #
    aux_bulk.ql_tendency_sedimentation .= FT(0)
    aux_bulk.qi_tendency_sedimentation .= FT(0)
    aux_bulk.qt_tendency_sedimentation .= FT(0)
    aux_bulk.θ_liq_ice_tendency_sedimentation .= FT(0)
    #
    aux_bulk.ql_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    aux_bulk.qi_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    aux_bulk.qt_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other updrafts
    aux_bulk.θ_liq_ice_tendency_sedimentation_other .= FT(0) # this is the sedimentation contribution to precip formation from the other up
    #
    aux_bulk.ql_tendency_acnv .= FT(0)
    aux_bulk.qi_tendency_acnv .= FT(0)
    aux_bulk.qi_tendency_acnv_dep .= FT(0)
    aux_bulk.qi_tendency_acnv_dep_is .= FT(0)
    aux_bulk.qi_tendency_acnv_dep_above .= FT(0)
    aux_bulk.qi_tendency_acnv_agg .= FT(0)
    aux_bulk.qi_tendency_acnv_agg_other .= FT(0) # this is the aggregation contribution to precip formation from the other updrafts
    aux_bulk.qi_tendency_acnv_agg_mix .= FT(0) # this is the aggregation tendency that goes into snow (i.e. not the one that goes into rain)
    aux_bulk.qi_tendency_acnv_thresh .= FT(0) # this is the aggregation tendency that goes into snow (i.e. not the one that goes into rain)
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
    aux_bulk.qi_tendency_melt .= FT(0)
    #
    #
    aux_tc.qs_tendency_accr_rai_sno .= FT(0) # this is edited in both aux_up and aux_en because the outcome is temperature dependent, so we zero it out here
    ## ================================================================================================================= ##



    # ============================================================================================================================================================================== # 
    # calculate and  q_vap_sat_liq, q_vap_sat_ice (maybe move to assing_thermo_aux?)
    if edmf.moisture_model isa NonEquilibriumMoisture
        # en has been updated, call before microphysics(). Use aux_en bc that's what microphysics uses
        # use aux_en.ts because it has the real ρ for the ts_up...
        @. aux_en.q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, aux_en.T, TD.air_density(thermo_params, aux_en.ts), TD.Liquid())
        @. aux_en.q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, aux_en.T, TD.air_density(thermo_params, aux_en.ts), TD.Ice())
    end

    # en was already done, now do up/bulk before caling compute_nonequilibrium_moisture_tendencies!() now that aux up/bulk have ben updated.
    # use aux bc that's what compute_nonequilibrium_moisture_tendencies!() uses
    if edmf.moisture_model isa NonEquilibriumMoisture
        @inbounds for i in 1:N_up
            # use aux_up[i].ts because it has the real ρ for the ts_up...
            @. aux_up[i].q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, aux_up[i].T, TD.air_density(thermo_params, aux_up[i].ts) , TD.Liquid())
            @. aux_up[i].q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, aux_up[i].T, TD.air_density(thermo_params, aux_up[i].ts) , TD.Ice())
        end
        # This is nonlinear so really there's not a good answer here... just never use this variable I suppose...
        # @. aux_bulk.q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, aux_bulk.T, ρ_c, TD.Liquid())
        # @. aux_bulk.q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, aux_bulk.T, ρ_c, TD.Ice())
    end
    # ======================================================== #

    # update massflux quantities [[ moved up for use w/ update_N_τ_termvel!() ]]
    # ∇c = CCO.DivergenceF2C()
    @inbounds for i in 1:N_up
        massflux_face_i = aux_up_f[i].massflux
        parent(massflux_face_i) .= 0
        # aux_up_i = aux_up[i]
        a_up = aux_up[i].area
        w_up_i = aux_up_f[i].w
        @. aux_up_f[i].massflux =
            ifelse(ᶠinterp_a(aux_bulk.area) > 0, ρ_f * ᶠinterp_a(a_up) * (w_up_i - toscalar(w_gm)), FT(0))
        @. massflux_c += Ic(aux_up_f[i].massflux)
    end


    # update N, τ, and termvel now that everything including w is updated.
    update_N_τ_termvel!(edmf, grid, state, param_set, thermo_params, microphys_params) # update and store N, τ, and termvel... This gets things set but also needs updraft set for efficiency for the NN, but updraft isn't set until after...
    microphysics!(edmf.en_thermo, grid, state, edmf, edmf.moisture_model, edmf.precip_model, edmf.cloud_sedimentation_model, edmf.rain_formation_model, edmf.snow_formation_model, Δt, param_set, use_fallback_tendency_limiters) # set env tendencies for microphysics

    # return back to origin state incase it needs to be zerod out for recalc
    parent(massflux) .= 0
    parent(massflux_c) .= 0
    
    #####
    ##### compute_updraft_closures
    #####
    #TODO - AJ add the non-equilibrium tendency computation here
    if edmf.moisture_model isa NonEquilibriumMoisture
        compute_nonequilibrium_moisture_tendencies!(grid, state, edmf, Δt, param_set, use_fallback_tendency_limiters)
        compute_other_microphysics_tendencies!(grid, state, edmf, Δt, param_set, use_fallback_tendency_limiters) # i think only in noneq case is ok... idk... these are tendencies that would be overwritten in equilibrium case, and ql_tendency_noneq/qi_tendency_noneq don't exist in equilibrium case
    end
    compute_cloud_condensate_sedimentation_tendencies!(grid, state, edmf, Δt, param_set) # not sure on the merits of doing it here vs at the end w/ precipitation tendencies... leaving here bc originally i had it in compute_nonequilibrium_moisture_tendencies!()



    # domain interaction microphysics
    # compute_domain_interaction_microphysics_tendencies!(grid, state, edmf, Δt, param_set, use_fallback_tendency_limiters) # this is the microphysics that happens between updrafts and environment (e.g. accretion, evaporation, etc)


    # ============================================================================================================================================================================== #

    #######
    ##### deferred grid mean calculations (done after microphysics has updated things..)
    #######
    if edmf.moisture_model isa EquilibriumMoisture
        @inbounds for k in real_center_indices(grid)
            a_bulk_c = aux_bulk.area[k]
            #####
            ##### update_GMV_diagnostics (can't call in Eq until after microphysics())
            #####
            aux_gm.q_liq[k] = (aux_bulk.area[k] * aux_bulk.q_liq[k] + (1 - aux_bulk.area[k]) * aux_en.q_liq[k])
            aux_gm.q_ice[k] = (aux_bulk.area[k] * aux_bulk.q_ice[k] + (1 - aux_bulk.area[k]) * aux_en.q_ice[k])
            aux_gm.T[k] = (aux_bulk.area[k] * aux_bulk.T[k] + (1 - aux_bulk.area[k]) * aux_en.T[k])
            aux_gm.buoy[k] = (aux_bulk.area[k] * aux_bulk.buoy[k] + (1 - aux_bulk.area[k]) * aux_en.buoy[k])

            if iszero(aux_gm.T[k]) # print aux_gm and aux_bulk and aux_en values for q_tot, q_liq, q_ice, T, buoy and area in the error message...
                error("aux_gm.T[$k] is 0. Other values are aux_gm.q_tot[k] = $(aux_gm.q_tot[k]), aux_gm.q_liq[k] = $(aux_gm.q_liq[k]); aux_gm.q_ice[k] = $(aux_gm.q_ice[k]); aux_gm.buoy[k] = $(aux_gm.buoy[k]); aux_bulk.area[k] = $(aux_bulk.area[k]); aux_bulk.T[k] = $(aux_bulk.T[k]); aux_en.area[k] = $(aux_en.area[k]); aux_en.T[k] = $(aux_en.T[k]); aux_en.q_tot[k] = $(aux_en.q_tot[k]); aux_en.q_liq[k] = $(aux_en.q_liq[k]); aux_en.q_ice[k] = $(aux_en.q_ice[k])")
            end

            has_condensate = TD.has_condensate(aux_bulk.q_liq[k] + aux_bulk.q_ice[k])
            aux_bulk.cloud_fraction[k] = if has_condensate && a_bulk_c > 1e-3
                1
            else
                0
            end
        end
    end
    

    # update massflux quantities
    # ∇c = CCO.DivergenceF2C()
    @inbounds for i in 1:N_up
        massflux_face_i = aux_up_f[i].massflux
        parent(massflux_face_i) .= 0
        # aux_up_i = aux_up[i]
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
    # ∇c = CCO.DivergenceF2C()
    # Ic = CCO.InterpolateF2C()
    b_exch = center_aux_turbconv(state).b_exch
    parent(b_exch) .= 0
    a_en = aux_en.area
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

    # ∇0_bcs = (; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # moved earlier for dTdz
    # If0 = CCO.InterpolateC2F(; ∇0_bcs...)
    C123 = CCG.Covariant123Vector

    uₕ_gm = grid_mean_uₕ(state)
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
        aux_tc.∂b∂z[k] = bg.∂b∂z # my addition

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
            tke = max(aux_en.tke[k], FT(0)),
            b_exch = b_exch[k],
        )

        ml = mixing_length(mix_len_params, param_set, ml_model)
        aux_tc.mls[k] = ml.min_len_ind
        aux_tc.mixing_length[k] = ml.mixing_length
        aux_tc.ml_ratio[k] = ml.ml_ratio

        KM[k] = c_m * ml.mixing_length * sqrt(max(aux_en.tke[k], 0))
        KH[k] = KM[k] / aux_tc.prandtl_nvec[k]
        KQ[k] = KH[k] / Le

        # aux_en_2m.tke.buoy[k] = -aux_en.area[k] * ρ_c[k] * KH[k] * bg.∂b∂z

        # 1. Calculate the normal buoyancy production [[ commenting this part out didn't make the updrafts come back...]]
        original_buoy = -aux_en.area[k] * ρ_c[k] * KH[k] * bg.∂b∂z
        # 2. Calculate a "buoyancy floor" that ramps up linearly with instability. This single line implements the ramp you described.
        min_∂b∂z = FT(1e-5)  # threshold for tke starting to be produced from buoyancy. above this value, we want the buoyancy to be unchanged.
        max_∂b∂z = FT(-1e-4) # benchmark value for very strong instability, at this ∂b∂z, we want buoyancy tke production to be at least 1e-2
        # buoy_target = FT(1e-2)  # target minimum buoyancy production rate at strong instability
        buoy_target = FT(1e-3)  # target minimum buoyancy production rate at strong instability
        buoy_floor = (buoy_target / (min_∂b∂z - max_∂b∂z )) * max(FT(0), min_∂b∂z - bg.∂b∂z)        # 3. Apply the floor to ensure a minimum production rate in strong instability
        aux_en_2m.tke.buoy[k] = (bg.∂b∂z < min_∂b∂z) ? max(original_buoy, buoy_floor) : original_buoy
    
    end

    #####
    ##### compute covariances tendencies
    #####
    tke_press = aux_en_2m.tke.press
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
                my_terminal_velocity(param_set, rain_type, edmf.precip_model.rain_terminal_velocity_scheme, ρ_c[k], prog_pr.q_rai[k])  # .* edmf.precip_model.rain_sedimentation_scaling_factor
            term_vel_snow[k] =
                my_terminal_velocity(param_set, snow_type, edmf.precip_model.snow_terminal_velocity_scheme, ρ_c[k], prog_pr.q_sno[k])  # .* edmf.precip_model.snow_sedimentation_scaling_factor
        end
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




function update_N_τ_termvel!(edmf::EDMFModel, grid::Grid, state::State, param_set::APS, thermo_params::TDPS, microphys_params::ACMP)

    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)

    FT = float_type(state)
    N_up = n_updrafts(edmf)

    ρ_c = prog_gm.ρ

    ts_env = center_aux_environment(state).ts

    massflux_c = aux_tc.massflux

    kc_toa = kc_top_of_atmos(grid)




     # ======================================================== #
    # update N_i, calculate sedimentation velocity [ TODO: Implement for N_l]
    # if edmf.moisture_model isa NonEquilibriumMoisture # we only use w to calculate N here, we aren't currently passing it into calculate_sedimentation_velocity() (we used to but stopped)
        w::CC.Fields.Field = Ic.(aux_en_f.w) # start w env since we use it first
    # end

    # ============================================================================================================================================================================================================================================================================== #

    # Stabilize our estimate of cloud top height against index jumps

    # ice  cloud top height
    # Environment
    cloud_top_ice_z_en = FT(0)
    cloud_top_ice_T_en = FT(Inf)
    k_cloud_top_ice_en = kc_surface(grid)
    N_i_cloud_top_ice_en = FT(0)
    @inbounds for k in real_center_indices(grid)
        if (TD.has_condensate(aux_en.q_ice[k]) && aux_en.area[k] > 1e-6) && ((S_k = TD.supersaturation(thermo_params,  TD.PhasePartition(thermo_params, aux_en.ts[k]), ρ_c[k], aux_en.T[k], TD.Ice())) > FT(0) )
            
            if edmf.moisture_model isa NonEquilibriumMoisture
                N_i_cloud_top_ice_en_here = get_INP_concentration(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, aux_en.ts[k]), aux_en.T[k], ρ_c[k], w[k])
            else
                N_i_cloud_top_ice_en_here = get_N_i_Cooper_curve(aux_en.T[k]; clamp_N=true) # the fallback
            end

            if N_i_cloud_top_ice_en_here ≥ N_i_cloud_top_ice_en
                N_i_cloud_top_ice_en = N_i_cloud_top_ice_en_here
                cloud_top_ice_z_en = max(cloud_top_ice_z_en, grid.zc[k].z)
                k_cloud_top_ice_en = k
                # cloud_top_ice_T_en = min(cloud_top_ice_T_en, aux_en.T[k])
                if k == kc_toa
                    T_top = aux_en.T[k]
                else
                    S_kp1 = TD.supersaturation(thermo_params, TD.PhasePartition(thermo_params, aux_en.ts[k+1]), ρ_c[k+1], aux_en.T[k+1], TD.Ice())
                    T_top = (S_kp1 < FT(0)) ? linear_interpolate_extrapolate(FT(0), (S_k, S_kp1), (aux_en.T[k], aux_en.T[k+1])) : aux_en.T[k]
                end
                cloud_top_ice_T_en = min(cloud_top_ice_T_en, T_top)
            end
        end
    end


    # updraft
    cloud_top_ice_z_ups = fill(FT(0), N_up)
    cloud_top_ice_T_ups = fill(FT(Inf), N_up)
    k_cloud_top_ice_ups = fill(kc_surface(grid), N_up)
    N_i_cloud_top_ice_ups = fill(FT(0), N_up)
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:N_up
            if (TD.has_condensate(aux_up[i].q_ice[k])) && (aux_up[i].area[k] > 1e-6) && ((S_k = TD.supersaturation(thermo_params, TD.PhasePartition(thermo_params, aux_up[i].ts[k]), ρ_c[k], aux_up[i].T[k], TD.Ice())) > FT(0))
                
                if edmf.moisture_model isa NonEquilibriumMoisture
                    # we know w_mean = 0 so we can calulate w_up 
                    w_up = -w[k] * aux_en.area[k] / (1-aux_up[i].area[k] + eps(FT))
                    N_i_cloud_top_ice_ups_here = get_INP_concentration(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, aux_up[i].ts[k]), aux_up[i].T[k], ρ_c[k], w_up)
                else
                    N_i_cloud_top_ice_ups_here = get_N_i_Cooper_curve(aux_up[i].T[k]; clamp_N=true) # the fallback
                end

                if N_i_cloud_top_ice_ups_here ≥ N_i_cloud_top_ice_ups[i]
                    N_i_cloud_top_ice_ups[i] = N_i_cloud_top_ice_ups_here
                    cloud_top_ice_z_ups[i] = max(cloud_top_ice_z_ups[i], grid.zc[k].z)
                    k_cloud_top_ice_ups[i] = k
                    # cloud_top_ice_T_ups[i] = min(cloud_top_ice_T_ups[i], aux_up[i].T[k])
                    if k == kc_toa
                        T_top = aux_up[i].T[k]
                    else
                        S_kp1 = TD.supersaturation(thermo_params, TD.PhasePartition(thermo_params, aux_up[i].ts[k+1]), ρ_c[k+1], aux_up[i].T[k+1], TD.Ice())
                        T_top = (S_kp1 < FT(0)) ? linear_interpolate_extrapolate(FT(0), (S_k, S_kp1), (aux_up[i].T[k], aux_up[i].T[k+1])) : aux_up[i].T[k]
                    end
                    cloud_top_ice_T_ups[i] = min(cloud_top_ice_T_ups[i], T_top)
                end
            end
        end
    end

    # we need to rank from highest to lowest, waterfalling any remaining area, let en be domain 0, then ups 1, 2, ..., N_up. We sort by height rather than z but idk...
    domains_zs_ks_Ts = Vector{Tuple{Int,FT,Cent{Int64},FT}}(undef, N_up+1)
    N_INP_cloud_top_ices = fill(FT(0), N_up+1) # updrafts + 1 env
    domains_zs_ks_Ts[1] = (0, cloud_top_ice_z_en, k_cloud_top_ice_en, cloud_top_ice_T_en)
    N_INP_cloud_top_ices[1] = N_i_cloud_top_ice_en
    @inbounds for i in 1:N_up
        domains_zs_ks_Ts[i+1] = (i, cloud_top_ice_z_ups[i], k_cloud_top_ice_ups[i], cloud_top_ice_T_ups[i])
        N_INP_cloud_top_ices[i+1] = N_i_cloud_top_ice_ups[i]
    end
    sort!(domains_zs_ks_Ts, by=x->x[2], rev=true) # sort by height, highest first, should be more robust than sorting by temperature (with inversions and all possibly existing)

    @inbounds for (i, (domain, _, k_top, T_top)) in enumerate(domains_zs_ks_Ts)

        if !isinf(T_top)
            if iszero(domain) # env
                N_i_cloud_top_ice_here = N_i_cloud_top_ice_en
            else # updraft
                N_i_cloud_top_ice_here = N_i_cloud_top_ice_ups[domain]
            end

            # if it's an updraft, let it just keep its internal properties unless env came before and was larger
            if !iszero(domain)
                # N_INP_cloud_top_ices[domain+1] = N_i_cloud_top_ice_here
                N_INP_cloud_top_ices[domain+1] = max(N_INP_cloud_top_ices[domain+1], N_INP_cloud_top_ices[1])
            else # if it's env, let it be a waterfalled sum of everything above it
                # area_remaining = FT(1)
                # for j in 1:(i-1) # loop over all updrafts higher than it
                #     a_up = aux_up[domains_zs_ks_Ts[j][1]].area[domains_zs_ks_Ts[j][3]]
                #     N_i_up = N_INP_cloud_top_ices[domains_zs_ks_Ts[j][1]+1]
                #     if N_i_up > N_i_cloud_top_ice_here
                #         a_count = min(area_remaining, a_up * 50) # assume anvil seeds INP everywhere by detrainment, and also that our updraft core is an undercount of the true area...
                #         N_INP_cloud_top_ices[domain+1] += (N_i_up * a_count) # allow up to full detrainment as the contribution
                #         area_remaining -= a_count
                #     end
                # end
                # N_INP_cloud_top_ices[domain+1] += N_i_cloud_top_ice_here * min(area_remaining, aux_en.area[k_top])
                # N_INP_cloud_top_ices[domain+1] = maximum(N_INP_cloud_top_ices) # max of all...
                N_INP_cloud_top_ices[domain+1] = N_INP_cloud_top_ices[domain+1] # I think giving in to the arbitrarily large updraft values allows too much boosting? idk...

            end
        end
    end

    # # print with probability 1/100
    # if rand() < 0.01
    #     @warn "N_INP_cloud_top_ices = $N_INP_cloud_top_ices"
    # end

    # # get T at cloud top z
    # N_i_cloud_top_ice = FT(NaN)
    # if !isinf(cloud_top_ice_T)
    #     if edmf.moisture_model isa NonEquilibriumMoisture
    #         # N_i_cloud_top_ice = get_N_i_raw(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, aux_en.ts[k_cloud_top_ice]), cloud_top_ice_T,  ρ_c[k_cloud_top_ice], w[k_cloud_top_ice])
    #         # N_i at cloud top doesn't matter, that could be just a smidgen above 0 for all we know, it's INP that matters as our upper bound.
    #         N_i_cloud_top_ice = get_INP_concentration(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, aux_en.ts[k_cloud_top_ice]), cloud_top_ice_T, ρ_c[k_cloud_top_ice], w[k_cloud_top_ice])
    #     else
    #         N_i_cloud_top_ice = get_N_i_Cooper_curve(cloud_top_ice_T; clamp_N=true) # the fallback
    #     end
    # end

    # ------------------------ #

    # Get the highest q_r for SIP ? Everywhere below the SIP zone needs the max ICNC raised...
    # ICNC_SIP_scaling_factors = similar(aux_gm.q_ice) # it needs its own vector so we can work up from the bottom, one factor won't cut it
    ICNC_SIP_scaling_factors = aux_gm.f_ice_mult
    ICNC_SIP_scaling_factors .= FT(1) # reset
    ICNC_SIP_scaling_factor = FT(1) # this is the running value, we only allow it to increase as we move towards sfc. this is because an INP explosion will have to be propagated downwards, a 10^3 growth permanently impacts everywhere below. (we could for example propose just turning off the INP upper bound or something idk.  +seeder-feeder and all that)
    
    if param_set.user_params.use_ice_mult
        @inbounds for k in Base.Iterators.reverse(real_center_indices(grid)) # have to use fully qualified Base.Iterators.reverse(), see Grid.jl for implementation. But we wanna go TOA to SFC since the factor increases w/ depth
            if (prog_pr.q_rai[k] > FT(0)) # drizzle drive both Hallet-Mossop and Droplet Shattering ICNC growth
                if edmf.moisture_model isa NonEquilibriumMoisture
                    N_INP = get_INP_concentration(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, aux_gm.ts[k]), aux_gm.T[k], ρ_c[k], w[k])
                else
                    N_INP = get_N_i_Cooper_curve(aux_gm.T[k]; clamp_N=true)
                end

                if (ice_mult_factor_candidate = get_ice_mult_factor_ICNC_max(param_set, N_INP, aux_gm.N_i[k], aux_gm.q_ice[k], prog_pr.q_rai[k], prog_pr.q_sno[k], aux_gm.T[k], ρ_c[k])) > ICNC_SIP_scaling_factor
                    ICNC_SIP_scaling_factor = ice_mult_factor_candidate
                end
                ICNC_SIP_scaling_factors[k] = ICNC_SIP_scaling_factor
            else
                ICNC_SIP_scaling_factors[k] = ICNC_SIP_scaling_factor
            end
        end
    end


    # ============================================================================================================================================================================================================================================================================== #

    # ============================================================================================================================================================================================================================================================================== #

    apply_massflux_boost = param_set.user_params.apply_massflux_N_i_boost # env only
    apply_sedimentation_boost = param_set.user_params.apply_sedimentation_N_i_boost # everywhere, default to false

    # ------------------------ #

    # env
    if edmf.moisture_model isa NonEquilibriumMoisture
        if edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale # can we do the entire vector at once for efficiency?
            # get_τs_and_Ns!(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition.(thermo_params, ts_env), aux_en.T, aux_en.p, TD.air_density.(thermo_params, ts_env), w, aux_en.area, aux_en.τ_liq, aux_en.τ_ice, aux_en.N_l, aux_en.N_i, ICNC_SIP_scaling_factors, prog_pr.q_sno, massflux_c, aux_en.dTdz, aux_en.term_vel_ice; N_INP_top = N_INP_cloud_top_ices[1], apply_massflux_boost=apply_massflux_boost, apply_sedimentation_boost=apply_sedimentation_boost)
            get_τs_and_Ns_and_N_i_no_boost!(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition.(thermo_params, ts_env), aux_en.T, aux_en.p, TD.air_density.(thermo_params, ts_env), w, aux_en.area, aux_en.τ_liq, aux_en.τ_ice, aux_en.N_l, aux_en.N_i, aux_en.N_i_no_boost, ICNC_SIP_scaling_factors, prog_pr.q_sno, massflux_c, aux_en.dTdz, aux_en.term_vel_ice; N_INP_top = N_INP_cloud_top_ices[1], apply_massflux_boost=apply_massflux_boost, apply_sedimentation_boost=apply_sedimentation_boost)
        end
    end # no NN for eq

    @inbounds for k in real_center_indices(grid)
        if aux_en.area[k] > 0
            if edmf.moisture_model isa NonEquilibriumMoisture
                # aux_en.N_i[k] = get_N_i(param_set, edmf.moisture_model.scheme, ts_env[k], w[k])
                # aux_en.N_l[k] = get_N_l(param_set, edmf.moisture_model.scheme, ts_env[k], w[k])
                # aux_en.N_l[k], aux_en.N_i[k] = get_Ns(param_set, edmf.moisture_model.scheme, ts_env[k], w[k])
                # aux_en.N_l[k], aux_en.N_i[k] = get_Ns(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_env[k]), aux_en.T[k], TD.air_density(thermo_params, ts_env[k]), w[k]) # cheaper bc T and ρ are stored in the ts.
                # aux_en.τ_liq[k], aux_en.τ_ice[k] = get_τs(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_env[k]), aux_en.T[k], aux_en.p[k], TD.air_density(thermo_params, ts_env[k]), w[k]) # cheaper bc T and ρ are stored in the ts.
                if !(edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale)
                    # aux_en.τ_liq[k], aux_en.τ_ice[k], aux_en.N_l[k], aux_en.N_i[k] = get_τs_and_Ns(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_env[k]), aux_en.T[k], aux_en.p[k], TD.air_density(thermo_params, ts_env[k]), w[k]; N_INP_top = N_INP_cloud_top_ices[1], f_ice_mult = ICNC_SIP_scaling_factors[k], q_sno = prog_pr.q_sno[k], massflux=massflux_c[k], dTdz=aux_en.dTdz[k], w_i = aux_en.term_vel_ice[k], apply_massflux_boost=apply_massflux_boost, apply_sedimentation_boost=apply_sedimentation_boost)
                    aux_en.τ_liq[k], aux_en.τ_ice[k], aux_en.N_l[k], aux_en.N_i[k], aux_en.N_i_no_boost[k] = get_τs_and_Ns_and_N_i_no_boost(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_env[k]), aux_en.T[k], aux_en.p[k], TD.air_density(thermo_params, ts_env[k]), w[k]; N_INP_top = N_INP_cloud_top_ices[1], f_ice_mult = ICNC_SIP_scaling_factors[k], q_sno = prog_pr.q_sno[k], massflux=massflux_c[k], dTdz=aux_en.dTdz[k], w_i = aux_en.term_vel_ice[k], apply_massflux_boost=apply_massflux_boost, apply_sedimentation_boost=apply_sedimentation_boost)
                end
            else
                if !(edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale)
                    # aux_en.N_l[k], aux_en.N_i[k] = get_Ns(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_env[k]), aux_en.T[k], TD.air_density(thermo_params, ts_env[k]), w[k]; N_INP_top = N_INP_cloud_top_ices[1], f_ice_mult = ICNC_SIP_scaling_factors[k], q_sno = prog_pr.q_sno[k], massflux=massflux_c[k], dTdz=aux_en.dTdz[k], w_i = aux_en.term_vel_ice[k], apply_massflux_boost=apply_massflux_boost, apply_sedimentation_boost=apply_sedimentation_boost)
                    aux_en.N_l[k], aux_en.N_i[k], aux_en.N_i_no_boost[k] = get_Ns_and_N_i_no_boost(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_env[k]), aux_en.T[k], TD.air_density(thermo_params, ts_env[k]), w[k]; N_INP_top = N_INP_cloud_top_ices[1], f_ice_mult = ICNC_SIP_scaling_factors[k], q_sno = prog_pr.q_sno[k], massflux=massflux_c[k], dTdz=aux_en.dTdz[k], w_i = aux_en.term_vel_ice[k], apply_massflux_boost=apply_massflux_boost, apply_sedimentation_boost=apply_sedimentation_boost)
                end
            end
            aux_en.r_l_mean[k] = r_from_qN(param_set, liq_type, aux_en.q_liq[k], aux_en.N_l[k]; monodisperse = true, ρ = ρ_c[k])
            aux_en.r_i_mean[k] = r_from_qN(param_set, ice_type, aux_en.q_ice[k], aux_en.N_i[k]; monodisperse = false, ρ = ρ_c[k])

            if edmf.cloud_sedimentation_model isa CloudSedimentationModel
                if edmf.cloud_sedimentation_model.grid_mean
                    error("Not impelemented yet") # TODO: Impelement this. term_vel_ice is stored in aux_gm here I believe
                else

                    ρ_l = CMP.ρ_cloud_liq(microphys_params)
                    aux_en.term_vel_liq[k] = calculate_sedimentation_velocity(
                        param_set,
                        q_effective_nan_N_safe(param_set, liq_type, aux_en.q_liq[k], aux_en.N_l[k]; monodisperse = true), # if N is NaN this just returns q, which should be fine? the Chen2022Type makes some assumption about what N is... which is probably trustworthy? assume no NaNs in the vector...
                        ρ_c[k], # air density
                        liq_type,
                        aux_en.N_l[k]; # broadcasting is using _first and not going through properly for some reason...
                        velo_scheme = edmf.cloud_sedimentation_model.liq_terminal_velocity_scheme,
                        # Dmax = edmf.cloud_sedimentation_model.liq_Dmax,
                        Dmax = FT(Inf), # testing [i think it is better]
                    ) .* edmf.cloud_sedimentation_model.liq_sedimentation_scaling_factor
                    
                    ρ_i = CMP.ρ_cloud_ice(microphys_params)
                    aux_en.term_vel_ice[k] = calculate_sedimentation_velocity(
                        param_set,
                        q_effective_nan_N_safe(param_set, ice_type, aux_en.q_ice[k], aux_en.N_i[k]; monodisperse = true), # if N is NaN this just returns q, which should be fine? the Chen2022Type makes some assumption about what N is... which is probably trustworthy? assume no NaNs in the vector...
                        ρ_c[k], # air density
                        ice_type,
                        aux_en.N_i[k]; # broadcasting is using _first and not going through properly for some reason...
                        velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme,
                        # Dmax = edmf.cloud_sedimentation_model.ice_Dmax,
                        Dmax = FT(Inf), # testing [ i think it is better ]
                    ) .* edmf.cloud_sedimentation_model.ice_sedimentation_scaling_factor
                end
            end
        else
            aux_en.N_l[k] = FT(0)
            aux_en.N_i[k] = FT(0)
            aux_en.r_l_mean[k] = r_from_qN(param_set, liq_type, FT(0), FT(0); monodisperse = true, ρ = ρ_c[k])
            aux_en.r_i_mean[k] = r_from_qN(param_set, ice_type, FT(0), FT(0); monodisperse = false, ρ = ρ_c[k])
            aux_en.term_vel_liq[k] = FT(0)
            aux_en.term_vel_ice[k] = FT(0)
            aux_en.τ_liq[k], aux_en.τ_ice[k] = FT(Inf), FT(Inf) # we can't store NaNs or Infs and have the regridder not complain..., but we don't want to waste compute calling the methods when the values wont be used...
        end
    end

    # -------------------- #

    # updraft
    @. aux_bulk.N_l = FT(0) # reset N_l
    @. aux_bulk.N_i = FT(0) # reset N_i
    @. aux_bulk.r_l_mean = FT(0) # reset r_l_mean
    @. aux_bulk.r_i_mean = FT(0) # reset r_i_mean
    @. aux_bulk.term_vel_liq = FT(0) # reset term_vel_liq
    @. aux_bulk.term_vel_ice = FT(0) # reset term_vel_ice

    if edmf.moisture_model isa NonEquilibriumMoisture
        @. aux_bulk.τ_liq = FT(0) # reset τ_liq [ store 1/τ here so that we can sum and then take the inverse at the end , so 1/Inf is 0 ]
        @. aux_bulk.τ_ice = FT(0) # reset τ_ice [ store 1/τ here so that we can sum and then take the inverse at the end , so 1/Inf is 0 ]
    end


    @inbounds for i in 1:N_up
        if edmf.moisture_model isa NonEquilibriumMoisture
            @. w = Ic(aux_up_f[i].w) # reuse same w allocation
        end
        ts_up = aux_up[i].ts # reuse ts_up allocation


        if edmf.moisture_model isa NonEquilibriumMoisture
            if edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale # can we do the entire vector at once for efficiency?
                get_τs_and_Ns!(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition.(thermo_params, ts_up), aux_up[i].T, aux_up[i].p, TD.air_density.(thermo_params, ts_up), w, aux_up[i].area, aux_up[i].τ_liq, aux_up[i].τ_ice, aux_up[i].N_l, aux_up[i].N_i, ICNC_SIP_scaling_factors, prog_pr.q_sno, massflux_c, aux_up[i].dTdz, aux_up[i].term_vel_ice; N_INP_top = N_INP_cloud_top_ices[i+1], apply_massflux_boost=false, apply_sedimentation_boost=apply_sedimentation_boost) # not sure if massflux boost should be 0 in the updraft
            end
        else
            if edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale # can we do the entire vector at once for efficiency?
                get_Ns!(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition.(thermo_params, ts_up), aux_up[i].T, aux_up[i].p, TD.air_density.(thermo_params, ts_up), w, aux_up[i].area, aux_up[i].τ_liq, aux_up[i].τ_ice, aux_up[i].N_l, aux_up[i].N_i, ICNC_SIP_scaling_factors, prog_pr.q_sno, massflux_c, aux_up[i].dTdz, aux_up[i].term_vel_ice; N_INP_top = N_INP_cloud_top_ices[i+1], apply_massflux_boost=false, apply_sedimentation_boost=apply_sedimentation_boost) # not sure if massflux boost should be 0 in the updraft
            end
        end

        @inbounds for k in real_center_indices(grid)
            if aux_up[i].area[k] > 0
                # calculate N_i
                if edmf.moisture_model isa NonEquilibriumMoisture
                    # aux_up[i].N_i[k] = get_N_i(param_set, edmf.moisture_model.scheme, ts_up, w[k])
                    # aux_up[i].N_l[k] = get_N_l(param_set, edmf.moisture_model.scheme, ts_up, w[k])
                    # aux_up[i].N_l[k], aux_up[i].N_i[k] = get_Ns(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_up[k]), aux_up[i].T[k], TD.air_density(thermo_params, ts_up[k]), w[k]) # cheaper bc T and ρ are stored in the ts.
                    # aux_up[i].τ_liq[k], aux_up[i].τ_ice[k] = get_τs(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_up[k]), aux_up[i].T[k], aux_up[i].p[k], TD.air_density(thermo_params, ts_up[k]), w[k]) # cheaper bc T and ρ are stored in the ts.
                    if !(edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale) # you could use the vector methods but probably just makes more allocations since a for loop is probably sufficient. not sure the SIMD speedup is meanigful.
                        aux_up[i].τ_liq[k], aux_up[i].τ_ice[k], aux_up[i].N_l[k], aux_up[i].N_i[k] = get_τs_and_Ns(param_set, microphys_params, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_up[k]), aux_up[i].T[k], aux_up[i].p[k], TD.air_density(thermo_params, ts_up[k]), w[k]; N_INP_top = N_INP_cloud_top_ices[i+1], f_ice_mult = ICNC_SIP_scaling_factors[k], q_sno = prog_pr.q_sno[k], massflux = massflux_c[k], dTdz = aux_up[i].dTdz[k], w_i = aux_up[i].term_vel_ice[k], apply_massflux_boost=false, apply_sedimentation_boost=apply_sedimentation_boost) # not sure if massflux boost should be 0 in the updraft
                    end
                else
                    # we need to come up with what we do for equilibrium... rn for base and others we just have no N if it's not predicted, but maybe that's not ideal idk... maybe there should be a fallback method for eq and for the relaationtimescales that don't predict N that just gets the default values based on the given q.
                    # we just need to get_Ns
                    if !(edmf.moisture_model.scheme isa NeuralNetworkRelaxationTimescale) # you could use the vector methods but probably just makes more allocations since a for loop is probably sufficient. not sure the SIMD speedup is meanigful.
                        aux_up[i].N_l[k], aux_up[i].N_i[k] = get_Ns(param_set, edmf.moisture_model.scheme, TD.PhasePartition(thermo_params, ts_up[k]), aux_up[i].T[k], TD.air_density(thermo_params, ts_up[k]), w[k]; N_INP_top = N_INP_cloud_top_ices[i+1], f_ice_mult = ICNC_SIP_scaling_factors[k], q_sno = prog_pr.q_sno[k], massflux = massflux_c[k], dTdz = aux_up[i].dTdz[k], w_i = aux_up[i].term_vel_ice[k], apply_massflux_boost=false, apply_sedimentation_boost=apply_sedimentation_boost) # not sure if massflux boost should be 0 in the updraft
                    end
                end

                aux_up[i].r_l_mean[k] = r_from_qN(param_set, liq_type, aux_up[i].q_liq[k], aux_up[i].N_l[k]; monodisperse = true, ρ = ρ_c[k])
                aux_up[i].r_i_mean[k] = r_from_qN(param_set, ice_type, aux_up[i].q_ice[k], aux_up[i].N_i[k]; monodisperse = false, ρ = ρ_c[k])

                # store N_i_no_boost
                aux_up[i].N_i_no_boost[k] = aux_up[i].N_i[k] # updrafts don't get massflux boost
                
                # calculate sedimentation velocity
                if edmf.cloud_sedimentation_model isa CloudSedimentationModel

                    # calculate terminal velocity
                    if edmf.cloud_sedimentation_model.grid_mean
                        error("Not impelemented yet")     # TODO: Impelement this. term_vel_ice is stored in aux_gm here I believe
                    else

                        # liquid terminal velocity
                        ρ_l = CMP.ρ_cloud_liq(microphys_params)
                        aux_up[i].term_vel_liq[k] = calculate_sedimentation_velocity(
                            param_set,
                            q_effective_nan_N_safe(param_set, liq_type, aux_up[i].q_liq[k], aux_up[i].N_l[k]; monodisperse = true), # if N is NaN this just returns q, which should be fine? the Chen2022Type makes some assumption about what N is... which is probably trustworthy? assume no NaNs in the vector...
                            ρ_c[k], # air density
                            liq_type,
                            aux_up[i].N_l[k]; # broadcasting is using _first and not going through properly for some reason...
                            velo_scheme = edmf.cloud_sedimentation_model.liq_terminal_velocity_scheme,
                            # Dmax = edmf.cloud_sedimentation_model.liq_Dmax,
                            Dmax = FT(Inf), # testing [ i think it is better ]
                        ) .* edmf.cloud_sedimentation_model.liq_sedimentation_scaling_factor
                            
                        # ice terminal veocity
                        ρ_i = CMP.ρ_cloud_ice(microphys_params)
                        aux_up[i].term_vel_ice[k] = calculate_sedimentation_velocity(
                            param_set,
                            q_effective_nan_N_safe(param_set, ice_type, aux_up[i].q_ice[k], aux_up[i].N_i[k]; monodisperse = true), # if N is NaN this just returns q, which should be fine? the Chen2022Type makes some assumption about what N is... which is probably trustworthy? assume no NaNs in the vector...
                            ρ_c[k], # air density
                            ice_type,
                            aux_up[i].N_i[k]; # broadcasting is using _first and not going through properly for some reason...
                            velo_scheme = edmf.cloud_sedimentation_model.ice_terminal_velocity_scheme,
                            # Dmax = edmf.cloud_sedimentation_model.ice_Dmax,
                            Dmax = FT(Inf), # testing [ i think it is better ]
                        ) .* edmf.cloud_sedimentation_model.ice_sedimentation_scaling_factor
                    end
                end

                # calculate bulk N_i and term_vel_ice and N_l and term_vel_liq
                if aux_bulk.area[k] > FT(0)

                    # liquid 
                    if !isnan(aux_up[i].N_l[k])
                        aux_bulk.N_l[k] += aux_up[i].N_l[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add N weighted by updraft fraction
                        # r_l needs to be weighted by total N, so multiply by N*area
                        aux_bulk.r_l_mean[k] += aux_up[i].r_l_mean[k] * aux_up[i].N_l[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add r_l weighted by updraft fraction
                    else
                        # is 0 right here? or something else?
                    end
                    aux_bulk.term_vel_liq[k] += aux_up[i].term_vel_liq[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add sedimentation velocity weighted by updraft fraction

                    # ice 
                    if !isnan(aux_up[i].N_i[k])
                        aux_bulk.N_i[k] += aux_up[i].N_i[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add N weighted by updraft fraction
                        aux_bulk.r_i_mean[k] += aux_up[i].r_i_mean[k] * aux_up[i].N_i[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add r_i weighted by updraft fraction

                    else
                        # is 0 right here? or something else?
                    end
                    aux_bulk.term_vel_ice[k] += aux_up[i].term_vel_ice[k] * (aux_up[i].area[k] / aux_bulk.area[k]) # add sedimentation velocity weighted by updraft fraction

                    if edmf.moisture_model isa NonEquilibriumMoisture
                        aux_bulk.τ_liq[k] += inv(aux_up[i].τ_liq[k]) * (aux_up[i].area[k] / aux_bulk.area[k]) # add τ_liq weighted by updraft fraction
                        aux_bulk.τ_ice[k] += inv(aux_up[i].τ_ice[k]) * (aux_up[i].area[k] / aux_bulk.area[k]) # add τ_ice weighted by updraft fraction
                    end

                else
                    # liquid
                    if !isnan(aux_up[i].N_l[k])
                        aux_bulk.N_l[k] += aux_up[i].N_l[k] * (1. / N_up) # add N weighted by updraft fraction
                        aux_bulk.r_l_mean[k] += aux_up[i].r_l_mean[k] * aux_up[i].N_l[k] * (1. / N_up) # add r_l weighted by updraft fraction
                    else
                        # is 0 right here? or something else?
                    end
                    aux_bulk.N_i[k] += aux_up[i].N_i[k] * (1. / N_up) # add N weighted by updraft fraction

                    # ice
                    if !isnan(aux_up[i].N_i[k])
                        aux_bulk.N_i[k] += aux_up[i].N_i[k] * (1. / N_up)
                        aux_bulk.r_i_mean[k] += aux_up[i].r_i_mean[k] * aux_up[i].N_i[k] * (1. / N_up) # add r_i weighted by updraft fraction

                    else
                        # is 0 right here? or something else?
                    end
                    aux_bulk.term_vel_ice[k] += aux_up[i].term_vel_ice[k] * (1. / N_up)

                    if edmf.moisture_model isa NonEquilibriumMoisture
                        aux_bulk.τ_liq[k] += inv(aux_up[i].τ_liq[k]) * (1. / N_up) # add τ_liq weighted by updraft fraction
                        aux_bulk.τ_ice[k] += inv(aux_up[i].τ_ice[k]) * (1. / N_up) # add τ_ice weighted by updraft fraction
                    end

                end
            else
                aux_up[i].N_l[k] = FT(0)
                aux_up[i].N_i[k] = FT(0)
                aux_up[i].r_l_mean[k] = r_from_qN(param_set, liq_type, FT(0), FT(0); monodisperse = true, ρ = ρ_c[k])
                aux_up[i].r_i_mean[k] = r_from_qN(param_set, ice_type, FT(0), FT(0); monodisperse = false, ρ = ρ_c[k])
                aux_up[i].term_vel_ice[k] = FT(0)
                aux_up[i].term_vel_liq[k] = FT(0)


                if edmf.moisture_model isa NonEquilibriumMoisture
                    aux_up[i].τ_liq[k] = FT(Inf) # reset τ_liq
                    aux_up[i].τ_ice[k] = FT(Inf) # reset τ_ice
                end

            end
            #
        end
    end
    if edmf.moisture_model isa NonEquilibriumMoisture
        @. aux_bulk.τ_liq = inv(aux_bulk.τ_liq) # undo inverse
        @. aux_bulk.τ_ice = inv(aux_bulk.τ_ice) # undo inverse
    end

    # to get r_mean we need to divide by the toal N_l, N_i
    @inbounds for k in real_center_indices(grid)
        aux_bulk.r_l_mean[k] = iszero(aux_bulk.N_l[k]) ? r_from_qN(param_set, liq_type, FT(0), FT(0); monodisperse = true, ρ = ρ_c[k]) : aux_bulk.r_l_mean[k] / aux_bulk.N_l[k] # if N_l is 0 then r_l_mean is the default minimum
        aux_bulk.r_i_mean[k] = iszero(aux_bulk.N_i[k]) ? r_from_qN(param_set, ice_type, FT(0), FT(0); monodisperse = false, ρ = ρ_c[k]) : aux_bulk.r_i_mean[k] / aux_bulk.N_i[k] # if N_i is 0 then r_i_mean is the default minimum

        aux_bulk.N_i_no_boost[k] = aux_bulk.N_i[k]
    end

    # -------------------- #

    # update grid mean N_i and term_vel_ice, being sensitive to what could be NaNs in one but not the other
    if (edmf.cloud_sedimentation_model isa CloudSedimentationModel) && !edmf.cloud_sedimentation_model.grid_mean
        @inbounds for k in real_center_indices(grid)
            if aux_bulk.area[k] > FT(0)
                if aux_en.area[k] > FT(0)
                    aux_gm.N_i[k] = (aux_bulk.N_i[k] * aux_bulk.area[k]) + (aux_en.N_i[k] * aux_en.area[k])
                    aux_gm.N_l[k] = (aux_bulk.N_l[k] * aux_bulk.area[k]) + (aux_en.N_l[k] * aux_en.area[k])
                    aux_gm.r_l_mean[k] = (aux_bulk.r_l_mean[k] * aux_bulk.N_l[k] * aux_bulk.area[k] + aux_en.r_l_mean[k] * aux_en.N_l[k] * aux_en.area[k]) / (aux_gm.N_l[k]) # area is 1
                    aux_gm.r_i_mean[k] = (aux_bulk.r_i_mean[k] * aux_bulk.N_i[k] * aux_bulk.area[k] + aux_en.r_i_mean[k] * aux_en.N_i[k] * aux_en.area[k]) / (aux_gm.N_i[k]) # area is 1

                    aux_gm.N_i_no_boost[k] = (aux_bulk.N_i_no_boost[k] * aux_bulk.area[k]) + (aux_en.N_i_no_boost[k] * aux_en.area[k]) # bulk N_i_no_boost is just N_i since bulk doesn't get massflux boost

                    # calculate grid mean term_vel_liq
                    if aux_gm.q_liq[k] > FT(0) # if q_liq is 0 then N_l is 0
                        aux_gm.term_vel_liq[k] = (aux_bulk.term_vel_liq[k] * aux_bulk.area[k] * aux_bulk.q_liq[k] + aux_en.term_vel_liq[k] * aux_en.area[k] * aux_en.q_liq[k]) / aux_gm.q_liq[k] # mass weighted, area sums and cancels out... as does density
                    else
                        aux_gm.term_vel_liq[k] = FT(0) # i think this is already true but just to be sure
                    end

                    # calculate grid mean term_vel_ice
                    if aux_gm.q_ice[k] > FT(0)
                        aux_gm.term_vel_ice[k] = (aux_bulk.term_vel_ice[k] * aux_bulk.area[k] * aux_bulk.q_ice[k] + aux_en.term_vel_ice[k] * aux_en.area[k] * aux_en.q_ice[k]) / aux_gm.q_ice[k] # mass weighted, area sums and cancels out... as does density
                    else
                        aux_gm.term_vel_ice[k] = FT(0) # i think this is already true but just to be sure
                    end

                else
                    aux_gm.N_i[k] = aux_bulk.N_i[k]
                    aux_gm.r_i_mean[k] = aux_bulk.r_i_mean[k]
                    aux_gm.term_vel_ice[k] = aux_bulk.term_vel_ice[k]
                    aux_gm.N_l[k] = aux_bulk.N_l[k]
                    aux_gm.r_l_mean[k] = aux_bulk.r_l_mean[k]
                    aux_gm.term_vel_liq[k] = aux_bulk.term_vel_liq[k]

                    aux_gm.N_i_no_boost[k] = aux_bulk.N_i_no_boost[k] # bulk N_i_no_boost is just N_i since bulk doesn't get massflux boost
                end
            else
                aux_gm.N_i[k] = aux_en.N_i[k]
                aux_gm.r_i_mean[k] = aux_en.r_i_mean[k]
                aux_gm.term_vel_ice[k] = aux_en.term_vel_ice[k]
                aux_gm.N_l[k] = aux_en.N_l[k]
                aux_gm.r_l_mean[k] = aux_en.r_l_mean[k]
                aux_gm.term_vel_liq[k] = aux_en.term_vel_liq[k]

                aux_gm.N_i_no_boost[k] = aux_en.N_i_no_boost[k]
            end
        end
    end


    # update gradients #    
    wvec = CC.Geometry.WVector
    ∇0_bcs = (; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    If0 = CCO.InterpolateC2F(; ∇0_bcs...)
    RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0))) # top = CCO.SetValue(T_toa)) # right biased
    ∇c = CCO.DivergenceF2C()
    # @. aux_en.dr_i_mean_dz = ∇c(wvec(If0(aux_en.r_i_mean)))
    # @. aux_en.dr_l_mean_dz = ∇c(wvec(If0(aux_en.r_l_mean)))
    @. aux_en.dN_i_dz = ∇c(wvec(RB(aux_en.N_i)))
    # @. aux_en.dN_l_dz = ∇c(wvec(If0(aux_en.N_l)))
    @. aux_en.dqidz = ∇c(wvec(RB(aux_en.q_ice))) # right biased divergence for sedimentation flux

    @inbounds for i in 1:N_up
    #     @. aux_up[i].dr_i_mean_dz = ∇c(wvec(If0(aux_up[i].r_i_mean)))
    #     @. aux_up[i].dr_l_mean_dz = ∇c(wvec(If0(aux_up[i].r_l_mean)))
        @. aux_up[i].dN_i_dz = ∇c(wvec(RB(aux_up[i].N_i)))
    #     @. aux_up[i].dN_l_dz = ∇c(wvec(If0(aux_up[i].N_l)))
        @. aux_up[i].dqidz = ∇c(wvec(RB(aux_up[i].q_ice))) # right biased divergence for sedimentation flux
    end



    # ======================================================== #

    return nothing
end