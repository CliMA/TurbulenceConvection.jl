function update_aux!(edmf::EDMF_PrognosticTKE, gm, grid::Grid, state::State, Case, param_set::APS, t::Real, Δt::Real)
    #####
    ##### Unpack common variables
    #####
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    up = edmf.UpdVar
    ρ0_f = face_ref_state(state).ρ0
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    α0_c = center_ref_state(state).α0
    g = CPP.grav(param_set)
    c_m = CPEDMF.c_m(param_set)
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    surface = Case.Sur
    obukhov_length = surface.obukhov_length
    FT = eltype(grid)
    prog_gm = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_tc = center_aux_turbconv(state)
    aux_bulk = center_aux_bulk(state)
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_en_unsat = aux_en.unsat
    aux_en_sat = aux_en.sat
    m_entr_detr = aux_tc.ϕ_temporary
    ∇m_entr_detr = aux_tc.ψ_temporary
    w_up_c = aux_tc.w_up_c
    w_en_c = aux_tc.w_en_c
    wvec = CC.Geometry.WVector
    max_area = edmf.max_area
    sde_model = edmf.sde_model

    #####
    ##### center variables
    #####

    @inbounds for k in real_center_indices(grid)
        #####
        ##### Set primitive variables
        #####
        @inbounds for i in 1:N_up
            if is_surface_center(grid, k)
                if prog_up[i].ρarea[k] / ρ0_c[k] >= edmf.minimum_area
                    aux_up[i].θ_liq_ice[k] = edmf.h_surface_bc[i]
                    aux_up[i].q_tot[k] = edmf.qt_surface_bc[i]
                    aux_up[i].area[k] = edmf.area_surface_bc[i]
                else
                    aux_up[i].θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
                    aux_up[i].q_tot[k] = prog_gm.q_tot[k]
                end
            else
                if prog_up[i].ρarea[k] / ρ0_c[k] >= edmf.minimum_area
                    aux_up[i].θ_liq_ice[k] = prog_up[i].ρaθ_liq_ice[k] / prog_up[i].ρarea[k]
                    aux_up[i].q_tot[k] = prog_up[i].ρaq_tot[k] / prog_up[i].ρarea[k]
                    aux_up[i].area[k] = prog_up[i].ρarea[k] / ρ0_c[k]
                else
                    aux_up[i].θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
                    aux_up[i].q_tot[k] = prog_gm.q_tot[k]
                    aux_up[i].area[k] = 0
                end
            end
        end

        #####
        ##### compute bulk
        #####
        aux_bulk.q_tot[k] = 0
        aux_bulk.θ_liq_ice[k] = 0
        aux_bulk.area[k] = sum(i -> aux_up[i].area[k], 1:N_up)
        if aux_bulk.area[k] > 0
            @inbounds for i in 1:N_up
                aux_bulk.q_tot[k] += aux_up[i].area[k] * aux_up[i].q_tot[k] / aux_bulk.area[k]
                aux_bulk.θ_liq_ice[k] += aux_up[i].area[k] * aux_up[i].θ_liq_ice[k] / aux_bulk.area[k]
            end
        else
            aux_bulk.q_tot[k] = prog_gm.q_tot[k]
            aux_bulk.θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
        end
        aux_en.area[k] = 1 - aux_bulk.area[k]
        aux_en.tke[k] = prog_en.ρatke[k] / (ρ0_c[k] * aux_en.area[k])
        aux_en.Hvar[k] = prog_en.ρaHvar[k] / (ρ0_c[k] * aux_en.area[k])
        aux_en.QTvar[k] = prog_en.ρaQTvar[k] / (ρ0_c[k] * aux_en.area[k])
        aux_en.HQTcov[k] = prog_en.ρaHQTcov[k] / (ρ0_c[k] * aux_en.area[k])

        #####
        ##### decompose_environment
        #####
        a_bulk_c = aux_bulk.area[k]
        val1 = 1 / (1 - a_bulk_c)
        val2 = a_bulk_c * val1
        aux_en.q_tot[k] = max(val1 * prog_gm.q_tot[k] - val2 * aux_bulk.q_tot[k], 0) #Yair - this is here to prevent negative QT
        aux_en.θ_liq_ice[k] = val1 * prog_gm.θ_liq_ice[k] - val2 * aux_bulk.θ_liq_ice[k]

        #####
        ##### saturation_adjustment and buoyancy
        #####
        ts_gm = thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        aux_tc.θ_virt[k] = TD.virtual_pottemp(ts_gm)

        ts_en = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])

        aux_en.T[k] = TD.air_temperature(ts_en)
        aux_en.θ_virt[k] = TD.virtual_pottemp(ts_en)
        aux_en.θ_dry[k] = TD.dry_pottemp(ts_en)
        aux_en.q_liq[k] = TD.liquid_specific_humidity(ts_en)
        aux_en.q_ice[k] = TD.ice_specific_humidity(ts_en)
        rho = TD.air_density(ts_en)
        aux_en.buoy[k] = buoyancy_c(param_set, ρ0_c[k], rho)

        # update_sat_unsat
        if TD.has_condensate(ts_en)
            aux_en.cloud_fraction[k] = 1.0
            aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts_en)
            aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts_en)
            aux_en_sat.T[k] = TD.air_temperature(ts_en)
            aux_en_sat.q_tot[k] = TD.total_specific_humidity(ts_en)
            aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(ts_en)
        else
            aux_en.cloud_fraction[k] = 0.0
            aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts_en)
            aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(ts_en)
            aux_en_unsat.q_tot[k] = TD.total_specific_humidity(ts_en)
        end

        aux_en.RH[k] = TD.relative_humidity(ts_en)

        @inbounds for i in 1:N_up
            if aux_up[i].area[k] < edmf.minimum_area && k > kc_surf && aux_up[i].area[k - 1] > 0.0
                qt = aux_up[i].q_tot[k - 1]
                h = aux_up[i].θ_liq_ice[k - 1]
                ts_up = thermo_state_pθq(param_set, p0_c[k], h, qt)
            else
                ts_up = thermo_state_pθq(param_set, p0_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k])
            end
            aux_up[i].q_liq[k] = TD.liquid_specific_humidity(ts_up)
            aux_up[i].q_ice[k] = TD.ice_specific_humidity(ts_up)
            aux_up[i].T[k] = TD.air_temperature(ts_up)
            ρ = TD.air_density(ts_up)
            aux_up[i].buoy[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
            aux_up[i].RH[k] = TD.relative_humidity(ts_up)
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
    @inbounds for i in 1:N_up
        a_up_bcs = (; bottom = CCO.SetValue(edmf.area_surface_bc[i]), top = CCO.Extrapolate())
        If = CCO.InterpolateC2F(; a_up_bcs...)
        a_min = edmf.minimum_area
        a_up = aux_up[i].area
        @. aux_up_f[i].w = ifelse(If(a_up) >= a_min, max(prog_up_f[i].ρaw / (ρ0_f * If(a_up)), 0), FT(0))
    end
    @inbounds for i in 1:N_up
        aux_up_f[i].w[kf_surf] = edmf.w_surface_bc[i]
    end

    parent(aux_tc_f.bulk.w) .= 0
    a_bulk_bcs = (; bottom = CCO.SetValue(sum(edmf.area_surface_bc)), top = CCO.Extrapolate())
    Ifb = CCO.InterpolateC2F(; a_bulk_bcs...)
    @inbounds for i in 1:N_up
        a_up = aux_up[i].area
        a_up_bcs = (; bottom = CCO.SetValue(edmf.area_surface_bc[i]), top = CCO.Extrapolate())
        Ifu = CCO.InterpolateC2F(; a_up_bcs...)
        @. aux_tc_f.bulk.w += ifelse(Ifb(aux_bulk.area) > 0, Ifu(a_up) * aux_up_f[i].w / Ifb(aux_bulk.area), FT(0))
    end
    # Assuming gm.W = 0!
    @. aux_en_f.w = -Ifb(aux_bulk.area) / (1 - Ifb(aux_bulk.area)) * aux_tc_f.bulk.w

    #####
    #####  diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, grid, state, Val(:Hvar), Val(:θ_liq_ice), Val(:θ_liq_ice))
    get_GMV_CoVar(edmf, grid, state, Val(:QTvar), Val(:q_tot), Val(:q_tot))
    get_GMV_CoVar(edmf, grid, state, Val(:HQTcov), Val(:θ_liq_ice), Val(:q_tot))

    # TODO - use this inversion in free_convection_windspeed and not compute zi twice
    edmf.zi = get_inversion(grid, state, param_set, surface.Ri_bulk_crit)

    update_surface(Case, grid, state, gm, t, param_set)
    update_forcing(Case, grid, state, gm, t, param_set)
    update_radiation(Case, grid, state, gm, t, param_set)

    compute_pressure_plume_spacing(edmf, param_set)

    #####
    ##### compute_updraft_closures
    #####

    Ic = CCO.InterpolateF2C()
    ∇c = CCO.DivergenceF2C()
    LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))

    @inbounds for i in 1:N_up
        # compute ∇m at cell centers
        a_up = aux_up[i].area
        w_up = aux_up_f[i].w
        w_en = aux_en_f.w
        w_gm = prog_gm_f.w
        @. m_entr_detr = a_up * (Ic(w_up) - Ic(w_gm))
        @. ∇m_entr_detr = ∇c(wvec(LB(m_entr_detr)))
        @. w_up_c = Ic(w_up)
        @. w_en_c = Ic(w_en)
        @inbounds for k in real_center_indices(grid)
            # entrainment

            q_cond_up = TD.condensate(TD.PhasePartition(aux_up[i].q_tot[k], aux_up[i].q_liq[k], aux_up[i].q_ice[k]))
            q_cond_en = TD.condensate(TD.PhasePartition(aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k]))

            if aux_up[i].area[k] > 0.0

                εδ_model_vars = GeneralizedEntr(;
                    q_cond_up = q_cond_up,
                    q_cond_en = q_cond_en,
                    w_up = w_up_c[k],
                    w_en = w_en_c[k],
                    b_up = aux_up[i].buoy[k],
                    b_en = aux_en.buoy[k],
                    tke = aux_en.tke[k],
                    dMdz = ∇m_entr_detr[k],
                    M = m_entr_detr[k],
                    a_up = aux_up[i].area[k],
                    a_en = aux_en.area[k],
                    R_up = edmf.pressure_plume_spacing[i],
                    RH_up = aux_up[i].RH[k],
                    RH_en = aux_en.RH[k],
                    max_area = max_area,
                    zc_i = grid.zc[k].z,
                    Δt = Δt,
                )

                er = entr_detr(param_set, εδ_model_vars, edmf.entr_closure)
                aux_up[i].entr_sc[k] = er.ε_dyn
                aux_up[i].detr_sc[k] = er.δ_dyn
                # stochastic closure
                stoch_ε = stochastic_closure(param_set, sde_model, Entrainment())
                stoch_δ = stochastic_closure(param_set, sde_model, Detrainment())
                aux_up[i].entr_sc[k] *= stoch_ε
                aux_up[i].detr_sc[k] *= stoch_δ

                aux_up[i].frac_turb_entr[k] = er.ε_turb
                aux_up[i].horiz_K_eddy[k] = er.K_ε
            else
                aux_up[i].entr_sc[k] = 0.0
                aux_up[i].detr_sc[k] = 0.0
                aux_up[i].frac_turb_entr[k] = 0.0
                aux_up[i].horiz_K_eddy[k] = 0.0
            end
        end
    end

    @inbounds for i in 1:N_up

        # pressure
        b_up = aux_up[i].buoy
        a_up = aux_up[i].area
        w_up = aux_up_f[i].w
        w_en = aux_en_f.w
        updraft_top = up.updraft_top[i]
        asp_ratio = 1.0

        b_bcs = (; bottom = CCO.SetValue(b_up[kc_surf]), top = CCO.SetValue(b_up[kc_toa]))
        a_bcs = (; bottom = CCO.SetValue(edmf.area_surface_bc[i]), top = CCO.SetValue(FT(0)))
        w_bcs = (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
        bcs = (; a_up = a_bcs, b_up = b_bcs, w_up = w_bcs)

        nh_press_buoy = aux_up_f[i].nh_pressure_b
        nh_press_adv = aux_up_f[i].nh_pressure_adv
        nh_press_drag = aux_up_f[i].nh_pressure_drag

        nh_press_buoy .= nh_pressure_buoy(param_set, a_up, b_up, ρ0_f, asp_ratio, bcs)
        nh_press_adv .= nh_pressure_adv(param_set, updraft_top, a_up, ρ0_f, w_up, bcs)
        nh_press_drag .= nh_pressure_drag(param_set, updraft_top, a_up, ρ0_f, w_up, w_en, bcs)
        @. aux_up_f[i].nh_pressure = nh_press_buoy + nh_press_adv + nh_press_drag
    end

    #####
    ##### compute_eddy_diffusivities_tke
    #####

    # Subdomain exchange term
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
        ε_turb = aux_up[i].frac_turb_entr
        @. b_exch +=
            a_up * Ic(w_up) * δ_dyn / a_en * (1 / 2 * (Ic(w_up) - Ic(w_en))^2 - tke_en) -
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
    u_gm = prog_gm.u
    v_gm = prog_gm.v
    w_en = aux_en_f.w
    # compute shear
    @. Shear² = (∇c(wvec(If0(u_gm))))^2 + (∇c(wvec(If0(v_gm))))^2 + (∇c(wvec(w_en)))^2

    q_tot_en = aux_en.q_tot
    θ_liq_ice_en = aux_en.θ_liq_ice
    θ_virt_en = aux_en.θ_virt
    @. ∂qt∂z = ∇c(wvec(If0(q_tot_en)))
    @. ∂θl∂z = ∇c(wvec(If0(θ_liq_ice_en)))
    @. ∂θv∂z = ∇c(wvec(If0(θ_virt_en)))

    # Second order approximation: Use dry and cloudy environmental fields.
    cf = aux_en.cloud_fraction
    shm = copy(cf)
    parent(shm) .= shrink_mask(vec(cf))

    # Since NaN*0 ≠ 0, we need to conditionally replace
    # our gradients by their default values.
    @. ∂qt∂z_sat = ∇c(wvec(If0(aux_en_sat.q_tot)))
    @. ∂θl∂z_sat = ∇c(wvec(If0(aux_en_sat.θ_liq_ice)))
    @. ∂θv∂z_unsat = ∇c(wvec(If0(aux_en_unsat.θ_virt)))
    for k in real_center_indices(grid)
        if shm[k] == 0
            ∂qt∂z_sat[k] = ∂qt∂z[k]
            ∂θl∂z_sat[k] = ∂θl∂z[k]
            ∂θv∂z_unsat[k] = ∂θv∂z[k]
        end
    end

    @inbounds for k in real_center_indices(grid)

        # buoyancy_gradients
        if edmf.bg_closure == BuoyGradMean()
            # First order approximation: Use environmental mean fields.
            ts_en = TD.PhaseEquil_pTq(param_set, p0_c[k], aux_en.T[k], aux_en.q_tot[k])
            bg_kwargs = (;
                t_sat = aux_en.T[k],
                qv_sat = TD.vapor_specific_humidity(ts_en),
                qt_sat = aux_en.q_tot[k],
                θ_sat = aux_en.θ_dry[k],
                θ_liq_ice_sat = aux_en.θ_liq_ice[k],
                ∂θv∂z_unsat = ∂θv∂z[k],
                ∂qt∂z_sat = ∂qt∂z[k],
                ∂θl∂z_sat = ∂θl∂z[k],
                p0 = p0_c[k],
                en_cld_frac = aux_en.cloud_fraction[k],
                alpha0 = α0_c[k],
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
                p0 = p0_c[k],
                en_cld_frac = aux_en.cloud_fraction[k],
                alpha0 = α0_c[k],
            )
            bg_model = EnvBuoyGrad(edmf.bg_closure; bg_kwargs...)
        end
        bg = buoyancy_gradients(param_set, bg_model)

        # Limiting stratification scale (Deardorff, 1976)
        # compute ∇Ri and Pr
        ∇_Ri = gradient_Richardson_number(bg.∂b∂z, Shear²[k], eps(0.0))
        aux_tc.prandtl_nvec[k] = turbulent_Prandtl_number(param_set, obukhov_length, ∇_Ri)

        ml_model = MinDisspLen(;
            z = grid.zc[k].z,
            obukhov_length = obukhov_length,
            tke_surf = aux_en.tke[kc_surf],
            ustar = surface.ustar,
            Pr = aux_tc.prandtl_nvec[k],
            p0 = p0_c[k],
            ∇b = bg,
            Shear² = Shear²[k],
            tke = aux_en.tke[k],
            b_exch = b_exch[k],
        )

        ml = mixing_length(param_set, ml_model)
        aux_tc.mls[k] = ml.min_len_ind
        aux_tc.mixing_length[k] = ml.mixing_length
        aux_tc.ml_ratio[k] = ml.ml_ratio

        KM[k] = c_m * ml.mixing_length * sqrt(max(aux_en.tke[k], 0.0))
        KH[k] = KM[k] / aux_tc.prandtl_nvec[k]

        aux_en_2m.tke.buoy[k] = -aux_en.area[k] * ρ0_c[k] * KH[k] * bg.∂b∂z
    end


    #####
    ##### compute covarinaces tendendies
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
    compute_covariance_shear(edmf, grid, state, gm, Val(:tke), Val(:w), Val(:w))
    compute_covariance_shear(edmf, grid, state, gm, Val(:Hvar), Val(:θ_liq_ice), Val(:θ_liq_ice))
    compute_covariance_shear(edmf, grid, state, gm, Val(:QTvar), Val(:q_tot), Val(:q_tot))
    compute_covariance_shear(edmf, grid, state, gm, Val(:HQTcov), Val(:θ_liq_ice), Val(:q_tot))
    compute_covariance_dissipation(edmf, grid, state, Val(:tke), param_set)
    compute_covariance_dissipation(edmf, grid, state, Val(:Hvar), param_set)
    compute_covariance_dissipation(edmf, grid, state, Val(:QTvar), param_set)
    compute_covariance_dissipation(edmf, grid, state, Val(:HQTcov), param_set)

    # TODO defined again in compute_covariance_shear and compute_covaraince
    @inbounds for k in real_center_indices(grid)
        aux_en_2m.tke.rain_src[k] = 0
        aux_en_2m.Hvar.rain_src[k] = ρ0_c[k] * aux_en.area[k] * 2 * aux_en.Hvar_rain_dt[k]
        aux_en_2m.QTvar.rain_src[k] = ρ0_c[k] * aux_en.area[k] * 2 * aux_en.QTvar_rain_dt[k]
        aux_en_2m.HQTcov.rain_src[k] = ρ0_c[k] * aux_en.area[k] * aux_en.HQTcov_rain_dt[k]
    end

    get_GMV_CoVar(edmf, grid, state, Val(:tke), Val(:w), Val(:w))

    compute_diffusive_fluxes(edmf, grid, state, gm, Case, param_set)


    # TODO: use dispatch
    if edmf.precip_model isa Clima1M
        # helper to calculate the rain velocity
        # TODO: assuming gm.W = 0
        # TODO: verify translation
        term_vel_rain = aux_tc.term_vel_rain
        term_vel_snow = aux_tc.term_vel_snow
        prog_pr = center_prog_precipitation(state)

        @inbounds for k in real_center_indices(grid)
            term_vel_rain[k] = CM1.terminal_velocity(param_set, CM1.RainType(), ρ0_c[k], prog_pr.q_rai[k])
            term_vel_snow[k] = CM1.terminal_velocity(param_set, CM1.SnowType(), ρ0_c[k], prog_pr.q_sno[k])
        end
    end

    @inbounds for i in 1:N_up
        up.updraft_top[i] = 0.0

        @inbounds for k in real_center_indices(grid)
            if aux_up[i].area[k] > 1e-3
                up.updraft_top[i] = max(up.updraft_top[i], grid.zc[k])
            end
        end
    end

    return nothing
end
