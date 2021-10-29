function update_aux!(edmf, gm, grid, state, Case, param_set, TS)
    #####
    ##### Unpack common variables
    #####
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    en_thermo = edmf.EnvThermo
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
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)

    #####
    ##### Set primitive variables
    #####

    @inbounds for i in 1:(up.n_updrafts)

        # at the surface
        if prog_up[i].ρarea[kc_surf] / ρ0_c[kc_surf] >= edmf.minimum_area
            aux_up[i].θ_liq_ice[kc_surf] = edmf.h_surface_bc[i]
            aux_up[i].q_tot[kc_surf] = edmf.qt_surface_bc[i]
            aux_up[i].area[kc_surf] = edmf.area_surface_bc[i]
            aux_up_f[i].w[kf_surf] = edmf.w_surface_bc[i]
        else
            aux_up[i].θ_liq_ice[kc_surf] = prog_gm.θ_liq_ice[kc_surf]
            aux_up[i].q_tot[kc_surf] = prog_gm.q_tot[kc_surf]
        end

        @inbounds for k in real_center_indices(grid)
            is_surface_center(grid, k) && continue
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

    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            anew_k = interpc2f(aux_up[i].area, grid, k; a_up_bcs...)
            if anew_k >= edmf.minimum_area
                aux_up_f[i].w[k] = max(prog_up_f[i].ρaw[k] / (ρ0_f[k] * anew_k), 0)
            else
                aux_up_f[i].w[k] = 0
            end
        end
    end

    for k in real_center_indices(grid)
        aux_tc.bulk.area[k] = sum(ntuple(i -> aux_up[i].area[k], up.n_updrafts))
    end

    #####
    ##### diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, grid, state, :Hvar, :θ_liq_ice)
    get_GMV_CoVar(edmf, grid, state, :QTvar, :q_tot)
    get_GMV_CoVar(edmf, grid, state, :HQTcov, :θ_liq_ice, :q_tot)
    GMV_third_m(edmf, grid, state, :Hvar, :θ_liq_ice, :H_third_m)
    GMV_third_m(edmf, grid, state, :QTvar, :q_tot, :QT_third_m)
    GMV_third_m(edmf, grid, state, :tke, :w, :W_third_m)

    @inbounds for k in real_center_indices(grid)
        θ_liq_ice = prog_gm.θ_liq_ice[k]
        q_tot = prog_gm.q_tot[k]
        ts = thermo_state_pθq(param_set, p0_c[k], θ_liq_ice, q_tot)
        aux_gm.RH[k] = TD.relative_humidity(ts)
    end

    #####
    ##### decompose_environment
    #####
    # (Find values of environmental variables by subtracting updraft values from grid mean values.
    # Make sure the "bulkvalues" of the updraft variables are updated first)
    # velocity (face indicies)
    @inbounds for k in real_face_indices(grid)
        aux_tc_f.bulk.w[k] = 0
        a_bulk_bcs = (; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        a_bulk_f = interpc2f(aux_tc.bulk.area, grid, k; a_bulk_bcs...)
        if a_bulk_f > 1.0e-20
            @inbounds for i in 1:(up.n_updrafts)
                a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
                a_up_f = interpc2f(aux_up[i].area, grid, k; a_up_bcs...)
                aux_tc_f.bulk.w[k] += a_up_f * aux_up_f[i].w[k] / a_bulk_f
            end
        end
        # Assuming gm.W = 0!
        aux_en_f.w[k] = -a_bulk_f / (1 - a_bulk_f) * aux_tc_f.bulk.w[k]
    end
    @inbounds for k in real_center_indices(grid)
        a_bulk_c = aux_tc.bulk.area[k]
        aux_tc.bulk.q_tot[k] = 0
        aux_tc.bulk.q_liq[k] = 0
        aux_tc.bulk.q_ice[k] = 0
        aux_tc.bulk.θ_liq_ice[k] = 0
        aux_tc.bulk.T[k] = 0
        aux_tc.bulk.RH[k] = 0
        aux_tc.bulk.buoy[k] = 0

        @inbounds for i in 1:(up.n_updrafts)
            if aux_up[i].area[k] > 0.0
                ts_up = thermo_state_pθq(param_set, p0_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k])
                aux_up[i].q_liq[k] = TD.liquid_specific_humidity(ts_up)
                aux_up[i].q_ice[k] = TD.ice_specific_humidity(ts_up)
                aux_up[i].T[k] = TD.air_temperature(ts_up)
                ρ = TD.air_density(ts_up)
                aux_up[i].buoy[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
                aux_up[i].RH[k] = TD.relative_humidity(ts_up)
            elseif k > kc_surf
                if aux_up[i].area[k - 1] > 0.0 && edmf.extrapolate_buoyancy
                    qt = aux_up[i].q_tot[k - 1]
                    h = aux_up[i].θ_liq_ice[k - 1]
                    ts_up = thermo_state_pθq(param_set, p0_c[k], h, qt)
                    ρ = TD.air_density(ts_up)
                    aux_up[i].buoy[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
                    aux_up[i].RH[k] = TD.relative_humidity(ts_up)
                else
                    aux_up[i].buoy[k] = aux_en.buoy[k]
                    aux_up[i].RH[k] = aux_en.RH[k]
                end
            else
                aux_up[i].buoy[k] = aux_en.buoy[k]
                aux_up[i].RH[k] = aux_en.RH[k]
            end
        end

        if a_bulk_c > 1.0e-20
            @inbounds for i in 1:(up.n_updrafts)
                aux_tc.bulk.q_tot[k] += aux_up[i].area[k] * aux_up[i].q_tot[k] / a_bulk_c
                aux_tc.bulk.q_liq[k] += aux_up[i].area[k] * aux_up[i].q_liq[k] / a_bulk_c
                aux_tc.bulk.q_ice[k] += aux_up[i].area[k] * aux_up[i].q_ice[k] / a_bulk_c
                aux_tc.bulk.θ_liq_ice[k] += aux_up[i].area[k] * aux_up[i].θ_liq_ice[k] / a_bulk_c
                aux_tc.bulk.T[k] += aux_up[i].area[k] * aux_up[i].T[k] / a_bulk_c
                aux_tc.bulk.RH[k] += aux_up[i].area[k] * aux_up[i].RH[k] / a_bulk_c
                aux_tc.bulk.buoy[k] += aux_up[i].area[k] * aux_up[i].buoy[k] / a_bulk_c
            end
        else
            aux_tc.bulk.q_tot[k] = prog_gm.q_tot[k]
            aux_tc.bulk.θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
            aux_tc.bulk.RH[k] = aux_gm.RH[k]  # TODO - here we are using previous timestep values
            aux_tc.bulk.T[k] = aux_gm.T[k]    # TODO - here we are using previous timestep values
        end
        if TD.has_condensate(aux_tc.bulk.q_liq[k] + aux_tc.bulk.q_ice[k]) && a_bulk_c > 1e-3
            up.cloud_fraction[k] = 1.0
        else
            up.cloud_fraction[k] = 0.0
        end

        val1 = 1 / (1 - a_bulk_c)
        val2 = a_bulk_c * val1

        aux_en.area[k] = 1 - a_bulk_c
        aux_en.q_tot[k] = max(val1 * prog_gm.q_tot[k] - val2 * aux_tc.bulk.q_tot[k], 0) #Yair - this is here to prevent negative QT
        aux_en.θ_liq_ice[k] = val1 * prog_gm.θ_liq_ice[k] - val2 * aux_tc.bulk.θ_liq_ice[k]

        #####
        ##### saturation_adjustment
        #####
        ts_en = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])

        aux_en.T[k] = TD.air_temperature(ts_en)
        aux_en.q_liq[k] = TD.liquid_specific_humidity(ts_en)
        aux_en.q_ice[k] = TD.ice_specific_humidity(ts_en)
        rho = TD.air_density(ts_en)
        aux_en.buoy[k] = buoyancy_c(param_set, ρ0_c[k], rho)

        update_sat_unsat(en_thermo, state, k, ts_en)
        aux_en.RH[k] = TD.relative_humidity(ts_en)

        #####
        ##### buoyancy
        #####
        aux_gm.buoy[k] = (1.0 - aux_tc.bulk.area[k]) * aux_en.buoy[k]
        @inbounds for i in 1:(up.n_updrafts)
            aux_gm.buoy[k] += aux_up[i].area[k] * aux_up[i].buoy[k]
        end
        @inbounds for i in 1:(up.n_updrafts)
            aux_up[i].buoy[k] -= aux_gm.buoy[k]
        end
        aux_en.buoy[k] -= aux_gm.buoy[k]
    end
    # TODO - use this inversion in free_convection_windspeed and not compute zi twice
    θ_ρ = center_field(grid)
    @inbounds for k in real_center_indices(grid)
        ts = thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        θ_ρ[k] = TD.virtual_pottemp(ts)
    end
    edmf.zi = get_inversion(param_set, θ_ρ, prog_gm.u, prog_gm.v, grid, surface.Ri_bulk_crit)

    update_surface(Case, grid, state, gm, TS, param_set)
    update_forcing(Case, grid, state, gm, TS, param_set)
    update_radiation(Case, grid, state, gm, TS, param_set)

    #####
    ##### update_GMV_diagnostics
    #####
    a_up_bulk = aux_tc.bulk.area
    @inbounds for k in real_center_indices(grid)
        # TODO: Can we use thermo states here? Is this consistent with one constructed
        #       from a thermo state?
        aux_gm.q_liq[k] = (a_up_bulk[k] * aux_tc.bulk.q_liq[k] + (1 - a_up_bulk[k]) * aux_en.q_liq[k])
        aux_gm.q_ice[k] = (a_up_bulk[k] * aux_tc.bulk.q_ice[k] + (1 - a_up_bulk[k]) * aux_en.q_ice[k])
    end
    compute_pressure_plume_spacing(edmf, param_set)

    #####
    ##### compute_updraft_closures
    #####

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            # entrainment
            if aux_up[i].area[k] > 0.0
                # compute ∇m at cell centers
                a_up_c = aux_up[i].area[k]
                w_up_c = interpf2c(aux_up_f[i].w, grid, k)
                w_gm_c = interpf2c(prog_gm_f.w, grid, k)
                m = a_up_c * (w_up_c - w_gm_c)
                a_up_cut = ccut_upwind(aux_up[i].area, grid, k)
                w_up_cut = daul_f2c_upwind(aux_up_f[i].w, grid, k)
                w_gm_cut = daul_f2c_upwind(prog_gm_f.w, grid, k)
                m_cut = a_up_cut .* (w_up_cut .- w_gm_cut)
                ∇m = FT(c∇_upwind(m_cut, grid, k; bottom = SetValue(0), top = FreeBoundary()))

                w_min = 0.001

                εδ_model = MoistureDeficitEntr(;
                    q_cond_up = TD.condensate(TD.PhasePartition(
                        aux_up[i].q_tot[k],
                        aux_up[i].q_liq[k],
                        aux_up[i].q_ice[k],
                    )),
                    q_cond_en = TD.condensate(TD.PhasePartition(aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k])),
                    w_up = interpf2c(aux_up_f[i].w, grid, k),
                    w_en = interpf2c(aux_en_f.w, grid, k),
                    b_up = aux_up[i].buoy[k],
                    b_en = aux_en.buoy[k],
                    tke = prog_en.tke[k],
                    dMdz = ∇m,
                    M = m,
                    a_up = aux_up[i].area[k],
                    a_en = aux_en.area[k],
                    R_up = edmf.pressure_plume_spacing[i],
                    RH_up = aux_up[i].RH[k],
                    RH_en = aux_en.RH[k],
                    max_area = edmf.max_area,
                )

                er = entr_detr(param_set, εδ_model)
                edmf.entr_sc[i, k] = er.ε_dyn
                edmf.detr_sc[i, k] = er.δ_dyn
                # stochastic closure
                sde_model = edmf.sde_model
                stoch_ε = stochastic_closure(param_set, sde_model, Entrainment())
                stoch_δ = stochastic_closure(param_set, sde_model, Detrainment())
                edmf.entr_sc[i, k] *= stoch_ε
                edmf.detr_sc[i, k] *= stoch_δ

                edmf.frac_turb_entr[i, k] = er.ε_turb
                edmf.horiz_K_eddy[i, k] = er.K_ε
            else
                edmf.entr_sc[i, k] = 0.0
                edmf.detr_sc[i, k] = 0.0
                edmf.frac_turb_entr[i, k] = 0.0
                edmf.horiz_K_eddy[i, k] = 0.0
            end
        end
    end

    @inbounds for k in real_face_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)

            # pressure
            a_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetValue(0))
            a_kfull = interpc2f(aux_up[i].area, grid, k; a_bcs...)
            if a_kfull > 0.0
                B = aux_up[i].buoy
                b_bcs = (; bottom = SetValue(B[kc_surf]), top = SetValue(B[kc_toa]))
                b_kfull = interpc2f(aux_up[i].buoy, grid, k; b_bcs...)
                w_cut = fcut(aux_up_f[i].w, grid, k)
                ∇w_up = f∇(w_cut, grid, k; bottom = SetValue(0), top = SetGradient(0))
                asp_ratio = 1.0
                nh_pressure_b, nh_pressure_adv, nh_pressure_drag = perturbation_pressure(
                    param_set,
                    up.updraft_top[i],
                    a_kfull,
                    b_kfull,
                    ρ0_f[k],
                    aux_up_f[i].w[k],
                    ∇w_up,
                    aux_en_f.w[k],
                    asp_ratio,
                )
            else
                nh_pressure_b = 0.0
                nh_pressure_adv = 0.0
                nh_pressure_drag = 0.0
            end
            edmf.nh_pressure_b[i, k] = nh_pressure_b
            edmf.nh_pressure_adv[i, k] = nh_pressure_adv
            edmf.nh_pressure_drag[i, k] = nh_pressure_drag
            edmf.nh_pressure[i, k] = nh_pressure_b + nh_pressure_adv + nh_pressure_drag
        end
    end

    #####
    ##### compute_eddy_diffusivities_tke
    #####

    @inbounds for k in real_center_indices(grid)

        # compute shear
        U_cut = ccut(prog_gm.u, grid, k)
        V_cut = ccut(prog_gm.v, grid, k)
        wc_en = interpf2c(aux_en_f.w, grid, k)
        wc_up = ntuple(up.n_updrafts) do i
            interpf2c(aux_up_f[i].w, grid, k)
        end
        w_dual = dual_faces(aux_en_f.w, grid, k)

        ∇U = c∇(U_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇V = c∇(V_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇w = ∇f2c(w_dual, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        Shear² = ∇U^2 + ∇V^2 + ∇w^2

        QT_cut = ccut(aux_en.q_tot, grid, k)
        ∂qt∂z = c∇(QT_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        θ_liq_ice_cut = ccut(aux_en.θ_liq_ice, grid, k)
        ∂θl∂z = c∇(θ_liq_ice_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))

        p0_cut = ccut(p0_c, grid, k)
        T_cut = ccut(aux_en.T, grid, k)
        QT_cut = ccut(aux_en.q_tot, grid, k)
        ts_cut = TD.PhaseEquil_pTq.(param_set, p0_cut, T_cut, QT_cut)
        θv_cut = TD.virtual_pottemp.(ts_cut)

        ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])
        θ = TD.dry_pottemp(ts)
        θv = TD.virtual_pottemp(ts)
        ∂θv∂z = c∇(θv_cut, grid, k; bottom = SetGradient(0), top = Extrapolate())

        # buoyancy_gradients
        if edmf.bg_closure == BuoyGradMean()
            # First order approximation: Use environmental mean fields.
            bg_kwargs = (;
                t_sat = aux_en.T[k],
                qv_sat = TD.vapor_specific_humidity(ts),
                qt_sat = aux_en.q_tot[k],
                θ_sat = θ,
                θ_liq_ice_sat = aux_en.θ_liq_ice[k],
                ∂θv∂z_unsat = ∂θv∂z,
                ∂qt∂z_sat = ∂qt∂z,
                ∂θl∂z_sat = ∂θl∂z,
                p0 = p0_c[k],
                en_cld_frac = aux_en.cloud_fraction[k],
                alpha0 = α0_c[k],
            )
            bg_model = EnvBuoyGrad(edmf.bg_closure; bg_kwargs...)

        elseif edmf.bg_closure == BuoyGradQuadratures()
            # Second order approximation: Use dry and cloudy environmental fields.
            cf_cut = ccut(aux_en.cloud_fraction, grid, k)
            QT_sat_cut = ccut(en_thermo.qt_sat, grid, k)
            ∂qt∂z_sat = c∇_vanishing_subdomain(
                QT_sat_cut,
                cf_cut,
                ∂qt∂z,
                grid,
                k;
                bottom = SetGradient(0),
                top = SetGradient(0),
            )
            θ_liq_ice_sat_cut = ccut(en_thermo.θ_liq_ice_sat, grid, k)
            ∂θl∂z_sat = c∇_vanishing_subdomain(
                θ_liq_ice_sat_cut,
                cf_cut,
                ∂θl∂z,
                grid,
                k;
                bottom = SetGradient(0),
                top = SetGradient(0),
            )
            θv_unsat_cut = ccut(en_thermo.θv_unsat, grid, k)
            ∂θv∂z_unsat = c∇_vanishing_subdomain(
                θv_unsat_cut,
                cf_cut,
                ∂θv∂z,
                grid,
                k;
                bottom = SetGradient(0),
                top = SetGradient(0),
            )

            bg_kwargs = (;
                t_sat = en_thermo.t_sat[k],
                qv_sat = en_thermo.qv_sat[k],
                qt_sat = en_thermo.qt_sat[k],
                θ_sat = en_thermo.θ_sat[k],
                θ_liq_ice_sat = en_thermo.θ_liq_ice_sat[k],
                ∂θv∂z_unsat = ∂θv∂z_unsat,
                ∂qt∂z_sat = ∂qt∂z_sat,
                ∂θl∂z_sat = ∂θl∂z_sat,
                p0 = p0_c[k],
                en_cld_frac = aux_en.cloud_fraction[k],
                alpha0 = α0_c[k],
            )
            bg_model = EnvBuoyGrad(edmf.bg_closure; bg_kwargs...)
        end
        bg = buoyancy_gradients(param_set, bg_model)

        # Limiting stratification scale (Deardorff, 1976)
        # compute ∇Ri and Pr
        ∇_Ri = gradient_Richardson_number(bg.∂b∂z, Shear², eps(0.0))
        edmf.prandtl_nvec[k] = turbulent_Prandtl_number(param_set, obukhov_length, ∇_Ri)

        ml_model = MinDisspLen(;
            z = grid.zc[k].z,
            obukhov_length = obukhov_length,
            tke_surf = prog_en.tke[kc_surf],
            ustar = surface.ustar,
            Pr = edmf.prandtl_nvec[k],
            p0 = p0_c[k],
            ∇b = bg,
            Shear² = Shear²,
            tke = prog_en.tke[k],
            a_en = (1 - aux_tc.bulk.area[k]),
            wc_en = wc_en,
            wc_up = Tuple(wc_up),
            a_up = ntuple(i -> aux_up[i].area[k], up.n_updrafts),
            ε_turb = ntuple(i -> edmf.frac_turb_entr[i, k], up.n_updrafts),
            δ_dyn = ntuple(i -> edmf.detr_sc[i, k], up.n_updrafts),
            N_up = up.n_updrafts,
        )

        ml = mixing_length(param_set, ml_model)
        edmf.mls[k] = ml.min_len_ind
        edmf.mixing_length[k] = ml.mixing_length
        edmf.ml_ratio[k] = ml.ml_ratio

        KM[k] = c_m * edmf.mixing_length[k] * sqrt(max(prog_en.tke[k], 0.0))
        KH[k] = KM[k] / edmf.prandtl_nvec[k]

        aux_en_2m.tke.buoy[k] = -ml_model.a_en * ρ0_c[k] * KH[k] * bg.∂b∂z
    end

    compute_covariance_entr(edmf, grid, state, :tke, :w)
    compute_covariance_shear(edmf, grid, state, gm, :tke, :w)
    compute_covariance_interdomain_src(edmf, grid, state, :tke, :w)

    #####
    ##### compute_tke_pressure
    #####
    @inbounds for k in real_center_indices(grid)
        aux_en_2m.tke.press[k] = 0.0
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(aux_up_f[i].w, grid, k)
            w_en_c = interpf2c(aux_en_f.w, grid, k)
            press_c = interpf2c(edmf.nh_pressure, grid, k, i)
            aux_en_2m.tke.press[k] += (w_en_c - w_up_c) * press_c
        end
    end

    compute_covariance_entr(edmf, grid, state, :Hvar, :θ_liq_ice)
    compute_covariance_entr(edmf, grid, state, :QTvar, :q_tot)
    compute_covariance_entr(edmf, grid, state, :HQTcov, :θ_liq_ice, :q_tot)
    compute_covariance_shear(edmf, grid, state, gm, :Hvar, :θ_liq_ice)
    compute_covariance_shear(edmf, grid, state, gm, :QTvar, :q_tot)
    compute_covariance_shear(edmf, grid, state, gm, :HQTcov, :θ_liq_ice, :q_tot)
    compute_covariance_interdomain_src(edmf, grid, state, :Hvar, :θ_liq_ice)
    compute_covariance_interdomain_src(edmf, grid, state, :QTvar, :q_tot)
    compute_covariance_interdomain_src(edmf, grid, state, :HQTcov, :θ_liq_ice, :q_tot)

    #####
    ##### compute_covariance_rain # need to update this one
    #####

    # TODO defined again in compute_covariance_shear and compute_covaraince
    ae = 1 .- aux_tc.bulk.area # area of environment
    @inbounds for k in real_center_indices(grid)
        aux_en_2m.tke.rain_src[k] = 0
        aux_en_2m.Hvar.rain_src[k] = ρ0_c[k] * ae[k] * 2 * en_thermo.Hvar_rain_dt[k]
        aux_en_2m.QTvar.rain_src[k] = ρ0_c[k] * ae[k] * 2 * en_thermo.QTvar_rain_dt[k]
        aux_en_2m.HQTcov.rain_src[k] = ρ0_c[k] * ae[k] * en_thermo.HQTcov_rain_dt[k]
    end

    get_GMV_CoVar(edmf, grid, state, :tke, :w)

    compute_diffusive_fluxes(edmf, grid, state, gm, Case, TS, param_set)
    update_cloud_frac(edmf, grid, state, gm)

    compute_covariance_dissipation(edmf, grid, state, :tke, param_set)
    compute_covariance_detr(edmf, grid, state, :tke)

    compute_covariance_dissipation(edmf, grid, state, :Hvar, param_set)
    compute_covariance_dissipation(edmf, grid, state, :QTvar, param_set)
    compute_covariance_dissipation(edmf, grid, state, :HQTcov, param_set)
    compute_covariance_detr(edmf, grid, state, :Hvar)
    compute_covariance_detr(edmf, grid, state, :QTvar)
    compute_covariance_detr(edmf, grid, state, :HQTcov)

end
