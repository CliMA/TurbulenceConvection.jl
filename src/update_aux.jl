function update_aux!(edmf, gm, grid, Case, ref_state, param_set, TS)
    #####
    ##### Unpack common variables
    #####
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    en_thermo = edmf.EnvThermo
    p0_c = ref_state.p0_half
    ρ0_c = ref_state.rho0_half
    α0_c = ref_state.alpha0_half
    g = CPP.grav(param_set)
    c_m = CPEDMF.c_m(param_set)
    KM = diffusivity_m(edmf).values
    KH = diffusivity_h(edmf).values
    surface = Case.Sur
    obukhov_length = surface.obukhov_length
    FT = eltype(grid)

    up.Area.bulkvalues .= up_sum(up.Area.values)

    #####
    ##### diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, up.Area, up.H, up.H, en.H, en.H, en.Hvar, gm.H.values, gm.H.values, gm.Hvar.values)
    get_GMV_CoVar(edmf, up.Area, up.QT, up.QT, en.QT, en.QT, en.QTvar, gm.QT.values, gm.QT.values, gm.QTvar.values)
    get_GMV_CoVar(edmf, up.Area, up.H, up.QT, en.H, en.QT, en.HQTcov, gm.H.values, gm.QT.values, gm.HQTcov.values)
    GMV_third_m(edmf, gm.H_third_m, en.Hvar, en.H, up.H)
    GMV_third_m(edmf, gm.QT_third_m, en.QTvar, en.QT, up.QT)
    GMV_third_m(edmf, gm.W_third_m, en.TKE, en.W, up.W)

    #####
    ##### decompose_environment
    #####
    # Find values of environmental variables by subtracting updraft values from grid mean values
    # whichvals used to check which substep we are on--correspondingly use "GMV.SomeVar.value" (last timestep value)
    # first make sure the "bulkvalues" of the updraft variables are updated
    @inbounds for k in real_face_indices(grid)
        up.W.bulkvalues[k] = 0
        a_bulk_bcs = (; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        a_bulk_f = interpc2f(up.Area.bulkvalues, grid, k; a_bulk_bcs...)
        if a_bulk_f > 1.0e-20
            @inbounds for i in xrange(up.n_updrafts)
                a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
                a_up_f = interpc2f(up.Area.values, grid, k, i; a_up_bcs...)
                up.W.bulkvalues[k] += a_up_f * up.W.values[i, k] / a_bulk_f
            end
        end
        # Assuming gm.W = 0!
        en.W.values[k] = -a_bulk_f / (1 - a_bulk_f) * up.W.bulkvalues[k]
    end

    @inbounds for k in real_center_indices(grid)
        a_bulk_c = up.Area.bulkvalues[k]
        up.QT.bulkvalues[k] = 0
        up.QL.bulkvalues[k] = 0
        up.H.bulkvalues[k] = 0
        up.T.bulkvalues[k] = 0
        up.B.bulkvalues[k] = 0
        up.RH.bulkvalues[k] = 0
        if a_bulk_c > 1.0e-20
            @inbounds for i in xrange(up.n_updrafts)
                up.QT.bulkvalues[k] += up.Area.values[i, k] * up.QT.values[i, k] / a_bulk_c
                up.QL.bulkvalues[k] += up.Area.values[i, k] * up.QL.values[i, k] / a_bulk_c
                up.H.bulkvalues[k] += up.Area.values[i, k] * up.H.values[i, k] / a_bulk_c
                up.T.bulkvalues[k] += up.Area.values[i, k] * up.T.values[i, k] / a_bulk_c
                up.RH.bulkvalues[k] += up.Area.values[i, k] * up.RH.values[i, k] / a_bulk_c
                up.B.bulkvalues[k] += up.Area.values[i, k] * up.B.values[i, k] / a_bulk_c
            end
        else
            up.QT.bulkvalues[k] = gm.QT.values[k]
            up.H.bulkvalues[k] = gm.H.values[k]
            up.RH.bulkvalues[k] = gm.RH.values[k]
            up.T.bulkvalues[k] = gm.T.values[k]
        end
        if up.QL.bulkvalues[k] > 1e-8 && a_bulk_c > 1e-3
            up.cloud_fraction[k] = 1.0
        else
            up.cloud_fraction[k] = 0.0
        end

        val1 = 1 / (1 - a_bulk_c)
        val2 = a_bulk_c * val1

        en.Area.values[k] = 1 - a_bulk_c
        en.QT.values[k] = max(val1 * gm.QT.values[k] - val2 * up.QT.bulkvalues[k], 0) #Yair - this is here to prevent negative QT
        en.H.values[k] = val1 * gm.H.values[k] - val2 * up.H.bulkvalues[k]

        #####
        ##### saturation_adjustment
        #####

        ts_en = TD.PhaseEquil_pθq(param_set, p0_c[k], en.H.values[k], en.QT.values[k])

        en.T.values[k] = TD.air_temperature(ts_en)
        en.QL.values[k] = TD.liquid_specific_humidity(ts_en)
        rho = TD.air_density(ts_en)
        en.B.values[k] = buoyancy_c(param_set, ρ0_c[k], rho)

        update_cloud_dry(en_thermo, k, en, ts_en)
        en.RH.values[k] = TD.relative_humidity(ts_en)

        #####
        ##### buoyancy
        #####

        @inbounds for i in xrange(up.n_updrafts)
            if up.Area.values[i, k] > 0.0
                ts_up = TD.PhaseEquil_pθq(param_set, p0_c[k], up.H.values[i, k], up.QT.values[i, k])
                up.QL.values[i, k] = TD.liquid_specific_humidity(ts_up)
                up.T.values[i, k] = TD.air_temperature(ts_up)
                ρ = TD.air_density(ts_up)
                up.B.values[i, k] = buoyancy_c(param_set, ρ0_c[k], ρ)
                up.RH.values[i, k] = TD.relative_humidity(ts_up)
            elseif k > kc_surf
                if up.Area.values[i, k - 1] > 0.0 && edmf.extrapolate_buoyancy
                    qt = up.QT.values[i, k - 1]
                    h = up.H.values[i, k - 1]
                    ts_up = TD.PhaseEquil_pθq(param_set, p0_c[k], h, qt)
                    ρ = TD.air_density(ts_up)
                    up.B.values[i, k] = buoyancy_c(param_set, ρ0_c[k], ρ)
                    up.RH.values[i, k] = TD.relative_humidity(ts_up)
                else
                    up.B.values[i, k] = en.B.values[k]
                    up.RH.values[i, k] = en.RH.values[k]
                end
            else
                up.B.values[i, k] = en.B.values[k]
                up.RH.values[i, k] = en.RH.values[k]
            end
        end

        gm.B.values[k] = (1.0 - up.Area.bulkvalues[k]) * en.B.values[k]
        @inbounds for i in xrange(up.n_updrafts)
            gm.B.values[k] += up.Area.values[i, k] * up.B.values[i, k]
        end
        @inbounds for i in xrange(up.n_updrafts)
            up.B.values[i, k] -= gm.B.values[k]
        end
        en.B.values[k] -= gm.B.values[k]
    end
    # TODO - use this inversion in free_convection_windspeed and not compute zi twice
    θ_ρ = center_field(grid)
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, ref_state.p0_half[k], gm.H.values[k], gm.QT.values[k])
        θ_ρ[k] = TD.virtual_pottemp(ts)
    end
    edmf.zi = get_inversion(param_set, θ_ρ, gm.U.values, gm.V.values, grid, surface.Ri_bulk_crit)

    update_surface(Case, gm, TS)
    update_forcing(Case, gm, TS)
    update_radiation(Case, gm, TS)

    #####
    ##### update_GMV_diagnostics
    #####
    a_up_bulk = up.Area.bulkvalues
    @inbounds for k in real_center_indices(grid)
        gm.QL.values[k] = (a_up_bulk[k] * up.QL.bulkvalues[k] + (1 - a_up_bulk[k]) * en.QL.values[k])
        gm.T.values[k] = (a_up_bulk[k] * up.T.bulkvalues[k] + (1 - a_up_bulk[k]) * en.T.values[k])
        gm.B.values[k] = (a_up_bulk[k] * up.B.bulkvalues[k] + (1 - a_up_bulk[k]) * en.B.values[k])
    end
    compute_pressure_plume_spacing(edmf)

    #####
    ##### compute_updraft_closures
    #####
    upd_cloud_diagnostics(up, ref_state) # TODO: should this be moved to compute_diagnostics! ?

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in xrange(up.n_updrafts)
            # entrainment
            if up.Area.values[i, k] > 0.0
                # compute ∇m at cell centers
                a_up_c = up.Area.values[i, k]
                w_up_c = interpf2c(up.W.values, grid, k, i)
                w_gm_c = interpf2c(gm.W.values, grid, k)
                m = a_up_c * (w_up_c - w_gm_c)
                a_up_cut = ccut_upwind(up.Area.values, grid, k, i)
                w_up_cut = daul_f2c_upwind(up.W.values, grid, k, i)
                w_gm_cut = daul_f2c_upwind(gm.W.values, grid, k)
                m_cut = a_up_cut .* (w_up_cut .- w_gm_cut)
                ∇m = FT(c∇_upwind(m_cut, grid, k; bottom = SetValue(0), top = FreeBoundary()))

                w_min = 0.001

                εδ_model = MoistureDeficitEntr(;
                    q_liq_up = up.QL.values[i, k],
                    q_liq_en = en.QL.values[k],
                    w_up = interpf2c(up.W.values, grid, k, i),
                    w_en = interpf2c(en.W.values, grid, k),
                    b_up = up.B.values[i, k],
                    b_en = en.B.values[k],
                    tke = en.TKE.values[k],
                    dMdz = ∇m,
                    M = m,
                    a_up = up.Area.values[i, k],
                    a_en = en.Area.values[k],
                    R_up = edmf.pressure_plume_spacing[i],
                    RH_up = up.RH.values[i, k],
                    RH_en = en.RH.values[k],
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
        @inbounds for i in xrange(up.n_updrafts)

            # pressure
            a_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetValue(0))
            a_kfull = interpc2f(up.Area.values, grid, k, i; a_bcs...)
            if a_kfull > 0.0
                B = up.B.values
                b_bcs = (; bottom = SetValue(B[i, kc_surf]), top = SetValue(B[i, kc_toa]))
                b_kfull = interpc2f(up.B.values, grid, k, i; b_bcs...)
                w_cut = fcut(up.W.values, grid, k, i)
                ∇w_up = f∇(w_cut, grid, k; bottom = SetValue(0), top = SetGradient(0))
                asp_ratio = 1.0
                nh_pressure_b, nh_pressure_adv, nh_pressure_drag = perturbation_pressure(
                    param_set,
                    up.updraft_top[i],
                    a_kfull,
                    b_kfull,
                    ref_state.rho0[k],
                    up.W.values[i, k],
                    ∇w_up,
                    en.W.values[k],
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
        U_cut = ccut(gm.U.values, grid, k)
        V_cut = ccut(gm.V.values, grid, k)
        wc_en = interpf2c(en.W.values, grid, k)
        wc_up = interpf2c.(Ref(up.W.values), Ref(grid), k, 1:(up.n_updrafts))
        w_dual = dual_faces(en.W.values, grid, k)

        ∇U = c∇(U_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇V = c∇(V_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇w = ∇f2c(w_dual, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        Shear² = ∇U^2 + ∇V^2 + ∇w^2

        QT_cut = ccut(en.QT.values, grid, k)
        ∂qt∂z = c∇(QT_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        THL_cut = ccut(en.H.values, grid, k)
        ∂θl∂z = c∇(THL_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        # buoyancy_gradients
        bg_model = Tan2018(;
            qt_dry = en_thermo.qt_dry[k],
            th_dry = en_thermo.th_dry[k],
            t_cloudy = en_thermo.t_cloudy[k],
            qv_cloudy = en_thermo.qv_cloudy[k],
            qt_cloudy = en_thermo.qt_cloudy[k],
            th_cloudy = en_thermo.th_cloudy[k],
            ∂qt∂z = ∂qt∂z,
            ∂θl∂z = ∂θl∂z,
            p0 = p0_c[k],
            en_cld_frac = en.cloud_fraction.values[k],
            alpha0 = α0_c[k],
        )
        bg = buoyancy_gradients(param_set, bg_model)

        # Limiting stratification scale (Deardorff, 1976)
        p0_cut = ccut(p0_c, grid, k)
        T_cut = ccut(en.T.values, grid, k)
        QT_cut = ccut(en.QT.values, grid, k)
        QL_cut = ccut(en.QL.values, grid, k)
        ts_cut = TD.PhaseEquil_pTq.(param_set, p0_cut, T_cut, QT_cut)
        thv_cut = TD.virtual_pottemp.(ts_cut)

        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], en.H.values[k], en.QT.values[k])
        θv = TD.virtual_pottemp(ts)
        ∂θv∂z = c∇(thv_cut, grid, k; bottom = SetGradient(0), top = Extrapolate())
        # compute ∇Ri and Pr
        ∇_Ri = gradient_Richardson_number(bg.∂b∂z_θl, bg.∂b∂z_qt, Shear², eps(0.0))
        edmf.prandtl_nvec[k] = turbulent_Prandtl_number(param_set, obukhov_length, ∇_Ri)

        ml_model = MinDisspLen(;
            z = grid.zc[k].z,
            obukhov_length = obukhov_length,
            tke_surf = en.TKE.values[kc_surf],
            ustar = surface.ustar,
            Pr = edmf.prandtl_nvec[k],
            p0 = p0_c[k],
            ∂b∂z_θl = bg.∂b∂z_θl,
            Shear² = Shear²,
            ∂b∂z_qt = bg.∂b∂z_qt,
            ∂θv∂z = ∂θv∂z,
            ∂qt∂z = ∂qt∂z,
            ∂θl∂z = ∂θl∂z,
            θv = θv,
            tke = en.TKE.values[k],
            a_en = (1 - up.Area.bulkvalues[k]),
            wc_en = wc_en,
            wc_up = Tuple(wc_up),
            a_up = Tuple(up.Area.values[:, k]),
            ε_turb = Tuple(edmf.frac_turb_entr[:, k]),
            δ_dyn = Tuple(edmf.detr_sc[:, k]),
            en_cld_frac = en.cloud_fraction.values[k],
            θ_li_en = en.H.values[k],
            ql_en = en.QL.values[k],
            qt_en = en.QT.values[k],
            T_en = en.T.values[k],
            N_up = up.n_updrafts,
        )

        ml = mixing_length(param_set, ml_model)
        edmf.mls[k] = ml.min_len_ind
        edmf.mixing_length[k] = ml.mixing_length
        edmf.ml_ratio[k] = ml.ml_ratio

        KM[k] = c_m * edmf.mixing_length[k] * sqrt(max(en.TKE.values[k], 0.0))
        KH[k] = KM[k] / edmf.prandtl_nvec[k]

        en.TKE.buoy[k] = -ml_model.a_en * ρ0_c[k] * KH[k] * (bg.∂b∂z_θl + bg.∂b∂z_qt)
    end

    compute_covariance_entr(edmf, en.TKE, up.W, up.W, en.W, en.W, gm.W, gm.W)
    compute_covariance_shear(edmf, gm, en.TKE, en.W.values, en.W.values)
    compute_covariance_interdomain_src(edmf, up.Area, up.W, up.W, en.W, en.W, en.TKE)

    #####
    ##### compute_tke_pressure
    #####
    @inbounds for k in real_center_indices(grid)
        en.TKE.press[k] = 0.0
        @inbounds for i in xrange(up.n_updrafts)
            w_up_c = interpf2c(up.W.values, grid, k, i)
            w_en_c = interpf2c(en.W.values, grid, k)
            press_c = interpf2c(edmf.nh_pressure, grid, k, i)
            en.TKE.press[k] += (w_en_c - w_up_c) * press_c
        end
    end

    compute_covariance_entr(edmf, en.Hvar, up.H, up.H, en.H, en.H, gm.H, gm.H)
    compute_covariance_entr(edmf, en.QTvar, up.QT, up.QT, en.QT, en.QT, gm.QT, gm.QT)
    compute_covariance_entr(edmf, en.HQTcov, up.H, up.QT, en.H, en.QT, gm.H, gm.QT)
    compute_covariance_shear(edmf, gm, en.Hvar, en.H.values, en.H.values)
    compute_covariance_shear(edmf, gm, en.QTvar, en.QT.values, en.QT.values)
    compute_covariance_shear(edmf, gm, en.HQTcov, en.H.values, en.QT.values)
    compute_covariance_interdomain_src(edmf, up.Area, up.H, up.H, en.H, en.H, en.Hvar)
    compute_covariance_interdomain_src(edmf, up.Area, up.QT, up.QT, en.QT, en.QT, en.QTvar)
    compute_covariance_interdomain_src(edmf, up.Area, up.H, up.QT, en.H, en.QT, en.HQTcov)

    #####
    ##### compute_covariance_rain # need to update this one
    #####

    # TODO defined again in compute_covariance_shear and compute_covaraince
    ae = 1 .- up.Area.bulkvalues # area of environment
    @inbounds for k in real_center_indices(grid)
        en.TKE.rain_src[k] = 0
        en.Hvar.rain_src[k] = ρ0_c[k] * ae[k] * 2 * en_thermo.Hvar_rain_dt[k]
        en.QTvar.rain_src[k] = ρ0_c[k] * ae[k] * 2 * en_thermo.QTvar_rain_dt[k]
        en.HQTcov.rain_src[k] = ρ0_c[k] * ae[k] * en_thermo.HQTcov_rain_dt[k]
    end

    reset_surface_covariance(edmf, gm, Case)

    compute_diffusive_fluxes(edmf, gm, Case, TS)
    update_cloud_frac(edmf, gm)
end
