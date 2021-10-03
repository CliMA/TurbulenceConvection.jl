function update_aux!(edmf, gm, grid, Case, ref_state, param_set, TS)
    kc_surf = kc_surface(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    en_thermo = edmf.EnvThermo

    #####
    ##### decompose_environment
    #####
    # Find values of environmental variables by subtracting updraft values from grid mean values
    # whichvals used to check which substep we are on--correspondingly use "GMV.SomeVar.value" (last timestep value)
    # first make sure the "bulkvalues" of the updraft variables are updated

    up.Area.bulkvalues .= up_sum(up.Area.values)

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
        ρ0_c = ref_state.rho0_half[k]
        p0_c = ref_state.p0_half[k]
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

        ts_en = TD.PhaseEquil_pθq(param_set, p0_c, en.H.values[k], en.QT.values[k])

        en.T.values[k] = TD.air_temperature(ts_en)
        en.QL.values[k] = TD.liquid_specific_humidity(ts_en)
        rho = TD.air_density(ts_en)
        en.B.values[k] = buoyancy_c(param_set, ρ0_c, rho)

        # TODO: can we pass `ts_en` here instead?
        update_cloud_dry(
            en_thermo,
            k,
            en,
            en.T.values[k],
            en.H.values[k],
            en.QT.values[k],
            en.QL.values[k],
            en.QT.values[k] - en.QL.values[k],
        )
        en.RH.values[k] = TD.relative_humidity(ts_en)

        #####
        ##### buoyancy
        #####

        @inbounds for i in xrange(up.n_updrafts)
            if up.Area.values[i, k] > 0.0
                ts_up = TD.PhaseEquil_pθq(param_set, p0_c, up.H.values[i, k], up.QT.values[i, k])
                up.QL.values[i, k] = TD.liquid_specific_humidity(ts_up)
                up.T.values[i, k] = TD.air_temperature(ts_up)
                ρ = TD.air_density(ts_up)
                up.B.values[i, k] = buoyancy_c(param_set, ρ0_c, ρ)
                up.RH.values[i, k] = TD.relative_humidity(ts_up)
            elseif k > kc_surf
                if up.Area.values[i, k - 1] > 0.0 && edmf.extrapolate_buoyancy
                    qt = up.QT.values[i, k - 1]
                    h = up.H.values[i, k - 1]
                    ts_up = TD.PhaseEquil_pθq(param_set, p0_c, h, qt)
                    ρ = TD.air_density(ts_up)
                    up.B.values[i, k] = buoyancy_c(param_set, ρ0_c, ρ)
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

    #####
    ##### diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, up.Area, up.H, up.H, en.H, en.H, en.Hvar, gm.H.values, gm.H.values, gm.Hvar.values)
    get_GMV_CoVar(edmf, up.Area, up.QT, up.QT, en.QT, en.QT, en.QTvar, gm.QT.values, gm.QT.values, gm.QTvar.values)
    get_GMV_CoVar(edmf, up.Area, up.H, up.QT, en.H, en.QT, en.HQTcov, gm.H.values, gm.QT.values, gm.HQTcov.values)
    GMV_third_m(edmf, gm.H_third_m, en.Hvar, en.H, up.H)
    GMV_third_m(edmf, gm.QT_third_m, en.QTvar, en.QT, up.QT)
    GMV_third_m(edmf, gm.W_third_m, en.TKE, en.W, up.W)
    compute_pressure_plume_spacing(edmf)
    compute_updraft_closures(edmf, gm, Case)
    compute_eddy_diffusivities_tke(edmf, gm, Case)

    #####
    ##### compute_covariance_rhs
    #####
    compute_covariance_entr(edmf, en.TKE, up.W, up.W, en.W, en.W, gm.W, gm.W)
    compute_covariance_shear(edmf, gm, en.TKE, en.W.values, en.W.values)
    compute_covariance_interdomain_src(edmf, up.Area, up.W, up.W, en.W, en.W, en.TKE)
    compute_tke_pressure(edmf)
    compute_covariance_entr(edmf, en.Hvar, up.H, up.H, en.H, en.H, gm.H, gm.H)
    compute_covariance_entr(edmf, en.QTvar, up.QT, up.QT, en.QT, en.QT, gm.QT, gm.QT)
    compute_covariance_entr(edmf, en.HQTcov, up.H, up.QT, en.H, en.QT, gm.H, gm.QT)
    compute_covariance_shear(edmf, gm, en.Hvar, en.H.values, en.H.values)
    compute_covariance_shear(edmf, gm, en.QTvar, en.QT.values, en.QT.values)
    compute_covariance_shear(edmf, gm, en.HQTcov, en.H.values, en.QT.values)
    compute_covariance_interdomain_src(edmf, up.Area, up.H, up.H, en.H, en.H, en.Hvar)
    compute_covariance_interdomain_src(edmf, up.Area, up.QT, up.QT, en.QT, en.QT, en.QTvar)
    compute_covariance_interdomain_src(edmf, up.Area, up.H, up.QT, en.H, en.QT, en.HQTcov)
    compute_covariance_rain(edmf, TS) # need to update this one
    reset_surface_covariance(edmf, gm, Case)

    compute_diffusive_fluxes(edmf, gm, Case, TS)
end
