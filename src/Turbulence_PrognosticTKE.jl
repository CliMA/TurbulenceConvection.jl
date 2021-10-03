
function initialize(
    edmf::EDMF_PrognosticTKE,
    Case::CasesBase,
    gm::GridMeanVariables,
    ref_state::ReferenceState,
    TS::TimeStepping,
)
    if Case.casename == "DryBubble"
        initialize_DryBubble(edmf, edmf.UpdVar, gm, ref_state)
    else
        initialize(edmf, edmf.UpdVar, gm)
    end
    param_set = parameter_set(edmf)
    grid = get_grid(edmf)
    update_aux!(edmf, gm, grid, Case, ref_state, param_set, TS)
    update_inversion(edmf, gm, Case.inversion_option)
    edmf.wstar = get_wstar(Case.Sur.bflux, edmf.zi)
    microphysics(edmf.EnvThermo, edmf.EnvVar, edmf.Rain, TS.dt)
    initialize_covariance(edmf, gm, Case)
    return
end

# Initialize the IO pertaining to this class
function initialize_io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)

    initialize_io(edmf.UpdVar, Stats)
    initialize_io(edmf.EnvVar, Stats)
    initialize_io(edmf.Rain, Stats)

    add_profile(Stats, "eddy_viscosity")
    add_profile(Stats, "eddy_diffusivity")
    add_profile(Stats, "entrainment_sc")
    add_profile(Stats, "detrainment_sc")
    add_profile(Stats, "nh_pressure")
    add_profile(Stats, "nh_pressure_adv")
    add_profile(Stats, "nh_pressure_drag")
    add_profile(Stats, "nh_pressure_b")
    add_profile(Stats, "asp_ratio")
    add_profile(Stats, "horiz_K_eddy")
    add_profile(Stats, "sorting_function")
    add_profile(Stats, "b_mix")
    add_ts(Stats, "rd")
    add_profile(Stats, "turbulent_entrainment")
    add_profile(Stats, "turbulent_entrainment_full")
    add_profile(Stats, "turbulent_entrainment_W")
    add_profile(Stats, "turbulent_entrainment_H")
    add_profile(Stats, "turbulent_entrainment_QT")
    add_profile(Stats, "massflux")
    add_profile(Stats, "massflux_h")
    add_profile(Stats, "massflux_qt")
    add_profile(Stats, "massflux_tendency_h")
    add_profile(Stats, "massflux_tendency_qt")
    add_profile(Stats, "diffusive_flux_h")
    add_profile(Stats, "diffusive_flux_u")
    add_profile(Stats, "diffusive_flux_v")
    add_profile(Stats, "diffusive_flux_qt")
    add_profile(Stats, "diffusive_tendency_h")
    add_profile(Stats, "diffusive_tendency_qt")
    add_profile(Stats, "total_flux_h")
    add_profile(Stats, "total_flux_qt")
    add_profile(Stats, "mixing_length")
    add_profile(Stats, "updraft_qt_precip")
    add_profile(Stats, "updraft_thetal_precip")
    # Diff mixing lengths: Ignacio
    add_profile(Stats, "ed_length_scheme")
    add_profile(Stats, "mixing_length_ratio")
    add_profile(Stats, "entdet_balance_length")
    add_profile(Stats, "interdomain_tke_t")
    add_profile(Stats, "tke_buoy")
    add_profile(Stats, "tke_dissipation")
    add_profile(Stats, "tke_entr_gain")
    add_profile(Stats, "tke_detr_loss")
    add_profile(Stats, "tke_shear")
    add_profile(Stats, "tke_pressure")
    add_profile(Stats, "tke_interdomain")
    add_profile(Stats, "Hvar_dissipation")
    add_profile(Stats, "QTvar_dissipation")
    add_profile(Stats, "HQTcov_dissipation")
    add_profile(Stats, "Hvar_entr_gain")
    add_profile(Stats, "QTvar_entr_gain")
    add_profile(Stats, "Hvar_detr_loss")
    add_profile(Stats, "QTvar_detr_loss")
    add_profile(Stats, "HQTcov_detr_loss")
    add_profile(Stats, "HQTcov_entr_gain")
    add_profile(Stats, "Hvar_shear")
    add_profile(Stats, "QTvar_shear")
    add_profile(Stats, "HQTcov_shear")
    add_profile(Stats, "Hvar_rain")
    add_profile(Stats, "QTvar_rain")
    add_profile(Stats, "HQTcov_rain")
    add_profile(Stats, "Hvar_interdomain")
    add_profile(Stats, "QTvar_interdomain")
    add_profile(Stats, "HQTcov_interdomain")
    return
end

function io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats, TS::TimeStepping)

    grid = get_grid(edmf)
    ref_state = reference_state(edmf)

    mean_nh_pressure = face_field(grid)
    mean_nh_pressure_adv = face_field(grid)
    mean_nh_pressure_drag = face_field(grid)
    mean_nh_pressure_b = face_field(grid)

    mean_asp_ratio = center_field(grid)
    mean_entr_sc = center_field(grid)
    mean_detr_sc = center_field(grid)
    massflux = center_field(grid)
    mf_h = center_field(grid)
    mf_qt = center_field(grid)
    mean_frac_turb_entr = center_field(grid)
    mean_horiz_K_eddy = center_field(grid)
    mean_sorting_function = center_field(grid)
    mean_b_mix = center_field(grid)

    io(edmf.UpdVar, Stats, ref_state)
    io(edmf.EnvVar, Stats, ref_state)
    io(edmf.Rain, Stats, ref_state, edmf.UpdThermo, edmf.EnvThermo, TS)

    a_up_bulk = edmf.UpdVar.Area.bulkvalues
    a_up = edmf.UpdVar.Area.values

    write_profile(Stats, "eddy_viscosity", diffusivity_m(edmf).values)
    write_profile(Stats, "eddy_diffusivity", diffusivity_h(edmf).values)
    write_ts(Stats, "rd", StatsBase.mean(edmf.pressure_plume_spacing))

    @inbounds for k in real_center_indices(grid)
        mf_h[k] = interpf2c(edmf.massflux_h, grid, k)
        mf_qt[k] = interpf2c(edmf.massflux_qt, grid, k)
        if edmf.UpdVar.Area.bulkvalues[k] > 0.0
            @inbounds for i in xrange(edmf.n_updrafts)
                massflux[k] += interpf2c(edmf.m, grid, k, i)
                mean_entr_sc[k] += a_up[i, k] * edmf.entr_sc[i, k] / a_up_bulk[k]
                mean_detr_sc[k] += a_up[i, k] * edmf.detr_sc[i, k] / a_up_bulk[k]
                mean_asp_ratio[k] += a_up[i, k] * edmf.asp_ratio[i, k] / a_up_bulk[k]
                mean_frac_turb_entr[k] += a_up[i, k] * edmf.frac_turb_entr[i, k] / a_up_bulk[k]
                mean_horiz_K_eddy[k] += a_up[i, k] * edmf.horiz_K_eddy[i, k] / a_up_bulk[k]
                mean_sorting_function[k] += a_up[i, k] * edmf.sorting_function[i, k] / a_up_bulk[k]
                mean_b_mix[k] += a_up[i, k] * edmf.b_mix[i, k] / a_up_bulk[k]
            end
        end
    end

    @inbounds for k in real_face_indices(grid)
        a_up_bulk_f =
            interpc2f(a_up_bulk, grid, k; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        if a_up_bulk_f > 0.0
            @inbounds for i in xrange(edmf.n_updrafts)
                a_up_f =
                    interpc2f(a_up, grid, k, i; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
                mean_nh_pressure[k] += a_up_f * edmf.nh_pressure[i, k] / a_up_bulk_f
                mean_nh_pressure_b[k] += a_up_f * edmf.nh_pressure_b[i, k] / a_up_bulk_f
                mean_nh_pressure_adv[k] += a_up_f * edmf.nh_pressure_adv[i, k] / a_up_bulk_f
                mean_nh_pressure_drag[k] += a_up_f * edmf.nh_pressure_drag[i, k] / a_up_bulk_f
            end
        end
    end

    write_profile(Stats, "turbulent_entrainment", mean_frac_turb_entr)
    write_profile(Stats, "horiz_K_eddy", mean_horiz_K_eddy)
    write_profile(Stats, "entrainment_sc", mean_entr_sc)
    write_profile(Stats, "detrainment_sc", mean_detr_sc)
    write_profile(Stats, "nh_pressure", mean_nh_pressure)
    write_profile(Stats, "nh_pressure_adv", mean_nh_pressure_adv)
    write_profile(Stats, "nh_pressure_drag", mean_nh_pressure_drag)
    write_profile(Stats, "nh_pressure_b", mean_nh_pressure_b)
    write_profile(Stats, "asp_ratio", mean_asp_ratio)
    write_profile(Stats, "massflux", massflux)
    write_profile(Stats, "massflux_h", mf_h)
    write_profile(Stats, "massflux_qt", mf_qt)
    write_profile(Stats, "massflux_tendency_h", edmf.massflux_tendency_h)
    write_profile(Stats, "massflux_tendency_qt", edmf.massflux_tendency_qt)
    write_profile(Stats, "diffusive_flux_h", edmf.diffusive_flux_h)
    write_profile(Stats, "diffusive_flux_qt", edmf.diffusive_flux_qt)
    write_profile(Stats, "diffusive_flux_u", edmf.diffusive_flux_u)
    write_profile(Stats, "diffusive_flux_v", edmf.diffusive_flux_v)
    write_profile(Stats, "diffusive_tendency_h", edmf.diffusive_tendency_h)
    write_profile(Stats, "diffusive_tendency_qt", edmf.diffusive_tendency_qt)
    write_profile(Stats, "total_flux_h", mf_h .+ edmf.diffusive_flux_h)
    write_profile(Stats, "total_flux_qt", mf_qt .+ edmf.diffusive_flux_qt)
    write_profile(Stats, "mixing_length", edmf.mixing_length)
    write_profile(Stats, "updraft_qt_precip", edmf.UpdThermo.prec_source_qt_tot)
    write_profile(Stats, "updraft_thetal_precip", edmf.UpdThermo.prec_source_h_tot)

    #Different mixing lengths : Ignacio
    write_profile(Stats, "ed_length_scheme", edmf.mls)
    write_profile(Stats, "mixing_length_ratio", edmf.ml_ratio)
    write_profile(Stats, "entdet_balance_length", edmf.l_entdet)
    write_profile(Stats, "interdomain_tke_t", edmf.b)
    compute_covariance_dissipation(edmf, edmf.EnvVar.TKE)
    write_profile(Stats, "tke_dissipation", edmf.EnvVar.TKE.dissipation)
    write_profile(Stats, "tke_entr_gain", edmf.EnvVar.TKE.entr_gain)
    compute_covariance_detr(edmf, edmf.EnvVar.TKE)
    write_profile(Stats, "tke_detr_loss", edmf.EnvVar.TKE.detr_loss)
    write_profile(Stats, "tke_shear", edmf.EnvVar.TKE.shear)
    write_profile(Stats, "tke_buoy", edmf.EnvVar.TKE.buoy)
    write_profile(Stats, "tke_pressure", edmf.EnvVar.TKE.press)
    write_profile(Stats, "tke_interdomain", edmf.EnvVar.TKE.interdomain)

    compute_covariance_dissipation(edmf, edmf.EnvVar.Hvar)
    write_profile(Stats, "Hvar_dissipation", edmf.EnvVar.Hvar.dissipation)
    compute_covariance_dissipation(edmf, edmf.EnvVar.QTvar)
    write_profile(Stats, "QTvar_dissipation", edmf.EnvVar.QTvar.dissipation)
    compute_covariance_dissipation(edmf, edmf.EnvVar.HQTcov)
    write_profile(Stats, "HQTcov_dissipation", edmf.EnvVar.HQTcov.dissipation)
    write_profile(Stats, "Hvar_entr_gain", edmf.EnvVar.Hvar.entr_gain)
    write_profile(Stats, "QTvar_entr_gain", edmf.EnvVar.QTvar.entr_gain)
    write_profile(Stats, "HQTcov_entr_gain", edmf.EnvVar.HQTcov.entr_gain)
    compute_covariance_detr(edmf, edmf.EnvVar.Hvar)
    compute_covariance_detr(edmf, edmf.EnvVar.QTvar)
    compute_covariance_detr(edmf, edmf.EnvVar.HQTcov)
    write_profile(Stats, "Hvar_detr_loss", edmf.EnvVar.Hvar.detr_loss)
    write_profile(Stats, "QTvar_detr_loss", edmf.EnvVar.QTvar.detr_loss)
    write_profile(Stats, "HQTcov_detr_loss", edmf.EnvVar.HQTcov.detr_loss)
    write_profile(Stats, "Hvar_shear", edmf.EnvVar.Hvar.shear)
    write_profile(Stats, "QTvar_shear", edmf.EnvVar.QTvar.shear)
    write_profile(Stats, "HQTcov_shear", edmf.EnvVar.HQTcov.shear)
    write_profile(Stats, "Hvar_rain", edmf.EnvVar.Hvar.rain_src)
    write_profile(Stats, "QTvar_rain", edmf.EnvVar.QTvar.rain_src)
    write_profile(Stats, "HQTcov_rain", edmf.EnvVar.HQTcov.rain_src)
    write_profile(Stats, "Hvar_interdomain", edmf.EnvVar.Hvar.interdomain)
    write_profile(Stats, "QTvar_interdomain", edmf.EnvVar.QTvar.interdomain)
    write_profile(Stats, "HQTcov_interdomain", edmf.EnvVar.HQTcov.interdomain)
    return
end

update_surface(Case::CasesBase, GMV::GridMeanVariables, TS::TimeStepping) =
    error("update_surface should be overloaded in case-specific methods (in Cases.jl)")
update_forcing(Case::CasesBase, GMV::GridMeanVariables, TS::TimeStepping) =
    error("update_forcing should be overloaded in case-specific methods (in Cases.jl)")
update_radiation(Case::CasesBase, GMV::GridMeanVariables, TS::TimeStepping) =
    error("update_radiation should be overloaded in case-specific methods (in Cases.jl)")

function update_cloud_frac(edmf::EDMF_PrognosticTKE, GMV::GridMeanVariables)
    grid = get_grid(edmf)
    # update grid-mean cloud fraction and cloud cover
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        edmf.EnvVar.Area.values[k] = 1 - edmf.UpdVar.Area.bulkvalues[k]
        GMV.cloud_fraction.values[k] =
            edmf.EnvVar.Area.values[k] * edmf.EnvVar.cloud_fraction.values[k] +
            edmf.UpdVar.Area.bulkvalues[k] * edmf.UpdVar.cloud_fraction[k]
    end
    GMV.cloud_cover = min(edmf.EnvVar.cloud_cover + sum(edmf.UpdVar.cloud_cover), 1)
end

function compute_gm_tendencies!(edmf::EDMF_PrognosticTKE, grid, Case, gm, ref_state, TS)
    gm.U.tendencies .= 0
    gm.V.tendencies .= 0
    gm.QT.tendencies .= 0
    gm.H.tendencies .= 0
    param_set = parameter_set(gm)
    ρ0_f = ref_state.rho0
    p0_c = ref_state.p0_half
    α0_c = ref_state.alpha0_half
    kf_surf = kf_surface(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    ae = 1 .- up.Area.bulkvalues # area of environment

    @inbounds for k in real_center_indices(grid)
        # Apply large-scale horizontal advection tendencies
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], gm.H.values[k], gm.QT.values[k])
        Π = TD.exner(ts)

        if Case.Fo.apply_coriolis
            gm.U.tendencies[k] -= Case.Fo.coriolis_param * (Case.Fo.vg[k] - gm.V.values[k])
            gm.V.tendencies[k] += Case.Fo.coriolis_param * (Case.Fo.ug[k] - gm.U.values[k])
        end
        if rad_type(Case.Rad) <: Union{RadiationDYCOMS_RF01, RadiationLES}
            gm.H.tendencies[k] += Case.Rad.dTdt[k] / Π
        end
        H_cut = ccut_downwind(gm.H.values, grid, k)
        q_tot_cut = ccut_downwind(gm.QT.values, grid, k)
        ∇H = c∇_downwind(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
        ∇q_tot = c∇_downwind(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))

        if force_type(Case.Fo) <: ForcingDYCOMS_RF01
            gm.QT.tendencies[k] += Case.Fo.dqtdt[k]
            # Apply large-scale subsidence tendencies
            gm.H.tendencies[k] -= ∇H * Case.Fo.subsidence[k]
            gm.QT.tendencies[k] -= ∇q_tot * Case.Fo.subsidence[k]
        end

        if force_type(Case.Fo) <: ForcingStandard
            if Case.Fo.apply_subsidence
                gm.H.tendencies[k] -= ∇H * Case.Fo.subsidence[k]
                gm.QT.tendencies[k] -= ∇q_tot * Case.Fo.subsidence[k]
            end
            gm.H.tendencies[k] += Case.Fo.dTdt[k] / Π
            gm.QT.tendencies[k] += Case.Fo.dqtdt[k]
        end

        if force_type(Case.Fo) <: ForcingLES
            H_horz_adv = Case.Fo.dtdt_hadv[k] / Π
            H_nudge = Case.Fo.dtdt_nudge[k] / Π
            H_fluc = Case.Fo.dtdt_fluc[k] / Π

            gm_U_nudge_k = (Case.Fo.u_nudge[k] - gm.U.values[k]) / Case.Fo.nudge_tau
            gm_V_nudge_k = (Case.Fo.v_nudge[k] - gm.V.values[k]) / Case.Fo.nudge_tau
            if Case.Fo.apply_subsidence
                # Apply large-scale subsidence tendencies
                gm_H_subsidence_k = -∇H * Case.Fo.subsidence[k]
                gm_QT_subsidence_k = -∇q_tot * Case.Fo.subsidence[k]
            else
                gm_H_subsidence_k = 0.0
                gm_QT_subsidence_k = 0.0
            end

            gm.H.tendencies[k] += H_horz_adv + H_nudge + H_fluc + gm_H_subsidence_k
            gm.QT.tendencies[k] +=
                Case.Fo.dqtdt_hadv[k] + Case.Fo.dqtdt_nudge[k] + gm_QT_subsidence_k + Case.Fo.dqtdt_fluc[k]

            gm.U.tendencies[k] += gm_U_nudge_k
            gm.V.tendencies[k] += gm_V_nudge_k
        end
        gm.QT.tendencies[k] +=
            edmf.UpdThermo.prec_source_qt_tot[k] +
            edmf.EnvThermo.prec_source_qt[k] +
            edmf.rainphysics.rain_evap_source_qt[k]
        gm.H.tendencies[k] +=
            edmf.UpdThermo.prec_source_h_tot[k] +
            edmf.EnvThermo.prec_source_h[k] +
            edmf.rainphysics.rain_evap_source_h[k]
    end

    edmf.massflux_h .= 0.0
    edmf.massflux_qt .= 0.0
    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in xrange(up.n_updrafts)
        edmf.m[i, kf_surf] = 0.0
        a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
        @inbounds for k in real_face_indices(grid)
            a_up = interpc2f(up.Area.values, grid, k, i; a_up_bcs...)
            a_en = interpc2f(ae, grid, k; a_up_bcs...)
            edmf.m[i, k] = ρ0_f[k] * a_up * a_en * (up.W.values[i, k] - en.W.values[k])
        end
    end

    @inbounds for k in real_face_indices(grid)
        edmf.massflux_h[k] = 0.0
        edmf.massflux_qt[k] = 0.0
        # We know that, since W = 0 at z = 0, m = 0 also, and
        # therefore θ_liq_ice / q_tot values do not matter
        m_bcs = (; bottom = SetValue(0), top = SetValue(0))
        h_en_f = interpc2f(en.H.values, grid, k; m_bcs...)
        qt_en_f = interpc2f(en.QT.values, grid, k; m_bcs...)
        @inbounds for i in xrange(up.n_updrafts)
            h_up_f = interpc2f(up.H.values, grid, k, i; m_bcs...)
            qt_up_f = interpc2f(up.QT.values, grid, k, i; m_bcs...)
            edmf.massflux_h[k] += edmf.m[i, k] * (h_up_f - h_en_f)
            edmf.massflux_qt[k] += edmf.m[i, k] * (qt_up_f - qt_en_f)
        end
    end

    # Compute the  mass flux tendencies
    # Adjust the values of the grid mean variables
    @inbounds for k in real_center_indices(grid)
        mf_tend_h_dual = dual_faces(edmf.massflux_h, grid, k)
        mf_tend_qt_dual = dual_faces(edmf.massflux_qt, grid, k)

        ∇mf_tend_h = ∇f2c(mf_tend_h_dual, grid, k)
        ∇mf_tend_qt = ∇f2c(mf_tend_qt_dual, grid, k)

        mf_tend_h = -∇mf_tend_h * α0_c[k]
        mf_tend_qt = -∇mf_tend_qt * α0_c[k]

        # Prepare the output
        edmf.massflux_tendency_h[k] = mf_tend_h
        edmf.massflux_tendency_qt[k] = mf_tend_qt
        gm.H.tendencies[k] += mf_tend_h
        gm.QT.tendencies[k] += mf_tend_qt
    end
end

function compute_diffusive_fluxes(edmf::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)
    grid = get_grid(edmf)
    ref_state = reference_state(edmf)
    param_set = parameter_set(edmf)
    edmf.ae .= 1 .- edmf.UpdVar.Area.bulkvalues # area of environment
    KM = diffusivity_m(edmf).values
    KH = diffusivity_h(edmf).values
    aeKM = edmf.ae .* KM
    aeKH = edmf.ae .* KH
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_surf = kf_surface(grid)
    aeKM_bcs = (; bottom = SetValue(aeKM[kc_surf]), top = SetValue(aeKM[kc_toa]))
    aeKH_bcs = (; bottom = SetValue(aeKH[kc_surf]), top = SetValue(aeKH[kc_toa]))

    @inbounds for k in real_face_indices(grid)
        edmf.rho_ae_KH[k] = interpc2f(aeKH, grid, k; aeKH_bcs...) * ref_state.rho0[k]
        edmf.rho_ae_KM[k] = interpc2f(aeKM, grid, k; aeKM_bcs...) * ref_state.rho0[k]
    end
    q_bc = -Case.Sur.rho_qtflux / edmf.rho_ae_KH[kf_surf]
    θ_liq_ice_bc = -Case.Sur.rho_hflux / edmf.rho_ae_KH[kf_surf]
    u_bc = -Case.Sur.rho_uflux / edmf.rho_ae_KM[kf_surf]
    v_bc = -Case.Sur.rho_vflux / edmf.rho_ae_KM[kf_surf]
    @inbounds for k in real_center_indices(grid)
        q_cut = ccut(edmf.EnvVar.QT.values, grid, k)
        ∇q_tot = c∇(q_cut, grid, k; bottom = SetGradient(q_bc), top = SetGradient(0))
        edmf.diffusive_flux_qt[k] = -0.5 * ref_state.rho0_half[k] * edmf.ae[k] * KH[k] * ∇q_tot

        θ_liq_ice_cut = ccut(edmf.EnvVar.H.values, grid, k)
        ∇θ_liq_ice = c∇(θ_liq_ice_cut, grid, k; bottom = SetGradient(θ_liq_ice_bc), top = SetGradient(0))
        edmf.diffusive_flux_h[k] = -0.5 * ref_state.rho0_half[k] * edmf.ae[k] * KH[k] * ∇θ_liq_ice

        u_cut = ccut(GMV.U.values, grid, k)
        ∇u = c∇(u_cut, grid, k; bottom = SetGradient(u_bc), top = SetGradient(0))
        edmf.diffusive_flux_u[k] = -0.5 * ref_state.rho0_half[k] * edmf.ae[k] * KM[k] * ∇u

        v_cut = ccut(GMV.V.values, grid, k)
        ∇v = c∇(v_cut, grid, k; bottom = SetGradient(v_bc), top = SetGradient(0))
        edmf.diffusive_flux_v[k] = -0.5 * ref_state.rho0_half[k] * edmf.ae[k] * KM[k] * ∇v
    end
    return
end

# Perform the update of the scheme
function update(edmf::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    gm = GMV
    up = edmf.UpdVar
    en = edmf.EnvVar
    grid = get_grid(edmf)
    ref_state = reference_state(edmf)
    param_set = parameter_set(edmf)
    up_thermo = edmf.UpdThermo
    en_thermo = edmf.EnvThermo

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.
    set_updraft_surface_bc(edmf, gm, Case)

    #####
    ##### diagnose_GMV_moments
    #####
    get_GMV_CoVar(edmf, up.Area, up.H, up.H, en.H, en.H, en.Hvar, gm.H.values, gm.H.values, gm.Hvar.values)
    get_GMV_CoVar(edmf, up.Area, up.QT, up.QT, en.QT, en.QT, en.QTvar, gm.QT.values, gm.QT.values, gm.QTvar.values)
    get_GMV_CoVar(edmf, up.Area, up.H, up.QT, en.H, en.QT, en.HQTcov, gm.H.values, gm.QT.values, gm.HQTcov.values)
    GMV_third_m(edmf, gm.H_third_m, en.Hvar, en.H, up.H)
    GMV_third_m(edmf, gm.QT_third_m, en.QTvar, en.QT, up.QT)
    GMV_third_m(edmf, gm.W_third_m, en.TKE, en.W, up.W)

    update_aux!(edmf, gm, grid, Case, ref_state, param_set, TS)
    update_surface(Case, gm, TS)
    update_forcing(Case, gm, TS)
    update_radiation(Case, gm, TS)
    update_GMV_diagnostics(edmf, gm)
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

    clear_precip_sources(up_thermo)
    microphysics(up_thermo, edmf.UpdVar, edmf.Rain, TS.dt) # causes division error in dry bubble first time step
    microphysics(en_thermo, edmf.EnvVar, edmf.Rain, TS.dt) # saturation adjustment + rain creation
    update_total_precip_sources(up_thermo)
    if edmf.Rain.rain_model == "clima_1m"
        solve_rain_evap(edmf.rainphysics, gm, TS, edmf.Rain.QR)
        # sum updraft and environment rain into bulk rain
        sum_subdomains_rain(edmf.Rain, up_thermo, en_thermo, TS)
    end
    # compute tendencies
    compute_gm_tendencies!(edmf, grid, Case, gm, edmf.ref_state, TS)

    # update
    solve_updraft(edmf, gm, TS)
    if edmf.Rain.rain_model == "clima_1m"
        solve_rain_fall(edmf.rainphysics, gm, TS, edmf.Rain.QR)
    end
    update_cloud_frac(edmf, gm)

    # ----------- TODO: move to compute_tendencies
    implicit_eqs = edmf.implicit_eqs
    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    implicit_eqs.A_θq_gm .= construct_tridiag_diffusion_gm(grid, TS.dt, edmf.rho_ae_KH, ref_state.rho0_half, edmf.ae)
    implicit_eqs.A_uv_gm .= construct_tridiag_diffusion_gm(grid, TS.dt, edmf.rho_ae_KM, ref_state.rho0_half, edmf.ae)
    Δzi = grid.Δzi
    @inbounds for k in real_center_indices(grid)
        implicit_eqs.b_u_gm[k] = GMV.U.values[k]
        implicit_eqs.b_v_gm[k] = GMV.V.values[k]
        implicit_eqs.b_q_tot_gm[k] = edmf.EnvVar.QT.values[k]
        implicit_eqs.b_θ_liq_ice_gm[k] = edmf.EnvVar.H.values[k]
        if is_surface_center(grid, k)
            implicit_eqs.b_u_gm[k] += TS.dt * Case.Sur.rho_uflux * Δzi * ref_state.alpha0_half[k] / edmf.ae[k]
            implicit_eqs.b_v_gm[k] += TS.dt * Case.Sur.rho_vflux * Δzi * ref_state.alpha0_half[k] / edmf.ae[k]
            implicit_eqs.b_q_tot_gm[k] += TS.dt * Case.Sur.rho_qtflux * Δzi * ref_state.alpha0_half[k] / edmf.ae[k]
            implicit_eqs.b_θ_liq_ice_gm[k] += TS.dt * Case.Sur.rho_hflux * Δzi * ref_state.alpha0_half[k] / edmf.ae[k]
        end
    end

    KM = diffusivity_m(edmf).values
    KH = diffusivity_h(edmf).values
    common_args = (
        grid,
        param_set,
        ref_state,
        TS,
        KM,
        KH,
        up.Area.bulkvalues,
        up.Area.values,
        up.W.values,
        en.W.values,
        en.TKE.values,
        edmf.n_updrafts,
        edmf.minimum_area,
        edmf.pressure_plume_spacing,
        edmf.frac_turb_entr,
        edmf.entr_sc,
        edmf.mixing_length,
    )

    implicit_eqs.A_TKE .= construct_tridiag_diffusion_en(common_args..., true)
    implicit_eqs.A_Hvar .= construct_tridiag_diffusion_en(common_args..., false)
    implicit_eqs.A_QTvar .= construct_tridiag_diffusion_en(common_args..., false)
    implicit_eqs.A_HQTcov .= construct_tridiag_diffusion_en(common_args..., false)

    implicit_eqs.b_TKE .= en_diffusion_tendencies(grid, ref_state, TS, up.Area.values, en.TKE)
    implicit_eqs.b_Hvar .= en_diffusion_tendencies(grid, ref_state, TS, up.Area.values, en.Hvar)
    implicit_eqs.b_QTvar .= en_diffusion_tendencies(grid, ref_state, TS, up.Area.values, en.QTvar)
    implicit_eqs.b_HQTcov .= en_diffusion_tendencies(grid, ref_state, TS, up.Area.values, en.HQTcov)
    # -----------

    # Update the grid mean variables with the tendency due to eddy diffusion
    # Solve tri-diagonal systems
    x_q_tot_gm = center_field(grid) # for tridiag solver
    x_θ_liq_ice_gm = center_field(grid) # for tridiag solver
    x_q_tot_gm .= tridiag_solve(implicit_eqs.b_q_tot_gm, implicit_eqs.A_θq_gm)
    x_θ_liq_ice_gm .= tridiag_solve(implicit_eqs.b_θ_liq_ice_gm, implicit_eqs.A_θq_gm)
    GMV.U.new .= tridiag_solve(implicit_eqs.b_u_gm, implicit_eqs.A_uv_gm)
    GMV.V.new .= tridiag_solve(implicit_eqs.b_v_gm, implicit_eqs.A_uv_gm)
    en.TKE.values .= tridiag_solve(implicit_eqs.b_TKE, implicit_eqs.A_TKE)
    en.Hvar.values .= tridiag_solve(implicit_eqs.b_Hvar, implicit_eqs.A_Hvar)
    en.QTvar.values .= tridiag_solve(implicit_eqs.b_QTvar, implicit_eqs.A_QTvar)
    en.HQTcov.values .= tridiag_solve(implicit_eqs.b_HQTcov, implicit_eqs.A_HQTcov)

    @inbounds for k in real_center_indices(grid)
        GMV.QT.new[k] = max(GMV.QT.values[k] + edmf.ae[k] * (x_q_tot_gm[k] - edmf.EnvVar.QT.values[k]), 0.0)
        GMV.H.new[k] = GMV.H.values[k] + edmf.ae[k] * (x_θ_liq_ice_gm[k] - edmf.EnvVar.H.values[k])
        # get the diffusive flux, TODO: move to diagnostics (callbacks?)
        edmf.diffusive_tendency_h[k] = (GMV.H.new[k] - GMV.H.values[k]) * TS.dti
        edmf.diffusive_tendency_qt[k] = (GMV.QT.new[k] - GMV.QT.values[k]) * TS.dti
    end

    # Filter solution, TODO: fuse with `zero_area_fraction_cleanup` and put into `filter_variables!`
    @inbounds for k in real_center_indices(grid)
        en.TKE.values[k] = max(en.TKE.values[k], 0.0)
        en.Hvar.values[k] = max(en.Hvar.values[k], 0.0)
        en.QTvar.values[k] = max(en.QTvar.values[k], 0.0)
        en.HQTcov.values[k] = max(en.HQTcov.values[k], -sqrt(en.Hvar.values[k] * en.QTvar.values[k]))
        en.HQTcov.values[k] = min(en.HQTcov.values[k], sqrt(en.Hvar.values[k] * en.QTvar.values[k]))
    end

    update_GMV_turbulence(edmf, GMV, Case, TS)
    # set values
    set_values_with_new(edmf.UpdVar)
    zero_area_fraction_cleanup(edmf, GMV)

    update(GMV, TS)
    return
end

function update_GMV_turbulence(edmf::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    @inbounds for k in real_center_indices(edmf.grid)
        GMV.H.tendencies[k] += (GMV.H.new[k] - GMV.H.values[k]) * TS.dti
        GMV.QT.tendencies[k] += (GMV.QT.new[k] - GMV.QT.values[k]) * TS.dti
        GMV.U.tendencies[k] += (GMV.U.new[k] - GMV.U.values[k]) * TS.dti
        GMV.V.tendencies[k] += (GMV.V.new[k] - GMV.V.values[k]) * TS.dti
    end

    return
end

function compute_eddy_diffusivities_tke(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables, Case::CasesBase)
    grid = get_grid(edmf)
    param_set = parameter_set(gm)
    g = CPP.grav(param_set)
    ref_state = reference_state(edmf)
    kc_surf = kc_surface(grid)
    ω_pr = CPEDMF.ω_pr(param_set)
    Pr_n = CPEDMF.Pr_n(param_set)
    c_m = CPEDMF.c_m(param_set)
    KM = diffusivity_m(edmf)
    KH = diffusivity_h(edmf)
    en_thermo = edmf.EnvThermo
    en = edmf.EnvVar
    up = edmf.UpdVar
    p0_c = ref_state.p0_half
    ρ0_c = ref_state.rho0_half
    α0_c = ref_state.alpha0_half
    obukhov_length = Case.Sur.obukhov_length

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
        edmf.prandtl_nvec[k] = turbulent_Prandtl_number(obukhov_length, ∇_Ri, Pr_n, ω_pr)

        ml_model = MinDisspLen(;
            z = grid.zc[k].z,
            obukhov_length = obukhov_length,
            tke_surf = en.TKE.values[kc_surf],
            ustar = Case.Sur.ustar,
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

        KM.values[k] = c_m * edmf.mixing_length[k] * sqrt(max(en.TKE.values[k], 0.0))
        KH.values[k] = KM.values[k] / edmf.prandtl_nvec[k]

        en.TKE.buoy[k] = -ml_model.a_en * ρ0_c[k] * KH.values[k] * (bg.∂b∂z_θl + bg.∂b∂z_qt)
    end
    return
end

function set_updraft_surface_bc(edmf::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    grid = get_grid(edmf)
    kc_surf = kc_surface(grid)

    Δzi = grid.Δzi
    zLL = grid.zc[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    alpha0LL = reference_state(edmf).alpha0_half[kc_surf]
    qt_var = get_surface_variance(Case.Sur.rho_qtflux * alpha0LL, Case.Sur.rho_qtflux * alpha0LL, ustar, zLL, oblength)
    h_var = get_surface_variance(Case.Sur.rho_hflux * alpha0LL, Case.Sur.rho_hflux * alpha0LL, ustar, zLL, oblength)

    if Case.Sur.bflux > 0.0
        a_total = edmf.surface_area
        edmf.entr_surface_bc = 2.0 * Δzi
        edmf.detr_surface_bc = 0.0
    else
        # a_total = edmf.surface_area
        a_total = edmf.minimum_area * 0.9
        edmf.entr_surface_bc = 0.0
        edmf.detr_surface_bc = 2.0 * Δzi
    end

    a_ = a_total / edmf.n_updrafts
    @inbounds for i in xrange(edmf.n_updrafts)
        surface_scalar_coeff = percentile_bounds_mean_norm(1.0 - a_total + i * a_, 1.0 - a_total + (i + 1) * a_, 1000)
        edmf.area_surface_bc[i] = a_
        edmf.w_surface_bc[i] = 0.0
        edmf.h_surface_bc[i] = (GMV.H.values[kc_surf] + surface_scalar_coeff * sqrt(h_var))
        edmf.qt_surface_bc[i] = (GMV.QT.values[kc_surf] + surface_scalar_coeff * sqrt(qt_var))
    end
    return
end

function reset_surface_covariance(edmf::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    flux1 = Case.Sur.rho_hflux
    flux2 = Case.Sur.rho_qtflux
    grid = get_grid(edmf)
    kc_surf = kc_surface(grid)
    ref_state = reference_state(edmf)
    zLL = grid.zc[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    alpha0LL = ref_state.alpha0_half[kc_surf]

    gm = GMV
    up = edmf.UpdVar
    en = edmf.EnvVar

    en.TKE.values[kc_surf] = get_surface_tke(Case.Sur.ustar, grid.zc[kc_surf], Case.Sur.obukhov_length)
    get_GMV_CoVar(edmf, up.Area, up.W, up.W, en.W, en.W, en.TKE, gm.W.values, gm.W.values, gm.TKE.values)

    en.Hvar.values[kc_surf] = get_surface_variance(flux1 * alpha0LL, flux1 * alpha0LL, ustar, zLL, oblength)
    en.QTvar.values[kc_surf] = get_surface_variance(flux2 * alpha0LL, flux2 * alpha0LL, ustar, zLL, oblength)
    en.HQTcov.values[kc_surf] = get_surface_variance(flux1 * alpha0LL, flux2 * alpha0LL, ustar, zLL, oblength)
    return
end

# Note: this assumes all variables are defined on half levels not full levels (i.e. phi, psi are not w)
# if covar_e.name is not "tke".
function get_GMV_CoVar(
    edmf::EDMF_PrognosticTKE,
    au::UpdraftVariable,
    ϕ_up::UpdraftVariable,
    ψ_up::UpdraftVariable,
    ϕ_en::EnvironmentVariable,
    ψ_en::EnvironmentVariable,
    covar_e::EnvironmentVariable_2m,
    ϕ_gm,
    ψ_gm,
    gmv_covar,
)

    grid = get_grid(edmf)
    ae = 1 .- au.bulkvalues
    is_tke = covar_e.name == "tke"
    tke_factor = is_tke ? 0.5 : 1

    if is_tke
        @inbounds for k in real_center_indices(grid)
            ϕ_en_dual = dual_faces(ϕ_en.values, grid, k)
            ϕ_gm_dual = dual_faces(ϕ_gm, grid, k)
            ψ_en_dual = dual_faces(ψ_en.values, grid, k)
            ψ_gm_dual = dual_faces(ψ_gm, grid, k)
            Δϕ_dual = ϕ_en_dual .- ϕ_gm_dual
            Δψ_dual = ψ_en_dual .- ψ_gm_dual
            Δϕ = interpf2c(Δϕ_dual, grid, k)
            Δψ = interpf2c(Δψ_dual, grid, k)

            gmv_covar[k] = tke_factor * ae[k] * Δϕ * Δψ + ae[k] * covar_e.values[k]
            @inbounds for i in xrange(edmf.n_updrafts)
                ϕ_up_dual = dual_faces(ϕ_up.values, grid, k, i)
                ϕ_gm_dual = dual_faces(ϕ_gm, grid, k)
                ψ_up_dual = dual_faces(ψ_up.values, grid, k, i)
                ψ_gm_dual = dual_faces(ψ_gm, grid, k)
                Δϕ_dual = ϕ_up_dual .- ϕ_gm_dual
                Δψ_dual = ψ_up_dual .- ψ_gm_dual
                Δϕ = interpf2c(Δϕ_dual, grid, k)
                Δψ = interpf2c(Δψ_dual, grid, k)
                gmv_covar[k] += tke_factor * au.values[i, k] * Δϕ * Δψ
            end
        end
    else

        @inbounds for k in real_center_indices(grid)
            Δϕ = ϕ_en.values[k] - ϕ_gm[k]
            Δψ = ψ_en.values[k] - ψ_gm[k]

            gmv_covar[k] = tke_factor * ae[k] * Δϕ * Δψ + ae[k] * covar_e.values[k]
            @inbounds for i in xrange(edmf.n_updrafts)
                Δϕ = ϕ_up.values[i, k] - ϕ_gm[k]
                Δψ = ψ_up.values[i, k] - ψ_gm[k]
                gmv_covar[k] += tke_factor * au.values[i, k] * Δϕ * Δψ
            end
        end
    end
    return
end

function compute_updraft_closures(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables, Case::CasesBase)
    grid = get_grid(edmf)
    FT = eltype(grid)
    param_set = parameter_set(edmf)
    ref_state = reference_state(edmf)
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar

    upd_cloud_diagnostics(up, ref_state)

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
    return
end

function compute_pressure_plume_spacing(edmf::EDMF_PrognosticTKE)

    param_set = parameter_set(edmf)
    H_up_min = CPEDMF.H_up_min(param_set)
    @inbounds for i in xrange(edmf.n_updrafts)
        edmf.pressure_plume_spacing[i] =
            max(edmf.aspect_ratio * edmf.UpdVar.updraft_top[i], H_up_min * edmf.aspect_ratio)
    end
    return
end

function zero_area_fraction_cleanup(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables)

    up = edmf.UpdVar
    en = edmf.EnvVar
    grid = get_grid(edmf)
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in xrange(up.n_updrafts)
            if up.Area.values[i, k] < edmf.minimum_area
                up.Area.values[i, k] = 0.0
                up.B.values[i, k] = gm.B.values[k]
                up.H.values[i, k] = gm.H.values[k]
                up.QT.values[i, k] = gm.QT.values[k]
                up.T.values[i, k] = gm.T.values[k]
                up.QL.values[i, k] = gm.QL.values[k]
            end
        end

        if sum(up.Area.values[:, k]) == 0.0
            en.B.values[k] = gm.B.values[k]
            en.H.values[k] = gm.H.values[k]
            en.QT.values[k] = gm.QT.values[k]
            en.T.values[k] = gm.T.values[k]
            en.QL.values[k] = gm.QL.values[k]
        end
    end

    @inbounds for k in real_face_indices(grid)
        @inbounds for i in xrange(up.n_updrafts)
            a_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetValue(0))
            a_up_f_i = interpc2f(up.Area.values, grid, k, i; a_bcs...)
            if a_up_f_i < edmf.minimum_area
                up.W.values[i, k] = gm.W.values[k]
            end
        end

        a_up_f = map(1:(up.n_updrafts)) do i
            a_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetValue(0))
            interpc2f(up.Area.values, grid, k, i; a_bcs...)
        end

        if sum(a_up_f) == 0.0
            en.W.values[k] = gm.W.values[k]
        end
    end

    return
end

function solve_updraft(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables, TS::TimeStepping)
    grid = get_grid(edmf)
    param_set = parameter_set(gm)
    ref_state = reference_state(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dti_ = 1.0 / TS.dt
    Δt = TS.dt

    up = edmf.UpdVar
    up_thermo = edmf.UpdThermo
    en = edmf.EnvVar
    a_up = up.Area.values
    a_up_new = up.Area.new
    w_up = up.W.values
    w_up_new = up.W.new
    ρ_0_c = ref_state.rho0_half
    ρ_0_f = ref_state.rho0
    au_lim = edmf.max_area

    @inbounds for i in xrange(up.n_updrafts)
        edmf.entr_sc[i, kc_surf] = edmf.entr_surface_bc
        edmf.detr_sc[i, kc_surf] = edmf.detr_surface_bc
        w_up_new[i, kf_surf] = edmf.w_surface_bc[i]
        a_up_new[i, kc_surf] = edmf.area_surface_bc[i]
    end

    # Solve for updraft area fraction
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in xrange(up.n_updrafts)
            is_surface_center(grid, k) && continue
            w_up_c = interpf2c(w_up, grid, k, i)
            adv = upwind_advection_area(ρ_0_c, a_up[i, :], w_up[i, :], grid, k)

            a_up_c = a_up[i, k]
            entr_term = a_up_c * w_up_c * (edmf.entr_sc[i, k])
            detr_term = a_up_c * w_up_c * (-edmf.detr_sc[i, k])
            a_tendencies = adv + entr_term + detr_term

            a_up_candidate = max(a_up_c + Δt * a_tendencies, 0)

            if a_up_candidate > au_lim
                a_up_candidate = au_lim
                a_up_div = a_up_c > 0.0 ? a_up_c : au_lim
                edmf.detr_sc[i, k] = (((au_lim - a_up_c) * dti_ - adv - entr_term) / (-a_up_div * w_up_c))
            end
            a_up_new[i, k] = a_up_candidate
        end
    end

    entr_w_c = edmf.entr_sc .+ edmf.frac_turb_entr
    detr_w_c = edmf.detr_sc .+ edmf.frac_turb_entr


    # Solve for updraft velocity
    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in xrange(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            anew_k = interpc2f(a_up_new, grid, k, i; a_up_bcs...)

            if anew_k >= edmf.minimum_area
                a_k = interpc2f(a_up, grid, k, i; a_up_bcs...)
                # We know that, since W = 0 at z = 0, these BCs should
                # not matter in the end:
                entr_w = interpc2f(entr_w_c, grid, k, i; bottom = SetValue(0), top = SetValue(0))
                detr_w = interpc2f(detr_w_c, grid, k, i; bottom = SetValue(0), top = SetValue(0))
                B_k = interpc2f(up.B.values, grid, k, i; bottom = SetValue(0), top = SetValue(0))

                adv = upwind_advection_velocity(ρ_0_f, a_up[i, :], w_up[i, :], grid, k; a_up_bcs)
                exch = (ρ_0_f[k] * a_k * w_up[i, k] * (entr_w * en.W.values[k] - detr_w * w_up[i, k]))
                buoy = ρ_0_f[k] * a_k * B_k
                w_tendencies = -adv + exch + buoy + edmf.nh_pressure[i, k]
                w_up_new[i, k] = (ρ_0_f[k] * a_k * w_up[i, k] + Δt * w_tendencies) / (ρ_0_f[k] * anew_k)

                w_up_new[i, k] = max(w_up_new[i, k], 0)
                # TODO: remove a_up_new from this loop.
                if w_up_new[i, k] <= 0.0
                    if !(k > size(a_up_new, 2))
                        a_up_new[i, k] = 0
                    end
                end
            else
                w_up_new[i, k] = 0
                if !(k > size(a_up_new, 2))
                    a_up_new[i, k] = 0
                end
            end
        end
    end

    # Solve for θ_liq_ice & q_tot
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in xrange(up.n_updrafts)
            if is_surface_center(grid, k)
                # at the surface
                if a_up_new[i, k] >= edmf.minimum_area
                    up.H.new[i, k] = edmf.h_surface_bc[i]
                    up.QT.new[i, k] = edmf.qt_surface_bc[i]
                else
                    up.H.new[i, k] = gm.H.values[k]
                    up.QT.new[i, k] = gm.QT.values[k]
                end
                continue
            end

            # write the discrete equations in form
            # c1 * phi_new[k] = c2 * phi[k] + c3 * phi[k-1] + c4 * ϕ_enntr
            if a_up_new[i, k] >= edmf.minimum_area
                w_up_c = interpf2c(w_up, grid, k, i)
                m_k = (ρ_0_c[k] * a_up[i, k] * w_up_c)

                adv = upwind_advection_scalar(ρ_0_c, a_up[i, :], w_up[i, :], up.H.values[i, :], grid, k)
                entr = entr_w_c[i, k] * en.H.values[k]
                detr = detr_w_c[i, k] * up.H.values[i, k]
                θ_liq_ice_rain = ρ_0_c[k] * up_thermo.prec_source_h[i, k]
                θ_liq_ice_tendencies = -adv + m_k * (entr - detr) + θ_liq_ice_rain
                up.H.new[i, k] =
                    (ρ_0_c[k] * a_up[i, k] * up.H.values[i, k] + Δt * θ_liq_ice_tendencies) /
                    (ρ_0_c[k] * a_up_new[i, k])

                adv = upwind_advection_scalar(ρ_0_c, a_up[i, :], w_up[i, :], up.QT.values[i, :], grid, k)
                entr = entr_w_c[i, k] * en.QT.values[k]
                detr = detr_w_c[i, k] * up.QT.values[i, k]
                q_tot_rain = ρ_0_c[k] * up_thermo.prec_source_qt[i, k]
                q_tot_tendencies = -adv + m_k * (entr - detr) + q_tot_rain
                up.QT.new[i, k] = max(
                    (ρ_0_c[k] * a_up[i, k] * up.QT.values[i, k] + Δt * q_tot_tendencies) / (ρ_0_c[k] * a_up_new[i, k]),
                    0.0,
                )

            else
                up.H.new[i, k] = gm.H.values[k]
                up.QT.new[i, k] = gm.QT.values[k]
            end
        end
    end

    return
end

function compute_tke_pressure(edmf::EDMF_PrognosticTKE)
    grid = get_grid(edmf)

    @inbounds for k in real_center_indices(grid)
        edmf.EnvVar.TKE.press[k] = 0.0
        @inbounds for i in xrange(edmf.n_updrafts)
            w_up_c = interpf2c(edmf.UpdVar.W.values, grid, k, i)
            w_en_c = interpf2c(edmf.EnvVar.W.values, grid, k)
            press_c = interpf2c(edmf.nh_pressure, grid, k, i)
            edmf.EnvVar.TKE.press[k] += (w_en_c - w_up_c) * press_c
        end
    end
    return
end

function update_GMV_diagnostics(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables)

    grid = get_grid(edmf)
    en = edmf.EnvVar
    up = edmf.UpdVar
    a_up_bulk = up.Area.bulkvalues
    @inbounds for k in real_center_indices(grid)
        gm.QL.values[k] = (a_up_bulk[k] * up.QL.bulkvalues[k] + (1 - a_up_bulk[k]) * en.QL.values[k])
        gm.T.values[k] = (a_up_bulk[k] * up.T.bulkvalues[k] + (1 - a_up_bulk[k]) * en.T.values[k])
        gm.B.values[k] = (a_up_bulk[k] * up.B.bulkvalues[k] + (1 - a_up_bulk[k]) * en.B.values[k])
    end

    return
end

function initialize_covariance(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables, Case::CasesBase)

    ws = edmf.wstar
    us = Case.Sur.ustar
    zs = edmf.zi
    grid = get_grid(edmf)
    kc_surf = kc_surface(grid)
    en = edmf.EnvVar

    reset_surface_covariance(edmf, gm, Case)

    # grid-mean tke closure over all but z:
    tke_gm(z) = ws * 1.3 * cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) * sqrt(max(1 - z / zs, 0))


    if ws > 0.0
        @inbounds for k in real_center_indices(grid)
            z = grid.zc[k]
            gm.TKE.values[k] = tke_gm(z)
        end
    end
    # TKE initialization from Beare et al, 2006
    if Case.casename == "GABLS"
        @inbounds for k in real_center_indices(grid)
            z = grid.zc[k]
            if (z <= 250.0)
                gm.TKE.values[k] = 0.4 * (1.0 - z / 250.0) * (1.0 - z / 250.0) * (1.0 - z / 250.0)
            end
        end

    elseif Case.casename == "Bomex"
        @inbounds for k in real_center_indices(grid)
            z = grid.zc[k]
            if (z <= 2500.0)
                gm.TKE.values[k] = 1.0 - z / 3000.0
            end
        end

    elseif Case.casename == "Soares" || Case.casename == "Nieuwstadt"
        @inbounds for k in real_center_indices(grid)
            z = grid.zc[k]
            if (z <= 1600.0)
                gm.TKE.values[k] = 0.1 * 1.46 * 1.46 * (1.0 - z / 1600.0)
            end
        end
    end

    if ws > 0.0
        @inbounds for k in real_center_indices(grid)
            z = grid.zc[k]
            # need to rethink of how to initilize the covarinace profiles - for now took the TKE profile
            gm.Hvar.values[k] = gm.Hvar.values[kc_surf] * tke_gm(z)
            gm.QTvar.values[k] = gm.QTvar.values[kc_surf] * tke_gm(z)
            gm.HQTcov.values[k] = gm.HQTcov.values[kc_surf] * tke_gm(z)
        end
    else
        @inbounds for k in real_center_indices(grid)
            gm.Hvar.values[k] = 0.0
            gm.QTvar.values[k] = 0.0
            gm.HQTcov.values[k] = 0.0
        end
    end

    # TKE initialization from Beare et al, 2006
    if Case.casename == "GABLS"
        @inbounds for k in real_center_indices(grid)
            z = grid.zc[k]
            if (z <= 250.0)
                gm.Hvar.values[k] = 0.4 * (1.0 - z / 250.0) * (1.0 - z / 250.0) * (1.0 - z / 250.0)
            end
            gm.QTvar.values[k] = 0.0
            gm.HQTcov.values[k] = 0.0
        end
    end

    @inbounds for k in real_center_indices(grid)
        en.TKE.values[k] = gm.TKE.values[k]
        en.Hvar.values[k] = gm.Hvar.values[k]
        en.QTvar.values[k] = gm.QTvar.values[k]
        en.HQTcov.values[k] = gm.HQTcov.values[k]
    end
    return
end

function compute_covariance_shear(
    edmf::EDMF_PrognosticTKE,
    gm::GridMeanVariables,
    Covar::EnvironmentVariable_2m,
    EnvVar1,
    EnvVar2,
)

    grid = get_grid(edmf)
    ae = 1 .- edmf.UpdVar.Area.bulkvalues
    KH = diffusivity_h(edmf).values
    rho0_half = reference_state(edmf).rho0_half
    is_tke = Covar.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    k_eddy = is_tke ? diffusivity_m(edmf).values : diffusivity_h(edmf).values

    if is_tke
        @inbounds for k in real_center_indices(grid)
            v_cut = ccut(gm.V.values, grid, k)
            ∇v = c∇(v_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            u_cut = ccut(gm.U.values, grid, k)
            ∇u = c∇(u_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_dual = dual_faces(EnvVar2, grid, k)
            var1_dual = dual_faces(EnvVar1, grid, k)

            ∇var2 = ∇f2c(var2_dual, grid, k; bottom = SetValue(0), top = SetGradient(0))
            ∇var1 = ∇f2c(var1_dual, grid, k; bottom = SetValue(0), top = SetGradient(0))

            Covar.shear[k] = tke_factor * 2 * (rho0_half[k] * ae[k] * k_eddy[k] * (∇var1 * ∇var2 + ∇u^2 + ∇v^2))
        end
    else
        @inbounds for k in real_center_indices(grid)
            # Defined correctly only for covariance between half-level variables.
            var1_cut = ccut(EnvVar1, grid, k)
            ∇var1 = c∇(var1_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_cut = ccut(EnvVar2, grid, k)
            ∇var2 = c∇(var2_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            Covar.shear[k] = tke_factor * 2 * (rho0_half[k] * ae[k] * k_eddy[k] * (∇var1 * ∇var2))
        end
    end
    return
end

function compute_covariance_interdomain_src(
    edmf::EDMF_PrognosticTKE,
    au::UpdraftVariable,
    ϕ_up::UpdraftVariable,
    ψ_up::UpdraftVariable,
    ϕ_en::EnvironmentVariable,
    ψ_en::EnvironmentVariable,
    Covar::EnvironmentVariable_2m,
)

    is_tke = Covar.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    grid = get_grid(edmf)
    if is_tke
        @inbounds for k in real_center_indices(grid)
            Covar.interdomain[k] = 0.0
            @inbounds for i in xrange(edmf.n_updrafts)
                Δϕ = interpf2c(ϕ_up.values, grid, k, i) - interpf2c(ϕ_en.values, grid, k)
                Δψ = interpf2c(ψ_up.values, grid, k, i) - interpf2c(ψ_en.values, grid, k)

                Covar.interdomain[k] += tke_factor * au.values[i, k] * (1.0 - au.values[i, k]) * Δϕ * Δψ
            end
        end
    else
        @inbounds for k in real_center_indices(grid)
            Covar.interdomain[k] = 0.0
            @inbounds for i in xrange(edmf.n_updrafts)
                Δϕ = ϕ_up.values[i, k] - ϕ_en.values[k]
                Δψ = ψ_up.values[i, k] - ψ_en.values[k]
                Covar.interdomain[k] += tke_factor * au.values[i, k] * (1.0 - au.values[i, k]) * Δϕ * Δψ
            end
        end
    end
    return
end

function compute_covariance_entr(
    edmf::EDMF_PrognosticTKE,
    Covar::EnvironmentVariable_2m,
    UpdVar1::UpdraftVariable,
    UpdVar2::UpdraftVariable,
    EnvVar1::EnvironmentVariable,
    EnvVar2::EnvironmentVariable,
    GmvVar1::VariablePrognostic,
    GmvVar2::VariablePrognostic,
)

    grid = get_grid(edmf)
    ref_state = reference_state(edmf)
    rho0_half = ref_state.rho0_half
    is_tke = Covar.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    a_up = edmf.UpdVar.Area.values

    @inbounds for k in real_center_indices(grid)
        Covar.entr_gain[k] = 0.0
        Covar.detr_loss[k] = 0.0
        @inbounds for i in xrange(edmf.n_updrafts)
            if edmf.UpdVar.Area.values[i, k] > edmf.minimum_area
                R_up = edmf.pressure_plume_spacing[i]
                updvar1 = is_tke ? interpf2c(UpdVar1.values, grid, k, i) : UpdVar1.values[i, k]
                updvar2 = is_tke ? interpf2c(UpdVar2.values, grid, k, i) : UpdVar2.values[i, k]
                envvar1 = is_tke ? interpf2c(EnvVar1.values, grid, k) : EnvVar1.values[k]
                envvar2 = is_tke ? interpf2c(EnvVar2.values, grid, k) : EnvVar2.values[k]
                gmvvar1 = is_tke ? interpf2c(GmvVar1.values, grid, k) : GmvVar1.values[k]
                gmvvar2 = is_tke ? interpf2c(GmvVar2.values, grid, k) : GmvVar2.values[k]

                eps_turb = edmf.frac_turb_entr[i, k]

                w_u = interpf2c(edmf.UpdVar.W.values, grid, k, i)
                dynamic_entr =
                    tke_factor *
                    rho0_half[k] *
                    a_up[i, k] *
                    abs(w_u) *
                    edmf.detr_sc[i, k] *
                    (updvar1 - envvar1) *
                    (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    rho0_half[k] *
                    a_up[i, k] *
                    abs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                Covar.entr_gain[k] += dynamic_entr + turbulent_entr
                Covar.detr_loss[k] +=
                    tke_factor *
                    rho0_half[k] *
                    a_up[i, k] *
                    abs(w_u) *
                    (edmf.entr_sc[i, k] + eps_turb) *
                    Covar.values[k]
            end
        end
    end

    return
end

function compute_covariance_detr(edmf::EDMF_PrognosticTKE, Covar::EnvironmentVariable_2m)
    grid = get_grid(edmf)
    up = edmf.UpdVar
    ρ0_c = reference_state(edmf).rho0_half
    @inbounds for k in real_center_indices(grid)
        Covar.detr_loss[k] = 0.0
        @inbounds for i in xrange(up.n_updrafts)
            w_up_c = interpf2c(up.W.values, grid, k, i)
            Covar.detr_loss[k] += up.Area.values[i, k] * abs(w_up_c) * edmf.entr_sc[i, k]
        end
        Covar.detr_loss[k] *= ρ0_c[k] * Covar.values[k]
    end
    return
end

function compute_covariance_rain(edmf::EDMF_PrognosticTKE, TS::TimeStepping)
    # TODO defined again in compute_covariance_shear and compute_covaraince
    grid = get_grid(edmf)
    up = edmf.UpdVar
    en = edmf.EnvVar
    en_thermo = edmf.EnvThermo
    ae = 1 .- up.Area.bulkvalues # area of environment
    ρ0_c = reference_state(edmf).rho0_half

    @inbounds for k in real_center_indices(grid)
        en.TKE.rain_src[k] = 0.0
        en.Hvar.rain_src[k] = ρ0_c[k] * ae[k] * 2.0 * en_thermo.Hvar_rain_dt[k]
        en.QTvar.rain_src[k] = ρ0_c[k] * ae[k] * 2.0 * en_thermo.QTvar_rain_dt[k]
        en.HQTcov.rain_src[k] = ρ0_c[k] * ae[k] * en_thermo.HQTcov_rain_dt[k]
    end

    return
end

function compute_covariance_dissipation(edmf::EDMF_PrognosticTKE, Covar::EnvironmentVariable_2m)
    grid = get_grid(edmf)
    param_set = parameter_set(edmf)
    en = edmf.EnvVar
    c_d = CPEDMF.c_d(param_set)
    ae = 1 .- edmf.UpdVar.Area.bulkvalues
    rho0_half = reference_state(edmf).rho0_half

    @inbounds for k in real_center_indices(grid)
        Covar.dissipation[k] = (
            rho0_half[k] * ae[k] * Covar.values[k] * max(en.TKE.values[k], 0)^0.5 / max(edmf.mixing_length[k], 1.0e-3) * c_d
        )
    end
    return
end

function en_diffusion_tendencies(grid::Grid, ref_state::ReferenceState, TS, a_up, covar)
    dti = TS.dti
    b = center_field(grid)
    ρ0_c = ref_state.rho0_half

    ae = 1 .- up_sum(a_up)

    kc_surf = kc_surface(grid)
    covar_surf = covar.values[kc_surf]

    @inbounds for k in real_center_indices(grid)
        if is_surface_center(grid, k)
            b[k] = covar_surf
        else
            b[k] = (
                ρ0_c[k] * ae[k] * covar.values[k] * dti +
                covar.press[k] +
                covar.buoy[k] +
                covar.shear[k] +
                covar.entr_gain[k] +
                covar.rain_src[k]
            )
        end
    end

    return b
end

function GMV_third_m(
    edmf::EDMF_PrognosticTKE,
    Gmv_third_m::VariableDiagnostic,
    env_covar::EnvironmentVariable_2m,
    env_mean::EnvironmentVariable,
    upd_mean::UpdraftVariable,
)

    grid = get_grid(edmf)
    up = edmf.UpdVar
    en = edmf.EnvVar
    ae = 1 .- up.Area.bulkvalues
    au = up.Area.values
    is_tke = env_covar.name == "tke"

    @inbounds for k in real_center_indices(grid)
        mean_en = is_tke ? interpf2c(env_mean.values, grid, k) : env_mean.values[k]
        GMVv_ = ae[k] * mean_en
        @inbounds for i in xrange(up.n_updrafts)
            mean_up = is_tke ? interpf2c(upd_mean.values, grid, k, i) : upd_mean.values[i, k]
            GMVv_ += au[i, k] * mean_up
        end

        # TODO: report bug: i used outside of scope.
        # This is only valid (assuming correct) for 1
        # updraft.
        i_last = last(xrange(up.n_updrafts))
        if is_tke
            w_bcs = (; bottom = SetValue(0), top = SetValue(0))
            w_en_dual = dual_faces(en.W.values, grid, k)
            ∇w_en = ∇f2c(w_en_dual, grid, k; w_bcs...)
            Envcov_ = -edmf.horiz_K_eddy[i_last, k] * ∇w_en
        else
            Envcov_ = env_covar.values[k]
        end

        Upd_cubed = 0.0
        GMVcov_ = ae[k] * (Envcov_ + (mean_en - GMVv_)^2)
        @inbounds for i in xrange(up.n_updrafts)
            mean_up = is_tke ? interpf2c(upd_mean.values, grid, k, i) : upd_mean.values[i, k]
            GMVcov_ += au[i, k] * (mean_up - GMVv_)^2
            Upd_cubed += au[i, k] * mean_up^3
        end

        if is_surface_center(grid, k)
            Gmv_third_m.values[k] = 0.0 # this is here as first value is biased with BC area fraction
        else
            Gmv_third_m.values[k] =
                Upd_cubed + ae[k] * (mean_en^3 + 3 * mean_en * Envcov_) - GMVv_^3 - 3 * GMVcov_ * GMVv_
        end
    end
    return
end


# Update the diagnosis of the inversion height, using the maximum temperature gradient method
function update_inversion(edmf::EDMF_PrognosticTKE, gm::GridMeanVariables, option)
    grid = edmf.grid
    θ_ρ = center_field(grid)
    ∇θ_liq_max = 0.0
    kc_surf = kc_surface(grid)
    param_set = parameter_set(gm)
    ref_state = edmf.ref_state
    p0_c = ref_state.p0_half

    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], gm.H.values[k], gm.QT.values[k])
        θ_ρ[k] = TD.virtual_pottemp(ts)
    end


    if option == "theta_rho"
        @inbounds for k in real_center_indices(grid)
            if θ_ρ[k] > θ_ρ[kc_surf]
                edmf.zi = grid.zc[k]
                break
            end
        end
    elseif option == "thetal_maxgrad"

        @inbounds for k in real_center_indices(grid)
            ∇θ_liq_cut = ccut_downwind(gm.H.values, grid, k)
            ∇θ_liq = c∇_downwind(∇θ_liq_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
            if ∇θ_liq > ∇θ_liq_max
                ∇θ_liq_max = ∇θ_liq
                edmf.zi = grid.zc[k]
            end
        end
    elseif option == "critical_Ri"
        edmf.zi = get_inversion(param_set, θ_ρ, gm.U.values, gm.V.values, grid, Ri_bulk_crit(edmf))

    else
        error("INVERSION HEIGHT OPTION NOT RECOGNIZED")
    end

    return
end
