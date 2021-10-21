
function initialize(edmf::EDMF_PrognosticTKE, grid, state, Case::CasesBase, gm::GridMeanVariables, TS::TimeStepping)
    initialize_covariance(edmf, grid, state, gm, Case)
    if Case.casename == "DryBubble"
        initialize_DryBubble(edmf, grid, state, edmf.UpdVar, gm)
    else
        initialize(edmf, grid, state, edmf.UpdVar, gm)
    end
    return
end

# Initialize the IO pertaining to this class
function initialize_io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)

    initialize_io(edmf.UpdVar, Stats)
    initialize_io(edmf.EnvVar, Stats)
    initialize_io(edmf.Rain, Stats)

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
    return
end

function io(edmf::EDMF_PrognosticTKE, grid, state, Stats::NetCDFIO_Stats, TS::TimeStepping, param_set)

    mean_nh_pressure = face_field(grid)
    mean_nh_pressure_adv = face_field(grid)
    mean_nh_pressure_drag = face_field(grid)
    mean_nh_pressure_b = face_field(grid)

    mean_asp_ratio = center_field(grid)
    mean_entr_sc = center_field(grid)
    mean_detr_sc = center_field(grid)
    massflux = center_field(grid)
    mean_frac_turb_entr = center_field(grid)
    mean_horiz_K_eddy = center_field(grid)
    mean_sorting_function = center_field(grid)
    mean_b_mix = center_field(grid)

    io(edmf.UpdVar, grid, state, Stats)
    io(edmf.EnvVar, grid, state, Stats)
    io(edmf.Rain, grid, state, Stats, edmf.UpdThermo, edmf.EnvThermo, TS)

    prog_up = center_prog_updrafts(state)
    aux_tc = center_aux_tc(state)
    a_up_bulk = aux_tc.bulk.area

    write_ts(Stats, "rd", StatsBase.mean(edmf.pressure_plume_spacing))

    @inbounds for k in real_center_indices(grid)
        if a_up_bulk[k] > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                massflux[k] += interpf2c(edmf.m, grid, k, i)
                mean_entr_sc[k] += prog_up[i].area[k] * edmf.entr_sc[i, k] / a_up_bulk[k]
                mean_detr_sc[k] += prog_up[i].area[k] * edmf.detr_sc[i, k] / a_up_bulk[k]
                mean_asp_ratio[k] += prog_up[i].area[k] * edmf.asp_ratio[i, k] / a_up_bulk[k]
                mean_frac_turb_entr[k] += prog_up[i].area[k] * edmf.frac_turb_entr[i, k] / a_up_bulk[k]
                mean_horiz_K_eddy[k] += prog_up[i].area[k] * edmf.horiz_K_eddy[i, k] / a_up_bulk[k]
                mean_sorting_function[k] += prog_up[i].area[k] * edmf.sorting_function[i, k] / a_up_bulk[k]
                mean_b_mix[k] += prog_up[i].area[k] * edmf.b_mix[i, k] / a_up_bulk[k]
            end
        end
    end

    @inbounds for k in real_face_indices(grid)
        a_up_bulk_f =
            interpc2f(a_up_bulk, grid, k; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        if a_up_bulk_f > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                a_up_f = interpc2f(
                    prog_up[i].area,
                    grid,
                    k;
                    bottom = SetValue(edmf.area_surface_bc[i]),
                    top = SetZeroGradient(),
                )
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
    write_profile(Stats, "massflux_h", edmf.massflux_h)
    write_profile(Stats, "massflux_qt", edmf.massflux_qt)
    write_profile(Stats, "massflux_tendency_h", edmf.massflux_tendency_h)
    write_profile(Stats, "massflux_tendency_qt", edmf.massflux_tendency_qt)
    write_profile(Stats, "diffusive_flux_h", edmf.diffusive_flux_h)
    write_profile(Stats, "diffusive_flux_qt", edmf.diffusive_flux_qt)
    write_profile(Stats, "diffusive_flux_u", edmf.diffusive_flux_u)
    write_profile(Stats, "diffusive_flux_v", edmf.diffusive_flux_v)
    write_profile(Stats, "diffusive_tendency_h", edmf.diffusive_tendency_h)
    write_profile(Stats, "diffusive_tendency_qt", edmf.diffusive_tendency_qt)
    write_profile(Stats, "total_flux_h", edmf.massflux_h .+ edmf.diffusive_flux_h)
    write_profile(Stats, "total_flux_qt", edmf.massflux_h .+ edmf.diffusive_flux_qt)
    write_profile(Stats, "mixing_length", edmf.mixing_length)
    write_profile(Stats, "updraft_qt_precip", edmf.UpdThermo.qt_tendency_rain_formation_tot)
    write_profile(Stats, "updraft_thetal_precip", edmf.UpdThermo.θ_liq_ice_tendency_rain_formation_tot)

    #Different mixing lengths : Ignacio
    write_profile(Stats, "ed_length_scheme", edmf.mls)
    write_profile(Stats, "mixing_length_ratio", edmf.ml_ratio)
    write_profile(Stats, "entdet_balance_length", edmf.l_entdet)

    compute_covariance_dissipation(edmf, grid, state, :tke, param_set)
    compute_covariance_detr(edmf, grid, state, :tke)

    compute_covariance_dissipation(edmf, grid, state, :Hvar, param_set)
    compute_covariance_dissipation(edmf, grid, state, :QTvar, param_set)
    compute_covariance_dissipation(edmf, grid, state, :HQTcov, param_set)
    compute_covariance_detr(edmf, grid, state, :Hvar)
    compute_covariance_detr(edmf, grid, state, :QTvar)
    compute_covariance_detr(edmf, grid, state, :HQTcov)
    return
end

#= These methods are to be overloaded by Cases.jl =#
function update_surface end
function update_forcing end
function update_radiation end

function update_cloud_frac(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables)
    # update grid-mean cloud fraction and cloud cover
    prog_up = center_prog_updrafts(state)
    aux_tc = center_aux_tc(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    a_up_bulk = aux_tc.bulk.area
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        aux_en.area[k] = 1 - a_up_bulk[k]
        aux_gm.cloud_fraction[k] =
            aux_en.area[k] * aux_en.cloud_fraction[k] + a_up_bulk[k] * edmf.UpdVar.cloud_fraction[k]
    end
    gm.cloud_cover = min(edmf.EnvVar.cloud_cover + sum(edmf.UpdVar.cloud_cover), 1)
end

function compute_gm_tendencies!(edmf::EDMF_PrognosticTKE, grid, state, Case, gm, TS)
    tendencies_gm = center_tendencies_grid_mean(state)
    parent(tendencies_gm.u) .= 0
    parent(tendencies_gm.v) .= 0
    parent(tendencies_gm.q_tot) .= 0
    parent(tendencies_gm.θ_liq_ice) .= 0
    param_set = parameter_set(gm)
    prog_gm = center_prog_grid_mean(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    ρ0_f = face_ref_state(state).ρ0
    p0_c = center_ref_state(state).p0
    α0_c = center_ref_state(state).α0
    kf_surf = kf_surface(grid)
    kc_surf = kc_surface(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area # area of environment

    @inbounds for k in real_center_indices(grid)
        # Apply large-scale horizontal advection tendencies
        ts = thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(ts)

        if Case.Fo.apply_coriolis
            tendencies_gm.u[k] -= Case.Fo.coriolis_param * (Case.Fo.vg[k] - prog_gm.v[k])
            tendencies_gm.v[k] += Case.Fo.coriolis_param * (Case.Fo.ug[k] - prog_gm.u[k])
        end
        if rad_type(Case.Rad) <: Union{RadiationDYCOMS_RF01, RadiationLES}
            tendencies_gm.θ_liq_ice[k] += Case.Rad.dTdt[k] / Π
        end
        H_cut = ccut_downwind(prog_gm.θ_liq_ice, grid, k)
        q_tot_cut = ccut_downwind(prog_gm.q_tot, grid, k)
        ∇H = c∇_downwind(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
        ∇q_tot = c∇_downwind(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))

        if force_type(Case.Fo) <: ForcingDYCOMS_RF01
            tendencies_gm.q_tot[k] += Case.Fo.dqtdt[k]
            # Apply large-scale subsidence tendencies
            tendencies_gm.θ_liq_ice[k] -= ∇H * Case.Fo.subsidence[k]
            tendencies_gm.q_tot[k] -= ∇q_tot * Case.Fo.subsidence[k]
        end

        if force_type(Case.Fo) <: ForcingStandard
            if Case.Fo.apply_subsidence
                tendencies_gm.θ_liq_ice[k] -= ∇H * Case.Fo.subsidence[k]
                tendencies_gm.q_tot[k] -= ∇q_tot * Case.Fo.subsidence[k]
            end
            tendencies_gm.θ_liq_ice[k] += Case.Fo.dTdt[k] / Π
            tendencies_gm.q_tot[k] += Case.Fo.dqtdt[k]
        end

        if force_type(Case.Fo) <: ForcingLES
            H_horz_adv = Case.Fo.dtdt_hadv[k] / Π
            H_nudge = Case.Fo.dtdt_nudge[k] / Π
            H_fluc = Case.Fo.dtdt_fluc[k] / Π

            gm_U_nudge_k = (Case.Fo.u_nudge[k] - prog_gm.u[k]) / Case.Fo.nudge_tau
            gm_V_nudge_k = (Case.Fo.v_nudge[k] - prog_gm.v[k]) / Case.Fo.nudge_tau
            if Case.Fo.apply_subsidence
                # Apply large-scale subsidence tendencies
                gm_H_subsidence_k = -∇H * Case.Fo.subsidence[k]
                gm_QT_subsidence_k = -∇q_tot * Case.Fo.subsidence[k]
            else
                gm_H_subsidence_k = 0.0
                gm_QT_subsidence_k = 0.0
            end

            tendencies_gm.θ_liq_ice[k] += H_horz_adv + H_nudge + H_fluc + gm_H_subsidence_k
            tendencies_gm.q_tot[k] +=
                Case.Fo.dqtdt_hadv[k] + Case.Fo.dqtdt_nudge[k] + gm_QT_subsidence_k + Case.Fo.dqtdt_fluc[k]

            tendencies_gm.u[k] += gm_U_nudge_k
            tendencies_gm.v[k] += gm_V_nudge_k
        end
        tendencies_gm.q_tot[k] +=
            edmf.UpdThermo.qt_tendency_rain_formation_tot[k] +
            edmf.EnvThermo.qt_tendency_rain_formation[k] +
            edmf.RainPhys.qt_tendency_rain_evap[k]
        tendencies_gm.θ_liq_ice[k] +=
            edmf.UpdThermo.θ_liq_ice_tendency_rain_formation_tot[k] +
            edmf.EnvThermo.θ_liq_ice_tendency_rain_formation[k] +
            edmf.RainPhys.θ_liq_ice_tendency_rain_evap[k]
    end

    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    edmf.massflux_h .= 0.0
    edmf.massflux_qt .= 0.0
    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in 1:(up.n_updrafts)
        edmf.m[i, kf_surf] = 0.0
        a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
        @inbounds for k in real_face_indices(grid)
            a_up = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
            a_en = interpc2f(ae, grid, k; a_up_bcs...)
            edmf.m[i, k] = ρ0_f[k] * a_up * a_en * (prog_up_f[i].w[k] - aux_en_f.w[k])
        end
    end

    @inbounds for k in real_face_indices(grid)
        edmf.massflux_h[k] = 0.0
        edmf.massflux_qt[k] = 0.0
        # We know that, since W = 0 at z = 0, m = 0 also, and
        # therefore θ_liq_ice / q_tot values do not matter
        m_bcs = (; bottom = SetValue(0), top = SetValue(0))
        h_en_f = interpc2f(aux_en.θ_liq_ice, grid, k; m_bcs...)
        qt_en_f = interpc2f(aux_en.q_tot, grid, k; m_bcs...)
        @inbounds for i in 1:(up.n_updrafts)
            h_up_f = interpc2f(prog_up[i].θ_liq_ice, grid, k; m_bcs...)
            qt_up_f = interpc2f(prog_up[i].q_tot, grid, k; m_bcs...)
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
        tendencies_gm.θ_liq_ice[k] += mf_tend_h
        tendencies_gm.q_tot[k] += mf_tend_qt
    end

    aeKHq_tot_bc = Case.Sur.rho_qtflux / aux_en.area[kc_surf]
    aeKHθ_liq_ice_bc = Case.Sur.rho_hflux / aux_en.area[kc_surf]
    aeKHu_bc = Case.Sur.rho_uflux / aux_en.area[kc_surf]
    aeKHv_bc = Case.Sur.rho_vflux / aux_en.area[kc_surf]
    @inbounds for k in real_center_indices(grid)
        aeKH_q_tot_cut = dual_faces(edmf.diffusive_flux_qt, grid, k)
        ∇aeKH_q_tot = ∇f2c(aeKH_q_tot_cut, grid, k; bottom = SetValue(aeKHq_tot_bc), top = SetValue(0))
        tendencies_gm.q_tot[k] += -α0_c[k] * ae[k] * ∇aeKH_q_tot

        aeKH_θ_liq_ice_cut = dual_faces(edmf.diffusive_flux_h, grid, k)
        ∇aeKH_θ_liq_ice = ∇f2c(aeKH_θ_liq_ice_cut, grid, k; bottom = SetValue(aeKHθ_liq_ice_bc), top = SetValue(0))
        tendencies_gm.θ_liq_ice[k] += -α0_c[k] * ae[k] * ∇aeKH_θ_liq_ice

        aeKM_u_cut = dual_faces(edmf.diffusive_flux_u, grid, k)
        ∇aeKM_u = ∇f2c(aeKM_u_cut, grid, k; bottom = SetValue(aeKHu_bc), top = SetValue(0))
        tendencies_gm.u[k] += -α0_c[k] * ae[k] * ∇aeKM_u

        aeKM_v_cut = dual_faces(edmf.diffusive_flux_v, grid, k)
        ∇aeKM_v = ∇f2c(aeKM_v_cut, grid, k; bottom = SetValue(aeKHv_bc), top = SetValue(0))
        tendencies_gm.v[k] += -α0_c[k] * ae[k] * ∇aeKM_v
    end
end

function compute_diffusive_fluxes(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    gm::GridMeanVariables,
    Case::CasesBase,
    TS::TimeStepping,
    param_set,
)
    ρ0_f = face_ref_state(state).ρ0
    aux_tc = center_aux_tc(state)
    aux_tc_f = face_aux_tc(state)
    aux_en = center_aux_environment(state)
    aux_en.area .= 1 .- aux_tc.bulk.area # area of environment
    KM = center_aux_tc(state).KM
    KH = center_aux_tc(state).KH
    aeKM = aux_en.area .* KM
    aeKH = aux_en.area .* KH
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_surf = kf_surface(grid)
    prog_gm = center_prog_grid_mean(state)
    aeKM_bcs = (; bottom = SetValue(aeKM[kc_surf]), top = SetValue(aeKM[kc_toa]))
    aeKH_bcs = (; bottom = SetValue(aeKH[kc_surf]), top = SetValue(aeKH[kc_toa]))

    @inbounds for k in real_face_indices(grid)
        aux_tc_f.ρ_ae_KH[k] = interpc2f(aeKH, grid, k; aeKH_bcs...) * ρ0_f[k]
        aux_tc_f.ρ_ae_KM[k] = interpc2f(aeKM, grid, k; aeKM_bcs...) * ρ0_f[k]
    end

    aeKHq_tot_bc = -Case.Sur.rho_qtflux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KH[kf_surf]
    aeKHθ_liq_ice_bc = -Case.Sur.rho_hflux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KH[kf_surf]
    aeKHu_bc = -Case.Sur.rho_uflux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KM[kf_surf]
    aeKHv_bc = -Case.Sur.rho_vflux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KM[kf_surf]

    @inbounds for k in real_face_indices(grid)
        q_dual = dual_centers(aux_en.q_tot, grid, k)
        ∇q_tot_f = ∇c2f(q_dual, grid, k; bottom = SetGradient(aeKHq_tot_bc), top = SetGradient(0))
        edmf.diffusive_flux_qt[k] = -aux_tc_f.ρ_ae_KH[k] * ∇q_tot_f

        θ_liq_ice_dual = dual_centers(aux_en.θ_liq_ice, grid, k)
        ∇θ_liq_ice_f = ∇c2f(θ_liq_ice_dual, grid, k; bottom = SetGradient(aeKHθ_liq_ice_bc), top = SetGradient(0))
        edmf.diffusive_flux_h[k] = -aux_tc_f.ρ_ae_KH[k] * ∇θ_liq_ice_f

        u_dual = dual_centers(prog_gm.u, grid, k)
        ∇u_f = ∇c2f(u_dual, grid, k; bottom = SetGradient(aeKHu_bc), top = SetGradient(0))
        edmf.diffusive_flux_u[k] = -aux_tc_f.ρ_ae_KM[k] * ∇u_f

        v_dual = dual_centers(prog_gm.v, grid, k)
        ∇v_f = ∇c2f(v_dual, grid, k; bottom = SetGradient(aeKHv_bc), top = SetGradient(0))
        edmf.diffusive_flux_v[k] = -aux_tc_f.ρ_ae_KM[k] * ∇v_f
    end
    return
end

# Perform the update of the scheme
function update(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    gm = gm
    up = edmf.UpdVar
    en = edmf.EnvVar
    param_set = parameter_set(gm)
    up_thermo = edmf.UpdThermo
    en_thermo = edmf.EnvThermo
    n_updrafts = up.n_updrafts
    prog_gm = center_prog_grid_mean(state)
    prog_en = center_prog_environment(state)
    prog_up_f = face_prog_updrafts(state)

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.
    set_updraft_surface_bc(edmf, grid, state, gm, Case)
    update_aux!(edmf, gm, grid, state, Case, param_set, TS)

    tendencies_gm = center_tendencies_grid_mean(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_en = center_tendencies_environment(state)
    tendencies_ra = center_tendencies_rain(state)
    parent(tendencies_gm) .= 0
    parent(tendencies_up) .= 0
    parent(tendencies_en) .= 0
    parent(tendencies_ra) .= 0
    compute_rain_formation_tendencies(up_thermo, grid, state, edmf.UpdVar, edmf.Rain, TS.dt, param_set) # causes division error in dry bubble first time step
    microphysics(en_thermo, grid, state, edmf.EnvVar, edmf.Rain, TS.dt, param_set) # saturation adjustment + rain creation
    if edmf.Rain.rain_model == "clima_1m"
        compute_rain_evap_tendencies(edmf.RainPhys, grid, state, gm, TS)
        compute_rain_advection_tendencies(edmf.RainPhys, grid, state, gm, TS)
    end

    # compute tendencies
    compute_gm_tendencies!(edmf, grid, state, Case, gm, TS)
    compute_updraft_tendencies(edmf, grid, state, gm, TS)
    # ----------- TODO: move to compute_tendencies
    implicit_eqs = edmf.implicit_eqs
    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    prog_up = center_prog_updrafts(state)

    common_args = (
        grid,
        param_set,
        state,
        TS,
        up.n_updrafts,
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

    implicit_eqs.b_TKE .= en_diffusion_tendencies(grid, state, TS, :tke, n_updrafts)
    implicit_eqs.b_Hvar .= en_diffusion_tendencies(grid, state, TS, :Hvar, n_updrafts)
    implicit_eqs.b_QTvar .= en_diffusion_tendencies(grid, state, TS, :QTvar, n_updrafts)
    implicit_eqs.b_HQTcov .= en_diffusion_tendencies(grid, state, TS, :HQTcov, n_updrafts)
    # -----------

    ###
    ### update
    ###
    update_updraft(edmf, grid, state, gm, TS)
    updraft_set_values(edmf, grid, state, gm)
    if edmf.Rain.rain_model == "clima_1m"
        update_rain(edmf.Rain, grid, state, up_thermo, en_thermo, edmf.RainPhys, TS)
    end

    parent(prog_en.tke) .= implicit_eqs.A_TKE \ implicit_eqs.b_TKE
    parent(prog_en.Hvar) .= implicit_eqs.A_Hvar \ implicit_eqs.b_Hvar
    parent(prog_en.QTvar) .= implicit_eqs.A_QTvar \ implicit_eqs.b_QTvar
    parent(prog_en.HQTcov) .= implicit_eqs.A_HQTcov \ implicit_eqs.b_HQTcov
    @inbounds for k in real_center_indices(grid)
        prog_gm.u[k] += tendencies_gm.u[k] * TS.dt
        prog_gm.v[k] += tendencies_gm.v[k] * TS.dt
        prog_gm.θ_liq_ice[k] += tendencies_gm.θ_liq_ice[k] * TS.dt
        prog_gm.q_tot[k] += tendencies_gm.q_tot[k] * TS.dt
    end

    ###
    ### set values
    ###
    @inbounds for k in real_center_indices(grid)
        prog_en.tke[k] = max(prog_en.tke[k], 0.0)
        prog_en.Hvar[k] = max(prog_en.Hvar[k], 0.0)
        prog_en.QTvar[k] = max(prog_en.QTvar[k], 0.0)
        prog_en.HQTcov[k] = max(prog_en.HQTcov[k], -sqrt(prog_en.Hvar[k] * prog_en.QTvar[k]))
        prog_en.HQTcov[k] = min(prog_en.HQTcov[k], sqrt(prog_en.Hvar[k] * prog_en.QTvar[k]))
    end
    return
end

function set_updraft_surface_bc(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables, Case::CasesBase)
    kc_surf = kc_surface(grid)

    Δzi = grid.Δzi
    zLL = grid.zc[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]
    prog_gm = center_prog_grid_mean(state)
    qt_var = get_surface_variance(Case.Sur.rho_qtflux * α0LL, Case.Sur.rho_qtflux * α0LL, ustar, zLL, oblength)
    h_var = get_surface_variance(Case.Sur.rho_hflux * α0LL, Case.Sur.rho_hflux * α0LL, ustar, zLL, oblength)

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
    @inbounds for i in 1:(edmf.n_updrafts)
        surface_scalar_coeff = percentile_bounds_mean_norm(1.0 - a_total + i * a_, 1.0 - a_total + (i + 1) * a_, 1000)
        edmf.area_surface_bc[i] = a_
        edmf.w_surface_bc[i] = 0.0
        edmf.h_surface_bc[i] = (prog_gm.θ_liq_ice[kc_surf] + surface_scalar_coeff * sqrt(h_var))
        edmf.qt_surface_bc[i] = (prog_gm.q_tot[kc_surf] + surface_scalar_coeff * sqrt(qt_var))
    end
    return
end

function reset_surface_covariance(edmf::EDMF_PrognosticTKE, grid, state, gm, Case::CasesBase)
    flux1 = Case.Sur.rho_hflux
    flux2 = Case.Sur.rho_qtflux
    kc_surf = kc_surface(grid)
    zLL = grid.zc[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]

    gm = gm
    up = edmf.UpdVar
    en = edmf.EnvVar
    prog_en = center_prog_environment(state)

    prog_en.tke[kc_surf] = get_surface_tke(Case.Sur.ustar, grid.zc[kc_surf], Case.Sur.obukhov_length)
    get_GMV_CoVar(edmf, grid, state, :tke, :w)

    prog_en.Hvar[kc_surf] = get_surface_variance(flux1 * α0LL, flux1 * α0LL, ustar, zLL, oblength)
    prog_en.QTvar[kc_surf] = get_surface_variance(flux2 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    prog_en.HQTcov[kc_surf] = get_surface_variance(flux1 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    return
end

# Note: this assumes all variables are defined on half levels not full levels (i.e. phi, psi are not w)
# if covar_e.name is not "tke".
function get_GMV_CoVar(edmf::EDMF_PrognosticTKE, grid, state, covar_sym::Symbol, ϕ_sym::Symbol, ψ_sym::Symbol = ϕ_sym)

    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area
    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    gmv_covar = getproperty(center_aux_grid_mean(state), covar_sym)
    covar_e = getproperty(center_prog_environment(state), covar_sym)
    prog_gm = is_tke ? prog_gm_f : prog_gm_c
    ϕ_gm = getproperty(prog_gm, ϕ_sym)
    ψ_gm = getproperty(prog_gm, ψ_sym)
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)

    if is_tke
        @inbounds for k in real_center_indices(grid)
            ϕ_en_dual = dual_faces(ϕ_en, grid, k)
            ϕ_gm_dual = dual_faces(ϕ_gm, grid, k)
            ψ_en_dual = dual_faces(ψ_en, grid, k)
            ψ_gm_dual = dual_faces(ψ_gm, grid, k)
            Δϕ_dual = ϕ_en_dual .- ϕ_gm_dual
            Δψ_dual = ψ_en_dual .- ψ_gm_dual
            Δϕ = interpf2c(Δϕ_dual, grid, k)
            Δψ = interpf2c(Δψ_dual, grid, k)

            gmv_covar[k] = tke_factor * ae[k] * Δϕ * Δψ + ae[k] * covar_e[k]
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up_var = getproperty(prog_up_f[i], ϕ_sym)
                ψ_up_var = getproperty(prog_up_f[i], ψ_sym)
                ϕ_up_dual = dual_faces(ϕ_up_var, grid, k)
                ϕ_gm_dual = dual_faces(ϕ_gm, grid, k)
                ψ_up_dual = dual_faces(ψ_up_var, grid, k)
                ψ_gm_dual = dual_faces(ψ_gm, grid, k)
                Δϕ_dual = ϕ_up_dual .- ϕ_gm_dual
                Δψ_dual = ψ_up_dual .- ψ_gm_dual
                Δϕ = interpf2c(Δϕ_dual, grid, k)
                Δψ = interpf2c(Δψ_dual, grid, k)
                gmv_covar[k] += tke_factor * prog_up[i].area[k] * Δϕ * Δψ
            end
        end
    else

        @inbounds for k in real_center_indices(grid)
            Δϕ = ϕ_en[k] - ϕ_gm[k]
            Δψ = ψ_en[k] - ψ_gm[k]

            gmv_covar[k] = tke_factor * ae[k] * Δϕ * Δψ + ae[k] * covar_e[k]
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up_var = getproperty(prog_up[i], ϕ_sym)
                ψ_up_var = getproperty(prog_up[i], ψ_sym)
                Δϕ = ϕ_up_var[k] - ϕ_gm[k]
                Δψ = ψ_up_var[k] - ψ_gm[k]
                gmv_covar[k] += tke_factor * prog_up[i].area[k] * Δϕ * Δψ
            end
        end
    end
    return
end

function compute_pressure_plume_spacing(edmf::EDMF_PrognosticTKE, param_set)

    H_up_min = CPEDMF.H_up_min(param_set)
    @inbounds for i in 1:(edmf.n_updrafts)
        edmf.pressure_plume_spacing[i] =
            max(edmf.aspect_ratio * edmf.UpdVar.updraft_top[i], H_up_min * edmf.aspect_ratio)
    end
    return
end

function compute_updraft_tendencies(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables, TS::TimeStepping)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dti_ = 1.0 / TS.dt
    Δt = TS.dt

    up = edmf.UpdVar
    up_thermo = edmf.UpdThermo
    en = edmf.EnvVar
    prog_up = center_prog_updrafts(state)
    aux_up = center_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    prog_up_f = face_prog_updrafts(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    ρ_0_c = center_ref_state(state).ρ0
    ρ_0_f = face_ref_state(state).ρ0
    au_lim = edmf.max_area

    @inbounds for i in 1:(up.n_updrafts)
        @inbounds for k in real_center_indices(grid)
            # tendencies_up[i].area[k] = 0
            tendencies_up[i].ρarea[k] = 0
            tendencies_up[i].θ_liq_ice[k] = 0
            tendencies_up[i].q_tot[k] = 0
        end
        @inbounds for k in real_face_indices(grid)
            tendencies_up_f[i].w[k] = 0
        end
    end

    @inbounds for i in 1:(up.n_updrafts)
        edmf.entr_sc[i, kc_surf] = edmf.entr_surface_bc
        edmf.detr_sc[i, kc_surf] = edmf.detr_surface_bc
    end

    # Solve for updraft area fraction
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            adv = upwind_advection_area(ρ_0_c, prog_up[i].area, prog_up_f[i].w, grid, k)

            ρa_up_c = ρ_0_c[k] * prog_up[i].area[k]
            entr_term = ρa_up_c * w_up_c * (edmf.entr_sc[i, k])
            detr_term = ρa_up_c * w_up_c * (-edmf.detr_sc[i, k])
            tendencies_up[i].ρarea[k] = (ρ_0_c[k] * adv + entr_term + detr_term)
            ρa_up_candidate = ρa_up_c + Δt * tendencies_up[i].ρarea[k]
            ρau_lim = ρ_0_c[k] * au_lim
            if ρa_up_candidate > ρau_lim
                ρa_up_div = ρa_up_c > 0.0 ? ρa_up_c : ρau_lim
                ρa_up_candidate = ρau_lim
                edmf.detr_sc[i, k] = (((ρau_lim - ρa_up_c) * dti_ - ρ_0_c[k] * adv - entr_term) / (-ρa_up_div * w_up_c))
                tendencies_up[i].ρarea[k] = (ρa_up_candidate - ρa_up_c) * dti_
            end
        end
    end

    entr_w_c = edmf.entr_sc .+ edmf.frac_turb_entr
    detr_w_c = edmf.detr_sc .+ edmf.frac_turb_entr

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            m_k = (ρ_0_c[k] * prog_up[i].area[k] * w_up_c)

            adv = upwind_advection_scalar(ρ_0_c, prog_up[i].area, prog_up_f[i].w, prog_up[i].θ_liq_ice, grid, k)
            entr = entr_w_c[i, k] * aux_en.θ_liq_ice[k]
            detr = detr_w_c[i, k] * prog_up[i].θ_liq_ice[k]
            rain = ρ_0_c[k] * up_thermo.θ_liq_ice_tendency_rain_formation[i, k]
            tendencies_up[i].θ_liq_ice[k] = -adv + m_k * (entr - detr) + rain

            adv = upwind_advection_scalar(ρ_0_c, prog_up[i].area, prog_up_f[i].w, prog_up[i].q_tot, grid, k)
            entr = entr_w_c[i, k] * aux_en.q_tot[k]
            detr = detr_w_c[i, k] * prog_up[i].q_tot[k]
            rain = ρ_0_c[k] * up_thermo.qt_tendency_rain_formation[i, k]
            tendencies_up[i].q_tot[k] = -adv + m_k * (entr - detr) + rain
        end
    end

    # Solve for updraft velocity
    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            a_k = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
            # We know that, since W = 0 at z = 0, these BCs should
            # not matter in the end:
            entr_w = interpc2f(entr_w_c, grid, k, i; bottom = SetValue(0), top = SetValue(0))
            detr_w = interpc2f(detr_w_c, grid, k, i; bottom = SetValue(0), top = SetValue(0))
            B_k = interpc2f(aux_up[i].buoy, grid, k; bottom = SetValue(0), top = SetValue(0))

            adv = upwind_advection_velocity(ρ_0_f, prog_up[i].area, prog_up_f[i].w, grid, k; a_up_bcs)
            exch = (ρ_0_f[k] * a_k * prog_up_f[i].w[k] * (entr_w * aux_en_f.w[k] - detr_w * prog_up_f[i].w[k]))
            buoy = ρ_0_f[k] * a_k * B_k
            tendencies_up_f[i].w[k] = -adv + exch + buoy + edmf.nh_pressure[i, k]
        end
    end
    return
end


function update_updraft(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables, TS::TimeStepping)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dti_ = 1.0 / TS.dt
    Δt = TS.dt

    up = edmf.UpdVar
    up_thermo = edmf.UpdThermo
    en = edmf.EnvVar
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    prog_up = center_prog_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    prog_up_f = face_prog_updrafts(state)
    ρ_0_c = center_ref_state(state).ρ0
    ρ_0_f = face_ref_state(state).ρ0

    @inbounds for i in 1:(up.n_updrafts)
        prog_up_f[i].w[kf_surf] = edmf.w_surface_bc[i]
        prog_up[i].area[kc_surf] = edmf.area_surface_bc[i]
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            ρa_up_c = ρ_0_c[k] * prog_up[i].area[k]
            # a_up_candidate = max(a_up_c + Δt * tendencies_up[i].area[k], 0)
            prog_up[i].ρarea[k] = max(ρa_up_c + Δt * tendencies_up[i].ρarea[k], 0)
        end
    end

    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            a_k = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
            prog_up_f[i].ρaw[k] = ρ_0_f[k] * a_k * prog_up_f[i].w[k] + Δt * tendencies_up_f[i].w[k]
        end
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            if is_surface_center(grid, k)
                continue
            end

            prog_up[i].ρaθ_liq_ice[k] =
                ρ_0_c[k] * prog_up[i].area[k] * prog_up[i].θ_liq_ice[k] + Δt * tendencies_up[i].θ_liq_ice[k]
            prog_up[i].ρaq_tot[k] = ρ_0_c[k] * prog_up[i].area[k] * prog_up[i].q_tot[k] + Δt * tendencies_up[i].q_tot[k]
        end
    end
    return
end

function updraft_set_values(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)

    up = edmf.UpdVar
    prog_up = center_prog_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    prog_up_f = face_prog_updrafts(state)
    ρ_0_c = center_ref_state(state).ρ0
    ρ_0_f = face_ref_state(state).ρ0

    # set values
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            prog_up[i].area[k] = prog_up[i].ρarea[k] / ρ_0_c[k]
        end
    end

    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            anew_k = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
            if anew_k >= edmf.minimum_area
                prog_up_f[i].w[k] = max(prog_up_f[i].ρaw[k] / (ρ_0_f[k] * anew_k), 0)
                if prog_up_f[i].w[k] <= 0.0
                    if !(k.i > length(vec(prog_up[i].area)))
                        prog_up[i].area[Cent(k.i)] = 0
                    end
                end
            else
                prog_up_f[i].w[k] = 0
                if !(k.i > length(vec(prog_up[i].area)))
                    prog_up[i].area[Cent(k.i)] = 0
                end
            end
        end
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            if is_surface_center(grid, k)
                # at the surface
                if prog_up[i].area[k] >= edmf.minimum_area
                    prog_up[i].θ_liq_ice[k] = edmf.h_surface_bc[i]
                    prog_up[i].q_tot[k] = edmf.qt_surface_bc[i]
                else
                    prog_up[i].θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
                    prog_up[i].q_tot[k] = prog_gm.q_tot[k]
                end
                continue
            end
            if prog_up[i].area[k] >= edmf.minimum_area
                prog_up[i].θ_liq_ice[k] = prog_up[i].ρaθ_liq_ice[k] / (ρ_0_c[k] * prog_up[i].area[k])
                prog_up[i].q_tot[k] = max(prog_up[i].ρaq_tot[k] / (ρ_0_c[k] * prog_up[i].area[k]), 0.0)
            else
                prog_up[i].θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
                prog_up[i].q_tot[k] = prog_gm.q_tot[k]
            end
        end
    end
    return
end

function initialize_covariance(edmf::EDMF_PrognosticTKE, grid, state, gm, Case::CasesBase)

    kc_surf = kc_surface(grid)
    en = edmf.EnvVar
    aux_gm = center_aux_grid_mean(state)
    prog_en = center_prog_environment(state)

    prog_en.tke .= aux_gm.tke

    reset_surface_covariance(edmf, grid, state, gm, Case)
    aux_gm.Hvar .= aux_gm.Hvar[kc_surf] .* aux_gm.tke
    aux_gm.QTvar .= aux_gm.QTvar[kc_surf] .* aux_gm.tke
    aux_gm.HQTcov .= aux_gm.HQTcov[kc_surf] .* aux_gm.tke

    prog_en.Hvar .= aux_gm.Hvar
    prog_en.QTvar .= aux_gm.QTvar
    prog_en.HQTcov .= aux_gm.HQTcov
    return
end

function compute_covariance_shear(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    gm::GridMeanVariables,
    covar_sym::Symbol,
    en_var1_sym::Symbol,
    en_var2_sym::Symbol = en_var1_sym,
)

    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area # area of environment
    aux_tc = center_aux_tc(state)
    ρ0_c = center_ref_state(state).ρ0
    prog_gm = center_prog_grid_mean(state)
    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    k_eddy = is_tke ? aux_tc.KM : aux_tc.KH
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)

    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    EnvVar1 = getproperty(aux_en, en_var1_sym)
    EnvVar2 = getproperty(aux_en, en_var2_sym)

    if is_tke
        @inbounds for k in real_center_indices(grid)
            v_cut = ccut(prog_gm.v, grid, k)
            ∇v = c∇(v_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            u_cut = ccut(prog_gm.u, grid, k)
            ∇u = c∇(u_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_dual = dual_faces(EnvVar2, grid, k)
            var1_dual = dual_faces(EnvVar1, grid, k)

            ∇var2 = ∇f2c(var2_dual, grid, k; bottom = SetValue(0), top = SetGradient(0))
            ∇var1 = ∇f2c(var1_dual, grid, k; bottom = SetValue(0), top = SetGradient(0))

            aux_covar.shear[k] = tke_factor * 2 * (ρ0_c[k] * ae[k] * k_eddy[k] * (∇var1 * ∇var2 + ∇u^2 + ∇v^2))
        end
    else
        @inbounds for k in real_center_indices(grid)
            # Defined correctly only for covariance between half-level variables.
            var1_cut = ccut(EnvVar1, grid, k)
            ∇var1 = c∇(var1_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_cut = ccut(EnvVar2, grid, k)
            ∇var2 = c∇(var2_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            aux_covar.shear[k] = tke_factor * 2 * (ρ0_c[k] * ae[k] * k_eddy[k] * (∇var1 * ∇var2))
        end
    end
    return
end

function compute_covariance_interdomain_src(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    Covar_sym::Symbol,
    ϕ_var::Symbol,
    ψ_var::Symbol = ϕ_var,
)

    is_tke = Covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    prog_up_c = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, Covar_sym)
    prog_up = is_tke ? prog_up_f : prog_up_c
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    ϕ_en = getproperty(aux_en, ϕ_var)
    ψ_en = getproperty(aux_en, ψ_var)

    if is_tke
        @inbounds for k in real_center_indices(grid)
            aux_covar.interdomain[k] = 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up = getproperty(prog_up[i], ϕ_var)
                ψ_up = getproperty(prog_up[i], ψ_var)
                Δϕ = interpf2c(ϕ_up, grid, k) - interpf2c(ϕ_en, grid, k)
                Δψ = interpf2c(ψ_up, grid, k) - interpf2c(ψ_en, grid, k)

                aux_covar.interdomain[k] += tke_factor * prog_up_c[i].area[k] * (1.0 - prog_up_c[i].area[k]) * Δϕ * Δψ
            end
        end
    else
        @inbounds for k in real_center_indices(grid)
            aux_covar.interdomain[k] = 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up = getproperty(prog_up[i], ϕ_var)
                ψ_up = getproperty(prog_up[i], ψ_var)
                Δϕ = ϕ_up[k] - ϕ_en[k]
                Δψ = ψ_up[k] - ψ_en[k]
                aux_covar.interdomain[k] += tke_factor * prog_up_c[i].area[k] * (1.0 - prog_up_c[i].area[k]) * Δϕ * Δψ
            end
        end
    end
    return
end

function compute_covariance_entr(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    covar_sym::Symbol,
    var1::Symbol,
    var2::Symbol = var1,
)

    ρ_0_c = center_ref_state(state).ρ0

    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    prog_up_c = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    prog_gm = is_tke ? prog_gm_f : prog_gm_c
    prog_up = is_tke ? prog_up_f : prog_up_c
    GmvVar1 = getproperty(prog_gm, var1)
    GmvVar2 = getproperty(prog_gm, var2)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    prog_en = center_prog_environment(state)
    prog_covar = getproperty(prog_en, covar_sym)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    EnvVar1 = getproperty(aux_en, var1)
    EnvVar2 = getproperty(aux_en, var2)

    @inbounds for k in real_center_indices(grid)
        aux_covar.entr_gain[k] = 0.0
        aux_covar.detr_loss[k] = 0.0
        @inbounds for i in 1:(edmf.n_updrafts)
            a_up = prog_up_c[i].area[k]
            if a_up > edmf.minimum_area
                R_up = edmf.pressure_plume_spacing[i]
                up_var_1 = getproperty(prog_up[i], var1)
                up_var_2 = getproperty(prog_up[i], var2)
                updvar1 = is_tke ? interpf2c(up_var_1, grid, k) : up_var_1[k]
                updvar2 = is_tke ? interpf2c(up_var_2, grid, k) : up_var_2[k]
                envvar1 = is_tke ? interpf2c(EnvVar1, grid, k) : EnvVar1[k]
                envvar2 = is_tke ? interpf2c(EnvVar2, grid, k) : EnvVar2[k]
                gmvvar1 = is_tke ? interpf2c(GmvVar1, grid, k) : GmvVar1[k]
                gmvvar2 = is_tke ? interpf2c(GmvVar2, grid, k) : GmvVar2[k]

                eps_turb = edmf.frac_turb_entr[i, k]

                w_u = interpf2c(prog_up_f[i].w, grid, k)
                dynamic_entr =
                    tke_factor *
                    ρ_0_c[k] *
                    a_up *
                    abs(w_u) *
                    edmf.detr_sc[i, k] *
                    (updvar1 - envvar1) *
                    (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    ρ_0_c[k] *
                    a_up *
                    abs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                aux_covar.entr_gain[k] += dynamic_entr + turbulent_entr
                aux_covar.detr_loss[k] +=
                    tke_factor * ρ_0_c[k] * a_up * abs(w_u) * (edmf.entr_sc[i, k] + eps_turb) * prog_covar[k]
            end
        end
    end

    return
end

function compute_covariance_detr(edmf::EDMF_PrognosticTKE, grid, state, covar_sym::Symbol)
    up = edmf.UpdVar
    ρ0_c = center_ref_state(state).ρ0
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)

    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    prog_en = center_prog_environment(state)
    prog_covar = getproperty(prog_en, covar_sym)

    @inbounds for k in real_center_indices(grid)
        aux_covar.detr_loss[k] = 0.0
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            aux_covar.detr_loss[k] += prog_up[i].area[k] * abs(w_up_c) * edmf.entr_sc[i, k]
        end
        aux_covar.detr_loss[k] *= ρ0_c[k] * prog_covar[k]
    end
    return
end

function compute_covariance_dissipation(edmf::EDMF_PrognosticTKE, grid, state, covar_sym::Symbol, param_set)
    en = edmf.EnvVar
    c_d = CPEDMF.c_d(param_set)
    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area
    ρ0_c = center_ref_state(state).ρ0
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    prog_covar = getproperty(prog_en, covar_sym)
    prog_en = center_prog_environment(state)

    @inbounds for k in real_center_indices(grid)
        aux_covar.dissipation[k] =
            (ρ0_c[k] * ae[k] * prog_covar[k] * max(prog_en.tke[k], 0)^0.5 / max(edmf.mixing_length[k], 1.0e-3) * c_d)
    end
    return
end

function en_diffusion_tendencies(grid::Grid, state, TS, covar_sym::Symbol, n_updrafts)
    dti = TS.dti
    b = center_field(grid)
    ρ0_c = center_ref_state(state).ρ0
    prog_up = center_prog_updrafts(state)
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    prog_covar = getproperty(prog_en, covar_sym)
    aux_covar = getproperty(aux_en_2m, covar_sym)

    ae = center_field(grid)

    @inbounds for k in real_center_indices(grid)
        ae[k] = 1 .- sum(ntuple(i -> prog_up[i].area[k], n_updrafts))
    end

    kc_surf = kc_surface(grid)
    covar_surf = prog_covar[kc_surf]

    @inbounds for k in real_center_indices(grid)
        if is_surface_center(grid, k)
            b[k] = covar_surf
        else
            b[k] = (
                ρ0_c[k] * ae[k] * prog_covar[k] * dti +
                aux_covar.press[k] +
                aux_covar.buoy[k] +
                aux_covar.shear[k] +
                aux_covar.entr_gain[k] +
                aux_covar.rain_src[k]
            )
        end
    end

    return b
end

function GMV_third_m(edmf::EDMF_PrognosticTKE, grid, state, covar_en_sym::Symbol, var::Symbol, gm_third_m_sym::Symbol)

    gm_third_m = getproperty(center_aux_grid_mean(state), gm_third_m_sym)

    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    is_tke = covar_en_sym == :tke
    covar_en = getproperty(center_prog_environment(state), covar_en_sym)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    var_en = getproperty(aux_en, var)

    @inbounds for k in real_center_indices(grid)
        mean_en = is_tke ? interpf2c(var_en, grid, k) : var_en[k]
        GMVv_ = ae[k] * mean_en
        @inbounds for i in 1:(up.n_updrafts)
            var_up = is_tke ? getproperty(prog_up_f[i], var) : getproperty(prog_up[i], var)
            mean_up = is_tke ? interpf2c(var_up, grid, k) : var_up[k]
            GMVv_ += prog_up[i].area[k] * mean_up
        end

        # TODO: report bug: i used outside of scope.
        # This is only valid (assuming correct) for 1
        # updraft.
        i_last = last(1:(up.n_updrafts))
        if is_tke
            w_bcs = (; bottom = SetValue(0), top = SetValue(0))
            w_en_dual = dual_faces(aux_en_f.w, grid, k)
            ∇w_en = ∇f2c(w_en_dual, grid, k; w_bcs...)
            Envcov_ = -edmf.horiz_K_eddy[i_last, k] * ∇w_en
        else
            Envcov_ = covar_en[k]
        end

        Upd_cubed = 0.0
        GMVcov_ = ae[k] * (Envcov_ + (mean_en - GMVv_)^2)
        @inbounds for i in 1:(up.n_updrafts)
            var_up = is_tke ? getproperty(prog_up_f[i], var) : getproperty(prog_up[i], var)
            mean_up = is_tke ? interpf2c(var_up, grid, k) : var_up[k]
            GMVcov_ += prog_up[i].area[k] * (mean_up - GMVv_)^2
            Upd_cubed += prog_up[i].area[k] * mean_up^3
        end

        if is_surface_center(grid, k)
            gm_third_m[k] = 0.0 # this is here as first value is biased with BC area fraction
        else
            gm_third_m[k] = Upd_cubed + ae[k] * (mean_en^3 + 3 * mean_en * Envcov_) - GMVv_^3 - 3 * GMVcov_ * GMVv_
        end
    end
    return
end
