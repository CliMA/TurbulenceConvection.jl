
#= These methods are to be overloaded by Cases.jl =#
function update_surface end
function update_forcing end
function update_radiation end

function update_cloud_frac(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables)
    # update grid-mean cloud fraction and cloud cover
    aux_tc = center_aux_turbconv(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    a_up_bulk = aux_tc.bulk.area
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        aux_gm.cloud_fraction[k] =
            aux_en.area[k] * aux_en.cloud_fraction[k] + a_up_bulk[k] * aux_tc.bulk.cloud_fraction[k]
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
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    ρ0_f = face_ref_state(state).ρ0
    p0_c = center_ref_state(state).p0
    α0_c = center_ref_state(state).α0
    kf_surf = kf_surface(grid)
    kc_surf = kc_surface(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_tc = center_aux_turbconv(state)

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
            H_horz_adv = Case.Fo.dTdt_hadv[k] / Π
            H_nudge = Case.Fo.dTdt_nudge[k] / Π
            H_fluc = Case.Fo.dTdt_fluc[k] / Π

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
            aux_bulk.qt_tendency_precip_formation[k] +
            aux_en.qt_tendency_precip_formation[k] +
            aux_tc.qt_tendency_precip_sinks[k]
        tendencies_gm.θ_liq_ice[k] +=
            aux_bulk.θ_liq_ice_tendency_precip_formation[k] +
            aux_en.θ_liq_ice_tendency_precip_formation[k] +
            aux_tc.θ_liq_ice_tendency_precip_sinks[k]
    end

    aux_up_f = face_aux_updrafts(state)
    edmf.massflux_h .= 0.0
    edmf.massflux_qt .= 0.0
    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in 1:(up.n_updrafts)
        edmf.m[i, kf_surf] = 0.0
        a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
        @inbounds for k in real_face_indices(grid)
            a_up = interpc2f(aux_up[i].area, grid, k; a_up_bcs...)
            a_en = interpc2f(aux_en.area, grid, k; a_up_bcs...)
            edmf.m[i, k] = ρ0_f[k] * a_up * a_en * (aux_up_f[i].w[k] - aux_en_f.w[k])
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
            h_up_f = interpc2f(aux_up[i].θ_liq_ice, grid, k; m_bcs...)
            qt_up_f = interpc2f(aux_up[i].q_tot, grid, k; m_bcs...)
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
    aeKMu_bc = Case.Sur.rho_uflux / aux_en.area[kc_surf]
    aeKMv_bc = Case.Sur.rho_vflux / aux_en.area[kc_surf]
    @inbounds for k in real_center_indices(grid)
        aeKH_q_tot_cut = dual_faces(edmf.diffusive_flux_qt, grid, k)
        ∇aeKH_q_tot = ∇f2c(aeKH_q_tot_cut, grid, k; bottom = SetValue(aeKHq_tot_bc), top = SetValue(0))
        tendencies_gm.q_tot[k] += -α0_c[k] * aux_en.area[k] * ∇aeKH_q_tot

        aeKH_θ_liq_ice_cut = dual_faces(edmf.diffusive_flux_h, grid, k)
        ∇aeKH_θ_liq_ice = ∇f2c(aeKH_θ_liq_ice_cut, grid, k; bottom = SetValue(aeKHθ_liq_ice_bc), top = SetValue(0))
        tendencies_gm.θ_liq_ice[k] += -α0_c[k] * aux_en.area[k] * ∇aeKH_θ_liq_ice

        aeKM_u_cut = dual_faces(edmf.diffusive_flux_u, grid, k)
        ∇aeKM_u = ∇f2c(aeKM_u_cut, grid, k; bottom = SetValue(aeKMu_bc), top = SetValue(0))
        tendencies_gm.u[k] += -α0_c[k] * aux_en.area[k] * ∇aeKM_u

        aeKM_v_cut = dual_faces(edmf.diffusive_flux_v, grid, k)
        ∇aeKM_v = ∇f2c(aeKM_v_cut, grid, k; bottom = SetValue(aeKMv_bc), top = SetValue(0))
        tendencies_gm.v[k] += -α0_c[k] * aux_en.area[k] * ∇aeKM_v
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
    aux_tc = center_aux_turbconv(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_en = center_aux_environment(state)
    aux_en.area .= 1 .- aux_tc.bulk.area # area of environment
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
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
    aeKMu_bc = -Case.Sur.rho_uflux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KM[kf_surf]
    aeKMv_bc = -Case.Sur.rho_vflux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KM[kf_surf]

    @inbounds for k in real_face_indices(grid)
        q_dual = dual_centers(aux_en.q_tot, grid, k)
        ∇q_tot_f = ∇c2f(q_dual, grid, k; bottom = SetGradient(aeKHq_tot_bc), top = SetGradient(0))
        edmf.diffusive_flux_qt[k] = -aux_tc_f.ρ_ae_KH[k] * ∇q_tot_f

        θ_liq_ice_dual = dual_centers(aux_en.θ_liq_ice, grid, k)
        ∇θ_liq_ice_f = ∇c2f(θ_liq_ice_dual, grid, k; bottom = SetGradient(aeKHθ_liq_ice_bc), top = SetGradient(0))
        edmf.diffusive_flux_h[k] = -aux_tc_f.ρ_ae_KH[k] * ∇θ_liq_ice_f

        u_dual = dual_centers(prog_gm.u, grid, k)
        ∇u_f = ∇c2f(u_dual, grid, k; bottom = SetGradient(aeKMu_bc), top = SetGradient(0))
        edmf.diffusive_flux_u[k] = -aux_tc_f.ρ_ae_KM[k] * ∇u_f

        v_dual = dual_centers(prog_gm.v, grid, k)
        ∇v_f = ∇c2f(v_dual, grid, k; bottom = SetGradient(aeKMv_bc), top = SetGradient(0))
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
    en_thermo = edmf.EnvThermo
    n_updrafts = up.n_updrafts
    prog_gm = center_prog_grid_mean(state)
    prog_en = center_prog_environment(state)

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.

    compute_updraft_surface_bc(edmf, grid, state, Case)
    update_aux!(edmf, gm, grid, state, Case, param_set, TS)

    tendencies_gm = center_tendencies_grid_mean(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    tendencies_en = center_tendencies_environment(state)
    tendencies_pr = center_tendencies_precipitation(state)
    parent(tendencies_gm) .= 0
    parent(tendencies_up) .= 0
    parent(tendencies_up_f) .= 0
    parent(tendencies_en) .= 0
    parent(tendencies_pr) .= 0
    compute_precipitation_formation_tendencies(grid, state, edmf.UpdVar, edmf.Precip, TS.dt, param_set) # causes division error in dry bubble first time step
    microphysics(en_thermo, grid, state, en, edmf.Precip, TS.dt, param_set) # saturation adjustment + rain creation
    if edmf.Precip.precipitation_model == "clima_1m"
        compute_precipitation_sink_tendencies(edmf.PrecipPhys, grid, state, gm, TS)
        compute_precipitation_advection_tendencies(edmf.PrecipPhys, grid, state, gm, TS)
    end

    # compute tendencies
    compute_gm_tendencies!(edmf, grid, state, Case, gm, TS)
    compute_updraft_tendencies(edmf, grid, state, gm)

    compute_en_tendencies!(edmf, grid, state, param_set, TS, :tke, :ρatke, n_updrafts)
    compute_en_tendencies!(edmf, grid, state, param_set, TS, :Hvar, :ρaHvar, n_updrafts)
    compute_en_tendencies!(edmf, grid, state, param_set, TS, :QTvar, :ρaQTvar, n_updrafts)
    compute_en_tendencies!(edmf, grid, state, param_set, TS, :HQTcov, :ρaHQTcov, n_updrafts)

    ###
    ### update (to be removed)
    ###

    Δt = TS.dt
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    @inbounds for k in real_center_indices(grid)
        prog_en.ρatke[k] += TS.dt * tendencies_en.ρatke[k]
        prog_en.ρaHvar[k] += TS.dt * tendencies_en.ρaHvar[k]
        prog_en.ρaQTvar[k] += TS.dt * tendencies_en.ρaQTvar[k]
        prog_en.ρaHQTcov[k] += TS.dt * tendencies_en.ρaHQTcov[k]
        prog_gm.u[k] += tendencies_gm.u[k] * TS.dt
        prog_gm.v[k] += tendencies_gm.v[k] * TS.dt
        prog_gm.θ_liq_ice[k] += tendencies_gm.θ_liq_ice[k] * TS.dt
        prog_gm.q_tot[k] += tendencies_gm.q_tot[k] * TS.dt
        @inbounds for i in 1:(up.n_updrafts)
            prog_up[i].ρarea[k] += Δt * tendencies_up[i].ρarea[k]
            prog_up[i].ρaθ_liq_ice[k] += Δt * tendencies_up[i].ρaθ_liq_ice[k]
            prog_up[i].ρaq_tot[k] += Δt * tendencies_up[i].ρaq_tot[k]
        end
        if edmf.Precip.precipitation_model == "clima_1m"
            prog_pr.q_rai[k] += tendencies_pr.q_rai[k] * TS.dt
            prog_pr.q_sno[k] += tendencies_pr.q_sno[k] * TS.dt
        end
    end

    @inbounds for k in real_face_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            prog_up_f[i].ρaw[k] += Δt * tendencies_up_f[i].ρaw[k]
        end
    end

    ###
    ### Filters
    ###
    set_edmf_surface_bc(edmf, grid, state, up, Case.Sur)
    filter_updraft_vars(edmf, grid, state, gm)
    @inbounds for k in real_center_indices(grid)
        prog_en.ρatke[k] = max(prog_en.ρatke[k], 0.0)
        prog_en.ρaHvar[k] = max(prog_en.ρaHvar[k], 0.0)
        prog_en.ρaQTvar[k] = max(prog_en.ρaQTvar[k], 0.0)
        prog_en.ρaHQTcov[k] = max(prog_en.ρaHQTcov[k], -sqrt(prog_en.ρaHvar[k] * prog_en.ρaQTvar[k]))
        prog_en.ρaHQTcov[k] = min(prog_en.ρaHQTcov[k], sqrt(prog_en.ρaHvar[k] * prog_en.ρaQTvar[k]))
    end

    return
end

function set_edmf_surface_bc(edmf::EDMF_PrognosticTKE, grid, state, up, surface)
    en = edmf.EnvVar
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0
    prog_up = center_prog_updrafts(state)
    prog_en = center_prog_environment(state)
    prog_up_f = face_prog_updrafts(state)
    aux_tc = center_aux_turbconv(state)
    @inbounds for i in 1:(up.n_updrafts)
        prog_up[i].ρarea[kc_surf] = ρ0_c[kc_surf] * edmf.area_surface_bc[i]
        prog_up[i].ρaθ_liq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * edmf.h_surface_bc[i]
        prog_up[i].ρaq_tot[kc_surf] = prog_up[i].ρarea[kc_surf] * edmf.qt_surface_bc[i]
        prog_up_f[i].ρaw[kf_surf] = ρ0_f[kf_surf] * edmf.w_surface_bc[i]
    end

    flux1 = surface.rho_hflux
    flux2 = surface.rho_qtflux
    zLL = grid.zc[kc_surf]
    ustar = surface.ustar
    oblength = surface.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]
    # TODO: is bulk even defined before this is called?
    ae = 1 .- aux_tc.bulk.area # area of environment

    ρ0_ae = ρ0_c[kc_surf] * ae[kc_surf]

    prog_en.ρatke[kc_surf] = ρ0_ae * get_surface_tke(surface.ustar, grid.zc[kc_surf], surface.obukhov_length)
    prog_en.ρaHvar[kc_surf] = ρ0_ae * get_surface_variance(flux1 * α0LL, flux1 * α0LL, ustar, zLL, oblength)
    prog_en.ρaQTvar[kc_surf] = ρ0_ae * get_surface_variance(flux2 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    prog_en.ρaHQTcov[kc_surf] = ρ0_ae * get_surface_variance(flux1 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
end


function compute_updraft_surface_bc(edmf::EDMF_PrognosticTKE, grid, state, Case::CasesBase)
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
        a_ = a_total / edmf.n_updrafts
        @inbounds for i in 1:(edmf.n_updrafts)
            surface_scalar_coeff =
                percentile_bounds_mean_norm(1.0 - a_total + i * a_, 1.0 - a_total + (i + 1) * a_, 1000)
            edmf.area_surface_bc[i] = a_
            edmf.w_surface_bc[i] = 0.0
            edmf.h_surface_bc[i] = prog_gm.θ_liq_ice[kc_surf] + surface_scalar_coeff * sqrt(h_var)
            edmf.qt_surface_bc[i] = prog_gm.q_tot[kc_surf] + surface_scalar_coeff * sqrt(qt_var)
        end
    else
        edmf.entr_surface_bc = 0.0
        edmf.detr_surface_bc = 2.0 * Δzi
        @inbounds for i in 1:(edmf.n_updrafts)
            edmf.area_surface_bc[i] = 0
            edmf.w_surface_bc[i] = 0.0
            edmf.h_surface_bc[i] = prog_gm.θ_liq_ice[kc_surf]
            edmf.qt_surface_bc[i] = prog_gm.q_tot[kc_surf]
        end
    end

    return
end

# Note: this assumes all variables are defined on half levels not full levels (i.e. phi, psi are not w)
# if covar_e.name is not "tke".
function get_GMV_CoVar(edmf::EDMF_PrognosticTKE, grid, state, covar_sym::Symbol, ϕ_sym::Symbol, ψ_sym::Symbol = ϕ_sym)

    aux_tc = center_aux_turbconv(state)
    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    aux_up = center_aux_updrafts(state)
    gmv_covar = getproperty(center_aux_grid_mean(state), covar_sym)
    covar_e = getproperty(center_aux_environment(state), covar_sym)
    prog_gm = is_tke ? prog_gm_f : prog_gm_c
    ϕ_gm = getproperty(prog_gm, ϕ_sym)
    ψ_gm = getproperty(prog_gm, ψ_sym)
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)

    if is_tke
        area_en = aux_en_c.area
        Ic = CCO.InterpolateF2C()
        @. gmv_covar = tke_factor * area_en * Ic(ϕ_en - ϕ_gm) * Ic(ψ_en - ψ_gm) + aux_en_c.area * covar_e
        @inbounds for i in 1:(edmf.n_updrafts)
            ϕ_up = getproperty(aux_up_f[i], ϕ_sym)
            ψ_up = getproperty(aux_up_f[i], ψ_sym)
            @. gmv_covar += tke_factor * aux_up[i].area * Ic(ϕ_up - ϕ_gm) * Ic(ψ_up - ψ_gm)
        end
    else

        @. gmv_covar = tke_factor * aux_en_c.area * (ϕ_en - ϕ_gm) * (ψ_en - ψ_gm) + aux_en_c.area * covar_e
        @inbounds for i in 1:(edmf.n_updrafts)
            ϕ_up = getproperty(aux_up[i], ϕ_sym)
            ψ_up = getproperty(aux_up[i], ψ_sym)
            @. gmv_covar += tke_factor * aux_up[i].area * (ϕ_up - ϕ_gm) * (ψ_up - ψ_gm)
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

function compute_updraft_tendencies(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)

    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_up = center_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_up_f = face_aux_updrafts(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0
    au_lim = edmf.max_area

    parent(tendencies_up) .= 0
    entr_turb_dyn = edmf.entr_sc .+ edmf.frac_turb_entr
    detr_turb_dyn = edmf.detr_sc .+ edmf.frac_turb_entr

    # Solve for updraft area fraction
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            w_up_c = interpf2c(aux_up_f[i].w, grid, k)
            m_k = (ρ0_c[k] * aux_up[i].area[k] * w_up_c)

            adv = upwind_advection_area(ρ0_c, aux_up[i].area, aux_up_f[i].w, grid, k)
            entr_term = m_k * entr_turb_dyn[i, k]
            detr_term = m_k * detr_turb_dyn[i, k]
            tendencies_up[i].ρarea[k] = (ρ0_c[k] * adv + entr_term - detr_term)

            adv = upwind_advection_scalar(ρ0_c, aux_up[i].area, aux_up_f[i].w, aux_up[i].θ_liq_ice, grid, k)
            entr = entr_term * aux_en.θ_liq_ice[k]
            detr = detr_term * aux_up[i].θ_liq_ice[k]
            rain = ρ0_c[k] * aux_up[i].θ_liq_ice_tendency_precip_formation[k]
            tendencies_up[i].ρaθ_liq_ice[k] = -adv + entr - detr + rain

            adv = upwind_advection_scalar(ρ0_c, aux_up[i].area, aux_up_f[i].w, aux_up[i].q_tot, grid, k)
            entr = entr_term * aux_en.q_tot[k]
            detr = detr_term * aux_up[i].q_tot[k]
            rain = ρ0_c[k] * aux_up[i].qt_tendency_precip_formation[k]
            tendencies_up[i].ρaq_tot[k] = -adv + entr - detr + rain
        end
    end

    # Solve for updraft velocity
    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            a_k = interpc2f(aux_up[i].area, grid, k; a_up_bcs...)
            # We know that, since W = 0 at z = 0, these BCs should
            # not matter in the end:
            entr_w = interpc2f(entr_turb_dyn, grid, k, i; bottom = SetValue(0), top = SetValue(0))
            detr_w = interpc2f(detr_turb_dyn, grid, k, i; bottom = SetValue(0), top = SetValue(0))
            B_k = interpc2f(aux_up[i].buoy, grid, k; bottom = SetValue(0), top = SetValue(0))

            adv = upwind_advection_velocity(ρ0_f, aux_up[i].area, aux_up_f[i].w, grid, k; a_up_bcs)
            exch = (ρ0_f[k] * a_k * aux_up_f[i].w[k] * (entr_w * aux_en_f.w[k] - detr_w * aux_up_f[i].w[k]))
            buoy = ρ0_f[k] * a_k * B_k
            tendencies_up_f[i].ρaw[k] = -adv + exch + buoy + edmf.nh_pressure[i, k]
        end
    end
    return
end

function filter_updraft_vars(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)

    up = edmf.UpdVar
    prog_up = center_prog_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0

    @inbounds for i in 1:(up.n_updrafts)
        prog_up[i].ρarea .= max.(prog_up[i].ρarea, 0)
        prog_up[i].ρaθ_liq_ice .= max.(prog_up[i].ρaθ_liq_ice, 0)
        prog_up[i].ρaq_tot .= max.(prog_up[i].ρaq_tot, 0)
        @inbounds for k in real_center_indices(grid)
            prog_up[i].ρarea[k] = min(prog_up[i].ρarea[k], ρ0_c[k] * edmf.max_area)
        end
    end

    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            prog_up_f[i].ρaw[k] = max(prog_up_f[i].ρaw[k], 0)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            ρa_f = interpc2f(prog_up[i].ρarea, grid, k; a_up_bcs...)
            if ρa_f < ρ0_f[k] * edmf.minimum_area
                prog_up_f[i].ρaw[k] = 0
            end
        end
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            prog_up[i].ρaq_tot[k] = max(prog_up[i].ρaq_tot[k], 0)
            ρaw_up_c = interpf2c(prog_up_f[i].ρaw, grid, k)
            if ρaw_up_c <= 0
                prog_up[i].ρarea[k] = 0
                prog_up[i].ρaθ_liq_ice[k] = 0
                prog_up[i].ρaq_tot[k] = 0
            end
            # this is needed to make sure Rico is unchanged.
            # TODO : look into it further to see why
            # a similar filtering of ρaθ_liq_ice breaks the simulation
            if prog_up[i].ρarea[k] / ρ0_c[k] < edmf.minimum_area
                prog_up[i].ρaq_tot[k] = 0
            end
        end
    end
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

    aux_tc = center_aux_turbconv(state)
    aux_tc = center_aux_turbconv(state)
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

            aux_covar.shear[k] =
                tke_factor * 2 * (ρ0_c[k] * aux_en_c.area[k] * k_eddy[k] * (∇var1 * ∇var2 + ∇u^2 + ∇v^2))
        end
    else
        @inbounds for k in real_center_indices(grid)
            # Defined correctly only for covariance between half-level variables.
            var1_cut = ccut(EnvVar1, grid, k)
            ∇var1 = c∇(var1_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_cut = ccut(EnvVar2, grid, k)
            ∇var2 = c∇(var2_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            aux_covar.shear[k] = tke_factor * 2 * (ρ0_c[k] * aux_en_c.area[k] * k_eddy[k] * (∇var1 * ∇var2))
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
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, Covar_sym)
    prog_up = is_tke ? aux_up_f : aux_up
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

                aux_covar.interdomain[k] += tke_factor * aux_up[i].area[k] * (1.0 - aux_up[i].area[k]) * Δϕ * Δψ
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
                aux_covar.interdomain[k] += tke_factor * aux_up[i].area[k] * (1.0 - aux_up[i].area[k]) * Δϕ * Δψ
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

    ρ0_c = center_ref_state(state).ρ0

    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    prog_gm = is_tke ? prog_gm_f : prog_gm_c
    prog_up = is_tke ? aux_up_f : aux_up
    GmvVar1 = getproperty(prog_gm, var1)
    GmvVar2 = getproperty(prog_gm, var2)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_en_c = center_aux_environment(state)
    covar = getproperty(aux_en_c, covar_sym)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    EnvVar1 = getproperty(aux_en, var1)
    EnvVar2 = getproperty(aux_en, var2)

    @inbounds for k in real_center_indices(grid)
        aux_covar.entr_gain[k] = 0.0
        aux_covar.detr_loss[k] = 0.0
        @inbounds for i in 1:(edmf.n_updrafts)
            a_up = aux_up[i].area[k]
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

                w_u = interpf2c(aux_up_f[i].w, grid, k)
                dynamic_entr =
                    tke_factor *
                    ρ0_c[k] *
                    a_up *
                    abs(w_u) *
                    edmf.detr_sc[i, k] *
                    (updvar1 - envvar1) *
                    (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    ρ0_c[k] *
                    a_up *
                    abs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                aux_covar.entr_gain[k] += dynamic_entr + turbulent_entr
                aux_covar.detr_loss[k] +=
                    tke_factor * ρ0_c[k] * a_up * abs(w_u) * (edmf.entr_sc[i, k] + eps_turb) * covar[k]
            end
        end
    end

    return
end

function compute_covariance_detr(edmf::EDMF_PrognosticTKE, grid, state, covar_sym::Symbol)
    up = edmf.UpdVar
    ρ0_c = center_ref_state(state).ρ0
    aux_up_f = face_aux_updrafts(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_up = center_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    covar = getproperty(aux_en, covar_sym)

    @inbounds for k in real_center_indices(grid)
        aux_covar.detr_loss[k] = 0.0
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(aux_up_f[i].w, grid, k)
            aux_covar.detr_loss[k] += aux_up[i].area[k] * abs(w_up_c) * edmf.entr_sc[i, k]
        end
        aux_covar.detr_loss[k] *= ρ0_c[k] * covar[k]
    end
    return
end

function compute_covariance_dissipation(edmf::EDMF_PrognosticTKE, grid, state, covar_sym::Symbol, param_set)
    en = edmf.EnvVar
    c_d = CPEDMF.c_d(param_set)
    aux_tc = center_aux_turbconv(state)
    ρ0_c = center_ref_state(state).ρ0
    prog_en = center_prog_environment(state)
    aux_en = center_aux_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    covar = getproperty(aux_en, covar_sym)

    @inbounds for k in real_center_indices(grid)
        aux_covar.dissipation[k] = (
            ρ0_c[k] * aux_en.area[k] * covar[k] * max(aux_en.tke[k], 0)^0.5 / max(aux_tc.mixing_length[k], 1.0e-3) *
            c_d
        )
    end
    return
end

function compute_en_tendencies!(edmf, grid::Grid, state, param_set, TS, covar_sym::Symbol, prog_sym::Symbol, n_updrafts)
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    tendencies_en = center_tendencies_environment(state)
    prog_covar = getproperty(tendencies_en, prog_sym)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    covar = getproperty(aux_en, covar_sym)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_up = center_aux_updrafts(state)
    w_en_f = face_aux_environment(state).w
    w_en = center_aux_environment(state).w
    c_d = CPEDMF.c_d(param_set)
    is_tke = covar_sym == :tke


    ρ_ae_K = face_aux_turbconv(state).ρ_ae_K
    ρ_ae_K∇ϕ = face_aux_turbconv(state).ρ_ae_K∇ϕ
    ∇ρ_ae_K∇ϕ = center_aux_turbconv(state).∇ρ_ae_K∇ϕ
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    aux_tc = center_aux_turbconv(state)
    ae = 1 .- aux_tc.bulk.area
    aeK = is_tke ? ae .* KM : ae .* KH
    aeK_bcs = (; bottom = SetValue(aeK[kc_surf]), top = SetValue(aeK[kc_toa]))
    prog_bcs = (; bottom = SetGradient(0), top = SetGradient(0))

    @inbounds for k in real_face_indices(grid)
        ρ_ae_K[k] = interpc2f(aeK, grid, k; aeK_bcs...) * ρ0_f[k]
        ϕ_dual = dual_centers(covar, grid, k)
        ρ_ae_K∇ϕ[k] = ρ_ae_K[k] * ∇c2f(ϕ_dual, grid, k; prog_bcs...)
    end
    @inbounds for k in real_center_indices(grid)
        ρ_ae_K∇ϕ_dual = dual_faces(ρ_ae_K∇ϕ, grid, k)
        ∇ρ_ae_K∇ϕ[k] = ∇f2c(ρ_ae_K∇ϕ_dual, grid, k)
    end

    mixing_length = aux_tc.mixing_length
    minimum_area = edmf.minimum_area
    frac_turb_entr = edmf.frac_turb_entr
    pressure_plume_spacing = edmf.pressure_plume_spacing
    entr_sc = edmf.entr_sc

    ρaew_en_ϕ = center_aux_turbconv(state).ρaew_en_ϕ
    FT = eltype(grid)

    @inbounds for k in real_center_indices(grid)
        w_en[k] = interpf2c(w_en_f, grid, k)
        ρaew_en_ϕ[k] = ρ0_c[k] * aux_en.area[k] * w_en[k] * covar[k]
    end

    kc_surf = kc_surface(grid)
    covar_surf = covar[kc_surf]

    @inbounds for k in real_center_indices(grid)
        D_env = sum(1:n_updrafts) do i
            if aux_up[i].area[k] > minimum_area
                turb_entr = frac_turb_entr[i, k]
                R_up = pressure_plume_spacing[i]
                w_up_c = interpf2c(aux_up_f[i].w, grid, k)
                D_env_i = ρ0_c[k] * aux_up[i].area[k] * w_up_c * (entr_sc[i, k] + turb_entr)
            else
                D_env_i = FT(0)
            end
            D_env_i
        end
        dissipation = ρ0_c[k] * aux_en.area[k] * c_d * sqrt(max(aux_en.tke[k], 0)) / max(mixing_length[k], 1)
        ρaew_en_ϕ_cut = ccut_downwind(ρaew_en_ϕ, grid, k)
        ∇ρaew_en_ϕ = c∇_downwind(ρaew_en_ϕ_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
        if is_surface_center(grid, k)
            prog_covar[k] = covar_surf
        else
            prog_covar[k] = (
                aux_covar.press[k] +
                aux_covar.buoy[k] +
                aux_covar.shear[k] +
                aux_covar.entr_gain[k] +
                aux_covar.rain_src[k] - D_env * covar[k] - dissipation * covar[k] - ∇ρaew_en_ϕ + ∇ρ_ae_K∇ϕ[k]
            )
        end
    end

    return nothing
end

function GMV_third_m(edmf::EDMF_PrognosticTKE, grid, state, covar_en_sym::Symbol, var::Symbol, gm_third_m_sym::Symbol)

    gm_third_m = getproperty(center_aux_grid_mean(state), gm_third_m_sym)

    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_tc = center_aux_turbconv(state)
    aux_up_f = face_aux_updrafts(state)
    is_tke = covar_en_sym == :tke
    covar_en = getproperty(center_aux_environment(state), covar_en_sym)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    aux_up = center_aux_updrafts(state)
    var_en = getproperty(aux_en, var)

    @inbounds for k in real_center_indices(grid)
        mean_en = is_tke ? interpf2c(var_en, grid, k) : var_en[k]
        GMVv_ = aux_en_c.area[k] * mean_en
        @inbounds for i in 1:(up.n_updrafts)
            var_up = is_tke ? getproperty(aux_up_f[i], var) : getproperty(aux_up[i], var)
            mean_up = is_tke ? interpf2c(var_up, grid, k) : var_up[k]
            GMVv_ += aux_up[i].area[k] * mean_up
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
        GMVcov_ = aux_en_c.area[k] * (Envcov_ + (mean_en - GMVv_)^2)
        @inbounds for i in 1:(up.n_updrafts)
            var_up = is_tke ? getproperty(aux_up_f[i], var) : getproperty(aux_up[i], var)
            mean_up = is_tke ? interpf2c(var_up, grid, k) : var_up[k]
            GMVcov_ += aux_up[i].area[k] * (mean_up - GMVv_)^2
            Upd_cubed += aux_up[i].area[k] * mean_up^3
        end

        if is_surface_center(grid, k)
            gm_third_m[k] = 0.0 # this is here as first value is biased with BC area fraction
        else
            gm_third_m[k] =
                Upd_cubed + aux_en_c.area[k] * (mean_en^3 + 3 * mean_en * Envcov_) - GMVv_^3 - 3 * GMVcov_ * GMVv_
        end
    end
    return
end
