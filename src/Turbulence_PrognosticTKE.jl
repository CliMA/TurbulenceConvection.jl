
#= These methods are to be overloaded by Cases.jl =#
function update_surface end
function update_forcing end
function update_radiation end

function update_cloud_frac(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables)
    # update grid-mean cloud fraction and cloud cover
    aux_bulk = center_aux_bulk(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    a_up_bulk = aux_bulk.area
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        aux_gm.cloud_fraction[k] = aux_en.area[k] * aux_en.cloud_fraction[k] + a_up_bulk[k] * aux_bulk.cloud_fraction[k]
    end
    gm.cloud_cover = min(edmf.EnvVar.cloud_cover + sum(edmf.UpdVar.cloud_cover), 1)
end

function compute_les_Γᵣ(z::ClimaCore.Geometry.ZPoint, τᵣ::Real = 24.0 * 3600.0, zᵢ::Real = 3000.0, zᵣ::Real = 3500.0)
    # returns height-dependent relaxation timescale from eqn. 9 in `Shen et al. 2021`
    if z < zᵢ
        return 0.0
    elseif zᵢ <= z <= zᵣ
        cos_arg = pi * ((z - zᵢ) / (zᵣ - zᵢ))
        return (0.5 / τᵣ) * (1 - cos(cos_arg))
    elseif z > zᵣ
        return (1 / τᵣ)
    end
end

function compute_gm_tendencies!(edmf::EDMF_PrognosticTKE, grid, state, Case, gm, TS)
    tendencies_gm = center_tendencies_grid_mean(state)
    param_set = parameter_set(gm)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_up_f = face_aux_updrafts(state)
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
            tendencies_gm.u[k] -= Case.Fo.coriolis_param * (aux_gm.vg[k] - prog_gm.v[k])
            tendencies_gm.v[k] += Case.Fo.coriolis_param * (aux_gm.ug[k] - prog_gm.u[k])
        end
        if rad_type(Case.Rad) <: Union{RadiationDYCOMS_RF01, RadiationLES}
            tendencies_gm.θ_liq_ice[k] += aux_gm.dTdt_rad[k] / Π
        end
        H_cut = ccut_downwind(prog_gm.θ_liq_ice, grid, k)
        q_tot_cut = ccut_downwind(prog_gm.q_tot, grid, k)
        ∇H = c∇_downwind(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
        ∇q_tot = c∇_downwind(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))

        if force_type(Case.Fo) <: ForcingDYCOMS_RF01
            tendencies_gm.q_tot[k] += aux_gm.dqtdt[k]
            # Apply large-scale subsidence tendencies
            tendencies_gm.θ_liq_ice[k] -= ∇H * aux_gm.subsidence[k]
            tendencies_gm.q_tot[k] -= ∇q_tot * aux_gm.subsidence[k]
        end

        if force_type(Case.Fo) <: ForcingStandard
            if Case.Fo.apply_subsidence
                tendencies_gm.θ_liq_ice[k] -= ∇H * aux_gm.subsidence[k]
                tendencies_gm.q_tot[k] -= ∇q_tot * aux_gm.subsidence[k]
            end
            tendencies_gm.θ_liq_ice[k] += aux_gm.dTdt[k] / Π
            tendencies_gm.q_tot[k] += aux_gm.dqtdt[k]
        end

        if force_type(Case.Fo) <: ForcingLES
            H_horz_adv = aux_gm.dTdt_hadv[k] / Π
            H_fluc = aux_gm.dTdt_fluc[k] / Π

            gm_U_nudge_k = (aux_gm.u_nudge[k] - prog_gm.u[k]) / Case.Fo.nudge_tau
            gm_V_nudge_k = (aux_gm.v_nudge[k] - prog_gm.v[k]) / Case.Fo.nudge_tau

            Γᵣ = compute_les_Γᵣ(grid.zc[k])
            if Γᵣ != 0
                tau_k = 1.0 / Γᵣ
                gm_H_nudge_k = (aux_gm.H_nudge[k] - prog_gm.θ_liq_ice[k]) / tau_k
                gm_q_tot_nudge_k = (aux_gm.qt_nudge[k] - prog_gm.q_tot[k]) / tau_k
            else
                gm_H_nudge_k = 0.0
                gm_q_tot_nudge_k = 0.0
            end

            if Case.Fo.apply_subsidence
                # Apply large-scale subsidence tendencies
                gm_H_subsidence_k = -∇H * aux_gm.subsidence[k]
                gm_QT_subsidence_k = -∇q_tot * aux_gm.subsidence[k]
            else
                gm_H_subsidence_k = 0.0
                gm_QT_subsidence_k = 0.0
            end

            tendencies_gm.θ_liq_ice[k] += H_horz_adv + gm_H_nudge_k + H_fluc + gm_H_subsidence_k
            tendencies_gm.q_tot[k] +=
                aux_gm.dqtdt_hadv[k] + gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k] + gm_QT_subsidence_k

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

    # TODO: we shouldn't need to call parent here
    parent(aux_tc_f.massflux_h) .= 0
    parent(aux_tc_f.massflux_qt) .= 0
    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in 1:(up.n_updrafts)
        aux_up_f[i].massflux[kf_surf] = 0.0
        a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
        @inbounds for k in real_face_indices(grid)
            a_up = interpc2f(aux_up[i].area, grid, k; a_up_bcs...)
            a_en = interpc2f(aux_en.area, grid, k; a_up_bcs...)
            aux_up_f[i].massflux[k] = ρ0_f[k] * a_up * a_en * (aux_up_f[i].w[k] - aux_en_f.w[k])
        end
    end

    @inbounds for k in real_face_indices(grid)
        aux_tc_f.massflux_h[k] = 0.0
        aux_tc_f.massflux_qt[k] = 0.0
        # We know that, since W = 0 at z = 0, m = 0 also, and
        # therefore θ_liq_ice / q_tot values do not matter
        m_bcs = (; bottom = SetValue(0), top = SetValue(0))
        h_en_f = interpc2f(aux_en.θ_liq_ice, grid, k; m_bcs...)
        qt_en_f = interpc2f(aux_en.q_tot, grid, k; m_bcs...)
        @inbounds for i in 1:(up.n_updrafts)
            h_up_f = interpc2f(aux_up[i].θ_liq_ice, grid, k; m_bcs...)
            qt_up_f = interpc2f(aux_up[i].q_tot, grid, k; m_bcs...)
            aux_tc_f.massflux_h[k] += aux_up_f[i].massflux[k] * (h_up_f - h_en_f)
            aux_tc_f.massflux_qt[k] += aux_up_f[i].massflux[k] * (qt_up_f - qt_en_f)
        end
    end

    # Compute the  mass flux tendencies
    # Adjust the values of the grid mean variables
    @inbounds for k in real_center_indices(grid)
        mf_tend_h_dual = dual_faces(aux_tc_f.massflux_h, grid, k)
        mf_tend_qt_dual = dual_faces(aux_tc_f.massflux_qt, grid, k)

        ∇mf_tend_h = ∇f2c(mf_tend_h_dual, grid, k)
        ∇mf_tend_qt = ∇f2c(mf_tend_qt_dual, grid, k)

        mf_tend_h = -∇mf_tend_h * α0_c[k]
        mf_tend_qt = -∇mf_tend_qt * α0_c[k]

        # Prepare the output
        aux_tc.massflux_tendency_h[k] = mf_tend_h
        aux_tc.massflux_tendency_qt[k] = mf_tend_qt
        tendencies_gm.θ_liq_ice[k] += mf_tend_h
        tendencies_gm.q_tot[k] += mf_tend_qt
    end

    aeKHq_tot_bc = Case.Sur.rho_qtflux / aux_en.area[kc_surf]
    aeKHθ_liq_ice_bc = Case.Sur.rho_hflux / aux_en.area[kc_surf]
    aeKMu_bc = Case.Sur.rho_uflux / aux_en.area[kc_surf]
    aeKMv_bc = Case.Sur.rho_vflux / aux_en.area[kc_surf]
    @inbounds for k in real_center_indices(grid)
        aeKH_q_tot_cut = dual_faces(aux_tc_f.diffusive_flux_qt, grid, k)
        ∇aeKH_q_tot = ∇f2c(aeKH_q_tot_cut, grid, k; bottom = SetValue(aeKHq_tot_bc), top = SetValue(0))
        tendencies_gm.q_tot[k] += -α0_c[k] * aux_en.area[k] * ∇aeKH_q_tot

        aeKH_θ_liq_ice_cut = dual_faces(aux_tc_f.diffusive_flux_h, grid, k)
        ∇aeKH_θ_liq_ice = ∇f2c(aeKH_θ_liq_ice_cut, grid, k; bottom = SetValue(aeKHθ_liq_ice_bc), top = SetValue(0))
        tendencies_gm.θ_liq_ice[k] += -α0_c[k] * aux_en.area[k] * ∇aeKH_θ_liq_ice

        aeKM_u_cut = dual_faces(aux_tc_f.diffusive_flux_u, grid, k)
        ∇aeKM_u = ∇f2c(aeKM_u_cut, grid, k; bottom = SetValue(aeKMu_bc), top = SetValue(0))
        tendencies_gm.u[k] += -α0_c[k] * aux_en.area[k] * ∇aeKM_u

        aeKM_v_cut = dual_faces(aux_tc_f.diffusive_flux_v, grid, k)
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
    aux_bulk = center_aux_bulk(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_en = center_aux_environment(state)
    aux_en.area .= 1 .- aux_bulk.area # area of environment
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
        aux_tc_f.diffusive_flux_qt[k] = -aux_tc_f.ρ_ae_KH[k] * ∇q_tot_f

        θ_liq_ice_dual = dual_centers(aux_en.θ_liq_ice, grid, k)
        ∇θ_liq_ice_f = ∇c2f(θ_liq_ice_dual, grid, k; bottom = SetGradient(aeKHθ_liq_ice_bc), top = SetGradient(0))
        aux_tc_f.diffusive_flux_h[k] = -aux_tc_f.ρ_ae_KH[k] * ∇θ_liq_ice_f

        u_dual = dual_centers(prog_gm.u, grid, k)
        ∇u_f = ∇c2f(u_dual, grid, k; bottom = SetGradient(aeKMu_bc), top = SetGradient(0))
        aux_tc_f.diffusive_flux_u[k] = -aux_tc_f.ρ_ae_KM[k] * ∇u_f

        v_dual = dual_centers(prog_gm.v, grid, k)
        ∇v_f = ∇c2f(v_dual, grid, k; bottom = SetGradient(aeKMv_bc), top = SetGradient(0))
        aux_tc_f.diffusive_flux_v[k] = -aux_tc_f.ρ_ae_KM[k] * ∇v_f
    end
    return
end

function affect_filter!(edmf, grid, state, gm, Case, TS)
    # TODO: figure out why this filter kills the DryBubble results if called at t = 0.
    if Case.casename == "DryBubble" && !(TS.t > 0)
        return nothing
    end
    prog_en = center_prog_environment(state)
    up = edmf.UpdVar
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
end

# Perform the update of the scheme
function step!(tendencies, prog, params, t)
    UnPack.@unpack edmf, grid, gm, case, aux, TS = params

    TS.t = t
    state = State(prog, aux, tendencies)

    affect_filter!(edmf, grid, state, gm, case, TS)

    gm = gm
    up = edmf.UpdVar
    en = edmf.EnvVar
    param_set = parameter_set(gm)
    en_thermo = edmf.EnvThermo

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.

    compute_updraft_surface_bc(edmf, grid, state, case)
    update_aux!(edmf, gm, grid, state, case, param_set, TS)

    parent(tendencies) .= 0

    # causes division error in dry bubble first time step
    compute_precipitation_formation_tendencies(grid, state, up, edmf.Precip, TS.dt, param_set)

    microphysics(en_thermo, grid, state, en, edmf.Precip, TS.dt, param_set) # saturation adjustment + rain creation
    if edmf.Precip.precipitation_model == "clima_1m"
        compute_precipitation_sink_tendencies(grid, state, gm, TS)
        compute_precipitation_advection_tendencies(grid, state, gm, TS)
    end

    # compute tendencies
    compute_gm_tendencies!(edmf, grid, state, case, gm, TS)
    compute_updraft_tendencies(edmf, grid, state, gm)

    compute_en_tendencies!(edmf, grid, state, param_set, TS, Val(:tke), Val(:ρatke))
    compute_en_tendencies!(edmf, grid, state, param_set, TS, Val(:Hvar), Val(:ρaHvar))
    compute_en_tendencies!(edmf, grid, state, param_set, TS, Val(:QTvar), Val(:ρaQTvar))
    compute_en_tendencies!(edmf, grid, state, param_set, TS, Val(:HQTcov), Val(:ρaHQTcov))

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
    aux_bulk = center_aux_bulk(state)
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
    ae = 1 .- aux_bulk.area # area of environment

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
    area_en = aux_en_c.area

    if is_tke
        Ic = CCO.InterpolateF2C()
        @. gmv_covar = tke_factor * area_en * Ic(ϕ_en - ϕ_gm) * Ic(ψ_en - ψ_gm) + area_en * covar_e
        @inbounds for i in 1:(edmf.n_updrafts)
            ϕ_up = getproperty(aux_up_f[i], ϕ_sym)
            ψ_up = getproperty(aux_up_f[i], ψ_sym)
            @. gmv_covar += tke_factor * aux_up[i].area * Ic(ϕ_up - ϕ_gm) * Ic(ψ_up - ψ_gm)
        end
    else
        @. gmv_covar = tke_factor * area_en * (ϕ_en - ϕ_gm) * (ψ_en - ψ_gm) + area_en * covar_e
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

function compute_updraft_tendencies(edmf::EDMF_PrognosticTKE{N_up}, grid, state, gm::GridMeanVariables) where {N_up}
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)

    aux_up = center_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_up_f = face_aux_updrafts(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0
    au_lim = edmf.max_area

    @inbounds for i in 1:N_up
        aux_up[i].entr_turb_dyn .= aux_up[i].entr_sc .+ aux_up[i].frac_turb_entr
        aux_up[i].detr_turb_dyn .= aux_up[i].detr_sc .+ aux_up[i].frac_turb_entr
    end

    UB = CCO.UpwindBiasedProductC2F(bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    Ic = CCO.InterpolateF2C()

    wvec = CC.Geometry.WVector
    ∂ = CCO.DivergenceF2C()
    w_bcs = (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))

    # Solve for updraft area fraction
    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        w_up = aux_up_f[i].w
        a_up = aux_up_i.area
        q_tot_up = aux_up_i.q_tot
        q_tot_en = aux_en.q_tot
        θ_liq_ice_en = aux_en.θ_liq_ice
        θ_liq_ice_up = aux_up_i.θ_liq_ice
        entr_turb_dyn = aux_up_i.entr_turb_dyn
        detr_turb_dyn = aux_up_i.detr_turb_dyn
        θ_liq_ice_tendency_precip_formation = aux_up_i.θ_liq_ice_tendency_precip_formation
        qt_tendency_precip_formation = aux_up_i.qt_tendency_precip_formation

        tends_ρarea = tendencies_up[i].ρarea
        tends_ρaθ_liq_ice = tendencies_up[i].ρaθ_liq_ice
        tends_ρaq_tot = tendencies_up[i].ρaq_tot

        @. tends_ρarea =
            -∂(wvec(LB(Ic(w_up) * ρ0_c * a_up))) + (ρ0_c * a_up * Ic(w_up) * entr_turb_dyn) -
            (ρ0_c * a_up * Ic(w_up) * detr_turb_dyn)

        @. tends_ρaθ_liq_ice =
            -∂(wvec(LB(Ic(w_up) * ρ0_c * a_up * θ_liq_ice_up))) +
            (ρ0_c * a_up * Ic(w_up) * entr_turb_dyn * θ_liq_ice_en) -
            (ρ0_c * a_up * Ic(w_up) * detr_turb_dyn * θ_liq_ice_up) + (ρ0_c * θ_liq_ice_tendency_precip_formation)

        @. tends_ρaq_tot =
            -∂(wvec(LB(Ic(w_up) * ρ0_c * a_up * q_tot_up))) + (ρ0_c * a_up * Ic(w_up) * entr_turb_dyn * q_tot_en) -
            (ρ0_c * a_up * Ic(w_up) * detr_turb_dyn * q_tot_up) + (ρ0_c * qt_tendency_precip_formation)

        tends_ρarea[kc_surf] = 0
        tends_ρaθ_liq_ice[kc_surf] = 0
        tends_ρaq_tot[kc_surf] = 0
    end

    # Solve for updraft velocity

    # We know that, since W = 0 at z = 0, BCs for entr, detr,
    # and buoyancy should not matter in the end
    zero_bcs = (; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    I0f = CCO.InterpolateC2F(; zero_bcs...)
    adv_bcs = (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    LB = CCO.LeftBiasedF2C(; bottom = CCO.SetValue(FT(0)))
    ∂ = CCO.DivergenceC2F(; adv_bcs...)

    @inbounds for i in 1:N_up
        a_up_bcs = (; bottom = CCO.SetValue(edmf.area_surface_bc[i]), top = CCO.SetValue(FT(0)))
        Iaf = CCO.InterpolateC2F(; a_up_bcs...)
        tends_ρaw = tendencies_up_f[i].ρaw
        nh_pressure = aux_up_f[i].nh_pressure
        a_up = aux_up[i].area
        w_up = aux_up_f[i].w
        w_en = aux_en_f.w
        entr_w = aux_up[i].entr_turb_dyn
        detr_w = aux_up[i].detr_turb_dyn
        buoy = aux_up[i].buoy

        @. tends_ρaw =
            -(∂(wvec(LB(Iaf(a_up) * ρ0_f * w_up * w_up)))) +
            (ρ0_f * Iaf(a_up) * w_up * (I0f(entr_w) * w_en - I0f(detr_w) * w_up)) +
            (ρ0_f * Iaf(a_up) * I0f(buoy)) +
            nh_pressure
        tends_ρaw[kf_surf] = 0
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
    wvec = CC.Geometry.WVector
    EnvVar1 = getproperty(aux_en, en_var1_sym)
    EnvVar2 = getproperty(aux_en, en_var2_sym)
    FT = eltype(grid)

    bcs = (; bottom = CCO.Extrapolate(), top = CCO.SetGradient(wvec(zero(FT))))
    If = CCO.InterpolateC2F(; bcs...)
    ∇c = CCO.DivergenceF2C()
    u = prog_gm.u
    v = prog_gm.v
    area_en = aux_en_c.area
    shear = aux_covar.shear

    if is_tke
        @. shear =
            tke_factor *
            2 *
            ρ0_c *
            area_en *
            k_eddy *
            (∇c(wvec(EnvVar1)) * ∇c(wvec(EnvVar2)) + (∇c(wvec(If(u))))^2 + (∇c(wvec(If(v))))^2)
    else
        @. shear = tke_factor * 2 * ρ0_c * area_en * k_eddy * ∇c(wvec(If(EnvVar1))) * ∇c(wvec(If(EnvVar2)))
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
    interdomain = getproperty(aux_en_2m, Covar_sym).interdomain
    prog_up = is_tke ? aux_up_f : aux_up
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    ϕ_en = getproperty(aux_en, ϕ_var)
    ψ_en = getproperty(aux_en, ψ_var)
    Ic = is_tke ? CCO.InterpolateF2C() : x -> x

    parent(interdomain) .= 0
    @inbounds for i in 1:(edmf.n_updrafts)
        ϕ_up = getproperty(prog_up[i], ϕ_var)
        ψ_up = getproperty(prog_up[i], ψ_var)
        a_up = aux_up[i].area
        @. interdomain += tke_factor * a_up * (1 - a_up) * (Ic(ϕ_up) - Ic(ϕ_en)) * (Ic(ψ_up) - Ic(ψ_en))
    end
    return
end

function compute_covariance_entr(
    edmf::EDMF_PrognosticTKE{N_up},
    grid,
    state,
    ::Val{covar_sym},
    v1::Val{var1},
    v2::Val{var2} = v1,
) where {covar_sym, var1, var2, N_up}

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
    entr_gain = aux_covar.entr_gain
    detr_loss = aux_covar.detr_loss

    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        frac_turb_entr = aux_up_i.frac_turb_entr
        detr_sc = aux_up_i.detr_sc
        entr_sc = aux_up_i.entr_sc
        w_up = aux_up_f[i].w
        prog_up_i = prog_up[i]
        @inbounds for k in real_center_indices(grid)
            entr_gain[k] = 0
            detr_loss[k] = 0
            a_up = aux_up_i.area[k]
            if a_up > edmf.minimum_area
                up_var_1 = getproperty(prog_up_i, var1)
                up_var_2 = getproperty(prog_up_i, var2)
                updvar1 = is_tke ? interpf2c(up_var_1, grid, k) : up_var_1[k]
                updvar2 = is_tke ? interpf2c(up_var_2, grid, k) : up_var_2[k]
                envvar1 = is_tke ? interpf2c(EnvVar1, grid, k) : EnvVar1[k]
                envvar2 = is_tke ? interpf2c(EnvVar2, grid, k) : EnvVar2[k]
                gmvvar1 = is_tke ? interpf2c(GmvVar1, grid, k) : GmvVar1[k]
                gmvvar2 = is_tke ? interpf2c(GmvVar2, grid, k) : GmvVar2[k]

                eps_turb = frac_turb_entr[k]

                w_u = interpf2c(w_up, grid, k)
                dynamic_entr =
                    tke_factor * ρ0_c[k] * a_up * abs(w_u) * detr_sc[k] * (updvar1 - envvar1) * (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    ρ0_c[k] *
                    a_up *
                    abs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                entr_gain[k] += dynamic_entr + turbulent_entr
                detr_loss[k] += tke_factor * ρ0_c[k] * a_up * abs(w_u) * (entr_sc[k] + eps_turb) * covar[k]
            end
        end
    end

    return
end

function compute_covariance_detr(edmf::EDMF_PrognosticTKE{N_up}, grid, state, covar_sym::Symbol) where {N_up}
    up = edmf.UpdVar
    ρ0_c = center_ref_state(state).ρ0
    aux_up_f = face_aux_updrafts(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_up = center_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    covar = getproperty(aux_en, covar_sym)
    detr_loss = aux_covar.detr_loss
    Ic = CCO.InterpolateF2C()

    parent(detr_loss) .= 0
    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        area_up = aux_up_i.area
        entr_sc = aux_up_i.entr_sc
        w_up = aux_up_f[i].w
        @. detr_loss += ρ0_c * covar * area_up * abs(Ic(w_up)) * entr_sc
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
    dissipation = aux_covar.dissipation
    area_en = aux_en.area
    tke_en = aux_en.tke
    mixing_length = aux_tc.mixing_length

    @. dissipation = ρ0_c * area_en * covar * max(tke_en, 0)^0.5 / max(mixing_length, 1.0e-3) * c_d
    return
end

function compute_en_tendencies!(
    edmf::EDMF_PrognosticTKE{N_up},
    grid::Grid,
    state,
    param_set,
    TS,
    ::Val{covar_sym},
    ::Val{prog_sym},
) where {N_up, covar_sym, prog_sym}
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
    c_d = CPEDMF.c_d(param_set)
    is_tke = covar_sym == :tke
    FT = eltype(grid)

    ρ_ae_K = face_aux_turbconv(state).ρ_ae_K
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    aux_tc = center_aux_turbconv(state)
    aux_bulk = center_aux_bulk(state)
    D_env = aux_bulk.D_env
    ae = 1 .- aux_bulk.area
    aeK = is_tke ? ae .* KM : ae .* KH

    press = aux_covar.press
    buoy = aux_covar.buoy
    shear = aux_covar.shear
    entr_gain = aux_covar.entr_gain
    rain_src = aux_covar.rain_src

    wvec = CC.Geometry.WVector
    aeK_bcs = (; bottom = CCO.SetValue(aeK[kc_surf]), top = CCO.SetValue(aeK[kc_toa]))
    prog_bcs = (; bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))

    If = CCO.InterpolateC2F(; aeK_bcs...)
    ∇f = CCO.GradientC2F(; prog_bcs...)
    ∇c = CCO.DivergenceF2C()

    mixing_length = aux_tc.mixing_length
    min_area = edmf.minimum_area
    pressure_plume_spacing = edmf.pressure_plume_spacing

    Ic = CCO.InterpolateF2C()
    area_en = aux_en.area
    tke_en = aux_en.tke

    parent(D_env) .= 0

    @inbounds for i in 1:N_up
        turb_entr = aux_up[i].frac_turb_entr
        entr_sc = aux_up[i].entr_sc
        w_up = aux_up_f[i].w
        a_up = aux_up[i].area
        # TODO: using `Int(bool) *` means that NaNs can propagate
        # into the solution. Could we somehow call `ifelse` instead?
        @. D_env += Int(a_up > min_area) * ρ0_c * a_up * Ic(w_up) * (entr_sc + turb_entr)
    end

    RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0)))
    @. prog_covar =
        press + buoy + shear + entr_gain + rain_src - D_env * covar -
        (ρ0_c * area_en * c_d * sqrt(max(tke_en, 0)) / max(mixing_length, 1)) * covar -
        ∇c(wvec(RB(ρ0_c * area_en * Ic(w_en_f) * covar))) + ∇c(ρ0_f * If(aeK) * ∇f(covar))

    prog_covar[kc_surf] = covar[kc_surf]

    return nothing
end

function GMV_third_m(
    edmf::EDMF_PrognosticTKE{N_up},
    grid,
    state,
    ::Val{covar_en_sym},
    ::Val{var},
    ::Val{gm_third_m_sym},
) where {N_up, covar_en_sym, var, gm_third_m_sym}

    gm_third_m = getproperty(center_aux_grid_mean(state), gm_third_m_sym)
    kc_surf = kc_surface(grid)

    aux_bulk = center_aux_bulk(state)
    aux_up_f = face_aux_updrafts(state)
    is_tke = covar_en_sym == :tke
    aux_en_c = center_aux_environment(state)
    covar_en = getproperty(aux_en_c, covar_en_sym)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    aux_up_c = center_aux_updrafts(state)
    aux_tc = center_aux_turbconv(state)
    ϕ_gm = aux_tc.ϕ_gm
    ϕ_gm_cov = aux_tc.ϕ_gm_cov
    ϕ_en_cov = aux_tc.ϕ_en_cov
    ϕ_up_cubed = aux_tc.ϕ_up_cubed
    aux_up = is_tke ? aux_up_f : aux_up_c
    var_en = getproperty(aux_en, var)
    area_en = aux_en_c.area
    Ic = is_tke ? CCO.InterpolateF2C() : x -> x
    wvec = CC.Geometry.WVector
    ∇c = CCO.DivergenceF2C()
    w_en = aux_en_f.w

    @. ϕ_gm = area_en * Ic(var_en)
    @inbounds for i in 1:N_up
        a_up = aux_up_c[i].area
        var_up = getproperty(aux_up[i], var)
        @. ϕ_gm += a_up * Ic(var_up)
    end

    if is_tke
        parent(ϕ_en_cov) .= 0
        @inbounds for i in 1:N_up
            horiz_K_eddy = aux_up_c[i].horiz_K_eddy
            a_up = aux_up_c[i].area
            a_bulk = aux_bulk.area
            @. ϕ_en_cov += -horiz_K_eddy * ∇c(wvec(w_en)) * a_up / a_bulk
        end
    else
        @. ϕ_en_cov = covar_en
    end

    parent(ϕ_up_cubed) .= 0
    @. ϕ_gm_cov = area_en * (ϕ_en_cov + (Ic(var_en) - ϕ_gm)^2)
    @inbounds for i in 1:N_up
        a_up = aux_up_c[i].area
        var_up = getproperty(aux_up[i], var)
        @. ϕ_gm_cov += a_up * (Ic(var_up) - ϕ_gm)^2
        @. ϕ_up_cubed += a_up * Ic(var_up)^3
    end

    @. gm_third_m = ϕ_up_cubed + area_en * (Ic(var_en)^3 + 3 * Ic(var_en) * ϕ_en_cov) - ϕ_gm^3 - 3 * ϕ_gm_cov * ϕ_gm

    gm_third_m[kc_surf] = 0
    return
end
