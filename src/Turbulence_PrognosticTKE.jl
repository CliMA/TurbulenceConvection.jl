
#= These methods are to be overloaded by Cases.jl =#
function update_surface end
function update_forcing end
function update_radiation end

function update_cloud_frac(edmf::EDMF_PrognosticTKE, grid::Grid, state::State, gm::GridMeanVariables)
    # update grid-mean cloud fraction and cloud cover
    aux_bulk = center_aux_bulk(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    a_up_bulk = aux_bulk.area
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        aux_gm.cloud_fraction[k] = aux_en.area[k] * aux_en.cloud_fraction[k] + a_up_bulk[k] * aux_bulk.cloud_fraction[k]
    end
end

function compute_les_Γᵣ(
    z::ClimaCore.Geometry.ZPoint,
    τᵣ::FT = 24.0 * 3600.0,
    zᵢ::FT = 3000.0,
    zᵣ::FT = 3500.0,
) where {FT <: Real}
    # returns height-dependent relaxation timescale from eqn. 9 in `Shen et al. 2021`
    if z < zᵢ
        return FT(0)
    elseif zᵢ <= z <= zᵣ
        cos_arg = pi * ((z - zᵢ) / (zᵣ - zᵢ))
        return (FT(0.5) / τᵣ) * (1 - cos(cos_arg))
    elseif z > zᵣ
        return (1 / τᵣ)
    end
end

function compute_gm_tendencies!(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    surf::SurfaceBase,
    radiation::RadiationBase,
    force::ForcingBase,
    gm::GridMeanVariables,
)
    N_up = n_updrafts(edmf)
    tendencies_gm = center_tendencies_grid_mean(state)
    kc_toa = kc_top_of_atmos(grid)
    FT = eltype(grid)
    param_set = parameter_set(gm)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    ∇θ_liq_ice_gm = center_aux_grid_mean(state).∇θ_liq_ice_gm
    ∇q_tot_gm = center_aux_grid_mean(state).∇q_tot_gm
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
    massflux_h = aux_tc_f.massflux_h
    massflux_qt = aux_tc_f.massflux_qt
    aux_tc = center_aux_turbconv(state)

    θ_liq_ice_gm_toa = prog_gm.θ_liq_ice[kc_toa]
    q_tot_gm_toa = prog_gm.q_tot[kc_toa]
    RBθ = CCO.RightBiasedC2F(; top = CCO.SetValue(θ_liq_ice_gm_toa))
    RBq = CCO.RightBiasedC2F(; top = CCO.SetValue(q_tot_gm_toa))
    wvec = CC.Geometry.WVector
    ∇c = CCO.DivergenceF2C()
    @. ∇θ_liq_ice_gm = ∇c(wvec(RBθ(prog_gm.θ_liq_ice)))
    @. ∇q_tot_gm = ∇c(wvec(RBq(prog_gm.q_tot)))

    @inbounds for k in real_center_indices(grid)
        # Apply large-scale horizontal advection tendencies
        ts = thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(ts)

        if force.apply_coriolis
            tendencies_gm.u[k] -= force.coriolis_param * (aux_gm.vg[k] - prog_gm.v[k])
            tendencies_gm.v[k] += force.coriolis_param * (aux_gm.ug[k] - prog_gm.u[k])
        end
        if rad_type(radiation) <: Union{RadiationDYCOMS_RF01, RadiationLES}
            tendencies_gm.θ_liq_ice[k] += aux_gm.dTdt_rad[k] / Π
        end

        if force_type(force) <: ForcingDYCOMS_RF01
            tendencies_gm.q_tot[k] += aux_gm.dqtdt[k]
            # Apply large-scale subsidence tendencies
            tendencies_gm.θ_liq_ice[k] -= ∇θ_liq_ice_gm[k] * aux_gm.subsidence[k]
            tendencies_gm.q_tot[k] -= ∇q_tot_gm[k] * aux_gm.subsidence[k]
        end

        if force_type(force) <: ForcingStandard
            if force.apply_subsidence
                tendencies_gm.θ_liq_ice[k] -= ∇θ_liq_ice_gm[k] * aux_gm.subsidence[k]
                tendencies_gm.q_tot[k] -= ∇q_tot_gm[k] * aux_gm.subsidence[k]
            end
            tendencies_gm.θ_liq_ice[k] += aux_gm.dTdt[k] / Π
            tendencies_gm.q_tot[k] += aux_gm.dqtdt[k]
        end

        if force_type(force) <: ForcingLES
            H_horz_adv = aux_gm.dTdt_hadv[k] / Π
            H_fluc = aux_gm.dTdt_fluc[k] / Π

            gm_U_nudge_k = (aux_gm.u_nudge[k] - prog_gm.u[k]) / force.nudge_tau
            gm_V_nudge_k = (aux_gm.v_nudge[k] - prog_gm.v[k]) / force.nudge_tau

            Γᵣ = compute_les_Γᵣ(grid.zc[k])
            if Γᵣ != 0
                tau_k = 1.0 / Γᵣ
                gm_H_nudge_k = (aux_gm.H_nudge[k] - prog_gm.θ_liq_ice[k]) / tau_k
                gm_q_tot_nudge_k = (aux_gm.qt_nudge[k] - prog_gm.q_tot[k]) / tau_k
            else
                gm_H_nudge_k = 0.0
                gm_q_tot_nudge_k = 0.0
            end

            if force.apply_subsidence
                # Apply large-scale subsidence tendencies
                gm_H_subsidence_k = -∇θ_liq_ice_gm[k] * aux_gm.subsidence[k]
                gm_QT_subsidence_k = -∇q_tot_gm[k] * aux_gm.subsidence[k]
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
    a_en = aux_en.area
    w_en = aux_en_f.w
    a_en_bcs = a_en_boundary_conditions(surf, edmf)
    Ifae = CCO.InterpolateC2F(; a_en_bcs...)
    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in 1:N_up
        a_up_bcs = a_up_boundary_conditions(surf, edmf, i)
        Ifau = CCO.InterpolateC2F(; a_up_bcs...)
        a_up = aux_up[i].area
        w_up = aux_up_f[i].w
        massflux = aux_up_f[i].massflux
        @. massflux = ρ0_f * Ifau(a_up) * Ifae(a_en) * (w_up - w_en)
        massflux[kf_surf] = 0
    end

    parent(massflux_h) .= 0
    parent(massflux_qt) .= 0
    If = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    θ_liq_ice_en = aux_en.θ_liq_ice
    q_tot_en = aux_en.q_tot
    @inbounds for i in 1:N_up
        aux_up_f_i = aux_up_f[i]
        aux_up_i = aux_up[i]
        θ_liq_ice_up = aux_up_i.θ_liq_ice
        q_tot_up = aux_up_i.q_tot
        massflux = aux_up_f[i].massflux
        # We know that, since W = 0 at z = 0, m = 0 also, and
        # therefore θ_liq_ice / q_tot values do not matter
        @. massflux_h += massflux * (If(θ_liq_ice_up) - If(θ_liq_ice_en))
        @. massflux_qt += massflux * (If(q_tot_up) - If(q_tot_en))
    end

    # Compute the  mass flux tendencies
    # Adjust the values of the grid mean variables
    # Prepare the output
    @. aux_tc.massflux_tendency_h = -∇c(wvec(massflux_h)) * α0_c
    @. aux_tc.massflux_tendency_qt = ∇c(wvec(massflux_qt))
    @. tendencies_gm.θ_liq_ice += -∇c(wvec(massflux_h)) * α0_c
    @. tendencies_gm.q_tot += -∇c(wvec(massflux_qt)) * α0_c

    a_en = aux_en.area
    aeKHq_tot_bc = surf.ρq_tot_flux / a_en[kc_surf]
    aeKHθ_liq_ice_bc = surf.ρθ_liq_ice_flux / a_en[kc_surf]
    aeKMu_bc = surf.ρu_flux / a_en[kc_surf]
    aeKMv_bc = surf.ρv_flux / a_en[kc_surf]

    tbc = (; top = CCO.SetValue(wvec(FT(0))))
    ∇aeKH_q_tot = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(aeKHq_tot_bc)), tbc...)
    ∇aeKH_θ_liq_ice = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(aeKHθ_liq_ice_bc)), tbc...)
    ∇aeKM_u = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(aeKMu_bc)), tbc...)
    ∇aeKM_v = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(aeKMv_bc)), tbc...)

    @. tendencies_gm.q_tot += -α0_c * a_en * ∇aeKH_q_tot(wvec(aux_tc_f.diffusive_flux_qt))
    @. tendencies_gm.θ_liq_ice += -α0_c * a_en * ∇aeKH_θ_liq_ice(wvec(aux_tc_f.diffusive_flux_h))
    @. tendencies_gm.u += -α0_c * a_en * ∇aeKM_u(wvec(aux_tc_f.diffusive_flux_u))
    @. tendencies_gm.v += -α0_c * a_en * ∇aeKM_v(wvec(aux_tc_f.diffusive_flux_v))

    return nothing
end

function compute_diffusive_fluxes(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    gm::GridMeanVariables,
    surf::SurfaceBase,
    param_set::APS,
)
    FT = eltype(grid)
    ρ0_f = face_ref_state(state).ρ0
    aux_bulk = center_aux_bulk(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_en = center_aux_environment(state)
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    aeKM = aux_en.area .* KM
    aeKH = aux_en.area .* KH
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_surf = kf_surface(grid)
    prog_gm = center_prog_grid_mean(state)
    IfKH = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKH[kc_surf]), top = CCO.SetValue(aeKH[kc_toa]))
    IfKM = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKM[kc_surf]), top = CCO.SetValue(aeKM[kc_toa]))

    @. aux_tc_f.ρ_ae_KH = IfKH(aeKH) * ρ0_f
    @. aux_tc_f.ρ_ae_KM = IfKM(aeKM) * ρ0_f

    aeKHq_tot_bc = -surf.ρq_tot_flux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KH[kf_surf]
    aeKHθ_liq_ice_bc = -surf.ρθ_liq_ice_flux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KH[kf_surf]
    aeKMu_bc = -surf.ρu_flux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KM[kf_surf]
    aeKMv_bc = -surf.ρv_flux / aux_en.area[kc_surf] / aux_tc_f.ρ_ae_KM[kf_surf]

    ∇q_tot_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKHq_tot_bc), top = CCO.SetDivergence(FT(0)))
    ∇θ_liq_ice_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKHθ_liq_ice_bc), top = CCO.SetDivergence(FT(0)))
    ∇u_gm = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKMu_bc), top = CCO.SetDivergence(FT(0)))
    ∇v_gm = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKMv_bc), top = CCO.SetDivergence(FT(0)))

    wvec = CC.Geometry.WVector
    @. aux_tc_f.diffusive_flux_qt = -aux_tc_f.ρ_ae_KH * ∇q_tot_en(wvec(aux_en.q_tot))
    @. aux_tc_f.diffusive_flux_h = -aux_tc_f.ρ_ae_KH * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice))
    @. aux_tc_f.diffusive_flux_u = -aux_tc_f.ρ_ae_KM * ∇u_gm(wvec(prog_gm.u))
    @. aux_tc_f.diffusive_flux_v = -aux_tc_f.ρ_ae_KM * ∇v_gm(wvec(prog_gm.v))

    return nothing
end

function affect_filter!(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    gm::GridMeanVariables,
    surf::SurfaceBase,
    casename::String,
    t::Real,
)
    # TODO: figure out why this filter kills the DryBubble results if called at t = 0.
    if casename == "DryBubble" && !(t > 0)
        return nothing
    end
    prog_en = center_prog_environment(state)
    ###
    ### Filters
    ###
    set_edmf_surface_bc(edmf, grid, state, surf)
    filter_updraft_vars(edmf, grid, state, surf, gm)

    @inbounds for k in real_center_indices(grid)
        prog_en.ρatke[k] = max(prog_en.ρatke[k], 0.0)
        prog_en.ρaHvar[k] = max(prog_en.ρaHvar[k], 0.0)
        prog_en.ρaQTvar[k] = max(prog_en.ρaQTvar[k], 0.0)
        prog_en.ρaHQTcov[k] = max(prog_en.ρaHQTcov[k], -sqrt(prog_en.ρaHvar[k] * prog_en.ρaQTvar[k]))
        prog_en.ρaHQTcov[k] = min(prog_en.ρaHQTcov[k], sqrt(prog_en.ρaHvar[k] * prog_en.ρaQTvar[k]))
    end
    return nothing
end

# Compute the sum of tendencies for the scheme
function ∑tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    UnPack.@unpack edmf, grid, gm, case, aux, TS = params

    state = State(prog, aux, tendencies)

    Δt = TS.dt
    param_set = parameter_set(gm)
    surf = get_surface(case.surf_params, grid, state, gm, t, param_set)
    force = case.Fo
    radiation = case.Rad
    en_thermo = edmf.en_thermo
    precip_model = edmf.precip_model

    affect_filter!(edmf, grid, state, gm, surf, case.casename, t)

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.

    update_aux!(edmf, gm, grid, state, case, param_set, t, Δt)

    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0

    # causes division error in dry bubble first time step
    compute_precipitation_formation_tendencies(grid, state, edmf, precip_model, Δt, param_set)

    microphysics(en_thermo, grid, state, precip_model, Δt, param_set)
    compute_precipitation_sink_tendencies(precip_model, grid, state, gm, Δt)
    compute_precipitation_advection_tendencies(precip_model, edmf, grid, state, gm)

    # compute tendencies
    compute_gm_tendencies!(edmf, grid, state, surf, radiation, force, gm)
    compute_updraft_tendencies(edmf, grid, state, gm, surf)

    compute_en_tendencies!(edmf, grid, state, param_set, Val(:tke), Val(:ρatke))
    compute_en_tendencies!(edmf, grid, state, param_set, Val(:Hvar), Val(:ρaHvar))
    compute_en_tendencies!(edmf, grid, state, param_set, Val(:QTvar), Val(:ρaQTvar))
    compute_en_tendencies!(edmf, grid, state, param_set, Val(:HQTcov), Val(:ρaHQTcov))

    return nothing
end

function set_edmf_surface_bc(edmf::EDMF_PrognosticTKE, grid::Grid, state::State, surf::SurfaceBase)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0
    prog_up = center_prog_updrafts(state)
    prog_en = center_prog_environment(state)
    prog_up_f = face_prog_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    @inbounds for i in 1:N_up
        θ_surf = θ_surface_bc(surf, grid, state, edmf, i)
        q_surf = q_surface_bc(surf, grid, state, edmf, i)
        a_surf = area_surface_bc(surf, edmf, i)
        prog_up[i].ρarea[kc_surf] = ρ0_c[kc_surf] * a_surf
        prog_up[i].ρaθ_liq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * θ_surf
        prog_up[i].ρaq_tot[kc_surf] = prog_up[i].ρarea[kc_surf] * q_surf
        prog_up_f[i].ρaw[kf_surf] = ρ0_f[kf_surf] * w_surface_bc(surf)
    end

    flux1 = surf.ρθ_liq_ice_flux
    flux2 = surf.ρq_tot_flux
    zLL = grid.zc[kc_surf].z
    ustar = surf.ustar
    oblength = surf.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]
    # TODO: is bulk even defined before this is called?
    ae = 1 .- aux_bulk.area # area of environment

    ρ0_ae = ρ0_c[kc_surf] * ae[kc_surf]

    prog_en.ρatke[kc_surf] = ρ0_ae * get_surface_tke(surf.ustar, zLL, surf.obukhov_length)
    prog_en.ρaHvar[kc_surf] = ρ0_ae * get_surface_variance(flux1 * α0LL, flux1 * α0LL, ustar, zLL, oblength)
    prog_en.ρaQTvar[kc_surf] = ρ0_ae * get_surface_variance(flux2 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    prog_en.ρaHQTcov[kc_surf] = ρ0_ae * get_surface_variance(flux1 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    return nothing
end

function surface_helper(surf::SurfaceBase, grid::Grid, state::State)
    kc_surf = kc_surface(grid)
    zLL = grid.zc[kc_surf].z
    ustar = surf.ustar
    oblength = surf.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]
    return (; ustar, zLL, oblength, α0LL)
end

function a_up_boundary_conditions(surf::SurfaceBase, edmf::EDMF_PrognosticTKE, i::Int)
    a_surf = area_surface_bc(surf, edmf, i)
    return (; bottom = CCO.SetValue(a_surf), top = CCO.Extrapolate())
end

function a_bulk_boundary_conditions(surf::SurfaceBase, edmf::EDMF_PrognosticTKE)
    N_up = n_updrafts(edmf)
    a_surf = sum(i -> area_surface_bc(surf, edmf, i), 1:N_up)
    return (; bottom = CCO.SetValue(a_surf), top = CCO.Extrapolate())
end

function a_en_boundary_conditions(surf::SurfaceBase, edmf::EDMF_PrognosticTKE)
    N_up = n_updrafts(edmf)
    a_surf = 1 - sum(i -> area_surface_bc(surf, edmf, i), 1:N_up)
    return (; bottom = CCO.SetValue(a_surf), top = CCO.Extrapolate())
end

function area_surface_bc(surf::SurfaceBase{FT}, edmf::EDMF_PrognosticTKE, i::Int)::FT where {FT}
    N_up = n_updrafts(edmf)
    return surf.bflux > 0 ? edmf.surface_area / N_up : FT(0)
end

function w_surface_bc(::SurfaceBase{FT})::FT where {FT}
    return FT(0)
end
function θ_surface_bc(surf::SurfaceBase{FT}, grid::Grid, state::State, edmf::EDMF_PrognosticTKE, i::Int)::FT where {FT}
    prog_gm = center_prog_grid_mean(state)
    kc_surf = kc_surface(grid)
    surf.bflux > 0 || return FT(0)
    a_total = edmf.surface_area
    a_ = area_surface_bc(surf, edmf, i)
    UnPack.@unpack ustar, zLL, oblength, α0LL = surface_helper(surf, grid, state)
    ρθ_liq_ice_flux = surf.ρθ_liq_ice_flux
    h_var = get_surface_variance(ρθ_liq_ice_flux * α0LL, ρθ_liq_ice_flux * α0LL, ustar, zLL, oblength)
    surface_scalar_coeff = percentile_bounds_mean_norm(1 - a_total + i * a_, 1 - a_total + (i + 1) * a_, 1000)
    return prog_gm.θ_liq_ice[kc_surf] + surface_scalar_coeff * sqrt(h_var)
end
function q_surface_bc(surf::SurfaceBase{FT}, grid::Grid, state::State, edmf::EDMF_PrognosticTKE, i::Int)::FT where {FT}
    prog_gm = center_prog_grid_mean(state)
    kc_surf = kc_surface(grid)
    surf.bflux > 0 || return prog_gm.q_tot[kc_surf]
    a_total = edmf.surface_area
    a_ = area_surface_bc(surf, edmf, i)
    UnPack.@unpack ustar, zLL, oblength, α0LL = surface_helper(surf, grid, state)
    ρq_tot_flux = surf.ρq_tot_flux
    qt_var = get_surface_variance(ρq_tot_flux * α0LL, ρq_tot_flux * α0LL, ustar, zLL, oblength)
    surface_scalar_coeff = percentile_bounds_mean_norm(1 - a_total + i * a_, 1 - a_total + (i + 1) * a_, 1000)
    return prog_gm.q_tot[kc_surf] + surface_scalar_coeff * sqrt(qt_var)
end

function get_GMV_CoVar(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_sym},
    ::Val{ψ_sym},
) where {covar_sym, ϕ_sym, ψ_sym}
    N_up = n_updrafts(edmf)
    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    aux_up_c = center_aux_updrafts(state)
    aux_up = is_tke ? aux_up_f : aux_up_c
    gmv_covar = getproperty(center_aux_grid_mean(state), covar_sym)
    covar_e = getproperty(center_aux_environment(state), covar_sym)
    prog_gm = is_tke ? prog_gm_f : prog_gm_c
    ϕ_gm = getproperty(prog_gm, ϕ_sym)
    ψ_gm = getproperty(prog_gm, ψ_sym)
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    area_en = aux_en_c.area

    Icd = is_tke ? CCO.InterpolateF2C() : x -> x
    @. gmv_covar = tke_factor * area_en * Icd(ϕ_en - ϕ_gm) * Icd(ψ_en - ψ_gm) + area_en * covar_e
    @inbounds for i in 1:N_up
        ϕ_up = getproperty(aux_up[i], ϕ_sym)
        ψ_up = getproperty(aux_up[i], ψ_sym)
        @. gmv_covar += tke_factor * aux_up_c[i].area * Icd(ϕ_up - ϕ_gm) * Icd(ψ_up - ψ_gm)
    end
    return nothing
end

function compute_updraft_top(grid::Grid{FT}, state::State, i::Int)::FT where {FT}
    aux_up = center_aux_updrafts(state)
    return z_findlast_center(k -> aux_up[i].area[k] > 1e-3, grid)
end

function compute_pressure_plume_spacing(
    edmf::EDMF_PrognosticTKE,
    grid::Grid{FT},
    state::State,
    param_set::APS,
    i::Int,
)::FT where {FT}
    H_up_min::FT = CPEDMF.H_up_min(param_set)
    updraft_top = compute_updraft_top(grid, state, i)
    aspect_ratio = edmf.aspect_ratio
    return max(aspect_ratio * updraft_top, H_up_min * aspect_ratio)
end

function compute_updraft_tendencies(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    gm::GridMeanVariables,
    surf::SurfaceBase,
)
    N_up = n_updrafts(edmf)
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
        aux_up_i = aux_up[i]
        @. aux_up_i.entr_turb_dyn = aux_up_i.entr_sc + aux_up_i.frac_turb_entr
        @. aux_up_i.detr_turb_dyn = aux_up_i.detr_sc + aux_up_i.frac_turb_entr
    end

    UB = CCO.UpwindBiasedProductC2F(bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    Ic = CCO.InterpolateF2C()

    wvec = CC.Geometry.WVector
    ∇c = CCO.DivergenceF2C()
    w_bcs = (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    LBF = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))

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
            -∇c(wvec(LBF(Ic(w_up) * ρ0_c * a_up))) + (ρ0_c * a_up * Ic(w_up) * entr_turb_dyn) -
            (ρ0_c * a_up * Ic(w_up) * detr_turb_dyn)

        @. tends_ρaθ_liq_ice =
            -∇c(wvec(LBF(Ic(w_up) * ρ0_c * a_up * θ_liq_ice_up))) +
            (ρ0_c * a_up * Ic(w_up) * entr_turb_dyn * θ_liq_ice_en) -
            (ρ0_c * a_up * Ic(w_up) * detr_turb_dyn * θ_liq_ice_up) + (ρ0_c * θ_liq_ice_tendency_precip_formation)

        @. tends_ρaq_tot =
            -∇c(wvec(LBF(Ic(w_up) * ρ0_c * a_up * q_tot_up))) + (ρ0_c * a_up * Ic(w_up) * entr_turb_dyn * q_tot_en) -
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
    LBC = CCO.LeftBiasedF2C(; bottom = CCO.SetValue(FT(0)))
    ∇f = CCO.DivergenceC2F(; adv_bcs...)

    @inbounds for i in 1:N_up
        a_up_bcs = a_up_boundary_conditions(surf, edmf, i)
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
            -(∇f(wvec(LBC(Iaf(a_up) * ρ0_f * w_up * w_up)))) +
            (ρ0_f * Iaf(a_up) * w_up * (I0f(entr_w) * w_en - I0f(detr_w) * w_up)) +
            (ρ0_f * Iaf(a_up) * I0f(buoy)) +
            nh_pressure
        tends_ρaw[kf_surf] = 0
    end

    return nothing
end

function filter_updraft_vars(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    surf::SurfaceBase,
    gm::GridMeanVariables,
)
    N_up = n_updrafts(edmf)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    N_up = n_updrafts(edmf)

    prog_up = center_prog_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0
    a_min = edmf.minimum_area
    a_max = edmf.max_area

    @inbounds for i in 1:N_up
        prog_up[i].ρarea .= max.(prog_up[i].ρarea, 0)
        prog_up[i].ρaθ_liq_ice .= max.(prog_up[i].ρaθ_liq_ice, 0)
        prog_up[i].ρaq_tot .= max.(prog_up[i].ρaq_tot, 0)
        @inbounds for k in real_center_indices(grid)
            prog_up[i].ρarea[k] = min(prog_up[i].ρarea[k], ρ0_c[k] * a_max)
        end
    end

    @inbounds for i in 1:N_up
        @. prog_up_f[i].ρaw = max.(prog_up_f[i].ρaw, 0)
        a_up_bcs = a_up_boundary_conditions(surf, edmf, i)
        If = CCO.InterpolateC2F(; a_up_bcs...)
        @. prog_up_f[i].ρaw = Int(If(prog_up[i].ρarea) >= ρ0_f * a_min) * prog_up_f[i].ρaw
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:N_up
            is_surface_center(grid, k) && continue
            prog_up[i].ρaq_tot[k] = max(prog_up[i].ρaq_tot[k], 0)
            # this is needed to make sure Rico is unchanged.
            # TODO : look into it further to see why
            # a similar filtering of ρaθ_liq_ice breaks the simulation
            if prog_up[i].ρarea[k] / ρ0_c[k] < a_min
                prog_up[i].ρaq_tot[k] = 0
            end
        end
    end

    Ic = CCO.InterpolateF2C()
    @inbounds for i in 1:N_up
        @. prog_up[i].ρarea = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρarea)
        @. prog_up[i].ρaθ_liq_ice = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρaθ_liq_ice)
        @. prog_up[i].ρaq_tot = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρaq_tot)
        θ_surf = θ_surface_bc(surf, grid, state, edmf, i)
        q_surf = q_surface_bc(surf, grid, state, edmf, i)
        a_surf = area_surface_bc(surf, edmf, i)
        prog_up[i].ρarea[kc_surf] = ρ0_c[kc_surf] * a_surf
        prog_up[i].ρaθ_liq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * θ_surf
        prog_up[i].ρaq_tot[kc_surf] = prog_up[i].ρarea[kc_surf] * q_surf
    end
    return nothing
end

function compute_covariance_shear(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    gm::GridMeanVariables,
    ::Val{covar_sym},
    ::Val{ϕ_en_sym},
    ::Val{ψ_en_sym},
) where {covar_sym, ϕ_en_sym, ψ_en_sym}

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
    ϕ_en = getproperty(aux_en, ϕ_en_sym)
    ψ_en = getproperty(aux_en, ψ_en_sym)
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
            (∇c(wvec(ϕ_en)) * ∇c(wvec(ψ_en)) + (∇c(wvec(If(u))))^2 + (∇c(wvec(If(v))))^2)
    else
        @. shear = tke_factor * 2 * ρ0_c * area_en * k_eddy * ∇c(wvec(If(ϕ_en))) * ∇c(wvec(If(ψ_en)))
    end
    return nothing
end

function compute_covariance_interdomain_src(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_sym},
    ::Val{ψ_sym},
) where {covar_sym, ϕ_sym, ψ_sym}
    N_up = n_updrafts(edmf)
    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_2m = center_aux_environment_2m(state)
    interdomain = getproperty(aux_en_2m, covar_sym).interdomain
    prog_up = is_tke ? aux_up_f : aux_up
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    Ic = is_tke ? CCO.InterpolateF2C() : x -> x

    parent(interdomain) .= 0
    @inbounds for i in 1:N_up
        ϕ_up = getproperty(prog_up[i], ϕ_sym)
        ψ_up = getproperty(prog_up[i], ψ_sym)
        a_up = aux_up[i].area
        @. interdomain += tke_factor * a_up * (1 - a_up) * (Ic(ϕ_up) - Ic(ϕ_en)) * (Ic(ψ_up) - Ic(ψ_en))
    end
    return nothing
end

function compute_covariance_entr(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_sym},
    ::Val{ψ_sym},
) where {covar_sym, ϕ_sym, ψ_sym}

    N_up = n_updrafts(edmf)
    ρ0_c = center_ref_state(state).ρ0
    FT = eltype(grid)
    is_tke = covar_sym == :tke
    tke_factor = is_tke ? 0.5 : 1
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    prog_gm = is_tke ? prog_gm_f : prog_gm_c
    prog_up = is_tke ? aux_up_f : aux_up
    ϕ_gm = getproperty(prog_gm, ϕ_sym)
    ψ_gm = getproperty(prog_gm, ψ_sym)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_en_c = center_aux_environment(state)
    covar = getproperty(aux_en_c, covar_sym)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    entr_gain = aux_covar.entr_gain
    detr_loss = aux_covar.detr_loss
    Ic = CCO.InterpolateF2C()
    Idc = is_tke ? Ic : x -> x
    # TODO: we shouldn't need `parent` call here:
    parent(entr_gain) .= 0
    parent(detr_loss) .= 0
    min_area = edmf.minimum_area

    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        frac_turb_entr = aux_up_i.frac_turb_entr
        eps_turb = frac_turb_entr
        detr_sc = aux_up_i.detr_sc
        entr_sc = aux_up_i.entr_sc
        w_up = aux_up_f[i].w
        prog_up_i = prog_up[i]
        ϕ_up = getproperty(prog_up_i, ϕ_sym)
        ψ_up = getproperty(prog_up_i, ψ_sym)

        a_up = aux_up_i.area

        @. entr_gain +=
            Int(a_up > min_area) *
            (tke_factor * ρ0_c * a_up * abs(Ic(w_up)) * detr_sc * (Idc(ϕ_up) - Idc(ϕ_en)) * (Idc(ψ_up) - Idc(ψ_en))) + (
                tke_factor *
                ρ0_c *
                a_up *
                abs(Ic(w_up)) *
                eps_turb *
                ((Idc(ϕ_en) - Idc(ϕ_gm)) * (Idc(ψ_up) - Idc(ψ_en)) + (Idc(ψ_en) - Idc(ψ_gm)) * (Idc(ϕ_up) - Idc(ϕ_en)))
            )

        @. detr_loss += Int(a_up > min_area) * tke_factor * ρ0_c * a_up * abs(Ic(w_up)) * (entr_sc + eps_turb) * covar

    end

    return nothing
end

function compute_covariance_dissipation(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    ::Val{covar_sym},
    param_set::APS,
) where {covar_sym}
    FT = eltype(grid)
    c_d::FT = CPEDMF.c_d(param_set)
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

    @. dissipation = ρ0_c * area_en * covar * max(tke_en, 0)^FT(0.5) / max(mixing_length, FT(1.0e-3)) * c_d
    return nothing
end

function compute_en_tendencies!(
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    param_set::APS,
    ::Val{covar_sym},
    ::Val{prog_sym},
) where {S, PS, covar_sym, prog_sym}
    N_up = n_updrafts(edmf)
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
    D_env = aux_tc.ϕ_temporary
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
    edmf::EDMF_PrognosticTKE,
    grid::Grid,
    state::State,
    ::Val{covar_en_sym},
    ::Val{var},
    ::Val{gm_third_m_sym},
) where {S, covar_en_sym, var, gm_third_m_sym}

    N_up = n_updrafts(edmf)
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
    return nothing
end
