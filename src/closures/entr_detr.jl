#### Entrainment-Detrainment kernels

function compute_turbulent_entrainment(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    c_γ = FT(CPEDMF.c_γ(param_set))

    ε_turb = if εδ_model_vars.w_up * εδ_model_vars.a_up > 0
        2 * c_γ * sqrt(max(εδ_model_vars.tke, 0)) / (εδ_model_vars.w_up * εδ_model_vars.H_up)
    else
        FT(0)
    end

    return ε_turb
end

function compute_inverse_timescale(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    Δb = εδ_model_vars.b_up - εδ_model_vars.b_en
    Δw = get_Δw(param_set, εδ_model_vars)
    c_λ = FT(CPEDMF.c_λ(param_set))

    l_1 = c_λ * abs(Δb / sqrt(εδ_model_vars.tke + 1e-8))
    l_2 = abs(Δb / Δw)
    l = SA.SVector(l_1, l_2)
    return lamb_smooth_minimum(l, FT(0.1), FT(0.0005))
end

function get_Δw(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    Δw = εδ_model_vars.w_up - εδ_model_vars.w_en
    Δw += copysign(FT(CPEDMF.w_min(param_set)), Δw)
    return Δw
end

function get_MdMdz(εδ_model_vars)
    FT = eltype(εδ_model_vars)
    MdMdz_ε = max(εδ_model_vars.dMdz / max(εδ_model_vars.M, eps(FT)), 0)
    MdMdz_δ = max(-εδ_model_vars.dMdz / max(εδ_model_vars.M, eps(FT)), 0)
    return MdMdz_ε, MdMdz_δ
end

function entrainment_length_scale(param_set, εδ_model_vars)
    Δw = get_Δw(param_set, εδ_model_vars)
    λ = compute_inverse_timescale(param_set, εδ_model_vars)
    return (λ / Δw)
end

"""
    entr_detr(param_set, εδ_model_vars, εδ_model_type)

Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, following
Cohen et al. (JAMES, 2020), given:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: type of non-dimensional model for entrainment/detrainment
"""
function entr_detr(param_set::APS, εδ_model_vars, εδ_model_type)
    FT = eltype(εδ_model_vars)
    dim_scale = entrainment_length_scale(param_set, εδ_model_vars)
    area_limiter = max_area_limiter(param_set, εδ_model_vars)

    c_div = FT(ICP.entrainment_massflux_div_factor(param_set))
    MdMdz_ε, MdMdz_δ = get_MdMdz(εδ_model_vars) .* c_div

    nondim_ε, nondim_δ = non_dimensional_function(param_set, εδ_model_vars, εδ_model_type)

    # dynamic entrainment / detrainment
    ε_dyn = dim_scale * nondim_ε + MdMdz_ε
    δ_dyn = dim_scale * (nondim_δ + area_limiter) + MdMdz_δ

    # turbulent entrainment
    ε_turb = compute_turbulent_entrainment(param_set, εδ_model_vars)

    return EntrDetr{FT}(ε_dyn, δ_dyn, ε_turb, nondim_ε, nondim_δ)
end

"""
    entr_detr_given_scales(param_set, εδ_model_vars, εδ_model_type)

Returns the dynamic entrainment and detrainment rates
given non dimenational values given by another model,
as well as the turbulent entrainment rate, following
Cohen et al. (JAMES, 2020), given:
 - `param_set`      :: parameter set
 - `dim_scale`      :: dimensional scale
 - `nondim_ε,`      :: non dimentional entrainment value
 - `nondim_δ,`      :: non dimentional detrainment value
 - `εδ_model_vars`  :: structure containing variables
"""
function entr_detr_given_scales(param_set::APS, dim_scale, nondim_ε, nondim_δ, εδ_model_vars)
    area_limiter = max_area_limiter(param_set, εδ_model_vars)

    c_div = FT(ICP.entrainment_massflux_div_factor(param_set))
    MdMdz_ε, MdMdz_δ = get_MdMdz(εδ_model_vars) .* c_div

    # dynamic entrainment / detrainment
    ε_dyn = dim_scale * nondim_ε + MdMdz_ε
    δ_dyn = dim_scale * (nondim_δ + area_limiter) + MdMdz_δ

    # turbulent entrainment
    ε_turb = compute_turbulent_entrainment(param_set, εδ_model_vars)

    return EntrDetr{FT}(ε_dyn, δ_dyn, ε_turb)
end

##### Compute entr detr

function compute_entr_detr!(
    state::State,
    grid::Grid,
    edmf::EDMF_PrognosticTKE,
    param_set::APS,
    surf::SurfaceBase,
    Δt::Real,
    ::AbstractEntrDetrModel,
)
    FT = eltype(grid)
    N_up = n_updrafts(edmf)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    w_up_c = aux_tc.w_up_c
    w_en_c = aux_tc.w_en_c
    m_entr_detr = aux_tc.ϕ_temporary
    ∇m_entr_detr = aux_tc.ψ_temporary
    wvec = CC.Geometry.WVector
    max_area = edmf.max_area
    plume_scale_height = map(1:N_up) do i
        compute_plume_scale_height(grid, state, param_set, i)
    end

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

                εδ_model_vars = GeneralizedEntr{FT}(;
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
                    H_up = plume_scale_height[i],
                    RH_up = aux_up[i].RH[k],
                    RH_en = aux_en.RH[k],
                    max_area = max_area,
                    zc_i = grid.zc[k].z,
                    Δt = Δt,
                    # non-dimensional entr/detr state
                    nondim_entr_sc = aux_up[i].nondim_entr_sc[k],
                    nondim_detr_sc = aux_up[i].nondim_detr_sc[k],
                )

                er = entr_detr(param_set, εδ_model_vars, edmf.entr_closure)
                aux_up[i].entr_sc[k] = er.ε_dyn
                aux_up[i].detr_sc[k] = er.δ_dyn
                aux_up[i].frac_turb_entr[k] = er.ε_turb
            else
                aux_up[i].entr_sc[k] = 0.0
                aux_up[i].detr_sc[k] = 0.0
                aux_up[i].frac_turb_entr[k] = 0.0
            end
        end
    end
end

function compute_entr_detr!(
    state::State,
    grid::Grid,
    edmf::EDMF_PrognosticTKE,
    param_set::APS,
    surf::SurfaceBase,
    Δt::Real,
    ::AbstractNonLocalEntrDetrModel,
)
    FT = eltype(grid)
    N_up = n_updrafts(edmf)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    w_up_c = aux_tc.w_up_c
    w_en_c = aux_tc.w_en_c
    m_entr_detr = aux_tc.ϕ_temporary
    ∇m_entr_detr = aux_tc.ψ_temporary
    wvec = CC.Geometry.WVector
    max_area = edmf.max_area
    plume_scale_height = map(1:N_up) do i
        compute_plume_scale_height(grid, state, param_set, i)
    end

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
                εδ_model_vars = GeneralizedEntr{FT}(;
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
                    H_up = plume_scale_height[i],
                    RH_up = aux_up[i].RH[k],
                    RH_en = aux_en.RH[k],
                    max_area = max_area,
                    zc_i = grid.zc[k].z,
                    Δt = Δt,
                    # non-dimensional entr/detr state
                    nondim_entr_sc = aux_up[i].nondim_entr_sc[k],
                    nondim_detr_sc = aux_up[i].nondim_detr_sc[k],
                )
                aux_up[i].ε_dim[k] = entrainment_length_scale(param_set, εδ_model_vars)
                Π = non_dimensional_groups(param_set, εδ_model_vars)
                aux_up[i].Π₁[k] = Π[1]
                aux_up[i].Π₂[k] = Π[2]
                aux_up[i].Π₃[k] = Π[3]
                aux_up[i].Π₄[k] = Π[4]
            else
                # TODO: is there a better way to do this?
                aux_up[i].Π₁[k] = 0
                aux_up[i].Π₂[k] = 0
                aux_up[i].Π₃[k] = 0
                aux_up[i].Π₄[k] = 0
            end
        end
        ε_dim = parent(aux_up[i].ε_dim)
        Π₁ = parent(aux_up[i].Π₁)
        Π₂ = parent(aux_up[i].Π₂)
        Π₃ = parent(aux_up[i].Π₃)
        Π₄ = parent(aux_up[i].Π₄)
        nondim_ε, nondim_δ = non_dimensional_function(param_set, Π₁, Π₂, Π₃, Π₄. εδ_model_type)
        ε_dyn, δ_dyn, ε_turb = entr_detr_given_scales(APS, ε_dim, nondim_ε, nondim_δ, εδ_model_vars)
        aux_up[i].entr_sc = entr_sc
        aux_up[i].detr_sc = detr_sc
        aux_up[i].frac_turb_entr = frac_turb_entr
    end
end
