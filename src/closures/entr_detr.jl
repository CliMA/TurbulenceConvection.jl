#### Entrainment-Detrainment kernels

function compute_turbulent_entrainment(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    c_γ = FT(ICP.turbulent_entrainment_factor(param_set))

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

    return EntrDetr{FT}(ε_dyn, δ_dyn, ε_turb)
end
