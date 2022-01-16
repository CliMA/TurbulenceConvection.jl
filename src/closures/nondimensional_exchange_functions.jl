#### Non-dimensional Entrainment-Detrainment functions

function max_area_limiter(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    γ_lim = FT(ICP.area_limiter_scale(param_set))
    β_lim = FT(ICP.area_limiter_power(param_set))
    logistic_term = (2 - 1 / (1 + exp(-γ_lim * (εδ_model_vars.max_area - εδ_model_vars.a_up))))
    return logistic_term^β_lim - 1
end

function non_dimensional_groups(param_set, εδ_model_vars)
    Δw = get_Δw(param_set, εδ_model_vars)
    Δb = εδ_model_vars.b_up - εδ_model_vars.b_en
    Π_1 = Δw^2 / (εδ_model_vars.zc_i * Δb)
    Π_2 = Δw^2 / εδ_model_vars.tke
    Π_3 = √(εδ_model_vars.a_up)
    Π_4 = εδ_model_vars.RH_up - εδ_model_vars.RH_en
    return SA.SVector(Π_1, Π_2, Π_3, Π_4)
end

"""
    non_dimensional_function(param_set, εδ_model_vars, ::MDEntr)

Returns the nondimensional entrainment and detrainment
functions following Cohen et al. (JAMES, 2020), given:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: MDEntr - Moisture deficit entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, ::MDEntr)
    FT = eltype(εδ_model_vars)
    Δw = get_Δw(param_set, εδ_model_vars)
    c_ε = FT(CPEDMF.c_ε(param_set))
    μ_0 = FT(CPEDMF.μ_0(param_set))
    β = FT(CPEDMF.β(param_set))
    χ = FT(CPEDMF.χ(param_set))
    c_δ = if !TD.has_condensate(εδ_model_vars.q_cond_up + εδ_model_vars.q_cond_en)
        FT(0)
    else
        FT(CPEDMF.c_δ(param_set))
    end

    Δb = εδ_model_vars.b_up - εδ_model_vars.b_en
    μ_ij = (χ - εδ_model_vars.a_up / (εδ_model_vars.a_up + εδ_model_vars.a_en)) * Δb / Δw
    exp_arg = μ_ij / μ_0
    D_ε = 1 / (1 + exp(-exp_arg))
    D_δ = 1 / (1 + exp(exp_arg))

    M_δ = (max((εδ_model_vars.RH_up)^β - (εδ_model_vars.RH_en)^β, 0))^(1 / β)
    M_ε = (max((εδ_model_vars.RH_en)^β - (εδ_model_vars.RH_up)^β, 0))^(1 / β)

    nondim_ε = (c_ε * D_ε + c_δ * M_ε)
    nondim_δ = (c_ε * D_δ + c_δ * M_δ)
    return nondim_ε, nondim_δ
end
