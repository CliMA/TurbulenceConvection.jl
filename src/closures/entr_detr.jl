#### Entrainment-Detrainment kernels

function compute_turbulent_entrainment(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    c_t = FT(CPEDMF.c_t(param_set))
    K_ε = εδ_model_vars.a_up * c_t * sqrt(max(εδ_model_vars.tke, 0)) * εδ_model_vars.R_up

    ε_turb = if εδ_model_vars.w_up * εδ_model_vars.a_up > 0
        (2 / εδ_model_vars.R_up^2) * K_ε / (εδ_model_vars.w_up * εδ_model_vars.a_up)
    else
        FT(0)
    end

    return (ε_turb, K_ε)
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

function get_c_δ(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    c_δ = if !TD.has_condensate(εδ_model_vars.q_cond_up + εδ_model_vars.q_cond_en)
        FT(0)
    else
        FT(CPEDMF.c_δ(param_set))
    end
    return c_δ
end

function get_MdMdz(εδ_model_vars)
    FT = eltype(εδ_model_vars)
    MdMdz_ε = max(εδ_model_vars.dMdz / max(εδ_model_vars.M, eps(FT)), 0)
    MdMdz_δ = max(-εδ_model_vars.dMdz / max(εδ_model_vars.M, eps(FT)), 0)
    return MdMdz_ε, MdMdz_δ
end

function dimensional_part(param_set, εδ_model_vars)
    Δw = get_Δw(param_set, εδ_model_vars)
    λ = compute_inverse_timescale(param_set, εδ_model_vars)
    return (λ / Δw)
end

function max_area_limiter(param_set, εδ_model_vars)
    FT = eltype(εδ_model_vars)
    γ_lim = FT(ICP.area_limiter_scale(param_set))
    β_lim = FT(ICP.area_limiter_power(param_set))
    logistic_term = (2 - 1 / (1 + exp(-γ_lim * (εδ_model_vars.max_area - εδ_model_vars.a_up))))
    area_limiter = dimensional_part(param_set, εδ_model_vars) * (logistic_term^β_lim - 1)
    return area_limiter
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
    entr_detr(param_set, εδ_model_vars, εδ_model_type::MDEntr)

Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, following
Cohen et al. (JAMES, 2020), given:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: MDEntr - Moisture deficit entrainment closure
"""
function entr_detr(param_set, εδ_model_vars, εδ_model_type::MDEntr)
    FT = eltype(εδ_model_vars)
    dim_scale = dimensional_part(param_set, εδ_model_vars)
    area_limiter = max_area_limiter(param_set, εδ_model_vars)

    # Moisture deficit closure
    Δw = get_Δw(param_set, εδ_model_vars)
    D_ε, D_δ, M_δ, M_ε = nondimensional_exchange_functions(param_set, Δw, εδ_model_vars)
    c_ε = FT(CPEDMF.c_ε(param_set))
    c_δ = get_c_δ(param_set, εδ_model_vars)
    nondim_ε = (c_ε * D_ε + c_δ * M_ε)
    nondim_δ = (c_ε * D_δ + c_δ * M_δ)
    c_div = FT(ICP.entrainment_massflux_div_factor(param_set))
    MdMdz_ε, MdMdz_δ = get_MdMdz(εδ_model_vars) .* c_div

    # dynamic entrainment / detrainment
    ε_dyn = dim_scale * nondim_ε
    δ_dyn = dim_scale * nondim_δ + area_limiter

    # turbulent entrainment
    ε_turb, K_ε = compute_turbulent_entrainment(param_set, εδ_model_vars)

    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end

"""
    entr_detr(param_set, εδ_model_vars, εδ_model_type::NNEntr)

Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, using a neural
network to predict the non-dimensional component of dynamical entrainment:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NNEntr - Neural network entrainment closure
"""
function entr_detr(param_set, εδ_model_vars, εδ_model_type::NNEntr)
    c_gen = ICP.c_gen(param_set)
    dim_scale = dimensional_part(param_set, εδ_model_vars)
    area_limiter = max_area_limiter(param_set, εδ_model_vars)

    # Neural network closure
    nn_arc = (4, 2, 2)  # (#inputs, #neurons, #outputs)
    nn_model = Flux.Chain(
        Flux.Dense(reshape(c_gen[1:8], nn_arc[2], nn_arc[1]), c_gen[9:10], sigmoid),
        Flux.Dense(reshape(c_gen[11:14], nn_arc[3], nn_arc[2]), c_gen[15:16], softplus),
    )

    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars)
    nondim_ε, nondim_δ = nn_model(nondim_groups)

    # dynamic entrainment / detrainment
    ε_dyn = dim_scale * nondim_ε
    δ_dyn = dim_scale * nondim_δ + area_limiter

    # turbulent entrainment
    ε_turb, K_ε = compute_turbulent_entrainment(param_set, εδ_model_vars)

    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end
