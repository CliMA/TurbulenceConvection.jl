#### Entrainment-Detrainment kernels
using Flux

function compute_turbulent_entrainment(param_set, εδ_model_vars)

    c_t = CPEDMF.c_t(param_set)
    K_ε = εδ_model_vars.a_up * c_t * sqrt(max(εδ_model_vars.tke, 0.0)) * εδ_model_vars.R_up
    if εδ_model_vars.w_up * εδ_model_vars.a_up > 0.0
        ε_turb = (2.0 / εδ_model_vars.R_up^2.0) * K_ε / (εδ_model_vars.w_up * εδ_model_vars.a_up)
    else
        ε_turb = 0.0
    end

    return (ε_turb, K_ε)
end

function compute_inverse_timescale(Δb, Δw, param_set, εδ_model_vars)

    c_λ = CPEDMF.c_λ(param_set)

    l_1 = c_λ * abs(Δb / sqrt(εδ_model_vars.tke + 1e-8))
    l_2 = abs(Δb / Δw)
    l = (l_1, l_2)
    return lamb_smooth_minimum(l, 0.1, 0.0005)

end

function preprocess_inputs(param_set, εδ_model_vars, Δw)

    w_min = CPEDMF.w_min(param_set)
    c_δ = CPEDMF.c_δ(param_set)

    # should be: c_δ = sign(condensate(ts_en) + condensate(ts_up[i])) * entr.c_δ
    if !TD.has_condensate(εδ_model_vars.q_cond_up + εδ_model_vars.q_cond_en)
        c_δ = 0.0
    end

    # clip vertical velocity
    if Δw < 0.0
        Δw -= w_min
    else
        Δw += w_min
    end

    MdMdz = max(εδ_model_vars.dMdz / max(εδ_model_vars.M, 1e-12), 0.0)
    MdMdz = max(-εδ_model_vars.dMdz / max(εδ_model_vars.M, 1e-12), 0.0)

    return (Δw, MdMdz, c_δ)
end

"""
    entr_detr(
        param_set,
        εδ_model_vars,
        εδ_model_type:MDEntr,
    )

Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, following
Cohen et al. (JAMES, 2020), given:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: MDEntr - Moisture deficit entrainment closure
"""
function entr_detr(param_set, εδ_model_vars, εδ_model_type::MDEntr)

    γ_lim = ICP.area_limiter_scale(param_set)
    β_lim = ICP.area_limiter_power(param_set)
    c_ε = CPEDMF.c_ε(param_set)
    c_t = CPEDMF.c_t(param_set)
    c_div = ICP.entrainment_massflux_div_factor(param_set)

    Δw = εδ_model_vars.w_up - εδ_model_vars.w_en
    Δb = εδ_model_vars.b_up - εδ_model_vars.b_en

    Δw, MdMdz, c_δ = preprocess_inputs(param_set, εδ_model_vars, Δw)

    D_ε, D_δ, M_δ, M_ε = nondimensional_exchange_functions(param_set, Δw, Δb, εδ_model_vars)

    λ = compute_inverse_timescale(Δb, Δw, param_set, εδ_model_vars)

    # turbulent entrainment
    ε_turb, K_ε = compute_turbulent_entrainment(param_set, εδ_model_vars)

    ε_dyn = λ / Δw * (c_ε * D_ε + c_δ * M_ε) + MdMdz * c_div
    logistic_term = (2.0 - 1.0 / (1 + exp(-γ_lim * (εδ_model_vars.max_area - εδ_model_vars.a_up))))
    max_area_limiter = λ / Δw * (logistic_term^β_lim - 1.0)
    δ_dyn = (λ / Δw * (c_ε * D_δ + c_δ * M_δ) + MdMdz * c_div) + max_area_limiter

    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end

"""
    entr_detr(
        param_set,
        εδ_model_vars,
        εδ_model_type::NNEntr,
    )

Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, using a neural
network to predict the non-dimensional component of dynamical entrainment:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NNEntr - Neural network entrainment closure
"""
function entr_detr(param_set, εδ_model_vars, εδ_model_type::NNEntr)

    γ_lim = ICP.area_limiter_scale(param_set)
    β_lim = ICP.area_limiter_power(param_set)
    c_ε = CPEDMF.c_ε(param_set)
    c_t = CPEDMF.c_t(param_set)
    c_gen = ICP.c_gen(param_set)
    c_div = ICP.entrainment_massflux_div_factor(param_set)

    Δw = εδ_model_vars.w_up - εδ_model_vars.w_en
    Δb = εδ_model_vars.b_up - εδ_model_vars.b_en

    Δw, MdMdz, c_δ = preprocess_inputs(param_set, εδ_model_vars, Δw)

    λ = compute_inverse_timescale(Δb, Δw, param_set, εδ_model_vars)

    # turbulent entrainment
    ε_turb, K_ε = compute_turbulent_entrainment(param_set, εδ_model_vars)

    logistic_term = (2.0 - 1.0 / (1 + exp(-γ_lim * (εδ_model_vars.max_area - εδ_model_vars.a_up))))
    max_area_limiter = λ / Δw * (logistic_term^β_lim - 1.0)

    nn_arc = (4, 2, 2) # number of inputs, number of neurons, number of outputs
    nn_model = Chain(
        Dense(reshape(c_gen[1:8], nn_arc[2], nn_arc[1]), c_gen[9:10], sigmoid),
        Dense(reshape(c_gen[11:14], nn_arc[3], nn_arc[2]), c_gen[15:16], sigmoid),
    )

    pi_1 = εδ_model_vars.updraft_top * Δb / Δw^2
    pi_2 = εδ_model_vars.tke / Δw^2
    non_dim_functions = nn_model([pi_1, pi_2, εδ_model_vars.a_up, εδ_model_vars.RH_up - εδ_model_vars.RH_en])
    ε_dyn = (λ / Δw) * non_dim_functions[1]
    δ_dyn = (λ / Δw) * non_dim_functions[2] + max_area_limiter

    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end
