#### Entrainment-Detrainment kernels

"""
    entr_detr(
        param_set,
        εδ_model,
    )
Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, following
Cohen et al. (JAMES, 2020), given:
 - `param_set`: parameter set
 - `εδ_model`: a [`MoistureDeficitEntr`](@ref)
"""
function entr_detr(param_set, εδ_model::MoistureDeficitEntr)

    l = zeros(2)
    c_ε = CPEDMF.c_ε(param_set)
    c_δ = CPEDMF.c_δ(param_set)
    c_t = CPEDMF.c_t(param_set)
    c_λ = CPEDMF.c_λ(param_set)
    w_min = CPEDMF.w_min(param_set)
    c_div = ICP.entrainment_massflux_div_factor(param_set)

    # should be: c_δ = sign(condensate(ts_en) + condensate(ts_up[i])) * entr.c_δ
    if (εδ_model.q_liq_up + εδ_model.q_liq_en) == 0.0
        c_δ = 0.0
    end

    Δw = εδ_model.w_up - εδ_model.w_en
    if Δw < 0.0
        Δw -= w_min
    else
        Δw += w_min
    end

    Δb = (εδ_model.b_up - εδ_model.b_en)

    D_ε, D_δ, M_δ, M_ε = nondimensional_exchange_functions(param_set, Δw, Δb, εδ_model)

    l_1 = c_λ * abs(Δb / sqrt(εδ_model.tke + 1e-8))
    l_2 = abs(Δb / Δw)
    l = (l_1, l_2)
    λ = lamb_smooth_minimum(l, 0.1, 0.0005)

    MdMdz = max(εδ_model.dMdz / max(εδ_model.M, 1e-12), 0.0)
    MdMdz = max(-εδ_model.dMdz / max(εδ_model.M, 1e-12), 0.0)

    # turbulent entrainment
    K_ε = εδ_model.a_up * c_t * sqrt(max(εδ_model.tke, 0.0)) * εδ_model.R_up
    if εδ_model.w_up * εδ_model.a_up > 0.0
        ε_turb = (2.0 / εδ_model.R_up^2.0) * K_ε / (εδ_model.w_up * εδ_model.a_up)
    else
        ε_turb = 0.0
    end

    ε_dyn = λ / Δw * (c_ε * D_ε + c_δ * M_ε) + MdMdz * c_div
    δ_dyn = λ / Δw * (c_ε * D_δ + c_δ * M_δ) + MdMdz * c_div

    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end
