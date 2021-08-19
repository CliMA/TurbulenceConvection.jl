"""
    nondimensional_exchange_functions(
        Δw,
        Δb,
        εδ_model,
    )

Returns the nondimensional entrainment and detrainment
functions following Cohen et al. (JAMES, 2020), given:
 - `Δw`, updraft - environment vertical velocity differnce
 - `Δb`, updraft - environment buoynacy differnce
 - `εδ_model`, entrainment detrainment model type
"""
function nondimensional_exchange_functions(param_set, Δw, Δb, εδ_model)

    μ_0 = CPEDMF.μ_0(param_set)
    μ = ICP.entrainment_sigma(param_set)
    β = CPEDMF.β(param_set)
    χ = CPEDMF.χ(param_set)


    μ_ij = (χ - εδ_model.a_up / (εδ_model.a_up + εδ_model.a_en)) * Δb / Δw
    D_ε = 1.0 / (1.0 + exp(-μ / μ_0 * μ_ij))
    D_δ = 1.0 / (1.0 + exp(μ / μ_0 * μ_ij))

    M_δ = (max((εδ_model.RH_up)^β - (εδ_model.RH_en)^β, 0.0))^(1.0 / β)
    M_ε = (max((εδ_model.RH_en)^β - (εδ_model.RH_up)^β, 0.0))^(1.0 / β)

    return D_ε, D_δ, M_δ, M_ε
end;
