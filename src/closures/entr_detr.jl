#### Entrainment-Detrainment kernels

"""
    entr_detr(
        param_set,
        w_min,
        β,
        c_δ,
        c_μ,
        c_μ0,
        χ_upd,
        c_λ,
        εδ_model,
    )

Returns the dynamic entrainment and detrainment rates,
as well as the turbulent entrainment rate, following
Cohen et al. (JAMES, 2020), given:
 - `param_set`: parameter set
 - `w_min`: minimum veritcal velocity
 - `β`: sorting power for moist mixing
 - `c_δ`: detrainment factor
 - `c_μ`: logisitc function scale
 - `c_μ0`: logisitc function timescale
 - `χ_upd`: updraft mixing fraction
 - `c_λ`: tke scale factor
 - `c_εt`: turbulent entrainment factor
 - `εδ_model`: a [`MoistureDeficitEntr`](@ref)
"""
function entr_detr(param_set, w_min, β, c_δ, c_μ, c_μ0, χ_upd, c_λ, c_εt, εδ_model::MoistureDeficitEntr)

    c_ε = CPEDMF.c_ε(param_set)
    c_div = ICP.entrainment_massflux_div_factor(param_set)

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

    D_ε, D_δ, M_δ, M_ε = nondimensional_exchange_functions(param_set, c_δ, c_μ, c_μ0, β, χ_upd, Δw, Δb, εδ_model)

    l_1 = c_λ * abs(Δb / sqrt(εδ_model.tke + 1e-8))
    l_2 = abs(Δb / Δw)
    l = (l_1, l_2)
    λ = lamb_smooth_minimum(l, 0.1, 0.0005)

    MdMdz = max(εδ_model.dMdz / max(εδ_model.M, 1e-12), 0.0)
    MdMdz = max(-εδ_model.dMdz / max(εδ_model.M, 1e-12), 0.0)



    # turbulent entrainment
    K_ε = εδ_model.a_up * c_εt * sqrt(max(εδ_model.tke, 0.0)) * εδ_model.R_up
    if εδ_model.w_up * εδ_model.a_up > 0.0
        ε_turb = (2.0 / εδ_model.R_up^2.0) * K_ε / (εδ_model.w_up * εδ_model.a_up)
    else
        ε_turb = 0.0
    end

    ε_dyn = λ / Δw * (c_ε * D_ε + c_δ * M_ε) + MdMdz * c_div
    δ_dyn = λ / Δw * (c_ε * D_δ + c_δ * M_δ) + MdMdz * c_div


    return EntrDetr(ε_dyn, δ_dyn, ε_turb, K_ε)
end
