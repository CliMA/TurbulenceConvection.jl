#### Entrainment-Detrainment kernels
"""
    entr_detr(
        param_set,
        w_min,
        β,
        c_δ,
        c_div,
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
 - `c_div`: divergence factor for bubble case (zero otherwise)
 - `c_μ`: logisitc function scale 
 - `c_μ0`: logisitc function timescale
 - `χ_upd`: updraft mixing fraction
 - `c_λ`: tke scale factor
 - `c_εt`: turbulent entrainment factor
 - `εδ_model`: entrainment detrainment model type
"""
function entr_detr(param_set, w_min, β, c_δ, c_div, c_μ, c_μ0, χ_upd, c_λ, c_εt, εδ_model)

    l = zeros(2)
    c_ε = CPEDMF.c_ε(param_set)

    if (εδ_model.ql_upd + εδ_model.ql_env) == 0.0
        c_δ = 0.0
    end

    Δw = εδ_model.w_upd - εδ_model.w_env
    if Δw < 0.0
        Δw -= w_min
    else
        Δw += w_min
    end

    Δb = (εδ_model.b_upd - εδ_model.b_env)

    D_ε, D_δ, M_δ, M_ε = nondimensional_exchange_functions(c_ε, c_δ, c_μ, c_μ0, β, χ_upd, Δw, Δb, εδ_model)

    l[1] = c_λ * fabs(Δb / sqrt(εδ_model.tke + 1e-8))
    l[2] = fabs(Δb / Δw)
    λ = lamb_smooth_minimum(l, 0.1, 0.0005)

    MdMdz = fmax(εδ_model.dMdz / fmax(εδ_model.M, 1e-12), 0.0)
    MdMdz = fmax(-εδ_model.dMdz / fmax(εδ_model.M, 1e-12), 0.0)



    # turbulent entrainment
    k_ε = εδ_model.a_upd * c_εt * sqrt(fmax(εδ_model.tke, 0.0)) * εδ_model.R_up
    if εδ_model.w_upd * εδ_model.a_upd > 0.0
        ε_turb = (2.0 / εδ_model.R_up^2.0) * k_ε / (εδ_model.w_upd * εδ_model.a_upd)
    else
        ε_turb = 0.0
    end

    ε_dyn = λ / Δw * (c_ε * D_ε + c_δ * M_ε) + MdMdz * c_div
    δ_dyn = λ / Δw * (c_ε * D_δ + c_δ * M_δ) + MdMdz * c_div

    return ε_dyn, δ_dyn, ε_turb, k_ε
end
