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

"""
    non_dimensional_function(param_set, εδ_model_vars, ::NNEntr)

Uses a fully connected neural network to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NNEntr - Neural network entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, ::NNEntr)
    c_gen = ICP.c_gen(param_set)

    # Neural network closure
    nn_arc = (4, 2, 2)  # (#inputs, #neurons, #outputs)
    nn_model = Flux.Chain(
        Flux.Dense(reshape(c_gen[1:8], nn_arc[2], nn_arc[1]), c_gen[9:10], Flux.sigmoid),
        Flux.Dense(reshape(c_gen[11:14], nn_arc[3], nn_arc[2]), c_gen[15:16], Flux.softplus),
    )

    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars)
    nondim_ε, nondim_δ = nn_model(nondim_groups)
    return nondim_ε, nondim_δ
end

"""
    non_dimensional_function(param_set, εδ_model_vars, ::LinearEntr)

Uses a simple linear model to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: LinearEntr - linear entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, ::LinearEntr)
    c_gen = ICP.c_gen(param_set)

    # Linear closure
    lin_arc = (4, 1)  # (#weights, #outputs)
    lin_model_ε = Flux.Dense(reshape(c_gen[1:4], lin_arc[2], lin_arc[1]), [c_gen[5]], Flux.relu)
    lin_model_δ = Flux.Dense(reshape(c_gen[6:9], lin_arc[2], lin_arc[1]), [c_gen[10]], Flux.relu)

    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars)
    nondim_ε = lin_model_ε(nondim_groups)[1]
    nondim_δ = lin_model_δ(nondim_groups)[1]
    return nondim_ε, nondim_δ
end

"""
    non_dimensional_function(param_set, εδ_model_vars, εδ_model_type::LogNormalScalingProcess)

Uses a LogNormal random variable to scale a deterministic process
to predict the non-dimensional components of dynamical entrainment/detrainment.

Arguments:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: LogNormalScalingProcess - Stochastic lognormal scaling
"""
function non_dimensional_function(param_set, εδ_model_vars, εδ_model_type::LogNormalScalingProcess)
    FT = eltype(εδ_model_vars)
    # model parameters
    mean_model = εδ_model_type.mean_model
    c_gen_stoch = ICP.c_gen_stoch(param_set)
    ε_σ² = c_gen_stoch[1]
    δ_σ² = c_gen_stoch[2]

    # Mean model closure
    ε_mean_nondim, δ_mean_nondim = non_dimensional_function(param_set, εδ_model_vars, mean_model)

    # lognormal scaling
    nondim_ε = ε_mean_nondim * lognormal_sampler(FT(1), ε_σ²)
    nondim_δ = δ_mean_nondim * lognormal_sampler(FT(1), δ_σ²)

    return nondim_ε, nondim_δ
end

function lognormal_sampler(m::FT, var::FT)::FT where {FT}
    μ = log(m^2 / √(m^2 + var))
    σ = √(log(1 + var / m^2))
    return rand(Distributions.LogNormal(μ, σ))
end

"""
    non_dimensional_function(param_set, εδ_model_vars, εδ_model_type::NoisyRelaxationProcess)

Uses a noisy relaxation process to predict the non-dimensional components 
of dynamical entrainment/detrainment. A deterministic closure is used as the
equilibrium mean function for the relaxation process.

Arguments:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NoisyRelaxationProcess - A noisy relaxation process closure
"""
function non_dimensional_function(param_set, εδ_model_vars, εδ_model_type::NoisyRelaxationProcess)
    # model parameters
    mean_model = εδ_model_type.mean_model
    c_gen_stoch = ICP.c_gen_stoch(param_set)
    ε_σ² = c_gen_stoch[1]
    δ_σ² = c_gen_stoch[2]
    ε_λ = c_gen_stoch[3]
    δ_λ = c_gen_stoch[4]

    # Mean model closure
    ε_mean_nondim, δ_mean_nondim = non_dimensional_function(param_set, εδ_model_vars, mean_model)

    # noisy relaxation process
    ε_u0 = εδ_model_vars.nondim_entr_sc
    δ_u0 = εδ_model_vars.nondim_detr_sc
    Δt = εδ_model_vars.Δt
    nondim_ε = noisy_relaxation_process(ε_mean_nondim, ε_λ, ε_σ², ε_u0, Δt)
    nondim_δ = noisy_relaxation_process(δ_mean_nondim, δ_λ, δ_σ², δ_u0, Δt)

    return nondim_ε, nondim_δ
end

""" 
    Solve a noisy relaxation process numerically 

In this formulation, the noise amplitude is scaled by the speed of 
reversion λ and the long-term mean μ, in addition to the variance σ² as is usual,

    `du = λ(μ - u)⋅dt + √(2λμσ²)⋅dW`

To ensure non-negativity, the solution is passed through a relu filter.
"""
function noisy_relaxation_process(μ::FT, λ::FT, σ²::FT, u0::FT, Δt::FT)::FT where {FT}
    f(u, p, t) = λ * (μ - u)        # mean-reverting process
    g(u, p, t) = √(2λ * μ * σ²)     # noise fluctuation
    tspan = (0.0, Δt)
    prob = SDE.SDEProblem(f, g, u0, tspan; save_start = false, saveat = last(tspan))
    sol = SDE.solve(prob, SDE.SOSRI())
    return Flux.relu(sol.u[end])
end
