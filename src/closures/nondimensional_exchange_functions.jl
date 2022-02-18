#### Non-dimensional Entrainment-Detrainment functions
function max_area_limiter(param_set, max_area, a_up)
    FT = eltype(a_up)
    γ_lim = FT(ICP.area_limiter_scale(param_set))
    β_lim = FT(ICP.area_limiter_power(param_set))
    logistic_term = (2 - 1 / (1 + exp(-γ_lim * (max_area - a_up))))
    return (logistic_term)^β_lim - 1
end

function non_dimensional_groups(param_set, εδ_model_vars, entr_pi_groups)
    Δw = get_Δw(param_set, εδ_model_vars.w_up, εδ_model_vars.w_en)
    Δb = εδ_model_vars.b_up - εδ_model_vars.b_en
    FT = eltype(εδ_model_vars.tke_en)
    Π₁ = (εδ_model_vars.zc_i * Δb) / (Δw^2 + εδ_model_vars.wstar^2)
    Π₁ = Π₁ / FT(478.298) # Normalize by max(abs(Π₁)) found in {Bomex, DYCOMS, TRMM}
    Π₂ = (εδ_model_vars.tke_gm - εδ_model_vars.a_en * εδ_model_vars.tke_en) / (εδ_model_vars.tke_gm + eps(FT))
    Π₃ = √(εδ_model_vars.a_up)
    Π₄ = εδ_model_vars.RH_up - εδ_model_vars.RH_en
    Π₅ = εδ_model_vars.zc_i / εδ_model_vars.H_up
    Π₆ = εδ_model_vars.zc_i / εδ_model_vars.ref_H_up
    # pi_groups = (Π₁, Π₂, Π₃, Π₄, Π₅, Π₆)[entr_pi_groups] # needed to avoid subarray
    return SA.SVector(Π₁, Π₂, Π₃, Π₄, Π₅, Π₆)[entr_pi_groups]
end

"""
    non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::MDEntr)

Returns the nondimensional entrainment and detrainment
functions following Cohen et al. (JAMES, 2020), given:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: MDEntr - Moisture deficit entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::MDEntr)
    FT = eltype(εδ_model_vars)
    Δw = get_Δw(param_set, εδ_model_vars.w_up, εδ_model_vars.w_en)
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
    non_dimensional_function!(nondim_ε ,nondim_δ ,param_set ,Π₁ ,Π₂ ,Π₃ ,Π₄, εδ_model::FNOEntr)

Uses a non local (Fourier) neural network to predict the fields of
    non-dimensional components of dynamical entrainment/detrainment.
 - `nondim_ε`   :: output - non dimensional entr from FNO, as column fields
 - `nondim_δ`   :: output - non dimensional detr from FNO, as column fields
 - `param_set`  :: input - parameter set
 - `Π₁,₂,₃,₄`   :: input - non dimensional groups, as column fields
 - `::FNOEntr ` a non-local entrainment-detrainment model type
"""
function non_dimensional_function!(
    nondim_ε::AbstractArray{FT}, # output
    nondim_δ::AbstractArray{FT}, # output
    param_set::APS,
    Π₁::AbstractArray{FT}, # input
    Π₂::AbstractArray{FT}, # input
    Π₃::AbstractArray{FT}, # input
    Π₄::AbstractArray{FT}, # input
    Π₅::AbstractArray{FT}, # input
    Π₆::AbstractArray{FT}, # input
    entr_pi_groups,
    εδ_model::FNOEntr,
) where {FT <: Real}

    # define the model
    Π = hcat(Π₁, Π₂, Π₃, Π₄, Π₅, Π₆)'
    Π = reshape(Π, (size(Π)..., 1))[entr_pi_groups,:,:]
    trafo = OF.FourierTransform(modes = (2,))
    model = OF.Chain(OF.SpectralKernelOperator(trafo, 4 => 2, Flux.relu),)

    # set the parameters
    c_fno = ICP.c_fno(param_set)
    index = 1
    for p in Flux.params(model)
        len_p = length(p)
        p_slice = p[:]
        if eltype(p_slice) <: Real
            p_slice .= c_fno[index:(index + len_p - 1)]
            index += len_p
        elseif eltype(p_slice) <: Complex
            c_fno_slice = c_fno[index:(index + len_p - 1)]
            p_slice .= c_fno_slice + c_fno[(index + len_p):(index + len_p * 2 - 1)] * im
            index += len_p * 2
        else
            error("Bad eltype in Flux params")
        end
    end

    output = model(Π)

    # we need a sigmoid that recieves output and produced non negative nondim_ε, nondim_δ
    nondim_ε .= output[1, :]
    nondim_δ .= output[2, :]
    return nothing
end

"""
    Count number of parameters in fully-connected NN model given tuple specifying architecture following
         the pattern: (#inputs, #neurons in L1, #neurons in L2, ...., #outputs)
"""
num_params_from_arc(nn_arc::NTuple) = sum([(nn_arc[i] * nn_arc[i + 1] + nn_arc[i + 1]) for i in 1:(length(nn_arc) - 1)])

"""

    Given network architecture and parameter vector, construct NN model and unpack parameters as weights and biases.
    - `nn_arc` :: tuple specifying network architecture
    - `nn_params` :: parameter vector containing weights and biases
    - `activation_function` :: activation function for hidden layers
    - `output_layer_activation_function` :: activation function for output layer
"""
function construct_fully_connected_nn(
    nn_arc::NTuple,
    parameters::AbstractArray{FT};
    activation_function::Flux.Function = Flux.sigmoid,
    output_layer_activation_function::Flux.Function = Flux.relu,
) where {FT <: Real}

    n_params_nn = num_params_from_arc(nn_arc)
    n_params_vect = length(parameters)
    if n_params_vect != n_params_nn
        error("Incorrect number of parameters ($n_params_vect) for requested NN architecture ($n_params_nn)!")
    end

    layers = []
    parameters_i = 1
    # unpack parameters in parameter vector into network
    for layer_i in 1:(length(nn_arc) - 1)
        if layer_i == length(nn_arc) - 1
            activation_function = output_layer_activation_function
        end
        layer_num_weights = nn_arc[layer_i] * nn_arc[layer_i + 1]
        layer = Flux.Dense(
            reshape(
                parameters[parameters_i:(parameters_i + layer_num_weights - 1)],
                nn_arc[layer_i + 1],
                nn_arc[layer_i],
            ),
            parameters[(parameters_i + layer_num_weights):(parameters_i + layer_num_weights + nn_arc[layer_i + 1] - 1)],
            activation_function,
        )
        parameters_i += layer_num_weights + nn_arc[layer_i + 1]
        push!(layers, layer)
    end

    return Flux.Chain(layers...)
end

"""
    non_dimensional_function!(nondim_ε ,nondim_δ ,param_set ,Π₁ ,Π₂ ,Π₃ ,Π₄, εδ_model::NNEntrNonlocal)

Uses a fully connected neural network to predict the non-dimensional components of dynamical entrainment/detrainment.
    non-dimensional components of dynamical entrainment/detrainment.
 - `nondim_ε`   :: output - non dimensional entr from FNO, as column fields
 - `nondim_δ`   :: output - non dimensional detr from FNO, as column fields
 - `param_set`  :: input - parameter set
 - `Π₁,₂,₃,₄`   :: input - non dimensional groups, as column fields
 - `::NNEntrNonlocal ` a non-local entrainment-detrainment model type
"""
function non_dimensional_function!(
    nondim_ε::AbstractArray{FT}, # output
    nondim_δ::AbstractArray{FT}, # output
    param_set::APS,
    Π₁::AbstractArray{FT}, # input
    Π₂::AbstractArray{FT}, # input
    Π₃::AbstractArray{FT}, # input
    Π₄::AbstractArray{FT}, # input
    Π₅::AbstractArray{FT}, # input
    Π₆::AbstractArray{FT}, # input
    entr_pi_groups,
    εδ_model::NNEntrNonlocal,
) where {FT <: Real}
    # neural network architecture
    nn_arc = (length(entr_pi_groups), 5, 4, 2) # (#inputs, #neurons in L1, #neurons in L2, ...., #outputs)
    c_gen = ICP.c_gen(param_set)
    Π = hcat(Π₁, Π₂, Π₃, Π₄, Π₅, Π₆)'
    nn_model = construct_fully_connected_nn(nn_arc, c_gen)
    output = nn_model(Π)

    nondim_ε .= output[1, :]
    nondim_δ .= output[2, :]

    return nothing
end

"""
    non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::NNEntr)

Uses a fully connected neural network to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NNEntr - Neural network entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::NNEntr)
    c_gen = ICP.c_gen(param_set)

    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars, entr_pi_groups)
    # neural network architecture
    nn_arc = (length(nondim_groups), 6, 4, 2) # (#inputs, #neurons in L1, #neurons in L2, ...., #outputs)
    nn_model = construct_fully_connected_nn(nn_arc, c_gen)

    nondim_ε, nondim_δ = nn_model(nondim_groups)
    return nondim_ε, nondim_δ
end

"""
    non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::LinearEntr)

Uses a simple linear model to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: LinearEntr - linear entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::LinearEntr)
    c_gen = ICP.c_gen(param_set)

    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars, entr_pi_groups)
    # Linear closure
    lin_arc = (length(nondim_groups), 1)  # (#weights, #outputs)
    lin_model_ε = Flux.Dense(reshape(c_gen[1:4], lin_arc[2], lin_arc[1]), [c_gen[5]], Flux.relu)
    lin_model_δ = Flux.Dense(reshape(c_gen[6:9], lin_arc[2], lin_arc[1]), [c_gen[10]], Flux.relu)

    nondim_ε = lin_model_ε(nondim_groups)[1]
    nondim_δ = lin_model_δ(nondim_groups)[1]
    return nondim_ε, nondim_δ
end


"""
    non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::RFEntr)

Uses a Random Feature model to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: RFEntr - basic RF entrainment closure
"""
function non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, ::RFEntr)
    # d=4 inputs, p=2 outputs, m random features
    nondim_groups = non_dimensional_groups(param_set, εδ_model_vars, entr_pi_groups)
    d = size(nondim_groups)[1]

    # Learnable and fixed parameters
    c_rf_fix = ICP.c_rf_fix(param_set)      # 2 x m x (1 + d), fix
    c_rf_fix = reshape(c_rf_fix, 2, :, 1 + d)
    m = size(c_rf_fix)[2]
    c_rf_opt = ICP.c_rf_opt(param_set)      # 2 x (m + 1 + d), learn
    c_rf_opt = reshape(c_rf_opt, 2, m + 1 + d)

    # Random Features
    scale_x_entr = (c_rf_opt[1, (m + 2):(m + d + 1)] .^ 2) .* nondim_groups
    scale_x_detr = (c_rf_opt[2, (m + 2):(m + d + 1)] .^ 2) .* nondim_groups
    f_entr = c_rf_opt[1, m + 1]^2 * sqrt(2) * cos.(c_rf_fix[1, :, 2:(d + 1)] * scale_x_entr + c_rf_fix[1, :, 1])
    f_detr = c_rf_opt[2, m + 1]^2 * sqrt(2) * cos.(c_rf_fix[2, :, 2:(d + 1)] * scale_x_detr + c_rf_fix[2, :, 1])

    # Square output for nonnegativity for prediction
    nondim_ε = sum(c_rf_opt[1, 1:m] .* f_entr) / m
    nondim_δ = sum(c_rf_opt[2, 1:m] .* f_detr) / m
    return nondim_ε^2, nondim_δ^2
end

"""
    non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, εδ_model_type::LogNormalScalingProcess)

Uses a LogNormal random variable to scale a deterministic process
to predict the non-dimensional components of dynamical entrainment/detrainment.

Arguments:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: LogNormalScalingProcess - Stochastic lognormal scaling
"""
function non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, εδ_model_type::LogNormalScalingProcess)
    FT = eltype(εδ_model_vars)
    # model parameters
    mean_model = εδ_model_type.mean_model
    c_gen_stoch = ICP.c_gen_stoch(param_set)
    ε_σ² = c_gen_stoch[1]
    δ_σ² = c_gen_stoch[2]

    # Mean model closure
    ε_mean_nondim, δ_mean_nondim = non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, mean_model)

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
    non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, εδ_model_type::NoisyRelaxationProcess)

Uses a noisy relaxation process to predict the non-dimensional components
of dynamical entrainment/detrainment. A deterministic closure is used as the
equilibrium mean function for the relaxation process.

Arguments:
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NoisyRelaxationProcess - A noisy relaxation process closure
"""
function non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, εδ_model_type::NoisyRelaxationProcess)
    # model parameters
    mean_model = εδ_model_type.mean_model
    c_gen_stoch = ICP.c_gen_stoch(param_set)
    ε_σ² = c_gen_stoch[1]
    δ_σ² = c_gen_stoch[2]
    ε_λ = c_gen_stoch[3]
    δ_λ = c_gen_stoch[4]

    # Mean model closure
    ε_mean_nondim, δ_mean_nondim = non_dimensional_function(param_set, εδ_model_vars, entr_pi_groups, mean_model)

    # noisy relaxation process
    ε_u0 = εδ_model_vars.ε_nondim
    δ_u0 = εδ_model_vars.δ_nondim
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
