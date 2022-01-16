import TurbulenceConvection
import Flux
const TC = TurbulenceConvection
const ICP = TC.ClimaParams

"""
    TC.non_dimensional_function(param_set, εδ_model_vars, ::NNEntr)

Uses a fully connected neural network to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: NNEntr - Neural network entrainment closure
"""
function TC.non_dimensional_function(param_set, εδ_model_vars, ::TC.NNEntr)
    c_gen = ICP.c_gen(param_set)

    # Neural network closure
    nn_arc = (4, 2, 2)  # (#inputs, #neurons, #outputs)
    nn_model = Flux.Chain(
        Flux.Dense(reshape(c_gen[1:8], nn_arc[2], nn_arc[1]), c_gen[9:10], Flux.sigmoid),
        Flux.Dense(reshape(c_gen[11:14], nn_arc[3], nn_arc[2]), c_gen[15:16], Flux.softplus),
    )

    nondim_groups = TC.non_dimensional_groups(param_set, εδ_model_vars)
    nondim_ε, nondim_δ = nn_model(nondim_groups)
    return nondim_ε, nondim_δ
end

"""
    TC.non_dimensional_function(param_set, εδ_model_vars, ::LinearEntr)

Uses a simple linear model to predict the non-dimensional components of dynamical entrainment/detrainment.
 - `param_set`      :: parameter set
 - `εδ_model_vars`  :: structure containing variables
 - `εδ_model_type`  :: LinearEntr - linear entrainment closure
"""
function TC.non_dimensional_function(param_set, εδ_model_vars, ::TC.LinearEntr)
    c_gen = ICP.c_gen(param_set)

    # Linear closure
    lin_arc = (4, 1)  # (#weights, #outputs)
    lin_model_ε = Flux.Dense(reshape(c_gen[1:4], lin_arc[2], lin_arc[1]), [c_gen[5]], Flux.relu)
    lin_model_δ = Flux.Dense(reshape(c_gen[6:9], lin_arc[2], lin_arc[1]), [c_gen[10]], Flux.relu)

    nondim_groups = TC.non_dimensional_groups(param_set, εδ_model_vars)
    nondim_ε = lin_model_ε(nondim_groups)[1]
    nondim_δ = lin_model_δ(nondim_groups)[1]
    return nondim_ε, nondim_δ
end
