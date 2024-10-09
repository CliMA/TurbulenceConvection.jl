
using Flux
using JLD2

function predict_τ(ρ,T,q, w, NN; FT=Float32, norm=1.)
    # normalize
    x_0 = FT.([ρ, T, q.liq, q.ice, w]) ./ FT.(norm)
    log_τ_liq, log_τ_ice = NN(x_0)
    return exp10(log_τ_liq), exp10(log_τ_ice)
end

NN_to_vec(NN::Flux.Chain) =  Flux.destructure(NN) # see Flux.destructure, returns vector (θ), and representation (re)
vec_to_NN(param_vector::Vector, re) = re(param_vector) # see Flux.destructure  
model_destructure_re_from_file(filename; var="re") = JLD2.load(filename, var)
model_destructure_re_to_file(re, filename) = JLD2.@save filename re

model_state_from_file(filename) =  JLD2.load(filename, "model_state")
model_state_to_file(model_state, filename) = JLD2.@save filename model_state
