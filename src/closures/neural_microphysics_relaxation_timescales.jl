
using Flux
using JLD2

const FTNN = Float32
# ρ_0, T_0, q_liq_0, q_ice_0, w_0 = FTNN(1.), FTNN(273.15), FTNN(1e-4), FTNN(1e-7), FTNN(1e-3) # characteristic values
ρ_0, T_0, q_liq_0, q_ice_0, w_0 = FTNN(1.0), FTNN(1), FTNN(1e-4), FTNN(1e-7), FTNN(1e-3) # characteristic values, shifted -- still divide liquid values despite log scale, moves us up towards 0....
x_0_characteristic = FTNN[ρ_0, T_0, q_liq_0, q_ice_0, w_0]
x_0_characteristic_tuple = Tuple(x_0_characteristic)

function prepare_for_NN(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT}
    # normalize

    ρ = ρ / norm[1]
    T = (T - 273.15) ./ norm[2]
    q_liq = q_liq / norm[3]
    q_ice = q_ice / norm[4]
    w = w / norm[5]

    # log of the condensate amounts since they have a large range
    q_liq = log10(FTT(1e-10) + q_liq)
    q_ice = log10(FTT(1e-10) + q_ice)

    return FTNN.([ρ, T, q_liq, q_ice, w])
end
prepare_for_NN(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, norm::AbstractVector{FTT}) where {FTT} = prepare_for_NN(ρ, T, q_liq, q_ice, w, Tuple(norm)) # default in method is already tuple for norm

function prepare_for_NN(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}; norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT}
    # normalize

    ρ = ρ ./ norm[1]
    T = (T .- 273.15) ./ norm[2]
    q_liq = q_liq ./ norm[3]
    q_ice = q_ice ./ norm[4]
    w = w ./ norm[5]

    # log of the condensate amounts since they have a large range
    q_liq = log10.(FTT(1e-10) .+ q_liq)
    q_ice = log10.(FTT(1e-10) .+ q_ice)

    return (FTNN.(ρ), FTNN.(T), FTNN.(q_liq), FTNN.(q_ice), FTNN.(w))
end

prepare_for_NN(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, norm::AbstractVector{FTT}) where {FTT} = prepare_for_NN(ρ, T, q_liq, q_ice, w, Tuple(norm))

   
function predict_τ(ρ::FTT, T::FTT, q::TD.PhasePartition{FTT}, w::FTT, NN::Flux.Chain, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT}
    # normalize

    x_0 = prepare_for_NN(ρ, T, q.liq, q.ice, w, norm)

    log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = NN(x_0)
    return exp10(log_τ_liq), exp10(log_τ_ice), exp10(log_N_liq), exp10(log_N_ice)

    # log_τ_liq, log_τ_ice, N_liq, N_ice = NN(x_0)
    # return exp10(log_τ_liq), exp10(log_τ_ice), N_liq, N_ice
end
predict_τ(ρ::FTT, T::FTT, q::TD.PhasePartition{FTT}, w::FTT, NN::Flux.Chain, norm::AbstractVector{FTT}) where {FTT} = predict_τ(ρ, T, q, w, NN, Tuple(norm)) # default in method is already tuple for norm

NN_to_vec(NN::Flux.Chain) = Flux.destructure(NN) # see Flux.destructure, returns vector (θ), and representation (re)
vec_to_NN(param_vector::Vector, re) = re(param_vector) # see Flux.destructure  
model_destructure_re_from_file(filename; var = "re") = JLD2.load(filename, var)


# model_destructure_re_from_file(filename; var = "re") = begin 
#     @info "Loading $var from $filename"
#     println("Loading $var from $filename")
#     @warn "Loading $var from $filename"
#     JLD2.load(filename, var)
# end

model_destructure_re_to_file(re, filename) = JLD2.@save filename re

model_state_from_file(filename) = JLD2.load(filename, "model_state")
model_state_to_file(model_state, filename) = JLD2.@save filename model_state
