
using Flux
using JLD2

FT = Float32
# ρ_0, T_0, q_liq_0, q_ice_0, w_0 = FT(1.), FT(273.15), FT(1e-4), FT(1e-7), FT(1e-3) # characteristic values
ρ_0, T_0, q_liq_0, q_ice_0, w_0 = FT(1.), FT(1), FT(1e-4), FT(1e-7), FT(1e-3) # characteristic values, shifted -- still divide liquid values despite log scale, moves us up towards 0....
x_0_characteristic = [ρ_0, T_0, q_liq_0, q_ice_0, w_0] 

function prepare_for_NN(ρ, T, q_liq, q_ice, w; FT=Float32, norm = x_0_characteristic) 
    # normalize

    ρ = ρ ./ norm[1]
    T = (T.-273.15) ./ norm[2]
    q_liq = q_liq ./ norm[3]
    q_ice = q_ice ./ norm[4]
    w = w ./ norm[5]

    # log of the condensate amounts since they have a large range
    q_liq = log10.(FT(1e-10) .+ q_liq)
    q_ice = log10.(FT(1e-10) .+ q_ice)

    return FT.([ρ, T, q_liq, q_ice, w])
end

function predict_τ(ρ, T, q, w, NN; FT = Float32, norm = x_0_characteristic)
    # normalize

    x_0 = prepare_for_NN(ρ, T, q.liq, q.ice, w; norm = norm)

    log_τ_liq, log_τ_ice = NN(x_0)
    return exp10(log_τ_liq), exp10(log_τ_ice)
end

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
