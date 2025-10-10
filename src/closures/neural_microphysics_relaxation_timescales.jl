
using Flux
using JLD2

#=
Julia has a bug where you can't both define a global at runtime (e.g. get_relaxation_timescale_type() and modify it later (e.g. in get_τ_helper()
because you can't reference global in get_τ_helper() because it will be parsed into the namespace without a type...

You can get away with it for neural_network because we do not modify it, but for the Cache, no such luck
=#


# const FTNN = Float32 # moved to top level. For pre-training NN, just define it yourself...
# const ρ_0, T_0, q_liq_0, q_ice_0, w_0 = FTNN(1.), FTNN(273.15), FTNN(1e-4), FTNN(1e-7), FTNN(1e-3) # characteristic values
const ρ_0, T_0, q_liq_0, q_ice_0, w_0 = FTNN(1.0), FTNN(1), FTNN(1e-4), FTNN(1e-7), FTNN(1e-2) # characteristic values, shifted -- still divide liquid values despite log scale, moves us up towards 0....
const x_0_characteristic = FTNN[ρ_0, T_0, q_liq_0, q_ice_0, w_0]
const x_0_characteristic_tuple = Tuple(x_0_characteristic)

# global neural_network::Flux.Chain{Tuple{Flux.Dense{typeof(NNlib.relu), Matrix{Float32}, Vector{Float32}}, Flux.Dense{typeof(NNlib.relu), Matrix{Float32}, Vector{Float32}}, Flux.Dense{typeof(NNlib.relu), Matrix{Float32}, Vector{Float32}}, Flux.Dense{typeof(identity), Matrix{Float32}, Vector{Float32}}}}
# global NN_cache::Pair{NTuple{5, Float64}, NTuple{4, FTNN}} # this is rough bc we don't always know that Float64 will match param_set... hence why a runtime definition would be nice...


function prepare_for_NN(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT <: Real}
    # normalize
    ρ /= norm[1]
    T = (T - FTT(273.15)) / norm[2] # not sure if is better than two atomic poerations
    q_liq /= norm[3]
    q_ice /= norm[4]

    # ρ = ρ / norm[1]
    # T = (T - 273.15) / norm[2]
    # q_liq = q_liq / norm[3]
    # q_ice = q_ice / norm[4]
    # w = w / norm[5]

    #= w can go up to about 1m/s. about 1 cm/s is enough to impact updraft significantly. so for the NN, we can maybe take raw w, but tanh limit it so it saturates around 2 m/s. However this can't scale with teh other inputs which are much less bounded
        - w = tanh(w / norm)  ::  asymptotes to 1, so you have to be careful to choose a scale commensurate with the other variables (which go up to about 10) M 1 
        - w = sign(w) * log10(1 + abs(w) / norm)  ::  never fully asymptotes. however it's about 30 ns
        -w = sign(x) * abs(x/norm)^(1/p) :: Also never fully asymptotes and actually keps growing more quickly, but you've got to pick a power as well which is arbitrary. about 5ns
    =#
    # w = tanh.(w ./ norm[5]) 
    w = sign(w) * log10(1 + abs(w) / norm[5])
    # w = sign(w) * abs(w ./ norm[5])^(1/2) # this is about 5ns faster than the tanh, but it doesn't saturate, so we have to be careful with the scale of w

    # log of the condensate amounts since they have a large range
    q_liq = log10(FTT(1e-10) + q_liq)
    q_ice = log10(FTT(1e-10) + q_ice)


    # return FTNN.([ρ, T, q_liq, q_ice, w])
    return SA.SVector{5,FTNN}(ρ, T, q_liq, q_ice, w) # use SVector for better performance. also if using static_strided_aray from svector, matching input to be svector has best performance. doesn't really impact Flux at all.
end
prepare_for_NN(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, norm::AbstractVector{FTT}) where {FTT <: Real} = prepare_for_NN(ρ, T, q_liq, q_ice, w, Tuple(norm)) # default in method is already tuple for norm

function prepare_for_NN!(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT <: Real}
    # normalize
    # in place versions -- it's a little dangerous..
    ρ ./= norm[1]
    T .-= FTT(273.15)
    T ./= norm[2]
    q_liq ./= norm[3]
    q_ice ./= norm[4]
    # w ./= norm[5]
    @. w = sign(w) * log10(1 + abs(w) / norm[5])

    # log of the condensate amounts since they have a large range
    @. q_liq = log10(FTT(1e-10) + q_liq)
    @. q_ice = log10(FTT(1e-10) + q_ice)

    
    ρ .= FTNN.(ρ) # convert to FTNN
    T .= FTNN.(T)
    q_liq .= FTNN.(q_liq)
    q_ice .= FTNN.(q_ice)
    w .= FTNN.(w)

    return ρ, T, q_liq, q_ice, w # return the modified vectors
end
function prepare_for_NN(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT <: Real}
    ρ = ρ ./ norm[1]
    T = (T .- FTT(273.15)) ./ norm[2]
    q_liq = q_liq ./ norm[3]
    q_ice = q_ice ./ norm[4]
    # w = w ./ norm[5]
    w = sign.(w) .* log10.(1 .+ abs.(w) ./ norm[5])

    @. q_liq = log10(FTT(1e-10) + q_liq)
    @. q_ice = log10(FTT(1e-10) + q_ice)
    return FTNN.(ρ), FTNN.(T), FTNN.(q_liq), FTNN.(q_ice), FTNN.(w)
end
prepare_for_NN!(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, norm::AbstractVector{FTT}) where {FTT <: Real} = prepare_for_NN!(ρ, T, q_liq, q_ice, w, Tuple(norm)) # default in method is already tuple for norm
prepare_for_NN(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, norm::AbstractVector{FTT}) where {FTT <: Real} = prepare_for_NN(ρ, T, q_liq, q_ice, w, Tuple(norm)) # default in method is already tuple for norm


# function predict_τ(ρ::FTT, T::FTT, q::TD.PhasePartition{FTT}, w::FTT, NN::Flux.Chain, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT}
# note, we could use a Vector instead of a SimpleChains.StrideArraysCore.StaticStrideArray  and it would be similarly fast, we can't directly pass an SVector though... rn our conversion is via to_static_strided_array().
function predict_τ(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, NN::Union{Flux.Chain, Tuple{SimpleChains.SimpleChain, SimpleChains.StrideArraysCore.StaticStrideArray }}, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT <: Real}
    # normalize

    x_0 = prepare_for_NN(ρ, T, q_liq, q_ice, w, norm)

    if (NN isa Flux.Chain)
        log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = NN(x_0)
    else
        log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = NN[1](x_0, NN[2]) # NN[1] is the SimpleChain, NN[2] is the parameters
    end

    # return exp10(log_τ_liq), exp10(log_τ_ice), exp10(log_N_liq), exp10(log_N_ice)
    return FTT(exp10(log_τ_liq)), FTT(exp10(log_τ_ice)), FTT(exp10(log_N_liq)), FTT(exp10(log_N_ice))

    # log_τ_liq, log_τ_ice, N_liq, N_ice = NN(x_0)
    # return exp10(log_τ_liq), exp10(log_τ_ice), N_liq, N_ice
end
# predict_τ(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, NN::Flux.Chain, norm::AbstractVector{FTT}) where {FTT} = predict_τ(ρ, T, q, w, NN, Tuple(norm)) # default in method is already tuple for norm
predict_τ(ρ::FTT, T::FTT, q_liq::FTT, q_ice::FTT, w::FTT, NN::Union{Flux.Chain, Tuple{SimpleChains.SimpleChain, SimpleChains.StrideArraysCore.StaticStrideArray }}, norm::AbstractVector{FTT}) where {FTT <: Real} = predict_τ(ρ, T, q_liq, q_ice, w, NN, Tuple(norm)) # default in method is already tuple for norm
# predict_τ(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, NN::Union{Flux.Chain, Tuple{SimpleChains.SimpleChain, SimpleChains.StrideArraysCore.StaticStrideArray }}, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT} = predict_τ.(ρ, T, q_liq, q_ice, w, NN, norm)

# this really should make a matrix call for efficiency, since you can already use the dot version yourself
function predict_τ(ρ::AbstractVector{FTT}, T::AbstractVector{FTT}, q_liq::AbstractVector{FTT}, q_ice::AbstractVector{FTT}, w::AbstractVector{FTT}, NN::Union{Flux.Chain, Tuple{SimpleChains.SimpleChain, SimpleChains.StrideArraysCore.StaticStrideArray }}, norm::NTuple{5,FTT} = x_0_characteristic_tuple) where {FTT <: Real}
    # normalize
    x_0 = prepare_for_NN(ρ, T, q_liq, q_ice, w, norm)

    if (NN isa Flux.Chain)
        # log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = eachcol(NN(reduce(hcat, x_0)'))
        # NN takes [n_in x n_obs] input and returns [n_out x n_obs] output. So we need to stack the input vectors on dims=1 to get [n_obs x n_in] input, and take eachrow() from the output
        log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = eachrow(NN(stack(x_0, dims=1))) # stack on dims_1 should eliminate the need for a transpose I think, which means we wont have to cast back to Matrix for simple chains either (SimpleChains doesn't accept a transpose)
    else
        # log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = eachcol(NN[1](reduce(hcat, x_0), NN[2])) # NN[1] is the SimpleChain, NN[2] is the parameters
        log_τ_liq, log_τ_ice, log_N_liq, log_N_ice = eachrow(NN[1](stack(x_0, dims=1), NN[2])) # NN[1] is the SimpleChain, NN[2] is the parameters
    end

    # @warn "status: log_τ_liq = $log_τ_liq, log_τ_ice = $log_τ_ice, log_N_liq = $log_N_liq, log_N_ice = $log_N_ice"

    # alias (i hope)
    τ_liq = log_τ_liq
    τ_ice = log_τ_ice
    N_liq = log_N_liq
    N_ice = log_N_ice

    # in place
    @. τ_liq  = FTT(exp10(τ_liq))
    @. τ_ice  = FTT(exp10(τ_ice))
    @. N_liq  = FTT(exp10(N_liq))
    @. N_ice  = FTT(exp10(N_ice))

    return τ_liq, τ_ice, N_liq, N_ice

end

NN_to_vec(NN::Flux.Chain) = Flux.destructure(NN) # see Flux.destructure, returns vector (θ), and representation (re)
vec_to_NN(param_vector::Vector, re) = re(param_vector) # see Flux.destructure  
model_destructure_re_from_file(filename::String; var = "re") = JLD2.load(filename, var)


function simple_chain_model_from_file(filename::String) #, param_vector::Union{AbstractVector{FTNN}, Nothing} = nothing)
    NN_simple_chain = JLD2.load(filename, "NN_simple_chain",) # don't need to load nn_simple_chain_params, we can do that when creating the param_set, don't need to load x_0_characteristic, we don't use it from the file it's just defined here.
    # if !isnothing(param_vector)
    #     nn_simple_chain_params = SimpleChains.init_params(NN_simple_chain) # this is not needed, we already have the params
    #     nn_simple_chain_params .= param_vector # should this be static or does it not matter?
    # end
    return NN_simple_chain
end


# model_destructure_re_from_file(filename; var = "re") = begin 
#     @info "Loading $var from $filename"
#     println("Loading $var from $filename")
#     @warn "Loading $var from $filename"
#     JLD2.load(filename, var)
# end

model_destructure_re_to_file(re, filename::String) = JLD2.@save filename re

model_state_from_file(filename::String) = JLD2.load(filename, "model_state")
model_state_to_file(model_state, filename::String) = JLD2.@save filename model_state
