module TurbulenceConvection

import ClimaCore as CC
import ClimaCore.Geometry as CCG
import ClimaCore.Geometry: ⊗
import ClimaCore.Operators as CCO
import LinearAlgebra as LA
import LinearAlgebra: ×
import DocStringExtensions
import StaticArrays as SA
import StatsBase
import Dierckx
import LambertW
import Thermodynamics as TD
import Distributions
import FastGaussQuadrature
import CloudMicrophysics as CM
import CloudMicrophysics.MicrophysicsNonEq as CMNe
import CloudMicrophysics.Microphysics0M as CM0
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.CommonTypes as CMT
import UnPack
import Random
import StochasticDiffEq as SDE
import Flux
import SimpleChains # much faster than Flux and its networks are isbits so they can go in the paramset etc...
import OperatorFlux as OF

# I think doing it in `export OPENBLAS_NUM_THREADS=1`  is better, cause if you can use BLAS threads, then this wont stop you. [ this won't run at import w/o an __init__() ]
# function __init__()
#     using LinearAlgebra
#     LinearAlgebra.BLAS.set_num_threads(1) # if you're on HPC this is essential to A\b not slowing down by 5 orders of magnitude from 1ms to 100s 
# end

const liq_type = CM.CommonTypes.LiquidType()
const ice_type = CM.CommonTypes.IceType()
const rain_type = CM.CommonTypes.RainType()
const snow_type = CM.CommonTypes.SnowType()

const Blk1MVel = CMT.Blk1MVelType() # for terminal velocity
const Chen2022Vel = CMT.Chen2022Type() # for terminal velocity

# const CMTWaterTypes = Union{CMT.AbstractCloudType, CMT.AbstractPrecipType}
const CMTWaterTypes = Union{CMT.LiquidType, CMT.IceType, CMT.RainType, CMT.SnowType} # more precise


function get_termvel_type(termvel_type::Symbol)
    if termvel_type === :Blk1MVel
        return Blk1MVel
    elseif termvel_type === :Chen2022Vel
        return Chen2022Vel
    else
        error("Unknown termvel_type: $termvel_type")
    end
end
get_termvel_type(termvel_type::String) = get_termvel_type(Symbol(termvel_type)) # for convenience   

# CMP = CM.Parameters
import CloudMicrophysics.Parameters as CMP
# CT = CMT
const ACMP = CMP.AbstractCloudMicrophysicsParameters #
# const TDP = TD.parameters
const TDPS = TD.Parameters.ThermodynamicsParameters

function safe_clamp(x, lo, hi)
    # current clamp doesn't safely work when min < max.  Because of its order of operations it returns x if lo < x < hi, then if x > hi, returns hi (bad if hi is lo), and if x < lo, returns lo (bad if lo is hi)
    return (lo < hi) ? clamp(x, lo, hi) : clamp(x, hi, lo)
end # safe_clamp
safe_clamp(x, lo_hi) =  safe_clamp(x, lo_hi[1], lo_hi[2]) # for convenience

resolve_nan(x::FT, val = FT(0.0)) where {FT} = isnan(x) ? FT(val) : x # replace nan w/ 0
function resolve_nan!(x::AbstractArray{FT}, val = FT(0.0)) where {FT}
    @inbounds for i in eachindex(x)
        x[i] = resolve_nan(x[i], FT(val))
    end
end

resolve_not_finite(x::FT, val = FT(NaN)) where {FT} = !isfinite(x) ? FT(val) : x # replace inf with NaN
function resolve_not_finite!(x::AbstractArray{FT}, val = FT(NaN)) where {FT}
    @inbounds for i in eachindex(x)
        x[i] = resolve_not_finite(x[i], FT(val))
    end
end 

# resolve_inf(x::FT; val::FT=FT(NaN)) where {FT} = isinf(x) ? val : x # replace inf with NaN

full_print(x) = show(IOContext(stdout, :limit => false), "text/plain", x)

"""
Maybe i'm dumb but i can't find the constructor to create a StrideArray from an existing array
StrideArraysCore.StrideArray(a)
and 
StrideArraysCore.StrideArray( svector(a) ) are close, but not quite the same type
"""

create_svector(x::AbstractVector) = SA.SVector{length(x)}(x) # create an SVector from an array, this is the same as SVector(x...) but faster [still slow bc length(x) is looked up at compile time]
create_svector(x::NTuple{N, FT}) where {FT, N} = SA.SVector{N, FT}(x) # create an SVector from a tuple, this is the same as SVector(x...) but faster [still slow bc length(x) is looked up at compile time]

# function to_strided_array(x::AbstractArray{FT}) where {FT}
#     out = StrideArray{FT}(undef, size(int))
#     out .= x
#     return out
# end

# function to_strided_array(x::NTuple{N, FT}) where {FT, N}
#     out = StrideArray{FT}(undef, N)
#     # out .= x # doesn't work, but this does
#     for i in eachindex(x)
#         out[i] = x[i]
#     end
#     return out
# end

# Note you could use StrideArrays.StaticInt or SimpleChains.StaticInt == SimpleChains.StrideArraysCore.StaticInt, but SimpleChains only makes StrideArraysCore available not StrideArrays, so we use SimpleChains.StaticInt: Both inherit from Static.StaticInt
# StrideArrays.StrideArray(svector) creates the right type of object but the values are wrong.... however the underlying .data correct so it must be the pointers that are broken
"""
Called on an SVector, this takes about 133 ns so about the same as just calling Vector. so we create the faster method below which takes 115ns
Calling NN_simiple_chain(input, ssv) takes about 90 ns
"""
function to_static_strided_array(x::AbstractArray{FT}) where {FT}  # this can be very slow on vectors... got  65.675 μs (201 allocations: 10.12 KiB)
    out = SimpleChains.StrideArray{FT}(undef, SimpleChains.StaticInt.(size(x)))
    out .= x
    return out
end

# 115 ns
function to_static_strided_array(x::SA.StaticVector{N, FT}) where {FT, N}
    out = SimpleChains.StrideArray{FT}(undef, SimpleChains.StaticInt(N))
    out .= x
    return out
end

function to_static_strided_array(x::NTuple{N, FT}) where {FT, N}
    out = SimpleChains.StrideArray{FT}(undef, SimpleChains.StaticInt(N))
    # out .= x # doesn't work, but this does
    for i in eachindex(x)
        out[i] = x[i]
    end
    return out
end

# """
# Vector(svec) takes like 138 ns...
# calling NN_simiple_chain(input, v) tkes about 87 ns
# """
# function fast_vec_from_svec(x::SVector{N, FT}) where {N, FT}
#     out = Vector{FT}(undef, N) # 125 ns
#     # out = similar(Vector{FT}, N) # this is faster than Vector{FT}(undef, length(x)) # 123 ns 
#     @inbounds for i in eachindex(x)
#         out[i] = x[i]
#     end
#     return out
# end 

@inline linear_interpolate_extrapolate(x::FT, xs::Tuple{FT,FT}, vs::Tuple{FT,FT}) where {FT} = (x - xs[1]) / (xs[2] - xs[1]) * (vs[2] - vs[1]) + vs[1] # fast linear interpolation/extrapolation
@inline positive_linear_interpolate_extrapolate(x::FT, xs::Tuple{FT,FT}, vs::Tuple{FT,FT}) where {FT} = max(0, linear_interpolate_extrapolate(x, xs, vs)) # fast linear interpolation/extrapolation



const FTNN = Float32 # Neural Network Type. declare as early as possible so it can be used both in closures and in types.jl for relaxation types


# NN_type = <> # Figure this out
# const neural_network = Ref{NN_type}() # empty init, fill in later...

#=
Julia has a bug where you can't both define a global at runtime (e.g. get_relaxation_timescale_type() and modify it later (e.g. in get_τ_helper()
because you can't reference global in get_τ_helper() because it will be parsed into the namespace without a type...

You can get away with it for neural_network because we do not modify it, but for the Cache, no such luck



I tried to put the definition in /home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/src/closures/neural_microphysics_relaxation_timescales.jl where both Flux and NNlib are avilable but it complained about
ERROR: LoadError: cannot set type for global TurbulenceConvection.NN_cache. It already has a value or is already set to a different type.
Maybe it parsed /home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/src/relaxation_timescale_types.jl first or something? I don't know... but it seems like a bug in the parsing and assignment logic...
=#
# using Flux
# global neural_network::Flux.Chain{Tuple{Flux.Dense{typeof(NNlib.relu), Matrix{Float32}, Vector{Float32}}, Flux.Dense{typeof(NNlib.relu), Matrix{Float32}, Vector{Float32}}, Flux.Dense{typeof(NNlib.relu), Matrix{Float32}, Vector{Float32}}, Flux.Dense{typeof(identity), Matrix{Float32}, Vector{Float32}}}}
# global NN_cache::Pair{Tuple{Float64,Float64,TD.PhasePartition{Float64},Float64} , NTuple{4, Float32}} # this is rough bc we don't always know that Float64 will match param_set... hence why a runtime definition would be nice...

# try using simplechains
# global neural_network::SimpleChain{Tuple{Static.StaticInt{5}}, Tuple{TurboDense{true, Static.StaticInt{10}, typeof(relu)}, TurboDense{true, Static.StaticInt{8}, typeof(relu)}, TurboDense{true, Static.StaticInt{4}, typeof(relu)}, TurboDense{true, Static.StaticInt{4}, typeof(identity)}}}
# global nn_simple_chain_params::StrideArraysCore.StaticStrideArray{Float32, 1, (1,), Tuple{Static.StaticInt{204}}, Tuple{Nothing}, Tuple{Static.StaticInt{1}}, 204}

# global neural_network::SimpleChains.SimpleChain{Tuple{SimpleChains.Static.StaticInt{5}}, Tuple{SimpleChains.TurboDense{true, SimpleChains.Static.StaticInt{10}, typeof(NNlib.relu)}, SimpleChains.TurboDense{true, SimpleChains.Static.StaticInt{8}, typeof(NNlib.relu)}, SimpleChains.TurboDense{true, SimpleChains.Static.StaticInt{4}, typeof(NNlib.relu)}, SimpleChains.TurboDense{true, SimpleChains.Static.StaticInt{4}, typeof(identity)}}}
# global nn_simple_chain_params::SimpleChains.StrideArraysCore.StaticStrideArray{Float32, 1, (1,), Tuple{SimpleChains.Static.StaticInt{204}}, Tuple{Nothing}, Tuple{SimpleChains.Static.StaticInt{1}}, 204}

include("Parameters.jl")
import .Parameters as TCP
const APS = TCP.AbstractTurbulenceConvectionParameters

up_sum(vals::AbstractArray) = reshape(sum(vals; dims = 1), size(vals, 2))

function parse_namelist(namelist, keys...; default = nothing, valid_options = nothing)
    param = namelist
    for k in keys
        if haskey(param, k)
            param = param[k]
            if !isnothing(valid_options) && !(param isa Dict)
                @assert param in valid_options
            end
        else
            if isnothing(default)
                error("No default value given for parameter (`$(join(keys, ", "))`).")
            else
                @info "Using default value, $default, for parameter (`$(join(keys, ", "))`)."
                return default
            end
        end
    end
    return param
end

Base.broadcastable(param_set::APS) = Ref(param_set)

#=
    debug_state(state, code_location::String)

A simple function for debugging the entire state,
specifically for when quantities that should remain
positive-definite become negative.

=#
function debug_state(state, code_location::String)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)

    prog_up = center_prog_updrafts(state)
    aux_up = center_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_up_f = face_aux_updrafts(state)

    prog_en = center_prog_environment(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)

    ######
    ###### Positive-definite variables
    ######

    vars_positive = [
        vec(prog_gm.ρθ_liq_ice),
        vec(prog_gm_f.w),
        vec(prog_up[1].ρarea),
        vec(prog_up[1].ρaθ_liq_ice),
        vec(prog_up_f[1].ρaw),
        vec(aux_en.area),
        vec(aux_en.θ_liq_ice),
    ]
    vars = vars_positive
    vars_conds = map(v -> any(v .< 0), vars)

    if any(vars_conds)
        @show code_location
        for (i, vc, v) in zip(1:length(vars), vars_conds, vars)
            vc || continue
            @show i, v
        end
        @show vars_conds
        error("Negative state for positive-definite field(s)")
    end

    ######
    ###### All listed variables
    ######
    vars = vars_positive
    vars_conds = map(v -> any(isnan.(v)) || any(isinf.(v)), vars)

    if any(vars_conds)
        @show code_location
        for (i, vc, v) in zip(1:length(vars), vars_conds, vars)
            vc || continue
            @show i, v
        end
        @show vars_conds
        error("Nan/Inf state for field(s)")
    end
end

const Ic = CCO.InterpolateF2C() # no bcs on F2C, this gets used all the time, so define it once here
const Ifx = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # C2F with extrapolate bcs, this gets used all the time, so define it once here
const ∇c = CCO.DivergenceF2C() 

include("Grid.jl")
include("dycore_api.jl")
include("diagnostics.jl")
include("Fields.jl")
include("types.jl")
include("name_aliases.jl")

# Domains we use [[ mostly for reweighting, and for threshold_acnv etc...]]
const Env = EnvDomain()
const Up = UpDomain()
const Bulk = BulkDomain()
#
const CloakUp = CloakUpDomain()
const CloakDown = CloakDownDomain()
const EnvRemaining = EnvRemainingDomain()

include("microphysics_coupling.jl")
include("turbulence_functions.jl")
include("utility_functions.jl")
include("variables.jl")
include("EDMF_Precipitation.jl")
include("closures/ql_qi_supersaturation_covariance_closure.jl") # for SGS quadrature w/ ql, qi
include("EDMF_Environment.jl")
include("EDMF_Updrafts.jl")
include("update_aux.jl")
include("EDMF_functions.jl")
include("thermodynamics.jl")
include("closures/perturbation_pressure.jl")
include("closures/entr_detr.jl")
include("closures/nondimensional_exchange_functions.jl")
include("closures/mixing_length.jl")
include("closures/buoyancy_gradients.jl")

include("closures/relaxation_timescales.jl")
include("closures/N_r_closures.jl") # testing different N/r distribution closures
include("closures/reweight_processes_for_grid.jl") # reweighting process rates to account for grid spacing...
include("closures/neural_microphysics_relaxation_timescales.jl") # testing different microphysics relaxation timescales
include("closures/korolev_mazin_2007.jl")
include("closures/morrison_milbrandt_2015_style.jl")
include("closures/morrison_milbrandt_2015_style_exponential_part_only.jl")
include("closures/sedimentation.jl")
include("closures/terminal_velocity.jl")

end
