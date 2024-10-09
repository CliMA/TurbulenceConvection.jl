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
import OperatorFlux as OF

const liq_type = CM.CommonTypes.LiquidType()
const ice_type = CM.CommonTypes.IceType()
const rain_type = CM.CommonTypes.RainType()
const snow_type = CM.CommonTypes.SnowType()

const Blk1MVel = CMT.Blk1MVelType() # for terminal velocity
const Chen2022Vel = CMT.Chen2022Type() # for terminal velocity


function get_termvel_type(termvel_type::Symbol)
    if termvel_type == :Blk1MVel
        return Blk1MVel
    elseif termvel_type == :Chen2022Vel
        return Chen2022Vel
    else
        error("Unknown termvel_type: $termvel_type")
    end
end

# CMP = CM.Parameters
import CloudMicrophysics.Parameters as CMP
# CT = CMT
const ACMP = CMP.AbstractCloudMicrophysicsParameters #
# const TDP = TD.parameters
const TDPS = TD.Parameters.ThermodynamicsParameters


"""
Because we pack our parameters in a named tuple and put symbols in Vals so they're isbits, we use a convenience get fcn.
Put in here so can use here and in driver
"""
function get_isbits_nt(
    named_tuple::NamedTuple,
    symbol::Symbol,
    default = nothing
)

    val = get(named_tuple, symbol, default)

    if isa(val, Val) # extract from it's Val wrapped prison
        val = typeof(val).parameters[1] # this is a hack to get the case name from the Val object
    end
    return val
end

include("Parameters.jl")
import .Parameters as TCP
const APS = TCP.AbstractTurbulenceConvectionParameters

up_sum(vals::AbstractArray) = reshape(sum(vals; dims = 1), size(vals, 2))

function parse_namelist(namelist, keys...; default = nothing, valid_options = nothing)
    param = namelist
    for k in keys
        if haskey(param, k)
            param = param[k]
            if valid_options ≠ nothing && !(param isa Dict)
                @assert param in valid_options
            end
        else
            if default == nothing
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

include("Grid.jl")
include("dycore_api.jl")
include("diagnostics.jl")
include("Fields.jl")
include("types.jl")
include("name_aliases.jl")

include("microphysics_coupling.jl")
include("turbulence_functions.jl")
include("utility_functions.jl")
include("variables.jl")
include("EDMF_Precipitation.jl")
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
include("closures/neural_microphysics_relaxation_timescales.jl") # testing different microphysics relaxation timescales
include("closures/korolev_mazin_2007.jl")
include("closures/morrison_milbrandt_2015_style.jl")
include("closures/sedimentation.jl")
include("closures/terminal_velocity.jl")

end
