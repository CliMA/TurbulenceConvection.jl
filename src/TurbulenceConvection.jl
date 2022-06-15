module TurbulenceConvection

import ClimaCore
import LinearAlgebra
import ClimaCore.Geometry: ⊗
import LinearAlgebra: ×
import DocStringExtensions
import StaticArrays
import StatsBase
import Dierckx
import LambertW
import Thermodynamics
import Distributions
import FastGaussQuadrature
import CloudMicrophysics
import UnPack
import Random
import StochasticDiffEq
import Flux
import OperatorFlux
const OF = OperatorFlux

const SDE = StochasticDiffEq
const CC = ClimaCore
const CCG = CC.Geometry
const CCO = CC.Operators
const SA = StaticArrays

const TD = Thermodynamics
const LA = LinearAlgebra

const CM = CloudMicrophysics
const CMNe = CloudMicrophysics.MicrophysicsNonEq
const CM0 = CloudMicrophysics.Microphysics0M
const CM1 = CloudMicrophysics.Microphysics1M
const liq_type = CM.CommonTypes.LiquidType()
const ice_type = CM.CommonTypes.IceType()
const rain_type = CM.CommonTypes.RainType()
const snow_type = CM.CommonTypes.SnowType()

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractEarthParameterSet

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

include("TurbulenceConvectionParameters.jl")
import .TurbulenceConvectionParameters
const TCP = TurbulenceConvectionParameters

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
        vec(prog_gm.ρe_tot),
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
include("Operators.jl")

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

thermodynamics_params(ps::APS) = ps.thermo_params
surface_fluxes_params(ps::APS) = ps.surf_flux_params
microphysics_params(ps::APS) = ps.microphys_params

end
