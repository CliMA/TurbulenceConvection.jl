module TurbulenceConvection

import ClimaCore
import LinearAlgebra
import DocStringExtensions
import StaticArrays
import StatsBase
import NCDatasets
import JSON
import Dierckx
import RootSolvers
import Statistics
import LambertW
import Thermodynamics
import Distributions
import FastGaussQuadrature
import CLIMAParameters
import OrdinaryDiffEq
import CloudMicrophysics
import SurfaceFluxes
import UnPack

const ODE = OrdinaryDiffEq
const CC = ClimaCore
const CCO = CC.Operators
const NC = NCDatasets
const SA = StaticArrays

const TD = Thermodynamics
const SF = SurfaceFluxes
const UF = SF.UniversalFunctions
const RS = RootSolvers

const CM = CloudMicrophysics
const CM0 = CloudMicrophysics.Microphysics_0M
const CM1 = CloudMicrophysics.Microphysics_1M
const liq_type = CM1.LiquidType()
const ice_type = CM1.IceType()
const rain_type = CM1.RainType()
const snow_type = CM1.SnowType()

const CP = CLIMAParameters
const CPP = CP.Planet
const CPSGS = CP.SubgridScale
const APS = CP.AbstractEarthParameterSet

const CPMP = CP.Atmos.Microphysics
const CPEDMF = CP.Atmos.EDMF

# For dispatching to inherited class
struct BaseCase end

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
                error("No default value given for parameter ($(join(keys, ", "))).")
            else
                @info "Using default value, $default, for parameter ($(join(keys, ", ")))."
                return default
            end
        end
    end
    return param
end

include("ClimaParams.jl")
import .ClimaParams
const ICP = ClimaParams # internal clima parameters

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
        vec(prog_gm.θ_liq_ice),
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

include("parameters.jl")

include("Grid.jl")
include("dycore_api.jl")
include("NetCDFIO.jl")
include("diagnostics.jl")
include("Fields.jl")
include("ReferenceState.jl")
include("types.jl")
include("name_aliases.jl")
include("Operators.jl")

include("variable_tags.jl")

include("microphysics_coupling.jl")
include("turbulence_functions.jl")
include("utility_functions.jl")
include("TimeStepping.jl")
include("Variables.jl")
include("EDMF_Precipitation.jl")
include("EDMF_Environment.jl")
include("EDMF_Updrafts.jl")
include("stochastic_closures.jl")
include("update_aux.jl")
include("Turbulence_PrognosticTKE.jl")
include("Forcing.jl")
include("Radiation.jl")
include("Surface.jl")
include("io.jl")
include("thermodynamics.jl")
include("surface_functions.jl")
include("closures/perturbation_pressure.jl")
include("closures/entr_detr.jl")
include("closures/nondimensional_exchange_functions.jl")
include("closures/mixing_length.jl")
include("closures/buoyancy_gradients.jl")

end
