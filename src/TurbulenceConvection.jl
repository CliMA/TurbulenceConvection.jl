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

const ODE = OrdinaryDiffEq
const CC = ClimaCore
const NC = NCDatasets
const SA = StaticArrays

const TD = Thermodynamics
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

include("python_primitives.jl")
include("parameters.jl")

include("thermodynamic_functions.jl")
import .TCThermodynamics
const TCTD = TCThermodynamics

include("Grid.jl")
include("dycore_api.jl")
include("NetCDFIO.jl")
include("diagnostics.jl")
include("Fields.jl")
include("ReferenceState.jl")
include("types.jl")
include("name_aliases.jl")
include("Operators.jl")

include("microphysics_functions.jl")
include("turbulence_functions.jl")
include("utility_functions.jl")
include("TimeStepping.jl")
include("Variables.jl")
include("EDMF_Rain.jl")
include("EDMF_Environment.jl")
include("EDMF_Updrafts.jl")
include("stochastic_closures.jl")
include("Turbulence_PrognosticTKE.jl")
include("Forcing.jl")
include("Radiation.jl")
include("Surface.jl")
include("surface_functions.jl")
include("closures/perturbation_pressure.jl")
include("closures/entr_detr.jl")
include("closures/nondimensional_exchange_functions.jl")
include("closures/mixing_length.jl")
include("closures/buoyancy_gradients.jl")

end
