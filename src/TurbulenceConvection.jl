module TurbulenceConvection

using StaticArrays
using StatsBase
using LinearAlgebra
import DocStringExtensions

import Dierckx
import Statistics
import LambertW
import Thermodynamics
import Distributions
import FastGaussQuadrature
import CLIMAParameters

const TD = Thermodynamics

const CP = CLIMAParameters
const CPP = CP.Planet
const APS = CP.AbstractEarthParameterSet

const CPMP = CP.Atmos.Microphysics
const CPEDMF = CP.Atmos.EDMF
const CPSGS = CP.Atmos.SubgridScale

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
include("Grid.jl")
include("Fields.jl")
include("Operators.jl")
include("NetCDFIO.jl")
include("ReferenceState.jl")
include("types.jl")

include("set_bcs.jl")
include("thermodynamic_functions.jl")
include("microphysics_functions.jl")
include("turbulence_functions.jl")
include("utility_functions.jl")
include("TimeStepping.jl")
include("Variables.jl")
include("EDMF_Rain.jl")
include("EDMF_Environment.jl")
include("EDMF_Updrafts.jl")
include("Turbulence_PrognosticTKE.jl")
include("Forcing.jl")
include("Radiation.jl")
include("forcing_functions.jl")
include("Surface.jl")
include("surface_functions.jl")
include("closures/perturbation_pressure.jl")
include("closures/entr_detr.jl")
include("closures/nondimensional_exchange_functions.jl")

end
