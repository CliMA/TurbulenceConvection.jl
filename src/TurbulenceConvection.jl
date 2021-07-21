module TurbulenceConvection

import Dierckx
using StatsBase
import Statistics
using PoissonRandom: pois_rand
import LambertW
using Distributions: Normal, quantile
using OffsetArrays
using FastGaussQuadrature: gausshermite
using LinearAlgebra

# For dispatching to inherited class
struct BaseCase end

up_sum(vals::OffsetArray) = off_arr(reshape(sum(vals; dims = 1), size(vals, 2)))

function parse_param(namelist, keys...; default = nothing)
    @assert default ≠ nothing
    param = namelist
    for k in keys
        if haskey(param, k)
            param = param[k]
        else
            @info "Using default value, $default, for parameter ($(join(keys, ", ")))."
            return default
        end
    end
    return param
end

include("python_primitives.jl")
include("parameters.jl")
include("Grid.jl")
include("NetCDFIO.jl")
include("ReferenceState.jl")
include("types.jl")

include("thermodynamic_functions.jl")
include("microphysics_functions.jl")
include("turbulence_functions.jl")
include("utility_functions.jl")
include("TimeStepping.jl")
include("Variables.jl")
include("EDMF_Rain.jl")
include("EDMF_Environment.jl")
include("EDMF_Updrafts.jl")
include("Turbulence.jl")
include("Turbulence_PrognosticTKE.jl")
include("Forcing.jl")
include("Radiation.jl")
include("forcing_functions.jl")
include("Surface.jl")
include("surface_functions.jl")

end
