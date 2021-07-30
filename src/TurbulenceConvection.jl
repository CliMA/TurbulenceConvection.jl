module TurbulenceConvection

import Dierckx
using StatsBase
import Statistics
using PoissonRandom: pois_rand
import LambertW
using Distributions: Normal, quantile
using OffsetArrays
using FastGaussQuadrature: gausshermite

# For dispatching to inherited class
struct BaseCase end

up_sum(vals::OffsetArray) = off_arr(reshape(sum(vals; dims=1), size(vals, 2)))

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
include("forcing_functions.jl")
include("Surface.jl")
include("surface_functions.jl")

function export_all(Case, Turb, GMV, TS, Stats)
    # if TS.i_iter == 1 && TS.isolate
    # if TS.i_iter == 10


    if TS.i_iter == 22
        println("------- inside export_all")
        @show Case.Sur.rho_hflux
        @show Case.Sur.rho_qtflux
        @show GMV.Gr.z_half[GMV.Gr.gw]
        @show Case.Sur.ustar
        @show Case.Sur.obukhov_length
        @show Turb.wstar
    # if TS.i_iter == 5759
        @show TS.i_iter
        open_files(Stats)
        write_simulation_time(Stats, TS.t)
        io(GMV, Stats)
        io(Case, Stats)
        io(Turb, Stats, TS)
        close_files(Stats)
        error("Done")
    end
end

end