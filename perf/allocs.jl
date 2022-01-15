import Pkg
Pkg.develop(path = ".")

# Track allocations in TC.jl plus all _direct_ dependencies
exhaustive = "exhaustive=true" in ARGS
@show exhaustive

import ReportMetrics
const RM = ReportMetrics

# Packages to monitor
import TurbulenceConvection
import CLIMAParameters
import ClimaCore
import CloudMicrophysics
import Dierckx
import Distributions
import FastGaussQuadrature
import Flux
import LinearAlgebra
import NCDatasets
import OrderedCollections
import OrdinaryDiffEq
import PoissonRandom
import Random
import RootSolvers
import SciMLBase
import StaticArrays
import Statistics
import StatsBase
import StochasticDiffEq
import SurfaceFluxes
import TerminalLoggers
import Thermodynamics

deps_to_monitor = [
    CLIMAParameters,
    ClimaCore,
    CloudMicrophysics,
    Dierckx,
    Distributions,
    FastGaussQuadrature,
    Flux,
    LinearAlgebra,
    NCDatasets,
    OrderedCollections,
    OrdinaryDiffEq,
    PoissonRandom,
    Random,
    RootSolvers,
    SciMLBase,
    StaticArrays,
    Statistics,
    StatsBase,
    StochasticDiffEq,
    SurfaceFluxes,
    TerminalLoggers,
    Thermodynamics,
]

deps_to_monitor = exhaustive ? deps_to_monitor : Module[]
# Packages to monitor

all_cases = [
    # "ARM_SGP",
    # "Bomex",
    # "DryBubble",
    # "DYCOMS_RF01",
    # "GABLS",
    # "GATE_III",
    # "life_cycle_Tan2018",
    # "Nieuwstadt",
    "Rico",
    # "Soares",
    # "SP",
    "TRMM_LBA",
    "LES_driven_SCM",
]

# only one case for exhaustive alloc analysis
all_cases = exhaustive ? ["Bomex"] : all_cases
n_unique_allocs = 40

for case in all_cases
    ENV["ALLOCATION_CASE_NAME"] = case
    if exhaustive
        run_cmd = `$(Base.julia_cmd()) --project=test/ --track-allocation=all perf/alloc_per_case_with_init_io.jl`
    else
        run_cmd = `$(Base.julia_cmd()) --project=test/ --track-allocation=all perf/alloc_per_case.jl`
    end
    RM.report_allocs(;
        job_name = case,
        run_cmd = run_cmd,
        dirs_to_monitor = [RM.mod_dir(TurbulenceConvection), RM.mod_dir.(deps_to_monitor)...],
        process_filename = function process_fn(fn)
            fn = "TurbulenceConvection.jl/" * last(split(fn, "turbulenceconvection-ci/"))
            fn = last(split(fn, "depot/cpu/packages/"))
            return fn
        end,
        n_unique_allocs = n_unique_allocs,
    )
end
