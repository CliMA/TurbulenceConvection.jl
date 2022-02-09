if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end
using Test
import TurbulenceConvection
const TC = TurbulenceConvection
const tc_dir = dirname(dirname(pathof(TC)))

include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
import .NameList
best_mse = all_best_mse["Bomex"]

@testset "Unit tests" begin
    # TODO: add unit tests.
end

@testset "Bomex" begin

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename, return_code = main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 4 * 3600,
        t_stop = 6 * 3600,
    )
    nothing
end

@testset "ClimaCore extensions" begin
    include("clima_core_extensions.jl")
end
