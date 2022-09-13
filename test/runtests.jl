using Test
using JSON
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

import Aqua
@testset "Quality Assurance" begin
    @test length(Aqua.detect_ambiguities(TurbulenceConvection; recursive = true)) ≤ 136
end

@testset "Bomex" begin

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    integrator, ds_tc_filenames, return_code = main(namelist)

    for ds_tc_filename in ds_tc_filenames
        computed_mse = compute_mse_wrapper(
            case_name,
            best_mse,
            ds_tc_filename;
            plot_comparison = true,
            t_start = 4 * 3600,
            t_stop = 6 * 3600,
        )
    end
    nothing

    # Test running from parsed namelist
    namelist_filename = "Output.Bomex.01/namelist_Bomex.in"
    read_namelist = JSON.parsefile(namelist_filename)
    read_namelist["meta"]["uuid"] = "02"
    read_namelist["time_stepping"]["t_max"] = read_namelist["time_stepping"]["dt_max"]
    main(read_namelist)
end

@testset "ClimaCore extensions" begin
    include("clima_core_extensions.jl")
end
