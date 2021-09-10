if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
using Test
import TurbulenceConvection
const TC = TurbulenceConvection
const tc_dir = dirname(dirname(pathof(TC)))

include(joinpath(tc_dir, "integration_tests", "utils", "main.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "compute_mse.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "mse_tables.jl"))
using .NameList
best_mse = all_best_mse["Bomex"]

@testset "Unit tests" begin
    # TODO: add unit tests.
end

@testset "Bomex" begin

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

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
