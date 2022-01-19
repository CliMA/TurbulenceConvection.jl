import Pkg
Pkg.develop(path = ".")
import TurbulenceConvection
using Test

const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
import .NameList

best_mse = all_best_mse["GABLS"]

case_name = "GABLS"
println("Running $case_name...")
namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "01"
ds_tc_filename, return_code = main(namelist)

computed_mse = compute_mse_wrapper(
    case_name,
    best_mse,
    ds_tc_filename;
    plot_comparison = true,
    t_start = 7 * 3600,
    t_stop = 9 * 3600,
)

open("computed_mse_$case_name.json", "w") do io
    JSON.print(io, computed_mse)
end

@testset "GABLS" begin
    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    include(joinpath(tc_dir, "post_processing", "post_run_tests.jl"))
    nothing
end
