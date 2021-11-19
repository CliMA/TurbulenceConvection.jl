if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
include(joinpath("utils", "mse_tables.jl"))
import .NameList

best_mse = all_best_mse["life_cycle_Tan2018"]

case_name = "life_cycle_Tan2018"
println("Running $case_name...")
namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "01"
ds_tc_filename, return_code = @time main(namelist; time_run = true)

computed_mse = compute_mse_wrapper(
    case_name,
    best_mse,
    ds_tc_filename;
    plot_comparison = true,
    t_start = 4 * 3600,
    t_stop = 6 * 3600,
)

open("computed_mse_$case_name.json", "w") do io
    JSON.print(io, computed_mse)
end

@testset "life_cycle_Tan2018" begin
    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    include(joinpath("utils", "post_run_tests.jl"))
    nothing
end
