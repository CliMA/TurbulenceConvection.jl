if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using Test

const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
import .NameList

case_name = "GATE_III"

best_mse = Dict()

println("Running $case_name...")
namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "01"
ds_tc_filename, return_code = @time main(namelist; time_run = true)

@testset "GATE_III" begin
    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    include(joinpath(tc_dir, "post_processing", "post_run_tests.jl"))
    nothing
end
