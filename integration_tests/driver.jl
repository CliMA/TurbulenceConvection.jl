import ArgParse

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--case"
        help = "Case to run"
        arg_type = String
        default = "Bomex"
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    case_name = parsed_args["case"]
end
if @isdefined case_name
    @info "Running $case_name..."
else
    case_name = "Bomex"
    @info "Running default case ($case_name)."
    @info "Set `case_name` if you'd like to run a different case"
end

import TurbulenceConvection
using Test

const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
import .NameList

best_mse = all_best_mse[case_name]

namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "01"
ds_tc_filename, return_code = main(namelist)

# Post-processing case kwargs
include(joinpath(tc_dir, "post_processing", "case_kwargs.jl"))

computed_mse = compute_mse_wrapper(case_name, best_mse, ds_tc_filename; case_kwargs[case_name]...)

open("computed_mse_$case_name.json", "w") do io
    JSON.print(io, computed_mse)
end

@testset "$case_name" begin
    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    include(joinpath(tc_dir, "post_processing", "post_run_tests.jl"))
    nothing
end
