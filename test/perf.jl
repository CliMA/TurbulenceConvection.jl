import TurbulenceConvection
using TurbulenceConvection
using Test

using Profile

include(joinpath("integration_tests", "utils", "Cases.jl"))
include(joinpath("integration_tests", "utils", "generate_namelist.jl"))
using .Cases
using .NameList

include(joinpath("integration_tests", "utils", "main.jl"))

function run_main(; time_run = false)
    namelist = default_namelist("Bomex")
    namelist["meta"]["uuid"] = "01"
    main(namelist; time_run = time_run)
end

run_main() # run first to compile
# run_main(;time_run=true)
