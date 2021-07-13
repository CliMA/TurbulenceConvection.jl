import TurbulenceConvection
using TurbulenceConvection
using Test
using Random

# Make deterministic:
Random.seed!(1234)

using Profile

include(joinpath("integration_tests", "utils", "Cases.jl"))
include(joinpath("integration_tests", "utils", "generate_paramlist.jl"))
include(joinpath("integration_tests", "utils", "generate_namelist.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("integration_tests", "utils", "main.jl"))

function run_main(; time_run=false)
    namelist = NameList.Bomex(default_namelist("Bomex"))
    paramlist = ParamList.Bomex(default_paramlist("Bomex"))
    namelist["meta"]["uuid"] = "01"
    main(namelist, paramlist; time_run=time_run)
end

run_main() # run first to compile
# run_main(;time_run=true)
