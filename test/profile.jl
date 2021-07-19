import TurbulenceConvection
using TurbulenceConvection
using Test

using Profile

include(joinpath("integration_tests", "utils", "Cases.jl"))
include(joinpath("integration_tests", "utils", "generate_paramlist.jl"))
include(joinpath("integration_tests", "utils", "generate_namelist.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("integration_tests", "utils", "main.jl"))

function run_main()
    namelist = NameList.GATE_III(default_namelist("GATE_III"))
    paramlist = ParamList.GATE_III(default_paramlist("GATE_III"))
    namelist["meta"]["uuid"] = "01"
    main(namelist, paramlist)
end

run_main() # run first to compile
@profile run_main()
# Profile.print()
Profile.print(; format = :flat, sortedby = :count)
