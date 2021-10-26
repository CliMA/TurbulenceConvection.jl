import TurbulenceConvection
using Test

using Profile

include(joinpath("integration_tests", "utils", "Cases.jl"))
include(joinpath("integration_tests", "utils", "generate_namelist.jl"))
using .Cases
import .NameList

include(joinpath("integration_tests", "utils", "main.jl"))

function run_main()
    namelist = default_namelist("GATE_III")
    namelist["meta"]["uuid"] = "01"
    main(namelist)
end

run_main() # run first to compile
@profile run_main()
# Profile.print()
Profile.print(; format = :flat, sortedby = :count)
