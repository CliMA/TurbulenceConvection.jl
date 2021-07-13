import TurbulenceConvection
using TurbulenceConvection
using Test
using Random

# Make deterministic:
Random.seed!(1234)

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

@testset "GATE_III" begin
    println("Running GATE_III...")
    namelist = NameList.GATE_III(default_namelist("GATE_III"))
    paramlist = ParamList.GATE_III(default_paramlist("GATE_III"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
