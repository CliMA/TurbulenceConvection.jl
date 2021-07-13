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

@testset "GABLS" begin
    println("Running GABLS...")
    namelist = NameList.GABLS(default_namelist("GABLS"))
    paramlist = ParamList.GABLS(default_paramlist("GABLS"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
