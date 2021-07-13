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

@testset "Rico" begin
    println("Running Rico...")
    namelist = NameList.Rico(default_namelist("Rico"))
    paramlist = ParamList.Rico(default_paramlist("Rico"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
