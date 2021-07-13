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

@testset "DYCOMS_RF01" begin
    println("Running DYCOMS_RF01...")
    namelist = NameList.DYCOMS_RF01(default_namelist("DYCOMS_RF01"))
    paramlist = ParamList.DYCOMS_RF01(default_paramlist("DYCOMS_RF01"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
