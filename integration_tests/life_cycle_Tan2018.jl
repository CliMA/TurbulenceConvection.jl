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

@testset "life_cycle_Tan2018" begin
    println("Running life_cycle_Tan2018...")
    namelist = NameList.life_cycle_Tan2018(default_namelist("life_cycle_Tan2018"))
    paramlist = ParamList.life_cycle_Tan2018(default_paramlist("life_cycle_Tan2018"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
