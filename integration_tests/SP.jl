if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
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

@testset "SP" begin
    println("Running SP...")
    namelist = NameList.SP(default_namelist("SP"))
    paramlist = ParamList.SP(default_paramlist("SP"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
