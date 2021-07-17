if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

@testset "LES_driven_SCM" begin
    println("Running LES_driven_SCM...")
    namelist = NameList.LES_driven_SCM(default_namelist("LES_driven_SCM"))
    paramlist = ParamList.LES_driven_SCM(default_paramlist("LES_driven_SCM"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
