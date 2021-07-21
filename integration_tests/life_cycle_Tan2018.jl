if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList
using .ParamList

@testset "life_cycle_Tan2018" begin
    println("Running life_cycle_Tan2018...")
    namelist = default_namelist("life_cycle_Tan2018")
    paramlist = default_paramlist("life_cycle_Tan2018")
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
