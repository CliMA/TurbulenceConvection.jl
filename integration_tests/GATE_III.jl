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

@testset "GATE_III" begin
    println("Running GATE_III...")
    namelist = NameList.GATE_III(default_namelist("GATE_III"))
    paramlist = ParamList.GATE_III(default_paramlist("GATE_III"))
    namelist["meta"]["uuid"] = "01"
    @time main(namelist, paramlist)
end
