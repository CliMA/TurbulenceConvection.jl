if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList

@testset "DryBubble" begin
    println("Running DryBubble...")
    namelist = default_namelist("DryBubble")
    namelist["meta"]["uuid"] = "01"
    @time main(namelist)
end
