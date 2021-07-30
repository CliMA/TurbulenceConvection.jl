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

@testset "SP" begin
    println("Running SP...")
    namelist = default_namelist("SP")
    namelist["meta"]["uuid"] = "01"
    @time main(namelist)
    print_artifact_file("SP")
end
