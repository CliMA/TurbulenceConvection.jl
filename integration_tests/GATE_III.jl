if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
import .NameList

# TODO - waiting for a better root solver
CLIMAParameters.Planet.T_freeze(::EarthParameterSet) = 100.0

@testset "GATE_III" begin
    println("Running GATE_III...")
    namelist = default_namelist("GATE_III")
    namelist["meta"]["uuid"] = "01"
    @time main(namelist)
end
