if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
import .NameList

println("Running GATE_III...")
namelist = default_namelist("GATE_III")
namelist["meta"]["uuid"] = "01"
ds_tc_filename = @time main(namelist; time_run = true)

@testset "GATE_III" begin end

include(joinpath("utils", "post_run_tests.jl"))
