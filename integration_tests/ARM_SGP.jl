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
include(joinpath("utils", "compute_mse.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

best_mse = OrderedDict()
best_mse["qt_mean"] = 5.1999754372401501e+00
best_mse["updraft_area"] = 4.7197978700467644e+01
best_mse["updraft_w"] = 6.6472456347075465e-02
best_mse["updraft_qt"] = 1.0454097540672662e+02
best_mse["updraft_thetal"] = 4.1837835135772195e+02
best_mse["v_mean"] = 8.5596723994911804e+02
best_mse["u_mean"] = 5.7620127283873725e+08
best_mse["tke_mean"] = 1.0150322806956440e+00

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "ARM_SGP.nc"), "r")

@testset "ARM_SGP" begin
    println("Running ARM_SGP...")
    namelist = NameList.ARM_SGP(default_namelist("ARM_SGP"))
    paramlist = ParamList.ARM_SGP(default_paramlist("ARM_SGP"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "ARM_SGP",
            best_mse,
            dirname(ds_filename);
            plot_comparison=true
        )
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end

