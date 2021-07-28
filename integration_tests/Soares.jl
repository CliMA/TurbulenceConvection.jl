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

best_mse = OrderedDict()
best_mse["qt_mean"] = 2.5843807856554074e-01
best_mse["updraft_area"] = 8.8039804794002623e+02
best_mse["updraft_w"] = 2.6632875471066519e+01
best_mse["updraft_qt"] = 1.0651734450298703e+01
best_mse["updraft_thetal"] = 2.1622611589959376e+01
best_mse["u_mean"] = 3.9550058876376834e+03
best_mse["tke_mean"] = 5.8683619102955831e+01

@testset "Soares" begin
    println("Running Soares...")
    namelist = default_namelist("Soares")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Soares.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Soares.nc"), "r") do ds_scampy
                compute_mse(
                    "Soares",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end
