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
best_mse["qt_mean"] = 4.8938409872334176e-01
best_mse["updraft_area"] = 1.9006975794173386e+03
best_mse["updraft_w"] = 4.9377041328364481e+02
best_mse["updraft_qt"] = 2.3571647470794897e+01
best_mse["updraft_thetal"] = 6.5599302555401209e+01
best_mse["v_mean"] = 1.0614722053614342e+02
best_mse["u_mean"] = 1.1381848101828112e+02
best_mse["tke_mean"] = 8.7760700153879964e+02

@testset "Rico" begin
    println("Running Rico...")
    namelist = default_namelist("Rico")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Rico.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Rico.nc"), "r") do ds_scampy
                compute_mse(
                    "Rico",
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
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end
