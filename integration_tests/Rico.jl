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
best_mse["qt_mean"] = 5.0707745192217146e-01
best_mse["updraft_area"] = 1.6283285530615315e+03
best_mse["updraft_w"] = 5.7196482018497807e+02
best_mse["updraft_qt"] = 2.4565269303654755e+01
best_mse["updraft_thetal"] = 6.5505090179439435e+01
best_mse["v_mean"] = 1.0611825566231954e+02
best_mse["u_mean"] = 1.1390377464452069e+02
best_mse["tke_mean"] = 8.9525123248795580e+02
best_mse["temperature_mean"] = 1.4787364980827869e-04
best_mse["ql_mean"] = 1.7937532436604911e+02
best_mse["thetal_mean"] = 1.4098220595492702e-04
best_mse["Hvar_mean"] = 1.6499437998682095e+03
best_mse["QTvar_mean"] = 7.3181087637066230e+02

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
    test_mse(computed_mse, best_mse, "temperature_mean")
    test_mse(computed_mse, best_mse, "ql_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    test_mse(computed_mse, best_mse, "Hvar_mean")
    test_mse(computed_mse, best_mse, "QTvar_mean")
    nothing
end
