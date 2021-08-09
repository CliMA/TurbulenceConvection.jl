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
best_mse["qt_mean"] = 4.9523533702786315e-01
best_mse["updraft_area"] = 1.7419532961183722e+03
best_mse["updraft_w"] = 5.3893010779565452e+02
best_mse["updraft_qt"] = 2.4530572533274270e+01
best_mse["updraft_thetal"] = 6.5509509302767754e+01
best_mse["v_mean"] = 1.0615326468150907e+02
best_mse["u_mean"] = 1.1393191702720974e+02
best_mse["tke_mean"] = 8.4090034387266689e+02
best_mse["temperature_mean"] = 1.4274314152822176e-04
best_mse["ql_mean"] = 2.1441653002617767e+02
best_mse["thetal_mean"] = 1.3609606969128344e-04
best_mse["Hvar_mean"] = 3.7305257502890240e+03
best_mse["QTvar_mean"] = 1.2241261056877631e+03

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
                    t_start = 22 * 3600,
                    t_stop = 24 * 3600,
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
