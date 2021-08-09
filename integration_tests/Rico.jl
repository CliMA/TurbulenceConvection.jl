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
best_mse["qt_mean"] = 4.2359420114784679e-01
best_mse["updraft_area"] = 1.6722742675049749e+03
best_mse["updraft_w"] = 2.8300484466673169e+02
best_mse["updraft_qt"] = 1.7904453430368505e+01
best_mse["updraft_thetal"] = 6.3573804437983789e+01
best_mse["v_mean"] = 1.0634096403180666e+02
best_mse["u_mean"] = 1.1434656097951681e+02
best_mse["tke_mean"] = 8.6824438228445820e+02
best_mse["temperature_mean"] = 1.2576042829475203e-04
best_mse["ql_mean"] = 8.9402683009257942e+01
best_mse["thetal_mean"] = 1.1896802334339157e-04
best_mse["Hvar_mean"] = 2.0228054123757493e+03
best_mse["QTvar_mean"] = 8.7748179858993274e+02

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
