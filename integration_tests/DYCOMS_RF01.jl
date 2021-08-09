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
best_mse["qt_mean"] = 1.5987989356173508e-02
best_mse["ql_mean"] = 5.9262542639081852e+00
best_mse["updraft_area"] = 2.3887443993964089e+02
best_mse["updraft_w"] = 4.3317357122037770e+00
best_mse["updraft_qt"] = 1.1983341056077728e+00
best_mse["updraft_thetal"] = 1.2739627176096105e+01
best_mse["v_mean"] = 3.9736878160027324e+01
best_mse["u_mean"] = 3.7037227811674803e+01
best_mse["tke_mean"] = 1.4955313626928097e+01
best_mse["temperature_mean"] = 1.6122516172306611e-05
best_mse["thetal_mean"] = 1.7630641025752975e-05
best_mse["Hvar_mean"] = 1.2452024821778659e+04
best_mse["QTvar_mean"] = 7.6626261862502452e+02


@testset "DYCOMS_RF01" begin
    println("Running DYCOMS_RF01...")
    namelist = default_namelist("DYCOMS_RF01")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "DYCOMS_RF01.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "DYCOMS_RF01.nc"), "r") do ds_scampy
                compute_mse(
                    "DYCOMS_RF01",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 2 * 3600,
                    t_stop = 4 * 3600,
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
    test_mse(computed_mse, best_mse, "thetal_mean")
    test_mse(computed_mse, best_mse, "Hvar_mean")
    test_mse(computed_mse, best_mse, "QTvar_mean")
    nothing
end
