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
best_mse["qt_mean"] = 3.9046069426558216e+00
best_mse["updraft_area"] = 2.8331099674735910e+04
best_mse["updraft_w"] = 6.5629758352987142e+02
best_mse["updraft_qt"] = 2.5012703977470046e+01
best_mse["updraft_thetal"] = 1.1072811439926345e+02
best_mse["v_mean"] = 2.9406070727525389e+02
best_mse["u_mean"] = 1.6903159369983448e+03
best_mse["tke_mean"] = 2.8961135369556791e+03
best_mse["temperature_mean"] = 8.3946608562868500e-04
best_mse["ql_mean"] = 4.3631392987274085e+03
best_mse["thetal_mean"] = 5.7967137448949413e-04
best_mse["Hvar_mean"] = 1.3320049142979460e+04
best_mse["QTvar_mean"] = 1.5870164338554601e+04

@testset "TRMM_LBA" begin
    println("Running TRMM_LBA...")
    namelist = default_namelist("TRMM_LBA")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_scampy
                compute_mse(
                    "TRMM_LBA",
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
    nothing
end
