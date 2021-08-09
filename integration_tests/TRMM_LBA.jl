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

# Note: temperatures in this case become extremely low.
CLIMAParameters.Planet.T_freeze(::EarthParameterSet) = 100.0

best_mse = OrderedDict()
best_mse["qt_mean"] = 1.7995063250630168e+00
best_mse["updraft_area"] = 3.1529411391564761e+04
best_mse["updraft_w"] = 1.0325943376733990e+03
best_mse["updraft_qt"] = 3.1156138180627952e+01
best_mse["updraft_thetal"] = 1.1001893832675252e+02
best_mse["v_mean"] = 2.9275592916609139e+02
best_mse["u_mean"] = 1.6873153880206969e+03
best_mse["tke_mean"] = 1.3627971247685127e+03
best_mse["temperature_mean"] = 6.8597416865395539e-04
best_mse["ql_mean"] = 1.0669889889408917e+03
best_mse["thetal_mean"] = 8.1757856795585861e-03
best_mse["Hvar_mean"] = 6.8922213953686887e+03
best_mse["QTvar_mean"] = 2.5975685843697461e+03

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
                    t_start = 4 * 3600,
                    t_stop = 6 * 3600,
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
