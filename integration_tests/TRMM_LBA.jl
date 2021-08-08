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
best_mse["qt_mean"] = 3.0609461714234958e+00
best_mse["updraft_area"] = 3.5992222800156262e+04
best_mse["updraft_w"] = 1.0100451430716917e+03
best_mse["updraft_qt"] = 2.7507151910663865e+01
best_mse["updraft_thetal"] = 1.1074070908480367e+02
best_mse["v_mean"] = 2.9438922611828133e+02
best_mse["u_mean"] = 1.6906200771006561e+03
best_mse["tke_mean"] = 2.5854592261914117e+03
best_mse["temperature_mean"] = 9.2883217582418853e-04
best_mse["ql_mean"] = 2.1864394985446243e+03
best_mse["thetal_mean"] = 6.4455936975381660e-04
best_mse["Hvar_mean"] = 1.2394474700136290e+04
best_mse["QTvar_mean"] = 9.2129571715440597e+03

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
