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
best_mse["qt_mean"] = 1.1113745694252045e-01
best_mse["updraft_area"] = 7.2538205557139236e+02
best_mse["updraft_w"] = 2.5943804836195852e+01
best_mse["updraft_qt"] = 4.0075397466347100e+00
best_mse["updraft_thetal"] = 2.1549266568293614e+01
best_mse["v_mean"] = 6.5506192741355221e+01
best_mse["u_mean"] = 5.3295866692805944e+01
best_mse["tke_mean"] = 3.9596900219869937e+01
best_mse["temperature_mean"] = 4.6114814306316929e-05
best_mse["ql_mean"] = 6.7039710795893939e+00
best_mse["thetal_mean"] = 4.6514684598969401e-05
best_mse["Hvar_mean"] = 2.5387994922244559e+02
best_mse["QTvar_mean"] = 4.5843881812484575e+01

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = default_namelist("Bomex")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Bomex.nc"), "r") do ds_scampy
                compute_mse(
                    "Bomex",
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
    test_mse(computed_mse, best_mse, "ql_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    nothing
end
