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
best_mse["qt_mean"] = 1.5987989288888820e-02
best_mse["ql_mean"] = 5.9262500323990217e+00
best_mse["updraft_area"] = 2.3887443639936907e+02
best_mse["updraft_w"] = 4.3317426538007693e+00
best_mse["updraft_qt"] = 1.1983341214387264e+00
best_mse["updraft_thetal"] = 1.2739627177179765e+01
best_mse["v_mean"] = 3.9736878190256498e+01
best_mse["u_mean"] = 3.7037227839730683e+01
best_mse["tke_mean"] = 1.4955295527017412e+01
best_mse["temperature_mean"] = 1.6122517253407575e-05
best_mse["thetal_mean"] = 1.7630642054888154e-05
best_mse["Hvar_mean"] = 1.2452017842800817e+04
best_mse["QTvar_mean"] = 7.6626218978684426e+02

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
