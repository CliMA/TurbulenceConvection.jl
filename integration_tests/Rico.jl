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
best_mse["qt_mean"] = 4.2745465446436787e-01
best_mse["updraft_area"] = 1.5810580213880007e+03
best_mse["updraft_w"] = 3.0293830804824944e+02
best_mse["updraft_qt"] = 1.8036840808532091e+01
best_mse["updraft_thetal"] = 6.3565923109974378e+01
best_mse["v_mean"] = 1.0637284782337001e+02
best_mse["u_mean"] = 1.1444354015737491e+02
best_mse["tke_mean"] = 8.6553999307149547e+02
best_mse["temperature_mean"] = 1.2540517259642333e-04
best_mse["ql_mean"] = 5.1325163153016852e+01
best_mse["thetal_mean"] = 1.2561593164473264e-04
best_mse["Hvar_mean"] = 1.4739277424672016e+03
best_mse["QTvar_mean"] = 7.1366176286414179e+02

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

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
