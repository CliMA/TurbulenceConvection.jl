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

best_mse["updraft_thetal"] = 5.0248694008424222e+00
best_mse["v_mean"] = 4.4630163253493214e+00
best_mse["u_mean"] = 9.6412500983148846e+00
best_mse["tke_mean"] = 2.4683562232142533e+00
best_mse["temperature_mean"] = 8.8678308948818863e-06
best_mse["thetal_mean"] = 8.7900622849842488e-06
best_mse["Hvar_mean"] = 1.2891667586703637e+01
best_mse["QTvar_mean"] = 4.4424755799649029e-01

@testset "GABLS" begin
    println("Running GABLS...")
    namelist = default_namelist("GABLS")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Gabls.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "GABLS.nc"), "r") do ds_scampy
                compute_mse(
                    "GABLS",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 7 * 3600,
                    t_stop = 9 * 3600,
                )
            end
        end
    end

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
