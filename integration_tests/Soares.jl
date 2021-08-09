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
best_mse["qt_mean"] = 3.2147414017234566e-01
best_mse["updraft_area"] = 7.5523687950242140e+02
best_mse["updraft_w"] = 3.5564080086151499e+01
best_mse["updraft_qt"] = 1.0702812362369221e+01
best_mse["updraft_thetal"] = 2.1534108878970788e+01
best_mse["u_mean"] = 3.9385484811475408e+03
best_mse["tke_mean"] = 5.2022934994314660e+01
best_mse["temperature_mean"] = 1.9333588833166558e-05
best_mse["thetal_mean"] = 1.9486304719984318e-05
best_mse["Hvar_mean"] = 2.6039266998435176e+02


@testset "Soares" begin
    println("Running Soares...")
    namelist = default_namelist("Soares")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Soares.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Soares.nc"), "r") do ds_scampy
                compute_mse(
                    "Soares",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 6 * 3600,
                    t_stop = 8 * 3600,
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    test_mse(computed_mse, best_mse, "temperature_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    test_mse(computed_mse, best_mse, "Hvar_mean")
    nothing
end
