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
best_mse["updraft_area"] = 5.9103576007452034e+02
best_mse["updraft_w"] = 2.6741082264407293e+01
best_mse["updraft_thetal"] = 3.0475272117346940e+01
best_mse["u_mean"] = 1.5245873874565478e+02
best_mse["tke_mean"] = 7.3707728452173882e+01
best_mse["temperature_mean"] = 1.2199914281992930e-05
best_mse["thetal_mean"] = 1.2067109520894632e-05
best_mse["Hvar_mean"] = 1.9140690462510639e+02

@testset "Nieuwstadt" begin
    println("Running Nieuwstadt...")
    namelist = default_namelist("Nieuwstadt")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Nieuwstadt.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Nieuwstadt.nc"), "r") do ds_scampy
                compute_mse(
                    "Nieuwstadt",
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

    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    test_mse(computed_mse, best_mse, "temperature_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    test_mse(computed_mse, best_mse, "Hvar_mean")
    nothing
end
