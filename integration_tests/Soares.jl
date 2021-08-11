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
best_mse["qt_mean"] = 1.9368730493269354e-01
best_mse["updraft_area"] = 4.8515131142491634e+02
best_mse["updraft_w"] = 3.1260485803642322e+01
best_mse["updraft_qt"] = 1.1137263873078020e+01
best_mse["updraft_thetal"] = 2.2395183773363378e+01
best_mse["u_mean"] = 7.2111916666642423e+02
best_mse["tke_mean"] = 4.7882808859338290e+01
best_mse["temperature_mean"] = 1.3818642368852120e-05
best_mse["thetal_mean"] = 1.4383770315305301e-05
best_mse["Hvar_mean"] = 2.3235278783484819e+02

@testset "Soares" begin
    println("Running Soares...")
    namelist = default_namelist("Soares")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds_tc
        Dataset(joinpath(PyCLES_output_dataset_path, "Soares.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Soares.nc"), "r") do ds_scampy
                compute_mse(
                    "Soares",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_tc = ds_tc,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 6 * 3600,
                    t_stop = 8 * 3600,
                )
            end
        end
    end

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
