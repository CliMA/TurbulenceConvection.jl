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
best_mse["updraft_area"] = 6.2270686101489321e+02
best_mse["updraft_w"] = 3.9197923601593011e+01
best_mse["updraft_thetal"] = 2.9892400636969811e+01
best_mse["u_mean"] = 7.6422297657642343e+02
best_mse["tke_mean"] = 6.9178066709527812e+01

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
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    print_artifact_file("Nieuwstadt")
    nothing
end
