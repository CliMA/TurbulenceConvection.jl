if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList
using .ParamList

best_mse = OrderedDict()

best_mse["updraft_thetal"] = 5.0418613105717345e+00
best_mse["v_mean"] = 5.2842359224499695e+00
best_mse["u_mean"] = 9.7164923816186999e+00
best_mse["tke_mean"] = 4.1067897074868842e+00

@testset "GABLS" begin
    println("Running GABLS...")
    namelist = default_namelist("GABLS")
    paramlist = default_paramlist("GABLS")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

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
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end
