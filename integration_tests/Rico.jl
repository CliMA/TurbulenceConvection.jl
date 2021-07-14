if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test
using Random

# Make deterministic:
Random.seed!(1234)

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

best_mse = OrderedDict()
best_mse["qt_mean"] = 4.0145821221362488e-01
best_mse["updraft_area"] = 4.6525267644195679e+02
best_mse["updraft_w"] = 2.2908548924319413e+02
best_mse["updraft_qt"] = 2.5661494340253533e+01
best_mse["updraft_thetal"] = 4.2409565376016479e+02
best_mse["v_mean"] = 3.6252873079755642e+04
best_mse["u_mean"] = 2.4726879986917771e+02
best_mse["tke_mean"] = 2.0708781629689062e+02

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "Rico.nc"), "r")

@testset "Rico" begin
    println("Running Rico...")
    namelist = NameList.Rico(default_namelist("Rico"))
    paramlist = ParamList.Rico(default_paramlist("Rico"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "Rico",
            best_mse,
            dirname(ds_filename);
            plot_comparison=true
        )
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end

