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
best_mse["qt_mean"] = 1.1424341230332884e-01
best_mse["updraft_area"] = 1.1765062251021085e+02
best_mse["updraft_w"] = 6.9115336258511732e+01
best_mse["updraft_qt"] = 6.6796467668487018e+00
best_mse["updraft_thetal"] = 7.2026005353910350e+01
best_mse["v_mean"] = 1.2520237439186344e+02
best_mse["u_mean"] = 2.2785293881438179e+03
best_mse["tke_mean"] = 4.1588123897851069e+01

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r")

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = NameList.Bomex(default_namelist("Bomex"))
    paramlist = ParamList.Bomex(default_paramlist("Bomex"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "Bomex",
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

