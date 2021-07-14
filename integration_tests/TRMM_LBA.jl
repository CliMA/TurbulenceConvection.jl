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
best_mse["qt_mean"] = 7.5561281422031961e+00
best_mse["updraft_area"] = 8.3276151294410829e+03
best_mse["updraft_w"] = 3.7727335913983225e+03
best_mse["updraft_qt"] = 6.2330579729990070e+01
best_mse["updraft_thetal"] = 1.9862147491681133e+03
best_mse["v_mean"] = 8.7742101712067142e+02
best_mse["u_mean"] = 1.0046223117017050e+03
best_mse["tke_mean"] = 6.1122707931388168e+03

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "TRMM_LBA.nc"), "r")

@testset "TRMM_LBA" begin
    println("Running TRMM_LBA...")
    namelist = NameList.TRMM_LBA(default_namelist("TRMM_LBA"))
    paramlist = ParamList.TRMM_LBA(default_paramlist("TRMM_LBA"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "TRMM_LBA",
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

