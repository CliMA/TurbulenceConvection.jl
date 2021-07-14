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
best_mse["qt_mean"] = 2.5102418088986439e-01
best_mse["updraft_area"] = 2.0267630391221073e+02
best_mse["updraft_w"] = 1.6282411683570519e+01
best_mse["updraft_qt"] = 2.1825912039857140e+01
best_mse["updraft_thetal"] = 6.2489182389073072e+01
best_mse["v_mean"] = 1.9736543073311742e+02
best_mse["u_mean"] = 5.3385203089809568e+02
best_mse["tke_mean"] = 1.7549619573373798e+02

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "Soares.nc"), "r")

@testset "Soares" begin
    println("Running Soares...")
    namelist = NameList.Soares(default_namelist("Soares"))
    paramlist = ParamList.Soares(default_paramlist("Soares"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "Soares",
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

