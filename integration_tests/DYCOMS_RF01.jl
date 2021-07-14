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
best_mse["qt_mean"] = 3.3444957379390586e-02
best_mse["updraft_area"] = 2.5801657716690602e+01
best_mse["updraft_w"] = 5.3015266139161090e+00
best_mse["updraft_qt"] = 1.9044631066492126e+00
best_mse["updraft_thetal"] = 4.6206573044439402e+01
best_mse["v_mean"] = 1.9191349011646031e+04
best_mse["u_mean"] = 7.3347017004521811e+04
best_mse["tke_mean"] = 2.2114155090507843e+01

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "DYCOMS_RF01.nc"), "r")

@testset "DYCOMS_RF01" begin
    println("Running DYCOMS_RF01...")
    namelist = NameList.DYCOMS_RF01(default_namelist("DYCOMS_RF01"))
    paramlist = ParamList.DYCOMS_RF01(default_paramlist("DYCOMS_RF01"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "DYCOMS_RF01",
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

