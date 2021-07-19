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
best_mse["qt_mean"] = 1.1424341230332884e-01
best_mse["updraft_area"] = 1.1765062251021085e+02
best_mse["updraft_w"] = 6.9115336258511732e+01
best_mse["updraft_qt"] = 6.6796467668487018e+00
best_mse["updraft_thetal"] = 7.2026005353910350e+01
best_mse["v_mean"] = 1.2520237439186344e+02
best_mse["u_mean"] = 2.2785293881438179e+03
best_mse["tke_mean"] = 4.1588123897851069e+01

@testset "Nieuwstadt" begin
    println("Running Nieuwstadt...")
    namelist = NameList.Nieuwstadt(default_namelist("Nieuwstadt"))
    paramlist = ParamList.Nieuwstadt(default_paramlist("Nieuwstadt"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Nieuwstadt.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Nieuwstadt.nc"), "r") do ds_scampy
                compute_mse(
                    "Nieuwstadt",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv=ds,
                    ds_scampy=ds_scampy,
                    ds_pycles=ds_pycles,
                    plot_comparison=true
                )
            end
        end
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

