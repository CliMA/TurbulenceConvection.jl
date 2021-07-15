if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

best_mse = OrderedDict()
best_mse["qt_mean"] = 1.0701019653530450e-01
best_mse["updraft_area"] = 6.9127918666220057e+02
best_mse["updraft_w"] = 8.3196631916073699e+01
best_mse["updraft_qt"] = 5.9300799798912323e+00
best_mse["updraft_thetal"] = 2.3063368634456552e+01
best_mse["v_mean"] = 1.2992861032456707e+02
best_mse["u_mean"] = 5.3536353327695572e+01
best_mse["tke_mean"] = 3.3193072857888232e+01

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = NameList.Bomex(default_namelist("Bomex"))
    paramlist = ParamList.Bomex(default_paramlist("Bomex"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Bomex.nc"), "r") do ds_scampy
                compute_mse(
                    "Bomex",
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

