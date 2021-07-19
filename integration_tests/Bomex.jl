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
best_mse["qt_mean"] = 8.2625377129815783e-02
best_mse["updraft_area"] = 7.0269263421038499e+02
best_mse["updraft_w"] = 8.7981763090171825e+01
best_mse["updraft_qt"] = 5.9526663775057314e+00
best_mse["updraft_thetal"] = 2.3059913728916317e+01
best_mse["v_mean"] = 1.2340566846095439e+02
best_mse["u_mean"] = 5.3487046322597514e+01
best_mse["tke_mean"] = 3.2160344619556398e+01

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
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
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
