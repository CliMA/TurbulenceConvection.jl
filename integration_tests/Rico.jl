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
best_mse["qt_mean"] = 4.0432375728696790e-01
best_mse["updraft_area"] = 1.8875614747005029e+03
best_mse["updraft_w"] = 4.7945553218441751e+02
best_mse["updraft_qt"] = 2.3512519674647539e+01
best_mse["updraft_thetal"] = 6.5601910609184642e+01
best_mse["v_mean"] = 1.0614382605973881e+02
best_mse["u_mean"] = 1.1383366513167900e+02
best_mse["tke_mean"] = 8.9241612675577824e+02

@testset "Rico" begin
    println("Running Rico...")
    namelist = NameList.Rico(default_namelist("Rico"))
    paramlist = ParamList.Rico(default_paramlist("Rico"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Rico.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Rico.nc"), "r") do ds_scampy
                compute_mse(
                    "Rico",
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

