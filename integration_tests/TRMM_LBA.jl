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
best_mse["qt_mean"] = 3.8877970266171977e+00
best_mse["updraft_area"] = 3.1640134647479666e+04
best_mse["updraft_w"] = 8.0040687916282764e+02
best_mse["updraft_qt"] = 1.1943729301590135e+01
best_mse["updraft_thetal"] = 1.1110312052943021e+02
best_mse["v_mean"] = 2.9891844806438138e+02
best_mse["u_mean"] = 1.6669521624928782e+03
best_mse["tke_mean"] = 2.6933890983728834e+03

@testset "TRMM_LBA" begin
    println("Running TRMM_LBA...")
    namelist = NameList.TRMM_LBA(default_namelist("TRMM_LBA"))
    paramlist = ParamList.TRMM_LBA(default_paramlist("TRMM_LBA"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_scampy
                compute_mse(
                    "TRMM_LBA",
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

