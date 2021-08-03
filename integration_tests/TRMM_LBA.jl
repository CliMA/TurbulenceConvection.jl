if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList

best_mse = OrderedDict()
best_mse["qt_mean"] = 3.9066401784413203e+00
best_mse["updraft_area"] = 2.8333333642947131e+04
best_mse["updraft_w"] = 6.5662893402745362e+02
best_mse["updraft_qt"] = 2.5017304239563906e+01
best_mse["updraft_thetal"] = 1.1072812435859557e+02
best_mse["v_mean"] = 2.9403959926399466e+02
best_mse["u_mean"] = 1.6903038672404982e+03
best_mse["tke_mean"] = 2.8972134698112682e+03
best_mse["temperature_mean"] = 8.3944025946955543e-04
best_mse["ql_mean"] = 4.3199525166349194e+03
best_mse["thetal_mean"] = 5.7974301596273598e-04
best_mse["Hvar_mean"] = 1.3419697125412104e+04
best_mse["QTvar_mean"] = 1.6238639092892809e+04

@testset "TRMM_LBA" begin
    println("Running TRMM_LBA...")
    namelist = default_namelist("TRMM_LBA")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_scampy
                compute_mse(
                    "TRMM_LBA",
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
    test_mse(computed_mse, best_mse, "temperature_mean")
    nothing
end
