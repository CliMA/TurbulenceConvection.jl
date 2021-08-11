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
best_mse["qt_mean"] = 9.8858939409301599e-02
best_mse["updraft_area"] = 7.2423270403289905e+02
best_mse["updraft_w"] = 2.7297489157558520e+01
best_mse["updraft_qt"] = 4.0686060603282321e+00
best_mse["updraft_thetal"] = 2.1547376374778047e+01
best_mse["v_mean"] = 6.5159453481189033e+01
best_mse["u_mean"] = 5.3292347667845092e+01
best_mse["tke_mean"] = 3.8707689709123862e+01
best_mse["temperature_mean"] = 4.0396510560887023e-05
best_mse["ql_mean"] = 3.9678285547654855e+00
best_mse["thetal_mean"] = 4.1130421765389724e-05
best_mse["Hvar_mean"] = 5.9914307236483694e+01
best_mse["QTvar_mean"] = 1.9732739289298404e+01

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = default_namelist("Bomex")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds_tc
        Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Bomex.nc"), "r") do ds_scampy
                compute_mse(
                    "Bomex",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_tc = ds_tc,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 4 * 3600,
                    t_stop = 6 * 3600,
                )
            end
        end
    end

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
