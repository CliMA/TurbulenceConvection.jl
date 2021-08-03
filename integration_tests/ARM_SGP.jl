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
best_mse["qt_mean"] = 5.8920349209467660e-01
best_mse["updraft_area"] = 2.9161664375922396e+02
best_mse["updraft_w"] = 9.2064966152387040e-01
best_mse["updraft_qt"] = 5.2856934357682029e+01
best_mse["updraft_thetal"] = 6.9072600122408147e+01
best_mse["u_mean"] = 8.7994360629094189e+01
best_mse["tke_mean"] = 2.8644243409080627e+00
best_mse["temperature_mean"] = 2.7496219242928010e-04
best_mse["ql_mean"] = 5.0850007235659434e-03
best_mse["thetal_mean"] = 2.7328104527915207e-04
best_mse["Hvar_mean"] = 2.2244310672644488e+01
best_mse["QTvar_mean"] = 2.6090373899600383e+01

@testset "ARM_SGP" begin
    println("Running ARM_SGP...")
    namelist = default_namelist("ARM_SGP")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "ARM_SGP.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "ARM_SGP.nc"), "r") do ds_scampy
                compute_mse(
                    "ARM_SGP",
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
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    test_mse(computed_mse, best_mse, "temperature_mean")
    test_mse(computed_mse, best_mse, "ql_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    test_mse(computed_mse, best_mse, "Hvar_mean")
    test_mse(computed_mse, best_mse, "QTvar_mean")
    nothing
end
