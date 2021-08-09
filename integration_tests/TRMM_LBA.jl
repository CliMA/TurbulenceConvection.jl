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
best_mse["qt_mean"] = 1.7430285968163504e+00
best_mse["updraft_area"] = 3.2545401197014067e+04
best_mse["updraft_w"] = 9.2113782722003907e+02
best_mse["updraft_qt"] = 3.0306781659879494e+01
best_mse["updraft_thetal"] = 1.1108867837608678e+02
best_mse["v_mean"] = 2.9291577336310036e+02
best_mse["u_mean"] = 1.6874237352911189e+03
best_mse["tke_mean"] = 1.5857027462501208e+03
best_mse["temperature_mean"] = 6.3088676948654218e-04
best_mse["ql_mean"] = 9.4722581733362244e+02
best_mse["thetal_mean"] = 4.1952296308447006e-04
best_mse["Hvar_mean"] = 8.1275664794363738e+03
best_mse["QTvar_mean"] = 2.7043370821860535e+03

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
                    t_start = 4 * 3600,
                    t_stop = 6 * 3600,
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
