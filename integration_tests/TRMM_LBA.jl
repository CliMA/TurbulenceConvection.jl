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

# Note: temperatures in this case become extremely low.
CLIMAParameters.Planet.T_freeze(::EarthParameterSet) = 100.0

best_mse = OrderedDict()
best_mse["qt_mean"] = 1.6060537595274345e+00
best_mse["updraft_area"] = 2.4326719919434730e+04
best_mse["updraft_w"] = 8.9781354744734904e+02
best_mse["updraft_qt"] = 2.7840588043213661e+01
best_mse["updraft_thetal"] = 1.1000857010242179e+02
best_mse["v_mean"] = 2.9255417406952165e+02
best_mse["u_mean"] = 1.6872488577305819e+03
best_mse["tke_mean"] = 1.5447663031746142e+03
best_mse["temperature_mean"] = 6.8027838562588429e-04
best_mse["ql_mean"] = 9.5570328764458566e+02
best_mse["thetal_mean"] = 8.1813107572736785e-03
best_mse["Hvar_mean"] = 6.7055563585233622e+03
best_mse["QTvar_mean"] = 2.6582568278836857e+03

@testset "TRMM_LBA" begin
    println("Running TRMM_LBA...")
    namelist = default_namelist("TRMM_LBA")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds_tc
        Dataset(joinpath(PyCLES_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "TRMM_LBA.nc"), "r") do ds_scampy
                compute_mse(
                    "TRMM_LBA",
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
