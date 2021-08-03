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
best_mse["qt_mean"] = 8.4924728795344143e-02
best_mse["updraft_area"] = 7.0195368138769777e+02
best_mse["updraft_w"] = 8.7703050153771670e+01
best_mse["updraft_qt"] = 5.9194060591321715e+00
best_mse["updraft_thetal"] = 2.2932547498112353e+01
best_mse["v_mean"] = 1.2350993295040250e+02
best_mse["u_mean"] = 5.3487467706554988e+01
best_mse["tke_mean"] = 3.2337733855931930e+01
best_mse["temperature_mean"] = 3.1767291118329973e-05
best_mse["ql_mean"] = 1.4250971600048450e+01
best_mse["thetal_mean"] = 3.2420072527836003e-05
best_mse["Hvar_mean"] = 6.2697869671981366e+01
best_mse["QTvar_mean"] = 2.2154880539077912e+01

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = default_namelist("Bomex")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

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
    test_mse(computed_mse, best_mse, "temperature_mean")
    test_mse(computed_mse, best_mse, "ql_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    nothing
end
