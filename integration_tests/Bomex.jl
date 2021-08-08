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
best_mse["qt_mean"] = 8.9644597356208541e-02
best_mse["updraft_area"] = 6.9327817412026047e+02
best_mse["updraft_w"] = 8.9057089858539797e+01
best_mse["updraft_qt"] = 5.9179688063528078e+00
best_mse["updraft_thetal"] = 2.2933264810554988e+01
best_mse["v_mean"] = 1.2348440297549568e+02
best_mse["u_mean"] = 5.3484778922733959e+01
best_mse["tke_mean"] = 3.2847767205446033e+01
best_mse["temperature_mean"] = 3.4874387520736297e-05
best_mse["ql_mean"] = 2.1344864368686135e+01
best_mse["thetal_mean"] = 3.5327188090292184e-05
best_mse["Hvar_mean"] = 1.7623971358198961e+02
best_mse["QTvar_mean"] = 4.0583986764371645e+01

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
