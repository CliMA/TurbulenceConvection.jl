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
best_mse["qt_mean"] = 1.6028016672114283e-02
best_mse["ql_mean"] = 5.9785694169946080e+00
best_mse["updraft_area"] = 2.3855918966457915e+02
best_mse["updraft_w"] = 4.3334326748232748e+00
best_mse["updraft_qt"] = 1.1941194458793654e+00
best_mse["updraft_thetal"] = 1.2740391229149941e+01
best_mse["v_mean"] = 3.9737974530858637e+01
best_mse["u_mean"] = 3.7038316285973458e+01
best_mse["tke_mean"] = 1.4893837285777673e+01
best_mse["temperature_mean"] = 1.6663829622999760e-05
best_mse["thetal_mean"] = 1.7779259864427720e-05
best_mse["Hvar_mean"] = 1.1587461866566149e+04
best_mse["QTvar_mean"] = 7.3381590434030170e+02

@testset "DYCOMS_RF01" begin
    case_name = "DYCOMS_RF01"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 2 * 3600,
        t_stop = 4 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
