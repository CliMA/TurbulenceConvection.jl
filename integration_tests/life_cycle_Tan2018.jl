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
best_mse["qt_mean"] = 5.2649432939457785e-03
best_mse["ql_mean"] = 8.3701044168155769e-01
best_mse["updraft_area"] = 7.0432678120046099e-01
best_mse["updraft_w"] = 5.8558890544083542e-01
best_mse["updraft_qt"] = 1.1615774393933681e-01
best_mse["updraft_thetal"] = 6.2885864461831401e-05
best_mse["v_mean"] = 2.4748668041280980e-01
best_mse["u_mean"] = 7.1361727014691585e-04
best_mse["tke_mean"] = 2.0664613598760551e-01
best_mse["temperature_mean"] = 2.5719775476904971e-06
best_mse["thetal_mean"] = 2.4566433106218294e-06
best_mse["Hvar_mean"] = 2.1515325662209743e+03
best_mse["QTvar_mean"] = 1.1458060992053790e+03

@testset "life_cycle_Tan2018" begin
    case_name = "life_cycle_Tan2018"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 4 * 3600,
        t_stop = 6 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
