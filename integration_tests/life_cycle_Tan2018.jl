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
best_mse["qt_mean"] = 5.2649429859732335e-03
best_mse["ql_mean"] = 8.3701130817214975e-01
best_mse["updraft_area"] = 7.0432677991444692e-01
best_mse["updraft_w"] = 5.8558890484998805e-01
best_mse["updraft_qt"] = 1.1615774284759468e-01
best_mse["updraft_thetal"] = 6.2885863821116435e-05
best_mse["v_mean"] = 2.4748668316753225e-01
best_mse["u_mean"] = 7.1361729071471190e-04
best_mse["tke_mean"] = 2.0664613818639557e-01
best_mse["temperature_mean"] = 2.5719773912461363e-06
best_mse["thetal_mean"] = 2.4566431564965449e-06
best_mse["Hvar_mean"] = 2.1515252048343000e+03
best_mse["QTvar_mean"] = 1.1458013034475746e+03

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
