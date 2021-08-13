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
best_mse["qt_mean"] = 9.7923185944396543e-02
best_mse["updraft_area"] = 6.9825342418932712e+02
best_mse["updraft_w"] = 3.2817058320329416e+01
best_mse["updraft_qt"] = 4.2756945036338720e+00
best_mse["updraft_thetal"] = 2.1546731002204076e+01
best_mse["v_mean"] = 6.8320914112603603e+01
best_mse["u_mean"] = 5.3308019185027945e+01
best_mse["tke_mean"] = 4.2619460317296351e+01
best_mse["temperature_mean"] = 4.2264960813453363e-05
best_mse["ql_mean"] = 6.1078690857394591e+00
best_mse["thetal_mean"] = 4.3075669150884208e-05
best_mse["Hvar_mean"] = 1.4320193838595969e+03
best_mse["QTvar_mean"] = 7.4620132655591669e+02

@testset "Bomex" begin
    case_name = "Bomex"
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
