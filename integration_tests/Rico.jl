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
best_mse["qt_mean"] = 3.5489212020759542e-01
best_mse["updraft_area"] = 1.9189473424865719e+03
best_mse["updraft_w"] = 1.6330042256516663e+02
best_mse["updraft_qt"] = 1.5365813368597401e+01
best_mse["updraft_thetal"] = 6.3601413658374895e+01
best_mse["v_mean"] = 1.0634033027478135e+02
best_mse["u_mean"] = 1.1440791317739802e+02
best_mse["tke_mean"] = 3.4396944042394813e+02
best_mse["temperature_mean"] = 1.7375385006813271e-04
best_mse["ql_mean"] = 1.2777774668393246e+02
best_mse["thetal_mean"] = 1.4818668975836064e-04
best_mse["Hvar_mean"] = 5.4593421490878145e+03
best_mse["QTvar_mean"] = 2.4103920470177113e+03

@testset "Rico" begin
    case_name = "Rico"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 22 * 3600,
        t_stop = 24 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
