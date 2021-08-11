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

best_mse["updraft_thetal"] = 5.0248694008424222e+00
best_mse["v_mean"] = 4.4630163253493214e+00
best_mse["u_mean"] = 9.6412500983148846e+00
best_mse["tke_mean"] = 2.4683562232142533e+00
best_mse["temperature_mean"] = 8.8678308948818863e-06
best_mse["thetal_mean"] = 8.7900622849842488e-06
best_mse["Hvar_mean"] = 1.2891667586703637e+01
best_mse["QTvar_mean"] = 4.4424755799649029e-01

@testset "GABLS" begin
    case_name = "GABLS"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 7 * 3600,
        t_stop = 9 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
