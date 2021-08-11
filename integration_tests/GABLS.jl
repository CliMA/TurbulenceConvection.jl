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

best_mse["updraft_thetal"] = 5.0248696023347037e+00
best_mse["v_mean"] = 4.4593534457868529e+00
best_mse["u_mean"] = 9.6414943665200035e+00
best_mse["tke_mean"] = 2.4674095133951375e+00
best_mse["temperature_mean"] = 8.8584843672667532e-06
best_mse["thetal_mean"] = 8.7856734759460943e-06
best_mse["Hvar_mean"] = 1.2892749042279126e+01
best_mse["QTvar_mean"] = 4.4456710317999498e-01

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
