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
best_mse["qt_mean"] = 3.2231237240917160e-01
best_mse["updraft_area"] = 2.0078340302396784e+03
best_mse["updraft_w"] = 3.3936433966850967e+02
best_mse["updraft_qt"] = 1.3509451158258978e+01
best_mse["updraft_thetal"] = 2.7683871839666477e+01
best_mse["u_mean"] = 8.7998547277817920e+01
best_mse["tke_mean"] = 6.9016396604551733e+02
best_mse["temperature_mean"] = 1.3318895186222515e-04
best_mse["ql_mean"] = 3.5021120774412930e+02
best_mse["thetal_mean"] = 1.3800056762098192e-04
best_mse["Hvar_mean"] = 6.2976082776935509e+03
best_mse["QTvar_mean"] = 4.3109900526010560e+03

@testset "ARM_SGP" begin
    case_name = "ARM_SGP"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 8 * 3600,
        t_stop = 11 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
