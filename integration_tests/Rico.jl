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
best_mse["qt_mean"] = 4.3104248538733420e-01
best_mse["updraft_area"] = 1.5789778929006447e+03
best_mse["updraft_w"] = 2.9474627527892022e+02
best_mse["updraft_qt"] = 1.7944807685758455e+01
best_mse["updraft_thetal"] = 6.3569098888788886e+01
best_mse["v_mean"] = 1.0627585162881577e+02
best_mse["u_mean"] = 1.1431212367609878e+02
best_mse["tke_mean"] = 8.8307805728796393e+02
best_mse["temperature_mean"] = 1.2372854790132712e-04
best_mse["ql_mean"] = 6.9240573720579249e+01
best_mse["thetal_mean"] = 1.2362373538586256e-04
best_mse["Hvar_mean"] = 1.5559806585826123e+03
best_mse["QTvar_mean"] = 7.3615953416203286e+02

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
