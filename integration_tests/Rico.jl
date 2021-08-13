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
best_mse["qt_mean"] = 3.6183738581707564e-01
best_mse["updraft_area"] = 1.9158389580482219e+03
best_mse["updraft_w"] = 1.7058718197542106e+02
best_mse["updraft_qt"] = 1.5449827839901584e+01
best_mse["updraft_thetal"] = 6.3602297256853468e+01
best_mse["v_mean"] = 1.0630514621668171e+02
best_mse["u_mean"] = 1.1443613156137732e+02
best_mse["tke_mean"] = 3.1910880264617197e+02
best_mse["temperature_mean"] = 1.8245777363727451e-04
best_mse["ql_mean"] = 2.2808514972615623e+02
best_mse["thetal_mean"] = 1.5458919462734390e-04
best_mse["Hvar_mean"] = 1.0744099006973258e+04
best_mse["QTvar_mean"] = 4.7454844707739849e+03

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
