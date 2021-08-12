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
best_mse["qt_mean"] = 1.0525043577985331e-01
best_mse["updraft_area"] = 6.9819789786346007e+02
best_mse["updraft_w"] = 3.2787315815313967e+01
best_mse["updraft_qt"] = 4.2909342180533088e+00
best_mse["updraft_thetal"] = 2.1546675837401239e+01
best_mse["v_mean"] = 6.8303124684318490e+01
best_mse["u_mean"] = 5.3308086120368444e+01
best_mse["tke_mean"] = 4.2971736553362518e+01
best_mse["temperature_mean"] = 5.0980133005003681e-05
best_mse["ql_mean"] = 6.6199873853720544e+02
best_mse["thetal_mean"] = 4.6078132934362555e-05
best_mse["Hvar_mean"] = 3.7127649327550598e+03
best_mse["QTvar_mean"] = 2.0267697132357714e+03

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
