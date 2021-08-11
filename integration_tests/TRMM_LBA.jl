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

# Note: temperatures in this case become extremely low.
CLIMAParameters.Planet.T_freeze(::EarthParameterSet) = 100.0

best_mse = OrderedDict()
best_mse["qt_mean"] = 2.1180570149373197e+00
best_mse["updraft_area"] = 2.2911123097125790e+04
best_mse["updraft_w"] = 9.9122631422628353e+02
best_mse["updraft_qt"] = 3.0750107437154107e+01
best_mse["updraft_thetal"] = 1.1001770046016014e+02
best_mse["v_mean"] = 2.9250578870751696e+02
best_mse["u_mean"] = 1.6873177041648653e+03
best_mse["tke_mean"] = 9.3810901861977175e+02
best_mse["temperature_mean"] = 8.1897112533893871e-04
best_mse["ql_mean"] = 7.3150520783875129e+02
best_mse["thetal_mean"] = 8.2746696205323478e-03
best_mse["Hvar_mean"] = 3.5185010182273427e+03
best_mse["QTvar_mean"] = 1.7745546315637387e+03

@testset "TRMM_LBA" begin
    case_name = "TRMM_LBA"
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
