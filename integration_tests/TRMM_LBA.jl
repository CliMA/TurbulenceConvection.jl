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
best_mse["qt_mean"] = 1.6060537595274345e+00
best_mse["updraft_area"] = 2.4326719919434730e+04
best_mse["updraft_w"] = 8.9781354744734904e+02
best_mse["updraft_qt"] = 2.7840588043213661e+01
best_mse["updraft_thetal"] = 1.1000857010242179e+02
best_mse["v_mean"] = 2.9255417406952165e+02
best_mse["u_mean"] = 1.6872488577305819e+03
best_mse["tke_mean"] = 1.5447663031746142e+03
best_mse["temperature_mean"] = 6.8027838562588429e-04
best_mse["ql_mean"] = 9.5570328764458566e+02
best_mse["thetal_mean"] = 8.1813107572736785e-03
best_mse["Hvar_mean"] = 6.7055563585233622e+03
best_mse["QTvar_mean"] = 2.6582568278836857e+03

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
