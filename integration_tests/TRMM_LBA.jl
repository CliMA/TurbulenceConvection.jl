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
best_mse["qt_mean"] = 2.1123818327762196e+00
best_mse["updraft_area"] = 2.2936100671363765e+04
best_mse["updraft_w"] = 9.7871834647382809e+02
best_mse["updraft_qt"] = 3.0602076184487512e+01
best_mse["updraft_thetal"] = 1.1001487658208509e+02
best_mse["v_mean"] = 2.9250821023986430e+02
best_mse["u_mean"] = 1.6873198920581733e+03
best_mse["tke_mean"] = 9.3648065006898219e+02
best_mse["temperature_mean"] = 8.1826613372337135e-04
best_mse["ql_mean"] = 7.2432912110535301e+02
best_mse["thetal_mean"] = 8.2748526945473598e-03
best_mse["Hvar_mean"] = 3.5701761092936244e+03
best_mse["QTvar_mean"] = 1.7865858561907439e+03

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
