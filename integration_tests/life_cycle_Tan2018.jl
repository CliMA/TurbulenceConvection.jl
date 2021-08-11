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
best_mse["qt_mean"] = 1.1173360639315133e-03
best_mse["ql_mean"] = 3.4722265591098894e-01
best_mse["updraft_area"] = 7.5356102070111355e-02
best_mse["updraft_w"] = 6.0687216986053602e-02
best_mse["updraft_qt"] = 1.3804511006887350e-02
best_mse["updraft_thetal"] = 7.7957427519555628e-06
best_mse["v_mean"] = 3.3600765939282476e-02
best_mse["u_mean"] = 2.8061427949272865e-04
best_mse["tke_mean"] = 5.5713881488773903e-02
best_mse["temperature_mean"] = 6.5814724055873443e-07
best_mse["thetal_mean"] = 5.1637429951845128e-07
best_mse["Hvar_mean"] = 2.6446943448553096e+01
best_mse["QTvar_mean"] = 9.2756222550545093e+00

@testset "life_cycle_Tan2018" begin
    case_name = "life_cycle_Tan2018"
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
