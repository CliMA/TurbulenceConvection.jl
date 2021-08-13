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
best_mse["qt_mean"] = 1.6511493474010635e-02
best_mse["ql_mean"] = 5.2388152471828944e+00
best_mse["updraft_area"] = 2.3937655332711455e+02
best_mse["updraft_w"] = 4.2950818025070792e+00
best_mse["updraft_qt"] = 1.1670622065147132e+00
best_mse["updraft_thetal"] = 1.2740701334370417e+01
best_mse["v_mean"] = 3.9746921720554738e+01
best_mse["u_mean"] = 3.7046560343557694e+01
best_mse["tke_mean"] = 1.4700070268032464e+01
best_mse["temperature_mean"] = 2.1532443068564967e-05
best_mse["thetal_mean"] = 2.2397858587086213e-05
best_mse["Hvar_mean"] = 8.2677316057727712e+03
best_mse["QTvar_mean"] = 6.0266525106233109e+02

@testset "DYCOMS_RF01" begin
    case_name = "DYCOMS_RF01"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 2 * 3600,
        t_stop = 4 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
