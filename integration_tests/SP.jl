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
best_mse["qt_mean"] = 3.5073036120908645e+00
best_mse["updraft_area"] = 3.9071035000833305e+00
best_mse["updraft_w"] = 9.3648209543530914e-01
best_mse["updraft_qt"] = 1.3868858638582151e+00
best_mse["updraft_thetal"] = 1.0515272505644875e-01
best_mse["v_mean"] = 4.6000262198139841e-01
best_mse["u_mean"] = 7.3724093159641325e-05
best_mse["tke_mean"] = 4.7833665690724869e-01
best_mse["temperature_mean"] = 6.8550010655096847e-07
best_mse["thetal_mean"] = 5.1299556224606434e-07
best_mse["Hvar_mean"] = 3.1719858004705035e+01
best_mse["QTvar_mean"] = 3.9762684967056670e+00

@testset "SP" begin
    case_name = "SP"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse =
        compute_mse_wrapper(case_name, best_mse, ds_tc_filename; plot_comparison = true, t_start = 0, t_stop = 2 * 3600)

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
