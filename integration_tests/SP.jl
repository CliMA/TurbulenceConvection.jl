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
best_mse["qt_mean"] = 3.5246740566216994e+00
best_mse["updraft_area"] = 3.4922216987359178e+00
best_mse["updraft_w"] = 1.0216753899526712e+00
best_mse["updraft_qt"] = 1.4126704196501478e+00
best_mse["updraft_thetal"] = 1.0513205076428579e-01
best_mse["v_mean"] = 4.6470417414900161e-01
best_mse["u_mean"] = 7.3604590399213739e-05
best_mse["tke_mean"] = 4.7848458015049938e-01
best_mse["temperature_mean"] = 6.9541909880873242e-07
best_mse["thetal_mean"] = 5.2370121754610367e-07
best_mse["Hvar_mean"] = 2.3999517246677664e+02
best_mse["QTvar_mean"] = 1.6280851822250700e+01

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
