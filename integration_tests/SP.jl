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
best_mse["qt_mean"] = 3.5073036121827599e+00
best_mse["updraft_area"] = 3.9071034998892715e+00
best_mse["updraft_w"] = 9.3648209541424188e-01
best_mse["updraft_qt"] = 1.3868858637940826e+00
best_mse["updraft_thetal"] = 1.0515272505644156e-01
best_mse["v_mean"] = 4.6000262200176228e-01
best_mse["u_mean"] = 7.3724093159961937e-05
best_mse["tke_mean"] = 4.7833665685836724e-01
best_mse["temperature_mean"] = 6.8550010657332516e-07
best_mse["thetal_mean"] = 5.1299556226337377e-07
best_mse["Hvar_mean"] = 3.1719859098824500e+01
best_mse["QTvar_mean"] = 3.9762684439302052e+00

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
