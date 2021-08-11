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
best_mse["updraft_area"] = 8.8107267096730493e+02
best_mse["updraft_w"] = 1.8973498180481789e+02
best_mse["updraft_thetal"] = 5.2905285445151222e-05
best_mse["u_mean"] = 2.4701556169535316e-27
best_mse["tke_mean"] = 4.0953307784558532e+00
best_mse["temperature_mean"] = 3.2275977622863615e-05
best_mse["thetal_mean"] = 2.5687103839396313e-05
best_mse["Hvar_mean"] = 3.2146754706000658e+00

@testset "DryBubble" begin
    case_name = "DryBubble"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse =
        compute_mse_wrapper(case_name, best_mse, ds_tc_filename; plot_comparison = true, t_start = 900, t_stop = 1000)

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
