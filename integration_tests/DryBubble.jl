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
best_mse["updraft_area"] = 6.8552893703976156e+02
best_mse["updraft_w"] = 1.6342412689086376e+02
best_mse["updraft_thetal"] = 3.9780037295736014e-05
best_mse["u_mean"] = 1.9502448099351233e-27
best_mse["tke_mean"] = 1.9987696066076023e+05
best_mse["temperature_mean"] = 3.2539779149902821e-05
best_mse["thetal_mean"] = 2.5848228458179228e-05
best_mse["Hvar_mean"] = 7.3771757968047757e+02

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
