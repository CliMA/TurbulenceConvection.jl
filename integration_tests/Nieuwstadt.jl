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
best_mse["updraft_area"] = 5.9567286602931904e+02
best_mse["updraft_w"] = 2.6450205443296568e+01
best_mse["updraft_thetal"] = 3.0475209174359087e+01
best_mse["u_mean"] = 1.5244498152508007e+02
best_mse["tke_mean"] = 7.3585026092564277e+01
best_mse["temperature_mean"] = 1.1872218655217143e-05
best_mse["thetal_mean"] = 1.2035241147904184e-05
best_mse["Hvar_mean"] = 1.8640506843913366e+02

@testset "Nieuwstadt" begin
    case_name = "Nieuwstadt"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)

    computed_mse = compute_mse_wrapper(
        case_name,
        best_mse,
        ds_tc_filename;
        plot_comparison = true,
        t_start = 6 * 3600,
        t_stop = 8 * 3600,
    )

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
