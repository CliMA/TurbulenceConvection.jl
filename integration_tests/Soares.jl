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
best_mse["qt_mean"] = 1.4966760733610701e-01
best_mse["updraft_area"] = 4.4752734228297282e+02
best_mse["updraft_w"] = 2.1338748581172648e+01
best_mse["updraft_qt"] = 1.1050242532811241e+01
best_mse["updraft_thetal"] = 2.2394553996747693e+01
best_mse["u_mean"] = 7.3068371493591656e+02
best_mse["tke_mean"] = 5.9234734475094214e+01
best_mse["temperature_mean"] = 1.1583915526092593e-05
best_mse["thetal_mean"] = 1.2172922371309679e-05
best_mse["Hvar_mean"] = 2.2601616515636465e+02

@testset "Soares" begin
    case_name = "Soares"
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
