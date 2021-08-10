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
best_mse["qt_mean"] = 7.2336971753194734e-04
best_mse["ql_mean"] = 3.0610499166624661e-01
best_mse["updraft_area"] = 1.0360077651771860e-01
best_mse["updraft_w"] = 3.1237334051283400e-02
best_mse["updraft_qt"] = 4.8044705803477854e-03
best_mse["updraft_thetal"] = 2.7425547386093370e-06
best_mse["v_mean"] = 4.0731966024459484e-02
best_mse["u_mean"] = 2.7156874610105889e-04
best_mse["tke_mean"] = 4.2551526373553514e-02
best_mse["temperature_mean"] = 4.7751904212962137e-07
best_mse["thetal_mean"] = 3.2930923191757882e-07
best_mse["Hvar_mean"] = 5.4910770906253944e+01
best_mse["QTvar_mean"] = 1.7379261184205511e+01

@testset "life_cycle_Tan2018" begin
    println("Running life_cycle_Tan2018...")
    namelist = default_namelist("life_cycle_Tan2018")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(SCAMPy_output_dataset_path, "life_cycle_Tan2018.nc"), "r") do ds_scampy
            compute_mse(
                "life_cycle_Tan2018",
                best_mse,
                joinpath(dirname(ds_filename), "comparison");
                ds_turb_conv = ds,
                ds_scampy = ds_scampy,
                plot_comparison = true,
                t_start = 4 * 3600,
                t_stop = 6 * 3600,
            )
        end
    end

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
