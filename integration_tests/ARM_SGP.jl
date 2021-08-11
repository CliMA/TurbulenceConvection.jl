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
best_mse["qt_mean"] = 3.5691847012203143e-01
best_mse["updraft_area"] = 1.9645743867371027e+03
best_mse["updraft_w"] = 3.5781735855735855e+02
best_mse["updraft_qt"] = 1.3457137434986310e+01
best_mse["updraft_thetal"] = 2.7680764603306077e+01
best_mse["u_mean"] = 8.7998547277817920e+01
best_mse["tke_mean"] = 6.3658705239459096e+02
best_mse["temperature_mean"] = 1.3818806435138376e-04
best_mse["ql_mean"] = 3.6220527567539864e+02
best_mse["thetal_mean"] = 1.4145191959567767e-04
best_mse["Hvar_mean"] = 1.5980068858564114e+03
best_mse["QTvar_mean"] = 3.6076370680937271e+02

@testset "ARM_SGP" begin
    println("Running ARM_SGP...")
    namelist = default_namelist("ARM_SGP")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "ARM_SGP.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "ARM_SGP.nc"), "r") do ds_scampy
                compute_mse(
                    "ARM_SGP",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 8 * 3600,
                    t_stop = 11 * 3600,
                )
            end
        end
    end

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end
