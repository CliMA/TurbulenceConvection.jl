if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList
using .ParamList

best_mse = OrderedDict()
best_mse["qt_mean"] = 3.8001673126194867e-02
best_mse["ql_mean"] = 8.4578478184878048e+00
best_mse["updraft_area"] = 2.2324280164372593e+02
best_mse["updraft_w"] = 3.4444820774461329e+00
best_mse["updraft_qt"] = 1.3810722649106608e+00
best_mse["updraft_thetal"] = 1.2761175148454377e+01
best_mse["v_mean"] = 4.0030908804636987e+01
best_mse["u_mean"] = 3.5747541514399444e+01
best_mse["tke_mean"] = 1.4604768548409767e+01


@testset "DYCOMS_RF01" begin
    println("Running DYCOMS_RF01...")
    namelist = NameList.DYCOMS_RF01(default_namelist("DYCOMS_RF01"))
    paramlist = ParamList.DYCOMS_RF01(default_paramlist("DYCOMS_RF01"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "DYCOMS_RF01.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "DYCOMS_RF01.nc"), "r") do ds_scampy
                compute_mse(
                    "DYCOMS_RF01",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv=ds,
                    ds_scampy=ds_scampy,
                    ds_pycles=ds_pycles,
                    plot_comparison=true
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end

