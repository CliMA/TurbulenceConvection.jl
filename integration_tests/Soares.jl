if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

best_mse = OrderedDict()
best_mse["qt_mean"] = 2.5379707482392216e-01
best_mse["updraft_area"] = 6.9134632696535311e+02
best_mse["updraft_w"] = 2.2088311400548207e+01
best_mse["updraft_qt"] = 1.0737238763474318e+01
best_mse["updraft_thetal"] = 2.2118711589015625e+01
best_mse["u_mean"] = 4.1797211980595075e+03
best_mse["tke_mean"] = 7.9767693285775977e+01

@testset "Soares" begin
    println("Running Soares...")
    namelist = NameList.Soares(default_namelist("Soares"))
    paramlist = ParamList.Soares(default_paramlist("Soares"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Soares.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Soares.nc"), "r") do ds_scampy
                compute_mse(
                    "Soares",
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

