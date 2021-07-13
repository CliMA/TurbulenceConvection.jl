import TurbulenceConvection
using TurbulenceConvection
using Test
using Random

# Make deterministic:
Random.seed!(1234)

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

best_mse = OrderedDict()
best_mse["qt_mean"] = 6.7598618147035705e-01
best_mse["updraft_fraction"] = 1.1151497561239513e+02
best_mse["updraft_w"] = 6.7164800402645284e+01
best_mse["updraft_qt"] = 1.4426382216507470e+01
best_mse["updraft_thetali"] = 1.0040573552118677e+02
best_mse["v_mean"] = 1.0040573552118677e+02
best_mse["u_mean"] = 1.0040573552118677e+02
best_mse["rho"] = 1.0040573552118677e+02
best_mse["tke_mean"] = 9.4610544693497943e+01
best_mse["env_qt2"] = 9.4610544693497943e+01

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r")
get_ds_filename(case) = "test/Output.{$case}.{01}/statsStats.{$case}.nc";

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = NameList.Bomex(default_namelist("Bomex"))
    paramlist = ParamList.Bomex(default_paramlist("Bomex"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)
    @show ds_filename

    # computed_mse = Dataset(ds_filename, "r") do ds
    #     compute_mse(
    #         ds,
    #         ds_pycles,
    #         "Bomex",
    #         best_mse,
    #         dirname(ds_filename);
    #         plot_comparison=true
    #     )
    # end

    # test_mse(computed_mse, best_mse, :q_tot_gm)
    # test_mse(computed_mse, best_mse, :a_up)
    # test_mse(computed_mse, best_mse, :w_up)
    # test_mse(computed_mse, best_mse, :q_tot_up)
    # test_mse(computed_mse, best_mse, :θ_liq_up)
    # test_mse(computed_mse, best_mse, :tke_en)
end

