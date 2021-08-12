if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
using Test
import TurbulenceConvection
const TC = TurbulenceConvection
const tc_dir = dirname(dirname(pathof(TC)))

include(joinpath(tc_dir, "integration_tests", "utils", "main.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "compute_mse.jl"))
using .NameList

@testset "Unit tests" begin
    # TODO: add unit tests.
end

@testset "Bomex" begin
    best_mse = OrderedDict()
    best_mse["qt_mean"] = 9.8858939409301599e-02
    best_mse["updraft_area"] = 7.2423270403289905e+02
    best_mse["updraft_w"] = 2.7297489157558520e+01
    best_mse["updraft_qt"] = 4.0686060603282321e+00
    best_mse["updraft_thetal"] = 2.1547376374778047e+01
    best_mse["v_mean"] = 6.5159453481189033e+01
    best_mse["u_mean"] = 5.3292347667845092e+01
    best_mse["tke_mean"] = 3.8707689709123862e+01
    best_mse["temperature_mean"] = 4.0396510560887023e-05
    best_mse["ql_mean"] = 3.9678285547654855e+00
    best_mse["thetal_mean"] = 4.1130421765389724e-05
    best_mse["Hvar_mean"] = 5.9914307236483694e+01
    best_mse["QTvar_mean"] = 1.9732739289298404e+01

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename = @time main(namelist)
    nothing
end
