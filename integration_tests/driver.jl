import ArgParse

run_from_command_line = abspath(PROGRAM_FILE) == @__FILE__
if run_from_command_line
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--case"
        help = "Case to run"
        default = "Bomex"
        # Command-line args to change default namelist params
        "--micro"          # Try other microphysics quadrature
        "--entr"           # Try other entr-detr models
        "--stoch_entr"     # Choose type of stochastic entr-detr model
        "--calibrate_io"   # Test that calibration IO passes regression tests
        "--skip_io"        # Test that skipping IO passes regression tests
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    case_name = parsed_args["case"]
end

#! format: off
micro        = run_from_command_line ? parsed_args["micro"] : nothing
entr         = run_from_command_line ? parsed_args["entr"] : nothing
stoch_entr   = run_from_command_line ? parsed_args["stoch_entr"] : nothing
calibrate_io = run_from_command_line ? parsed_args["calibrate_io"] : nothing
skip_io      = run_from_command_line ? parsed_args["skip_io"] : nothing

calibrate_io ≠ nothing && (calibrate_io = parse(Bool, calibrate_io))
skip_io ≠ nothing      && (skip_io = parse(Bool, skip_io))
#! format: on

if @isdefined case_name
    @info "Running $case_name..."
else
    case_name = "Bomex"
    @info "Running default case ($case_name)."
    @info "Set `case_name` if you'd like to run a different case"
end

import TurbulenceConvection
using Test

const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
import .NameList

namelist = NameList.default_namelist(case_name)

uuid_suffix(s::String) = "_" * s
uuid_suffix(::Nothing) = ""
suffix = uuid_suffix(micro)
suffix *= uuid_suffix(entr)
suffix *= uuid_suffix(stoch_entr)
calibrate_io ≠ nothing && (suffix *= "_calibrate_io_$calibrate_io")
skip_io ≠ nothing && (suffix *= "_skip_io_$skip_io")

namelist["meta"]["uuid"] = "01$suffix"

#! format: off
micro        ≠ nothing && (namelist["thermodynamics"]["quadrature_type"] = micro)
entr         ≠ nothing && (namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = entr)
stoch_entr   ≠ nothing && (namelist["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] = stoch_entr)
calibrate_io ≠ nothing && (namelist["stats_io"]["calibrate_io"] = calibrate_io)
skip_io      ≠ nothing && (namelist["stats_io"]["skip"] = skip_io)
#! format: on

ds_tc_filename, return_code = main(namelist)

# Post-processing case kwargs
include(joinpath(tc_dir, "post_processing", "case_kwargs.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
best_mse = all_best_mse[case_name]

# TODO: Remove this and compare only the prognostic state
#       we may need a more generic / flexible version of
#       compute_mse_wrapper.
calibrate_io ≠ nothing && exit()
skip_io ≠ nothing && exit()

computed_mse = compute_mse_wrapper(case_name, best_mse, ds_tc_filename; case_kwargs[case_name]...)

# We'll need to add mse tables if we want
# to use regression tests on experiments
# that deviate from defaults.
micro ≠ nothing && exit()
entr ≠ nothing && exit()
stoch_entr ≠ nothing && exit()

open("computed_mse_$case_name.json", "w") do io
    JSON.print(io, computed_mse)
end

@testset "$case_name" begin
    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    @testset "Post-run tests" begin
        isnan_or_inf(x) = isnan(x) || isinf(x)
        NC.Dataset(ds_tc_filename, "r") do ds
            profile = ds.group["profiles"]
            @test !any(isnan_or_inf.(Array(profile["qt_mean"])))
            @test !any(isnan_or_inf.(Array(profile["updraft_area"])))
            @test !any(isnan_or_inf.(Array(profile["updraft_w"])))
            @test !any(isnan_or_inf.(Array(profile["updraft_qt"])))
            @test !any(isnan_or_inf.(Array(profile["updraft_thetal"])))
            @test !any(isnan_or_inf.(Array(profile["u_mean"])))
            @test !any(isnan_or_inf.(Array(profile["tke_mean"])))
            @test !any(isnan_or_inf.(Array(profile["temperature_mean"])))
            @test !any(isnan_or_inf.(Array(profile["ql_mean"])))
            @test !any(isnan_or_inf.(Array(profile["qi_mean"])))
            @test !any(isnan_or_inf.(Array(profile["thetal_mean"])))
            @test !any(isnan_or_inf.(Array(profile["Hvar_mean"])))
            @test !any(isnan_or_inf.(Array(profile["QTvar_mean"])))
            @test !any(isnan_or_inf.(Array(profile["v_mean"])))
            @test !any(isnan_or_inf.(Array(profile["qr_mean"])))
        end
        nothing
    end

    @testset "Simulation completion" begin
        # Test that the simulation has actually finished,
        # and not aborted early.
        @test !(return_code == :simulation_aborted)
        @test return_code == :success
    end
    nothing
end
