include("cli_options.jl")

(s, parsed_args) = parse_commandline()

suffix = parsed_args["suffix"]
if !@isdefined case_name
    case_name = parsed_args["case"]
end

@info "Running $case_name..."
@info "`suffix`: `$(suffix)`"
@info "See `cli_options.jl` changing defaults"

import TurbulenceConvection
using Test

const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
import .NameList

namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "01$suffix"

#! format: off
!isnothing(parsed_args["sgs"]) && (namelist["thermodynamics"]["sgs"] = parsed_args["sgs"])
!isnothing(parsed_args["quad_type"]) && (namelist["thermodynamics"]["quadrature_type"] = parsed_args["quad_type"])
!isnothing(parsed_args["entr"]) && (namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = parsed_args["entr"])
!isnothing(parsed_args["stoch_entr"]) && (namelist["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] = parsed_args["stoch_entr"])
!isnothing(parsed_args["t_max"]) && (namelist["time_stepping"]["t_max"] = parsed_args["t_max"])
!isnothing(parsed_args["adapt_dt"]) && (namelist["time_stepping"]["adapt_dt"] = parsed_args["adapt_dt"])
!isnothing(parsed_args["dt"]) && (namelist["time_stepping"]["dt_min"] = parsed_args["dt"])
!isnothing(parsed_args["calibrate_io"]) && (namelist["stats_io"]["calibrate_io"] = parsed_args["calibrate_io"])
!isnothing(parsed_args["stretch_grid"]) && (namelist["grid"]["stretch"]["flag"] = parsed_args["stretch_grid"])
!isnothing(parsed_args["skip_io"]) && (namelist["stats_io"]["skip"] = parsed_args["skip_io"])
!isnothing(parsed_args["n_up"]) && (namelist["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = parsed_args["n_up"])
#! format: on

ds_tc_filename, return_code = main(namelist)

# Post-processing case kwargs
include(joinpath(tc_dir, "post_processing", "case_kwargs.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
best_mse = all_best_mse[case_name]

parsed_args["skip_post_proc"] && exit()

computed_mse = compute_mse_wrapper(case_name, best_mse, ds_tc_filename; case_kwargs[case_name]...)

parsed_args["skip_tests"] && exit()

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

# Cleanup
# https://github.com/JuliaLang/Pkg.jl/issues/3014
rm("LocalPreferences.toml"; force = true)
