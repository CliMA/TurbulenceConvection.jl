include("cli_options.jl")

if !@isdefined parsed_args
    (s, parsed_args) = parse_commandline()
end

suffix = parsed_args["suffix"]
if !@isdefined case_name
    case_name = parsed_args["case"]
end

import ClimaCore
if parsed_args["trunc_field_type_print"]
    ClimaCore.Fields.truncate_printing_field_types() = true
end

@info "Running $case_name..."
@info "`suffix`: `$(suffix)`"
@info "See `cli_options.jl` changing defaults"

import TurbulenceConvection
using Test

const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "driver", "main.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "overwrite_namelist.jl"))
import .NameList

namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "01$suffix"
overwrite_namelist!(namelist, parsed_args)

integrator, ds_tc_filenames, return_code = main(namelist);

parsed_args["broken_tests"] && exit()

@testset "Simulation completion" begin
    # Test that the simulation has actually finished,
    # and not aborted early.
    @test !(return_code === :simulation_aborted)
    @test return_code === :success
end

# Post-processing case kwargs
include(joinpath(tc_dir, "integration_tests", "sphere_utils.jl"))
include(joinpath(tc_dir, "post_processing", "case_kwargs.jl"))
include(joinpath(tc_dir, "post_processing", "compute_mse.jl"))
include(joinpath(tc_dir, "post_processing", "mse_tables.jl"))
best_mse = get(all_best_mse, namelist["meta"]["casename"], all_best_mse["missing_table_placeholder"]) # default to empty if mse table isn't there ("NA" is what we use for empty to pass test_mse()

parsed_args["skip_post_proc"] && exit()

plot_dir = joinpath(dirname(first(ds_tc_filenames)), "comparison")
if parsed_args["config"] == "sphere"
    plot_profiles(integrator.sol.u[end], plot_dir)
    test_zero_horizontal_variance(integrator.sol.u[end])
end

for ds_tc_filename in ds_tc_filenames
    computed_mse =
        compute_mse_wrapper(case_name, best_mse, ds_tc_filename; case_kwargs[namelist["meta"]["casename"]]..., plot_dir)

    parsed_args["skip_tests"] && exit()

    if parsed_args["config"] == "column"
        open("computed_mse_$case_name.json", "w") do io
            JSON.print(io, computed_mse)
        end
    end

    @testset "$case_name" begin
        for k in keys(best_mse)
            test_mse(computed_mse, best_mse, k)
        end
        @testset "Post-run tests" begin
            isnan_inf_or_filled(x) = isnan(x) || isinf(x) || x ≈ NC.fillvalue(typeof(x))
            NC.Dataset(ds_tc_filename, "r") do ds
                profile = ds.group["profiles"]
                @test !any(isnan_inf_or_filled.(Array(profile["qt_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["updraft_area"])))
                @test !any(isnan_inf_or_filled.(Array(profile["updraft_w"])))
                @test !any(isnan_inf_or_filled.(Array(profile["updraft_qt"])))
                @test !any(isnan_inf_or_filled.(Array(profile["updraft_thetal"])))
                @test !any(isnan_inf_or_filled.(Array(profile["u_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["tke_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["temperature_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["ql_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["qi_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["thetal_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["Hvar_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["QTvar_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["v_mean"])))
                @test !any(isnan_inf_or_filled.(Array(profile["qr_mean"])))
                timeseries = ds.group["timeseries"]
                @test !any(isnan_inf_or_filled.(Array(timeseries["lwp_mean"])))
            end
            nothing
        end
        nothing
    end
end
# Cleanup
# https://github.com/JuliaLang/Pkg.jl/issues/3014
rm("LocalPreferences.toml"; force = true)
