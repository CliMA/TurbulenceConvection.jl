include("cli_options.jl")

(s, parsed_args) = parse_commandline()

suffix = parsed_args["suffix"]
if !@isdefined case_name
    case_name = parsed_args["case"]
end

parsed_args["trunc_field_type_print"] && include("trunc_field_type_print.jl")

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
overwrite_namelist_map = Dict(
"sgs"                     => (nl, pa, key) -> (nl["thermodynamics"]["sgs"] = pa[key]),
"quad_type"               => (nl, pa, key) -> (nl["thermodynamics"]["quadrature_type"] = pa[key]),
"entr"                    => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = pa[key]),
"stoch_entr"              => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] = pa[key]),
"t_max"                   => (nl, pa, key) -> (nl["time_stepping"]["t_max"] = pa[key]),
"adapt_dt"                => (nl, pa, key) -> (nl["time_stepping"]["adapt_dt"] = pa[key]),
"dt"                      => (nl, pa, key) -> (nl["time_stepping"]["dt_min"] = pa[key]),
"calibrate_io"            => (nl, pa, key) -> (nl["stats_io"]["calibrate_io"] = pa[key]),
"stretch_grid"            => (nl, pa, key) -> (nl["grid"]["stretch"]["flag"] = pa[key]),
"skip_io"                 => (nl, pa, key) -> (nl["stats_io"]["skip"] = pa[key]),
"n_up"                    => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = pa[key]),
"moisture_model"          => (nl, pa, key) -> (nl["thermodynamics"]["moisture_model"] = pa[key]),
"precipitation_model"     => (nl, pa, key) -> (nl["microphysics"]["precipitation_model"] = pa[key]),
"precip_fraction_model"   => (nl, pa, key) -> (nl["microphysics"]["precip_fraction_model"] = pa[key]),
"prescribed_precip_frac_value" => (nl, pa, key) -> (nl["microphysics"]["prescribed_precip_frac_value"] = pa[key]),
"precip_fraction_limiter" => (nl, pa, key) -> (nl["microphysics"]["precip_fraction_limiter"] = pa[key]),
"thermo_covariance_model" => (nl, pa, key) -> (nl["thermodynamics"]["thermo_covariance_model"] = pa[key]),
"energy_var"              => (nl, pa, key) -> (nl["energy_var"] = pa[key]),
)
no_overwrites = (
    "case", # default_namelist already overwrites namelist["meta"]["casename"]
    "skip_post_proc",
    "skip_tests",
    "broken_tests",
    "trunc_field_type_print",
    "suffix",
)
#! format: on
for key in keys(overwrite_namelist_map)
    if !isnothing(parsed_args[key])
        @warn "Parameter `$key` overwriting namelist"
        overwrite_namelist_map[key](namelist, parsed_args, key)
    end
end

# Error check overwrites:
# A tuple of strings for all CL arguments that do _not_
overwrite_list = map(collect(keys(overwrite_namelist_map))) do key
    (key, !isnothing(parsed_args[key]))
end
filter!(x -> x[2], overwrite_list)
cl_list = map(collect(keys(parsed_args))) do key
    (key, !isnothing(parsed_args[key]) && !(key in no_overwrites))
end
filter!(x -> x[2], cl_list)
if length(overwrite_list) ≠ length(cl_list)
    error(
        string(
            "A prescribed CL argument is not overwriting the namelist.",
            "It seems that a CL argument was added, and the `no_overwrites`",
            "or `overwrite_namelist_map` must be updated.",
        ),
    )
end

ds_tc_filename, return_code = main(namelist)

parsed_args["broken_tests"] && exit()

@testset "Simulation completion" begin
    # Test that the simulation has actually finished,
    # and not aborted early.
    @test !(return_code == :simulation_aborted)
    @test return_code == :success
end

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

# Cleanup
# https://github.com/JuliaLang/Pkg.jl/issues/3014
rm("LocalPreferences.toml"; force = true)
