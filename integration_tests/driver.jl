include("cli_options.jl")

(s, parsed_args) = parse_commandline()

suffix = parsed_args["suffix"]
if !@isdefined case_name
    case_name = parsed_args["case"]
end

parsed_args["trunc_stack_traces"] && include("truncate_stack_traces.jl")

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
"sgs"          => (nl, pa, key) -> (nl["thermodynamics"]["sgs"] = pa[key]),
"quad_type"    => (nl, pa, key) -> (nl["thermodynamics"]["quadrature_type"] = pa[key]),
"entr"         => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = pa[key]),
"stoch_entr"   => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] = pa[key]),
"t_max"        => (nl, pa, key) -> (nl["time_stepping"]["t_max"] = pa[key]),
"adapt_dt"     => (nl, pa, key) -> (nl["time_stepping"]["adapt_dt"] = pa[key]),
"dt"           => (nl, pa, key) -> (nl["time_stepping"]["dt_min"] = pa[key]),
"calibrate_io" => (nl, pa, key) -> (nl["stats_io"]["calibrate_io"] = pa[key]),
"stretch_grid" => (nl, pa, key) -> (nl["grid"]["stretch"]["flag"] = pa[key]),
"skip_io"      => (nl, pa, key) -> (nl["stats_io"]["skip"] = pa[key]),
"n_up"         => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = pa[key]),
)
#! format: on
for key in keys(overwrite_namelist_map)
    if !isnothing(parsed_args[key])
        @warn "Parameter `$key` overwriting namelist"
        overwrite_namelist_map[key](namelist, parsed_args, key)
    end
end

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
