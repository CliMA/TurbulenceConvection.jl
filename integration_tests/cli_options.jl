import ArgParse

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--job_id"
        help = "A unique job identifier (currently only supported in flame graphs)"
        default = "UNDEFINED"
        "--case"
        help = "Case to run"
        default = "Bomex"
        "--sgs"            # SGS environmental teratment (mean or quadratures)
        arg_type = String
        "--quad_type"      # SGS quadrature type (lognormal or Gaussian)
        arg_type = String
        "--entr"           # Try other entr-detr models
        arg_type = String
        "--entr_dim_scale" # Specify dimensional scale for entrainment
        arg_type = String
        "--detr_dim_scale" # Specify dimensional scale for detrainment
        arg_type = String
        "--area_limiter_power" # Specify area limiter power
        arg_type = Float64
        "--nn_ent_biases" # Specify whether NN parameter vector contains biases
        arg_type = Bool
        "--stoch_entr"     # Choose type of stochastic entr-detr model
        arg_type = String
        "--t_max"          # Simulation time to run to
        arg_type = Float64
        "--adapt_dt"       # use adaptive timestepping
        arg_type = Bool
        "--dt"             # Specify model time step (when not using adaptive dt)
        arg_type = Float64
        "--dt_max"         # Specify maximum model time step (when using adaptive dt)
        arg_type = Float64
        "--calibrate_io"   # Test that calibration IO passes regression tests
        arg_type = Bool
        default = false
        # TODO: Improve this name, it's confusing!!
        "--stretch_grid"   # Test stretched grid option
        arg_type = Bool
        "--skip_io"        # Test that skipping IO passes regression tests
        arg_type = Bool
        default = false
        "--skip_post_proc" # Skip post processing
        arg_type = Bool
        default = false
        "--skip_tests"     # Skip regression tests
        arg_type = Bool
        default = false
        "--broken_tests"     # Tests are broken, skip tests
        arg_type = Bool
        default = false
        "--suffix"         # A suffix for the artifact folder
        arg_type = String
        default = ""
        "--n_up"           # Number of updrafts
        arg_type = Int
        "--moisture_model"
        help = "Cloud condensate formation model [`equilibrium` (default), `nonequilibrium`]"
        arg_type = String
        "--precipitation_model"
        help = "Precipitation model [`None` (default), `cutoff` or `clima_1m`]"
        arg_type = String
        "--rain_formation_scheme"
        help = "Rain autoconversion and accretion scheme [`clima_1m_default` (default) , `KK2000`, `B1994`, `TC1980`, `LD2004`]"
        arg_type = String
        "--prescribed_Nd"
        help = "Prescribed cloud droplet number concentration. Valid when rain_formation_scheme is `KK2000`, `B1994`, `TC1980` or `LD2004`]"
        arg_type = Float64
        "--precip_fraction_model"
        help = "Assumed (constant with height) precipitation fraction [`prescribed` (default), `cloud_cover`]"
        arg_type = String
        "--prescribed_precip_frac_value"
        help = "Value of the precipitation fraction, if prescribed."
        arg_type = Float64
        "--precip_fraction_limiter"
        help = "Minimum precipitation fraction, if diagnostic."
        arg_type = Float64
        "--thermo_covariance_model"
        help = "The type of equation for the sgs covariance [`prognostic`, `diagnostic` (default)]"
        arg_type = String
        "--float_type"
        arg_type = String
        default = "Float64"
        "--config"
        help = "Spatial configuration [`sphere`, `column` (default)]"
        arg_type = String
        default = "column"
        "--set_src_seed"
        help = "Set random seeds for reproducible results per column"
        arg_type = Bool
        default = false
        "--test_duals"
        help = "Test that we can use Duals through ∑tendencies to ForwardDiff through the model"
        arg_type = Bool
        default = false
        "--trunc_field_type_print"
        help = "Set to `true` to truncate printing of ClimaCore `Field` types"
        arg_type = Bool
        default = true
        "--acnv_scaling"
        help = "Coefficient multiplying the autoconversion rate. Default = 1"
        arg_type = Float64
        "--accr_scaling"
        help = "Coefficient multiplying the accretion rate. Default = 1"
        arg_type = Float64
        "--evap_scaling"
        help = "Coefficient multiplying the evaporation rate. Default = 1"
        arg_type = Float64
        "--depsub_scaling"
        help = "Coefficient multiplying the deposition sublimation rate. Default = 1"
        arg_type = Float64
        "--melt_scaling"
        help = "Coefficient multiplying the melting rate. Default = 1"
        arg_type = Float64
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end


"""
    print_repl_script(str::String)

Generate a block of code to run a particular
buildkite job given the `command:` string.

Example:
```julia
include("integration_tests/cli_options.jl")
print_repl_script("...")
```
"""
function print_repl_script(str)
    ib = """"""
    ib *= """\n"""
    ib *= """using Revise; include("integration_tests/cli_options.jl");\n"""
    ib *= """\n"""
    ib *= """(s, parsed_args) = parse_commandline();\n"""
    parsed_args = parsed_args_from_command_line_flags(str)
    for (flag, val) in parsed_args
        if val isa AbstractString
            ib *= "parsed_args[\"$flag\"] = \"$val\";\n"
        else
            ib *= "parsed_args[\"$flag\"] = $val;\n"
        end
    end
    ib *= """\n"""
    ib *= """include("integration_tests/driver.jl")\n"""
    println(ib)
end

function parsed_args_from_command_line_flags(str, parsed_args = Dict())
    s = str
    s = last(split(s, ".jl"))
    s = strip(s)
    parsed_args_list = split(s, " ")
    @assert iseven(length(parsed_args_list))
    parsed_arg_pairs = map(1:2:(length(parsed_args_list) - 1)) do i
        Pair(parsed_args_list[i], strip(parsed_args_list[i + 1], '\"'))
    end
    function parse_arg(val)
        for T in (Bool, Int, Float32, Float64)
            try
                return parse(T, val)
            catch
            end
        end
        return val # string
    end
    for (flag, val) in parsed_arg_pairs
        parsed_args[replace(flag, "--" => "")] = parse_arg(val)
    end
    return parsed_args
end
