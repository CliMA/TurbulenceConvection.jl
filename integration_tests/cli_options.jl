import ArgParse

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--case"
        help = "Case to run"
        default = "Bomex"
        "--sgs"            # SGS environmental teratment (mean or quadratures)
        arg_type = String
        "--quad_type"      # SGS quadrature type (lognormal or Gaussian)
        arg_type = String
        "--entr"           # Try other entr-detr models
        arg_type = String
        "--stoch_entr"     # Choose type of stochastic entr-detr model
        arg_type = String
        "--output_root"
        arg_type = String
        default = "./"
        "--t_max"          # Simulation time to run to
        arg_type = Float64
        "--adapt_dt"       # use adaptive timestepping
        arg_type = Bool
        "--dt"             # Specify model time step (when not using adaptive dt)
        arg_type = Float64
        "--calibrate_io"   # Test that calibration IO passes regression tests
        arg_type = Bool
        default = false
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
        "--suffix"         # A suffix for the artifact folder
        arg_type = String
        default = ""
        "--n_up"           # Number of updrafts
        arg_type = Int
        "--moisture_model" # Moisture model (equilibrium or non-equilibrium)
        arg_type = String
        "--trunc_stack_traces"
        help = "Set to `true` to truncate printing of ClimaCore `Field`s"
        arg_type = Bool
        default = true
        "--tau_cond_evap" #Until we figure out how to do it from main file
        arg_type = Float64
        "--tau_sub_dep" #Until we figure out how to do it from main file
        arg_type = Float64
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end
