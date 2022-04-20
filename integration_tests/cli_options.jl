import ArgParse

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--case"
        help = "Case to run"
        default = "Bomex"
        "--micro"          # Try other microphysics quadrature
        default = "log-normal"
        "--entr"           # Try other entr-detr models
        arg_type = String
        default = "moisture_deficit"
        "--stoch_entr"     # Choose type of stochastic entr-detr model
        arg_type = String
        "--t_max"          # Simulation time to run to
        arg_type = Float64
        "--calibrate_io"   # Test that calibration IO passes regression tests
        arg_type = Bool
        default = false
        "--stretch_grid"   # Test stretched grid option
        arg_type = Bool
        default = false
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
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end
