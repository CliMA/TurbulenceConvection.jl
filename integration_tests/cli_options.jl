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
        "--moisture_model" # Moisture model (equilibrium or non-equilibrium)
        arg_type = String
        "--precipitation_model" # Precipitation model (None, cutoff or clima_1m)
        arg_type = String
        "--precip_fraction_model" # Precipitation model (prescribed or cloud_cover)
        arg_type = String
        "--prescribed_precip_frac_value" # Value of the precipitation fraction, if prescribed
        arg_type = Float64
        "--precip_fraction_limiter" # Minimum precipitation fraction, if diagnostic
        arg_type = Float64
        "--thermo_covariance_model" # covariance model (prognostic or diagnostic)
        arg_type = String
        "--config"
        help = "Spatial configuration [`sphere` (default), `column`]"
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
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end
