function overwrite_namelist!(namelist, parsed_args)
    #! format: off
    overwrite_namelist_map = Dict(
    "sgs"                          => (nl, pa, key) -> (nl["thermodynamics"]["sgs"] = pa[key]),
    "quad_type"                    => (nl, pa, key) -> (nl["thermodynamics"]["quadrature_type"] = pa[key]),
    "entr"                         => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = pa[key]),
    "ml_entr"                      => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["ml_entrainment"] = pa[key]),
    "entr_dim_scale"               => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] = pa[key]),
    "detr_dim_scale"               => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["detr_dim_scale"] = pa[key]),
    "nn_ent_biases"                => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["nn_ent_biases"] = pa[key]),
    "stoch_entr"                   => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] = pa[key]),
    "t_max"                        => (nl, pa, key) -> (nl["time_stepping"]["t_max"] = pa[key]),
    "adapt_dt"                     => (nl, pa, key) -> (nl["time_stepping"]["adapt_dt"] = pa[key]),
    "dt"                           => (nl, pa, key) -> (nl["time_stepping"]["dt_min"] = pa[key]),
    "dt_max"                       => (nl, pa, key) -> (nl["time_stepping"]["dt_max"] = pa[key]),
    "calibrate_io"                 => (nl, pa, key) -> (nl["stats_io"]["calibrate_io"] = pa[key]),
    "stretch_grid"                 => (nl, pa, key) -> (nl["grid"]["stretch"]["flag"] = pa[key]),
    "skip_io"                      => (nl, pa, key) -> (nl["stats_io"]["skip"] = pa[key]),
    "n_up"                         => (nl, pa, key) -> (nl["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = pa[key]),
    "moisture_model"               => (nl, pa, key) -> (nl["thermodynamics"]["moisture_model"] = pa[key]),
    "precipitation_model"          => (nl, pa, key) -> (nl["microphysics"]["precipitation_model"] = pa[key]),
    "rain_formation_scheme"        => (nl, pa, key) -> (nl["microphysics"]["rain_formation_scheme"] = pa[key]),
    "prescribed_Nd"                => (nl, pa, key) -> (nl["microphysics"]["prescribed_Nd"] = pa[key]),
    "precip_fraction_model"        => (nl, pa, key) -> (nl["microphysics"]["precip_fraction_model"] = pa[key]),
    "prescribed_precip_frac_value" => (nl, pa, key) -> (nl["microphysics"]["prescribed_precip_frac_value"] = pa[key]),
    "precip_fraction_limiter"      => (nl, pa, key) -> (nl["microphysics"]["precip_fraction_limiter"] = pa[key]),
    "acnv_scaling"                 => (nl, pa, key) -> (nl["microphysics"]["microph_scaling_acnv"] = pa[key]),
    "accr_scaling"                 => (nl, pa, key) -> (nl["microphysics"]["microph_scaling_accr"] = pa[key]),
    "evap_scaling"                 => (nl, pa, key) -> (nl["microphysics"]["microph_scaling"] = pa[key]),
    "depsub_scaling"               => (nl, pa, key) -> (nl["microphysics"]["microph_scaling_dep_sub"] = pa[key]),
    "melt_scaling"                 => (nl, pa, key) -> (nl["microphysics"]["microph_scaling_melt"] = pa[key]),
    "thermo_covariance_model"      => (nl, pa, key) -> (nl["thermodynamics"]["thermo_covariance_model"] = pa[key]),
    "config"                       => (nl, pa, key) -> (nl["config"] = pa[key]),
    "set_src_seed"                 => (nl, pa, key) -> (nl["set_src_seed"] = pa[key]),
    "test_duals"                   => (nl, pa, key) -> (nl["test_duals"] = pa[key]),
    "float_type"                   => (nl, pa, key) -> (nl["float_type"] = pa[key]),
    )
    no_overwrites = (
        "job_id",
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
            @warn "Parameter `$key` overwriting namelist, `$key` = $(parsed_args[key])"
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
    return nothing
end
