"""
    name_aliases()

Returns a `Dict` containing:
 - a key (`String`), which we consider to be our core variable name
 - values (`Tuple` of `String`s), which we consider to be aliases of our core variable name
"""
function name_aliases()
    dict = Dict(
        "zc" => ("z_half",),
        "zf" => ("z",),
        "α0_c" => ("alpha0_half",),
        "α0_f" => ("alpha0",),
        "p0_c" => ("p0_half",),
        "p0_f" => ("p0",),
        "ρ0_c" => ("rho0_half",),
        "ρ0_f" => ("rho0",),
        "updraft_area" => ("updraft_fraction",),
        "updraft_thetal" => ("updraft_thetali",),
        "thetal_mean" => ("thetali_mean",),
    )
    return dict
end

"""
    get_nc_data(ds::NCDatasets.Dataset, var::String)

Returns the data for variable `var`, and also tries it's aliases
defined in `name_aliases`, in the `ds::NCDatasets.Dataset`.
"""
function get_nc_data(ds, var::String)
    dict = name_aliases()
    key_options = haskey(dict, var) ? (var, dict[var]...) : (var,)

    for key in key_options
        if haskey(ds, key)
            return ds[key]
        else
            for group_option in ["profiles", "reference", "timeseries"]
                haskey(ds.group, group_option) || continue
                if haskey(ds.group[group_option], key)
                    return ds.group[group_option][key]
                end
            end
        end
    end
    return nothing
end
