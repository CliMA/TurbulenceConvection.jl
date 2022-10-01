# compute the mean of the values between two percentiles (0 to 1) for a standard normal distribution
# this gives the surface scalar coefficients for 1 to n-1 updrafts when using n updrafts
function percentile_bounds_mean_norm(
    low_percentile::FT,
    high_percentile::FT,
    n_samples::Int,
    set_src_seed,
) where {FT <: Real}
    D = Distributions
    set_src_seed && Random.seed!(123)
    ∑samples = zero(FT)
    xp_low = D.quantile(D.Normal(), low_percentile)
    xp_high = D.quantile(D.Normal(), high_percentile)
    n_filtered_samples = 0
    @inbounds for i in 1:n_samples
        x = rand(D.Normal())
        if xp_low < x < xp_high
            ∑samples += x
            n_filtered_samples += 1
        end
    end
    return ∑samples / n_filtered_samples
end

function logistic(x, slope, mid)
    return 1 / (1 + exp(-slope * (x - mid)))
end

# lambert_2_over_e(::Type{FT}) where {FT} = FT(LambertW.lambertw(FT(2) / FT(MathConstants.e)))
lambert_2_over_e(::Type{FT}) where {FT} = FT(0.46305551336554884) # since we can evaluate

function lamb_smooth_minimum(l::SA.SVector, lower_bound::FT, upper_bound::FT) where {FT}
    x_min = minimum(l)
    λ_0 = max(x_min * lower_bound / lambert_2_over_e(FT), upper_bound)

    num = sum(l_i -> l_i * exp(-(l_i - x_min) / λ_0), l)
    den = sum(l_i -> exp(-(l_i - x_min) / λ_0), l)
    smin = num / den
    return smin
end

mean_nc_data(data, group, var, imin, imax) = StatsBase.mean(data.group[group][var][:][:, imin:imax], dims = 2)[:]
init_nc_data(data, group, var) = data.group[group][var][:][:, 1]

#= Simple linear interpolation function, wrapping Dierckx =#
function pyinterp(x, xp, fp)
    spl = Dierckx.Spline1D(xp, fp; k = 1)
    return spl(vec(x))
end

"""
    compare(a::Field, a::Field)
    compare(a::FieldVector, b::FieldVector)

Recursively compare two identically structured
`Field`s, or `FieldVector`s, with potentially different
data. If `!(maximum(abs.(parent(a) .- parent(b))) == 0.0)`
for any single field, then `compare` will print out the fields
and `err`. This can be helpful for debugging where and why
two `Field`/`FieldVector`s are different.
"""
function compare(a::FV, b::FV, pn0 = "") where {FV <: CC.Fields.FieldVector}
    for pn in propertynames(a)
        pa = getproperty(a, pn)
        pb = getproperty(b, pn)
        compare(pa, pb, "$pn0.$pn")
    end
end

function compare(a::F, b::F, pn0 = "") where {F <: CC.Fields.Field}
    if isempty(propertynames(a))
        err = abs.(parent(a) .- parent(b))
        if !(maximum(err) == 0.0)
            println("--- Comparing field $pn0")
            @show a
            @show b
            @show err
        end
    else
        for pn in propertynames(a)
            pa = getproperty(a, pn)
            pb = getproperty(b, pn)
            compare(pa, pb, "$pn0.$pn")
        end
    end
end
### Utility functions for LES_driven_SCM cases

"""
    get_LES_library

Hierarchical dictionary of available cfSite LES simulations, following similar
forcing to that described in [Shen2022](@cite), but including additional sites.
The following cfsites are available across listed models, months,
and experiments.
"""
function get_LES_library()
    LES_library = get_shallow_LES_library()
    deep_sites = (collect(30:33)..., collect(66:70)..., 82, 92, 94, 96, 99, 100)

    append!(LES_library["HadGEM2-A"]["07"]["cfsite_numbers"], deep_sites)
    append!(LES_library["HadGEM2-A"]["01"]["cfsite_numbers"], deep_sites)
    sites_04 = deepcopy(setdiff(deep_sites, [32, 92, 94]))
    append!(LES_library["HadGEM2-A"]["04"]["cfsite_numbers"], sites_04)
    sites_10 = deepcopy(setdiff(deep_sites, [94, 100]))
    append!(LES_library["HadGEM2-A"]["10"]["cfsite_numbers"], sites_10)

    LES_library_full = deepcopy(LES_library)
    for model in keys(LES_library_full)
        for month in keys(LES_library_full[model])
            LES_library_full[model][month]["cfsite_numbers"] = Dict()
            for cfsite_number in LES_library[model][month]["cfsite_numbers"]
                cfsite_number_str = string(cfsite_number, pad = 2)
                LES_library_full[model][month]["cfsite_numbers"][cfsite_number_str] = if cfsite_number >= 30
                    "deep"
                else
                    "shallow"
                end
            end
        end
    end
    return LES_library_full
end

"""
    get_shallow_LES_library

Hierarchical dictionary of available LES simulations described in [Shen2022](@cite).
The following cfsites are available across listed models, months,
and experiments.
"""
function get_shallow_LES_library()
    LES_library = Dict("HadGEM2-A" => Dict(), "CNRM-CM5" => Dict(), "CNRM-CM6-1" => Dict())
    Shen_et_al_sites = collect(2:15)
    append!(Shen_et_al_sites, collect(17:23))

    # HadGEM2-A model (76 AMIP-AMIP4K pairs)
    LES_library["HadGEM2-A"]["10"] = Dict()
    LES_library["HadGEM2-A"]["10"]["cfsite_numbers"] = Shen_et_al_sites
    LES_library["HadGEM2-A"]["07"] = Dict()
    LES_library["HadGEM2-A"]["07"]["cfsite_numbers"] = deepcopy(Shen_et_al_sites)
    LES_library["HadGEM2-A"]["04"] = Dict()
    LES_library["HadGEM2-A"]["04"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [15, 17, 18])
    LES_library["HadGEM2-A"]["01"] = Dict()
    LES_library["HadGEM2-A"]["01"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [15, 17, 18, 19, 20])

    # CNRM-CM5 model (59 AMIP-AMIP4K pairs)
    LES_library["CNRM-CM5"]["10"] = Dict()
    LES_library["CNRM-CM5"]["10"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [15, 22, 23])
    LES_library["CNRM-CM5"]["07"] = Dict()
    LES_library["CNRM-CM5"]["07"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [13, 14, 15, 18])
    LES_library["CNRM-CM5"]["04"] = Dict()
    LES_library["CNRM-CM5"]["04"]["cfsite_numbers"] =
        setdiff(Shen_et_al_sites, [11, 12, 13, 14, 15, 17, 18, 21, 22, 23])
    LES_library["CNRM-CM5"]["01"] = Dict()
    LES_library["CNRM-CM5"]["01"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [14, 15, 17, 18, 19, 20, 21, 22, 23])

    # CNRM-CM6-1 model (69 AMIP-AMIP4K pairs)
    LES_library["CNRM-CM6-1"]["10"] = Dict()
    LES_library["CNRM-CM6-1"]["10"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [22, 23])
    LES_library["CNRM-CM6-1"]["07"] = Dict()
    LES_library["CNRM-CM6-1"]["07"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [12, 13, 14, 15, 17])
    LES_library["CNRM-CM6-1"]["04"] = Dict()
    LES_library["CNRM-CM6-1"]["04"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [13, 14, 15])
    LES_library["CNRM-CM6-1"]["01"] = Dict()
    LES_library["CNRM-CM6-1"]["01"]["cfsite_numbers"] = setdiff(Shen_et_al_sites, [14, 15, 21, 22, 23])

    for month in ["01", "04", "07", "10"]
        LES_library["HadGEM2-A"][month]["experiments"] = ["amip", "amip4K"]
        LES_library["CNRM-CM5"][month]["experiments"] = ["amip", "amip4K"]
        LES_library["CNRM-CM6-1"][month]["experiments"] = ["amip", "amip4K"]
    end
    return LES_library
end


"""
    parse_les_path(les_path)
Given path to LES stats file, return cfsite_number, forcing_model, month, and experiment from filename.
    LES filename should follow pattern: `Stats.cfsite<SITE-NUMBER>_<FORCING-MODEL>_<EXPERIMENT>_2004-2008.<MONTH>.nc`
Inputs:
 - les_path - path to les simulation containing stats folder
Outputs:
 - cfsite_number  :: cfsite number
 - forcing_model :: {"HadGEM2-A", "CNRM-CM5", "CNRM-CM6-1", "IPSL-CM6A-LR"} - name of climate model used for forcing. Currently, only "HadGEM2-A" simulations are available reliably.
 - month :: {1, 4, 7, 10} - month of simulation.
 - experiment :: {"amip", "amip4K"} - experiment from which LES was forced.
"""
function parse_les_path(les_path)
    fname = basename(les_path)
    fname_split = split(fname, ('.', '_'))
    forcing_model = fname_split[3]
    experiment = fname_split[4]
    month = parse(Int64, fname_split[6])
    cfsite_number = parse(Int64, replace(fname_split[2], "cfsite" => ""))
    return (cfsite_number, forcing_model, month, experiment)
end

function valid_lespath(les_path)
    cfsite_number, forcing_model, month, experiment = parse_les_path(les_path)
    month = string(month, pad = 2)
    cfsite_number = string(cfsite_number, pad = 2)
    LES_library = get_LES_library()
    @assert forcing_model in keys(LES_library) "Forcing model $(forcing_model) not valid."
    @assert month in keys(LES_library[forcing_model]) "Month $(month) not available for $(forcing_model)."
    @assert cfsite_number in keys(LES_library[forcing_model][month]["cfsite_numbers"]) "cfSite $(cfsite_number) not found for $(forcing_model), month $(month)."
    @assert experiment in LES_library[forcing_model][month]["experiments"]
end
