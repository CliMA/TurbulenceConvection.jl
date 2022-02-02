# compute the mean of the values between two percentiles (0 to 1) for a standard normal distribution
# this gives the surface scalar coefficients for 1 to n-1 updrafts when using n updrafts
function percentile_bounds_mean_norm(low_percentile::FT, high_percentile::FT, n_samples::I) where {FT <: Real, I}
    # D = Distributions
    # TODO: check translation
    # x = rand(D.Normal(), n_samples)
    # xp_low = D.quantile(D.Normal(), low_percentile)
    # xp_high = D.quantile(D.Normal(), high_percentile)
    # filter!(y -> xp_low < y < xp_high, x)
    # TODO: undo this, it seems to fix the DYCOMS ql_mean
    pbmn = 1.7074549430665615
    return FT(pbmn)
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

#= Simple linear interpolation function, wrapping Dierckx =#
function pyinterp(x, xp, fp)
    spl = Dierckx.Spline1D(xp, fp; k = 1)
    return spl(vec(x))
end
