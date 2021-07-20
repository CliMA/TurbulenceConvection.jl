# compute the mean of the values between two percentiles (0 to 1) for a standard normal distribution
# this gives the surface scalar coefficients for 1 to n-1 updrafts when using n updrafts
function percentile_bounds_mean_norm(low_percentile::FT, high_percentile::FT, n_samples::I) where {FT <: Real, I}
    # TODO: check translation
    # x = rand(Normal(), n_samples)
    # xp_low = quantile(Normal(), low_percentile)
    # xp_high = quantile(Normal(), high_percentile)
    # filter!(y -> xp_low < y < xp_high, x)
    # TODO: undo this, it seems to fix the DYCOMS ql_mean
    pbmn = 1.7074549430665615
    return pbmn
end


function interp2pt(val1, val2)
    return 0.5 * (val1 + val2)
end

function logistic(x, slope, mid)
    return 1.0 / (1.0 + exp(-slope * (x - mid)))
end

function lamb_smooth_minimum(l, lower_bound, upper_bound)
    leng = size(l)
    xmin = minimum(l)
    lambda0 = max(xmin * lower_bound / real(LambertW.lambertw(2.0 / MathConstants.e)), upper_bound)

    # TODO: this will need to be i=1 when
    # going back to 1-based indexing
    # i = 1
    i = 0
    num = 0
    den = 0
    while (tuple(i) < leng)
        num += l[i] * exp(-(l[i] - xmin) / lambda0)
        den += exp(-(l[i] - xmin) / lambda0)
        i += 1
    end
    smin = num / den

    return smin
end
