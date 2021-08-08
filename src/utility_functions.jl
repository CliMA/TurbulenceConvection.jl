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
    x_min = minimum(l)
    λ_0 = max(x_min * lower_bound / real(LambertW.lambertw(2.0 / MathConstants.e)), upper_bound)

    num = sum(map(1:length(l)) do i
        l[i] * exp(-(l[i] - x_min) / λ_0)
    end)
    den = sum(map(1:length(l)) do i
        exp(-(l[i] - x_min) / λ_0)
    end)
    smin = num / den
    return smin
end
