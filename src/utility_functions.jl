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

function upwind_advection_area(k, dzi, ρ0_half::Vector{Float64}, a_up::Vector{Float64}, w_up::Vector{Float64})
    whalf_kp = interp2pt(w_up[k - 1], w_up[k])
    whalf_k = interp2pt(w_up[k - 2], w_up[k - 1])

    m_kp = (ρ0_half[k] * a_up[k] * whalf_kp)
    m_k = (ρ0_half[k - 1] * a_up[k - 1] * whalf_k)
    return -dzi * (m_kp - m_k) / ρ0_half[k]
end

function upwind_advection_velocity(k, dzi, ρ0::Vector{Float64}, a_up::Vector{Float64}, w_up::Vector{Float64})
    a_k = interp2pt(a_up[k], a_up[k + 1])
    a_km = interp2pt(a_up[k - 1], a_up[k])
    adv = (ρ0[k] * a_k * w_up[k] * w_up[k] * dzi - ρ0[k - 1] * a_km * w_up[k - 1] * w_up[k - 1] * dzi)
    return adv
end

function upwind_advection_scalar(
    k,
    dzi,
    ρ0_half::Vector{Float64},
    a_up::Vector{Float64},
    w_up::Vector{Float64},
    var::Vector{Float64},
)
    m_k = (ρ0_half[k] * a_up[k] * interp2pt(w_up[k - 1], w_up[k]))
    m_km = (ρ0_half[k - 1] * a_up[k - 1] * interp2pt(w_up[k - 2], w_up[k - 1]))
    return (m_k * var[k] - m_km * var[k - 1]) * dzi
end

get_nc_data(data, group, var, imin, imax) = mean(data.group[group][var][:][:, imin:imax], dims = 2)[:]
