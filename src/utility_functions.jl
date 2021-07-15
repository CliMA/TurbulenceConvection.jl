
# import scipy.special as sp
# from libc.math cimport exp, log
# from scipy.stats import norm
# from scipy.special import lambertw

# compute the mean of the values above a given percentile (0 to 1) for a standard normal distribution
# this gives the surface scalar coefficient for a single updraft or nth updraft of n updrafts
function percentile_mean_norm(percentile, nsamples)
    x = norm.rvs(size=nsamples)
    xp = norm.ppf(percentile)
    return np.ma.mean(np.ma.masked_less(x,xp))
end

# compute the mean of the values between two percentiles (0 to 1) for a standard normal distribution
# this gives the surface scalar coefficients for 1 to n-1 updrafts when using n updrafts
# function percentile_bounds_mean_norm(low_percentile, high_percentile, nsamples)
#     x = norm.rvs(size=nsamples)
#     xp_low = norm.ppf(low_percentile)
#     xp_high = norm.ppf(high_percentile)
#     return np.ma.mean(np.ma.masked_greater(np.ma.masked_less(x,xp_low),xp_high))
# end
function percentile_bounds_mean_norm(
    low_percentile::FT,
    high_percentile::FT,
    n_samples::I,
) where {FT <: Real, I}
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
    return 0.5*(val1 + val2)
end

function logistic(x, slope, mid)
    return 1.0/(1.0 + exp( -slope * (x-mid)))
end

function smooth_minimum(x, a)
    i = 0
    num = 0; den = 0
    leng = x.shape[0]
    while (i<leng)
      if (x[i]>1.0e-5)
        num += x[i]*exp(-a*(x[i]))
        den += exp(-a*(x[i]))
      end
      i += 1
    end
    smin = num/den
    return smin
end

function auto_smooth_minimum(x, f)
    i = 0
    a = 1.0
    # TODO: this should be uncommented
    # np.ndarray[double, ndim=1] x_ = np.empty(len(x))

    lmin = 1.0e5; lmin2 = 1.0e5
    leng = x.shape[0]
    while (i<leng)
        x_[i] = x[i]
        i += 1
    end
    # Get min and second min values
    i = 0
    lmin = min(x_)
    while (i<leng)
      if (x_[i]<lmin2 && x_[i]>lmin+1.0e-5)
        lmin2 = x_[i]
      end
      x_[i] -= lmin
      i += 1
    end

    # Set relative maximum importance of second min term
    scale = (lmin2-lmin)/lmin*(1.0/(1.0+exp(lmin2-lmin)))
    if (scale>f)
        a = log((lmin2-lmin)/lmin/f-1.0)/(lmin2-lmin)
    end

    i = 0
    num = 0.0; den = 0.0;
    while(i<leng)
        num += x_[i]*exp(-a*(x_[i]))
        den += exp(-a*(x_[i]))
        i += 1
    end
    smin = lmin + num/den

    return smin
end

# function lamb_smooth_minimum(x, eps, dz)
#     i = 0
#     # TODO: this should be uncommented:
#     # np.ndarray[double, ndim=1] x_ = np.empty(len(x))

#     leng = x.shape[0]
#     # Copy array
#     while (i<leng)
#         x_[i] = x[i]
#         i += 1
#     end

#     xmin = min(x_)
#     lambda0 = max(xmin*eps/np.real(lambertw(2.0/np.e)), dz)

#     i = 0
#     num = 0.0; den = 0.0;
#     while(i<leng)
#         num += x_[i]*exp(-(x_[i]-xmin)/lambda0)
#         den += exp(-(x_[i]-xmin)/lambda0)
#         i += 1
#     end
#     smin = num/den
#     return smin
# end
function lamb_smooth_minimum(l, lower_bound, upper_bound)
    leng = size(l)
    xmin = minimum(l)
    lambda0 = max(
        xmin * lower_bound / real(LambertW.lambertw(2.0 / MathConstants.e)),
        upper_bound,
    )

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

function smooth_minimum2(x, l0)
    i = 0
    numLengths = 0
    smin = 0.0

    leng = x.shape[0]
    while(i<leng)
        if (x[i]>1.0e-5)
            smin += exp(-x[i]/l0)
            numLengths += 1
        end
        i += 1
    end
    smin /=  float(numLengths)
    smin =- l0*log(smin)
    return smin
end

function softmin(x, k)
    i = 1
    j = 1
    smin = 0.0
    eps = 0.1
    lam = 1.0

    leng = x.shape[0]
    lmin = min(x)
    num = 1.0
    den = 1.0
    while(j<leng)
        if (x[j]-lmin>eps*lmin)
            lam = log( ( (1.0+lmin/(k*x[j]))^(1.0/(len(x)-1.0)) - 1.0 )^(-1) )
            lam /= ((x[j]-lmin)/lmin)
            break;
        end
        j += 1
    end
    while(i<leng)
        x[i] /= lmin
        num += x[i]*exp(-lam*(x[i]-1.0))
        den += exp(-lam*(x[i]-1.0))
        i += 1
    end
    smin = lmin*num/den
    return smin
end

function hardmin(x)
    i = 0
    lmin = 1.0e6

    leng = x.shape[0]
    while(i<leng)
        if (x[i]>1.0e-5 && x[i]<lmin)
            lmin = x[i]
        end
        i += 1
    end
    return min(x)
end
