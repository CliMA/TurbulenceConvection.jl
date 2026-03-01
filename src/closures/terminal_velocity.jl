#=
Test reimplimenting terminval_velocity so we don't have to update CloudMicrophysics.jld2
We mostly want truncated terminal velocities at some maximum radius

taken from https://github.com/CliMA/CloudMicrophysics.jl/blob/v0.14.0/src/Microphysics1M.jl
=#


# CMP = CM.Parameters
# CT = CMT
# CO = CM.Common
# const ACMP = CMP.AbstractCloudMicrophysicsParameters #
# const ACMP = CMP.AbstractCloudMicrophysicsParameters #




# ============================================================================= #
"""
    my_Chen2022ice_coeffs(prs, ρ_i)

 - prs - set with model parameters
 - ρ_i - cloud ice density

Returns the coefficients from Appendix B, Table B3 in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
needed for snow and ice terminal velocity
"""
function my_Chen2022ice_coeffs(prs::ACMP, ρ_i::FT) where {FT <: Real}

    # These appear to be from table B2 for small ice crystals ( < 625 μm)
    # Technically for snow we should integrate 62.5 to to 625 μm, and switch to go to infinity for larger ice crystals (snow)
    As_1::FT = CMP.As_coeff_1_Ch2022(prs)
    As_2::FT = CMP.As_coeff_2_Ch2022(prs)
    As_3::FT = CMP.As_coeff_3_Ch2022(prs)
    Bs_1::FT = CMP.Bs_coeff_1_Ch2022(prs)
    Bs_2::FT = CMP.Bs_coeff_2_Ch2022(prs)
    Bs_3::FT = CMP.Bs_coeff_3_Ch2022(prs)
    Cs_1::FT = CMP.Cs_coeff_1_Ch2022(prs)
    Cs_2::FT = CMP.Cs_coeff_2_Ch2022(prs)
    Cs_3::FT = CMP.Cs_coeff_3_Ch2022(prs)
    Cs_4::FT = CMP.Cs_coeff_4_Ch2022(prs)
    Es_1::FT = CMP.Es_coeff_1_Ch2022(prs)
    Es_2::FT = CMP.Es_coeff_2_Ch2022(prs)
    Es_3::FT = CMP.Es_coeff_3_Ch2022(prs)
    Fs_1::FT = CMP.Fs_coeff_1_Ch2022(prs)
    Fs_2::FT = CMP.Fs_coeff_2_Ch2022(prs)
    Fs_3::FT = CMP.Fs_coeff_3_Ch2022(prs)
    Gs_1::FT = CMP.Gs_coeff_1_Ch2022(prs)
    Gs_2::FT = CMP.Gs_coeff_2_Ch2022(prs)
    Gs_3::FT = CMP.Gs_coeff_3_Ch2022(prs)

    As = As_1 * (log(ρ_i))^2 − As_2 * log(ρ_i) - As_3
    Bs = FT(1) / (Bs_1 + Bs_2 * log(ρ_i) + Bs_3 / sqrt(ρ_i))
    Cs = Cs_1 + Cs_2 * exp(Cs_3 * ρ_i) + Cs_4 * sqrt(ρ_i)
    Es = Es_1 - Es_2 * (log(ρ_i))^2 + Es_3 * sqrt(ρ_i)
    Fs = -exp(Fs_1 - Fs_2 * (log(ρ_i))^2 + Fs_3 * log(ρ_i))
    Gs = FT(1) / (Gs_1 + Gs_2 / (log(ρ_i)) - Gs_3 * log(ρ_i) / ρ_i)

    return (As, Bs, Cs, Es, Fs, Gs)
end


"""
    my_Chen2022snow_coeffs(prs, ρ_i)

 - prs - set with model parameters
 - ρ_i - cloud ice density

Returns the coefficients from Appendix B, Table B3 in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
needed for snow and ice terminal velocity

    Copied Table B3 for small ice
    Added Table B5 for large ie, 
"""
function my_Chen2022snow_coeffs(prs::ACMP, ρ_i::FT) where {FT <: Real}

    # These appear to be from table B2 for small ice crystals ( < 625 μm)
    # Technically for snow we should integrate 62.5 to to 625 μm, and switch to go to infinity for larger ice crystals (snow)

    # small ice  (Table B3)
    As_s::FT, Bs_s::FT, Cs_s::FT, Es_s::FT, Fs_s, Gs_s::FT = my_Chen2022ice_coeffs(prs, ρ_i)

    # Large Ice (Tables B5)
    As_1l::FT, As_2l::FT, As_3l::FT = 0.475897, 0.00231270, 1.12293
    Bs_1l::FT, Bs_2l::FT, Bs_3l::FT = 2.56289, 0.00513504, 0.608459
    Cs_1l::FT, Cs_2l::FT, Cs_3l::FT = 0.756064, 0.935922, 1.70952
    Es_1l::FT, Es_2l::FT, Es_3l::FT = 0.00639847, 0.00906454, 0.108232
    Fs_1l::FT, Fs_2l::FT, Fs_3l::FT = 0.515453, 0.0725042, 1.86810e19
    Gs_1l::FT, Gs_2l::FT, Gs_3l::FT = 2.65236, 0.00158269, 259.935
    Hs_1l::FT, Hs_2l::FT, Hs_3l::FT = 0.346044, 7.17829e-11, 1.24394e20

    As_l = -As_1l - As_2l * log(ρ_i) + As_3l * ρ_i^(-3 / 2)
    Bs_l = exp(-Bs_1l - Bs_2l * log(ρ_i)^2 + Bs_3l * log(ρ_i))
    Cs_l = exp(-Cs_1l + Cs_2l / log(ρ_i) - Cs_3l / ρ_i)
    Es_l = Es_1l + Es_2l * log(ρ_i) * sqrt(ρ_i) - Es_3l * sqrt(ρ_i)
    Fs_l = Fs_1l - Fs_2l * log(ρ_i) - Fs_3l * exp(-ρ_i)
    Gs_l = (Gs_1l + Gs_2l * log(ρ_i) * sqrt(ρ_i) + Gs_3l / sqrt(ρ_i))^(-1)
    Hs_l = -Hs_1l - Hs_2l * ρ_i^(5 / 2) - Hs_3l * exp(-ρ_i)

    return (As_s, Bs_s, Cs_s, Es_s, Fs_s, Gs_s), (As_l, Bs_l, Cs_l, Es_l, Fs_l, Gs_l, Hs_l) # small and large particles 
end

"""
    my_Chen2022_vel_coeffs(prs, precip_type, ρ)

 - prs - set with free parameters
 - precip_type - type for ice, rain or snow
 - ρ - air density

Returns the coefficients from Appendix B in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function my_Chen2022_vel_coeffs(prs::ACMP, ::CMT.RainType, ρ::FT) where {FT <: Real}

    ρ0::FT = CMP.q_coeff_rain_Ch2022(prs)
    a1::FT = CMP.a1_coeff_rain_Ch2022(prs)
    a2::FT = CMP.a2_coeff_rain_Ch2022(prs)
    a3::FT = CMP.a3_coeff_rain_Ch2022(prs)
    a3_pow::FT = CMP.a3_pow_coeff_rain_Ch2022(prs)
    b1::FT = CMP.b1_coeff_rain_Ch2022(prs)
    b2::FT = CMP.b2_coeff_rain_Ch2022(prs)
    b3::FT = CMP.b3_coeff_rain_Ch2022(prs)
    b_ρ::FT = CMP.b_rho_coeff_rain_Ch2022(prs)
    c1::FT = CMP.c1_coeff_rain_Ch2022(prs)
    c2::FT = CMP.c2_coeff_rain_Ch2022(prs)
    c3::FT = CMP.c3_coeff_rain_Ch2022(prs)

    q = exp(ρ0 * ρ)
    ai = (a1 * q, a2 * q, a3 * q * ρ^a3_pow)
    bi = (b1 - b_ρ * ρ, b2 - b_ρ * ρ, b3 - b_ρ * ρ)
    ci = (c1, c2, c3)

    # unit conversions
    aiu = ai .* 1000 .^ bi .* (2 .^ bi) # mm^-bi to m^-bi
    ciu = ci .* 1000 .* 2

    return (aiu, bi, ciu)
end

function my_Chen2022_vel_coeffs(prs::ACMP, ::CMT.IceType, ρ::FT) where {FT <: Real}

    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    _As::FT, _Bs::FT, _Cs::FT, _Es::FT, _Fs::FT, _Gs::FT = my_Chen2022ice_coeffs(prs, ρ_i)

    ai = (_Es * ρ^_As, _Fs * ρ^_As)
    bi = (_Bs + ρ * _Cs, _Bs + ρ * _Cs)
    ci = (FT(0), _Gs)
    # unit conversions [ mm to m and diameter to radius]
    aiu = ai .* 1000 .^ bi .* (2 .^ bi)
    ciu = ci .* 1000 .* 2

    return (aiu, bi, ciu)
end


function my_Chen2022_vel_coeffs(prs::ACMP, ::CMT.SnowType, ρ::FT) where {FT <: Real}

    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    (_As_s::FT, _Bs_s::FT, _Cs_s::FT, _Es_s::FT, _Fs_s::FT, _Gs_s::FT),
    (_As_l::FT, _Bs_l::FT, _Cs_l::FT, _Es_l::FT, _Fs_l::FT, _Gs_l::FT, _Hs_l::FT) = my_Chen2022snow_coeffs(prs, ρ_i)

    # == small ================================================== # Table B2
    ai_s = (_Es_s * ρ^_As_s, _Fs_s * ρ^_As_s)
    bi_s = (_Bs_s + ρ * _Cs_s, _Bs_s + ρ * _Cs_s)
    ci_s = (FT(0), _Gs_s)
    # unit conversions
    aiu_s = ai_s .* 1000 .^ bi_s .* (2 .^ bi_s)
    ciu_s = ci_s .* 1000 .* 2
    # =========================================================== #


    # == large ================================================== # Table B4
    ai_l = (_Bs_l * ρ^_As_l, _Es_l * ρ^_As_l * exp(_Hs_l * ρ))
    bi_l = (_Cs_l, _Fs_l)
    ci_l = (FT(0), _Gs_l)
    # unit conversions [ *2 for diameteter to radius...]
    aiu_l = ai_l .* 1000 .^ bi_l .* (2 .^ bi_l) # mm^-bi to m^-bi
    ciu_l = ci_l .* 1000 .* 2 # /mm to /m
    # =========================================================== #

    return (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l)
end


"""
    my_Chen2022_vel_add(a, b, c, λ, k)

 - a, b, c, - free parameters defined in Chen etl al 2022
 - λ - size distribution parameter
 - k - size distribution moment for which we compute the bulk fall speed

Returns the addends of the bulk fall speed of rain or ice particles
following Chen et al 2022 DOI: 10.1016/j.atmosres.2022.106171 in [m/s].
We are assuming exponential size distribution and hence μ=0.
"""
function my_Chen2022_vel_add(
    a::FT,
    b::FT,
    c::FT,
    λ::FT,
    k::FT;
    Dmax::FT = Inf,
    Dmin::FT = 0.0,
    μ::FT = FT(0),
) where {FT <: Real}
    δ = μ + k + FT(1)
    # return a * λ^δ * CM1.SF.gamma(b + δ) / (λ + c)^(b + δ) / CM1.SF.gamma(δ)

    upper_limit_num = Dmax * (λ + c) # Dmax should not be 0
    lower_limit_num = iszero(Dmin) ? FT(0) : Dmin * (λ + c) # avoid 0 * inf

    upper_limit_den = Dmax * λ
    lower_limit_den = iszero(Dmin) ? FT(0) : Dmin * λ

    return a * λ^δ / (λ + c)^(b + δ) * (CM1.SF.gamma(b + δ, upper_limit_num) - CM1.SF.gamma(b + δ, lower_limit_num)) /
           (CM1.SF.gamma(δ, upper_limit_den) - CM1.SF.gamma(δ, lower_limit_den))
    # this should be the result of truncating at Dmax, where we transition from rain to snow for example...
end




"""
# think this Γ(..., Dmax * (_λ+...)) - Γ(..., Dmin * (_λ+...)) is right for integrating between bounds

#= Before, λ became λ+c bc e^-λD from size distro and e^-cD from terminal velocity were combined
# for small and large ice we use the same ansatz form, so is the same outcome?
_λ_num = ??? λ+c still...?

I think the only difference in the formula is that instead of deriving λ alone, it's more complicated for snow
but in the exponential is still just only λ and c?
=#
"""
function my_Chen2022_vel_add_sno(
    t::FT,
    b::FT,
    aec::FT,
    mec::FT,
    κ::FT,
    k::FT,
    c::FT,
    λ::FT,
    Dmin::FT,
    Dmax::FT;
    μ::FT = FT(0),
) where {FT}

    k += μ # since we dont use δ here

    upper_limit_num = Dmax * (λ + c) # Dmax should not be 0
    lower_limit_num = iszero(Dmin) ? FT(0) : Dmin * (λ + c) # avoid 0 * inf
    upper_limit_den = Dmax * λ
    lower_limit_den = iszero(Dmin) ? FT(0) : Dmin * λ

    return t * (
        CM1.SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1, upper_limit_num) -
        CM1.SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1, lower_limit_num)
    ) / (CM1.SF.gamma(k + 1, upper_limit_den) - CM1.SF.gamma(k + 1, lower_limit_den))


    # Checking here
    # return my_Chen2022_vel_add(t, b, c, λ, k; Dmax=Dmax, Dmin=Dmin) # treat same as small ice one

end
# This is for \int_0^∞, but we want transitions
# my_Chen2022_vel_add_sno(t, b, aec, mec, κ, k; Dmax=Dmax) =
# t * CM1.SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1) /
# CM1.SF.gamma(k + 1)

# ============================================================================= #

"""
    nm_int(prs, precip, velo_scheme, ρ, q_)

 - `prs` - abstract set with Earth parameters
 - `precip` - a type for ice, rain or snow
 - `velo_scheme` - type for terminal velocity parameterization
 - `ρ` - air density
 - `q_` - rain or snow specific humidity

 Returns  ∫n(r)m(r) dr assuming n(r) = N_0 * D^μ e^(-D*λ), m(r) = χ_m * m_0 * (r/r_0)^(m_e + Δ_m)

 Should be:
    ∫ n(r)m(r) dr = -χ_m m_0 n_0 (1/(r_0*λ))^(m_e + Δ_m)) Γ(m_e + Δ_m + 1)
        though in principle we put bounds on the integral so it's not infinite...so it's really
    -χ_m m_0 n_0 (1/(r_0*λ))^(m_e + Δ_m)) [Γ(m_e + Δ_m + 1, Dmax*λ) - Γ(m_e + Δ_m + 1, Dmin*λ)]

    q = ∫ n(r)m(r) dr / ρ

 """
function int_nm_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    # _n0::FT = isnothing(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, Nt, precip) # use wanted Nt If given
    _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, precip, Nt, μ) # use wanted Nt If given
    λ = lambda(prs, precip, q, ρ, Nt, μ, Dmin, Dmax)

    return int_nm_dr(prs, precip, _n0, λ, ρ, μ; Dmin = Dmin, Dmax = Dmax)
end



function int_nm_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    n0::FT,
    λ::FT,
    ρ::FT,
    μ::FT;
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    _r0::FT = r0(prs, precip)
    _m0::FT = m0(prs, precip)
    _me::FT = me(prs, precip)
    _Δm::FT = Δm(prs, precip)
    _χm::FT = χm(prs, precip)


    upper_limit = Dmax * λ # Dmax should not be zero...
    lower_limit = iszero(Dmin) ? FT(0) : Dmin * λ # avoid 0*Inf


    μ_tot = μ + _me + _Δm
    parts = (
        _χm,
        _m0,
        n0,
        (1 / (_r0))^(_me + _Δm),
        λ^-(μ_tot + FT(1)), # this is the same as λ^(-μ) * λ^(-1)
        CM1.SF.gamma(μ_tot + FT(1), upper_limit) - CM1.SF.gamma(μ_tot + FT(1), lower_limit),
    )

    if any(iszero, parts)
        return FT(0) # avoid 0*Inf -- you could get th wrong answer if it's really close to but not quite 0, but at least you won't get NaN
    else
        out = -reduce(*, parts)
        return isnan(out) ? FT(0) : out # in case the inf was introduced later, default to 0... (is that bad should probably only happen with very small inputs so probably is fine.)
    end
    # return -_χm * _m0 * n0 * (CM1.SF.gamma(μ + FT(1), upper_limit) - CM1.SF.gamma(μ + FT(1), lower_limit))
end

int_q(
    prs::ACMP,
    precip::CMTWaterTypes,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN),
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT} = int_nm_dr(prs, precip, q, ρ, μ; Nt = Nt, Dmin = Dmin, Dmax = Dmax) / ρ
int_q(prs::ACMP, precip::CMTWaterTypes, n0::FT, λ::FT, ρ::FT, μ::FT; Dmin::FT = FT(0), Dmax::FT = FT(Inf)) where {FT} =
    int_nm_dr(prs, precip, n0, λ, ρ, μ; Dmin = Dmin, Dmax = Dmax) / ρ




"""
integrate n(r) from Dmin to Dmax
    n(r) = N_0 D^μ e^(-D*λ)

    ∫n(r) dr = -N_0 / λ^(μ+1) * (Γ(μ+1, Dmax*λ) - Γ(μ+1, Dmin*λ))
"""
function int_n_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, precip, Nt, μ; Dmin = Dmin, Dmax = Dmax) # use wanted Nt If given
    λ::FT = lambda(prs, precip, q, ρ, Nt, μ, Dmin, Dmax)
    return int_n_dr(prs, precip, _n0, λ, ρ, μ; Dmin = Dmin, Dmax = Dmax)
end

function int_n_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    n0::FT,
    λ::FT,
    ρ::FT,
    μ::FT;
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    _r0::FT = r0(prs, precip)

    upper_limit = Dmax * λ # Dmax should not be zero...
    lower_limit = iszero(Dmin) ? FT(0) : Dmin * λ # avoid 0*Inf


    parts = (
        n0,
        λ^-(μ + FT(1)), # this is the same as λ^(-μ) * λ^(-1)
        CM1.SF.gamma(μ + FT(1), upper_limit) - CM1.SF.gamma(μ + FT(1), lower_limit),
    )


    if any(iszero, parts)
        return FT(0) # avoid 0*Inf -- you could get th wrong answer if it's really close to but not quite 0
    else
        out = -reduce(*, parts)
        return isnan(out) ? FT(0) : out # in case the inf was introduced later
    end
end
n_int(
    prs::ACMP,
    precip::CMTWaterTypes,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN),
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT} = int_n_dr(prs, precip, q, ρ, μ; Nt = Nt, Dmin = Dmin, Dmax = Dmax) / ρ


"""
Integrate n(r) * r
    n(D) = N_0 D^μ e^(-D*λ)

    ∫n(r) r dr = ∫ N_0 D^(μ+1) e^(-D*λ) dr = N_0 / λ^(μ+2) * (Γ(μ+2, Dmax*λ) - Γ(μ+2, Dmin*λ))

"""
function int_nr_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    _n0 = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, precip, Nt, μ; Dmin = Dmin, Dmax = Dmax) # use wanted Nt If given
    λ = lambda(prs, precip, q, ρ, Nt, μ, Dmin, Dmax)

    return int_nr_dr(prs, precip, _n0, λ, ρ, μ; Dmin = Dmin, Dmax = Dmax)
end

function int_nr_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    n0::FT,
    λ::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    upper_limit = Dmax * λ # Dmax should not be zero...
    lower_limit = iszero(Dmin) ? FT(0) : Dmin * λ # avoid 0*Inf
    μ_tot = μ + FT(1) # we add 1 to the moment for x r
    parts = (
        n0,
        λ^-(μ + FT(1)), # this is the same as λ^(-μ) * λ^(-1)
        CM1.SF.gamma(μ_tot + FT(1), upper_limit) - CM1.SF.gamma(μ_tot + FT(1), lower_limit),
    )

    if any(iszero, parts)
        return FT(0) # avoid 0*Inf -- you could get th wrong answer if it's really close to but not quite 0
    else
        out = -reduce(*, parts)
        return isnan(out) ? FT(0) : out # in case the inf was introduced later
    end
end


"""
Integrate n(r)a(r) dr
    n(r) = N_0 D^μ e^(-D*λ)
    a(r) = χ_a a_0 (r/r_0)^(a_e + Δ_a)

    ∫n(r)a(r) dr = -χ_a a_0 N_0 / λ^(μ + a_e + Δ_a + 1) * (Γ(μ + a_e + Δ_a + 1, Dmax*λ) - Γ(μ + a_e + Δ_a + 1, Dmin*λ))
"""
function int_na_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, precip, Nt, μ; Dmin = Dmin, Dmax = Dmax) # use wanted Nt If given
    λ = lambda(prs, precip, q, ρ, Nt, μ, Dmin, Dmax)
    return int_na_dr(prs, precip, _n0, λ, ρ, μ; Dmin = Dmin, Dmax = Dmax)
end

function int_na_dr(
    prs::ACMP,
    precip::CMTWaterTypes,
    n0::FT,
    λ::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    _r0::FT = r0(prs, precip)
    _a0::FT = a0(prs, precip)
    _ae::FT = ae(prs, precip)
    _Δa::FT = Δa(prs, precip)
    _χa::FT = χa(prs, precip)


    upper_limit = Dmax * λ # Dmax should not be zero...
    lower_limit = iszero(Dmin) ? FT(0) : Dmin * λ # avoid 0*Inf

    μ_tot = μ + _ae + _Δa
    parts = (
        _χa,
        _a0,
        n0,
        (1 / (_r0))^(_ae + _Δa),
        λ^-(μ_tot + FT(1)), # this is the same as λ^(-μ) * λ^(-1)
        CM1.SF.gamma(μ_tot + FT(1), upper_limit) - CM1.SF.gamma(μ_tot + FT(1), lower_limit),
    )

    if any(iszero, parts)
        return FT(0) # avoid 0*Inf -- you could get th wrong answer if it's really close to but not quite 0
    else
        out = -reduce(*, parts)
        return isnan(out) ? FT(0) : out # in case the inf was introduced later
    end

end



"""
    Numerator of Eq 20 in Chen et al 2022 for the kth moment.
        = ∫_Dmin^Dmax V(D) D^k n(D) dD =  n0 ∑ -a / (λ + c)^(b + δ) * (CM1.SF.gamma(b + δ, upper_limit) - CM1.SF.gamma(b + δ, lower_limit))
        V(D) = ϕ^κ ∑(a_i * D_i^b_i * e^(-c_i * D))
        n(D) = n_0 D^μ e^(-λD)

        We ignore the ϕ^κ part, and just use the sum of the a_i * D_i^b_i * e^(-c_i * D) part [i.e. aspect ratio = 1]

        Then the integral should be as in Eq 20 in Chen et al 2022 but only The top part [ note they had a negative sign that cancelled out ]

        Helper for int_nav_dr
"""
function int__v_Dk_n__dD(
    a::FT,
    b::FT,
    c::FT,
    n0::FT,
    λ::FT,
    k::FT;
    Dmax::FT = FT(Inf),
    Dmin::FT = FT(0),
    μ::FT = FT(0),
) where {FT <: Real}
    # defaults to μ = 0,  Exponential (Marshall-Palmer) instead of Gamma distribution
    δ = μ + k + 1

    upper_limit = Dmax * (λ + c) # Dmax should not be 0
    lower_limit = iszero(Dmin) ? FT(0) : Dmin * (λ + c) # avoid 0 * inf

    return -a * n0 / (λ + c)^(b + δ) * (CM1.SF.gamma(b + δ, upper_limit) - CM1.SF.gamma(b + δ, lower_limit))
    # this should be the result of truncating at Dmax, where we transition from rain to snow for example...
end


"""
∫n(r)a(r)v(r) dr 
n(r) = N_0 D^μ e^(-D*λ)
a(r) = χ_a a_0 (r/r_0)^a_e = cross section
v(r) = χ_v v_0 (r/r_0)^(v_e + Δ_v) = velocity

∫ n(r)a(r)m(r) dr = -χ_m m_0 χ_a m_a n_0 (1/(r_0*λ))^(m_e + Δ_m + a_e + Δ_a)) Γ(m_e + Δ_m + a_e + Δ_a + 1)
    though in principle we put bounds on the integral so it's not infinite... so it's really
-χ_m m_0 χ_a m_a n_0 (1/(r_0*λ))^(m_e + Δ_m + a_e + Δ_a)) [Γ(m_e + Δ_m + m_a + Δ_a 1, Dmax*λ) - Γ(m_e + Δ_m + m_a + Δ_a + 1, Dmin*λ)]
"""
function int_nav_dr(
    # prs::ACMP,
    param_set::APS,
    precip::CMTWaterTypes,
    velo_scheme::CMT.Blk1MVelType,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    error("int_nav_dr not implemented yet for Blk1MVelType -- this fcn needs to bupdated with v(r) instead of a(r)")
end



"""
    Unlike Eq 5,20 in Chen et al 2022 [ https://doi.org/10.1016/j.atmosres.2022.106171 ], here we're not going for a mass-weighted average of the terminal velocities.
    Instead we just want the one integral ∫ n(r)a(r)v(r) dr

        Note we just also use ql instead of having m(r) for liquid in the integral like in accretion_snow_rain() and don't have m for ice at all. I'm not sure why that's ok.

    If we plug in a(r) = π (D/2)^2 things are simplest... 
    We really have a(r) = χ_a a_0 (r/r_0)^(a_e + Δ_a)

    Then we have v(r) = χ_v v_0 (r/r_0)^(v_e + Δ_v) = χ_v v_0 (D/r_0)^(v_e + Δ_v)

    ∫ n(r)a(r)v(r) dr = π/4 ∫ n(r)D^2v(r) dr so you just have k = 2


    From Eq 20 this means [ using V =  ϕ^κ ∑_i (a_i * D^b_i * e^(-c_i*D)) from Eq 19 ]
        ∫ n(r)a(r)v(r) dr = π/4 ϕ^κ n0 ∑_i (a_i * λ^δ Γ(b_i + δ) / (λ+c_i)^(b_i+δ)
"""
function int_nav_dr(
    # prs::ACMP,
    param_set::APS,
    precip::CMT.IceType,
    velo_scheme::CMT.Chen2022Type,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    D_transition::FT = FT((0.625e-3) / 2), # .625mm is the transition from small to large ice crystals in Chen paper
) where {FT <: Real}

    int::FT = FT(0)
    if q > FT(0)

        prs = TCP.microphysics_params(param_set)

        if isnan(μ)
            μ = μ_from_qN(param_set, precip, q, Nt; ρ = ρ)
        end

        # coefficients from Appendix B from Chen et. al. 2022
        (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l) = my_Chen2022_vel_coeffs(prs, CM.CommonTypes.SnowType(), ρ)

        # ρ_i::FT = CMP.ρ_cloud_ice(prs)
        # Ok so here's the thing.... we want q to fit between Dmin and Dmax...
        # so we need a lambda that works for that....
        # _λ::FT = lambda(param_set, precip, q, ρ, Nt, μ, Dmin, Dmax)
        _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, precip, Nt, μ; Dmin = Dmin, Dmax = Dmax) # use wanted Nt If given
        _λ::FT = lambda(prs, precip, q, ρ, Nt, μ; _n0 = _n0, Dmin = Dmin, Dmax = Dmax)

        a0c::FT = a0(prs, precip) * χa(prs, precip)
        aec::FT = ae(prs, precip) + Δa(prs, precip)
        _r0::FT = r0(prs, precip)


        if Dmin < D_transition
            if Dmax <= D_transition
                regions = ((Dmin, Dmax),)
                abcs = ((aiu_s, bi_s, ciu_s),)
            else
                regions = ((Dmin, D_transition), (D_transition, Dmax))
                abcs = ((aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l))
            end
        else
            regions = ((Dmin, Dmax),)
            abcs = ((aiu_l, bi_l, ciu_l),)
        end

        # Here we assume μ = 0, otherewise we would add it to k... (we had δ = μ + k + 1)

        Dmax = min(Dmax, FT(.625e-3) / 2) # v is I think poorly defined past that point... for Chen.. grows too big, so we upper bound it.


        # Here the integrals should be additive, so we can just sum them up without any weighting
        # k = 2 #  k = 2 here bc a ∝ r^2
        k = aec # should be a ∝ r^2
        # for ((Dmin, Dmax), (aiu, bi, ciu)) in zip(regions, abcs)
        #     # eq 20 from Chen et al 2022
        #     # int += sum(int__v_Dk_n__dD.(aiu, bi, ciu, _n0, _λ, k; Dmin = Dmin, Dmax = Dmax, μ=μ)) 
        #     int += sum(((aiu, bi, ciu) -> int__v_Dk_n__dD(aiu, bi, ciu, _n0, _λ, k; Dmin=Dmin, Dmax=Dmax, μ=μ)).(aiu, bi, ciu))
        # end
        # @inbounds for ((Dmin, Dmax), (aiu, bi, ciu)) in zip(regions, abcs) # zip allocates for some reason
        @inbounds for i in eachindex(regions)
            Dmin, Dmax = regions[i]
            aiu, bi, ciu = abcs[i]
            # Eq 20 from Chen et al 2022
            @inbounds for (aiu_j, bi_j, ciu_j) in zip(aiu, bi, ciu)
                int += int__v_Dk_n__dD(aiu_j, bi_j, ciu_j, _n0, _λ, k; Dmin = Dmin, Dmax = Dmax, μ = μ)
            end
        end

        # int *= π /4 # we need to multiply by π/4 to get the right answer for π(D/2)^2 = π/4 D^2 [CloudMicrophysics.j actually uses radius]
        int *= (a0c * _r0^-(aec)) * χv(prs, precip)

        int = max(zero(FT), int) # n, a, and v are all positive, so the integral should be positive
    end
    return resolve_nan(int)
end


function int_nav_dr(
    # prs::ACMP,
    param_set::APS,
    precip::CMT.RainType,
    velo_scheme::CMT.Chen2022Type,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}

    int = FT(0)
    if q > FT(0)

        prs = TCP.microphysics_params(param_set)
        aec = ae(prs, precip) + Δa(prs, precip)


        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = my_Chen2022_vel_coeffs(prs, precip, ρ)
        # size distribution parameter
        _λ::FT = lambda(param_set, precip, q, ρ, Nt, Dmin, Dmax; μ = μ)
        _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(param_set, q, ρ, precip, Nt; Dmin = Dmin, Dmax = Dmax, μ = μ)

        # what is n0 here?
        _χa::FT = χa(prs, precip)
        _a0::FT = a0(prs, precip)
        _r0::FT = r0(prs, precip)

        # eq 20 from Chen et al 2022
        # k = 2 here bc a ∝ r^2
        k = aec # should be right...
        int = sum(int__v_Dk_n__dD.(aiu, bi, ciu, _n0, _λ, k; Dmin = Dmin, Dmax = Dmax))
        # int *=  FT(π)/4 
        int *= _χa * _a0 * _r0^(-aec) # a0 is the cross section, r0 is the reference radius, so this is the same as π(D/2)^2 = π/4 D^2 [CloudMicrophysics.j actually uses radius]

        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        int = max(FT(0), int) # n, a, and v are all positive, so the integral should be positive
    end
    return resolve_nan(int)
end

function int_nav_dr(
    prs::ACMP,
    precip::CMT.SnowType,
    velo_scheme::CMT.Chen2022Type,
    q::FT,
    ρ::FT,
    μ::FT;
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    D_transition::FT = FT((0.625e-3) / 2), # .625mm is the transition from small to large ice crystals in Chen paper
) where {FT <: Real}
    error("Not implemented yet")
    # you would just look to terminal_velocity() for inspiration but we don't need this atm... so
end





"""
    terminal_velocity(prs, precip, velo_scheme, ρ, q_)

 - `prs` - abstract set with Earth parameters
 - `precip` - a type for ice, rain or snow
 - `velo_scheme` - type for terminal velocity parameterization
 - `ρ` - air density
 - `q_` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of particles.
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical power-law relations for `velo_scheme == Blk1MVelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171, for `velo_scheme == Chen2022Type`
"""
function my_terminal_velocity(
    # prs::ACMP,
    param_set::APS,
    precip::CMT.AbstractPrecipType,
    velo_scheme::CMT.Blk1MVelType,
    ρ::FT,
    q_::FT;
    Dmin::FT = FT(0.0), # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = FT(Inf), # not implemented yet
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    μ::FT = FT(NaN),
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        prs = TCP.microphysics_params(param_set)

        if isnan(μ)
            μ = μ_from_qN(param_set, precip, q_, Nt; ρ = ρ) # maybe we could have a joint get_μ_λ_from_qN() function
        end

        _r0::FT = r0(prs, precip)
        _me::FT = me(prs, precip)
        _Δm::FT = Δm(prs, precip)
        _χm::FT = χm(prs, precip)
        _χv::FT = χv(prs, precip)
        _v0::FT = v0(prs, ρ, precip)
        _ve::FT = ve(prs, precip)
        _Δv::FT = Δv(prs, precip)
        _λ::FT = lambda(param_set, precip, q_, ρ, Nt, Dmin, Dmax; μ = μ)

        upper_limit = Dmax * _λ # Dmax should not be zero...
        lower_limit = iszero(Dmin) ? FT(0) : Dmin * _λ # avoid 0*Inf

        fall_w =
            _χv *
            _v0 *
            (_λ * _r0)^(-_ve - _Δv) *
            (
                CM1.SF.gamma(_me + _ve + _Δm + _Δv + FT(1), upper_limit) -
                CM1.SF.gamma(_me + _Δm + _Δv + FT(1), lower_limit)
            ) / (CM1.SF.gamma(_me + _Δm + FT(1), upper_limit) - CM1.SF.gamma(_me + _Δm + FT(1), lower_limit))
        # CM1.SF.gamma(_me + _ve + _Δm + _Δv + FT(1)) / CM1.SF.gamma(_me + _Δm + FT(1))
    end

    return resolve_nan(fall_w)
end

function my_terminal_velocity(
    # prs::ACMP,
    param_set::APS,
    precip::CMT.RainType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = FT(0), # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = FT(Inf),
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    μ::FT = FT(NaN),
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)
        prs = TCP.microphysics_params(param_set)

        if isnan(μ)
            μ = μ_from_qN(param_set, precip, q_, Nt; ρ = ρ)
        end

        mec::FT = me(prs, precip) + Δm(prs, precip)
        k = mec

        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = my_Chen2022_vel_coeffs(prs, precip, ρ)
        # size distribution parameter
        _λ::FT = lambda(param_set, precip, q_, ρ, Nt, Dmin, Dmax; μ = μ)
        # eq 20 from Chen et al 2022
        fall_w = sum(my_Chen2022_vel_add.(aiu, bi, ciu, _λ, k; Dmin = Dmin, Dmax = Dmax, μ = μ))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)

        fall_w *= χv(prs, precip) # scaling_factor
    end
    return resolve_nan(fall_w)
end

"""
For liquid we will copy rain, but if Nt is not given, we will enforce a large N (250 / cm^3) bc the assumed rain params are not valid for small rain drops.
"""
function my_terminal_velocity(
    # prs::ACMP,
    param_set::APS,
    precip::CMT.LiquidType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = FT(0.0), # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = FT(Inf),
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
) where {FT <: Real}
    error(
        "my_terminal_velocity for LiquidType not implemented yet, use my_terminal_velocity for RainType instead, but be sure to pass in your own `N` as the default will be far too small for rain.",
    )
end


function my_terminal_velocity(
    # prs::ACMP,
    param_set::APS,
    precip::CMT.IceType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = FT(0), # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = FT((0.625e-3) / 2), # 625 μm but we're really in r units now
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
    D_transition::FT = FT((0.625e-3) / 2), # .625mm is the transition from small to large ice crystals in Chen paper
    μ::FT = FT(NaN),
) where {FT <: Real}
    fall_w::FT = FT(0)
    # if q_ > FT(0)
    if q_ > eps(FT) # idk if this matters but it makes the plots look prettier lol, also prevents N, q divergence problems... also once we shrink back to tiny particles, the speed should be slow... we maybe just guessed N wrong.


        #=
            This seems to get out of hand when <r> gets too large, surpassing what even the snow paramterization would produce.
            I think this is because the hard sphere assumption should have broken down but hasn't, and maybe Chen is too lax about it.

            So, if <r> > r_ice_snow, we'll re-index on r_ice_snow being true (we just reset lambda)
        =#


        prs = TCP.microphysics_params(param_set)

        ρ_i::FT = CMP.ρ_cloud_ice(prs)

        # Ok so here's the thing.... we want q to fit between Dmin and Dmax...
        # so we need a lambda that works for that....

        if isnan(μ)
            μ = μ_from_qN(param_set, precip, q_, Nt; ρ = ρ)
        end
        _λ::FT = lambda(param_set, precip, q_, ρ, Nt, Dmin, Dmax; μ = μ)

        microphys_params = TCP.microphysics_params(param_set)
        # r_is = CMP.r_ice_snow(microphys_params)
        # if _λ < (μ + 1)/r_is # too limiting, we def see w_i increases after <r_is>
        #     _λ = (μ + 1)/r_is # let's stop Chen from going overboard with the sedimentation rates... [[ idk if this would be good, we gotta go past r_is according to the LES so idk... but the w response can't be so extreme...]]
        # end
        # r_th = get_r_cond_precip(param_set, precip) * FT(param_set.user_params.r_ice_snow_threshold_scaling_factor)
        # if _λ < (μ + 1)/r_th # I think with the limits we put in on max r in the integration, we can go all the way up to r_th now... and maybe even not need a limit...
        #     _λ = (μ + 1)/r_th # let's stop Chen from going overboard with the sedimentation rates... [[ idk if this would be good, we gotta go past r_is according to the LES so idk... but the w response can't be so extreme...]]
        # end

        m0c::FT = m0(prs, precip) * χm(prs, precip)
        a0c::FT = a0(prs, precip) * χa(prs, precip)
        mec::FT = me(prs, precip) + Δm(prs, precip)
        aec::FT = ae(prs, precip) + Δa(prs, precip)
        _r0::FT = r0(prs, precip)

        # ================================================================= #

        # # coefficients from Appendix B from Chen et. al. 2022
        # aiu, bi, ciu = my_Chen2022_vel_coeffs(prs, precip, ρ)

        # # eq 20 from Chen et al 2022
        # fall_w = sum(my_Chen2022_vel_add.(aiu, bi, ciu, _λ, 3; Dmin=Dmin, Dmax=Dmax))
        # # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        # fall_w = max(FT(0), fall_w)

        # ================================================================= #
        # k = 3 # mass weighted
        k = mec
        # k += μ # add the size distribution exponent to the moment [[ if you don't do this, you need to add it in the mass weights ]]

        # coefficients from Appendix B from Chen et. al. 2022
        (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l) = my_Chen2022_vel_coeffs(prs, CM.CommonTypes.SnowType(), ρ)

        # local mass_weights::SA.MVector{2, FT}
        # mass_weights = SA.@MVector [FT(0), FT(0)]
        if Dmin < D_transition
            if Dmax <= D_transition
                regions = ((Dmin, Dmax),)
                abcs = ((aiu_s, bi_s, ciu_s),)
            else
                regions = ((Dmin, D_transition), (D_transition, Dmax))
                abcs = ((aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l))
            end
        else
            regions = ((Dmin, Dmax),)
            abcs = ((aiu_l, bi_l, ciu_l),)
        end

        # for (i, ((Dmin, Dmax), (aiu, bi, ciu))) in enumerate(zip(regions, abcs)) # # this allocates for some reaonwe basically need to sum the integral as before but over all regions
        #     # are the mass weights really necessary? It's just an additive integral, no?
        #     mass_weights[i] = _λ^-(μ + k + 1) * (-CM1.SF.gamma(μ + k + 1, Dmax * _λ) + CM1.SF.gamma(μ + k + 1, Dmin * _λ)) # missing constants from the integral (n_0, 4/3, π, etc) but those are all the same and cancel out

        #     # eq 20 from Chen et al 2022
        #     fall_w += sum(my_Chen2022_vel_add.(aiu, bi, ciu, _λ, k; Dmin = Dmin, Dmax = Dmax, μ=μ)) .* mass_weights[i] # I think it's supposed to be k not FT(3) here
        # end
        total_weight = FT(0)
        @inbounds for i in eachindex(regions)
            Dmin, Dmax = regions[i]
            aiu, bi, ciu = abcs[i]

            mass_weight = _λ^-(μ + k + 1) * (-CM1.SF.gamma(μ + k + 1, Dmax * _λ) + CM1.SF.gamma(μ + k + 1, Dmin * _λ)) # missing constants from the integral (n_0, 4/3, π, etc) but those are all the same and cancel out
            total_weight += mass_weight

            @inbounds for (aiu_j, bi_j, ciu_j) in zip(aiu, bi, ciu)
                fall_w += my_Chen2022_vel_add(aiu_j, bi_j, ciu_j, _λ, k; Dmin = Dmin, Dmax = Dmax, μ = μ) .* mass_weight
            end
        end

        if !iszero(total_weight)
            fall_w /= total_weight # normalize by total mass
        end


        # ================================================================= #
        fall_w *= χv(prs, precip) # scaling_factor

    end


    return resolve_nan(fall_w)
end




"""
Adjusted for snow...

We go to the from Dmin to Dmax, with a break at the transition point if necessary

Is this smooth?
Or are we supposed to integrate both the whole way and the gamma functions just decay as necessary? idk...

"""
# function my_terminal_velocity(
#     # prs::ACMP,
#     param_set::APS,
#     precip::CMT.SnowType,
#     velo_scheme::CMT.Chen2022Type,
#     ρ::FT,
#     q_::FT;
#     Dmin::FT = FT(0), # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
#     Dmax::FT = FT(Inf), # not really needed for snow? maybe a Dmin? but that wouldn't make a huge diff i dont think idk...
#     # Nt::Union{FT, Nothing} = nothing,
#     Nt::FT = FT(NaN), # testing for type stability, use NaN instead of nothing
#     D_transition::FT = FT((.625e-3)/2), # .625mm is the transition from small to large ice crystals in Chen paper
#     μ::FT = FT(NaN), # μ is the exponent in the PSD, default to 0 for marshall palmer
# ) where {FT <: Real}
#     fall_w = FT(0)
#     if q_ > FT(0)

#         prs = TCP.microphysics_params(param_set)

#         if isnan(μ)
#             μ = μ_from_qN(param_set, precip, q_, Nt; ρ = ρ)
#         end

#         # short circuit cause greater than .625 mm is giving negative values for some reason... so we'll just ignore that I guess...
#         # return my_terminal_velocity(prs, ice_type, velo_scheme, ρ, q_; Dmin=Dmin, Dmax=Dmax, Nt=Nt) # backup for now bc large ice is broken...

#         _r0::FT = r0(prs, precip)
#         _λ::FT = lambda(param_set, precip, q_, ρ, Nt, Dmin, Dmax; μ = μ)
#         m0c::FT = m0(prs, precip) * χm(prs, precip)
#         a0c::FT = a0(prs, precip) * χa(prs, precip)
#         mec::FT = me(prs, precip) + Δm(prs, precip)
#         aec::FT = ae(prs, precip) + Δa(prs, precip)

#         ρ_i::FT = CMP.ρ_cloud_ice(prs)

#         # D_transition::FT = 0.625e-3 # .625mm is the transition from small to large ice crystals in Chen paper

#         local mass_weights::SA.MVector{2, FT}

#         # coefficients from Appendix B from Chen et. al. 2022
#         (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l) = my_Chen2022_vel_coeffs(prs, precip, ρ)

#         # κ = FT(-1 / 3) #oblate
#         # κ = FT(1 / 3) # oblate (anna says this is right) [[ lol i think she maybe is wright lol...]]
#         κ = FT(-1 / 6) # prolate (anna says this is right) [[ lol i think she maybe is wright lol...]]

#         # k = 3 # mass weighted
#         k = mec # should still be mass weight
#         # κ = FT(0)

#         @assert(Dmin <= Dmax, "Dmin must be less than Dmax, got Dmin = $Dmin and Dmax = $Dmax")

#         # Is this smooth?
#         # or are we supposed to integrate both the whole way and the gamma functions just decay as necessary? idk...
#         mass_weights = SA.@MVector [FT(0), FT(0)]
#         if Dmin < D_transition
#             if Dmax <= D_transition
#                 regions = ((Dmin, Dmax),)
#                 abcs = ((aiu_s, bi_s, ciu_s),)
#             else
#                 regions = ((Dmin, D_transition), (D_transition, Dmax))
#                 abcs = ((aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l))
#             end
#         else
#             regions = ((Dmin, Dmax),)
#             abcs = ((aiu_l, bi_l, ciu_l),)
#         end

#         tmp = _λ^(k + 1) * ((16 * a0c^3 * ρ_i^2) / (9 * π * m0c^2 * _r0^(3 * aec - 2 * mec)))^κ # I got this from Anna for aspect ratio... not fuly sure if it's right tho

#         fall_w = FT(0)
#         for (i, ((Dmin, Dmax), (aiu, bi, ciu))) in enumerate(zip(regions, abcs)) # we basically need to sum the integral as before but over all regions
#             ci_pow = (2 .* ciu .+ _λ) .^ (.-(3 .* aec .* κ .- 2 .* mec .* κ .+ bi .+ k .+ 1))

#             # ti = tmp .* aiu .* FT(2) .^ bi .* ci_pow
#             ti = tmp .* aiu .* ci_pow

#             # ti = _λ^(k + 1) .* aiu;

#             upper_limit = Dmax * _λ # Dmax should not be zero...
#             lower_limit = iszero(Dmin) ? FT(0) : Dmin * _λ # avoid 0*Inf

#             # k = 3 for mass, μ = 0 for exponential size distribution instead of gamma
#             # k = 3 means we're already mass weighting, however we need to properly mass weight in between regions. For that we use the mass_weights
#             # mass_weights[i] = _λ^-(k + 1) * (-CM1.SF.gamma(k + 1, upper_limit) + CM1.SF.gamma(k + 1, lower_limit)) # missing constants from the integral (n_0, 4/3, π, etc) but those are all the same and cancel out
#             mass_weights[i] = _λ^-(μ + k + 1) * (-CM1.SF.gamma(μ + k + 1, upper_limit) + CM1.SF.gamma(μ + k + 1, lower_limit)) # missing constants from the integral (n_0, 4/3, π, etc) but those are all the same and cancel out

#             # In general, if we're cutting off at a (Dmin, Dmax) but the exponential size distribution wasn't derived with those in mind, will it work out for terminal velocity? might get unrealistic results?
#             # total mass q would be ∼  _λ^-(k+1) * (-Γ(k+1, ∞) + Γ(k+1, 0)) = Γ(k+1) / λ^(k+1)... so if we wanna limit ourselves to Dmin, Dmax, we are applying  the speed from that region to the entire dist...
#             # for single moment... maybe

#             # fall_w += resolve_nan(sum(my_Chen2022_vel_add_sno.(ti, bi, aec, mec, κ, k, ciu, _λ, Dmin, Dmax)))

#             fall_w += sum(my_Chen2022_vel_add_sno.(ti, bi, aec, mec, κ, k, ciu, _λ, Dmin, Dmax; μ=μ)) * mass_weights[i]  # need to weight by mass

#         end

#         if (total_weight = sum(mass_weights)) != 0
#             fall_w /= total_weight # normalize by total mass
#         end

#         fall_w *= χv(prs, precip) # scaling_factor
#         # fall_w = max(FT(0), fall_w)
#     end
#     return resolve_nan(fall_w)
# end



"""
    my_terminal_velocity(param_set, precip, velo_scheme, ρ, q_; ...)

Calculates the mass-weighted average terminal velocity for a snow
particle size distribution (PSD) using the Chen et al. 2022 parameterization.

This function integrates the velocity over the PSD, splitting the integral
into "small" and "large" particle regimes at `D_transition`. It also
applies a physics-based aspect ratio correction and is written to be
numerically stable by scaling the integrals.

---
### Derivation from First Principles
---

**1. The Goal: Mass-Weighted Average Velocity**

We want to find V_avg:
V_avg = ∫ V(R) ⋅ m(R) ⋅ N(R) dR ] / [ ∫ m(R) ⋅ N(R) dR 
The integrals run from R_min to R_max. This function calculates the
total numerator `fall_w` and total denominator `total_weight` by
summing the integrals over different regions.

**2. The Components**

The integral is over R (particle radius in meters).

* **PSD, N(R):** A Gamma distribution:
    N(R) dR ∝ R^μ ⋅ e^(-λR) dR
* **Mass, m(R):** A power law:
    m(R) ∝ R^k (where `k = mec`)
* **Velocity, V(R):** The Chen 2022 formula:
    V(D_mm) ≈ Φ^κ ⋅ Σ aᵢ D_mm^bᵢ e^(-cᵢ D_mm)

**3. Unit Conversion**

The integral is over radius in meters (R_m), but the paper's
coefficients (aᵢ, bᵢ, cᵢ) are for diameter in millimeters (D_mm).

* Conversion: D_mm = 2 ⋅ R_mm = 2 ⋅ (1000 ⋅ R_m) = 2000 ⋅ R_m
* Substituting this into V:
    V(R_m) ≈ Φ(R_m)^κ ⋅ Σ (aᵢ ⋅ 2000^bᵢ) ⋅ R_m^bᵢ ⋅ e^(-(cᵢ ⋅ 2000) R_m)
* The `my_Chen2022_vel_coeffs` function correctly calculates these
    new coefficients for R_m:
    * `aiu = aᵢ ⋅ 2000^bᵢ`
    * `ciu = cᵢ ⋅ 2000`
* This function uses `aiu` and `ciu` directly. No further `* 2`
    or `* 2^bᵢ` conversions are needed.

**4. Aspect Ratio (Φ^κ)**

* The aspect ratio Φ is also a power law of radius: Φ(R) = Φ₀ ⋅ (R/r₀)^α
* The `my_aspect_ratio_coeffs` helpers compute `Φ₀` (as `tmp_phys`),
    `α` (as `α_shape`), and `κ` for each shape.
* Physics:
    * Prolate (cigar): Φ > 1, requires κ = -1/6
    * Oblate (squashed): Φ < 1, requires κ = +1/3
    * Spherical: Φ = 1, requires κ = 0

**5. The Numerator Integral (Deriving `ti` and `ci_pow`)**

We must solve the numerator integral for each region `i` and term `j`:
Numᵢⱼ = ∫ Vᵢⱼ(R) ⋅ m(R) ⋅ N(R) dR
Numᵢⱼ = ∫ [ (Φ₀^κ (R/r₀)^(ακ)) ⋅ (a_{iu,j} R^bⱼ e^(-c_{iu,j} R)) ] ⋅ (R^k) ⋅ (R^μ e^(-λR)) dR

Grouping terms (and absorbing r₀ into constants):
Numᵢⱼ = Cᵢⱼ ⋅ ∫ R^(ακ + bⱼ + k + μ) ⋅ e^(-(λ + c_{iu,j}) R) dR
where Cᵢⱼ = (Φ₀^κ ⋅ a_{iu,j} ⋅ ...)

This integral is solved analytically. To do this robustly, we
calculate the average velocity V_avgᵢⱼ and multiply it by the mass
in the region Mᵢ.
`fall_w = Σ (V_avg_i ⋅ M_i)`
`total_weight = Σ M_i`

* **Mass (Denominator), `mass_weights[i]`:**
    Mᵢ = ∫ m(R)N(R) dR ∝ ∫ R^(k+μ) e^(-λR) dR
    The solution is ∝ λ^-(k+μ+1) ⋅ [Γ(k+μ+1, λR_min) - Γ(k+μ+1, λR_max)]
    The `mass_weights[i]` line calculates exactly this. It is the
    *unscaled* mass in the region, which has a natural λ dependency.

* **Avg. Velocity (Numerator), `sum(my_Chen..._sno.(...))`:**
    This function must calculate V_avgᵢ = Numᵢ / Denomᵢ.
    V_avgᵢ ∝ [ (Φ₀^κ a_{iu,j}) ⋅ γ⁻ᵝ ] / [ λ^-(k+μ+1) ]
    where:
    * γ = λ + c_{iu,j}
    * β = ακ + bⱼ + k + μ + 1
    
    To make this calculation, we pass `ti` as an argument. `ti`
    contains all the prefactors.
    `ti = (Φ₀^κ ⋅ a_{iu,j}) ⋅ (λ + c_{iu,j})⁻ᵝ ⋅ λ^(k+μ+1)`

* **Code Implementation:**
    * `tmp = (tmp_phys)^κ ⋅ _λ^(k + μ + 1)`
        (Aspect ratio ⋅ λ-scaling term)
    * `ci_pow = (ciu .+ _λ) .^ (.-(α_shape .* κ .+ bi .+ k .+ μ .+ 1))`
        (This is (λ + c_{iu,j})⁻ᵝ)
    * `ti = tmp .* aiu .* ci_pow`
        (This is the full prefactor)
"""
function my_terminal_velocity(
    param_set::APS,
    precip::CMT.SnowType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    Nt::FT = FT(NaN),
    D_transition::FT = FT((.625e-3) / 2),
    μ::FT = FT(NaN),
    shape::Val = Val(:Prolate), # Default to Prolate (matches original tmp formula)
) where {FT <: Real}
    fall_w::FT = FT(0)
    if q_ > FT(0)

        prs = TCP.microphysics_params(param_set)

        if isnan(μ)
            μ = μ_from_qN(param_set, precip, q_, Nt; ρ = ρ)
        end

        _r0::FT = r0(prs, precip)
        _λ::FT = lambda(param_set, precip, q_, ρ, Nt, Dmin, Dmax; μ = μ)
        m0c::FT = m0(prs, precip) * χm(prs, precip)
        a0c::FT = a0(prs, precip) * χa(prs, precip)
        mec::FT = me(prs, precip) + Δm(prs, precip)
        aec::FT = ae(prs, precip) + Δa(prs, precip)

        ρ_i::FT = CMP.ρ_cloud_ice(prs)

        # local mass_weights::SA.MVector{2, FT}

        (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l) = my_Chen2022_vel_coeffs(prs, precip, ρ)

        # k = mass weighted
        k = mec

        # Get shape-dependent physics parameters from our new helper
        (; tmp_phys, α_shape, κ) = my_aspect_ratio_coeffs(shape, m0c, mec, a0c, aec, ρ_i, _r0)

        @assert(Dmin <= Dmax, "Dmin must be less than Dmax, got Dmin = $Dmin and Dmax = $Dmax")

        if Dmin < D_transition
            if Dmax <= D_transition
                regions = ((Dmin, Dmax),)
                abcs = ((aiu_s, bi_s, ciu_s),)
            else
                regions = ((Dmin, D_transition), (D_transition, Dmax))
                abcs = ((aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l))
            end
        else
            regions = ((Dmin, Dmax),)
            abcs = ((aiu_l, bi_l, ciu_l),)
        end

        # This is the (Aspect Ratio) * (λ-scaling) term. This is physically consistent and accounts for μ
        tmp = (tmp_phys)^κ * _λ^(k + μ + 1)

        # for (i, ((Dmin, Dmax), (aiu, bi, ciu))) in enumerate(zip(regions, abcs))
        #     # This is the (λ + c_iu)^-β term
        #     # This is physically consistent and accounts for μ
        #     ci_pow = (ciu .+ _λ) .^ (.-(α_shape .* κ .+ bi .+ k .+ μ .+ 1))

        #     # This is the full prefactor for the avg. velocity calculation [[ unit conversions and rolling 2^b into a iu is handled in my_Chen2022_vel_coeffs ]]
        #     ti = tmp .* aiu .* ci_pow

        #     upper_limit = Dmax * _λ
        #     lower_limit = iszero(Dmin) ? FT(0) : Dmin * _λ

        #     # This is the unscaled mass in the region (the denominator of the weighted average for this region)
        #     mass_weights[i] = _λ^-(μ + k + 1) * (-CM1.SF.gamma(μ + k + 1, upper_limit) + CM1.SF.gamma(μ + k + 1, lower_limit))

        #     # V_avg_i = sum(my_Chen..._sno.(...))
        #     # fall_w += V_avg_i * mass_weights[i]
        #     # This computes the total numerator of the integral
        #     fall_w += sum(my_Chen2022_vel_add_sno.(ti, bi, aec, mec, κ, k, ciu, _λ, Dmin, Dmax; μ = μ)) * mass_weights[i]
        # end

        total_weight::FT = FT(0)
        @inbounds for i in eachindex(regions)
            Dmin, Dmax = regions[i]
            aiu, bi, ciu = abcs[i]

            mass_weight = _λ^-(μ + k + 1) * (-CM1.SF.gamma(μ + k + 1, Dmax * _λ) + CM1.SF.gamma(μ + k + 1, Dmin * _λ))
            total_weight += mass_weight
            @inbounds for (aiu_j, bi_j, ciu_j) in zip(aiu, bi, ciu)
                ci_pow = (ciu_j + _λ)^(-(α_shape * κ + bi_j + k + μ + 1))
                ti = tmp * aiu_j * ci_pow
                fall_w += my_Chen2022_vel_add_sno(ti, bi_j, aec, mec, κ, k, ciu_j, _λ, Dmin, Dmax; μ = μ) * mass_weight
            end
        end

        if !iszero(total_weight)
            fall_w /= total_weight # normalize by total mass
        end

        fall_w *= χv(prs, precip) # scaling_factor
    end
    return resolve_nan(fall_w)
end




"""
    my_aspect_ratio_coeffs(shape, m0c, mec, a0c, aec, ρ_i, _r0)

Returns coefficients for the aspect ratio physics based on particle shape.

Arguments:
- `shape`: Val{:Spherical}, Val{:Prolate}, or Val{:Oblate}
- `m0c`, `mec`: Mass-radius power law params (m(R) ∝ m0c * (R/r0)^mec)
- `a0c`, `aec`: Area-radius power law params (a(R) ∝ a0c * (R/r0)^aec)
- `ρ_i`: Particle density
- `_r0`: Reference radius

Returns a NamedTuple:
- `tmp_phys`: The core aspect ratio formula (Φ₀)
- `α_shape`: The exponent α in the power law Φ(R) = Φ₀ * (R/r0)^α
- `κ`: The exponent κ in the velocity formula V ∝ Φ^κ
"""
function my_aspect_ratio_coeffs(
    ::Val{:Spherical},
    m0c::FT,
    mec::FT,
    a0c::FT,
    aec::FT,
    ρ_i::FT,
    _r0::FT,
) where {FT <: Real}
    # For a sphere, Φ = 1 and κ = 0.
    tmp_phys = FT(1)
    α_shape = FT(0)
    κ = FT(0)
    return (; tmp_phys, α_shape, κ)
end

function my_aspect_ratio_coeffs(
    ::Val{:Prolate},
    m0c::FT,
    mec::FT,
    a0c::FT,
    aec::FT,
    ρ_i::FT,
    _r0::FT,
) where {FT <: Real}
    # Prolate (cigar) shape: Φ(R) = 16 * ρᵢ² * a(R)³ / (9 * π * m(R)²)
    # Φ > 1, so κ must be negative.
    # κ = -1/6
    α_shape = 3 * aec - 2 * mec
    tmp_phys = (16 * a0c^3 * ρ_i^2) / (9 * π * m0c^2 * _r0^(α_shape))
    κ = FT(-1 / 6)
    return (; tmp_phys, α_shape, κ)
end

function my_aspect_ratio_coeffs(::Val{:Oblate}, m0c::FT, mec::FT, a0c::FT, aec::FT, ρ_i::FT, _r0::FT) where {FT <: Real}
    # Oblate (squashed) shape: Φ(R) = 3 * sqrt(π) * m(R) / (4 * ρᵢ * a(R)^(3/2))
    # Φ < 1, so κ must be positive.
    # κ = +1/3
    α_shape = mec - FT(1.5) * aec
    tmp_phys = (3 * sqrt(FT(π)) / 4 / ρ_i * m0c / (a0c)^FT(1.5) / _r0^(α_shape))
    κ = FT(1 / 3)
    return (; tmp_phys, α_shape, κ)
end
