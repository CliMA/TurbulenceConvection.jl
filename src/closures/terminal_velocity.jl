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

import SpecialFunctions as SF


Γ_lower(a, x) = SF.gamma(a) - SF.gamma(a, x) # lower incomplete gamma function, ∫_0^x t^(a-1) e^(-t) dt instead of ∫_0^∞
resolve_nan(x::FT; val = 0.0) where {FT} = isnan(x) ? FT(val) : x # replace nan w/ 0


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
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end

function my_Chen2022_vel_coeffs(prs::ACMP, ::CMT.IceType, ρ::FT) where {FT <: Real}

    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    _As, _Bs, _Cs, _Es, _Fs, _Gs = my_Chen2022ice_coeffs(prs, ρ_i)

    ai = (_Es * ρ^_As, _Fs * ρ^_As)
    bi = (_Bs + ρ * _Cs, _Bs + ρ * _Cs)
    ci = (FT(0), _Gs)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end


function my_Chen2022_vel_coeffs(prs::ACMP, ::CMT.SnowType, ρ::FT) where {FT <: Real}

    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    (_As_s, _Bs_s, _Cs_s, _Es_s, _Fs_s, _Gs_s), (_As_l, _Bs_l, _Cs_l, _Es_l, _Fs_l, _Gs_l, _Hs_l) =
        my_Chen2022snow_coeffs(prs, ρ_i)

    # == small ================================================== # Table B2
    ai_s = (_Es_s * ρ^_As_s, _Fs_s * ρ^_As_s)
    bi_s = (_Bs_s + ρ * _Cs_s, _Bs_s + ρ * _Cs_s)
    ci_s = (FT(0), _Gs_s)
    # unit conversions
    aiu_s = ai_s .* 1000 .^ bi_s
    ciu_s = ci_s .* 1000
    # =========================================================== #


    # == large ================================================== # Table B4
    ai_l = (_Bs_l * ρ^_As_l, _Es_l * ρ^_As_l * exp(_Hs_l * ρ))
    bi_l = (_Cs_l, _Fs_l)
    ci_l = (FT(0), _Gs_l)
    # unit conversions
    aiu_l = ai_l .* 1000 .^ bi_l # mm^-bi to m^-bi
    ciu_l = ci_l .* 1000 # /mm to /m
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
function my_Chen2022_vel_add(a::FT, b::FT, c::FT, λ::FT, k::Int; Dmax::FT = Inf, Dmin::FT = 0.0) where {FT <: Real}
    μ = 0 # Exponantial instaed of gamma distribution
    δ = μ + k + 1
    # return a * λ^δ * SF.gamma(b + δ) / (λ + c)^(b + δ) / SF.gamma(δ)
    return a * λ^δ / (λ + c)^(b + δ) * (SF.gamma(b + δ, Dmax * (λ + c)) - SF.gamma(b + δ, Dmin * (λ + c))) /
           (SF.gamma(δ, Dmax * λ) - SF.gamma(δ, Dmin * λ))
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
    k::Int,
    c::FT,
    λ::FT,
    Dmin::FT,
    Dmax::FT,
) where {FT}

    return t * (
        SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1, Dmax * (λ + c)) -
        SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1, Dmin * (λ + c))
    ) / (SF.gamma(k + 1, Dmax * λ) - SF.gamma(k + 1, Dmin * λ))


    # Checking here
    # return my_Chen2022_vel_add(t, b, c, λ, k; Dmax=Dmax, Dmin=Dmin) # treat same as small ice one

end
# This is for \int_0^∞, but we want transitions
# my_Chen2022_vel_add_sno(t, b, aec, mec, κ, k; Dmax=Dmax) =
# t * SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1) /
# SF.gamma(k + 1)

# ============================================================================= #

"""
    n0(prs, q_sno, ρ, snow_type)

 - `prs` - abstract set with Earth parameters
 - `q_sno` -  snow specific humidity
 - `ρ` - air density
 - `type` - type for dispatch

Returns the intercept parameter of the assumed Marshall-Palmer distribution
"""
function n0(prs::ACMP, q_sno::FT, ρ::FT, ::CMT.SnowType) where {FT <: Real}

    _ν_sno::FT = CMP.ν_sno(prs)
    _μ_sno::FT = CMP.μ_sno(prs)

    # TODO               this max should be replaced by
    #                    limiting inside a PhasePartition struct for
    #                    precipitation (once it is implemented)
    return _μ_sno * (ρ * max(0, q_sno))^_ν_sno
end
n0(prs::ACMP, ::Any, ::Any, ::CMT.IceType) = CMP.n0_ice(prs)
n0(prs::ACMP, ::Any, ::Any, ::CMT.RainType) = CMP.n0_rai(prs)


"""
if you think you know what the total number concentration N_t should be, you can use this
For marshall palmer, we take gamma dist N_0 D^μ exp(-λD) and set μ=0, leaving N_0 exp(-λD)

the total mass is based on a bunch of constants below... and n_0 is just set.
However, integrating we can show ∫(PSD) = Γ(μ+1) / λ^(μ+1) N_0 = N_t  so N_0 = N_t λ^(μ+1) / Γ(μ+1)

For μ = 0 this is just N_0 =  N_t λ

If we take λ as written in https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/#Assumed-particle-size-relationships Eq 7 (which accounts for the size <-> mass relationship assumptions)
then we can take N_t = N_0 Γ(μ+1) / λ^(μ+1) = N_0 / λ, plug in our desired N_t, and solve for N_0

I'm not sure if we can use other μ, bc I am not certain on how λ was derived, but i assume it assumes μ=0
"""
function n0(
    prs::ACMP,
    q::FT,
    ρ::FT,
    Nt::FT,
    precip::Union{CMT.IceType, CMT.RainType, CMT.SnowType}, # not sure if have a liqtype
) where {FT <: Real}
    # _n0::FT = n0(prs, q, ρ, precip)
    _r0::FT = r0(prs, precip)
    _m0::FT = m0(prs, precip)
    _me::FT = me(prs, precip)
    _Δm::FT = Δm(prs, precip)
    _χm::FT = χm(prs, precip)

    n0::FT = FT(0)
    if q > FT(0)
        E::FT = FT(1 / (_me + _Δm + 1))
        λ_no_n0::FT = (_χm * _m0 * SF.gamma(_me + _Δm + FT(1)) / ρ / q / _r0^(_me + _Δm))^E # We've divided λ by n_0 ^ FT(1 / (_me + _Δm + 1)) = n_0^E so that λ = λ_no_n0 * n_0^E

        # call FT(1 / (_me + _Δm + 1)) E
        #  So we have  N_t = n_0 / λ = n_0 / (λ_no_n0 * n_0^E) = n_0^(1-E) / λ_no_n0
        # Then, n_0^(1-E) = N_t λ_no_n0
        # n_0 = (N_t λ_no_n0)^(1 / (1-E))

        n0 = (Nt * λ_no_n0)^(1 / (1 - E))

        return n0
    end

    return n0
end


"""
    v0(prs, ρ, rain_type)

 - `prs` - abstract set with Earth parameters
 - `ρ` - air density
 - `type` - type for dispatch

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
function v0(prs::ACMP, ρ::FT, ::CMT.RainType) where {FT <: Real}

    _ρ_cloud_liq::FT = CMP.ρ_cloud_liq(prs)
    _C_drag::FT = CMP.C_drag(prs)
    _grav::FT = CMP.grav(prs)
    _r0_rai::FT = CMP.r0_rai(prs)

    return sqrt(FT(8 / 3) / _C_drag * (_ρ_cloud_liq / ρ - FT(1)) * _grav * _r0_rai)
end
v0(prs::ACMP, ::Any, ::CMT.SnowType) = CMP.v0_sno(prs)

# Other ice/rain/snow parameters to dispatch over
a_vent(prs::ACMP, ::CMT.RainType) = CMP.a_vent_rai(prs)
b_vent(prs::ACMP, ::CMT.RainType) = CMP.b_vent_rai(prs)
a_vent(prs::ACMP, ::CMT.SnowType) = CMP.a_vent_sno(prs)
b_vent(prs::ACMP, ::CMT.SnowType) = CMP.b_vent_sno(prs)

r0(prs::ACMP, ::CMT.IceType) = CMP.r0_ice(prs)
m0(prs::ACMP, ::CMT.IceType) = CMP.m0_ice(prs)
me(prs::ACMP, ::CMT.IceType) = CMP.me_ice(prs)
χm(prs::ACMP, ::CMT.IceType) = CMP.χm_ice(prs)
Δm(prs::ACMP, ::CMT.IceType) = CMP.Δm_ice(prs)

r0(prs::ACMP, ::CMT.RainType) = CMP.r0_rai(prs)
m0(prs::ACMP, ::CMT.RainType) = CMP.m0_rai(prs)
me(prs::ACMP, ::CMT.RainType) = CMP.me_rai(prs)
a0(prs::ACMP, ::CMT.RainType) = CMP.a0_rai(prs)
ae(prs::ACMP, ::CMT.RainType) = CMP.ae_rai(prs)
ve(prs::ACMP, ::CMT.RainType) = CMP.ve_rai(prs)
χm(prs::ACMP, ::CMT.RainType) = CMP.χm_rai(prs)
Δm(prs::ACMP, ::CMT.RainType) = CMP.Δm_rai(prs)
χa(prs::ACMP, ::CMT.RainType) = CMP.χa_rai(prs)
Δa(prs::ACMP, ::CMT.RainType) = CMP.Δa_rai(prs)
χv(prs::ACMP, ::CMT.RainType) = CMP.χv_rai(prs)
Δv(prs::ACMP, ::CMT.RainType) = CMP.Δv_rai(prs)

r0(prs::ACMP, ::CMT.SnowType) = CMP.r0_sno(prs)
m0(prs::ACMP, ::CMT.SnowType) = CMP.m0_sno(prs)
me(prs::ACMP, ::CMT.SnowType) = CMP.me_sno(prs)
a0(prs::ACMP, ::CMT.SnowType) = CMP.a0_sno(prs)
ae(prs::ACMP, ::CMT.SnowType) = CMP.ae_sno(prs)
ve(prs::ACMP, ::CMT.SnowType) = CMP.ve_sno(prs)
χm(prs::ACMP, ::CMT.SnowType) = CMP.χm_sno(prs)
Δm(prs::ACMP, ::CMT.SnowType) = CMP.Δm_sno(prs)
χa(prs::ACMP, ::CMT.SnowType) = CMP.χa_sno(prs)
Δa(prs::ACMP, ::CMT.SnowType) = CMP.Δa_sno(prs)
χv(prs::ACMP, ::CMT.SnowType) = CMP.χv_sno(prs)
Δv(prs::ACMP, ::CMT.SnowType) = CMP.Δv_sno(prs)

"""
    lambda(q, ρ, n0, m0, me, r0, χm, Δm)

 - `prs` - set with free parameters
 - `precip` - a type for cloud ice, rain or snow
 - `q` - specific humidity of rain, ice or snow
 - `ρ` - air density

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(
    prs::ACMP,
    precip::Union{CMT.IceType, CMT.RainType, CMT.SnowType},
    q::FT,
    ρ::FT,
    # Nt::Union{FT, Nothing},
    Nt::FT, # testing for type stability, use NaN instead of nothing
    Dmin::FT,
    Dmax::FT,
) where {FT <: Real}

    # _n0::FT = isnothing(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, Nt, precip) # use wanted Nt If given
    _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, Nt, precip) # use wanted Nt If given
    _r0::FT = r0(prs, precip)
    _m0::FT = m0(prs, precip)
    _me::FT = me(prs, precip)
    _Δm::FT = Δm(prs, precip)
    _χm::FT = χm(prs, precip)


    λ::FT = FT(0)

    if q > FT(0)
        λ = (_χm * _m0 * _n0 * SF.gamma(_me + _Δm + FT(1)) / ρ / q / _r0^(_me + _Δm))^FT(1 / (_me + _Δm + 1))
    end

    # this would be scaling N to get the right q, but no point bc doesnt affect terminal velocity anyway (we wouldn't recompute λ, and the mass weighting wouldn't change...)
    # if Dmin > 0 || Dmax < Inf
    # _n0 = SF.Gamma(μ+1) / (SF.Gamma(μ+1, D_max * λ) - SF.Gamma(μ+1, D_min * λ))
    # end

    return λ
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
    prs::ACMP,
    precip::CMT.AbstractPrecipType,
    velo_scheme::CMT.Blk1MVelType,
    ρ::FT,
    q_::FT;
    Dmin::FT = 0.0, # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = Inf, # not implemented yet
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = NaN, # testing for type stability, use NaN instead of nothing
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        _r0::FT = r0(prs, precip)
        _me::FT = me(prs, precip)
        _Δm::FT = Δm(prs, precip)
        _χm::FT = χm(prs, precip)
        _χv::FT = χv(prs, precip)
        _v0::FT = v0(prs, ρ, precip)
        _ve::FT = ve(prs, precip)
        _Δv::FT = Δv(prs, precip)
        _λ::FT = lambda(prs, precip, q_, ρ, Nt, Dmin, Dmax)

        fall_w =
            _χv * _v0 * (_λ * _r0)^(-_ve - _Δv) * SF.gamma(_me + _ve + _Δm + _Δv + FT(1)) / SF.gamma(_me + _Δm + FT(1))
    end

    return resolve_nan(fall_w)
end
function my_terminal_velocity(
    prs::ACMP,
    precip::CMT.RainType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = 0.0, # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = Inf,
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = NaN, # testing for type stability, use NaN instead of nothing
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = my_Chen2022_vel_coeffs(prs, precip, ρ)
        # size distribution parameter
        _λ::FT = lambda(prs, precip, q_, ρ, Nt, Dmin, Dmax)

        # eq 20 from Chen et al 2022
        fall_w = sum(my_Chen2022_vel_add.(aiu, bi, ciu, _λ, 3; Dmin = Dmin, Dmax = Dmax))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return resolve_nan(fall_w)
end


function my_terminal_velocity(
    prs::ACMP,
    precip::CMT.IceType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = 0.0, # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = Inf,
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = NaN, # testing for type stability, use NaN instead of nothing
    D_transition::FT = 0.625e-3, # .625mm is the transition from small to large ice crystals in Chen paper
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)


        ρ_i::FT = CMP.ρ_cloud_ice(prs)

        # Ok so here's the thing.... we want q to fit between Dmin and Dmax...
        # so we need a lambda that works for that....

        _λ::FT = lambda(prs, precip, q_, ρ, Nt, Dmin, Dmax)

        # ================================================================= #

        # # coefficients from Appendix B from Chen et. al. 2022
        # aiu, bi, ciu = my_Chen2022_vel_coeffs(prs, precip, ρ)

        # # eq 20 from Chen et al 2022
        # fall_w = sum(my_Chen2022_vel_add.(aiu, bi, ciu, _λ, 3; Dmin=Dmin, Dmax=Dmax))
        # # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        # fall_w = max(FT(0), fall_w)

        # ================================================================= #
        k = 3 # mass weighted
        # coefficients from Appendix B from Chen et. al. 2022
        (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l) = my_Chen2022_vel_coeffs(prs, CM.CommonTypes.SnowType(), ρ)

        local mass_weights::SA.MVector{2, FT}
        mass_weights = SA.@MVector [FT(0), FT(0)]
        if Dmin < D_transition
            if Dmax <= D_transition
                regions = [(Dmin, Dmax)]
                abcs = [(aiu_s, bi_s, ciu_s)]
            else
                regions = [(Dmin, D_transition), (D_transition, Dmax)]
                abcs = [(aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l)]
            end
        else
            regions = [(Dmin, Dmax)]
            abcs = [(aiu_l, bi_l, ciu_l)]
        end

        fall_w = FT(0)
        for (i, ((Dmin, Dmax), (aiu, bi, ciu))) in enumerate(zip(regions, abcs)) # we basically need to sum the integral as before but over all regions

            mass_weights[i] = _λ^-(k + 1) * (-SF.gamma(k + 1, Dmax * _λ) + SF.gamma(k + 1, Dmin * _λ)) # missing constants from the integral (n_0, 4/3, π, etc) but those are all the same and cancel out


            # eq 20 from Chen et al 2022
            fall_w += sum(my_Chen2022_vel_add.(aiu, bi, ciu, _λ, 3; Dmin = Dmin, Dmax = Dmax)) .* mass_weights[i]
        end

        if (total_weight = sum(mass_weights)) != 0
            fall_w /= total_weight # normalize by total mass
        end


        # ================================================================= #

    end
    return resolve_nan(fall_w)
end




"""
Adjusted for snow...

We go to the from Dmin to Dmax, with a break at the transition point if necessary

Is this smooth?
Or are we supposed to integrate both the whole way and the gamma functions just decay as necessary? idk...

"""
function my_terminal_velocity(
    prs::ACMP,
    precip::CMT.SnowType,
    velo_scheme::CMT.Chen2022Type,
    ρ::FT,
    q_::FT;
    Dmin::FT = 0.0, # Default to 0, but in CloudMicrophysics.jl technically it's 62.5μm...
    Dmax::FT = Inf, # not really needed for snow? maybe a Dmin? but that wouldn't make a huge diff i dont think idk...
    # Nt::Union{FT, Nothing} = nothing,
    Nt::FT = NaN, # testing for type stability, use NaN instead of nothing
    D_transition::FT = Inf, # .625mm is the transition from small to large ice crystals in Chen paper
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)


        # short circuit cause greater than .625 mm is giving negative values for some reason... so we'll just ignore that I guess...
        # return my_terminal_velocity(prs, ice_type, velo_scheme, ρ, q_; Dmin=Dmin, Dmax=Dmax, Nt=Nt) # backup for now bc large ice is broken...

        _r0::FT = r0(prs, precip)
        _λ::FT = lambda(prs, precip, q_, ρ, Nt, Dmin, Dmax)
        m0c::FT = m0(prs, precip) * χm(prs, precip)
        a0c::FT = a0(prs, precip) * χa(prs, precip)
        mec::FT = me(prs, precip) + Δm(prs, precip)
        aec::FT = ae(prs, precip) + Δa(prs, precip)

        ρ_i::FT = CMP.ρ_cloud_ice(prs)

        # D_transition::FT = 0.625e-3 # .625mm is the transition from small to large ice crystals in Chen paper

        local mass_weights::SA.MVector{2, FT}

        # coefficients from Appendix B from Chen et. al. 2022
        (aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l) = my_Chen2022_vel_coeffs(prs, precip, ρ)

        κ = FT(-1 / 3) #oblate
        # κ = FT(1 / 3) #oblate (anna says this is right)

        k = 3 # mass weighted

        @assert(Dmin <= Dmax, "Dmin must be less than Dmax, got Dmin = $Dmin and Dmax = $Dmax")

        # Is this smooth?
        # or are we supposed to integrate both the whole way and the gamma functions just decay as necessary? idk...
        mass_weights = SA.@MVector [FT(0), FT(0)]
        if Dmin < D_transition
            if Dmax <= D_transition
                regions = [(Dmin, Dmax)]
                abcs = [(aiu_s, bi_s, ciu_s)]
            else
                regions = [(Dmin, D_transition), (D_transition, Dmax)]
                abcs = [(aiu_s, bi_s, ciu_s), (aiu_l, bi_l, ciu_l)]
            end
        else
            regions = [(Dmin, Dmax)]
            abcs = [(aiu_l, bi_l, ciu_l)]
        end

        tmp = _λ^(k + 1) * ((16 * a0c^3 * ρ_i^2) / (9 * π * m0c^2 * _r0^(3 * aec - 2 * mec)))^κ
        fall_w = FT(0)
        for (i, ((Dmin, Dmax), (aiu, bi, ciu))) in enumerate(zip(regions, abcs)) # we basically need to sum the integral as before but over all regions
            ci_pow = (2 .* ciu .+ _λ) .^ (.-(3 .* aec .* κ .- 2 .* mec .* κ .+ bi .+ k .+ 1))

            ti = tmp .* aiu .* FT(2) .^ bi .* ci_pow

            # k = 3 for mass, μ = 0 for exponential size distribution instead of gamma
            mass_weights[i] = _λ^-(k + 1) * (-SF.gamma(k + 1, Dmax * _λ) + SF.gamma(k + 1, Dmin * _λ)) # missing constants from the integral (n_0, 4/3, π, etc) but those are all the same and cancel out
            # In general, if we're cutting off at a (Dmin, Dmax) but the exponential size distribution wasn't derived with those in mind, will it work out for terminal velocity? might get unrealistic results?
            # total mass q would be ∼  _λ^-(k+1) * (-Γ(k+1, ∞) + Γ(k+1, 0)) = Γ(k+1) / λ^(k+1)... so if we wanna limit ourselves to Dmin, Dmax, we are applying  the speed from that region to the entire dist...
            # for single moment... maybe

            # fall_w += resolve_nan(sum(my_Chen2022_vel_add_sno.(ti, bi, aec, mec, κ, k, ciu, _λ, Dmin, Dmax)))

            fall_w += sum(my_Chen2022_vel_add_sno.(ti, bi, aec, mec, κ, k, ciu, _λ, Dmin, Dmax)) * mass_weights[i]  # need to weight by mass

        end

        if (total_weight = sum(mass_weights)) != 0
            fall_w /= total_weight # normalize by total mass
        end

        # fall_w = max(FT(0), fall_w)
    end
    return resolve_nan(fall_w)
end
