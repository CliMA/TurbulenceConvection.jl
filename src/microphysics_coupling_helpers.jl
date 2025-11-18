#=
Useful microphysics definitions
=#



Γ_lower(a, x) = CM1.SF.gamma(a) - CM1.SF.gamma(a, x) # lower incomplete gamma function, ∫_0^x t^(a-1) e^(-t) dt instead of ∫_0^∞
# Γ_lower(a, x) = CM1.SF.gamma_inc(a, x)[1] * CM1.SF.gamma(a) # lower incomplete gamma function, ∫_0^x t^(a-1) e^(-t) dt instead of ∫_0^∞

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

function n0(prs::ACMP, q_liq::FT, ρ::FT, ::CMT.LiquidType) where {FT <: Real}
    _μ_liq::FT = FT(5)
    N = FT(200) * 1e6
    # average 5 micron droplets
    _λ_liq::FT = (_μ_liq + 1) / (5e-6)
   # we want total N to be something like 200 / cm^3
    n0 = N * _λ_liq^(_μ_liq + 1) / CM1.SF.gamma(_μ_liq + 1)
    return n0
end


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
    precip::CMTWaterTypes,
    Nt::FT,
    μ::FT # μ is the exponent in the PSD, default to 0 for marshall palmer
    ;
    Dmin::FT = FT(0), # Dmin is the lower limit of the PSD, default to 0
    Dmax::FT = FT(Inf), # Dmax is the upper limit of the PSD, default to Inf    
) where {FT <: Real}
    # _n0::FT = n0(prs, q, ρ, precip)
    _r0::FT = r0(prs, precip)
    _m0::FT = m0(prs, precip)
    _me::FT = me(prs, precip)
    _Δm::FT = Δm(prs, precip)
    _χm::FT = χm(prs, precip)

    _n0::FT = FT(0)
    if q > FT(0)

        if iszero(Dmin) && isinf(Dmax)

            if isnan(μ)
                error("Got μ = NaN, please provide an override value")
            end

            if iszero(μ)
                E = FT(1 / (_me + _Δm + 1))
                λ_no_n0 = (_χm * _m0 * CM1.SF.gamma(_me + _Δm + FT(1)) / ρ / q / _r0^(_me + _Δm))^E # We've divided λ by n_0 ^ FT(1 / (_me + _Δm + 1)) = n_0^E so that λ = λ_no_n0 * n_0^E

                # Call FT(1 / (_me + _Δm + 1)) `E`
                #  So we have  N_t = n_0 / λ = n_0 / (λ_no_n0 * n_0^E) = n_0^(1-E) / λ_no_n0
                # Then, n_0^(1-E) = N_t λ_no_n0, and :
                _n0 = (Nt * λ_no_n0)^(1 / (1 - E))
            else
                E = FT(1 / (_me + _Δm + μ + 1)) # this is the exponent in the PSD, so we can use it to scale n_0
                λ_no_n0 = (_χm * _m0 * CM1.SF.gamma(_me + _Δm + μ + FT(1)) / ρ / q / _r0^(_me + _Δm))^E # We've divided λ by n_0 ^ E so that λ = λ_no_n0 * n_0^E
                # w/ loggama you could do something like    log_λ_no_n0 = E * (log(_χm * _m0) + CM1.SF.loggamma(_me + _Δm + μ + FT(1)) - log(ρ) - log(q) - (_me + _Δm) * log(_r0))

                # Call FT(1 / (_me + _Δm + μ + 1)) `E`
                #  So we have  N_t = n_0 *  Γ(μ+1) / λ^(μ+1) = n_0 * Γ(μ+1) / (λ_no_n0 * n_0^E)^(μ+1) = n_0^(1-E) * Γ(μ+1) / λ_no_n0^(μ+1)
                # Then, n_0^(1-E) = N_t λ_no_n0^(μ+1) / Γ(μ+1), and:
                _n0 = (Nt * λ_no_n0^(μ + FT(1)) / CM1.SF.gamma(μ + FT(1)))^(1 / (1 - E*(μ + FT(1))))
            end
        else
            if !iszero(Dmin) && !isinf(Dmax) # use numerical solver...
                _n0, _ = get_n0_lambda(prs, precip, q, ρ, Nt, μ; Dmin=Dmin, Dmax=Dmax)
            end
        end
    end

    if !isfinite(_n0) || (iszero(_n0) && !iszero(Nt) && !iszero(q)) # if we wanted a nonzero Nt and q, but got n0 = 0 or Inf, something is wrong
        error("Got n0 = $_n0; q = $q; ρ = $ρ; Nt = $Nt; μ = $μ; Dmin = $Dmin; Dmax = $Dmax; precip = $precip;")
    end

    if isinf(_n0) # can happen if q becomes very small for example but you have a diagnosed N from INP that is large...
        _n0 = floatmax(FT) # just set to a large number (still could lead to problems down the road I guess? idk, testing w/o for a second)
    end

    return _n0
end

function n0(param_set::APS, q::FT, ρ::FT, precip::CMTWaterTypes, Nt::FT; μ::FT = FT(NaN), Dmin::FT = FT(0), Dmax::FT = FT(Inf)) where {FT <: Real}

    microphys_params = TCP.microphysics_params(param_set)

    if !isnan(Nt) # aerosol inside each particle...
        r_min::FT = param_set.user_params.particle_min_radius # this is the minimum radius we want to consider, if r < r_min, we return 0
        χm::FT = get_χm(param_set, precip) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
        q_r_min = particle_mass(microphys_params, precip, r_min, χm)
        q += Nt * q_r_min
    end
    if isnan(μ)
        μ = μ_from_qN(param_set, precip, q, Nt; ρ=ρ)
    end

    return n0(microphys_params, q, ρ, precip, Nt, μ; Dmin=Dmin, Dmax=Dmax)
end

# We know N = n0 Γ(μ + 1) / λ^(μ+1), so we just invert this
function n0_from_Nλ(N::FT, λ::FT; μ::FT = FT(0)) where{FT <: Real}
    return N * λ^(μ+1) / CM1.SF.gamma(μ + 1)
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
v0(prs::ACMP, ::FT, ::CMT.SnowType) where{FT <: Real} = CMP.v0_sno(prs)

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

# my ice additions ----------------------------------------------- #
# For these, you might want to fall back to snow but they need to work with r0 for ice etc.... just assume regular
a0(prs::ACMP, ice_type::CMT.IceType) = eltype(prs)(π) * r0(prs, ice_type)^2 # [ should just be π * r0^2 ]
v0(prs::ACMP, ρ::FT, ::CMT.IceType) where {FT <: Real} = error("`v0` not implemented for ice yet. Try using Chen parameterization") # (we will be using Chen2022 mostly anyway)

χa(prs::ACMP, ::CMT.IceType) = eltype(prs)(1) 
ae(prs::ACMP, ::CMT.IceType) = eltype(prs)(2)
Δa(prs::ACMP, ::CMT.IceType) = eltype(prs)(0) # [ should just be 0 ]

χv(prs::ACMP, ::CMT.IceType) = eltype(prs)(1) #
ve(prs::ACMP, ::CMT.IceType) = error("`ve` not implemented for ice yet. Try using Chen parameterization") # (we will be using Chen2022 mostly anyway)
Δv(prs::ACMP, ::CMT.IceType) = eltype(prs)(0) #

# my liq additions ------------------------------------------------------------ #
r0(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(10e-6) # 10 micron mean.
function m0(prs::ACMP, liq_type::CMT.LiquidType)
    FT = eltype(prs)
    return FT(4/3) * FT(π) * CMP.ρ_cloud_liq(prs) * r0(prs, liq_type)^3
end
χm(prs::ACMP, ::CMT.LiquidType) = error("I put this in param_set.user_params, use that one...")
me(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(3) # this is the exponent in the PSD, so we can use it to scale n_0
Δm(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(0) #
χa(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(1) #
a0(prs::ACMP, liq_type::CMT.LiquidType) = eltype(prs)(π) * r0(prs, liq_type)^2
ae(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(2) # this is the exponent in the PSD, so we can use it to scale n_0
Δa(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(0) 
χv(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(1) #
ve(prs::ACMP, ::CMT.LiquidType) = error("`ve` not implemented for liquid yet. Try using Chen parameterization") # (we will be using Chen2022 mostly anyway)
Δv(prs::ACMP, ::CMT.LiquidType) = eltype(prs)(0) #



# ------------------------------------------------------------ #


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


Derivation:


"""
function lambda(
    prs::ACMP,
    precip::CMTWaterTypes, # not sure if have a liqtype
    q::FT,
    ρ::FT,
    # Nt::Union{FT, Nothing},
    Nt::FT, # testing for type stability, use NaN instead of nothing
    μ::FT,
    ;
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    _n0::FT = FT(NaN), # allow avoiding repeated caclulations of _n0, since a lot of fcns use both n0 and lambda,
    _χm::FT = FT(NaN),
) where {FT <: Real}



    λ::FT = FT(0)
    if q > FT(0)

        if isnan(μ)
            μ = FT(0) # default to marshall palmer
        end

        if !iszero(Dmin) || !isinf(Dmax) # use numerical solver...
            _n0, λ = get_n0_lambda(prs, precip, q, ρ, Nt, μ; Dmin=Dmin, Dmax=Dmax)
        else
            # _n0::FT = isnothing(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, Nt, precip) # use wanted Nt If given
            if isnan(_n0)
                _n0::FT = isnan(Nt) ? n0(prs, q, ρ, precip) : n0(prs, q, ρ, precip, Nt, μ; Dmin=Dmin, Dmax=Dmax) # use wanted Nt If given
            end





            _r0::FT = r0(prs, precip)
            _m0::FT = m0(prs, precip)
            _me::FT = me(prs, precip)
            _Δm::FT = Δm(prs, precip)
            # _χm::FT = χm(prs, precip) # technically we change this in user_params... 
            _χm::FT = isnan(_χm) ? χm(prs, precip) : _χm

            λ = (_χm * _m0 * _n0 * CM1.SF.gamma(_me + _Δm + μ + FT(1)) / ρ / q / _r0^(_me + _Δm))^FT(1 / (_me + _Δm + μ + 1))
            # w/ loggama you could do something like     λ = exp((1 / (_me + _Δm + μ + 1)) * (log(_χm * _m0 * _n0) + CM1.SF.loggamma(_me + _Δm + μ + FT(1)) - log(ρ) - log(q) - (_me + _Δm) * log(_r0)))

            if !isfinite(_n0) || (iszero(_n0) && !iszero(Nt) && !iszero(q)) # if we wanted a nonzero Nt and q, but got n0 = 0 or Inf, something is wrong
                @error("Got n0 = $_n0; λ = $λ; q = $q; ρ = $ρ; Nt = $Nt; μ = $μ; Dmin = $Dmin; Dmax = $Dmax; precip = $precip; _r0 = $_r0; _m0 = $_m0; _me = $_me; _Δm = $_Δm; _χm = $_χm;")
            end

        end
    end

    # this would be scaling N to get the right q, but no point bc doesnt affect terminal velocity anyway (we wouldn't recompute λ, and the mass weighting wouldn't change...)
    # if Dmin > 0 || Dmax < Inf
    # upper_limit = Dmax * λ # Dmax should not be zero...
    # lower_limit = iszero(Dmin) ? FT(0) : Dmin * λ # avoid 0*Inf
    # _n0 = CM1.SF.Gamma(μ+1) / (CM1.SF.Gamma(μ+1, upper_limit) - CM1.SF.Gamma(μ+1, lower_limit))
    # end

    if isinf(λ) # can happen if q becomes very small for example but you have a diagnosed N from INP that is large...
        λ = floatmax(FT) # just set to a large number (still could lead to problems down the road I guess? idk, testing w/o for a second)
    end

    return λ
end

lambda(prs::ACMP, precip::CMTWaterTypes, q::FT, ρ::FT, Nt::FT, μ::FT, Dmin::FT, Dmax::FT; _n0::FT = FT(NaN), _χm::FT = FT(NaN)) where {FT <: Real} = lambda(prs, precip, q, ρ, Nt, μ; Dmin=Dmin, Dmax=Dmax, _n0=_n0, _χm=_χm)

function lambda(param_set::APS, precip::CMTWaterTypes, q::FT, ρ::FT, Nt::FT; Dmin::FT = FT(0), Dmax::FT = FT(Inf), _n0::FT = FT(NaN), μ::FT = FT(NaN), _χm::FT = FT(NaN)) where {FT <: Real}
    microphys_params = TCP.microphysics_params(param_set)

    if !isnan(Nt) # aerosol inside each particle...
        r_min::FT = param_set.user_params.particle_min_radius # this is the minimum radius we want to consider, if r < r_min, we return 0
        # _χm::FT = get_χm(param_set, precip) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
        _χm::FT = isnan(_χm) ? get_χm(param_set, precip) : _χm
        q_r_min = particle_mass(microphys_params, precip, r_min, _χm)
        q += Nt*q_r_min
    end
    if isnan(μ)
        μ = μ_from_qN(param_set, precip, q, Nt; ρ=ρ)
    end

    _χm = isnan(_χm) ? get_χm(param_set, precip) : _χm

    return lambda(microphys_params, precip, q, ρ, Nt, μ; Dmin=Dmin, Dmax=Dmax, _n0=_n0, _χm=_χm)
end

lambda(param_set::APS, precip::CMTWaterTypes, q::FT, ρ::FT, Nt::FT, Dmin::FT, Dmax::FT; _n0::FT = FT(NaN), μ::FT = FT(NaN), _χm::FT = FT(NaN)) where {FT <: Real} = lambda(param_set, precip, q, ρ, Nt; μ=μ, Dmin=Dmin, Dmax=Dmax, _n0=_n0, _χm=_χm)















"""
    solve_lambda_truncated_gamma(q_cloud, N_cloud, D_th; μ=2.0, d=3.0, ρ=1.0, χ_m=π*1000/6)

Solve for the shape parameter `λ` and normalization `n₀` of a truncated gamma distribution

    n(D) = n₀ D^μ exp(-λ D),   D ∈ [0, D_th]

such that it reproduces the given number concentration `N_cloud` and mass mixing ratio `q_cloud`.

# Arguments
- `q_cloud`: total mass mixing ratio [kg/kg]
- `N_cloud`: number concentration [1/m³]
- `D_th`: truncation diameter (upper bound of integration) [m]

# Keyword Arguments
- `μ`: shape parameter of the gamma distribution (default: 2.0)
- `d`: power of D in the mass-size relation m(D) ∝ D^d (default: 3.0)
- `ρ`: air density [kg/m³] (default: 1.0)
- `χ_m`: geometric factor in m(D) = χ_m D^d (default: π*1000/6 for spherical water)

# Returns
Tuple `(λ, n₀)` or `nothing` if no valid solution exists.

If `D_th` is `Inf`, uses the closed-form full gamma solution.
"""

function get_n0_lambda(
    prs::ACMP, # won't work for liquid
    q_type::CMTWaterTypes,
    q::FT,
    ρ::FT,
    N::FT,
    μ::FT,
    ;
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    add_dry_aerosol_mass::Bool=false, # if true, we add the dry aerosol mass to q before solving
    particle_min_radius::FT = FT(0.2e-6)
) where {FT <: Real}


    if isnan(μ)
        error("μ is NaN, please provide an override value!")
    end

    if add_dry_aerosol_mass
        if !(q_type isa CMT.LiquidType)
            q += mass(prs, q_type, particle_min_radius, N, get_χm(prs, q_type); monodisperse=true)
        else
            q += mass(prs, q_type, particle_min_radius, N, FT(1); monodisperse=true) # liquid needs param_set to get χm
        end
    end
    
    if iszero(Dmin) && !isinf(Dmax)
        
        # return solve_n0_lambda_upper_truncated_gamma(prs, q, q_type, N, Dmax, μ, _me + _Δm, ρ, _χm * _m0 * _r0^(-(_me + _Δm)))
        # This is just too unreliable right now, it has essentially a 0 percent success rate...
        # The real problem is that as soon as cbrt(q/c_mass) hits about 5 times Dmax, the solutions vanish rapidly.
        # So even if we could solve it, we'd rapidly face problems with N and λ not being stable across cells/timesteps...
        # This is probably why they did that `threshold` style solve.
        # we really should take our `w` from the full distribution I think then, as depending on how badly low `N` is, the part under 62.5 could be totally fraudulent...
        # It sems even now that our N are reasonable, 62.5 just doesn't cut it

        # _n0 = n0(prs, q, ρ, q_type, N, μ; Dmin=Dmin, Dmax=Dmax) # this is the n0 we want to use
        _n0 = isnan(N) ? n0(prs, q, ρ, q_type) : n0(prs, q, ρ, q_type, N, μ; Dmin=Dmin, Dmax=Dmax)
        λ = lambda(prs, q_type, q, ρ, N, μ; _n0=_n0, Dmin=FT(0), Dmax=FT(Inf)) # this is the λ we want to use

        # we could check feasibility
        # is_feasible, q_max = check_N_q_μ_feasibility(q, N, μ, _χm * _m0 * _r0^(-(_me + _Δm)), _me + _Δm; Dmin=FT(0), Dmax=Dmax)
        # if not feasible, find the smalles Dmax that is feasible
        # if !is_feasible
        #     Dmax = minimum_rmax_for_solution(q, N, μ, _χm * _m0 * _r0^(-(_me + _Δm)), _me + _Δm)
        #     return ...
        # also `truncated_gamma_solver()` might be easier to use or more correct than `solve_n0_lambda_upper_truncated_gamma()`? not sure yet...
        # end
        return _n0, λ
    # elseif !iszero(Dmin) && isinf(Dmax)
    # elseif !iszero(Dmin) && !isinf(Dmax)
    else

        # @warn "Dmin = $Dmin; Dmax = $Dmax; q = $q; N = $N; μ = $μ; ρ = $ρ"

        # Normal Closed form solution...
        # _n0 = n0(prs, q, ρ, q_type, N, μ)
        _n0 = isnan(N) ? n0(prs, q, ρ, q_type) : n0(prs, q, ρ, q_type, N, μ; Dmin=Dmin, Dmax=Dmax)
        λ = lambda(prs, q_type, q, ρ, N, μ; _n0=_n0, Dmin=FT(0), Dmax=Dmax)
        return _n0, λ
    end
end

function get_n0_lambda(param_set::APS, q_type::CMTWaterTypes, q::FT, ρ::FT, N::FT; μ::FT = FT(NaN), Dmin::FT = FT(0), Dmax::FT = FT(Inf), add_dry_aerosol_mass::Bool=false) where {FT <: Real}
    if isnan(μ)
        μ = μ_from_qN(param_set, q_type, q, N; ρ=ρ)
    end
    microphys_params = TCP.microphysics_params(param_set)
    return get_n0_lambda(microphys_params, q_type, q, ρ, N, μ; Dmin=Dmin, Dmax=Dmax, add_dry_aerosol_mass=add_dry_aerosol_mass)
end


"""
We are fitting a truncated gamma distribution:

    n(D) = n₀ D^μ exp(-λ D)     for D ∈ [0, Dmax]

where each particle has mass:

    m(D) = c_mass D^me 

We require the distribution to satisfy:

    N = ∫₀ᴰᵐᵃˣ n(D) dD                  (number concentration)
    q = ρ ∫₀ᴰᵐᵃˣ n(D) m(D) dD          (mass mixing ratio)

Substituting m(D):

    q = ρ * c_mass * ∫₀ᴰᵐᵃˣ n₀ D^(μ + me) exp(-λ D) dD

Let:

    I_μ   = ∫₀ᴰᵐᵃˣ D^μ       exp(-λ D) dD
    I_μme = ∫₀ᴰᵐᵃˣ D^(μ + me) exp(-λ D) dD

Then:

    N = n₀ * I_μ
    q = ρ * c_mass * n₀ * I_μme

Eliminating n₀ gives:

    q / (ρ * N * c_mass) = I_μme / I_μ

This defines the residual function used to solve for λ:

    residual(λ) = λ^(-me) * Γ(μ + me + 1, λ Dmax) / Γ(μ + 1, λ Dmax) - target

where:

    target = q / (ρ * N * c_mass)

is the normalized mean mass per droplet (modulo D^me scaling).
"""
function solve_n0_lambda_upper_truncated_gamma( # assumes Dmin = 0, so you can bound yourself accordingly
    prs::ACMP,
    q::FT,
    q_type::CMTWaterTypes,
    N::FT,
    Dmax::FT,
    μ::FT,
    me::FT,
    ρ::FT,
    c_mass::FT, # should be, mass = c_mass * D^me, so c_mass = χm * m0 * r0^(-me) # this is the mass prefactor, so we can use it to scale the mass of the particles
    ;
    tol::FT = FT(1e-6),
    maxiter::Int = 1000
) where {FT <: Real}

    # note, Dmax should be > 0

    # target = q / (ρ * N * c_mass)
    target = q * ρ / (N * c_mass)

    if isinf(Dmax)
        # Closed form solution for the full gamma distribution
        return get_n0_lambda(prs, q_type, q, ρ, N, μ; Dmax=Dmax) # fallback only
    end
        
    function residual(λ)
        # if λ <= FT(0)
        #     return FT(NaN)
        # end
        x = λ * Dmax
        num = Γ_lower(μ + me + 1, x)
        den = Γ_lower(μ + 1, x)
        if iszero(den)  # probably means very small x, just fall back to x = 0 [ so we can't divide by 0... what do we do?]
            #=  we take the limit...
                As x = λ * Dmax → 0, the lower incomplete gamma functions 
                γ(μ + 1, x) and γ(μ + me + 1, x) both → 0, so their ratio becomes 
                indeterminate (0/0). But we can resolve this by using the series 
                expansion for small x:
                    γ(a, x) ≈ x^a / a   as x → 0
                Thus,
                    γ(μ + me + 1, x) / γ(μ + 1, x)  ≈  (x^(μ + me + 1) / (μ + me + 1)) / (x^(μ + 1) / (μ + 1))  = x^me * (μ + 1) / (μ + me + 1)
                Since x = λ * Dmax, we have:
                    x^me * λ^(-me) = (λ * Dmax)^me * λ^(-me) = Dmax^me
                So the full fallback approximation becomes:
                    λ^(-me) * γ(μ + me + 1, x) / γ(μ + 1, x) ≈ Dmax^me * (μ + 1) / (μ + me + 1)
                This gives a smooth and finite limit as x → 0 and avoids division by zero.
            =#
            num_o_den = Dmax^(me) * (μ+1) / (μ + me + 1) # this is the limit of the ratio as x -> 0, so we can use it to avoid division by 0
        elseif isnan(den)  # probably means very big x, just fall back to x = Inf |  use a fallback... at very large λ, just use infinity

            num = CM1.SF.gamma(μ + me + 1)
            den = CM1.SF.gamma(μ + 1)
            num_o_den = num / den
        else
            num_o_den = num / den
        end

        λ_to_neg_me = λ^(-me)

        if isinf(λ_to_neg_me)
            return floatmax(FT)
        else
            return λ_to_neg_me * num_o_den - target
        end
    end

    λ_lower = FT(0) # lower bound is 0

    n0_upper = n0(prs, q, ρ, q_type, N, μ; Dmin=FT(0), Dmax=FT(Inf))
    λ_upper = lambda(prs, q_type, q, ρ, N, μ; _n0=n0_upper, Dmin=FT(0), Dmax=FT(Inf)) # upper bound is the closed form solution for the full gamma distribution


    methods_to_try = (TD.RS.RegulaFalsiMethod, TD.RS.BrentsMethod, TD.RS.SecantMethod, TD.RS.BisectionMethod) # False is allegedly the fastest? Secant can converge out of bounds, Bisection is slow...
    for method in methods_to_try
        # Use RootSolvers.jl to find λ [[ it's already loaded in Thermodynamics ]]
        sol = TD.RS.find_zero(
            residual,
            method(λ_lower, λ_upper),
            TD.RS.CompactSolution(),
            TD.RS.ResidualTolerance(tol * target),
            maxiter
        )
        if sol.converged
            converged = true
            λ_est = sol.root
            x = λ_est * Dmax
            I_μ = λ_est^(-μ - 1) * Γ_lower(μ + 1, x)
            n₀ = N / I_μ
            return n₀, λ_est # The main return!!!
        end
    end

    error("No solution found")

    # Fallback if no method converged, to just using Dmax = Inf
    return get_n0_lambda(prs, q_type, q, ρ, N, μ; Dmax=FT(Inf)) # fallback only

end

"""
    verify_truncated_gamma_fit(λ, n₀, Dmax; μ=2.0, me=3.0, ρ=1.0, c_mass=1)

Returns `(q_model, N_model)` implied by a truncated gamma distribution with parameters `(λ, n₀)`.

The size distribution is assumed to follow:

    n(D) = n₀ D^μ e^{-λ D},   D ∈ [0, Dmax]

We compute:
  - N_model = total number concentration
  - q_model = mass mixing ratio [kg/kg]

The droplet mass is assumed to scale with diameter as:

    m(D) = c_mass D^me

where `c_mass` is a constant prefactor, and `me` is the mass exponent (typically 3 for spheres).

Internally, this evaluates:

    N_model = n₀ * I_μ
    q_model = 1/ρ * c_mass * n₀ * I_μd

where:
    I_μ  = λ^(-μ - 1)       * γ(μ + 1, λ Dmax)
    I_μd = λ^(-μ - me - 1)  * γ(μ + me + 1, λ Dmax)

If `Dmax` is infinite, the complete gamma function is used instead of the incomplete one.
"""

function verify_truncated_gamma_fit(λ::FT, n₀::FT, Dmax::FT;
    μ::FT = FT(0),
    me::FT = FT(3),
    ρ::FT = FT(1),
    c_mass::FT = FT(1)
    ) where {FT <: Real}

    if isinf(Dmax)
        I_μ  = λ^(-μ - 1)       * CM1.SF.gamma(μ + 1)
        I_μd = λ^(-μ - me - 1)   * CM1.SF.gamma(μ + me + 1)
    else
        x = λ * Dmax
        I_μ  = λ^(-μ - 1)       * Γ_lower(μ + 1, x)
        I_μd = λ^(-μ - me - 1)  * Γ_lower(μ + me + 1, x)
    end

    N_model = n₀ * I_μ
    q_model = (1/ρ) * n₀ * c_mass * I_μd
    return (q_model, N_model)
end


"""
    Check whether a solution is possible at all for the desired q, N, μ
    returns (is_feasible::Bool, qmax::FT)
"""
function check_N_q_μ_feasibility(
    q::FT,
    N::FT,
    μ::FT,
    c_mass::FT, # should be equal to χm * m0 * r0 ^(-(m_e + Δm))
    me::FT; # should be (me + Δm) combined
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf)
) where {FT <: Real}

    if iszero(Dmin)
        if isinf(Dmax)
            return true, Inf
        else
            C = N * c_mass # C * χm * m₀ * r₀^(-m_e)
            q_max = C * (μ + 1) / (μ + me + 1) * Dmax^me # asymptotes to C * r_max^me = N * mass_of_one_particle at μ=∞ This is because at μ=∞, we have a delta spike, so we can place it anywhere and capture anything apparently, just by scaling it up and down.
            return FT(0) < q < q_max, q_max
            # implement here
        end
    else
        error("Not implemented yet")
    end
end

function check_N_q_μ_feasibility(
    param_set::APS,
    q::FT,
    q_type::CMTWaterTypes,
    N::FT,
    μ::FT
    ;
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT <: Real}
    prs = TCP.microphysics_params(param_set)
    _χm = get_χm(param_set, q_type)


    _m0::FT = m0(prs, q_type)
    _r0::FT = r0(prs, q_type)
    _me::FT = me(prs, q_type)
    _Δm::FT = Δm(prs, q_type)
    c_mass = _χm * _m0 * _r0^(-(_me + _Δm)) # should be equal to χm * m0 * r0 ^(-(m_e + Δm))
    return check_N_q_μ_feasibility(q, N, μ, c_mass, _me; Dmin=Dmin, Dmax=Dmax)
end



function minimum_rmax_for_solution(q::FT, N::FT, μ::FT, c_mass::FT, me::FT) where {FT <: Real}
    factor = q / (N * c_mass) * (μ + me + 1) / (μ + 1)
    return factor^(1 / me)
end

# might be a simpler solver than the one above
function truncated_gamma_solver(
    q::FT,
    N::FT, 
    ρ::FT,
    r_max::FT,
    c_mass::FT, # should be χm * m0 * r0^(-(me + Δm))
    ;
    μ::FT = FT(0),
    # χm::FT = 1.0,
    # m₀::FT = 1.0,
    # r₀::FT = 1.0,
    me::FT = FT(3), # should be me + Δm
    tol::FT = FT(1e-6),
    maxiter::Int = 100,
    λ_lower::FT = FT(1e-5), # we see some potentially spurious tiny roots we should exclude...
    λ_upper::FT = FT(1e10), # need to set this up soon to use the old lambda perhaps?
    ) where {FT <: Real}

    # ν = m_e + Δm
    const_term = N * c_mass; # c_mass = χm * m₀ * r₀^(-(me + Δm))

    function lower_incomplete_gamma(a, x)
        P, _ = CM1.SF.gamma_inc(a, x)
        return P * CM1.SF.gamma(a)
    end

    function residual(λ)
        x = λ * r_max
        γ1 = lower_incomplete_gamma(μ + 1, x)
        γ2 = lower_incomplete_gamma(μ + me + 1, x)
        q_pred = const_term * (γ2 / (λ^me * γ1))
        return q_pred - q
    end

    @info "(residual(λ_lower), residual(λ_upper)) = $((residual(λ_lower), residual(λ_upper)))"

    try
        λ_guess = 5.0 / r_max
        # λ = Roots.find_zero(residual, (1e-6, 100*λ_guess), verbose=true)
        # λ = Roots.find_zero(residual, (λ_lower, λ_upper), verbose=true)
        λ = TD.RS.find_zero(
            residual,
            TD.RS.SecantMethod(λ_lower, λ_upper),
            TD.RS.VerboseSolution(),
            TD.RS.ResidualTolerance(tol * q),
            maxiter
        )

        @info "root should be : $λ"
    catch e
        @warn e
        @info "no root prediciton"
    end
    println()
    

    # ---------------------------------------------------------------------- #
    # n0_upper = n0(prs, q, ρ, N, q_type; μ=μ, Dmin=FT(0), Dmax=FT(Inf))
    # λ_upper = lambda(prs, q_type, q, ρ, N; _n0=n0_upper, μ=μ, Dmin=FT(0), Dmax=FT(Inf)) # upper bound is the closed form solution for the full gamma distribution

    methods_to_try = (TD.RS.RegulaFalsiMethod, TD.RS.BrentsMethod, TD.RS.SecantMethod, TD.RS.BisectionMethod) # False is allegedly the fastest? Secant can converge out of bounds, Bisection is slow...
    for method in methods_to_try
        # Use RootSolvers.jl to find λ [[ it's already loaded in Thermodynamics ]]
        @info "Trying method: $method"
        sol = TD.RS.find_zero(
            residual,
            method(λ_lower, λ_upper),
            # method(FT(1e-6), FT(100*5.0) / r_max),
            TD.RS.VerboseSolution(),
            TD.RS.ResidualTolerance(tol * q),
            maxiter
        )
        print(sol)

        if sol.converged
            λ = sol.root
            γ1 = lower_incomplete_gamma(μ + 1, λ * r_max)
            n₀ = N * λ^(μ + 1) / γ1
            return λ, n₀
        else
            #
        end
    end
    # ---------------------------------------------------------------------- #

    error("no sol'n found")
end


