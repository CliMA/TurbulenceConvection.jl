
# maybe add an oversider dispatch fcn here?

# function NR_monodisperse(N::FT, q::FT; r_0::FT = 0.2 * 1e-6, ρ_w = 1000) where {FT}
#     """
#     constant N
#     """
#     r = (q / ((4 / 3) * π * ρ_w * N))^(1 / 3)
#     r = max(r, 0.2 * 1e-6) # bound to be at least ~micron size...something like kohler crit radius (e.g. if the cloud is new, diagnosed r is 0 for selecting tau in principle)
#     return (; N, r)
# end


# function NR_fixed_radius(r::FT, q::FT; ρ_w = 1000) where {FT}
#     """
#     constant r
#     """
#     N = q / ((4 / 3) * π * r^3 * ρ_w)
#     return (; N, r)
# end




# function NR_inhomogeneous_mixing_liquid(
#     param_set::APS,
#     N::FT,
#     p::FT,
#     q_liq::FT,
#     p_LCL::FT,
#     T_LCL::FT,
#     q_LCL::FT,
# ) where {FT}
#     """
#     # add some DOI references

#     # Diagnose q_l from adiabatic parcel ascent from LCL to this altitude
#     # Calculate r from this q_l and our assumed N
#     # Calculate real n given our true q_l and this r -- implication is we're subadiabatic and n is depleted by entrainment.

#     In a mixed-phase cloud then, WBF would deplete q_l even further leading to reduced n -- an interesting idea considering the Tan et al style pocketing idea of mixed phase clouds, maybe this is better than reducing r homogeneously.

#     We have no similar method for ice because we have no expectation of an N -- INP are continually activated as a function of temperature.
#     The ice N/R relationship maybe is best left as being uniform  since earlier particles would have grown more but there's less of them (don't know relative scaling)... or some other data informed parameterization...
#     A potential avenue for thinking about equilibrium adiabat q_l, q_i is Korolev, Field 2006 https://doi.org/10.1175/2007JAS2355.1, but then it would be horrible to reconvert that back to some representation of anything for ice...

#     What we'd really want is a lagrangian parcel model that evolves the ice size distribution that we can compare to -- couldn't find obs of ice size dist as fcn of height...



#     what happens if the previous LCL is above where we are now? is the calculation order correct to resolve this? (then the cloud base might be at the toa... not ideal...)

#     """
#     # theta_li is conserved so we just need to back out states based on theta_li and our two pressures
#     if p_LCL < p # LCL is above current condensation, so new cloud... use default monodisperse relation
#         return NR_monodisperse(N, q_liq)
#     else
#         q_LCL = TD.PhasePartition(q_LCL) # assume all vapor at cloud base (or should I read in the real value?) would have to make q_pt an argument rather than q_LCL or only accept ts for example...
#         θ_liq_ice = TD.liquid_ice_pottemp_given_pressure(param_set, T_LCL, p_LCL, q_LCL)
#         # @show(p, θ_liq_ice)
#         q_adiabatic = PhaseEquil_pθq_q(param_set, p, θ_liq_ice, q_LCL.tot)# can we ascend the adiabat here without having to do a full solve? -- also how do I partition liquid and ice in this framework? assume one phase?e
#         ρ_w = FT(1)
#         r_adiabatic = (q_adiabatic.liq / (4 / 3 * π * ρ_w * N))^(1 / 3)
#         r_adiabatic = max(r_adiabatic, 0.2 * 1e-6) # bound to be at least ~micron size...something like kohler crit radius (e.g. if the cloud is new, diagnosed r is 0 for selecting tau in principle) (also if q_liq is 0 so is r_liq so rly need some bounds)
#         N = q_liq / (4 / 3 * π * r_adiabatic^3 * ρ_w)
#         return (; N, r = r_adiabatic)
#     end
# end

# function NR_inhomogeneous_mixing_liquid(
#     param_set::APS,
#     N::FT,
#     p::FT,
#     q_liq::FT,
#     ts_LCL::TD.ThermodynamicState,
# ) where {FT}
#     p_LCL = TD.air_pressure(param_set, ts_LCL)
#     T_LCL = TD.air_temperature(param_set, ts_LCL)
#     q_LCL = TD.total_specific_humidity(param_set, ts_LCL)
#     return NR_inhomogeneous_mixing_liquid(param_set, N, p, q_liq, p_LCL, T_LCL, q_LCL)
# end

# function NR_inhomogeneous_mixing_liquid(
#     param_set::APS,
#     N::FT,
#     ts::TD.ThermodynamicState,
#     q_liq::FT,
#     ts_LCL::TD.ThermodynamicState,
# ) where {FT}
#     p = TD.air_pressure(param_set, ts)
#     p_LCL = TD.air_pressure(param_set, ts_LCL)
#     T_LCL = TD.air_temperature(param_set, ts_LCL)
#     q_LCL = TD.total_specific_humidity(param_set, ts_LCL)
#     return NR_inhomogeneous_mixing_liquid(param_set, N, p, q_liq, p_LCL, T_LCL, q_LCL)
# end


# """
# Phase equil but the output is always liquid, pretend ice doesnt exist even below freezing...
# """
# function PhaseEquil_pθq_q(
#     # param_set::APS,
#     param_set::APS, # a thermo_params to pass in
#     p::FT,
#     θ_liq_ice::FT,
#     q_tot::FT,
#     maxiter::IT = nothing,
#     relative_temperature_tol::FTT = nothing,
#     ::Type{sat_adjust_method} = TD.RS.SecantMethod,
#     T_guess::Union{FT, Nothing} = nothing,
# ) where {FT <: Real, IT <: TD.ITERTYPE, FTT <: TD.TOLTYPE(FT), sat_adjust_method}
#     maxiter === nothing && (maxiter = 50)
#     relative_temperature_tol === nothing && (relative_temperature_tol = FT(1e-4))
#     phase_type = TD.PhaseEquil{FT} # equil or non equil
#     q_tot_safe = clamp(q_tot, FT(0), FT(1))
#     T = saturation_adjustment_given_pθq( # use my local version... force liq_frac to 1
#         sat_adjust_method,
#         param_set,
#         p,
#         θ_liq_ice,
#         q_tot_safe,
#         phase_type,
#         maxiter,
#         relative_temperature_tol,
#         T_guess,
#     ) # finds T based on liq/ice potential temperature and q_tot ... afterwards it's up to you to back out the phase partitioning I guess
#     q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type) # our version of this uses liq_frac always 1, and in principle I guess T from saturation adjustment didn't care about what the phase partition looked like? (that doesnt feel right but...) edit -- replaced with my own version whew
#     return q_pt
#     # ρ = TD.air_density(param_set, T, p, q_pt)
#     # e_int = TD.internal_energy(param_set, T, q_pt)
#     # return TD.PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T) # for our output we want the 
# end


# function PhasePartition_equil_given_p(
#     param_set::APS,
#     T::FT,
#     p::FT,
#     q_tot::FT,
#     ::Type{phase_type},
# ) where {FT <: Real, phase_type <: TD.ThermodynamicState}

#     q_v_sat = TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
#     # _liquid_frac = liquid_fraction(param_set, T, phase_type;ramp=false)
#     _liquid_frac = 1 #force liquid for our adiabtic parcel test...
#     q_c = q_tot - q_v_sat
#     q_liq = _liquid_frac * q_c
#     q_ice = (1 - _liquid_frac) * q_c
#     return TD.PhasePartition(q_tot, q_liq, q_ice)
# end


# # adjust this one to make sure it uses our local PhasePartition_equil_given_p which always has liq_frac = 1
# function saturation_adjustment_given_pθq(
#     ::Type{sat_adjust_method},
#     param_set::APS,
#     p::FT,
#     θ_liq_ice::FT,
#     q_tot::FT,
#     ::Type{phase_type},
#     maxiter::Int,
#     relative_temperature_tol::FT,
#     T_guess::Union{FT, Nothing} = nothing,
# ) where {FT <: Real, sat_adjust_method, phase_type <: TD.PhaseEquil}
#     tol = TD.RS.RelativeSolutionTolerance(relative_temperature_tol)
#     T_min::FT = TD.TP.T_min(param_set)
#     T_freeze::FT = TD.TP.T_freeze(param_set)
#     cp_d::FT = TD.TP.cp_d(param_set)
#     cp_v::FT = TD.TP.cp_v(param_set)
#     air_temp(q) = TD.air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
#     function θ_liq_ice_closure(T)
#         # use local version of PhasePartition_equil_given_p which always has liq_frac = 1
#         q_pt = PhasePartition_equil_given_p(param_set, T, oftype(T, p), oftype(T, q_tot), phase_type)
#         return TD.liquid_ice_pottemp_given_pressure(param_set, T, oftype(T, p), q_pt)
#     end
#     q_vap_sat(T) = TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
#     T_1 = max(T_min, air_temp(TD.PhasePartition(q_tot))) # Assume all vapor
#     q_v_sat_1 = q_vap_sat(T_1)
#     unsaturated = q_tot <= q_v_sat_1
#     if unsaturated && T_1 > T_min
#         return T_1
#     end
#     temperature_tol = T_freeze * relative_temperature_tol
#     θ_liq_ice_upper = θ_liq_ice_closure(T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
#     θ_liq_ice_lower = θ_liq_ice_closure(T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
#     if θ_liq_ice_lower < θ_liq_ice < θ_liq_ice_upper
#         return T_freeze
#     end
#     roots(T) = oftype(T, θ_liq_ice) - θ_liq_ice_closure(T)
#     sol = TD.RS.find_zero(
#         roots,
#         TD.sa_numerical_method_pθq(sat_adjust_method, param_set, p, θ_liq_ice, q_tot, phase_type, T_guess),
#         TD.RS.CompactSolution(),
#         tol,
#         maxiter,
#     )
#     if !sol.converged
#         if TD.print_warning()
#             TD.KA.@print("-----------------------------------------\n")
#             TD.KA.@print("maxiter reached in saturation_adjustment_given_pθq:\n")
#             TD.print_numerical_method(sat_adjust_method)
#             TD.print_T_guess(sat_adjust_method, T_guess)
#             TD.KA.@print(", p=", p)
#             TD.KA.@print(", θ_liq_ice=", θ_liq_ice)
#             TD.KA.@print(", q_tot=", q_tot)
#             TD.KA.@print(", T=", sol.root)
#             TD.KA.@print(", maxiter=", maxiter)
#             TD.KA.@print(", tol=", tol.tol, "\n")
#         end
#         if TD.error_on_non_convergence()
#             error("Failed to converge with printed set of inputs.")
#         end
#     end
#     return sol.root
# end


# function cloud_base(aux, grid, ts, mode)
#     """
#     adopted from driver/compute_diagnostics.jl

#     really just need to map from from cloud base p,T,q to local p,T,q

#     either runs in each updraft or in the environment so doesn't need to worry about that structure
#     """


#     k_end = collect(real_center_indices(grid))[end]

#     if mode === :env
#         cloud_base = zc_toa(grid).z
#         cloud_base_ts = ts[k_end] # i think it goes bottom to top but need to check...
#         @inbounds for k in real_center_indices(grid)
#             if aux.area[k] > 1e-3
#                 if TD.has_condensate(aux.q_liq[k] + aux.q_ice[k])
#                     if cloud_base > grid.zc[k].z
#                         cloud_base = grid.zc[k].z
#                         cloud_base_ts = ts[k]
#                     end
#                 end
#             end
#         end
#         # Note definition of cloud cover : each updraft is associated with
#         # a cloud cover equal to the maximum area fraction of the updraft
#         # where ql > 0. Each updraft is assumed to have maximum overlap with
#         # respect to itup (i.e. no consideration of tilting due to shear)
#         # while the updraft classes are assumed to have no overlap at all.
#         # Thus total updraft cover is the sum of each updraft's cover

#     elseif mode === :up
#         cloud_base = zc_toa(grid).z
#         cloud_base_ts = ts[k_end]
#         @inbounds for k in real_center_indices(grid)
#             if TD.has_condensate(aux.q_liq[k] + aux.q_ice[k]) && aux.area[k] > 1e-6
#                 if cloud_base > grid.zc[k].z
#                     cloud_base = grid.zc[k].z
#                     cloud_base_ts = ts[k]
#                 end
#             end
#         end
#     end

#     return (; cloud_base_ts, cloud_base)

# end


# ======================================================================================================================== #

# CloudMicrophysics.jl:Microphysics1M tie-ins

function get_r_cond_precip(param_set::APS, q_type::CMTWaterTypes)
    FT = eltype(param_set)

    if (q_type isa CMT.LiquidType) || (q_type isa CMT.RainType)
        return FT(param_set.user_params.r_liq_rain)
    elseif (q_type isa CMT.IceType) || (q_type isa CMT.SnowType)
        return FT(CMP.r_ice_snow(TCP.microphysics_params(param_set)))
    else
        error("Unknown water type for r_threshold: $q_type")
    end
end
get_r_cond_precip(microphys_params::ACMP, q_type::Union{CMT.IceType, CMT.SnowType}) = CMP.r_ice_snow(microphys_params)

function get_autoconversion_timescale(microphys_params::ACMP, q_type::CMTWaterTypes)
    # FT = eltype(microphys_params)
    if (q_type isa CMT.LiquidType) || (q_type isa CMT.RainType)
        return CMP.τ_acnv_rai(microphys_params)
    elseif (q_type isa CMT.IceType) || (q_type isa CMT.SnowType)
        return CMP.τ_acnv_sno(microphys_params)
    else
        error("Unknown water type for autoconversion timescale: $q_type")
    end
end
get_autoconversion_timescale(param_set::APS, q_type::CMTWaterTypes) = get_autoconversion_timescale(TCP.microphysics_params(param_set), q_type)

function get_χm(param_set::APS, q_type::CMTWaterTypes)
    FT = eltype(param_set)
    if (q_type isa CMT.LiquidType)
        microphys_params = TCP.microphysics_params(param_set)
        return FT(param_set.user_params.χm_liq)
    else
        microphys_params = TCP.microphysics_params(param_set)
        return FT(CM1.χm(microphys_params, q_type))
    end
end
get_χm(microphys_params::ACMP, q_type::Union{CMT.IceType, CMT.AbstractPrecipType}) = χm(microphys_params, q_type)

function get_χa(param_set::APS, q_type::CMTWaterTypes)
    microphys_params = TCP.microphysics_params(param_set)
    return χa(microphys_params, q_type)
end
get_χa(microphys_params::ACMP, q_type::Union{CMT.AbstractPrecipType}) = χa(microphys_params, q_type)

# ------------------------------- #
function particle_mass(microphys_params::ACMP, q_type::CMTWaterTypes, r::FT, χm::FT) where FT
    """
    Mass of a single droplet.
    Not defined in CM1 for liquid
    """
    _r0::FT = r0(microphys_params, q_type)
    _m0::FT = m0(microphys_params, q_type)
    _me::FT = me(microphys_params, q_type)
    _Δm::FT = Δm(microphys_params, q_type)
    return χm * _m0 * (r/_r0)^(_me + _Δm)
end
particle_mass(microphys_params::ACMP, q_type::Union{CMT.IceType, CMT.AbstractPrecipType}, r::FT) where FT = particle_mass(microphys_params, q_type, r, CM1.χm(microphys_params, q_type))
particle_mass(param_set::APS, q_type::CMTWaterTypes, r::FT) where FT = particle_mass(TCP.microphysics_params(param_set), q_type, r, get_χm(param_set, q_type))

# multiple droplets
mass(microphys_params::ACMP, q_type::CMTWaterTypes, r::FT, N::FT, χm::FT; monodisperse::Bool=true) where FT = monodisperse ? N * particle_mass(microphys_params, q_type, r, χm) : error("Without making assumptions about the distribution, we cannot calculate the mass. If you mean to use gamma distribution, needs param_set, use q_from_rN")
mass(microphys_params::ACMP, q_type::Union{CMT.IceType, CMT.AbstractPrecipType}, r::FT, N::FT; monodisperse::Bool=true) where FT = monodisperse ? N * particle_mass(microphys_params, q_type, r, CM1.χm(microphys_params, q_type)) : error("Without making assumptions about the distribution, we cannot calculate the mass. If you mean to use gamma distribution, needs param_set, use q_from_rN")
mass(param_set::APS, q_type::CMTWaterTypes, r::FT, N::FT; monodisperse::Bool=true) where FT = monodisperse ? N * particle_mass(TCP.microphysics_params(param_set), q_type, r, get_χm(param_set, q_type)) : q_from_rN(param_set, q_type, r, N; monodisperse=false, ρ=FT(1), μ=FT(0))

# ------------------------------- #
function particle_radius_from_mass(microphys_params::ACMP, q_type::CMTWaterTypes, m::FT, χm::FT; ) where FT
    """
    Radius of a single droplet given its mass.
    m(r) = χm * m0 * (r/r0)^(me + Δm) --> r(m) = r0 * (m / (χm * m0))^(1 / (me + Δm))
    """


    _m0::FT = m0(microphys_params, q_type)
    _me::FT = me(microphys_params, q_type)
    _Δm::FT = Δm(microphys_params, q_type)
    _r0::FT = r0(microphys_params, q_type)

    return _r0 * (m / (χm * _m0))^(1 / (_me + _Δm))
end
particle_radius_from_mass(microphys_params::ACMP, q_type::Union{CMT.IceType, CMT.AbstractPrecipType}, m::FT; ) where FT = particle_radius_from_mass(microphys_params, q_type, m, CM1.χm(microphys_params, q_type))
particle_radius_from_mass(param_set::APS, q_type::CMTWaterTypes, m::FT; ) where FT = particle_radius_from_mass(TCP.microphysics_params(param_set), q_type, m, get_χm(param_set, q_type))

# multiple droplets
function radius_from_mass(microphys_params::ACMP, q_type::CMTWaterTypes, q::FT, N::FT, χm::FT; monodisperse::Bool = true, ρ::FT=FT(1), Dmin::FT = FT(0), Dmax::FT = FT(Inf), μ::FT=FT(0)) where FT
    if monodisperse
        return particle_radius_from_mass(microphys_params, q_type, q / N, χm)
    else
        if iszero(Dmin) && isinf(Dmax)
            # we know <q> = ...
            _m0::FT = m0(microphys_params, q_type)
            _me::FT = me(microphys_params, q_type)
            _Δm::FT = Δm(microphys_params, q_type)
            _r0::FT = r0(microphys_params, q_type)

            factor = (χm * _m0 / ρ) * _r0^(-(_me + _Δm))
            G_num = CM1.SF.gamma(μ + _me + _Δm + 1)
            G_den = CM1.SF.gamma(μ + 1)

            λ = (factor * G_num / G_den / (q / N))^(1 / (_me + _Δm))

            return (μ + 1) / λ
            error("haven't implemented this yet. Need to find the closed form sol'n for q, then solve for n0, λ, then <r>")
        else
            error("Without making assumptions about the distribution, we cannot calculate the mass of a single droplet")
        end
    end
end
radius_from_mass(microphys_params::ACMP, q::FT, N::FT, χm::FT; monodisperse::Bool = true) where FT = monodisperse ? particle_radius_from_mass(microphys_params, CMT.LiquidType(), q / N, χm) : error("Without making assumptions about the distribution, we cannot calculate the mass of a single droplet")
radius_from_mass(param_set::APS, q_type::CMTWaterTypes, q::FT, N::FT; monodisperse::Bool = true) where FT = particle_radius_from_mass(TCP.microphysics_params(param_set), q_type, q/N, get_χm(param_set, q_type))
# ------------------------------- #


function particle_area(microphys_params::ACMP, q_type::CMTWaterTypes, r::FT, χa::FT) where FT
    """
    Area of a single droplet.
    Not defined in CM1 for liquid or ice
    """

    _a0::FT = a0(microphys_params, q_type)
    _r0::FT = r0(microphys_params, q_type)
    _ae::FT = ae(microphys_params, q_type)
    _Δa::FT = Δa(microphys_params, q_type)
    return χa * _a0 * (r/_r0)^(_Δa + _ae)
end
particle_area(microphys_params::ACMP, q_type::CMT.AbstractPrecipType, r::FT) where FT = particle_area(microphys_params, q_type, r, CM1.χa(microphys_params, q_type))
particle_area(param_set::APS, q_type::CMTWaterTypes, r::FT) where FT = particle_area(TCP.microphysics_params(param_set), q_type, r, get_χa(param_set, q_type))

function area(microphys_params::ACMP, q_type::CMTWaterTypes, r::FT, N::FT, χa::FT; ρ::FT=FT(1), Dmin::FT = FT(0), Dmax::FT = FT(Inf), μ::FT=FT(0), monodisperse::Bool = true) where FT
    """
    Integrates the distribution ∫ n(r) a(r) dr
    """

    if monodisperse
        # If monodisperse, then N is the number of particles of radius r
        return N * particle_area(microphys_params, q_type, r, χa)
    else
        # Integrating the marshall palmer n(r) = n0 exp(-λr) , N = ∫ n(r) dr = n0/λ. But then we still don't know n0 and λ because we dont know q, there are, in fact, infinitely many solutions -- for each λ, a different q

        # strictly speaking, if we integrate 0 to inf, we know <r> = (μ+1)/(λ), N = n0/λ, so we do have both n0 and λ. If Dmin,Dmax are not 0,inf though, it becomes unsolvable
        if iszero(Dmin) && isinf(Dmax)
            # using <r> = (μ+1)/(λ)
            λ = (μ + 1) / r
            n0 = N * λ # N = n0 / λ
            q = int_nm_dr(microphys_params, q_type, n0, λ, ρ, μ; Dmin = Dmin, Dmax = Dmax)

        else
            # If Dmin or Dmax are not 0 or Inf, we cannot calculate the area without more information about the distribution.

            error("Without making assumptions about the distribution, we cannot calculate a from N and r")
        end
    end
end
area(microphys_params::ACMP, q_type::CMT.AbstractPrecipType, r::FT, N::FT;  monodisperse::Bool=true) where FT = area(microphys_params, q_type, r, N, CM1.χa(microphys_params, q_type); monodisperse = monodisperse)
area(param_set::APS, q_type::CMTWaterTypes, r::FT, N::FT; monodisperse::Bool=true) where FT = area(TCP.microphysics_params(param_set), q_type, r, N, get_χa(param_set, q_type); monodisperse = monodisperse)

# TODO: area(microphys_params::ACMP, q_type::CMTWaterTypes, r::FT, q::FT, χa::FT) which can be non-monodisperse...

# ------------------------------- #




# ======================================================================================================================== #

function q_is(
    param_set::APS,
    q_type::CMT.IceType,
    )
    """
    The mass of a droplet of the ice/snow threshold radius
    """
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    r_is = CMP.r_ice_snow(microphys_params)
    return q_from_r(param_set, q_type, r_is)
end


function q_lr(
    param_set::APS,
    q_type::CMT.LiquidType,
    )
    """
    No trivial equivalent for liquid because r_lr is in user_params, and we kind of added it ourselves.
    """
    FT = eltype(param_set)
    r_lr::FT = param_set.user_params.r_liq_rain # There is no r_lr for liquid...instead we just set it in user_params.
    return q_from_r(param_set, q_type, r_lr)

end

function r_ice_acnv(
    param_set::APS,
    r_acnv_scaling_factor = eltype(param_set)(NaN),
    )
    """
    r_ice_acnv = r_is * r_acnv_scaling_factor
    """
    FT = eltype(param_set)
    particle_min_radius::FT = param_set.user_params.particle_min_radius
    r_acnv_scaling_factor = isnan(r_acnv_scaling_factor) ? FT(param_set.user_params.r_ice_acnv_scaling_factor) : r_acnv_scaling_factor
    # return CMP.r_ice_snow(microphys_params) * FT(r_acnv_scaling_factor) # r_acnv_scaling_factor can become small enough that it's not allowed
    return particle_min_radius + (CMP.r_ice_snow(TCP.microphysics_params(param_set)) - particle_min_radius) * FT(r_acnv_scaling_factor)
end

function r_liq_acnv(
    param_set::APS,
    r_acnv_scaling_factor = eltype(param_set)(1),
    )
    """
    r_liq_acnv = r_lr * r_acnv_scaling_factor
    """
    FT = eltype(param_set)
    particle_min_radius::FT = param_set.user_params.particle_min_radius
    # return param_set.user_params.r_liq_rain * FT(r_acnv_scaling_factor) # r_acnv_scaling_factor can become small enough that it's not allowed
    return particle_min_radius + (param_set.user_params.r_liq_rain - particle_min_radius) * FT(r_acnv_scaling_factor)
end

function r_acnv(
    param_set::APS,
    q_type::CMTWaterTypes,
    r_acnv_scaling_factor = eltype(param_set)(NaN),
)
    """
    The acnv radius for the given q_type
    """
    FT = eltype(param_set)
    # microphys_params::ACMP = TCP.microphysics_params(param_set)
    if q_type isa CMT.IceType
        r_acnv_scaling_factor::FT = isnan(r_acnv_scaling_factor) ? FT(param_set.user_params.r_ice_acnv_scaling_factor) : FT(r_acnv_scaling_factor)
        return r_ice_acnv(param_set, FT(r_acnv_scaling_factor))
    elseif q_type isa CMT.LiquidType
        r_acnv_scaling_factor = isnan(r_acnv_scaling_factor) ? FT(1) : FT(r_acnv_scaling_factor)
        return r_liq_acnv(param_set, FT(r_acnv_scaling_factor))
    else
        error("Unknown q_type: $q_type")
    end
end



function q_acnv_0(
    param_set::APS,
    ice_type::CMT.IceType,
    r_acnv_scaling_factor = 1,
    )
    """
    The mass of a droplet of the acnv radius
    """
    FT = eltype(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    r_acnv = r_ice_acnv(param_set, FT(r_acnv_scaling_factor))
    return q_from_r(param_set, ice_type, r_acnv) # this is the q_acnv_0 at the given r_is * r_acnv_scaling_factor, with a minimum radius of r_min
end


function q_acnv_0(
    param_set::APS,
    liq_type::CMT.LiquidType,
    r_acnv_scaling_factor = 1,
)
    """
    The mass of a droplet of the acnv radius
    """
    FT = eltype(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    r_acnv = r_liq_acnv(param_set, FT(r_acnv_scaling_factor))
    return q_from_r(param_set, liq_type, r_acnv) # this is the q_acnv_0 at the given r_is * r_acnv_scaling_factor, with a minimum radius of r_min
end



"""
Get the implied N threshold based on r_is, and q
"""
function get_N_threshold(
    param_set::APS, # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.IceType,
    ;
    N = eltype(param_set)(NaN),
)
    FT = eltype(param_set)
    if isnan(N)
        # Let N be the number of droplets of radius r_is that means q = q_threshold
        microphys_params::ACMP = TCP.microphysics_params(param_set)
        q_is_here::FT = q_is(param_set, q_type)
        q_threshold::FT = CMP.q_ice_threshold(microphys_params)
        N = q_threshold / q_is_here # the implied number concentration at threshold for the largest possible droplets... so N and threshold are linked... ( i guess that makes sense? idk...)
    end

    return FT(N)

end


function get_q_threshold( # symmetric method for ice
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    ::CMT.LiquidType,
    ;
    N = eltype(param_set)(NaN),
) 
    # right now we don't have the threshold dependent on N, we just use relaxatin to a fixed thresh. We don't even have r_liq_rai so his is really just a fallback method 
    FT = eltype(param_set)
    return CMP.q_liq_threshold(TCP.microphysics_params(param_set))::FT
end


# function get_q_threshold(param_set::APS, q_type::CMT.LiquidType) # shorthand bc we don't have a fancy q_liq_threshold so we can save some arguments [[ is now redundant with the method above]]
#     return CMP.q_liq_threshold(TCP.microphysics_params(param_set))::eltype(param_set)
# end


"""
Get the implied q threshold based on r_is, and q 
"""
function get_q_threshold(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.IceType,
    ;
    N = eltype(param_set)(NaN),
) 

    FT = eltype(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)

    if isnan(N)
        return CMP.q_ice_threshold(microphys_params)::FT
    else
        N_thresh::FT = get_N_threshold(param_set, q_type; N = FT(N)) # if N was nothing, get_N_threshold() makes N_thresh * q_is = q_threshold (we could add an if block check for that to save fcn calls...)
        q_is_here::FT = q_is(param_set, q_type)
        return N_thresh * q_is_here  # q_is is the volume of a single ice crystal, so N * q_is = q_threshold
    end
end



"""
Get the implied q threshold based on r_is*r_acnv_scaling_factor, and q. This is a threshold that is based on the current q but with an r factor, rather than just being a fixed N_is * q_is value.
"""
function get_q_threshold_acnv(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.IceType,
    ;
    N = eltype(param_set)(NaN),
    assume_N_is::Bool = true,
    r_acnv_scaling_factor = eltype(param_set)(1.0),
)

    FT = eltype(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)

    # get N_is * q_acnv_0, if it's greater than q_ice_threshold, then use that, otherwise use q_ice_threshold

    q_ice_threshold::FT = CMP.q_ice_threshold(microphys_params)

    if isnan(N) && !assume_N_is # if N is NaN, but we don't assume we're at N_is, then we need to use the q_ice_threshold
        # since we assume we're at N_is, the thresh is just q * r_is_acnv_scaling_factor^3
        return q_ice_threshold # just take this as given

    else
        # We can't use `ice_dep_acnv_scaling_factor`  bc we allow for an r_min that impacts N, we have a separate scaling factor here.
        N_thresh::FT = get_N_threshold(param_set, q_type; N = FT(N)) # either take the given N or use the threshold N that assumes we're already at r_is
        q_acnv_0_here::FT = q_acnv_0(param_set, q_type, FT(r_acnv_scaling_factor)) # this is the q_acnv_0 at the given r_is * r_acnv_scaling_factor
        q_thresh_acnv::FT = N_thresh * q_acnv_0_here # q_is is the volume of a single ice crystal, so N * q_is = q_threshold
    end

    return max(q_thresh_acnv, q_ice_threshold) # make sure we don't get too low of a threshold... (this is a bit arbitrary but I think it makes sense to not let the threshold go below the q_ice threshold)
end



# copy of CM1.conv_q_liq_to_q_rai that uses our local get_q_threshold
function my_conv_q_liq_to_q_rai(
    param_set::APS,   # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition{FT};
    N::FT = FT(NaN),
) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    return max(0, q.liq - get_q_threshold(param_set, CMT.LiquidType(); N=N)) / CMP.τ_acnv_rai(microphys_params)
end

# copy of CM1.conv_q_ice_to_q_sno_no_supersat that uses our local get_q_threshold if we supply N otherwise just uses the default q_ice_threshold
function my_conv_q_ice_to_q_sno_no_supersat(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    qi::FT
    ;
    N::FT = FT(NaN),
    τ::FT = FT(NaN),
    assume_N_is::Bool = true,
    r_acnv_scaling_factor::FT = FT(1)
) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    if isnan(τ)
        τ = CMP.τ_acnv_sno(microphys_params)
    end
    return max(0, qi - get_q_threshold_acnv(param_set, CMT.IceType(); N=N, assume_N_is = assume_N_is, r_acnv_scaling_factor = r_acnv_scaling_factor)) / τ
end

my_conv_q_ice_to_q_sno_no_supersat(param_set::APS, q::TD.PhasePartition{FT}; N = FT(NaN), τ = FT(NaN), assume_N_is = true, r_acnv_scaling_factor = FT(1)) where {FT} = 
    my_conv_q_ice_to_q_sno_no_supersat(param_set, q.ice; N=FT(N), τ=FT(τ), assume_N_is=assume_N_is, r_acnv_scaling_factor=FT(r_acnv_scaling_factor))


function my_conv_q_ice_to_q_sno_no_supersat_simple(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    qi::FT,
    qi_threshold::FT,
    τ::FT = FT(NaN),
) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    if isnan(τ)
        τ = CMP.τ_acnv_sno(microphys_params)
    end
    return max(0, qi - qi_threshold) / τ
end

my_conv_q_ice_to_q_sno_no_supersat_simple(param_set::APS, q::TD.PhasePartition{FT}, qi_threshold::FT, τ::FT = FT(NaN)) where {FT} =  my_conv_q_ice_to_q_sno_no_supersat_simple(param_set, q.ice, qi_threshold, τ)


"""
Deposition-based autoconversion

M2005 leaves off the integral part and keeps only the part from growth at the threshold radius
Growth is slow enough I guess that growth above the thresh is much smaller than entire particles passing the thresh?

"""
function my_conv_q_ice_to_q_sno(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition{FT},
    T::FT,
    p::FT;
    N::FT = FT(NaN),
    # τ::FT = FT(NaN), # This doesn't get used, if we have N we just use that...
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    μ::FT = FT(NaN),
    thresh_only::Bool = false
) where {FT}

    acnv_rate = FT(0) 
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)

    ρ::FT = TD.air_density(thermo_params, T, p, q)
    _S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())


    if (q.ice > FT(0) && _S > FT(0))

        if isnan(μ)
            μ = μ_from_qN(param_set, ice_type, q.ice, N; ρ=ρ)
        end

        _G::FT = CM.Common.G_func(microphys_params, T, TD.Ice())


        _r_ice_snow::FT = CMP.r_ice_snow(microphys_params)
        # _n0::FT = n0(microphys_params, FT(0), ρ, ice_type)
        _n0::FT = isnan(N) ? n0(microphys_params, q.ice, ρ, ice_type) : n0(param_set, q.ice, ρ, ice_type, N; μ=μ, Dmin=Dmin, Dmax=Dmax) # use wanted N If given
        _me::FT = me(microphys_params, ice_type)
        _Δm::FT = Δm(microphys_params, ice_type)
        # _λ::FT = lambda(microphys_params, CT.IceType(), q.ice, ρ)
        _λ::FT = lambda(param_set, ice_type, q.ice, ρ, N, Dmin, Dmax; μ = μ, _n0 = _n0)

        # note term 1 is particles exactly at r_ice_snow, crossing over. term 2 is particles already larger than r_ice_snow.
        if thresh_only # only particles that are growing past r_ice_snow, not particles that were already larger.
            acnv_rate =
                4 * FT(π) * _S * _G * _n0 / ρ *
                (
                    _r_ice_snow^(μ + FT(2)) * exp(-_λ * _r_ice_snow) / (_me + _Δm)
                )
        else  # our default... include everything
            # acnv_rate =
            #     4 * FT(π) * _S * _G * _n0 / ρ *
            #     exp(-_λ * _r_ice_snow) *
            #     (
            #         _r_ice_snow^FT(2) / (_me + _Δm) +
            #         (_r_ice_snow * _λ + FT(1)) / _λ^FT(2)
            #     )

            # updated μ aware version
            acnv_rate =
                4 * FT(π) * _S * _G * _n0 / ρ *
                (
                    _r_ice_snow^(μ + FT(2)) * exp(-_λ * _r_ice_snow) / (_me + _Δm) +
                    CM1.SF.gamma(μ + FT(2), _λ * _r_ice_snow) / _λ^(μ + FT(2))
                )
        end

    end



    # return acnv_rate
    return resolve_nan(acnv_rate, FT(0)) # if we get a NaN, just return 0... (n0 and lambda can blow up for very small q which map to no conv)
end

"""
    Given a process acting on a gamma distribution, this function returns the fraction of the mass added to (or taken away from) the distribution that goes to particles **above** or **below** a specified threshold radius `r_th`.

    e.g. `DUM` definition at https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L2976 and https://github.com/NCAR/icar/blob/7fe80b4d86e24acae9925cf15f092bb6490ae228/src/physics/mp_morrison.f90#L2965 


    Compute the fraction of a process' action that goes to particles above or below
    a threshold radius r_th for a gamma size distribution.

    This docstring contains a full derivation so readers can check the math without
    referring to external notes.

    Overview
    --------
    We assume a gamma-shaped number distribution in radius:
        n(r) = n0 * r^μ * exp(-λ r),    r ≥ 0
    normalized so that ∫_0^∞ n(r) dr = N.

    A process acts on each particle with a weight proportional to r^k (k∈ℝ).
    Typical physical meanings of k:
    • k = 0 : number-weighted (each particle contributes equally)
    • k = 1 : deposition-like weighting when τ(r) ∝ 1/r (dm/dt ∝ r/τ ∝ r)
    • k = 2 : surface-area weighting (∝ r^2)
    • k = 3 : mass weighting for spheres (m ∝ r^3)

    We want the fraction of the *total* process (i.e. the integrated contribution)
    that goes to particles **above** a threshold r_th (or below, the complement).

    Derivation (general)
    --------------------
    Let the per-particle contribution be proportional to r^k. The total contribution
    over all radii is proportional to
        I_total = ∫_0^∞ r^k n(r) dr.

    The contribution from the tail r ≥ r_th is
        I_tail = ∫_{r_th}^∞ r^k n(r) dr.

    Hence the fraction **above** the threshold is
        F_>(r_th) = I_tail / I_total.

    Substitute the gamma form:
        I_total = n0 ∫_0^∞ r^{μ+k} e^{-λ r} dr
                = n0 * λ^{-(μ+k+1)} * Γ(μ+k+1).

        I_tail  = n0 ∫_{r_th}^∞ r^{μ+k} e^{-λ r} dr
                = n0 * λ^{-(μ+k+1)} * Γ(μ+k+1, ξ)

    where ξ ≡ λ r_th and Γ(s,x) is the upper incomplete gamma.
    Cancel the common prefactor n0 λ^{-(μ+k+1)} to obtain the closed ratio:

        F_>(r_th) = Γ(μ + k + 1, ξ) / Γ(μ + k + 1)   = Q(μ+k+1, ξ)

    and the fraction **below** is

        F_<(r_th) = 1 - F_>(r_th) = γ(μ+k+1, ξ) / Γ(μ+k+1)

    where γ(s,x) is the lower incomplete gamma and Q is the regularized upper gamma.

    Special case: μ = 0 and integer k
    ---------------------------------
    If μ = 0 and p ≡ k is a nonnegative integer (p∈ℕ), the incomplete gamma
    simplifies to a finite sum:

        Γ(p+1, ξ) = p! * e^{-ξ} * ∑_{i=0}^{p} ξ^i / i!

    Hence

        F_>(r_th) = e^{-ξ} * ∑_{i=0}^{p} ξ^i / i!,
        F_<(r_th) = 1 - e^{-ξ} * ∑_{i=0}^{p} ξ^i / i!.

    Examples (common k, μ=0):
    • k=0 (number):   F_> = e^{-ξ}
    • k=1 (deposition): F_> = e^{-ξ}(1 + ξ)
    • k=2 (area):     F_> = e^{-ξ}(1 + ξ + ξ^2/2)
    • k=3 (mass):     F_> = e^{-ξ}(1 + ξ + ξ^2/2 + ξ^3/6)

    How to obtain λ and n0 from bulk q and N
    ----------------------------------------
    To evaluate n(r_th) or compute moments you need λ and n0. For a specified
    mass law m(r) (for example spherical particles: m(r) = (4π/3) ρ_p r^3),
    use the definitions:

    N = ∫_0^∞ n(r) dr = n0 * λ^{-(μ+1)} * Γ(μ+1)
    M = q * ρ_a = ∫_0^∞ n(r) m(r) dr.

    If m(r) = χ_m m0 (r / r0)^p (general power-law mass), then

    M = n0 * χ_m m0 r0^{-p} * λ^{-(μ+p+1)} * Γ(μ + p + 1).

    Solve these two equations for λ and n0:

    λ^{p} = [ N * χ_m m0 r0^{-p} * Γ(μ + p + 1) ] / [ M * Γ(μ+1) ]

    so

    λ = ( [ N * χ_m m0 r0^{-p} * Γ(μ + p + 1) ] / [ M * Γ(μ+1) ] )^(1/p)

    and

    n0 = N * λ^{μ+1} / Γ(μ+1).

    (For the spherical case p=3 and χ_m m0 = 4π/3 ρ_p r0^3 this reduces to the
    λ^3 expression commonly used.)

    Numerical implementation notes
    ------------------------------
    • Use `SpecialFunctions.jl`:
        - `CM1.SF.gamma(s)`  → Γ(s)
        - `CM1.SF.gamma(s, x)` → Γ(s, x) (upper incomplete gamma)
    To get the regularized value use Γ(s, x) / Γ(s) or use `gammainc` if you prefer.

    • For μ = 0 and integer k you can use the finite-sum closed form (faster,
    exact and numerically stable).

    • For non-integer μ+k, use the incomplete gamma ratio:
        F_> = CM1.SF.gamma(μ+k+1, ξ) / CM1.SF.gamma(μ+k+1)

    • Units:
        - r in meters (m)
        - λ in m⁻¹
        - n(r) in #·m⁻³·m⁻¹ (number per m³ per m radius) — in code we treat n(r) as #/m³ per unit radius
        - q in kg/kg, ρ_a in kg/m³ → M = q * ρ_a is condensed mass per m³ (kg/m³)

    • Edge cases & robustness:
    - If μ+k+1 ≤ 0 the gamma integrals may diverge; enforce physically-plausible μ
        or use a finite `Dmax` cutoff in integrals.
    - For very large ξ (λ r_th), `e^{-ξ}` underflows; use the gamma ratio form
        which is numerically stable in both tails.
    - If the user provides non-integer k but μ=0 you can still compute F via the
        incomplete gamma (no closed form).

    Return value / API
    ------------------
    The function returns a single `FT` value (fraction). Keyword `below::Bool=true`
    controls whether the returned value is the fraction **below** (default) or **above**
    the threshold. If you call the SpecialFunctions gamma routines in the function,
    ensure `using SpecialFunctions` (or call via your `CM1.SF` alias).

    References
    ----------
    - Incomplete gamma identities and finite-sum forms: standard special-function texts.
    - Implementation note: use Γ(s,x)/Γ(s) for Q(s,x) (regularized upper gamma). [gamma_inc() works just fine actually]

"""
function get_process_threshold_fraction(
    param_set::APS,
    q_type::CMTWaterTypes,
    q::FT,
    N::FT,
    ρ_a::FT,
    r_th::FT;
    k::Int = Int(1), # default to radius weighting... for something like deposition which goes as dm/dt = S/τ where τ = (4πDr)⁻¹ ∝ r, so the mass contribution is linear in r. Then it's just a matter of doing the integral r^k n(r) and getting the fraction out.
    Dmax = FT(Inf),
    return_below::Bool = false, # default to returning fraction above.
    μ = FT(NaN),
    add_dry_aerosol_mass::Bool = false, # if true, add the dry aerosol mass to q when calculating λ and n0
) where {FT}

    microphys_params = TCP.microphysics_params(param_set)

    # μ is the shape parameter of the gamma distribution
    if isnan(μ)
        μ = μ_from_qN(param_set, q_type, q, N; ρ=ρ_a)
    end

    if add_dry_aerosol_mass
        q += mass(param_set, q_type, param_set.user_params.particle_min_radius, N; monodisperse=true) / ρ_a # add the dry aerosol mass to q when calculating λ and n0
    end

    # get the number distribution parameters
    _, λ = get_n0_lambda(microphys_params, q_type, q, ρ_a, N, μ; Dmax=Dmax, add_dry_aerosol_mass = false)

    ξ = λ*r_th

    if iszero(μ)
        # closed-form solution for integer k ≥ 0]

        if isone(k) # shortcut, see DUM in https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L2976
            F_below = one(FT) - exp(-ξ)*(1 + ξ)
            return return_below ? F_below : one(FT) - F_below
        end

        sum_val = zero(FT)
        for i in 0:k
            sum_val += ξ^i / factorial(i)
        end
        F_target = exp(-ξ) * sum_val
        return return_below ? one(FT) - F_target : F_target
    else
        # μ ≠ 0: use upper incomplete gamma function from SpecialFunctions.jl
        s = μ + k + 1
        # Γ_s = CM1.SF.gamma(s)
        # Γ_upper = CM1.SF.gamma(s, ξ)   # upper incomplete gamma
        # return return_below ? one(FT) - (Γ_upper / Γ_s) : (Γ_upper / Γ_s)
        P, Q = CM1.SF.gamma_inc(s, ξ) # regularized lower and upper incomplete gamma functions
        return return_below ? P : Q  # use the regularized lower gamma function for below
    end
end


"""
    Calculate the fraction of deposition that goes to snow directly from the sub_dep tendency.
    If τ_sub_dep is not a direct function of N and r (e.g. base, NN), then this is most useful, otherwise I think this should be similar to `my_conv_q_ice_to_q_sno()`
    One big difference though is that this does not include particles crossing r_is, only the growth of ice particles already above r_is.
    Also whatever is in `G` etc may not be exactly the same.

    The biggest other difference is this is fixed to the limiter we are using...
    By that I mean that if we're using the MM2015-EPA limiter, `my_conv_q_ice_to_q_sno()` does not respect that. this version does.
"""
function my_conv_q_ice_to_q_sno_by_fraction(
    param_set::APS,
    qi_tendency_sub_dep::FT,
    q::FT,
    N::FT,
    ρ_a::FT,
    add_dry_aerosol_mass::Bool = true,
) where {FT}
    acnv_rate = FT(0)
    if qi_tendency_sub_dep > FT(0)
        microphys_params = TCP.microphysics_params(param_set)
        r_is = CMP.r_ice_snow(microphys_params)
        fraction_to_snow = get_process_threshold_fraction(param_set, ice_type, q, N, ρ_a, r_is; k=1, return_below=false, add_dry_aerosol_mass=add_dry_aerosol_mass) # radius weighting for deposition
        acnv_rate = fraction_to_snow * qi_tendency_sub_dep
    end
    return resolve_nan(acnv_rate, FT(0)) # if we get a NaN, just return 0... (n0 and lambda can blow up for very small q which map to no conv)
end





"""
    get_fraction_below_or_above_thresh_from_qN(
        param_set, q_type, q, N, ρ_a, r_th;
        return_N=false, below_thresh=true, return_fraction=true, μ=NaN, Dmax=Inf
    )

Compute the fraction of mass or number of particles below or above a threshold radius `r_th`
for a gamma-distributed particle population.

Derivation:
1. Particle number distribution: n(r) = n0 * r^μ * exp(-λ r)
2. Particle mass: m(r) ∝ r^(me + Δm)
3. Total number: N = ∫_0^∞ n(r) dr
4. Total mass: q = (1/ρ_a) ∫_0^∞ n(r) m(r) dr
5. Fraction below threshold r_th:
   - Number fraction: f_N = ∫_0^{r_th} n(r) dr / N
   - Mass fraction:   f_q = ∫_0^{r_th} n(r) m(r) dr / q
6. Using gamma function identities:
   - f_N = γ(μ + 1, λ * r_th) / Γ(μ + 1)
   - f_q = γ(μ + p + 1, λ * r_th) / Γ(μ + p + 1)
   where γ(s,x) is the lower incomplete gamma function and p = me + Δm.


    Fraction is computed as:

        f_q(r_th) = γ(μ + me + Δm + 1, λ * r_th) / Γ(μ + me + Δm + 1) =  1 - Γ(μ + me + Δm + 1, λ * r_th) / Γ(μ + me + Δm + 1)
        f_N(r_th) = γ(μ + 1, λ * r_th) / Γ(μ + 1) =  1 - Γ(μ + 1, λ * r_th) / Γ(μ + 1)  # lower incomplete fraction

    where Γ(s, x) is the upper incomplete gamma function, γ is the lower incomplete gamma function,  and Γ the gamma function.

Parameters
----------
- `param_set::APS` : microphysics parameter set
- `q_type::CMTWaterTypes` : particle type (ice/liquid)
- `q::FT` : mass mixing ratio (kg/kg)
- `N::FT` : number concentration (#/m³)
- `ρ_a::FT` : air density (kg/m³)
- `r_th::FT` : threshold radius (m)
- `return_N::Bool` : compute number fraction if true, mass fraction if false
- `below_thresh::Bool` : fraction below threshold if true
- `return_fraction::Bool` : returns fraction if true, absolute q/N if false
- `μ::FT` : gamma shape parameter (computed if NaN)
- `Dmax::FT` : optional max particle diameter

Returns
-------
- Fraction (dimensionless) or absolute value (q*f or N*f)
"""

function get_fraction_below_or_above_thresh_from_qN(
    param_set::APS,
    q_type::CMTWaterTypes,
    q::FT,
    N::FT,
    ρ_a::FT,
    r_th::FT;
    return_N::Bool=false,
    below_thresh::Bool=true,
    return_fraction::Bool=true,
    μ::FT=FT(NaN),
    Dmax::FT=FT(Inf)
) where {FT}

    if !isfinite(q) || !isfinite(N)
        error("received invalid, non-finite q or N: q=$q; N=$N")
    end

    microphys_params = TCP.microphysics_params(param_set)

    # Determine μ if not provided
    if isnan(μ)
        μ = μ_from_qN(param_set, q_type, q, N; ρ=ρ_a)
    end

    # Get number distribution parameter λ
    _, λ = get_n0_lambda(microphys_params, q_type, q, ρ_a, N, μ; Dmax=Dmax)

    # Exponent k: number fraction -> 0, mass fraction -> me + Δm
    _me = me(microphys_params, q_type)
    _Δm = Δm(microphys_params, q_type)
    k = return_N ? FT(0) : (_me + _Δm)

    s = μ + k + one(FT)
    ξ = λ * r_th

    # get regularized incomplete gamma functions
    P, Q = CM1.SF.gamma_inc(s, ξ)

    fraction = below_thresh ? FT(P) : FT(Q)

    if !isfinite(fraction)
        error("computed non-finite fraction: fraction=$fraction; params: μ=$μ; k=$k; s=$s; ξ=$ξ; λ=$λ; r_th=$r_th")
    end


    if return_fraction
        # return fraction
        return fraction
    else
        return return_N ? (N * fraction) : (q * fraction)
    end
end














"""
Compute the mass crossing rate at threshold radius r_is for a gamma
particle distribution, given a known bulk deposition tendency.

Derivation
----------
Let the particle mass scale as
    m(r) = χ_m * m_0 * (r/r_0)^(m_e + Δm)

and let n(r) be a gamma distribution normalized to number concentration N:
    n(r) = n0 * r^μ * exp(-λ r),    n0 = N λ^(μ+1) / Γ(μ+1)

The deposition growth rate for a single particle is dm/dt = qi / tau_bulk (per bulk timescale).

Then, for particles crossing r_th:

    dq/dt |_cross = n(r_th) * m(r_th) * dr/dt

where dr/dt is obtained from the chain rule:

    dr/dt = (dm/dt) / (dm/dr)
          = (qi / tau_bulk) / (d/d r (m(r)))
          = (qi / tau_bulk) / ( (m_e + Δm) * m(r_th) / r_th )

Thus, the final expression:

    dq/dt |_cross = n(r_th) * m(r_th) * (qi / tau_bulk) / ((m_e + Δm) * m(r_th) / r_th)
                  = n(r_th) * (r_th / (m_e + Δm)) * (qi / tau_bulk)

This computes the mixing ratio tendency (kg/kg/s) crossing the threshold radius.

Arguments
---------
- tau_bulk :: FT   : total bulk deposition timescale (s)
- qi       :: FT   : bulk condensed mixing ratio (kg/kg)
- N        :: FT   : number concentration (#/m^3)
- mu       :: FT   : gamma distribution shape parameter μ
- r_th     :: FT   : threshold radius (m)
- chi_m    :: FT   : mass prefactor χ_m
- m0       :: FT   : reference particle mass m_0
- r0       :: FT   : reference radius r_0
- me       :: FT   : mass exponent m_e
- Δm       :: FT   : mass exponent increment Δm
- rho_a    :: FT   : air density (kg/m^3)
Keyword
---------
- Dmax     :: FT = Inf : optional upper cutoff for gamma distribution

Returns
-------
- dqdt_cross :: FT : mixing ratio tendency (kg/kg/s) crossing r_th,
                     proportional to supersaturation S if unknown

Notes
-----
- If you know supersaturation S, multiply dqdt_cross by S.
- n(r) is normalized to N.
"""
function my_conv_q_ice_to_q_sno_at_r_is_given_qi_tendency_sub_dep_old(
    param_set::APS,
    qi_tendency_sub_dep::FT,
    q::FT,
    N::FT,
    ρ_a::FT;
    Dmax::FT = FT(Inf),
    μ::FT = FT(NaN)
) where {FT}

    # condensed mass per unit volume

    acnv_rate = FT(0)

    if (q > FT(0)) && (qi_tendency_sub_dep > FT(0))
        microphys_params = TCP.microphysics_params(param_set)
        if isnan(μ)
            μ = μ_from_qN(param_set, ice_type, q, N; ρ=ρ_a)
        end

        # get the number distribution parameters
        n0, λ = get_n0_lambda(microphys_params, ice_type, q, ρ_a, N, μ; Dmax=Dmax)

        
        _r0::FT = r0(microphys_params, ice_type)
        _m0::FT = m0(microphys_params, ice_type)
        _me::FT = me(microphys_params, ice_type)
        _Δm::FT = Δm(microphys_params, ice_type)
        _χm::FT = χm(microphys_params, ice_type)

        r_is::FT = CMP.r_ice_snow(microphys_params)

        # particle mass at threshold
        m_r_is = _χm * _m0 * (r_is/_r0)^(_me + _Δm)

        # number at threshold
        n_r_is = n0 * r_is^μ * exp(-λ*r_is)

        # total condensed mass per volume
        M_cond = q * ρ_a

        # crossing flux
        acnv_rate = qi_tendency_sub_dep * (n_r_is * m_r_is * r_is) / M_cond
    end

    return acnv_rate
end


"""
Compute the mixing-ratio tendency (dq/dt, kg/kg/s) of ice mass crossing the
threshold radius r_is given a known bulk deposition tendency qi_tendency_sub_dep.

Implements the diffusion-limited growth derivation for a gamma distribution,
where the crossing fraction is

    frac_cross = (1/p) * (ξ^(μ+2) * exp(-ξ)) / Γ(μ+2),

with ξ = λ * r_is and p = m_e + Δm. The returned rate is

    dqdt_cross = qi_tendency_sub_dep * frac_cross

Arguments
---------
- param_set :: APS
- qi_tendency_sub_dep :: FT   : bulk mixing-ratio tendency (kg/kg/s)
- q  :: FT                    : condensed mixing ratio (kg/kg)
- N  :: FT                    : number concentration (#/m^3)
- ρ_a:: FT                    : air density (kg/m^3)

Keyword Arguments
-----------------
- Dmax :: FT = Inf            : optional cutoff
- μ    :: FT = NaN            : optional gamma shape (computed if NaN)

Returns
-------
- acnv_rate :: FT : crossing mixing-ratio tendency (kg/kg/s)
"""
function my_conv_q_ice_to_q_sno_at_r_is_given_qi_tendency_sub_dep(
    param_set::APS,
    qi_tendency_sub_dep::FT,
    q::FT,
    N::FT,
    ρ_a::FT;
    Dmax::FT = FT(Inf),
    μ::FT = FT(NaN)
) where {FT}

    acnv_rate = zero(FT)

    if (q > zero(FT)) && (qi_tendency_sub_dep > zero(FT))

        microphys_params = TCP.microphysics_params(param_set)

        # shape parameter
        if isnan(μ)
            μ = μ_from_qN(param_set, ice_type, q, N; ρ=ρ_a)
        end

        # get gamma distribution parameters
        _, λ = get_n0_lambda(microphys_params, ice_type, q, ρ_a, N, μ; Dmax=Dmax)

        # mass–radius law parameters
        _me = me(microphys_params, ice_type)
        _Δm = Δm(microphys_params, ice_type)
        p = _me + _Δm

        # threshold radius
        r_is = CMP.r_ice_snow(microphys_params)
        ξ = λ * r_is

        # closed-form crossing fraction
        frac_cross = (one(FT) / p) * (ξ^(μ + 2) * exp(-ξ)) / CM1.SF.gamma(μ + 2)

        # crossing tendency
        acnv_rate = qi_tendency_sub_dep * frac_cross
    end

    return resolve_nan(acnv_rate, zero(FT)) # if we get a NaN, just return 0... (n0 and lambda can blow up for very small q which map to no conv)
end



function invert__my_conv_q_ice_to_q_sno_at_r_is_given_qi_tendency_sub_dep__to_N(
    param_set::APS,
    qi_tendency_sub_dep::FT,
    acnv_rate::FT,
    q::FT,
    ρ_a::FT;
    Dmax::FT = FT(Inf),
    μ::FT = FT(0),
    branch::Int = 0
) where {FT}

    error("I'm not sure if this form is completely up to date -- see the python version for fully correct form - it's not really usable though except for exact inversion because the range of valid acnv_rates is so small")

    N = FT(NaN)

    if (q > zero(FT)) && (qi_tendency_sub_dep > zero(FT)) && (acnv_rate > zero(FT))

        microphys_params = TCP.microphysics_params(param_set)

        # shape parameter must be zero here
        @assert μ == 0 "This closed-form solver is valid only for μ = 0."

        # --- mass–radius parameters
        _me  = me(microphys_params, ice_type)
        _Δm  = Δm(microphys_params, ice_type)
        p    = _me + _Δm

        # --- observed crossing fraction
        f_obs = acnv_rate / qi_tendency_sub_dep

        # Existence check: max possible is 4/e^2
        if f_obs * p > 4 / exp(2)
            return FT(NaN)
        end

        # --- coefficients from microphys_params
        _m0::FT  = m0(microphys_params, ice_type)
        _χm::FT  = χm(microphys_params, ice_type)
        _r0::FT  = r0(microphys_params, ice_type)

        # Prefactor C
        C = _m0 * _χm * _r0^(-p) * CM1.SF.gamma(p + 1)

        # Threshold radius
        r_is = CMP.r_ice_snow(microphys_params)

        # Constant A
        A = r_is * (C / (ρ_a * q))^(one(FT) / p)

        # Lambert W argument
        s = -FT(0.5) * sqrt(f_obs * p)

        # Evaluate Lambert W on chosen branch
        W = LambertW.lambertw(branch, s)

        # Solve for N
        z = -2 * W
        y = z / A
        if y > zero(FT)
            N = y^p
        end
    end

    return N
end
















"""
Separate fcn for type stability.
    As long as we employ a separate acnv threshold, knowing what that thresh is is important for the limiter.
    Otherwise you can get into a scenario where if τ is very fast, you build up q and then suddenly demolish it all at once because of a fast τ and long enough timestep... this looks like oscillatory `pistoning` and is bad.
   
    You could try to separate out the part of the threshold that is q based even if we have the fallback q_ice_threshold...
    I guess you could argue we shouldn't limit a q based threshold bc it will just decay to 0, but really it's an exponential decay to 0 so it feels like we should stick to our `you can't do everything in one timestep` rule.
    and if the N is from T we really shouldn't go below it in 1 timestep ever.


    To make it simpler and reduce code duplication, although it's slower, we'll make this the main fcn and have my_conv_q_ice_to_q_sno_thresh() just only return the source
"""
function get_thresh_and_my_conv_q_ice_to_q_sno_thresh(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    qi::FT,
    p::FT;
    N::FT = FT(NaN),
    τ::FT = FT(NaN),  # consider making this just one timestep?
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    # v_ice::FT = FT(0),
    # v_snow::FT = FT(0),
    ice_acnv_power::FT = FT(2),
    r_acnv_scaling_factor::FT = FT(1),
    w::FT = FT(0), # this is the updraft speed, if we have it, we can use it to adjust the source term
    dwdz::FT = FT(0), # this is the sedimentation rate of ice, if we have it, we can use it to adjust the source term
    use_sed::Bool = false, # if true, we use dqidt_sed to adjust the source term, otherwise we just use the source term as is
    ice_type::CMT.IceType = ice_type,
) where {FT}

    # if r_acnv_scaling_factor = 1, we will just get q_thresh = q at most, in any case where N is calculated from q (adjust_ice_N = true), unless q_ice_threshold is not 0.
    
    microphys_params::ACMP = TCP.microphysics_params(param_set)

    if isnan(τ)
        τ = CMP.τ_acnv_sno(microphys_params)
    end

    S_qs::FT = FT(0)
    q_thresh::FT = FT(0) # We use this in limiters. if q = 0, S_qi is 0, so it doesn't matter too much but let's just also say 0.

    if qi > FT(0)
        
        #=
        Integrating from r_ice_snow to Dmax w/ int_nm_dr() doesnt work bc we cant actually create the distribution with q and N below r_ice_snow, only going to infinity (bc incomplete gamma fcns) so we end up getting wild guesses for the distributions that still mostly seem to scale w/ q...
        
        
        Instead we do something like the threshold (q-q_thresh) / τ, but to assess what is crossing the size thresh we scale it by (r/r_is)**3 (maybe the power could be calibrateable but i think it's 2 for area for coll and one for itself being of size to cross r_is)
        
        Note that if we don't/can't predict N, we just assume we can either:
            - assume we're always at r_is and just get 0 agg driven acnv (or use r_acnv scaled down from r_is))
            - assume a fixed thresh.
            The fixed thresh seems a bit subpar but beats adding another parameter I suppose.
        =#
        if isnan(N) 
            # we could also use r_acnv here w/
            # [ assume_N_is means we're always at r_is which means q_thresh is just q but we use q_acnv here... (maybe look at switching that idk, might be out of date now) ]
            # [ not assuming N_is means q_thresh is just q_ice_threshold ]
            q_thresh = get_q_threshold_acnv(param_set, ice_type; N=N, assume_N_is = false, r_acnv_scaling_factor = r_acnv_scaling_factor) 
            S_qs = my_conv_q_ice_to_q_sno_no_supersat_simple(param_set, qi, q_thresh, τ) # this is the rate of change of q_ice with respect to N, so we need to divide by τ to get the rate of change of N
        else
            q_thresh = get_q_threshold_acnv(param_set, ice_type; N=N, assume_N_is = false, r_acnv_scaling_factor = r_acnv_scaling_factor)
            S_qs = my_conv_q_ice_to_q_sno_no_supersat_simple(param_set, qi, q_thresh, τ) # this is the rate of change of q_ice with respect to N, so we need to divide by τ to get the rate of change of N
        end

        q_0 = FT(1e-7)

        if !use_sed
            if S_qs > FT(0) # avoid 0^neg power = Inf... (S_qs is 0 or positive)
                S_qs *= (S_qs * τ / q_0)^(ice_acnv_power-1) # we want q_0 ((q-q_thresh)/q_0)^p_acnv / τ and we have S_qs = (q-q_thresh)/τ . (S_qs * τ / q_0) ^ (p_acnv-1) = ((q - q_thresh)/q_0)^(p_acnv-1) so we'll multiply and get  (q - q_thresh)^(p_acnv) / (q_0^(p_acnv-1) τ) = q_0 * ((q - q_thresh)/q_0)^(p_acnv) / τ
            end
        else
            # here, the rate comes from sedimentation collisions and the size of the underlying particle. We will take the rate as ((q-q_thresh)/q_0)^p_acnv
            K = max(dwdz, FT(0)) + w/10 # max (dwdz,0) bc positive means more collisions... some background from w to smooth things in.
            S_qs = qi * K / τ * (qi/q_thresh)^(ice_acnv_power) # we want to scale the sedimentation rate by the size of the particles crossing the threshold, so we use dqidt_sed as the rate and scale it by (dqidt_sed * τ / q_0)^(ice_acnv_power-1)
        end
    end


    return S_qs, q_thresh # return the source term and the threshold
end

get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set::APS, q::TD.PhasePartition{FT}; N::FT = FT(NaN), τ::FT = FT(NaN), Dmin::FT = FT(0), Dmax::FT = FT(Inf), ice_acnv_power::FT = FT(2), r_acnv_scaling_factor::FT = FT(1), w::FT = FT(0), dwdz::FT = FT(0), use_sed::Bool = false) where {FT} =
    get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set, q.ice; N=N, τ=τ, Dmin=Dmin, Dmax=Dmax, ice_acnv_power=ice_acnv_power, r_acnv_scaling_factor=r_acnv_scaling_factor, w=w, dwdz=dwdz, use_sed=use_sed)


"""
    To make it simpler and reduce code duplication, although it's slower, we'll make this the secondary fcn and only return the source
    
    Having duplicate code in my_conv_q_ice_to_q_sno_thresh() and get_thresh_and_my_conv_q_ice_to_q_sno_thresh() would increase the risk of making a mistake trying to keep them in sync, especially if we make changes in my_conv_q_ice_to_q_sno_no_supersat() or get_q_threshold_acnv()
"""
function my_conv_q_ice_to_q_sno_thresh(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    qi::FT,
    N::FT = FT(NaN),
    τ::FT = FT(NaN),  # consider making this just one timestep?
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    # v_ice::FT = FT(0),
    # v_snow::FT = FT(0),
    ice_acnv_power::FT = FT(2),
    r_acnv_scaling_factor::FT = FT(1),
) where {FT}

    S_qs, _ = get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set, qi; N=N, τ=τ, Dmin=Dmin, Dmax=Dmax, ice_acnv_power = ice_acnv_power, r_acnv_scaling_factor = r_acnv_scaling_factor)
    return S_qs
end

my_conv_q_ice_to_q_sno_thresh(param_set::APS, q::TD.PhasePartition{FT}; N::FT = FT(NaN), τ::FT = FT(NaN), Dmin::FT = FT(0), Dmax::FT = FT(Inf), ice_acnv_power::FT = FT(2), r_acnv_scaling_factor::FT = FT(1)) where {FT} =
    my_conv_q_ice_to_q_sno_thresh(param_set, q.ice; N=N, τ=τ, Dmin=Dmin, Dmax=Dmax, ice_acnv_power=ice_acnv_power, r_acnv_scaling_factor=r_acnv_scaling_factor)


radius_scaling_factor_from_mass_scaling_factor(χm::FT) where {FT} = inv(cbrt(χm)) # 1/f_r^3 = χm --> f_r = 1 / cbrt(χm) # this is the factor by which we scale the mean radius to get the mean radius for the ice crystals, so that we can use it in the N_i and N_l calculations
mass_scaling_factor_from_radius_scaling_factor(f_r::FT) where {FT} =  inv(f_r^3) # χm = 1/f_r^3


"""
Get an estimate of Gamma shape parameter μ from q and N


We have a subsaturation booster... if we are subsaturated..., we assume r is larger because the small droplets are gone.

<r> = <r>_0 * f, where f = α[1-(1-S)]^β

Chat GPT wants H = max(0, 1-S), f = 1 + αH^β
"""
function μ_from_qN(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.IceType,
    q_i::FT,
    N_i::FT;
    ρ::FT = FT(1),
    ice_type::CMT.IceType = ice_type,
    max_μ::FT = FT(3.0), # this is the maximum μ we allow, if we get a larger one, we just return max_μ
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    # S::FT = FT(1) # supersaturation boost factor, default to 1 (no boost)
) where {FT}

    if !isfinite(N_i) || !isfinite(q_i)
        return FT(0)
    end

    return FT(0)  # I think we could try to move towards mu > 0 as <r> grows, but it's too hard to do in a principled way. I think we'd also need N_INP etc... idk

    # new_N_i = <adjust ice N somehow>
    # # compare monodisperse radii to get what we think the true μ is
    # r_i = ri_from_qN(param_set, q_i, N_i; ice_type = ice_type, monodisperse = true, Dmin=Dmin, Dmax=Dmax)
    # r_new = ri_from_qN(param_set, q_i, new_N_i; ice_type = ice_type, monodisperse = true, Dmin=Dmin, Dmax=Dmax)

    # if r_new > r_i # probably means one had N=0 or something and returned min_radius instead of infty ... just go with μ = 0 (you could also make an argument in that case to go to μ_max, but μ=0 is a more pure sol'n, especially for very small q. )
    #     return FT(0)
    #     # error("Got r_new = $r_new > r_i = $r_i; from inputs q_i = $q_i; N_i = $N_i; ρ = $ρ; new_N_i = $new_N_i; r_i = $r_i; r_new = $r_new; ")
    # end
    # #=
    #     since N can only go up, then r_new would be smaller... not bigger.
    #     if the new N is much bigger (r much smaller), we want to be closer to monodisperse (μ = ∞)
    #     if the new N is unchanged, we want to be closer to the original state (μ = 0)
    # =#

    # μ = r_i / r_new - FT(1)
    # # μ = r_new / r_i - 1 # if the new r is much bigger, we're probably also getting more monodisperse, in the other limit, we go to 0 [[ this should be positve by construction, new_N_i can only go up ]]
    # μ = min(μ, max_μ) # I saw a μ = 66. and knew we had to add this XD, thought i'd done seen it all

end

# For liquid this should be something like 10 really but... we'll just assume monodisperse I guess...
function μ_from_qN(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.LiquidType,
    q_l::FT,
    N_l::FT;
    ρ::FT = FT(1),
    ice_type::CMT.IceType = ice_type,
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    ηmax::FT = FT(0.577), # from MG08
) where {FT}

    if !isfinite(N_l)
        return FT(0)
    end

    #=
        Morrison 2005 suggests
         μ = p_c = \frac{1}{m_l} \frac{L_v q_c}{\left(g-\gamma_s c_p\right)}

         and they use mixing_legth = 30m
         We have mixing length but I don't wanna pass it in here, too complex

         We would also need to calculate the lapse rate, etc...
    =#

    η = FT(0.0005714) * N_l/FT(1e6) + FT(0.2714) # from Morrison & Gettelman, 2008 derived from Martin et al, 1994
    η = min(η, ηmax)
    μ =  1/η^2 -1
    return clamp(μ, FT(0), FT(10)) #  Would probably be fine just running 0s but it should be ok...

    # return FT(0)
end


function μ_from_qN(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.SnowType,
    q_i::FT, 
    N_i::FT;
    ρ::FT = FT(1),
    ice_type::CMT.IceType = ice_type,
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT}
    # return FT(NaN) # for type completenes, but these are not defined
    # return FT(0) # I think 0 is the same as not defined bc it adds nothing... snow has it's own mu thingy wrapped up in n0 or m0, afaik rain has nothing
    
    # NOTE: They explicitly use Mashall-Palmer, so μ = 0
    # It's unfortunate they used the name CMP.μ_sno(TCP.microphysics_params(param_set)), it is NOT size disribution shape parameter μ, just a coefficient to set the intercept...
    return FT(0) # 
end

function μ_from_qN(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMT.RainType,
    q_l::FT, 
    N_l::FT;
    ρ::FT = FT(1),
    ice_type::CMT.IceType = ice_type,
    liq_type::CMT.LiquidType = liq_type,
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT}
    # return FT(NaN) # for type completenes, but these are not defined
    # return FT(0) # I think 0 is the same as not defined bc it adds nothing... snow has it's own mu thingy wrapped up in n0 or m0, afaik rain has nothing
    return μ_from_qN(param_set, liq_type, q_l, N_l; ρ=ρ, ice_type=ice_type, Dmin=Dmin, Dmax=Dmax) # we could defer to the liquid one but it does crazy things... idk why...
end


function get_T_top_from_N_i(
    param_set::APS,
    relaxation_timescale::AbstractRelaxationTimescaleType,
    q_i::FT,
    N_i::FT;
    ρ::FT = FT(1),
    ice_type::CMT.IceType = ice_type,
) where {FT}
    return FT(NaN)
end

function get_T_top_from_N_i(
    param_set::APS,
    relaxation_timescale::INP_Aware_Timescale,
    q_i::FT;
    ρ::FT = FT(1),
    ice_type::CMT.IceType = ice_type,
    μ::FT = FT(0),
    monodisperse::Bool = false,
) where {FT}

    if iszero(q_i)
        error("q_i is 0")
    end

    microphys_params = TCP.microphysics_params(param_set)

    # r_acnv_scaling_factor = param_set.user_params.r_ice_acnv_scaling_factor # this MUST be less than 1!!!
    # _q_acnv_0_ = q_acnv_0(param_set, ice_type, r_acnv_scaling_factor) # this is the volume of a single ice crystal at the acnv radius, so N * q_acnv_0 = q_threshold
    
    r_is = CMP.r_ice_snow(microphys_params)
    _q_is = q_is(param_set, ice_type) # this is the volume of a single ice crystal at the r_is, so N * q_is = q_threshold

    if monodisperse
        N_top = q_i * ρ / _q_is
    else
        λ = (μ + FT(1)) / r_is # this is the scaling factor for the radius, so we can use it to get n0
        _χm::FT = get_χm(param_set, ice_type) # this is the mass scaling factor, so we can use it to get n0
        _m0::FT = m0(microphys_params, ice_type) # this is the mass of the ice crystal at the acnv radius
        _r0::FT = r0(microphys_params, ice_type) # this is the radius of the ice crystal at the acnv radius
        _me::FT = me(microphys_params, ice_type) # this is the exponent for the mass of the ice crystal
        _Δm::FT = Δm(microphys_params, ice_type) # this is the exponent for the mass of the ice crystal
        N_top = (q_i * ρ) / (_χm * _m0 * _r0^(-(_me + _Δm)) * (CM1.SF.gamma(μ + _me + _Δm + FT(1)) / CM1.SF.gamma(μ + FT(1))) * λ^(-(_me + _Δm))) # this is the number concentration of ice crystals at the acnv radius
    end
    # Now we would just need to get T_top. For the exponential and powerlaw T-scalings it's easy... otherwise it's hard.
    T = get_T_from_N_i(param_set, relaxation_timescale, N_top)

    return isnan(T) ? get_T_top_from_N_i_Cooper_curve(N_top) : T
end

get_T_top_from_N_i_Cooper_curve(N::FT) where {FT} =  FT(273.15) - (FT(1.0)/FT(0.304)) * log(N / (FT(0.005)*FT(1000.0))) # technically cooper is clamped but idk that it matters so much...
function get_N_i_Cooper_curve(T::FT; clamp_N::Bool = false) where{FT} # see https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L2348, they also have fletcher, we should prolly just pick one...
    N = FT(0.005) * exp(FT(0.304) * (FT(273.15) - T)) * FT(1000.0)  # L⁻¹ → m⁻³
    if clamp_N
        N = clamp(N, FT(0), FT(500e3)) # why did we divide by ρ here? (maybe it's a density depending concentration? idk.. this way goes up at lower ρ though...)
    end
    return N
end

"""
Ignores any contribution from dρ/dz...
"""
function get_d_N_i_Cooper_curve_dT(T::FT) where {FT}
    N = get_N_i_Cooper_curve(T)
    return -FT(0.304) * N
end
 get_d_N_i_Cooper_curve_dz(T::FT, dTdz::FT) where {FT} = get_d_N_i_Cooper_curve_dT(T) * dTdz

function get_T_from_N_i(param_set::APS, relaxation_timescale::Union{ExponentialTScalingIceRelaxationTimescale, ExponentialTScalingIceRawRelaxationTimescale, GeometricLiqExponentialTScalingIceRelaxationTimescale}, N::FT) where {FT}
    if iszero(N)
        return FT(NaN)
    end
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    # N = FT(c_1 * exp(c_2 * (T - T_fr))), so we invert this to solve for T
    return T_fr + (1/c_2) * log(N/c_1)
end
function get_T_from_N_i(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, N::FT) where {FT}
    if iszero(N)
        return FT(NaN)
    end
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_4i
    c_2 = relaxation_timescale.c_5i
    # N = FT(c_1 * exp(c_2 * (T - T_fr))), so we invert this to solve for T
    return T_fr + (1/c_2) * log(N/c_1)
end

function get_T_from_N_i(param_set::APS, relaxation_timescale::Union{PowerlawTScalingIceRelaxationTimescale, GeometricLiqPowerlawTScalingIceRelaxationTimescale}, N::FT) where {FT}
    if iszero(N)
        return FT(NaN)
    end
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    # N_i = (10^c_1) * (-(T - T_fr))^c_2 so invert this to get T
    return T_fr - (10^c_1 / N)^(1/c_2)
end







"""
    N_from_qr(param_set, q_type, q, r; monodisperse=true, μ=0, ρ=1,
              Dmin=0, Dmax=Inf, add_dry_aerosol_mass=false)

Compute number concentration `N` from a mass mixing ratio `q` and scale radius `r`
for a gamma particle size distribution with shape `μ` and scale `λ = (μ+1)/r`.

# Assumptions / notation
- Number distribution (full-domain, r ≥ 0):

        n(r) = N * λ^(μ+1) / Γ(μ+1) * r^μ * exp(-λ r)

  where `N` is the total number concentration when the domain is [0, ∞).

- Particle mass law:

        m(r) = _χm * _m0 * _r0^(-s_m) * r^(s_m)

  where `s_m = _me + _Δm`.

- Total mass per volume:

        M = q * ρ

  (mixing ratio times air density).

# Derivation (no truncation, r_min = 0)
1. Number normalization (r ∈ [0, ∞)):

        ∫ n(r) dr = N

   by construction.

2. Mass integral:

        M = ∫ m(r) n(r) dr
          = (_χm * _m0 * _r0^(-s_m)) * N * λ^(μ+1)/Γ(μ+1) * ∫_0^∞ r^(μ+s_m) e^(-λ r) dr

   Using

        ∫_0^∞ r^α e^(-λ r) dr = Γ(α+1) / λ^(α+1),   α = μ + s_m

   we get

        M = (_χm * _m0 * _r0^(-s_m)) * N * (Γ(μ+s_m+1)/Γ(μ+1)) * λ^(-s_m)

3. Solve for N:

        N = (q * ρ) / (_χm * _m0 * _r0^(-s_m) * (Γ(μ+s_m+1)/Γ(μ+1)) * λ^(-s_m))

# Truncated lower cutoff (r ∈ [r_min, ∞))
Let `s1 = μ+1`, `s2 = μ+s_m+1`, `x = λ * r_min`.

- Unnormalized integrals:

        ∫_{r_min}^∞ r^(s1-1) e^(-λ r) dr = Γ(s1, x) / λ^s1
        ∫_{r_min}^∞ r^(s2-1) e^(-λ r) dr = Γ(s2, x) / λ^s2

- Normalization requires using the upper regularized gamma `Q(s,x) = Γ(s,x)/Γ(s)`.

  Carrying through the algebra gives:

        M = (_χm * _m0 * _r0^(-s_m)) * N * (Γ(s2)/Γ(s1)) * (Q2/Q1) * λ^(-s_m)

- Solving for N:

        N = (q * ρ) / (_χm * _m0 * _r0^(-s_m) * (Γ(s2)/Γ(s1)) * (Q2/Q1) * λ^(-s_m))

Numerically, this is evaluated with `loggamma` to avoid overflow/underflow.
"""
function N_from_qr(param_set::APS, q_type::CMTWaterTypes, q::FT, r::FT; monodisperse::Bool=true, μ::FT=FT(0), ρ::FT=FT(1), Dmin::FT=FT(0), Dmax::FT=FT(Inf), add_dry_aerosol_mass::Bool=false) where {FT}
    microphys_params = TCP.microphysics_params(param_set)
    if monodisperse
        # For monodisperse case, we can directly use the formulas
        N = q / particle_mass(param_set, q_type, r)
    else

        if (iszero(Dmin) && isinf(Dmax))
            λ = (μ + FT(1)) / r # this is the scaling factor for the radius, so we can use it to get n0
            _χm::FT = get_χm(param_set, ice_type) # this is the mass scaling factor, so we can use it to get n0
            _m0::FT = m0(microphys_params, ice_type) # this is the mass of the ice crystal at the acnv radius
            _r0::FT = r0(microphys_params, ice_type) # this is the radius of the ice crystal at the acnv radius
            _me::FT = me(microphys_params, ice_type) # this is the exponent for the mass of the ice crystal
            _Δm::FT = Δm(microphys_params, ice_type) # this is the exponent for the mass of the ice crystal

            Dmin = (add_dry_aerosol_mass ? max(param_set.user_params.particle_min_radius, Dmin) : Dmin)
            if iszero(Dmin)
                # N = (q * ρ) / (_χm * _m0 * _r0^(-(_me + _Δm)) * (CM1.SF.gamma(μ + _me + _Δm + FT(1)) / CM1.SF.gamma(μ + FT(1))) * λ^(-(_me + _Δm))) # [works!] this is the number concentration of ice crystals at the acnv radius
                γratio = exp(CM1.SF.loggamma(μ + _me + _Δm + 1) - CM1.SF.loggamma(μ + 1)) # might be more numerically stable
                # N = (q * ρ) / (_χm * _m0 * _r0^(-(_me + _Δm)) * γratio * λ^(-(_me + _Δm)))
            else
                # -------- with aerosol core --------
                # _, Q1 = CM1.SF.gamma_inc(μ + 1, λ * Dmin)
                # _, Q2 = CM1.SF.gamma_inc(μ + _me + _Δm + 1, λ * Dmin)
                # γratio = exp(CM1.SF.loggamma(μ + 1) - CM1.SF.loggamma(μ + _me + _Δm + 1)) * (Q1 / Q2) # might be more numerically stable
                γratio = exp(CM1.SF.loggamma(μ + _me + _Δm + 1, λ * Dmin) - CM1.SF.loggamma(μ + 1, λ * Dmin))
            end
            N = (q * ρ) / (_χm * _m0 * _r0^(-(_me + _Δm)) * γratio * λ^(-(_me + _Δm)))
        else
            error("N_from_qr not implemented for non-monodisperse with Dmin/Dmax limits")
        end
    end

    if !isfinite(N)
        N = eps(FT)
    end

    return N
end


"
This much times N_below_r_is is how much we should expect to see added... N/N_below_r_is = f_boost., upto the mass flux factor.

For simplicity we always apply the boost, so N_i w/o boost = N_i / f_boost
"
function f_boost_MF(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    dNINP_dz::FT,
    ;
    ρ::FT = FT(1.0), # density of the air, if we have it, we can use it to adjust the source term
    S_i::FT = FT(0), # supersaturation
    massflux::FT = FT(0), # if we have a massflux, we can use it to adjust the source term, (combines area and velocity, velocity is prolly more important for penetration and downdraft generation but area matters.
    NINP_top_over_N_INP::FT = FT(1),
    massflux_min::FT = FT(0.05),
    massflux_0::FT = FT(0.05),
) where {FT}
    

    #=
        Even in the cases w/ no massflux, N_i+N_s (N_INP) still diffuses downwards by sedimentation. Sure we can get a higher boost N_i instead of N_s due to massflux but we still get an overall boost.
        Now that we apply the boost and then split N_s off, we should force some minimum level of this effect always.
        Otherwise, we risk N_i falling to unrealistically low levels (essentially back to NINP, but then still having N_s taken out)
    =# 
    massflux = max(massflux, massflux_min * ρ)


    if isfinite(NINP_top_over_N_INP)
        # Only boost if NINP < (0.75 NINP_top)  =>  dNINP_dz < (0.75 * dNINP_dz_top)  =>  inv_dINP_dz > (1/.75) * inv_dINP_dz_top
        dINP_dz_top = dNINP_dz * NINP_top_over_N_INP
        # min_inv_dINP_dz = (1/.75) * (1/dINP_dz_top)
    else
        # min_inv_dINP_dz = one(FT) # a classic cloud top value?
        dINP_dz_top = dNINP_dz # ? why
    end

       
    c_MF::FT = param_set.user_params.massflux_N_i_boost_factor
    dNINP_dz = max(dNINP_dz, eps(FT)) # avoid div by 0


    α = (massflux/ρ) / massflux_0 # normalize to saturate at  massflux_0 MF, yielding c_MF boost.
    r = (dINP_dz_top) / max(abs(dNINP_dz), eps(FT))      # ratio = N_top / N
    # target = max(c_MF * clamp((0.6*dINP_dz_top)/(dNINP_dz), one(FT), FT(3)), one(FT))
    # cap = max(FT(0.6) * (dINP_dz_top/dNINP_dz), one(FT))
    # f_boost = min(one(FT) + (target - one(FT))*α, cap)

    # let's only go halfway to top value, so we dont get stationary dN/dz. so instead of the full α * dINP_dz_top/dNINP_dz to go from N --> N_top, we want to go from N --> (N + αN_top)/2
    # Thus, the target is now (N + αN_top)/2 / N = 0.5 + 0.5 * (αN_top/N) = 0.5 + 0.5 * α(dINP_dz_top/dNINP_dz)

    # I feel like f_r should be bigger and f_max smaller?
    # f_r = FT(0.7) # the ratio of the max ratio we can boost to relative to the cloud top value.
    # f_max = FT(0.7) # how far we want to allow going towards the top value (rn is 0.6 N_top)
    f_r::FT = param_set.user_params.massflux_N_i_boost_max_ratio # the ratio of the max ratio we can boost to relative to the cloud top value.
    f_max::FT = param_set.user_params.massflux_N_i_boost_progress_fraction # how far we want to allow going towards the top value (rn is 0.6 N_top)
    target = max(one(FT), (one(FT)*(1-f_max) + f_max * f_r * clamp(r, one(FT), FT(6))))
    cap = max(one(FT), (one(FT)*(1-f_max) + f_max * f_r * r))
    f_boost = min(one(FT) + c_MF * α * (target - one(FT)), cap)


    # RH control, since it relies on moving down ice. note in paper variance should matter. Goes from fully off at S_i = -0.05 to fully on at S_i = 0.05
    # variance seems to be about 5%, so let's try off by S_i = -.05, fully on by 0.05
    f_boost = one(FT) + clamp(FT(0.5 + 10.0*S_i), FT(0.0), FT(1.0)) * (f_boost - one(FT))

    f_boost = max(f_boost, one(FT)) # never reduce N_i

    return f_boost
end


"""

    If subsaturated, we can assume none of the existing ice is being lost to growth based autoconversion (PITOSN threshold is still a possibility)
    
    If subsaturated, we don't really know what to do. We should assume we have large-ish particles so that sedimentatin is happening, and that we're not losing things to acnv by growth.
    So we will just err towards r_th... (maybe halfway between Ni_from_INP_qi_qs and N_th?)
"""
function get_Ni_from_INP_qi_qs(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    N_i::FT,
    N_INP::FT, # I think these should be equal?
    q_i::FT,
    q_s::FT;
    ρ::FT = FT(1.0), # density of the air, if we have it, we can use it to adjust the source term
    ice_type::CMT.IceType = ice_type,
    monodisperse::Bool = true, # if true, we assume all ice crystals
    μ::FT = FT(NaN),
    add_dry_aerosol_mass::Bool = false,
) where {FT}

    if monodisperse
        error("monodisperse not implemented")
    else

        if q_i > FT(0)

            r_min::FT = param_set.user_params.particle_min_radius
            if add_dry_aerosol_mass
                q_i += mass(param_set, ice_type, r_min, N_INP; monodisperse = true) 
            end

            if isnan(μ)
                μ = μ_from_qN(param_set, ice_type, q_i, N_i; ρ=ρ, ice_type=ice_type, Dmin=0, Dmax=Inf)
            end
            
            r_i_acnv = r_ice_acnv(param_set, param_set.user_params.r_ice_acnv_scaling_factor) # this is the radius of the ice crystal at the acnv radius
            r_thresh = get_r_cond_precip(param_set, ice_type) * FT(param_set.user_params.r_ice_snow_threshold_scaling_factor)
            N_i_acnv = N_from_qr(param_set, ice_type, q_i, r_i_acnv; monodisperse = false, μ=μ, ρ=ρ) # this is the minimum N_i we should have, since we can't have particles smaller than r_acnv
            # N_thresh = N_from_qr(param_set, ice_type, q_i, r_thresh; monodisperse = false, μ=μ, ρ=ρ)

            if q_s > FT(0)

                q_here = min(q_i + q_s, q_from_rN(param_set, ice_type, r_thresh, N_INP; monodisperse = monodisperse, μ=μ, ρ=ρ)) # we don't want to go above the threshold, so we limit q here to the threshold value

                N_above_r_is = get_fraction_below_or_above_thresh_from_qN(param_set, ice_type, q_here, N_INP, ρ, r_thresh; return_N = true, below_thresh = false, return_fraction=false, μ=μ, Dmax=FT(Inf))
                NI_below_r_is = N_INP - N_above_r_is # use r_thresh vs r_is here...

                # NI_below_r_is = clamp(NI_below_r_is, N_thresh, N_i_acnv) # we can't go below the acnv number, or above the INP number
                NI_below_r_is = min(NI_below_r_is, N_i_acnv) # Force q to at least mean N_i_acnv exists, but don't enforce threshold radius.

                # ensure N_s supports the N we predict
                r_i = r_from_qN(param_set, ice_type, q_i, NI_below_r_is; monodisperse = false, μ=μ, ρ=ρ) # rought r_i based on current N below r_is
                if !isfinite(r_i)
                    @warn "Got non-finite r_i = $r_i from q_i = $q_i; NI_below_r_is = $NI_below_r_is; N_INP = $N_INP; μ = $μ; ρ = $ρ"
                end
                # r_s_i = r_thresh # snow should be at least r_thresh [[ the radius snow would have if it was still just ice]]
                r_s_i = clamp(2.5 * r_i, r_thresh, 3r_thresh) # [[ empirical fit from data in RF09 , goes from about 2x at cloud top to 3x at bottom... but then r_i disappears ]]
                Ns_est = N_from_qr(param_set, ice_type, q_s, r_s_i; monodisperse = false, μ=μ, ρ=ρ) # this is how many crystals we would expect to see in the snow
                NI_below_r_is += max((N_INP - NI_below_r_is) - Ns_est, FT(0))

                # if rand() < 1e-6
                #     @warn "I feel like as qs grows, we prevent qi growth by forcing N to go over to snow rather than allowing r to grow.. i'm not sure which param is controlling this... maybe we can change it."
                #     @warn "We also want to allow larger <r>, not more <r>... reducing N_s will increase N... but still..."
                #     @warn "Maybe our acnv/production ratio peaks at too low an <r>? but we have calibrated params for that."
                #     @warn "Our wi also maybe is peaking at too low an <r>? idk..."
                #     @warn "Or it could have just been our too-low midddle from using no_boost... meant that not enoguh production reached middle levels but sed had to remain high. also we didnt have vert_adv right before so maybe that artificially amde no boost tau look more optimal by reducing growth with height..."
                # end

                return NI_below_r_is

            else
                # return N_i # no snow, so everything is ice
                # return clamp(N_INP, N_i_acnv, N_thresh) # we can't go below the acnv number, or above the INP number
                return min(N_INP, N_i_acnv) # Force q to at least mean N_i_acnv exists, but don't enforce threshold radius.
            end
            
        else # no ice, so really should just be 0. Leave unchanged in case we do that limiting elsewhere? idk
            # return N_i
            return FT(0)
        end
    end
end

"""

We have a subsaturation booster... if we are subsaturated..., we assume r is larger because the small droplets are gone.

<r> = <r>_0 * f, where f = α[1-(1-S)]^β

Chat GPT wants H = max(0, 1-S), f = 1 + αH^β
"""
function adjust_ice_N(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    N_i::FT, 
    N_INP::FT,
    q_i::FT,
    ;
    ρ::FT = FT(1.0), # density of the air, if we have it, we can use it to adjust the source term
    S_i::FT = FT(1), # supersaturation
    # S_l::FT = FT(1), # maybe use q_l instead bc in steady state this is right at 0...
    q_l::FT = FT(0),
    q_s::FT = FT(0),
    monodisperse::Bool = true, # if true, we assume all ice crystals are the same size, otherwise we assume a distribution of sizes
    ice_type::CMT.IceType = ice_type,
    decrease_N_if_subsaturated::Bool = false,
    N_INP_top::FT = FT(NaN), # we can't go past this value I guess?
    massflux::FT = FT(0), # if we have a massflux, we can use it to adjust the source term, (combines area and velocity, velocity is prolly more important for penetration and downdraft generation but area matters.
    dNINP_dz::FT = FT(0), # if we have a dINP/dz, we can use it to adjust the source term, negative means more INP below us
    w_i::FT = FT(0), # this is the terminal velocity, from the last iter though...
    N_i_from_INP::Bool = false, # if true, N_i was derived from N_INP and should be adjusted accordingly,
    apply_massflux_boost::Bool = false, # turn off in updraft
    apply_sedimentation_boost::Bool = false, # seems to be too strong and create waves of downward moving ice pulses, i think it should probably not exceed 0.1, since anything falling fast should be snow...
) where {FT}

    if (q_i > FT(0)) || (N_i > FT(0))

        if isinf(N_i)
            N_i = floatmax(FT)
        end

        if iszero(q_i)
            return FT(0) # no ice, no N_i
        end

        N_i_in = N_i # save input value for debugging

        microphys_params::ACMP = TCP.microphysics_params(param_set)
        r_acnv_scaling_factor::FT = param_set.user_params.r_ice_acnv_scaling_factor # this MUST be less than 1!!!
        r_thresh::FT = get_r_cond_precip(param_set, ice_type) * FT(param_set.user_params.r_ice_snow_threshold_scaling_factor)
        r_is = CMP.r_ice_snow(microphys_params)

        # r_thresh *= (80/70) # we really just need this to stop wi from growing too fast, but no reason we need to limit this ourselves right?... # since we aren't enforcing r_thresh at the lower end but are just letting N_INP rock, I think we're ok...
        # [[ if we do not expand this range, then by only ever setting r to maybe 5% over r_thresh, \tau_acnv_thresh can be slower... but keep in mind <r> gets diagnosed at the same time as the tendency so you don't wanna let it go too far overboard because the diagnosed value will stick since it's diagnosed and not prognosed. the LES looks like maybe 80 micron peak and activating around 70 micron.. this thresh scaling allows activation w/o having to boost N_i or anything ]]

        if monodisperse
            μ = FT(NaN) # not used in monodisperse case
            _q_acnv_0_ = q_acnv_0(param_set, ice_type, r_acnv_scaling_factor) # this is the volume of a single ice crystal at the acnv radius, so N * q_acnv_0 = q_threshold
            N_i_acnv = (q_i * ρ) / _q_acnv_0_ # make sure existing ice contributes to the size, particularly w/ sedimentation bringing in new ice
        
            q_thresh = q_from_r(param_set, ice_type, r_thresh)
            N_thresh = (q_i * ρ) / q_thresh
        else

            r_i_acnv = r_ice_acnv(param_set, FT(r_acnv_scaling_factor)) # this is the radius of the ice crystal at the acnv radius
            μ = μ_from_qN(param_set, ice_type, q_i, N_i; ρ=ρ) # this is the factor by which we scale the mean radius to get the mean radius for the ice crystals, so that we can use it in the N_i and N_l calculations
            N_i_acnv = N_from_qr(param_set, ice_type, q_i, r_i_acnv; monodisperse = false, μ=μ, ρ=ρ)

            # if decrease_N_if_subsaturated
            #     if S_i < FT(0)
            #         r_thresh = r_thresh + (r_thresh - 1.05 * r_thresh)/(0 - -0.05) * clamp(S_i, -0.05, 0) # Let r_thresh be 5\% bigger by S = -5%, meaning N_thresh can be inferred smaller and <r> above the true r_thresh [[ feels like we don't need this since N doesn't seem to shrink any faster in subsaturated conditions, if anything it's slower lol. might be tied to acnv and <r> ]]
            #         # decrease N_i if subsat (5% decrease at 15 percent subsat -- supposed to mean only large particles survive...)
            #         # N_i = N_i - (N_i - 0.95N_i) / (0 - -0.15) * clamp(S_i, -0.15, 0) # allow a slight decrease in N to account for sublimation, up to 5% at -5% subsaturation... (why would subsat matter? i guess a continuous fcn is good though)
            #     end
            # end
            N_thresh = N_from_qr(param_set, ice_type, q_i, r_thresh; monodisperse = false, μ=μ, ρ=ρ)
        end

        # if N_i < N_thresh # Boosts N_i if supersaturated [ we dont want to do this actually -- it looks like <r> grows if liquid is available ... also only doing this under N_thresh is a little weird]
            # N_i = N_i + (1.05N_i - N_i)/(0.2 - 0.) * clamp(S_i, 0, 0.2) # allow a slight boost to account for import from elsewhere, up to 5% at 0.2 supersaturation... (why would supersat matter? i guess a continuous fcn is good though)
        # end

        # if we have liquid, we seem to skew to higher <r>, which implies lower N... (not sure how the lower N makes sense, maybe it's a gamma artifact? accr pushes the entire dist right)
        # if  (q_l > FT(0)) # in steady state S_l might just be fixed near 0 due to our integration method... maybe we should use liquid presence but the scaling is hard idk...
        #     # decrease N_i, by 5% at 5 percent subsat (liq adjustment is very fast....)
        #     # N_i = N_i - (N_i - 0.95N_i) / (0 - -0.05) * clamp(S_l, -0.05, 0)
        #     # let's go for a 15% decrease by q_l = 2e-4
        #     N_i = N_i - (N_i - 0.85*N_i) / (FT(2e-4) - FT(0)) * clamp(q_l, FT(0), FT(2e-4))
        # end

        # Added 12/01/2025 -- I think we're are too susceptible to a low N_INP_top... We should never go below N_thresh/8 I think. If we have ice it came from somewhere. I think this also feeds into wi terminal velocity explosions and crashes.
        if !isnan(N_INP_top)
            N_INP_top = max(N_INP_top, N_thresh/8)
        end

    


        # ----------------------------------------------------------------- #
        # if param_set.user_params.apply_massflux_N_i_boost && (N_i > FT(1e-6))
        if (N_i > FT(1e-6)) # prevent near 0 problems



            # N_below_r_is = get_Ni_from_INP_qi_qs(param_set, N_i, N_INP, q_i, q_s; ρ=ρ, ice_type=ice_type, μ=μ, monodisperse=monodisperse, add_dry_aerosol_mass=true)
            # N_below_r_is = N_INP # try doing the boost to all INP, then figuring out what is ice and what is snow...
            
            #=
            
                When we are subsaturated, we should err towards being at r_thresh... since INP are no longer being lost to dep_acnv, we don't need to reduce as much
                At the same time, we know the particles are large enough to be sedimenting
                So, erring towards r_thresh is a good compromise.

                Note if we're below r_thresh (above N_thresh), we don't necessarily have to go down? But we should def go up? idk...
                Aim to hit r_th at -5% subsat

                However, just staying at r_th can lead to runaway q_growth... So we'll assume slightly above r_th so forced PITOSN acnv happens.
                Additionally, don't allow growth far beyond N_INP... wherever the particles came from should be nearby though it could be downdraft boosted etc. this just should stop runaway growth by further enforcing PITOSN and conversion to snow which is reasonable since theyre sedimenting particles anyway...
            =# 
            # if S_i < FT(0)
            #     N_thresh_adjusted = N_from_qr(param_set, ice_type, q_i, 1.1 * r_thresh; monodisperse = false, μ=μ, ρ=ρ) # aim for 10% above r_thresh so we can still have acnv
            #     if N_below_r_is < N_thresh_adjusted
            #         # Bring N up towards N_th if it's below it, since we should err towards r_thresh, do not reduce N if it's above N_th
            #         # N_below_r_is += (N_thresh_adjusted - N_below_r_is) * (-clamp(S_i, -0.05, 0) / 0.05)
            #         N_below_r_is += (N_thresh_adjusted - N_below_r_is) * (clamp(S_i, -0.05, 0) / 0.05 + 1)

            #         N_below_r_is = min(N_below_r_is, 2.5*N_INP) # don't let it exceed 2.5*INP
            #     else
            #         # do nothing, do not reduce, since we know subsat can lead to increased N_i due to reduced acnv. but we do not want to effect an increase
            #     end
            # end

            # if N_i_from_INP
            #     N_i = N_below_r_is # just set it directly, since we assume N_i was derived from N_INP
            # else
            #     N_i *= (N_below_r_is / N_INP)
            # end
            # N_INP = N_below_r_is


            if apply_massflux_boost
                massflux_N_i_boost_factor = f_boost_MF(param_set, dNINP_dz, ρ=ρ, S_i=S_i, massflux=massflux, NINP_top_over_N_INP = (N_INP_top/N_INP) ) # recalc in case we want to change it later

                N_INP_adjusted = N_INP * massflux_N_i_boost_factor
                if N_i_from_INP
                    N_i = N_INP_adjusted
                else # still try it
                    N_i *= massflux_N_i_boost_factor
                    # N_i = N_INP_adjusted
                end
            else
                N_INP_adjusted = N_INP
                N_i = min(N_i, N_INP) # just make sure we don't exceed INP if we aren't applying the boost
            end


            N_below_r_is = get_Ni_from_INP_qi_qs(param_set, N_i, N_INP_adjusted, q_i, q_s; ρ=ρ, ice_type=ice_type, μ=μ, monodisperse=monodisperse, add_dry_aerosol_mass=true)
            if N_i_from_INP
                N_i = N_below_r_is
            else
                N_i *= (N_below_r_is / N_INP_adjusted) # adjust N_i to match the new N_INP_adjusted
            end
            N_INP_adjusted = N_below_r_is


            # apply the subsat boost towards r_thresh again...
            #=
                The smaller N_INP_adjusted is, the further we should be from r_thresh, so we should err to higher r, since that's more time for PITOSN to have acted.
                # So we assume in the growth region,  we max out at r_thresh. below that, we very quickly growth to r_th and keep growing from there.
                # Really it is distance that matters, and not even the RH at all... since even 1% RH is huge relative to qi often... but this is an ok proxy.

                We should actually clamp ourselves to be between N_thresh and N_thresh_adjusted?
                At 0% subsat, we should be at N_INP_adjusted. Very rapidly this should go to N_thresh, and then grow up to N_thresh_adjusted
            =#
            if S_i < FT(0)
                # N_thresh_adjusted = N_from_qr(param_set, ice_type, q_i, 1.1 * r_thresh; monodisperse = false, μ=μ, ρ=ρ) # aim for 10% above r_thresh so we can still have acnv

                # N_thresh_adjusted should increase as N_INP_adjusted becomes smaller... let's have it asymptote to N_th as N_INP_adjusted --> N_th, and to r = 1.2 * r_th as N_INP_adjusted --> 0
                # r_thresh_adjusted = r_thresh + (1.2*r_thresh - r_thresh) * (1 - clamp(N_INP_adjusted / N_thresh, 0, 1)) # so at N_INP_adjusted = N_th, r_thresh_adjusted = r_thresh, at N_INP_adjusted = 0, r_thresh_adjusted = 1.2 * r_thresh
                r_thresh_adjusted = 1.08*r_thresh + (1.2*r_thresh - 1.08*r_thresh) * (1 - clamp(N_INP_adjusted / N_thresh, 0, 1)) # so at N_INP_adjusted = N_th, r_thresh_adjusted = r_thresh, at N_INP_adjusted = 0, r_thresh_adjusted = 1.2 * r_thresh [[ starts above r_thresh so we at least have some pitosn which can reduce qi and push us towards more pitosn as N shrinks (or if it's temp controlled at least we'll have some pitosn by 25% subsat -- it it's not strong enough we can still jump straight to 1.2 r_thresh but that might be too sharp given the big S_i gradients we have and the assumption that more native N_i really does slow pitosn...]]

                N_thresh_adjusted = N_from_qr(param_set, ice_type, q_i, r_thresh_adjusted; monodisperse = false, μ=μ, ρ=ρ)

                if N_INP_adjusted < N_thresh
                    if N_INP_adjusted < N_thresh_adjusted
                        # Bring N up towards N_th if it's below it, since we should err towards r_thresh, do not reduce N if it's above N_th
                        # N_INP_adjusted_here = N_INP_adjusted + (N_thresh_adjusted - N_INP_adjusted) * (-clamp(S_i, -0.05, 0) / 0.05) # by -5% subsat, hit r_thresh(adjusted). If we're already above r_thresh, 
                        # N_INP_adjusted_here = N_INP_adjusted + (N_thresh_adjusted - N_INP_adjusted) * (1 + clamp(S_i, -0.05, 0) / 0.05) # hit r_thresh (adjusted ) at 0% subsaturation, tapers to no change at -5% subsat

                        # We want to hit hit N_thresh at 0%, and N_thresh_adjusted at -5% subsat, with the understanding that PITOSN should take over after that.
                        N_INP_adjusted_here = N_thresh + (N_thresh_adjusted - N_thresh) * (-clamp(S_i, -0.05, 0) / 0.05)

                        N_INP_adjusted = min(N_INP_adjusted, 10*N_INP) # don't let it exceed 10*INP
                        if N_i_from_INP
                            N_i = N_INP_adjusted_here
                        else
                            N_i *= (N_INP_adjusted_here / N_INP_adjusted) # adjust N_i to match the new N_INP_adjusted
                        end
                        N_INP_adjusted = N_INP_adjusted_here
                    else
                        # r is between r_thresh and 1.1 r_thresh, so do nothing rn
                    end
                else
                    # do reduce but not as strongly as we boosted when N_INP_adjusted < N_thresh.
                    # We know we're in a region with a high capacity for INP generation, so we should be more cautious about reducing N too much.
                    # so instead of tapering off by 5% subsat, we'll extend it to 10% subsat
                    # we want to hit the original N_INP_adjusted at 0% subsat, N_thresh at -12.5% subsat, and N_thresh_adjusted at -25% subsat
                    if S_i > FT(-0.125)
                        N_INP_adjusted_here = N_INP_adjusted + (N_thresh - N_INP_adjusted) * (-clamp(S_i, -0.125, 0) / 0.125)
                    else
                        N_INP_adjusted_here = N_thresh + (N_thresh_adjusted - N_thresh) * ((-clamp(S_i, -0.25, -0.125) - 0.125) / 0.125)
                    end
                    #
                    
                    # Even when RH is super low, at cloud top this is still resulting in very high <r>. I think in scenarios when lambda is high and larger, falling INP arent advected up there, we're overly skewing to high <r>/...
                    # In narrow dry layers like RF10, we don't let <r> grow too much either since it takes time for size sorting to happen and high NINP means large amounts of ice probably came from saturation nearby.
                    # We can bound the reduction in N_INP_adjusted, but this might make it hard in calibration if N_INP is universally overpowering...
                    # N_INP_adjusted_here = max(N_INP_adjusted_here, N_INP_adjusted/2) # added 10/23 1028 am, /10 seems to be too strong, /2 is ok....could use raw N_INP but that feels too overpowering...


                    # Gradually reduce enforcement from full (N_INP_adjusted/2) at cloud top # [[ 10/31 to replace the line above that i think is too agressive low down and prevents autoconv in RF01 -- gradually go from full enforcement (N_INP_adjusted/2))]] to none at N_INP_adjusted = N_INP_top/10
                    if !isnan(N_INP_top)
                        w = clamp((N_INP_adjusted - N_INP_top/FT(10)) / (N_INP_top - N_INP_top/FT(10)), FT(0), one(FT))
                        N_limit = (one(FT) - w)*N_INP_adjusted_here + w*(N_INP_adjusted/2)
                        N_INP_adjusted_here = max(N_INP_adjusted_here, N_limit)
                    else
                        # Not sure
                    end

                    #
                    if N_i_from_INP
                        N_i = N_INP_adjusted_here
                    else
                        N_i *= (N_INP_adjusted_here / N_INP_adjusted) # adjust N_i to match the new N_INP_adjusted
                    end
                    N_INP_adjusted = N_INP_adjusted_here
                end
            else
                N_INP_adjusted = max(N_INP_adjusted, N_thresh/8) # added 12/01/2025 to prevent underflow probs
            end


            if apply_sedimentation_boost
                # sedimentation boost factor -- maybe 20 percent by .5 m/s? sedimentation only reaches as far as it can before acnv, so shorter at high RH, and at subsat it's longer except you could be losing particles.
                # maybe scale opposite of RH? idk.
                sedimentation_N_i_boost_factor::FT = FT(param_set.user_params.sedimentation_N_i_boost_factor)
                # sedimentation_N_i_boost_factor = FT(0)
                sedimentation_N_i_boost_factor *= clamp(w_i / FT(0.5), FT(0), FT(1)) # so at 0.5 m/s we get the full amount, at 0 m/s we get none
                # rh scaling [max at S_i = 0, by 8 percent supersat, reduce to 10 percent.]
                # sedimentation_N_i_boost_factor *= clamp(FT(1.0 - ((1.0 - 0.1) / 0.08) * (S_i - 1.0)), FT(0.1), FT(1.0)) % S_i already has the minus one...
                sedimentation_N_i_boost_factor *= clamp(FT(1 - (1.0 - 0.1)/0.08 * S_i), FT(0.1), FT(1.0))
                N_INP_adjusted *= (one(FT) + sedimentation_N_i_boost_factor)
                N_INP_adjusted = max(N_INP_adjusted, N_thresh/8) # if we have a negative massflux boost or massflux, don't allow that to deplete all N.... (allow r = 2 * r_thresh if needed)

                if N_i_from_INP
                    N_i *= (one(FT) + sedimentation_N_i_boost_factor)
                    N_i = max(N_i, N_thresh/8) # if we have a negative massflux boost or massflux, don't allow that to deplete all N.... (allow r = 2 * r_thresh if needed)
                else # still try it
                    N_i *= (one(FT) + sedimentation_N_i_boost_factor)
                    N_i = max(N_i, N_thresh/8) # if we have a negative massflux boost or massflux, don't allow that to deplete all N.... (allow r = 2 * r_thresh if needed)
                end

                if !isnan(N_INP_top)
                    N_INP_adjusted = min(N_INP_adjusted, N_INP_top) # don't let it go above N_INP_top w/o ice_mult...
                    N_i = min(N_i, N_INP_top) # don't let it go above N_INP_top
                end
            else
                # already done above
            end

            if !isfinite(N_INP_adjusted) || !isfinite(N_i)
                @warn "Got non-finite N_INP_adjusted = $N_INP_adjusted or N_i = $N_i from inputs N_INP = $N_INP; N_i_in = $N_i_in; q_i = $q_i; ρ = $ρ; r_is = $r_is; μ = $μ; massflux = $massflux; w_i = $w_i; S_i = $S_i; dNINP_dz = $dNINP_dz; N_INP_top = $N_INP_top; N_i_from_INP = $N_i_from_INP; q_l = $q_l; q_s = $q_s; apply_massflux_boost = $apply_massflux_boost; apply_sedimentation_boost = $apply_sedimentation_boost"
            end

        else
            N_INP_adjusted = max(N_INP, N_thresh/8) # changed to `N_thresh/8` 12/01/2025 to prevent underflow probs
        end
        # ----------------------------------------------------------------- #


        # N_bounds = (N_thresh, N_i_acnv) 
        N_bounds = (N_thresh/8, N_i_acnv) # allow r to grow up to 2x r_thresh. otherwise, just setting at r_thresh all the time doesn't lead to any PITOSN and hampers calibration... if N is inferred from q then it shouldn't be at r_thresh anyway.


        if !isnan(N_INP_top)
            # N_bounds = min.(N_bounds, N_INP_top)
            # we used to assume N_INP_top ? N_INP. This isn't strictly true, e.g. an inversion, but also because  f_mult_ice can increase INP upper bound
            N_bounds = min.(N_bounds, max(N_INP_top, N_INP_adjusted)) # This works but maybe it's giving too much leeway to increase N to keep <r> within bounds. We see a lot of N < N_INP_top but > N_INP where <r> = r_acnv... I don't understand why we get those values...
            # I think this is just always N_INP_top... N_INP_top > N_INP always presumably.

            # N_bounds = min.(N_bounds, N_INP) #  This works except we know N goes above N_INP a little bit... maybe something like  N_bounds = min.(N_bounds, 1.1 N_INP) is better or something idk...
            N_bounds = min.(N_bounds, min(1.2 * N_INP_adjusted, N_INP_top)) # Try this.. the higher the sedimentation, the lower the exported N_INP (larger particles), so we just will keep the same old N_INP

            # this line let's N_INP be very low, wich violates thee whole point of N_thresh... no? as in it breaks the N_thresh upper bound... at the same time I think we wanted to avoid just having some fixed value of <r> at the top... r_raw_max does put a bound on r_thresh that is 1.05 * r_thresh below though.

        
            # why do we bother with q_top at all?
            # N_INP sets an upper bound on N, which is a lower bound on <r>, but not an upper bound.
            # To upper bound <r>, we are using N relative to N_top
            # Note that the real starting <r> tends to be around r_is / 2 ≈ r_acnv, I'm not sure why. The max <q> seems to be around q_is * N_top though
            # For simplicity, we'll just assume we start at r = r_acnv and grow the upper bound from there..., it's a mixed bound, lower q_top -> smaller max <r> -> bigger min N, so looser in <r> but tighter in N... but i feel like it's better as it should enforce starting @ N_acnv.
            
            
            # I think this is bad bc any low down supersaturation could have its max r limited by some random initial remote q far above it even if tons of things happened in between... while any surrounding subsaturated  air would have no such limit... doesn't even make sense if SIP exists too.
            # q_top = q_from_rN(param_set, ice_type, r_i_acnv, N_INP_top; monodisperse = monodisperse, ρ=ρ, μ=μ) # this is the volume of a single ice crystal at the acnv radius, so N * q_top = q_threshold
            # if decrease_N_if_subsaturated && (S < FT(0))
            #     q_top = min(q_top, q_i) # allow q to turn back towards 0..., which constrain <r> to be smaller? (at q=0. it would be back to r_i_acnv)
            # end
            q_top = q_i


            
            # if N_i_from_INP # then we know the N we diagnosed was wrong and we should fully switch to this one...? unless it's the eTgi one... if it's somethign like geometric we dont know...
            # r_raw_max = r_from_qN(param_set, ice_type, q_top, N_INP_adjusted; monodisperse = monodisperse, ρ=ρ, μ=μ) # r is bounded by available INP from raw, if INP are available, <r> cannot grow.
            # r_raw_max = clamp(r_raw_max, r_i_acnv, r_thresh * FT(75/70)) # let the other 5 percent come from subsaturation if needed
            # N_raw_min = N_from_qr(param_set, ice_type, q_i, r_raw_max; monodisperse = monodisperse, ρ=ρ, μ=μ)
            # N_bounds = max.(N_bounds, N_raw_min)

            # N_thresh_adjusted = N_from_qr(param_set, ice_type, q_top, r_thresh * FT(75/70); monodisperse = monodisperse, ρ=ρ, μ=μ) # let the other 5 percent come from subsaturation if needed
            # N_bounds = max.(N_bounds, max(N_thresh_adjusted, min(N_i_acnv, N_INP_adjusted))) # I think we just need to prevent going to N_thresh if we have INP available..., this can allow unfettered growth of N_i once you hit r_thresh_adjusted...

            N_thresh_adjusted = N_from_qr(param_set, ice_type, q_top, r_thresh * FT(1.3); monodisperse = monodisperse, ρ=ρ, μ=μ) # let the other 5 percent come from subsaturation if needed
            N_bounds = max.(N_bounds, max(N_thresh_adjusted, min(N_i_acnv, N_INP_adjusted))) # bigger buffer so no unfettered growth...


            # don't let N_i grow w/o bound, cut it off at N_INP_adjusted, dont go to large N_INP if there's not q_i to support it.
            # we know N_i_acnv > N_thresh_adjusted, but if N_INP_adjusted was smaller, take it unless it's too far below N_thresh_adjusted.
            # with this system, we can't have runaway q_i,N_i becaue N_INP_adjusted will constrain us. if N_thresh_adjusted grows too large, there's nothing to stop it...

            # end

        else
            # I wish we could diagnose N_INP_top... but we can only guess it from q...
            #=
                We can surmise that if r < r_acnv if N_INP * r_acnv > q...

                Note though that we're tied to local q_i.
                This makes the bounding more local, we could try to estimate T_top based on assuming current q = q_top = N_INP(T_top) * q_is.
                To do that, we'd need to pass in relaxation_timescale etc and so forth in order to use `get_T_top_from_N_i()` etc.
                But that N_INP_top is always greater than the current N_i, so it is not a true upper bound.

                This is the same as if N_INP_top was provided if q_i = q_i_top... this is generally not a bad guess.

                They will diverge as q_i goes to 0, this will turn back towards r_max = 0 (r_acnv), even if we're otherwise fine.
                The advantage of knowing N_INP_top is you don't have to turn back per-se. In real clouds, <r> often doesn't turn back or barely does just at the very bottom in sublimation.
                I think you maybe in sublimation, even if you knew N_INP_top, you'd want to maybe allow a turnaround?
            =#

            # N_bounds = min.(N_bounds, N_INP) #  This works except we know N goes above N_INP a little bit... maybe something like  N_bounds = min.(N_bounds, 1.1 N_INP) is better or something idk...
            N_bounds = min.(N_bounds, 1.2 * N_INP_adjusted) # Try this.. the higher the sedimentation, the lower the exported N_INP (larger particles), so we just will keep the same old N_INP

            q_top = q_i # assumption

            # r_raw_max = r_from_qN(param_set, ice_type, q_top, N_INP_adjusted; monodisperse = monodisperse, ρ=ρ, μ=μ) # Here, we just check using the current INP. Do not let <r> grow beyond unless all INP are activated.
            # # because we are using q_i which is less than q_top in principle, r_raw_max will be smaller. as q goes to 0, r_raw_max will also go to 0. (r_i_acnv).
            # # So we have a battle between the local q and the N_INP. This is good, however it may push us back towards N_acnv when in reality that occurs very late in the process, and more through threshold acnv than evaporation sometimes...
            # r_raw_max = clamp(r_raw_max, r_i_acnv, r_thresh * FT(75/70)) # let the other 5 percent come from subsaturation if needed
            # N_raw_min = N_from_qr(param_set, ice_type, q_i, r_raw_max; monodisperse = monodisperse, ρ=ρ, μ=μ) # this would be a round trip except we clamped.

            N_thresh_adjusted = N_from_qr(param_set, ice_type, q_top, r_thresh * FT(2); monodisperse = monodisperse, ρ=ρ, μ=μ) # let the other 5 percent come from subsaturation if needed
            N_raw_min = max(N_thresh_adjusted, min(N_i_acnv, N_INP_adjusted)) # I think we just need to prevent going to N_thresh if we have IN

            N_bounds = max.(N_bounds, N_raw_min)

        end

        # Added 10/23/2025 -- I think we really do need some control on N outside of just the bounds from r_acnv and r_thresh, since those don't account for INP availability at all.
        # Since the incoming NINP accounts for ice mult we're ok there, and if q gets too large most new sub/dep goes straight to snow so there are limits to how big N can become...  I saw N_i = 7000+ in an updraft... ridiculous.
        # Let's bound ourselves to below N_INP_top for sure..., if we don't have it maybe 5xN_INP_adjusted? idk...
        N_below_r_is_top = (isnan(N_INP_top) || iszero(N_INP_top)) ? FT(Inf) : max(get_Ni_from_INP_qi_qs(param_set, N_i, N_INP_top, q_i, q_s; ρ=ρ, ice_type=ice_type, μ=μ, monodisperse=monodisperse, add_dry_aerosol_mass=true), N_thresh/8) # added min bound to prevent 0 division probs
        N_bounds = min.(N_bounds, N_below_r_is_top, (S_i < FT(0)) ? (8 * N_INP_adjusted) : (2 * N_INP_adjusted)) # not sure if should use N_INP_adjusted here or N_INP directly. I guess N_INP_adjusted is probably better since it includes MF boost, but it can also be too small with sedimentation... maybe only use it if S_i > 0? idk. it really is only a runaway growth stopper. we'll enforce a stronger limit in supersat when <r> isn't so big so terminal velocity is less important...? but also in growth regions with mf boost we should be limited by inp fairly strongly so maybe both can be more tightly bounded? idk.


        N_i = safe_clamp(N_i, N_bounds[1], N_bounds[2]) # clamp N_i to the bounds

        # technically we should go from N_INP_top * q_is and not allow q to go past that... but that should be emergent not enforced...

        # You should start @ N_raw and then move into desired range i suppose?
        # So, if N_raw is very high relative to N_top, move toward r_is, if N_raw is very low relative to N_top, move toward r_thresh.
    end

    if !isfinite(N_i) || (iszero(N_i) && (q_i > 0))
        @warn "Got non-finite or zero N_i = $N_i from inputs N_i_in = $N_i_in; N_INP = $N_INP; q_i = $q_i; ρ = $ρ; S_i = $S_i; q_l = $q_l; q_s = $q_s; monodisperse = $monodisperse; ice_type = $ice_type; decrease_N_if_subsaturated = $decrease_N_if_subsaturated; N_INP_top = $N_INP_top; massflux = $massflux; w_i = $w_i; N_i_from_INP = $N_i_from_INP"
    end



    return N_i
end
adjust_ice_N(param_set::APS, N_i::FT, N_INP::FT, q::TD.PhasePartition{FT}; ρ::FT = FT(1.0), S::FT = FT(1), monodisperse::Bool=true, ice_type::CMT.IceType = ice_type, decrease_N_if_subsaturated::Bool=true, N_INP_top::FT=FT(NaN), q_s::FT = FT(0), massflux::FT=FT(0), w_i::FT=FT(0), N_i_from_INP::Bool=false) where {FT} = adjust_ice_N(param_set, N_i, N_INP, q.ice; ρ=ρ, S=S, monodisperse=monodisperse, ice_type=ice_type, decrease_N_if_subsaturated=decrease_N_if_subsaturated, N_INP_top=N_INP_top, q_s=q_s, massflux=massflux, w_i=w_i, N_i_from_INP=N_i_from_INP)

    
    
function adjust_liq_N(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    N_l::FT, 
    q_l::FT; 
    ρ::FT = FT(1.0), # density of the air, if we have it, we can use it to adjust the source term
    S::FT = FT(1), # supersaturation
    monodisperse::Bool = true, # if true, we assume all liquid droplets are the same size, otherwise we assume a distribution of sizes
    liq_type::CMT.LiquidType = liq_type, # this is the type of liquid we're dealing with, e.g. rain or cloud liquid
    decrease_N_if_subsaturated::Bool = false
) where {FT}


    N_CCN = param_set.user_params.N_CCN
    N_l *= N_CCN/(100e6)
    
    # r_lr = param_set.user_params.r_liq_rain # There is no r_lr for liquid...instead we just set it in user_params.
    # ρ_l = CMP.ρ_cloud_liq(microphys_params)

    if monodisperse
        _q_lr_ = q_lr(param_set, liq_type) # this is the volume of a single liquid droplet, so N * q_lr = q_threshold
        # N_l = max(N_l, q_l / _q_lr_) # make sure existing liquid contributes to the size, particularly w/ sedimentation bringing in new liquid

        # really this is a highly peaked gamma distribution with μ of order 5-10... so <r> shouldn't be liq/rain... but maybe half that, which is 1/8 the volume. see https://ntrs.nasa.gov/api/citations/20220009419/downloads/SSantosJAMESLimitationsReprint.pdf
        N_l = max(N_l, q_l / (_q_lr_/8))


    else

        microphys_params::ACMP = TCP.microphysics_params(param_set)
        μ = μ_from_qN(param_set, liq_type, q_l, N_l; ρ=ρ) # this is the factor by which we scale the mean radius to get the mean radius for the ice crystals, so that we can use it in the N_i and N_l calculations

        # make sure <r> < (r_lr/2)
        r_lr::FT = param_set.user_params.r_liq_rain # There is no r_lr for liquid...instead we just set it in user_params.
        λ = (μ + FT(1)) / r_lr

        _χm::FT = get_χm(param_set, liq_type) 
        _m0::FT = m0(microphys_params, liq_type)
        _r0::FT = r0(microphys_params, liq_type) 
        _me::FT = me(microphys_params, liq_type) 
        _Δm::FT = Δm(microphys_params, liq_type) 
        N_lr = (q_l * ρ) / (_χm * _m0 * _r0^(-(_me + _Δm)) * (CM1.SF.gamma(μ + _me + _Δm + FT(1)) / CM1.SF.gamma(μ + FT(1))) * λ^(-(_me + _Δm))) # this is the number concentration of ice crystals at the ice-snow radius
        N_l = max(N_l, N_lr) # make sure we don't get too many liquid droplets
    end

    return N_l
end
adjust_liq_N(param_set::APS, N_l::FT, q::TD.PhasePartition{FT}; ρ::FT = FT(1.0), S::FT = FT(1), monodisperse::Bool=true, liq_type::CMT.LiquidType = liq_type, decrease_N_if_subsaturated::Bool = false) where {FT} = adjust_liq_N(param_set, N_l, q.liq; ρ=ρ, S=S, monodisperse=monodisperse, liq_type=liq_type, decrease_N_if_subsaturated=decrease_N_if_subsaturated)


"""
Calculate r from q, N, and ρ.
mass_scaling_factor exists so you can choose dfiferent estimates of <r> based on the q and N you're using.
For example in gamma distribution, <r> = 1/2λ, q = π ρ N_0 / λ^4 and N = N_0 / λ so r = ( q / ( 8 π N ρ))^(1/3)

We are currently using q = 4/3 π r^3 N ρ, so r = (q / (4/3 π N ρ))^(1/3)

The difference in factor is then  FT(((4/3) / 8)^(1/3)) ≈ 0.55, that is, our r estimates will be double the size they should be, and our timescales will be about half of what they should be.
This sounds small but it doubles any sources... if you rely on having a good prediction of N... this could result in say deposition sources being 2x too large, which isn't insignificant.

Could try `TCP.TurbulenceConvectionParameters{FT}` over APS

"""

function r_from_qN(param_set::APS, q_type::CMTWaterTypes, q::FT, N::FT; monodisperse::Bool=true, μ::FT=FT(NaN), ρ::FT=FT(1), Dmin::FT=FT(0), Dmax::FT=FT(Inf)) where {FT}
    """
    Calculate r from q, N
    """
    r_min::FT = FT(param_set.user_params.particle_min_radius) # this is the minimum radius we want to consider, if r < r_min, we return 0
    if iszero(N) || isinf(N) || iszero(q)
        return FT(r_min) # N/N will give NaN...
        # even if N = 0 and q is not 0 (in which case N shouldn't be 0), we need some fix...
        # if N = Inf, also bad but don't let r go below r_min, 1/Nr = 1/Inf will yield timescale = 0
        # If N = 0, we'll just return r_min, then 1/Nr = 1/0 will yield timescale = Inf
    end



    microphys_params::ACMP = TCP.microphysics_params(param_set)
    _χm::FT = get_χm(param_set, q_type) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
    q_r_min = particle_mass(microphys_params, q_type, r_min, _χm)

    if (q < 0) || (N < 0)
        # return FT(r_min)
        error("Got negative q = $q or q_r_min = $q_r_min; inputs were q_type = $q_type; N = $N; r_min = $r_min; ρ = $ρ; μ = $μ")
    end

    if isnan(N)
        _n0 = n0(microphys_params, q, ρ, q_type)
        # λ = CM1.lambda(microphys_params, q_type, q, ρ)
        λ = lambda(microphys_params, q_type, q, ρ,  N, μ, Dmin, Dmax; _n0 = _n0, _χm=_χm) # we need λ to get N from n0
        N = _n0 / λ # needed for q_r_min
        if isnan(μ) # if μ is NaN, we need to calculate it
            μ = μ_from_qN(param_set, q_type, q + N*q_r_min, N; ρ=ρ) # this is a test to get the μ value for the given N and q, so we can use it to scale the mean radius
        end
        if (_n0 * CM1.SF.gamma(μ + 1) / λ^(μ + 1)) < zero(FT)
            @error "Got negative N = $(_n0 * CM1.SF.gamma(μ + 1) / λ^(μ + 1)) from inputs q = $q; N = $N; Dmin = $Dmin; Dmax = $Dmax; μ = $μ; ρ = $ρ; _n0 = $_n0; λ = $λ"
        elseif !isfinite(_n0 * CM1.SF.gamma(μ + 1) / λ^(μ + 1))
            @error "Got non-finite N = $(_n0 * CM1.SF.gamma(μ + 1) / λ^(μ + 1)) from inputs q = $q; N = $N; Dmin = $Dmin; Dmax = $Dmax; μ = $μ; ρ = $ρ; _n0 = $_n0; λ = $λ"
        end
        N = _n0 * CM1.SF.gamma(μ + 1) / λ^(μ + 1) # this is the number concentration of particles at the ice-snow radius
    end


    if monodisperse
        return particle_radius_from_mass(microphys_params, q_type, q/N + q_r_min, _χm) # we add the q_r_min to get the effective q, so we can use the same r_min for all types
    else
        # error("r_from_qN for non-monodisperse not implemented yet. but with q and N it can be.")
        # we want ∫ r n(r) dr / ∫ n(r) dr = ∫ r n(r) dr / N = <r>... We here ignore Dmin, Dmax though...
        if isnan(μ) # if μ is NaN, we need to calculate it
            # μ = FT(0) # we don't really have infrastructure to calculate μ in lambda and other functions right now bc it relies on user params, maybe i should look into fixing that...
            μ = μ_from_qN(param_set, q_type, q + N*q_r_min, N; ρ=ρ) # this is a test to get the μ value for the given N and q, so we can use it to scale the mean radius
        end
        λ = lambda(microphys_params, q_type, q + N*q_r_min, ρ, N, μ, Dmin, Dmax)
        if iszero(Dmin) && isinf(Dmax)
            r = (μ + 1) / (λ) # see paper derivation ( this is true from 0 to inf, idk about otherwise...)
        else
            @error "how did we end up here? inputs were q = $q; N = $N; Dmin = $Dmin; Dmax = $Dmax; μ = $μ; ρ = $ρ"
            r = int_nr_dr(microphys_params, q_type, q + N*q_r_min, ρ, μ; Nt=N, Dmin=Dmin, Dmax=Dmax) / int_n_dr(microphys_params, q_type, q + N*q_r_min, ρ, μ; Nt=N, Dmin=Dmin, Dmax=Dmax) # this is the integral of r n(r) dr / ∫ n(r) dr = <r>... we ignore Dmin, Dmax though...
        end

        if !isfinite(r)
            if N < zero(FT)
                error("Got non-finite r = $r from inputs q = $q; q_type = $q_type; N = $N; Dmin = $Dmin; Dmax = $Dmax; μ = $μ; ρ = $ρ")
            else
                error("Got non-finite r = $r from inputs q = $q; q_type = $q_type; N = $N; Dmin = $Dmin; Dmax = $Dmax; μ = $μ; ρ = $ρ")
            end
            r = FT(r_min)
        end

        return r
    end
end

ri_from_qN(param_set::APS, q::FT, N::FT; monodisperse::Bool=true, ice_type::CMT.IceType = ice_type, μ::FT=FT(NaN), ρ::FT=FT(1), Dmin::FT=FT(0), Dmax::FT=FT(Inf)) where {FT} = r_from_qN(param_set, ice_type, q, N; monodisperse = monodisperse, μ=μ, ρ=ρ, Dmin=Dmin, Dmax=Dmax) # shorthand for ice type
rl_from_qN(param_set::APS, q::FT, N::FT; monodisperse::Bool=true, liq_type::CMT.LiquidType = liq_type, μ::FT=FT(NaN), ρ::FT=FT(1), Dmin::FT=FT(0), Dmax::FT=FT(Inf)) where {FT} = r_from_qN(param_set, liq_type, q, N; monodisperse = monodisperse, μ=μ, ρ=ρ, Dmin=Dmin, Dmax=Dmax) # shorthand for liquid type

function r_from_q(param_set::APS, q_type::CMTWaterTypes, q::FT) where {FT}
    """
    Calculate r from q
    """
    r_min::FT = param_set.user_params.particle_min_radius # this is the minimum radius we want to consider, if r < r_min, we return 0
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    χm::FT = get_χm(param_set, q_type) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
    q_r_min = particle_mass(microphys_params, q_type, r_min, χm) # this is the volume of a single droplet with radius r_min, so we add it to q to get the effective q
    return particle_radius_from_mass(microphys_params, q_type, q+q_r_min, χm) # we use q=0 to get the mean radius for the given N
end
ri_from_q(param_set::APS, q::FT; ice_type::CMT.IceType = ice_type) where {FT} = r_from_q(param_set, ice_type, q) # shorthand for ice type
rl_from_q(param_set::APS, q::FT; liq_type::CMT.LiquidType = liq_type) where {FT} = r_from_q(param_set, liq_type, q) # shorthand for liquid type

# -------------------------- #

function q_from_rN(param_set::APS, q_type::CMTWaterTypes, r::FT, N::FT; monodisperse::Bool=true, μ::FT=FT(0), ρ::FT=FT(1), _χm::FT=FT(NaN)) where {FT}
    """
    Calculate q from r, N
    """
    if iszero(N)
        return FT(0) # N/N will give NaN...
    elseif isinf(N)
        error("N is infinite, cannot calculate q from r and N")
    else
        r_min::FT = param_set.user_params.particle_min_radius # this is the minimum radius we want to consider, if r < r_min, we return 0
        microphys_params::ACMP = TCP.microphysics_params(param_set)
        # _χm::FT = get_χm(param_set, q_type) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
        _χm::FT = isnan(_χm) ? get_χm(param_set, q_type) : _χm
        if r > r_min
            if monodisperse
                return mass(microphys_params, q_type, r, N, _χm; monodisperse=true) - mass(microphys_params, q_type, r_min, N, _χm; monodisperse=true) # we subtract the r_min radius sphere volume to get the effective q
            else
                # _χm = get_χm(param_set, q_type) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
                # _χm::FT = isnan(_χm) ? get_χm(param_set, q_type) : _χm
                _m0::FT = m0(microphys_params, q_type) # this is the mass of the droplet at the acnv radius
                _r0::FT = r0(microphys_params, q_type) # this is the radius of the droplet at the acnv radius
                _me::FT = me(microphys_params, q_type)
                _Δm::FT = Δm(microphys_params, q_type)
                _λ = (μ+1) / r
                _n0 = n0_from_Nλ(N, _λ; μ=μ)

                q = _χm * _m0 * _n0 * CM1.SF.gamma(μ + _me + _Δm + FT(1)) / ( _r0^(_me + _Δm) * _λ^(μ + _me + _Δm + FT(1)) )
                q = max(FT(0), q - mass(microphys_params, q_type, r_min, N, _χm; monodisperse=true))
                q /= ρ  # 1/ρ to go from /m^3 to kg/kg
                return q

                # error("q_from_rN for non-monodisperse not implemented yet. but with r and N and μ perhaps it can be.") # i'm not sure how easy it is to go from r, N to q... is that enough to define λ?
            end
        else
            # consider erroring for r < r_min
            return FT(0)
        end
    end
end
qi_from_rN(param_set::APS, r::FT, N::FT; monodisperse::Bool=true, ice_type::CMT.IceType = ice_type, μ::FT=FT(0), ρ::FT=FT(1), _χm::FT = FT(NaN)) where {FT} = q_from_rN(param_set, ice_type, r, N; monodisperse = monodisperse, μ=μ, ρ=ρ, _χm=_χm) # shorthand for ice type
ql_from_rN(param_set::APS, r::FT, N::FT; monodisperse::Bool=true, liq_type::CMT.LiquidType = liq_type, μ::FT=FT(0), ρ::FT=FT(1), _χm::FT = FT(NaN)) where {FT} = q_from_rN(param_set, liq_type, r, N; monodisperse = monodisperse, μ=μ, ρ=ρ, _χm=_χm) # shorthand for liquid type

function q_from_r(param_set::APS, q_type::CMTWaterTypes, r::FT) where {FT}
    """
    Calculate q from r
    """
    r_min::FT = param_set.user_params.particle_min_radius # this is the minimum radius we want to consider, if r < r_min, we return 0
    if r > r_min
        microphys_params::ACMP = TCP.microphysics_params(param_set)
        χm::FT = get_χm(param_set, q_type) # this is the mass scaling factor for the mass diameter relationship, so we can use it to scale the mean radius
        return particle_mass(microphys_params, q_type, r, χm) - particle_mass(microphys_params, q_type, r_min, χm) # we subtract the r_min radius sphere volume to get the effective q
    else
        return FT(0) # consider erroring for r < r_min
    end
end
qi_from_r(param_set::APS, r::FT; ice_type::CMT.IceType = ice_type) where {FT} = q_from_r(param_set, ice_type, r) # shorthand for ice type
ql_from_r(param_set::APS, r::FT; liq_type::CMT.LiquidType = liq_type) where {FT} = q_from_r(param_set, liq_type, r) # shorthand for liquid type


    

"""
We have no way of passing in a minimum aerosol size into λ, N etc in terminal_velocity...
Instead, we calculate what q would be if there was a r_min radius sphere in the center of every drop.
q0 = liquid in one droplet

but then q and N won't match no? you'll have one giant drop... which will have a very high terminal velocity... so adusting N really is the right call...

You would then pass this q into sedimentation to get a more realistic speed for extreme N
"""
function q_effective(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q_type::CMTWaterTypes,
    q::FT,
    N::FT,
    r_min::FT;
    monodisperse::Bool = true, # if true, we assume the q is for a monodisperse distribution, otherwise we assume it's for a polydisperse distribution
) where {FT}
    """
    effective q for a given N, r, ρ that basically shoves the r_min sphere into the center of every drop
    """
    # ρ_eff = ρ * mass_scaling_factor # we scale the density by the mass_scaling_factor to get the effective density for the mass diameter relationship (we assume same power, but different scaling)
    # return q + FT(4/3) * FT(π) * ρ_eff * N * (r_min^3) 
    return q + mass(param_set, q_type, r_min, N; monodisperse = monodisperse) # this is the volume of a single droplet with radius r_min, so we add it to q to get the effective q
end


#= 
q if there were no r_min sized aerosols in the middle of every drop, but we rinstead replaced that volume with liquid...
Useful for things that only care about droplet size but estimate it from a q and N only (like my sedimentation terminal velocity implementation)
I suppose you could try to figure how how they're doing that q <--> N mapping, but when we supply our own N, things can get wonky for very large/small N so hopefully this helps.
=#
# q_effective_nan_N_safe(q::FT, N::FT, ρ::FT, r_min::FT) where{FT} = isnan(N) ? q : q_effective(q, N, ρ, r_min)
# q_effective_nan_N_safe(param_set::APS, q::FT, N::FT, ρ::FT) where{FT} = isnan(N) ? q :  q_effective(q, N, ρ, FT(param_set.user_params.particle_min_radius))

q_effective_nan_N_safe(param_set::APS, q_type::CMTWaterTypes, q::FT, N::FT; monodisperse::Bool=true) where{FT} = isnan(N) ? q : q_effective(param_set, q_type, q, N, FT(param_set.user_params.particle_min_radius); monodisperse = monodisperse)


