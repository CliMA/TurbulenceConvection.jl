"""
Get timescales for relaxation of the supersaturation


Throughout, mean_r_factor exists so you can choose dfiferent estimates of <r> based on the q and N you're using.
For example in gamma distribution, <r> = 1/2λ, q = N_0 Γ(5) / λ^4 , and N = N_0/λ so <r> = (32q/N)^1/3. Then to go from r = (q/(4/3πρN))^(1/3) to  (32q/N)^1/3 we need to multiply by (24/(πρ))^(1/3)
"""

# N can be arbitrarily small so D * N_l could yield 0, times r = Inf could yield NaN... always group (N * r) to avoid this, can't have Inf and 0 here bc of condition in r_from_qN()

# ======================================================================================================================================== #
# D_func(T::FT, p::FT) where {FT} = FT((2.11 * 1e-5) * (T / FT(273.15))^1.94 * (p / 101325))::FT  # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
D_func(T::FT, p::FT) where {FT} = FT((2.11 * 1e-5) * (T / FT(273.15))^1.94 * (FT(101325)/p))::FT  # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226) [e.g. 12.1 in https://apps.dtic.mil/sti/tr/pdf/ADA440352.pdf]
# ======================================================================================================================================== #

"""
Calculate the relaxation timescale τ. 
Assumes that N is adjusted for q.  
If there is a raw N_raw value, we assume that that anything in N_raw - N is at particle_min_radius
"""
function calculate_τ(D::FT, N::FT, r::FT, ρ::FT; 
    N_raw::FT = FT(N),
    particle_min_radius::FT = FT(0.2e-6),
    scaling_factor::FT = FT(1)
) where {FT}
    # return inv(FT(4) * FT(π) * D * N * r * ρ) * scaling_factor
    return inv(FT(4) * FT(π) * D * (N * r + max(FT(0), N_raw - N) * particle_min_radius) * ρ) * scaling_factor
end
# ======================================================================================================================================== #

function get_τs(param_set::APS, microphys_params::ACMP, relaxation_timescale::AbstractRelaxationTimescaleType, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT = FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false, apply_supersaturation_scaling=true) where {FT}
    τ_liq, τ_ice = get_τ_helper(param_set, microphys_params, relaxation_timescale, q, T, p, ρ, w; N_l = N_l, N_i_raw = N_i_raw, N_i_adjusted = N_i_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)

    # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster... change to be dt related) [used to have FT(1)]
    τ_liq = clamp(τ_liq, relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    τ_ice = clamp(τ_ice, relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)

    # if we're below critical subsat, tau_ice can go to in regardless of max_τ_ice
    # r0_activation = FT(10e-6) # 10 μm particle post activation, matches https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L1089
    # m0_activation = particle_mass(param_set, ice_type, r0_activation)

    if apply_supersaturation_scaling # need to do to Ns too
        thermo_params::TDPS = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        if S_i > zero(FT)
            # we use FT(Inf) instead of q_ice_threshold so that we are never forced to 1 :: This is useful because ??
            # This is bad because even when we have plenty of q, the timescale can get decimated....
            # Also isn't N already adjusted? Why do we need to readjust based on supersaturation again? At least for INP_Aware_Timescales?
            r0_activation = FT(10e-6) # 10 microns
            m0_activation = particle_mass(param_set, ice_type, r0_activation)
            # scaling_factor = INP_supersaturation_scaling_factor(param_set, S_i, q.ice, FT(Inf)) # never saturation q_ice_threshold so that we are never forced to 1 [[ is this a mistake? ]]
            scaling_factor = INP_supersaturation_scaling_factor(param_set, S_i, q.ice, m0_activation) # never saturation q_ice_threshold so that we are never forced to 1 [[ is this a mistake? ]]
            scaling_factor = iszero(scaling_factor) ? FT(Inf) : 1/scaling_factor
        else
            scaling_factor = one(FT)
        end
        τ_ice = clamp(τ_ice, relaxation_timescale.args.min_τ_ice * scaling_factor, relaxation_timescale.args.max_τ_ice * scaling_factor) # this can go below the lower limit as written, meaning with S_i > 0 can be slower than with < 0...
    end


    if isnan(τ_liq) || isnan(τ_ice)  # chasing down NaNs
        error("Timescale calculation failed. Got τ_liq: $τ_liq, τ_ice: $τ_ice for inputs q = $q; T = $T; p = $p; ρ = $ρ; w = $w; relaxation_timescale = $relaxation_timescale; N_l = $N_l; N_i_raw = $N_i_raw; N_i_adjusted = $N_i_adjusted; N_INP_top = $N_INP_top; f_ice_mult = $f_ice_mult; q_sno = $q_sno; massflux = $massflux; dTdz = $dTdz; w_i = $w_i; apply_massflux_boost = $apply_massflux_boost; apply_sedimentation_boost = $apply_sedimentation_boost")
    end

    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end
get_τs(param_set::APS, relaxation_timescale::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT = FT(0), w_i::FT = FT(0), apply_massflux_boost::Bool = false, apply_sedimentation_boost::Bool = false) where {FT} = get_τs(param_set, TCP.microphysics_params(param_set), relaxation_timescale, q, T, p, ρ, w; N_l = N_l, N_i_raw = N_i_raw, N_i_adjusted = N_i_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)


# the symbol methods now use namelist...
# function get_τs(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::Symbol, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT) where {FT}
#     return get_τs(param_set, microphys_params, get_relaxation_timescale_type(relaxation_timescale_type, param_set), q, T, p, ρ, w)
# end
# function get_τs(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT} # should not use this if possible bc it constructs another relaxation tiemscale object...
#     return get_τs(param_set, microphys_params, get_relaxation_timescale_type(relaxation_timescale_type, param_set), ts, w)
# end
# get_τs(param_set::APS, relaxation_timescale_type::Symbol, ts::TD.ThermodynamicState, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT) where {FT} = get_τs(param_set, TCP.microphysics_params(param_set), relaxation_timescale_type, q, T, p, ρ, w)
# get_τs(param_set::APS, relaxation_timescale_type::Symbol, ts::TD.ThermodynamicState, w::FT) where{FT} = get_τs(param_set, TCP.microphysics_params(param_set), relaxation_timescale_type, ts, w)


function get_τs(param_set::APS, microphys_params::ACMP, relaxation_timescale::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT = FT(0), w_i::FT = FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_τs(param_set, microphys_params, relaxation_timescale, q, T, p, ρ, w; N_l = N_l, N_i_raw = N_i_raw, N_i_adjusted = N_i_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end
get_τs(param_set::APS, relaxation_timescale::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT = FT(0), w_i::FT = FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT} = get_τs(param_set, TCP.microphysics_params(param_set), relaxation_timescale, ts, w; N_l = N_l, N_i_raw = N_i_raw, N_i_adjusted = N_i_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)

# ---------------------------------------------------------------------------------------------------------------------------------------- #

# :Base
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT=FT(NaN), N_i_raw::FT=FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT = FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    return relaxation_timescale.τ_liq, relaxation_timescale.τ_ice
end

# :exponential_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::ExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT = FT(0), w_i::FT = FT(0), apply_massflux_boost::Bool = false, apply_sedimentation_boost::Bool = false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    # N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w)

    if isnan(N_i_raw) || isnan(N_i_adjusted)
        (; N_i_raw, N_i_adjusted) = get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    end

    r_i = ri_from_qN(param_set, q.ice, N_i_adjusted; monodisperse = false, μ=FT(NaN), ρ=ρ)
    τ_ice = calculate_τ(D, N_i_adjusted, r_i, ρ; N_raw = N_i_raw, particle_min_radius = param_set.user_params.particle_min_radius, scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :exponential_T_scaling_ice_raw
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::ExponentialTScalingIceRawRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    τ_ice = FT(1 / (c_1 * exp(c_2 * (T - T_fr))))
    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :powerlaw_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::PowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)
    if T >= T_fr
        τ_ice = FT(Inf)
    else
        # N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w)
        if isnan(N_i_raw) || isnan(N_i_adjusted)
            (N_i_raw, N_i_adjusted) = get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        end
        r_i = ri_from_qN(param_set, q.ice, N_i_adjusted; monodisperse = false, μ=FT(NaN), ρ=ρ)
        # τ_ice = FT(1 / (4 * π * D * (N_i * r_i))) * relaxation_timescale.τ_sub_dep_scaling_factor
        τ_ice = calculate_τ(D, N_i_adjusted, r_i, ρ; N_raw = N_i_raw, particle_min_radius = param_set.user_params.particle_min_radius, scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    end
    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :exponential_times_powerlaw_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::ExponentialTimesPowerlawScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    error("NotImplmentedError: This relaxation_timescale functionality has not been implemented yet")
end

# :geometric_liq__geometric_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::GeometricLiqGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)

    if isnan(N_l)
        N_l = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    end
    r_l = rl_from_qN(param_set, q.liq, N_l; monodisperse = true)
    # τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) * relaxation_timescale.τ_cond_evap_scaling_factor
    τ_liq = calculate_τ(D, N_l, r_l, ρ; scaling_factor = relaxation_timescale.τ_cond_evap_scaling_factor)

    if isnan(N_i_raw) || isnan(N_i_adjusted)
        N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    else
        N_i = N_i_adjusted
    end
    r_i = ri_from_qN(param_set, q.ice, N_i; monodisperse = false, μ=FT(NaN), ρ=ρ)
    # τ_ice = FT(1 / (4 * π * D * (N_i * r_i))) * relaxation_timescale.τ_sub_dep_scaling_factor
    τ_ice = calculate_τ(D, N_i, r_i, ρ; scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :geometric_liq__exponential_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::GeometricLiqExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)

    if isnan(N_l)
        N_l = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    end
    r_l = rl_from_qN(param_set, q.liq, N_l; monodisperse = true)
    # τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) * relaxation_timescale.τ_cond_evap_scaling_factor
    τ_liq = calculate_τ(D, N_l, r_l, ρ; scaling_factor = relaxation_timescale.τ_cond_evap_scaling_factor)

    if isnan(N_i_raw) || isnan(N_i_adjusted)
        # N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w)
        (N_i_raw, N_i_adjusted) = get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    end
    r_i = ri_from_qN(param_set, q.ice, N_i_adjusted; monodisperse = false, μ=FT(NaN), ρ=ρ)
    # τ_ice = FT(1 / (4 * π * D * (N_i * r_i))) * relaxation_timescale.τ_sub_dep_scaling_factor
    τ_ice = calculate_τ(D, N_i_adjusted, r_i, ρ; N_raw = N_i_raw, particle_min_radius = param_set.user_params.particle_min_radius, scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :geometric_liq__powerlaw_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::GeometricLiqPowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)

    if isnan(N_l)
        N_l = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    end
    r_l = rl_from_qN(param_set, q.liq, N_l; monodisperse = true)
    # τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) * relaxation_timescale.τ_cond_evap_scaling_factor 
    τ_liq = calculate_τ(D, N_l, r_l, ρ; scaling_factor = relaxation_timescale.τ_cond_evap_scaling_factor)

    if T >= T_fr
        τ_ice = FT(Inf)
    else
        if isnan(N_i_raw) || isnan(N_i_adjusted)
            # N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w)
            (N_i_raw, N_i_adjusted) = get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        end
        r_i = ri_from_qN(param_set, q.ice, N_i_adjusted; monodisperse = false, μ=FT(NaN), ρ=ρ)
        # τ_ice = FT(1 / (4 * π * D * (N_i * r_i))) * relaxation_timescale.τ_sub_dep_scaling_factor
        τ_ice = calculate_τ(D, N_i_adjusted, r_i, ρ; N_raw = N_i_raw, particle_min_radius = param_set.user_params.particle_min_radius, scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    end
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

"""
    [4πDNr]⁻¹ was what I'd been using but I believe we do need to use ρ?
    N is already in m⁻³, so we need to work in those units.

    (from the amount of vapor perspective including ρ doesn't make a difference, but that vapor is all going to N, and how much air for each droplet is dependent on ρ_a
    
    The amount of water vapor available then is ρ

    (q_v - q_v_sat) * ρ_a =   (ρ_v - ρ_v_sat)/ρ_a * ρ_a = (ρ_v - ρ_v_sat) which is physical units...
"""


# :geometric_liq__exponential_T_scaling_and_geometric_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)

    if isnan(N_l)
        N_l = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    end
    r_l = rl_from_qN(param_set, q.liq, N_l; monodisperse = true)
    # τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) * relaxation_timescale.τ_cond_evap_scaling_factor
    τ_liq = calculate_τ(D, N_l, r_l, ρ; scaling_factor = relaxation_timescale.τ_cond_evap_scaling_factor)

    # N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, q_sno = q_sno, massflux = massflux, w_i = w_i)
    # r_i = ri_from_qN(param_set, q.ice, N_i; monodisperse = false, μ=FT(NaN), ρ=ρ)
    # τ_ice = 1 / (4 * π * D * (N_i * r_i)) * relaxation_timescale.τ_sub_dep_scaling_factor
    # τ_ice = calculate_τ(D, N_i, r_i, ρ; scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)


    #=  Testing having just a fallback growth rate that just uses c_4 directly rather than adding something to the timescale regardless..
        This can help decouple the N(T)-τ(T) relationship a little bit during activation.... and remote the artificial temperature contribution later as we move to lower N

        I think ordinarily c_3 would get suppressed low so that it didn't contribute at low N, and then growth at cloud top would be forced to be sluggish with no N to start and growing essentially only by q which we've found doesn't work in time for sedimentation. We couldn't raise c_3 because then the timescale would be too fast down low.
        Also c_3 and c_4 were redundant, but we needed c_3 so c_4 could match the other exponential T scaling models. This way, we hopefully don't need _c_3 at all.
    =#
    
    # c_4 = relaxation_timescale.c_4i
    # c_5 = relaxation_timescale.c_5i
    # N_i_T = c_4 * exp(c_5 * (T - T_fr))
    # N_i_T = max(N_i_T - N_i, FT(0))  # ensure N_i_T is non-negative

    # N_i_T = clamp(N_i_T, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    # r_i_T = param_set.user_params.particle_min_radius
    # τ_ice = 1 / (4 * π * D * ((N_i * r_i) + (N_i_T * r_i_T))) * relaxation_timescale.τ_sub_dep_scaling_factor
    if isnan(N_i_raw) || isnan(N_i_adjusted)
        # N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w)
        (N_i_raw, N_i_adjusted) = get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    end
    r_i = ri_from_qN(param_set, q.ice, N_i_adjusted; monodisperse = false, μ=FT(NaN), ρ=ρ)

    # τ_ice = calculate_τ(D, N_i + N_i_T, ((N_i * r_i) + (N_i_T * r_i_T))/(N_i + N_i_T), ρ; scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    # Here we have to be more careful because the value el,gTi returns is teh minimum of N(q) and N(T) so N(T) is what we want to pass as raw,
    τ_ice = calculate_τ(D, N_i_adjusted, r_i, ρ; N_raw = N_i_raw, particle_min_radius = param_set.user_params.particle_min_radius, scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)

    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :linear_combination
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::LinearCombinationRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)


    if isnan(N_l)
        N_l = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    end
    r_l = rl_from_qN(param_set, q.liq, N_l; monodisperse = true)
    # τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) * relaxation_timescale.τ_cond_evap_scaling_factor 
    τ_liq = calculate_τ(D, N_l, r_l, ρ; scaling_factor = relaxation_timescale.τ_cond_evap_scaling_factor)

    if isnan(N_i_raw) || isnan(N_i_adjusted)
        N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    else
        N_i = N_i_adjusted
    end
    r_i = ri_from_qN(param_set, q.ice, N_i; monodisperse = false, μ=FT(NaN), ρ=ρ)
    # τ_ice = FT(1 / (4 * π * D * (N_i * r_i))) * relaxation_timescale.τ_sub_dep_scaling_factor
    τ_ice = calculate_τ(D, N_i, r_i, ρ; scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)

    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :linear_combination_with_w
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::LinearCombinationWithWRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)
    #
    if isnan(N_l)
        N_l = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    end
    r_l = rl_from_qN(param_set, q.liq, N_l; monodisperse = true)
    # τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) * relaxation_timescale.τ_cond_evap_scaling_factor 
    τ_liq = calculate_τ(D, N_l, r_l, ρ; scaling_factor = relaxation_timescale.τ_cond_evap_scaling_factor)

    #
    if isnan(N_i_raw) || isnan(N_i_adjusted)
        N_i = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    else
        N_i = N_i_adjusted
    end
    r_i = ri_from_qN(param_set, q.ice, N_i; monodisperse = false, μ=FT(NaN), ρ=ρ)
    # τ_ice = FT(1 / (4 * π * D * (N_i * r_i))) * relaxation_timescale.τ_sub_dep_scaling_factor
    τ_ice = calculate_τ(D, N_i, r_i, ρ; scaling_factor = relaxation_timescale.τ_sub_dep_scaling_factor)
    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end

# :neural_network & :neural_network_no_weights & :neural_network_random_init & neural_network_pca_noise
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_l::FT = FT(NaN), N_i_raw::FT = FT(NaN), N_i_adjusted::FT = FT(NaN), N_INP_top::FT = FT(NaN), f_ice_mult::FT = FT(1), q_sno::FT = FT(0), massflux::FT = FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    # convert params to static_strided_array as simplechains doesn't accept svectors directly [if it comes to it we can set the static_strided_array as a global, it more than doubles the call time of the NN)]
    τ_liq, τ_ice, _, _ = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?

    return (; τ_liq = τ_liq, τ_ice = τ_ice)
end



# ======================================================================================================================================== #




# ======================================================================================================================================== #


function get_N_i(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_i = get_N_i_helper(param_set, relaxation_timescale, q, T, ρ, w)
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = (relaxation_timescale isa INP_Aware_Timescale)
        N_i = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    end
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return N_i
end
# ======================================================================================================================================== #

    
function get_N_i_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_i = get_N_i_helper(param_set, relaxation_timescale, q, T, ρ, w)
    N_i_no_boost = N_i
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = (relaxation_timescale isa INP_Aware_Timescale)
        N_i = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        if apply_massflux_boost
            N_i_no_boost = adjust_ice_N(param_set, N_i_no_boost, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = false, apply_sedimentation_boost = false)
        else
            N_i_no_boost = N_i
        end
    end
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_no_boost = clamp(N_i_no_boost, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return (; N_i, N_i_no_boost)
end
# ======================================================================================================================================== #

function get_N_i_raw_and_adjusted(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_i = get_N_i_helper(param_set, relaxation_timescale, q, T, ρ, w)
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = (relaxation_timescale isa INP_Aware_Timescale)
        N_i_adjusted = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    else
        N_i_adjusted = N_i
    end
    N_i_raw = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_adjusted = clamp(N_i_adjusted, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return (; N_i_raw, N_i_adjusted)
end
# ======================================================================================================================================== #

function get_N_i_raw_and_adjusted_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_i = get_N_i_helper(param_set, relaxation_timescale, q, T, ρ, w)
    N_i_no_boost = N_i
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = (relaxation_timescale isa INP_Aware_Timescale)
        N_i_adjusted = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        if apply_massflux_boost
            N_i_no_boost = adjust_ice_N(param_set, N_i_no_boost, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = false, apply_sedimentation_boost = false)
        else
            N_i_no_boost = N_i
        end
    else
        N_i_adjusted = N_i
    end
    N_i_raw = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_adjusted = clamp(N_i_adjusted, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_no_boost = clamp(N_i_no_boost, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return (; N_i_raw, N_i_adjusted, N_i_no_boost)
end
# ======================================================================================================================================== #

function get_N_i_raw(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_i = get_N_i_helper(param_set, relaxation_timescale, q, T, ρ, w)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return N_i
end
# ======================================================================================================================================== #

function get_N_l(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_l = get_N_l_helper(param_set, relaxation_timescale, q, T, ρ, w)
    if get_adjust_liq_N(relaxation_timescale)
        N_l = adjust_liq_N(param_set, N_l, q.liq; ρ = ρ, monodisperse = true, decrease_N_if_subsaturated = false)
    end
    N_l = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    return N_l
end
# ======================================================================================================================================== #

function get_N_l_raw_and_adjusted(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_l = get_N_l_helper(param_set, relaxation_timescale, q, T, ρ, w)
    if get_adjust_liq_N(relaxation_timescale)
        N_l_adjusted = adjust_liq_N(param_set, N_l, q.liq; ρ = ρ, monodisperse = true, decrease_N_if_subsaturated = false)
    end
    N_l_raw = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    N_l_adjusted = clamp(N_l_adjusted, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)

    return (; N_l_raw, N_l_adjusted)
end
# ======================================================================================================================================== #

function get_N_l_raw(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_l = get_N_l_helper(param_set, relaxation_timescale, q, T, ρ, w)
    N_l = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    return N_l
end
# ======================================================================================================================================== #



# --------------------------------------- #

function get_N_i(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

function get_N_i_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_i_and_N_i_no_boost(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end


function get_N_i_raw_and_adjusted(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

function get_N_i_raw_and_adjusted_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_i_raw_and_adjusted_and_N_i_no_boost(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

function get_N_i_raw(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_i_raw(param_set, relaxation_timescale, q, T, ρ, w)
end 

get_N_l(param_set::APS, relaxation_timescale::RelaxToEquilibrium, ts::TD.ThermodynamicState, w::FT) where {FT} = FT(NaN)
get_N_l_raw_and_adjusted(param_set::APS, relaxation_timescale::RelaxToEquilibrium, ts::TD.ThermodynamicState, w::FT) where {FT} = (N_l_raw = FT(NaN), N_l_adjusted = FT(NaN))
get_N_l_raw(param_set::APS, relaxation_timescale::RelaxToEquilibrium, ts::TD.ThermodynamicState, w::FT) where {FT} = FT(NaN)

function get_N_l(param_set::APS, relaxation_timescale::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
end

function get_N_l_raw_and_adjusted(param_set::APS, relaxation_timescale::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_l_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w)
end

function get_N_l_raw(param_set::APS, relaxation_timescale::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_l_raw(param_set, relaxation_timescale, q, T, ρ, w)
end

get_N_i(param_set::APS, Ni::FT, ts::TD.ThermodynamicState, w::FT) where {FT} = Ni
get_N_l(param_set::APS, Nl::FT, ts::TD.ThermodynamicState, w::FT) where {FT} = Nl
get_N_i(param_set::APS, Ni::FT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = Ni
get_N_l(param_set::APS, Nl::FT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = Nl

# ======================================================================================================================================== #

function get_N_i_raw(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_4 = relaxation_timescale.c_4i
    c_5 = relaxation_timescale.c_5i
    N_i_T = c_4 * exp(c_5 * (T - T_fr))
    N_i_T = clamp(N_i_T, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return N_i_T
end

function get_N_i_raw_and_adjusted(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    N_i = get_N_i_raw(param_set, relaxation_timescale, q, T, ρ, w)
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = (relaxation_timescale isa INP_Aware_Timescale)
        N_i_adjusted = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    else
        N_i_adjusted = N_i
    end
    N_i_raw = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_adjusted = clamp(N_i_adjusted, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return (; N_i_raw, N_i_adjusted)
end

function get_N_i_raw_and_adjusted_and_N_i_no_boost(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    N_i = get_N_i_raw(param_set, relaxation_timescale, q, T, ρ, w)
    N_i_no_boost = N_i
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = (relaxation_timescale isa INP_Aware_Timescale)
        N_i_adjusted = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        if apply_massflux_boost
            N_i_no_boost = adjust_ice_N(param_set, N_i_no_boost, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = false, apply_sedimentation_boost = false)
        else
            N_i_no_boost = N_i
        end
    else
        N_i_adjusted = N_i
    end
    N_i_raw = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_adjusted = clamp(N_i_adjusted, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_no_boost = clamp(N_i_no_boost, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    return (; N_i_raw, N_i_adjusted, N_i_no_boost)
end

# ======================================================================================================================================== #

function get_INP_concentration(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT, S_i::FT; apply_supersaturation_scaling::Bool=true) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    # return FT(NaN) # for ansatz where we don't actually predict INP... needed for say ice nucleation.
    N_INP = get_N_i_Cooper_curve(T; clamp_N=true)
    if apply_supersaturation_scaling
        r0_activation = FT(10e-6) # 10 μm particle post activation, matches https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L1089
        m0_activation = particle_mass(param_set, ice_type, r0_activation)
        N_INP *= INP_supersaturation_scaling_factor(param_set, S_i, q.ice, N_INP * m0_activation)
    end
    return N_INP
end
function get_INP_concentration(param_set::APS, relaxation_timescale::INP_Aware_Timescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT, S_i::FT; apply_supersaturation_scaling::Bool=true) where {FT}
    N_INP = get_N_i_raw(param_set, relaxation_timescale, q, T, ρ, w)
    if apply_supersaturation_scaling
        r0_activation = FT(10e-6) # 10 μm particle post activation, matches
        m0_activation = particle_mass(param_set, ice_type, r0_activation)
        N_INP *= INP_supersaturation_scaling_factor(param_set, S_i, q.ice, N_INP * m0_activation)
        if !isfinite(N_INP)
            error("N_INP = $N_INP not finite, N_INP_raw was $(get_N_i_raw(param_set, relaxation_timescale, q, T, ρ, w)); scaling_factor was $(INP_supersaturation_scaling_factor(param_set, S_i, q.ice, N_INP * m0_activation)); inputs were T=$T, ρ=$ρ, w=$w, S_i=$S_i, q_ice=$(q.ice)")
        end
    end
    return N_INP
end

"""
    INP_supersaturation_scaling_factor(param_set, S_i, q_i, q_ice_threshold; kwargs...)

Returns a scaling factor [0, 1] for INP concentration by blending two physical constraints:
1. **Supersaturation Activation**: INPs activate as `S_i` exceeds `S_ice_min_activation` (from `param_set`).
2. **Ice Feedback Masking**: If cloud ice (`q_i`) is already present, the dependence on `S_i` is reduced. 
   As `q_i` approaches `q_ice_threshold`, the factor is forced to 1.0 regardless of supersaturation.

# Arguments
- `param_set`: Struct containing `user_params.S_ice_min_activation`.
- `S_i`: Current ice supersaturation ratio (or supersaturation, depending on config).
- `q_i`: Current cloud ice specific humidity.
- `q_ice_threshold`: The ice concentration at which the scaling factor is forced to 1.0.

# Example Behavior
(Assuming `S_min = 0.1` and `q_thresh = 1e-4`)

| S_i  | q_i    | Result | Status              |
|:---- |:------ |:------ |:------------------- |
| 0.00 | 0.0    | 0.00   | Inactive            |
| 0.12 | 0.0    | 1.00   | Active (by S_i)     |
| 0.00 | 5e-5   | 0.50   | Lifted (by Ice)     |
| 0.00 | 1e-4   | 1.00   | Masked (Full Ice)   |
"""
function INP_supersaturation_scaling_factor(
    param_set::APS, 
    S_i::FT,
    q_i::FT,
    q_ice_threshold::FT; 
    smooth_fcn::Bool=true, 
    cutoff::FT = FT(0.9),
    transition_width_param::FT = FT(9.0)
    ) where {FT}

    S_min = FT(param_set.user_params.S_ice_min_activation)

    # --- 1. S_i Activation Scaling ---
    # SHORTCUT: If activation threshold is 0, use a step function
    if S_min ≤ zero(FT)
        # Safety check: Parameter cannot be negative
        if (S_min < zero(FT))
            error("S_ice_min_activation cannot be negative")
        end
        
        # Step function: 1.0 if supersaturated, 0.0 if subsaturated
        scaling_factor = (S_i ≥ zero(FT)) ? one(FT) : zero(FT)

    else
        # STANDARD: Calculate smooth transition or linear clamp

        # --- 1. S_i Activation Scaling ---
        S_start = cutoff * S_min
        S_width = S_min - S_start 

        scaling_factor = FT(0)

        if !smooth_fcn
            scaling_factor = clamp((S_i - S_start) / S_width, FT(0), FT(1))
        else
            S_mid = (S_start + S_min) / FT(2)
            k_steep = transition_width_param / S_width
            scaling_factor = FT(1) / (FT(1) + exp(-k_steep * (S_i - S_mid)))
        end
    end

    # --- 2. q_i Masking Factor ---
    q_mask = FT(0)

    # CHECK 1: If we exceeded the threshold (and threshold is not Inf)
    if q_i >= q_ice_threshold
        q_mask = FT(1)
        
    # CHECK 2: Only calculate smooth curve if we have ice AND a valid finite threshold
    elseif q_i > FT(0) && isfinite(q_ice_threshold)
        if !smooth_fcn
            q_mask = q_i / q_ice_threshold
        else
            k_steep_q = transition_width_param / q_ice_threshold
            q_mid = q_ice_threshold / FT(2)
            q_mask = FT(1) / (FT(1) + exp(-k_steep_q * (q_i - q_mid)))
        end
    end

    # --- 3. Blend ---
    scaling_factor = scaling_factor + (FT(1) - scaling_factor) * q_mask

    return min(scaling_factor, FT(1))
end


# ======================================================================================================================================== #
# Gradients
# ---------------- #

function get_dNINP_dz(param_set::APS, relaxation_timescale::RTT, T::FT, dTdz::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    dNdT =  get_dNINP_dT(param_set, relaxation_timescale, T)
    dNdz = dNdT * dTdz
    return dNdz
end
get_dNINP_dz(T::FT, dTdz::FT) where {FT} = get_dNINP_dT(T) * dTdz

# Base form using Cooper curve
function get_dNINP_dT(T::FT) where {FT}
    return get_d_N_i_Cooper_curve_dT(T)
end
get_dNINP_dT(param_set::APS, relaxation_timescale::RTT, T::FT) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}} = get_dNINP_dT(T) # Fallback, definitely should be true for NN, but LC and Base are weakpoints, and for arbitrary NNs is less defined.


# exponential
function get_dNINP_dT(param_set::APS, relaxation_timescale::Union{ExponentialTScalingIceRelaxationTimescale, GeometricLiqExponentialTScalingIceRelaxationTimescale}, T::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    # N_i_T = c_1 * exp(c_2 * (T - T_fr))
    # dNdT = c_2 * N_i_T
    return c_2 * c_1 * exp(c_2 * (T - T_fr))
end

function get_dNINP_dT(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, T::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_4 = relaxation_timescale.c_4i
    c_5 = relaxation_timescale.c_5i
    # N_i_T = c_4 * exp(c_5 * (T - T_fr))
    # dNdT = c_5 * N_i_T
    return c_5 * c_4 * exp(c_5 * (T - T_fr))
end

# powerlaw
function get_dNINP_dT(param_set::APS, relaxation_timescale::Union{PowerlawTScalingIceRelaxationTimescale, GeometricLiqPowerlawTScalingIceRelaxationTimescale}, T::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    
    # if T >= T_fr
    #     N_i = FT(0) # cant do powerlaw otherwise...
    # else
    #     N_i = (10^c_1) * (-(T - T_fr))^c_2
    # end
    
    if T >= T_fr
        return FT(0)
    else
        return -(c_2) * (10^c_1) * (-(T - T_fr))^(c_2 - 1)
    end
end
# ======================================================================================================================================== #


function get_ice_mult_ICNC_max(param_set::APS, N_INP::FT, N_i::FT, q_i::FT, q_r::FT, q_s::FT, T::FT, ρ::FT) where {FT}



    microphys_params = TCP.microphysics_params(param_set)
    # get <r_rain>, assume μ = 0
    μ_rai = FT(0)
    λ_r = CM1.lambda(microphys_params, rain_type, q_r, ρ)
    r_rain = (μ_rai + FT(1)) / λ_r

    # Droplet Shattering
    # Assume an enhancement of (r_rai / 400 μm) * 10^3 * K_DS, where K_DS is a gaussian that peaks at 0.1 at 257.5 K and is approimately 0 at 273.15K as per https://acp.copernicus.org/articles/21/15115/2021/#&gid=1&pid=1
    # N_DS_New = N_r * f_DS * K_DS # this wouldn't necessarily be the final enhancement, so we use empirical evidence
    K_DS = 0.1 * exp(-0.5 * ((T - 257.5) / 7.5)^2)
    f_DS = FT(1) + K_DS * (r_rain / FT(300e-6)) * 1e3 # based on https://www.pnas.org/doi/10.1073/pnas.2021387118

    f_DS /= 100 # reduce weight relative to HM, because we don't see any DS in M2005

    # Hallet-Mossop
    f_HM = begin
        T_fr = TCP.T_freeze(param_set)
        T_C = T - T_fr
        if (T_C < FT(-8)) || (T_C > FT(-3))
            f_HM = FT(0)
        elseif FT(-8) <= T_C <= FT(-5)
            f_HM = (T_C + 8) / 3
        elseif FT(-5) < T_C <= FT(-3)
            f_HM = (T_C + 3) / -2
        end
    end
    K_HM = FT(3.5e8) # Number of fragments per collision.
    # Because we don't track ICNC, rather than trying to track the collision rate directly, we'll assume we can generate an enhancement max based on available q_r, q_s, q_i
    # Assume the enhancement can reach 1e3 if we have enough liquid to do it. Enough means, assuming each rain droplet makes 3.5e8 fragments, we can get to N_INP = 10^3 * N_INP

    f_HM *= min(10^3, (K_HM * q_r * ρ + N_INP) / N_INP)
    ICNC_max = N_INP * max(f_HM, f_DS)

    # To avoid suddent jumps in N_i from using this, we should note that inferred N_i multiplication shouldn't be high if q_i and existing N_i are not high. So taking N_i from the last timestep, we don't let ICNC_max be more than 2x that, assuming 2*N_i is more than N_INP
    # This also helps avoid bad feedbacks w/ threshold acnv, where we may drastically swing τ_i by changing N_i too much, which feeds back with rain coming and going, which suddenly change the multiplication factor for everywhere underneath by up to 1000x...
    if !isnan(N_i) && (N_i > 0)

        r_acnv_scaling_factor = param_set.user_params.r_ice_acnv_scaling_factor # this MUST be less than 1!!!
        r_thresh = get_r_cond_precip(param_set, ice_type) * FT(param_set.user_params.r_ice_snow_threshold_scaling_factor)
        r_i_acnv = r_ice_acnv(param_set, FT(r_acnv_scaling_factor)) # this is the radius of the ice crystal at the acnv radius
        μ = μ_from_qN(param_set, ice_type, q_i, N_i; ρ=ρ) # this is the factor by which we scale the mean radius to get the mean radius for the ice crystals, so that we can use it in the N_i and N_l calculations
        N_i_acnv = N_from_qr(param_set, ice_type, q_i, r_i_acnv; monodisperse = false, μ=μ, ρ=ρ)
        N_thresh = N_from_qr(param_set, ice_type, q_i, r_thresh; monodisperse = false, μ=μ, ρ=ρ)
        N_i = clamp(N_i, N_thresh, N_i_acnv) # since this is the old N_i, we still. clamp it... multiplication shouldn't be able to raise N_i above N_i_acnv either, since it's just some fragments leaving while the rest remain...

        ICNC_max = min(ICNC_max, max( FT(1 + 1/3600) * N_i, FT(1+1/3600) * N_INP))
        # ICNC_max = clamp(ICNC_max, N_thresh, N_i_acnv) # don't allow ICNC_max to be outside physical bounds either [[ i think this enforces being at N_thresh too strongly maybe? like if mult isn't large enough we stay at thresh and never pitosn... ]]
        ICNC_max = clamp(ICNC_max, N_thresh/8, N_i_acnv) # allow to go to N_thresh /8 (r_thresh*2) to enable pitosn to do acnv.
    else
        # @warn "Got here from inputs N_INP = $N_INP; N_i = $N_i; q_i = $q_i; q_r = $q_r; q_s = $q_s; T = $T; ρ = $ρ; ICNC_max = $ICNC_max; f_HM = $f_HM; f_DS = $f_DS; r_rain = $r_rain; K_DS = $K_DS; f_DS = $f_DS"
        # e.g. Equilibrium N_i  is nan
        # ICNC_max = N_INP # no mult if we don't estimate N_i, bc we'll just get N from q.
    end

    return ICNC_max
end
get_ice_mult_factor_ICNC_max(param_set::APS, N_INP::FT, N_i::FT, q_i::FT, q_r::FT, q_s::FT, T::FT, ρ::FT) where {FT} = max(FT(1), get_ice_mult_ICNC_max(param_set, N_INP, N_i, q_i, q_r, q_s, T, ρ) / N_INP)


# ======================================================================================================================================== #

get_N_i_helper(param_set::APS, relaxation_timescale::RelaxToEquilibrium, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = get_N_i_Cooper_curve(T; clamp_N=true)
get_N_l_helper(param_set::APS, relaxation_timescale::RelaxToEquilibrium, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# :Base
# get_N_i_helper(param_set::APS, relaxation_timescale::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)   
get_N_i_helper(param_set::APS, relaxation_timescale::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = get_N_i_Cooper_curve(T; clamp_N=true) # trial using this to force PITOSN
# get_N_l_helper(param_set::APS, relaxation_timescale::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)


# :exponential_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale::ExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)

    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    return FT(c_1 * exp(c_2 * (T - T_fr)))
end
get_N_l_helper(param_set::APS, relaxation_timescale::ExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# :exponential_T_scaling_ice_raw
get_N_i_helper(param_set::APS, relaxation_timescale::ExponentialTScalingIceRawRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)
get_N_l_helper(param_set::APS, relaxation_timescale::ExponentialTScalingIceRawRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# :powerlaw_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale::PowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)

    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i

    if T >= T_fr
        N_i = FT(0) # cant do powerlaw otherwise...
    else
        N_i = (10^c_1) * (-(T - T_fr))^c_2
    end

    return N_i
end
get_N_l_helper(param_set::APS, relaxation_timescale::PowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT} =  FT(NaN)

# :exponential_times_powerlaw_scaling_ice
get_N_i_helper(param_set::APS, relaxation_timescale::ExponentialTimesPowerlawScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = error("NotImplmentedError: This relaxation_timescale functionality has not been implemented yet")
get_N_l_helper(param_set::APS, relaxation_timescale::ExponentialTimesPowerlawScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = error("NotImplmentedError: This relaxation_timescale functionality has not been implemented yet")


# :geometric_liq__geometric_ice
function get_N_i_helper(param_set::APS, relaxation_timescale::GeometricLiqGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    c_3 = relaxation_timescale.c_3i
    c_4 = relaxation_timescale.c_4i
    return min( FT(c_1 * (q.ice/1e-7)^(c_2) + 10^c_3), FT(10^c_3 + 10^c_4))
end
function get_N_l_helper(param_set::APS, relaxation_timescale::GeometricLiqGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale.c_1l
    c_2 = relaxation_timescale.c_2l
    c_3 = relaxation_timescale.c_3l
    c_4 = relaxation_timescale.c_4l
    return min(FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3), FT(10^c_3 + 10^c_4))
end

# :geometric_liq__exponential_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    return FT(c_1 * exp(c_2 * (T - T_fr)))
end
function get_N_l_helper(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale.c_1l
    c_2 = relaxation_timescale.c_2l
    c_3 = relaxation_timescale.c_3l
    c_4 = relaxation_timescale.c_4l
    return min( FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3) , FT(10^c_3 + 10^c_4))
end

# :geometric_liq__powerlaw_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale::GeometricLiqPowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    if T >= T_fr
        N_i = FT(0) # cant do powerlaw otherwise...
    else
        N_i = (10^c_1) * (-(T - T_fr))^c_2
    end
    return N_i
end
function get_N_l_helper(param_set::APS, relaxation_timescale::GeometricLiqPowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale.c_1l
    c_2 = relaxation_timescale.c_2l
    c_3 = relaxation_timescale.c_3l
    c_4 = relaxation_timescale.c_4l
    return min(FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3), FT(10^c_3 + 10^c_4))
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
# function get_N_i_helper(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}
#     T_fr = TCP.T_freeze(param_set)
#     c_1 = relaxation_timescale.c_1i
#     c_2 = relaxation_timescale.c_2i
#     c_3 = relaxation_timescale.c_3i
#     c_4 = relaxation_timescale.c_4i
#     return (c_1 * (q.ice/1e-7)^(c_2) + c_3) * exp(c_4 * (T - T_fr))
# end

function get_N_i_helper(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    # c_3 = relaxation_timescale.c_3i
    c_4 = relaxation_timescale.c_4i
    c_5 = relaxation_timescale.c_5i
    N_i_T = c_4 * exp(c_5 * (T - T_fr))
    N_i_q = c_1 * (q.ice/1e-7)^(c_2) #+ (10^c_3) * N_i_T
    return min(N_i_q, N_i_T)
end

function get_N_l_helper(param_set::APS, relaxation_timescale::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale.c_1l
    c_2 = relaxation_timescale.c_2l
    c_3 = relaxation_timescale.c_3l
    c_4 = relaxation_timescale.c_4l
    return min( FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3) , FT(10^c_3 + 10^c_4))
end

# :linear_combination
function get_N_i_helper(param_set::APS, relaxation_timescale::LinearCombinationRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    c_3 = relaxation_timescale.c_3i
    N_i = exp(c_1 + c_2 * (T - T_fr) + c_3 * q.ice / 1e-7)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice/2, relaxation_timescale.args.max_N_ice*2) # preclamp bc linearcombo sometimes goes to inf...
    return N_i

end
function get_N_l_helper(param_set::APS, relaxation_timescale::LinearCombinationRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale.c_1l
    c_2 = relaxation_timescale.c_2l
    c_3 = relaxation_timescale.c_3l
    return exp(c_1 + c_2 * (T - T_fr) + c_3 * q.liq )
end

# :linear_combination_with_w
function get_N_i_helper(param_set::APS, relaxation_timescale::LinearCombinationWithWRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    w_0 = FT(1e-3) # 1 mm/s
    c_1 = relaxation_timescale.c_1i
    c_2 = relaxation_timescale.c_2i
    c_3 = relaxation_timescale.c_3i
    c_4 = relaxation_timescale.c_4i
    N_i = exp(c_1 + c_2 * (T - T_fr) + c_3 * (q.ice / 1e-7) + c_4 * w / w_0)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice/2, relaxation_timescale.args.max_N_ice*2) # preclamp bc linearcombo sometimes goes to inf...
    return N_i
end
function get_N_l_helper(param_set::APS, relaxation_timescale::LinearCombinationWithWRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    w_0 = FT(1e-3) # 1 mm/s
    c_1 = relaxation_timescale.c_1l
    c_2 = relaxation_timescale.c_2l
    c_3 = relaxation_timescale.c_3l
    c_4 = relaxation_timescale.c_4l
    return exp(c_1 + c_2 * (T - T_fr) + c_3 * q.liq + c_4 * w / w_0)
end

# :neural_network & :neural_network_no_weights & :neural_network_random_init & neural_network_pca_noise
function get_N_i_helper(param_set::APS, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale
    _, _, _, N_i = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?
    
    return N_i
end

function get_N_l_helper(param_set::APS, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale
    _, _, N_l, _ = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?
    
    return N_l
end

# # :raymond_ice_test
# get_N_i_helper(param_set::APS, relaxation_timescale::RaymondIceTestRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)
# get_N_l_helper(param_set::APS, relaxation_timescale::RaymondIceTestRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# ======================================================================================================================================== #


function get_Ns(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT = FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_l::FT = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    N_i::FT = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    return (; N_liq = N_l, N_ice = N_i)
end

function get_Ns_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT = FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    N_l::FT = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    N_i::FT = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    if apply_massflux_boost
        N_i_no_boost = get_N_i(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = FT(0), dTdz = dTdz, w_i = w_i, apply_massflux_boost = false, apply_sedimentation_boost = apply_sedimentation_boost)
    else 
        N_i_no_boost = N_i
    end
    return (; N_liq = N_l, N_ice = N_i, N_ice_no_boost = N_i_no_boost)
end

function get_Ns_raw_and_adjusted(param_set::APS, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT = FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    (; N_i_raw, N_i_adjusted) = get_N_i_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    N_l::FT = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    return (; N_liq = N_l, N_ice_raw = N_i_raw, N_ice_adjusted = N_i_adjusted)
end

function get_Ns_raw_and_adjusted_and_N_i_no_boost(param_set::APS, relaxation_timescale::AbstractRelaxationTimescaleType, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT = FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    (; N_i_raw, N_i_adjusted, N_i_no_boost) = get_N_i_raw_and_adjusted_and_N_i_no_boost(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    N_l::FT = get_N_l(param_set, relaxation_timescale, q, T, ρ, w)
    return (; N_liq = N_l, N_ice_raw = N_i_raw, N_ice_adjusted = N_i_adjusted, N_ice_no_boost = N_i_no_boost)
end

function get_Ns(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_Ns(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

function get_Ns_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_Ns_and_N_i_no_boost(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

function get_Ns_raw_and_adjusted(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_Ns_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

function get_Ns_raw_and_adjusted_and_N_i_no_boost(param_set::APS, relaxation_timescale::RTT, ts::TD.ThermodynamicState, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_Ns_raw_and_adjusted_and_N_i_no_boost(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
end

# --------------------------------------------------------------------------- #

function get_τs_and_Ns(param_set::APS, microphys_params::ACMP, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    # Ns = get_Ns(param_set,                   relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    Ns = get_Ns_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    if apply_massflux_boost || apply_sedimentation_boost
        Ns_no_boost = get_Ns_raw_and_adjusted(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = false, apply_sedimentation_boost = apply_sedimentation_boost) # No massflux boost for τs
    else
        Ns_no_boost = Ns
    end
    # Don't use boosted values for τs
    # τs = get_τs(param_set, microphys_params, relaxation_timescale, q, T, p, ρ, w; N_l = Ns_no_boost.N_liq, N_i_raw = Ns_no_boost.N_ice_raw, N_i_adjusted = Ns_no_boost.N_ice_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    τs = get_τs(param_set, microphys_params, relaxation_timescale, q, T, p, ρ, w; N_l = Ns_no_boost.N_liq, N_i_raw = Ns_no_boost.N_ice_adjusted, N_i_adjusted = Ns_no_boost.N_ice_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost) # added 2026-01-14 :: I think our raw shouldn't be pure N_ice_raw anymore... it should be adjusted if adjust_ice_N is on... Now we have activation we dont need a shadow realm of tiny particles esp since we have supersaturation gates...
    return (; τ_liq = τs.τ_liq, τ_ice = τs.τ_ice, N_liq = Ns.N_liq, N_ice = Ns.N_ice_adjusted)

    # Ns = get_Ns(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    # return (; τ_liq = τs.τ_liq, τ_ice = τs.τ_ice, N_liq = Ns.N_liq, N_ice = Ns.N_ice)
end


function get_τs_and_Ns_and_N_i_no_boost(param_set::APS, microphys_params::ACMP, relaxation_timescale::RTT, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false, use_boost_for_τ::Bool=true) where {FT, RTT <: Union{RelaxToEquilibrium, AbstractRelaxationTimescaleType}}
    # Ns = get_Ns(param_set,                   relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    Ns = get_Ns_raw_and_adjusted_and_N_i_no_boost(param_set, relaxation_timescale, q, T, ρ, w; N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)

    # Don't use boosted values for τs, also the raw value should just be w/o the boost? (we assume the remaining particles are in snow... so we don't use N_i_raw....
    N_i_adjusted = use_boost_for_τ ? Ns.N_ice_adjusted : Ns.N_ice_no_boost
    τs = get_τs(param_set, microphys_params, relaxation_timescale, q, T, p, ρ, w; N_l = Ns.N_liq, N_i_raw = Ns.N_ice_adjusted, N_i_adjusted = N_i_adjusted, N_INP_top = N_INP_top, f_ice_mult = f_ice_mult, q_sno = q_sno, massflux = massflux, dTdz = dTdz, w_i = w_i, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    return (; τ_liq = τs.τ_liq, τ_ice = τs.τ_ice, N_liq = Ns.N_liq, N_ice = Ns.N_ice_adjusted, N_ice_no_boost = Ns.N_ice_no_boost)
end




# -- add fallback for NN that will only call the network once.... -- #
function get_Ns(param_set::APS, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale
    _, _, N_l, N_i = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?


    if get_adjust_liq_N(relaxation_timescale)
        # microphys_params::ACMP = TCP.microphysics_params(param_set)
        N_l = adjust_liq_N(param_set, N_l, q.liq; ρ = ρ, monodisperse = true, decrease_N_if_subsaturated = false)
    end
    if get_adjust_ice_N(relaxation_timescale)
        # microphys_params = TCP.microphysics_params(param_set) # hopefully the compiler optimizes this out if it's loaded earlier
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = false
        N_i = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    end

    # because this bypasses the normal get_N_i(, get_N_l(), and get_τs() we need to clamp the values here
    N_l = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    
    return (; N_liq = N_l, N_ice = N_i)
end

function get_Ns_and_N_i_no_boost(param_set::APS, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT = FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale
    _, _, N_l, N_i = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?


    if get_adjust_liq_N(relaxation_timescale)
        # microphys_params::ACMP = TCP.microphysics_params(param_set)
        N_l = adjust_liq_N(param_set, N_l, q.liq; ρ = ρ, monodisperse = true, decrease_N_if_subsaturated = false)
    end
    N_i_no_boost = N_i
    
    if get_adjust_ice_N(relaxation_timescale)
        # microphys_params = TCP.microphysics_params(param_set) # hopefully the compiler optimizes this out if it's loaded earlier
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz) # should this be limited to only when we have ice supersat?
        N_i_from_INP = false
        N_i = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        if apply_massflux_boost
            N_i_no_boost = adjust_ice_N(param_set, N_i_no_boost, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = FT(0), dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = false, apply_sedimentation_boost = apply_sedimentation_boost) # no massflux boost for no boost value
        else
            N_i_no_boost = N_i
        end
    end

    # because this bypasses the normal get_N_i(, get_N_l(), and get_τs() we need to clamp the values here
    N_l = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_no_boost = clamp(N_i_no_boost, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    
    return (; N_liq = N_l, N_ice = N_i, N_ice_no_boost = N_i_no_boost)
end


function get_τs_and_Ns(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale
    τ_liq, τ_ice, N_l, N_i = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?

    if get_adjust_liq_N(relaxation_timescale)
        N_l = adjust_liq_N(param_set, N_l, q.liq; ρ = ρ, monodisperse = true, decrease_N_if_subsaturated = false)
    end
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz)
        N_i_from_INP = false
        N_i = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, w_i = w_i, dNINP_dz = dNINP_dz, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    end
    
    # because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here
    τ_liq = clamp(τ_liq, relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    τ_ice = clamp(τ_ice, relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)
    N_l = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)

    return (; τ_liq = τ_liq, τ_ice = τ_ice, N_liq = N_l, N_ice = N_i)
end

function get_τs_and_Ns_and_N_i_no_boost(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT; N_INP_top::FT = FT(NaN), f_ice_mult::FT=FT(1), q_sno::FT=FT(0), massflux::FT=FT(0), dTdz::FT=FT(0), w_i::FT=FT(0), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where {FT}
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale
    τ_liq, τ_ice, N_l, N_i = predict_τ(ρ, T, q.liq, q.ice, w, (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?

    if get_adjust_liq_N(relaxation_timescale)
        N_l = adjust_liq_N(param_set, N_l, q.liq; ρ = ρ, monodisperse = true, decrease_N_if_subsaturated = false)
    end
    N_i_no_boost = N_i
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        N_INP = get_INP_concentration(param_set, relaxation_timescale, q, T, ρ, w, S_i)
        dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T, dTdz)
        N_i_from_INP = false
        N_i = adjust_ice_N(param_set, N_i, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = massflux, w_i = w_i, dNINP_dz = dNINP_dz, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
        if apply_massflux_boost
            N_i_no_boost = adjust_ice_N(param_set, N_i_no_boost, N_INP * f_ice_mult, q.ice; ρ = ρ, S_i = S_i, monodisperse = false, decrease_N_if_subsaturated = true, N_INP_top = N_INP_top * f_ice_mult, q_l = q.liq, q_s = q_sno, massflux = FT(0), w_i = w_i, dNINP_dz = dNINP_dz, N_i_from_INP = N_i_from_INP, apply_massflux_boost = false, apply_sedimentation_boost = apply_sedimentation_boost) # no massflux boost for no boost value
        else
            N_i_no_boost = N_i
        end
    end
    
    # because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here
    τ_liq = clamp(τ_liq, relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    τ_ice = clamp(τ_ice, relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)
    N_l = clamp(N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    N_i = clamp(N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    N_i_no_boost = clamp(N_i_no_boost, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)

    return (; τ_liq = τ_liq, τ_ice = τ_ice, N_liq = N_l, N_ice = N_i, N_ice_no_boost = N_i_no_boost)
end

# methods for dispatch because kwargs are not broadcast over which is a pain for views below... (bool is fine to keep, just not arrays we wannna broadcast)
adjust_ice_N_no_kwargs(param_set::APS, N_i::FT, N_INP::FT, q_ice::FT, ρ::FT, S_i::FT, q_liq::FT, q_sno::FT, massflux::FT, dNINP_dz::FT, w_i::FT, N_INP_top::FT; monodisperse::Bool, decrease_N_if_subsaturated::Bool, N_i_from_INP::Bool = false, apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false) where{FT} = adjust_ice_N(param_set, N_i, N_INP, q_ice; ρ = ρ, S_i = S_i, monodisperse = monodisperse, decrease_N_if_subsaturated = decrease_N_if_subsaturated, N_INP_top = N_INP_top, q_l = q_liq, q_s = q_sno, massflux = massflux, dNINP_dz = dNINP_dz, w_i = w_i, N_i_from_INP = N_i_from_INP, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
adjust_liq_N_no_kwargs(param_set::APS, N_l::FT, q_liq::FT, ρ::FT; monodisperse::Bool, decrease_N_if_subsaturated::Bool) where{FT} = adjust_liq_N(param_set, N_l, q_liq; ρ = ρ, monodisperse = monodisperse, decrease_N_if_subsaturated = decrease_N_if_subsaturated)

function get_τs_and_Ns!(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::CC.Fields.Field, T::CC.Fields.Field, p::CC.Fields.Field, ρ::CC.Fields.Field, w::CC.Fields.Field, area::CC.Fields.Field, τ_liq::CC.Fields.Field, τ_ice::CC.Fields.Field, N_liq::CC.Fields.Field, N_ice::CC.Fields.Field, f_ice_mult::CC.Fields.Field, q_sno::CC.Fields.Field, massflux::CC.Fields.Field, dTdz::CC.Fields.Field, w_i::CC.Fields.Field; N_INP_top = eltype(param_set)(NaN), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false)
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale

    # NOTE: parent(q) conveniently becomes a n x 3 Matrix (3 for [tot, liq, ice]), so q_liq = parent(q)[..., 2] and q_ice = parent(q)[..., 3]

    valid_inds = vec(parent(area) .> 0) # get the valid indices where area is greater than 0, convert to vec from nx1 matrix
    _τ_liq, _τ_ice, _N_l, _N_i = predict_τ(parent(ρ)[valid_inds], parent(T)[valid_inds], parent(q)[valid_inds, 2], parent(q)[valid_inds, 3], parent(w)[valid_inds], (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?
    
    # Because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here
    # _τ_liq .= clamp(_τ_liq, relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    # _τ_ice .= clamp(_τ_ice, relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)
    # _N_l .= clamp(_N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    # _N_i .= 


    parent(τ_liq)[valid_inds] .= _τ_liq
    parent(τ_ice)[valid_inds] .= _τ_ice
    parent(N_liq)[valid_inds] .= _N_l
    parent(N_ice)[valid_inds] .= _N_i

    if get_adjust_liq_N(relaxation_timescale) || get_adjust_ice_N(relaxation_timescale)
        q_p      = parent(q)
        ρ_p      = parent(ρ)
        T_p      = parent(T)
        w_p      = parent(w)
        dTdz_p   = parent(dTdz)
        N_liq_p  = parent(N_liq)
        N_ice_p  = parent(N_ice)
        f_mult_p = parent(f_ice_mult)
        q_sno_p  = parent(q_sno)
        mf_p     = parent(massflux)
        w_i_p    = parent(w_i)
    end


    # Because  this bypasses get_N_l(), get_N_i(), we need to apply adjust_liq/ice_N() here.
    if get_adjust_liq_N(relaxation_timescale)
        # @view(parent(N_liq)[valid_inds]) .= adjust_liq_N_no_kwargs.(param_set, @view(parent(N_liq)[valid_inds]), @view(parent(q)[valid_inds, 2]), @view(parent(ρ)[valid_inds]); monodisperse=true, decrease_N_if_subsaturated=false)
        for i in eachindex(valid_inds)
            valid_inds[i] || continue
            N_liq_p[i] = adjust_liq_N_no_kwargs(param_set, N_liq_p[i], q_p[i, 2], ρ_p[i]; monodisperse=true, decrease_N_if_subsaturated=false)
        end
    end
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        # FT = eltype(param_set)
        # qs = TD.PhasePartition{FT}[TD.PhasePartition(parent(q)[i, 1], parent(q)[i, 2], parent(q)[i, 3]) for i in eachindex(valid_inds) if valid_inds[i]]
        # S_i = TD.supersaturation.(thermo_params, qs, @view(parent(ρ)[valid_inds]), @view(parent(T)[valid_inds]), TD.Ice())
        # N_INP = get_INP_concentration.(param_set, relaxation_timescale, qs, @view(parent(T)[valid_inds]), @view(parent(ρ)[valid_inds]), @view(parent(w)[valid_inds]), S_i)
        # dNINP_dz = get_dNINP_dz.(param_set, relaxation_timescale, @view(parent(T)[valid_inds]), @view(parent(dTdz)[valid_inds]))
        # @view(parent(N_ice)[valid_inds]) .= adjust_ice_N_no_kwargs.(param_set,
        #     @view(parent(N_ice)[valid_inds]),
        #     N_INP .* @view(parent(f_ice_mult)[valid_inds]),
        #     @view(parent(q)[valid_inds, 3]),
        #     @view(parent(ρ)[valid_inds]),
        #     S_i,
        #     @view(parent(q)[valid_inds, 2]), # q_l for boosting <r>
        #     @view(parent(q_sno)[valid_inds]),
        #     @view(parent(massflux)[valid_inds]),
        #     dNINP_dz,
        #     @view(parent(w_i)[valid_inds]),
        #     N_INP_top .* @view(parent(f_ice_mult)[valid_inds]),
        #     ; monodisperse=false, decrease_N_if_subsaturated=true, N_i_from_INP = false, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)

        q_p      = parent(q)
        ρ_p      = parent(ρ)
        T_p      = parent(T)
        w_p      = parent(w)
        dTdz_p   = parent(dTdz)
        N_ice_p  = parent(N_ice)
        f_mult_p = parent(f_ice_mult)
        q_sno_p  = parent(q_sno)
        mf_p     = parent(massflux)
        w_i_p    = parent(w_i)
        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue

            qs = TD.PhasePartition(q_p[i,1], q_p[i,2], q_p[i,3])
            S_i = TD.supersaturation(thermo_params, qs, ρ_p[i], T_p[i], TD.Ice())
            N_INP = get_INP_concentration(param_set, relaxation_timescale, qs, T_p[i], ρ_p[i], w_p[i], S_i)
            dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T_p[i], dTdz_p[i])
            N_ice_p[i] = adjust_ice_N_no_kwargs(
                param_set,
                N_ice_p[i],
                N_INP * f_mult_p[i],
                q_p[i, 3],
                ρ_p[i],
                S_i,
                q_p[i, 2],
                q_sno_p[i],
                mf_p[i],
                dNINP_dz,
                w_i_p[i],
                N_INP_top * f_mult_p[i];
                monodisperse=false,
                decrease_N_if_subsaturated=true,
                N_i_from_INP=false,
                apply_massflux_boost=apply_massflux_boost,
                apply_sedimentation_boost=apply_sedimentation_boost,
            )
        end

    end


    # Because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here [[ (clamping the static array didn't work) -- e.g. clamp(_N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice) return an error bc of the stride arrays ]]
    clamp!(@view(parent(τ_liq)[valid_inds]), relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    clamp!(@view(parent(τ_ice)[valid_inds]), relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)
    clamp!(@view(parent(N_liq)[valid_inds]), relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    clamp!(@view(parent(N_ice)[valid_inds]), relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)

    return 
end


function get_τs_and_Ns_and_N_i_no_boost!(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::CC.Fields.Field, T::CC.Fields.Field, p::CC.Fields.Field, ρ::CC.Fields.Field, w::CC.Fields.Field, area::CC.Fields.Field, τ_liq::CC.Fields.Field, τ_ice::CC.Fields.Field, N_liq::CC.Fields.Field, N_ice::CC.Fields.Field, N_ice_no_boost::CC.Fields.Field, f_ice_mult::CC.Fields.Field, q_sno::CC.Fields.Field, massflux::CC.Fields.Field, dTdz::CC.Fields.Field, w_i::CC.Fields.Field; N_INP_top = eltype(param_set)(NaN), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false)
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale

    # NOTE: parent(q) conveniently becomes a n x 3 Matrix (3 for [tot, liq, ice]), so q_liq = parent(q)[..., 2] and q_ice = parent(q)[..., 3]

    valid_inds = vec(parent(area) .> 0) # get the valid indices where area is greater than 0, convert to vec from nx1 matrix
    _τ_liq, _τ_ice, _N_l, _N_i = predict_τ(parent(ρ)[valid_inds], parent(T)[valid_inds], parent(q)[valid_inds, 2], parent(q)[valid_inds, 3], parent(w)[valid_inds], (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?
    
    # Because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here
    # _τ_liq .= clamp(_τ_liq, relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    # _τ_ice .= clamp(_τ_ice, relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)
    # _N_l .= clamp(_N_l, relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    # _N_i .= clamp(_N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)


    parent(τ_liq)[valid_inds] .= _τ_liq
    parent(τ_ice)[valid_inds] .= _τ_ice
    parent(N_liq)[valid_inds] .= _N_l
    parent(N_ice)[valid_inds] .= _N_i

    if get_adjust_liq_N(relaxation_timescale) || get_adjust_ice_N(relaxation_timescale)
        q_p      = parent(q)
        ρ_p      = parent(ρ)
        T_p      = parent(T)
        w_p      = parent(w)
        dTdz_p   = parent(dTdz)
        N_liq_p  = parent(N_liq)
        N_ice_p  = parent(N_ice)
        f_mult_p = parent(f_ice_mult)
        q_sno_p  = parent(q_sno)
        mf_p     = parent(massflux)
        w_i_p    = parent(w_i)
    end

    # Because  this bypasses get_N_l(), get_N_i(), we need to apply adjust_liq/ice_N() here.
    if get_adjust_liq_N(relaxation_timescale)
        # @view(parent(N_liq)[valid_inds]) .= adjust_liq_N_no_kwargs.(param_set, @view(parent(N_liq)[valid_inds]), @view(parent(q)[valid_inds, 2]), @view(parent(ρ)[valid_inds]); monodisperse=true, decrease_N_if_subsaturated=false)
        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue
            parent(N_liq)[i] = adjust_liq_N_no_kwargs(param_set, N_liq_p[i], q_p[i, 2], ρ_p[i]; monodisperse=true, decrease_N_if_subsaturated=false)
        end
    end
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        # FT = eltype(param_set)
        # qs = TD.PhasePartition{FT}[TD.PhasePartition(parent(q)[i, 1], parent(q)[i, 2], parent(q)[i, 3]) for i in eachindex(valid_inds) if valid_inds[i]]
        # S_i = TD.supersaturation.(thermo_params, qs, @view(parent(ρ)[valid_inds]), @view(parent(T)[valid_inds]), TD.Ice())
        # N_INP = get_INP_concentration.(param_set, relaxation_timescale, qs, @view(parent(T)[valid_inds]), @view(parent(ρ)[valid_inds]), @view(parent(w)[valid_inds]), S_i)
        # dNINP_dz = get_dNINP_dz.(param_set, relaxation_timescale, @view(parent(T)[valid_inds]), @view(parent(dTdz)[valid_inds]))

        # @view(parent(N_ice)[valid_inds]) .= adjust_ice_N_no_kwargs.(param_set,
        #     @view(parent(N_ice)[valid_inds]),
        #     N_INP .* @view(parent(f_ice_mult)[valid_inds]),
        #     @view(parent(q)[valid_inds, 3]),
        #     @view(parent(ρ)[valid_inds]),
        #     S_i,
        #     @view(parent(q)[valid_inds, 2]), # q_l for boosting <r>
        #     @view(parent(q_sno)[valid_inds]),
        #     @view(parent(massflux)[valid_inds]),
        #     dNINP_dz,
        #     @view(parent(w_i)[valid_inds]),
        #     N_INP_top .* @view(parent(f_ice_mult)[valid_inds]),
        #     ; monodisperse=false, decrease_N_if_subsaturated=true, N_i_from_INP = false, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)

        # if apply_massflux_boost
        #     @view(parent(N_ice_no_boost)[valid_inds]) .= adjust_ice_N_no_kwargs.(param_set,
        #         @view(parent(N_ice)[valid_inds]),
        #         N_INP .* @view(parent(f_ice_mult)[valid_inds]),
        #         @view(parent(q)[valid_inds, 3]),
        #         @view(parent(ρ)[valid_inds]),
        #         S_i,
        #         @view(parent(q)[valid_inds, 2]), # q_l for boosting <r>
        #         @view(parent(q_sno)[valid_inds]),
        #         @view(parent(massflux)[valid_inds]),
        #         dNINP_dz,
        #         @view(parent(w_i)[valid_inds]),
        #         N_INP_top .* @view(parent(f_ice_mult)[valid_inds]),
        #         ; monodisperse=false, decrease_N_if_subsaturated=true, N_i_from_INP = false, apply_massflux_boost = false, apply_sedimentation_boost = apply_sedimentation_boost)
        # else
        #     @view(parent(N_ice_no_boost)[valid_inds]) .= @view(parent(N_ice)[valid_inds])
        # end

        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue

            qs = TD.PhasePartition(q_p[i,1], q_p[i,2], q_p[i,3])
            S_i = TD.supersaturation(thermo_params, qs, ρ_p[i], T_p[i], TD.Ice())
            N_INP = get_INP_concentration(param_set, relaxation_timescale, qs, T_p[i], ρ_p[i], w_p[i], S_i)
            dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T_p[i], dTdz_p[i])
            N_ice_p[i] = adjust_ice_N_no_kwargs(
                param_set,
                N_ice_p[i],
                N_INP * f_mult_p[i],
                q_p[i, 3],
                ρ_p[i],
                S_i,
                q_p[i, 2],
                q_sno_p[i],
                mf_p[i],
                dNINP_dz,
                w_i_p[i],
                N_INP_top * f_mult_p[i];
                monodisperse=false,
                decrease_N_if_subsaturated=true,
                N_i_from_INP=false,
                apply_massflux_boost=apply_massflux_boost,
                apply_sedimentation_boost=apply_sedimentation_boost,
            )

            if apply_massflux_boost
                parent(N_ice_no_boost)[i] = adjust_ice_N_no_kwargs(
                    param_set,
                    N_ice_p[i],
                    N_INP * f_mult_p[i],
                    q_p[i, 3],
                    ρ_p[i],
                    S_i,
                    q_p[i, 2],
                    q_sno_p[i],
                    mf_p[i],
                    dNINP_dz,
                    w_i_p[i],
                    N_INP_top * f_mult_p[i];
                    monodisperse=false,
                    decrease_N_if_subsaturated=true,
                    N_i_from_INP=false,
                    apply_massflux_boost=false,
                    apply_sedimentation_boost=apply_sedimentation_boost,
                )
            else
                parent(N_ice_no_boost)[i] = N_ice_p[i]
            end
        end


    end


    # Because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here [[ (clamping the static array didn't work) -- e.g. clamp(_N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice) return an error bc of the stride arrays ]]
    clamp!(@view(parent(τ_liq)[valid_inds]), relaxation_timescale.args.min_τ_liq, relaxation_timescale.args.max_τ_liq)
    clamp!(@view(parent(τ_ice)[valid_inds]), relaxation_timescale.args.min_τ_ice, relaxation_timescale.args.max_τ_ice)
    clamp!(@view(parent(N_liq)[valid_inds]), relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    clamp!(@view(parent(N_ice)[valid_inds]), relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    clamp!(@view(parent(N_ice_no_boost)[valid_inds]), relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)

    return 
end


function get_Ns!(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::CC.Fields.Field, T::CC.Fields.Field, p::CC.Fields.Field, ρ::CC.Fields.Field, w::CC.Fields.Field, area::CC.Fields.Field, τ_liq::CC.Fields.Field, τ_ice::CC.Fields.Field, N_liq::CC.Fields.Field, N_ice::CC.Fields.Field, f_ice_mult::CC.Fields.Field, q_sno::CC.Fields.Field, massflux::CC.Fields.Field, dTdz::CC.Fields.Field, w_i::CC.Fields.Field; N_INP_top = eltype(param_set)(NaN), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false)
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale

    # NOTE: parent(q) conveniently becomes a n x 3 Matrix (3 for [tot, liq, ice]), so q_liq = parent(q)[..., 2] and q_ice = parent(q)[..., 3]

    valid_inds = vec(parent(area) .> 0) # get the valid indices where area is greater than 0, convert to vec from nx1 matrix
    _, _, _N_l, _N_i = predict_τ(parent(ρ)[valid_inds], parent(T)[valid_inds], parent(q)[valid_inds, 2], parent(q)[valid_inds, 3], parent(w)[valid_inds], (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?

    parent(N_liq)[valid_inds] .= _N_l
    parent(N_ice)[valid_inds] .= _N_i

    if get_adjust_liq_N(relaxation_timescale) || get_adjust_ice_N(relaxation_timescale)
        q_p      = parent(q)
        ρ_p      = parent(ρ)
        T_p      = parent(T)
        w_p      = parent(w)
        dTdz_p   = parent(dTdz)
        N_liq_p  = parent(N_liq)
        N_ice_p  = parent(N_ice)
        f_mult_p = parent(f_ice_mult)
        q_sno_p  = parent(q_sno)
        mf_p     = parent(massflux)
        w_i_p    = parent(w_i)
    end

    # Because  this bypasses get_N_l(), get_N_i(), we need to apply adjust_liq/ice_N() here.
    if get_adjust_liq_N(relaxation_timescale)
        # @view(parent(N_liq)[valid_inds]) .= adjust_liq_N_no_kwargs.(param_set, @view(parent(N_liq)[valid_inds]), @view(parent(q)[valid_inds, 2]), @view(parent(ρ)[valid_inds]); monodisperse=true, decrease_N_if_subsaturated=false)
        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue
            N_liq_p[i] = adjust_liq_N_no_kwargs(param_set, N_liq_p[i], q_p[i, 2], ρ_p[i]; monodisperse=true, decrease_N_if_subsaturated=false)
        end
    end
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        FT = eltype(param_set)
    #     qs = TD.PhasePartition{FT}[TD.PhasePartition(parent(q)[i, 1], parent(q)[i, 2], parent(q)[i, 3]) for i in eachindex(valid_inds) if valid_inds[i]] # this allocates a matrix
    #     S_i = TD.supersaturation.(thermo_params, qs, @view(parent(ρ)[valid_inds]), @view(parent(T)[valid_inds]), TD.Ice()) # this allocation seems hard to get around while maintaining only matrix ops
    #     N_INP = get_INP_concentration.(param_set, relaxation_timescale, qs, @view(parent(T)[valid_inds]), @view(parent(ρ)[valid_inds]), @view(parent(w)[valid_inds]), S_i) # this allocation seems hard to get around while maintaining only matrix ops
    #     dNINP_dz = get_dNINP_dz.(param_set, relaxation_timescale, @view(parent(T)[valid_inds]), @view(parent(dTdz)[valid_inds]))
    #     @view(parent(N_ice)[valid_inds]) .= adjust_ice_N_no_kwargs.(param_set,
    #         @view(parent(N_ice)[valid_inds]),
    #         N_INP .* @view(parent(f_ice_mult)[valid_inds]),
    #         @view(parent(q)[valid_inds, 3]),
    #         @view(parent(ρ)[valid_inds]),
    #         S_i,
    #         @view(parent(q)[valid_inds, 2]), # q_l for boosting <r>
    #         @view(parent(q_sno)[valid_inds]),
    #         @view(parent(massflux)[valid_inds]),
    #         dNINP_dz,
    #         @view(parent(w_i)[valid_inds]),
    #         N_INP_top .* @view(parent(f_ice_mult)[valid_inds]),
    #         ; monodisperse=false, decrease_N_if_subsaturated=true, N_i_from_INP = false, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)

        q_p      = parent(q)
        ρ_p      = parent(ρ)
        T_p      = parent(T)
        w_p      = parent(w)
        dTdz_p   = parent(dTdz)
        N_ice_p  = parent(N_ice)
        f_mult_p = parent(f_ice_mult)
        q_sno_p  = parent(q_sno)
        mf_p     = parent(massflux)
        w_i_p    = parent(w_i)
        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue

            qs = TD.PhasePartition(q_p[i,1], q_p[i,2], q_p[i,3])
            S_i = TD.supersaturation(thermo_params, qs, ρ_p[i], T_p[i], TD.Ice())
            N_INP = get_INP_concentration(param_set, relaxation_timescale, qs, T_p[i], ρ_p[i], w_p[i], S_i)
            dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T_p[i], dTdz_p[i])
            N_ice_p[i] = adjust_ice_N_no_kwargs(
                param_set,
                N_ice_p[i],
                N_INP * f_mult_p[i],
                q_p[i, 3],
                ρ_p[i],
                S_i,
                q_p[i, 2],
                q_sno_p[i],
                mf_p[i],
                dNINP_dz,
                w_i_p[i],
                N_INP_top * f_mult_p[i];
                monodisperse=false,
                decrease_N_if_subsaturated=true,
                N_i_from_INP=false,
                apply_massflux_boost=apply_massflux_boost,
                apply_sedimentation_boost=apply_sedimentation_boost,
            )
        end
    end

    # Because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here [[ (clamping the static array didn't work) -- e.g. clamp(_N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice) return an error bc of the stride arrays ]]
    clamp!(@view(parent(N_liq)[valid_inds]), relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    clamp!(@view(parent(N_ice)[valid_inds]), relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)

    return 
end



function get_Ns_and_N_i_no_boost!(param_set::APS, microphys_params::ACMP, relaxation_timescale::NeuralNetworkRelaxationTimescale, q::CC.Fields.Field, T::CC.Fields.Field, p::CC.Fields.Field, ρ::CC.Fields.Field, w::CC.Fields.Field, area::CC.Fields.Field, τ_liq::CC.Fields.Field, τ_ice::CC.Fields.Field, N_liq::CC.Fields.Field, N_ice::CC.Fields.Field, N_ice_no_boost::CC.Fields.Field, f_ice_mult::CC.Fields.Field, q_sno::CC.Fields.Field, massflux::CC.Fields.Field, dTdz::CC.Fields.Field, w_i::CC.Fields.Field; N_INP_top = eltype(param_set)(NaN), apply_massflux_boost::Bool=false, apply_sedimentation_boost::Bool=false)
    # model_x_0_characteristic = relaxation_timescale.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale

    # NOTE: parent(q) conveniently becomes a n x 3 Matrix (3 for [tot, liq, ice]), so q_liq = parent(q)[..., 2] and q_ice = parent(q)[..., 3]

    valid_inds = vec(parent(area) .> 0) # get the valid indices where area is greater than 0, convert to vec from nx1 matrix
    _, _, _N_l, _N_i = predict_τ(parent(ρ)[valid_inds], parent(T)[valid_inds], parent(q)[valid_inds, 2], parent(q)[valid_inds, 3], parent(w)[valid_inds], (relaxation_timescale.neural_network, to_static_strided_array(relaxation_timescale.neural_network_params)), relaxation_timescale.model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?

    parent(N_liq)[valid_inds] .= _N_l
    parent(N_ice)[valid_inds] .= _N_i

    if get_adjust_liq_N(relaxation_timescale) || get_adjust_ice_N(relaxation_timescale)
        q_p      = parent(q)
        ρ_p      = parent(ρ)
        T_p      = parent(T)
        w_p      = parent(w)
        dTdz_p   = parent(dTdz)
        N_liq_p  = parent(N_liq)
        N_ice_p  = parent(N_ice)
        f_mult_p = parent(f_ice_mult)
        q_sno_p  = parent(q_sno)
        mf_p     = parent(massflux)
        w_i_p    = parent(w_i)
    end

    # Because  this bypasses get_N_l(), get_N_i(), we need to apply adjust_liq/ice_N() here.
    if get_adjust_liq_N(relaxation_timescale)
        # @view(parent(N_liq)[valid_inds]) .= adjust_liq_N_no_kwargs.(param_set, @view(parent(N_liq)[valid_inds]), @view(parent(q)[valid_inds, 2]), @view(parent(ρ)[valid_inds]); monodisperse=true, decrease_N_if_subsaturated=false)
        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue
            N_liq_p[i] = adjust_liq_N_no_kwargs(param_set, N_liq_p[i], q_p[i, 2], ρ_p[i]; monodisperse=true, decrease_N_if_subsaturated=false)
        end
    end
    if get_adjust_ice_N(relaxation_timescale)
        thermo_params = TCP.thermodynamics_params(param_set)
        # FT = eltype(param_set)
        # qs = TD.PhasePartition{FT}[TD.PhasePartition(parent(q)[i, 1], parent(q)[i, 2], parent(q)[i, 3]) for i in eachindex(valid_inds) if valid_inds[i]]
        # S_i = TD.supersaturation.(thermo_params, qs, @view(parent(ρ)[valid_inds]), @view(parent(T)[valid_inds]), TD.Ice())
        # N_INP = get_INP_concentration.(param_set, relaxation_timescale, qs, @view(parent(T)[valid_inds]), @view(parent(ρ)[valid_inds]), @view(parent(w)[valid_inds]))
        # dNINP_dz = get_dNINP_dz.(param_set, relaxation_timescale, @view(parent(T)[valid_inds]), @view(parent(dTdz)[valid_inds]))
        # @view(parent(N_ice)[valid_inds]) .= adjust_ice_N_no_kwargs.(param_set,
        #     @view(parent(N_ice)[valid_inds]),
        #     N_INP .* @view(parent(f_ice_mult)[valid_inds]),
        #     @view(parent(q)[valid_inds, 3]),
        #     @view(parent(ρ)[valid_inds]),
        #     S_i,
        #     @view(parent(q)[valid_inds, 2]), # q_l for boosting <r>
        #     @view(parent(q_sno)[valid_inds]),
        #     @view(parent(massflux)[valid_inds]),
        #     dNINP_dz,
        #     @view(parent(w_i)[valid_inds]),
        #     N_INP_top .* @view(parent(f_ice_mult)[valid_inds]),
        #     ; monodisperse=false, decrease_N_if_subsaturated=true, N_i_from_INP = false, apply_massflux_boost = apply_massflux_boost, apply_sedimentation_boost = apply_sedimentation_boost)
    
    
    
        # if apply_massflux_boost
        #     @view(parent(N_ice_no_boost)[valid_inds]) .= adjust_ice_N_no_kwargs.(param_set,
        #         @view(parent(N_ice)[valid_inds]),
        #         N_INP .* @view(parent(f_ice_mult)[valid_inds]),
        #         @view(parent(q)[valid_inds, 3]),
        #         @view(parent(ρ)[valid_inds]),
        #         S_i,
        #         @view(parent(q)[valid_inds, 2]), # q_l for boosting <r>
        #         @view(parent(q_sno)[valid_inds]),
        #         @view(parent(massflux)[valid_inds]),
        #         dNINP_dz,
        #         @view(parent(w_i)[valid_inds]),
        #         N_INP_top .* @view(parent(f_ice_mult)[valid_inds]),
        #         ; monodisperse=false, decrease_N_if_subsaturated=true, N_i_from_INP = false, apply_massflux_boost = false, apply_sedimentation_boost = apply_sedimentation_boost)
        # else
        #     @view(parent(N_ice_no_boost)[valid_inds]) .= @view(parent(N_ice)[valid_inds])
        # end

        @inbounds for i in eachindex(valid_inds)
            valid_inds[i] || continue

            qs = TD.PhasePartition(q_p[i,1], q_p[i,2], q_p[i,3])
            S_i = TD.supersaturation(thermo_params, qs, ρ_p[i], T_p[i], TD.Ice())
            N_INP = get_INP_concentration(param_set, relaxation_timescale, qs, T_p[i], ρ_p[i], w_p[i], S_i)
            dNINP_dz = get_dNINP_dz(param_set, relaxation_timescale, T_p[i], dTdz_p[i])
            N_ice_p[i] = adjust_ice_N_no_kwargs(
                param_set,
                N_ice_p[i],
                N_INP * f_mult_p[i],
                q_p[i, 3],
                ρ_p[i],
                S_i,
                q_p[i, 2],
                q_sno_p[i],
                mf_p[i],
                dNINP_dz,
                w_i_p[i],
                N_INP_top * f_mult_p[i];
                monodisperse=false,
                decrease_N_if_subsaturated=true,
                N_i_from_INP=false,
                apply_massflux_boost=apply_massflux_boost,
                apply_sedimentation_boost=apply_sedimentation_boost,
            )

            if apply_massflux_boost
                parent(N_ice_no_boost)[i] = adjust_ice_N_no_kwargs(
                    param_set,
                    N_ice_p[i],
                    N_INP * f_mult_p[i],
                    q_p[i, 3],
                    ρ_p[i],
                    S_i,
                    q_p[i, 2],
                    q_sno_p[i],
                    mf_p[i],
                    dNINP_dz,
                    w_i_p[i],
                    N_INP_top * f_mult_p[i];
                    monodisperse=false,
                    decrease_N_if_subsaturated=true,
                    N_i_from_INP=false,
                    apply_massflux_boost=false,
                    apply_sedimentation_boost=apply_sedimentation_boost,
                )
            else
                parent(N_ice_no_boost)[i] = N_ice_p[i]
            end
        end
    
    end


    # Because this bypasses get_N_l(), get_N_i(), and get_τs() we need to clamp the values here [[ (clamping the static array didn't work) -- e.g. clamp(_N_i, relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice) return an error bc of the stride arrays ]]
    clamp!(@view(parent(N_liq)[valid_inds]), relaxation_timescale.args.min_N_liq, relaxation_timescale.args.max_N_liq)
    clamp!(@view(parent(N_ice)[valid_inds]), relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)
    clamp!(@view(parent(N_ice_no_boost)[valid_inds]), relaxation_timescale.args.min_N_ice, relaxation_timescale.args.max_N_ice)

    return 
end


# ======================================================================================================================================== #
# ---------------------------------------------------------------------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------------------------------------------------------------------- #
