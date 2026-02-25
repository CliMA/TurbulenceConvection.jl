"""
    condensate_qt_SD(qt::FT, qt_var::FT, q_sat::FT, qc_real::FT, σ_max::FT = FT(1)) where {FT}

Estimates the condensate-weighted mean total water (`qt_cond`) and calculates 
how many standard deviations that conditional mean is shifted from the grid mean (`condensate_qt_SD`).

The shift is derived from the exact analytical integration of the condensate-weighted 
1st moment of a Gaussian distribution.
"""
function condensate_qt_SD(
    qt::FT,
    qt_var::FT,
    q_sat::FT,
    # qc_real::FT, # Kept for API compatibility, but unused in the pure statistical shift
    σ_max::FT = FT(1),
) where {FT}

    if qt_var <= eps(FT)
        return qt, zero(FT)
    end

    qt_std = sqrt(qt_var)
    z_star = (q_sat - qt) / qt_std

    # Asymptotic limit for highly subsaturated clouds (evaporating/non-equilibrium).
    # As z* -> ∞, the analytical shift cleanly converges to exactly z*.
    # We catch this early to avoid 0/0 floating-point underflow (NaNs) in the math below.
    if z_star > FT(3)
        condensate_qt_SD = min(z_star, σ_max)
        qt_cond = qt + condensate_qt_SD * qt_std
        return qt_cond, condensate_qt_SD
    end

    dist = Distributions.Normal(zero(FT), one(FT))
    
    # Use ccdf for better numerical precision in the tail
    prob_sat = Distributions.ccdf(dist, z_star) 
    phi_z = Distributions.pdf(dist, z_star)        

    # Theoretical equilibrium condensate
    qc_sat_adjust = qt_std * max(eps(FT), phi_z - z_star * prob_sat)

    # 1. Pure statistical condensate-weighted mean shift
    shift = qt_var * prob_sat / qc_sat_adjust
    qt_cond = qt + shift

    # 2. Normalized shift (capped by σ_max to prevent "cooked SD" in out-of-equilibrium tails)
    condensate_qt_SD = min(shift / qt_std, σ_max)

    return qt_cond, condensate_qt_SD
end


function partition_condensate_into_sgs_fractions(
    q::FT,
    supersaturated_fraction::FT;
    max_boost_factor::FT = FT(1.0),
) where {FT}
    # partition assuming all liq into supersaturated fraction, not to exceed max_boost_factor
    q_supersat = iszero(supersaturated_fraction) ? (max_boost_factor * q) : min(q / supersaturated_fraction, max_boost_factor * q)
    q_subsat = isone(supersaturated_fraction) ? zero(FT) : (q - (q_supersat * supersaturated_fraction)) / (1 - supersaturated_fraction)

    return (; supersaturated_fraction, q_supersat, q_subsat)
end

"""
    Partitions cloud liquid into supersaturated and subsaturated regions
    This is a quick fcn for when you dont want to make up an assumed PDF esp given no saturation adjustment.
"""
function partition_condensate_into_sgs_fractions(
    thermo_params::TDPS,
    ::CMT.LiquidType,
    ts::TD.ThermodynamicState{FT},
    # θ_liq_ice::FT,
    qt_var::FT,
    h_var::FT,
    h_qt_cov::FT,
    ;
    supersaturated_fraction_liq::FT = FT(NaN),
    max_boost_factor::FT = FT(2.0),
) where {FT}
    q = TD.PhasePartition(ts)
    if isnan(supersaturated_fraction_liq)
        T = TD.air_temperature(thermo_params, ts)
        ρ = TD.air_density(thermo_params, ts)
        p = TD.air_pressure(thermo_params, ts)
        q_vap = TD.vapor_specific_humidity(thermo_params, ts)
        q_sat_l = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
        dq_sat_dT_l = TD.∂q_vap_sat_∂T(thermo_params, one(FT), T, q_sat_l)  # Calculate derivative
        supersaturated_fraction_liq = sgs_saturated_fraction(thermo_params, q_vap, q_sat_l, qt_var, h_var, h_qt_cov, dq_sat_dT_l, p)
    end
    return partition_condensate_into_sgs_fractions(TD.Liquid(), q.liq, supersaturated_fraction_liq; max_boost_factor)
end
partition_condensate_into_sgs_fractions(::CMT.LiquidType, ql::FT, supersaturated_fraction_liq::FT; max_boost_factor::FT = FT(2.0)) where {FT} = partition_condensate_into_sgs_fractions(ql, supersaturated_fraction_liq; max_boost_factor)



"""
Like for liquid but τ_i is much slower than τ_l so our max boost factor default is 1.1 (could be calibrated)
"""
function partition_condensate_into_sgs_fractions(
    thermo_params::TDPS,
    ::CMT.IceType,
    ts::TD.ThermodynamicState{FT},
    # θ_liq_ice::FT,
    qt_var::FT,
    h_var::FT,
    h_qt_cov::FT,
    ;
    supersaturated_fraction_ice::FT = FT(NaN),
    max_boost_factor::FT = FT(1.1),
) where {FT}
    q = TD.PhasePartition(ts)
    if isnan(supersaturated_fraction_ice)
        T = TD.air_temperature(thermo_params, ts)
        ρ = TD.air_density(thermo_params, ts)
        p = TD.air_pressure(thermo_params, ts)
        q_vap = TD.vapor_specific_humidity(thermo_params, ts)
        q_sat = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
        dq_sat_dT = TD.∂q_vap_sat_∂T(thermo_params, zero(FT), T, q_sat)  # Calculate derivative
        supersaturated_fraction_ice = sgs_saturated_fraction(thermo_params, q_vap, q_sat, qt_var, h_var, h_qt_cov, dq_sat_dT, p) # pass q_vap here to get supersat fraction (pass q_tot to get cloud fraction (under saturation adjustment to form ql))
    end
    return partition_condensate_into_sgs_fractions(TD.Ice(), q.ice, supersaturated_fraction_ice; max_boost_factor)
end
partition_condensate_into_sgs_fractions(::CMT.IceType, qi::FT, supersaturated_fraction_ice::FT; max_boost_factor::FT = FT(1.1)) where {FT} = partition_condensate_into_sgs_fractions(qi, supersaturated_fraction_ice; max_boost_factor)

# ======================================================================================================================================================== #
# ======================================================================================================================================================== #

"""
    get_qc_stats_from_moments(...)

Calculates turbulent fluxes using a 'Physical Blending' approach for Ice:

1. LIQUID FLUX (Path A):
   Derived strictly from the Liquid PDF (Thermodynamic Phase Lock).
   w'ql' = min( w'qc_pot, (w'qc_pot / qc_pdf_pot) * ql_true )

2. ICE FLUX (Path B):
   Calculated as: Flux = Transport_Velocity * Ice_Variance
   
   - Transport Velocity (V_trans): Derived from Liquid PDF (w'qc_pot / sigma_pot).
   - Ice Variance: Blends two regimes based on Liquid Fraction (f_liq):
       * Active (f_liq ~ 1): Scales with Liquid Potential CV (Phase Locked).
       * Passive (f_liq ~ 0): Scales with Total Water CV (Tracer).
   
   This captures the transition from 'Active' mixed-phase turbulence to 'Passive' cirrus mixing.
"""
function get_qc_stats_from_moments(
    param_set::APS,
    ts::TD.ThermodynamicState,
    qt_var::FT,
    h_var::FT,
    h_qt_cov::FT,
    w_qt_cov::FT, 
    w_h_cov::FT, 
    ql_true::FT = FT(NaN),
    qi_true::FT = FT(NaN)
) where {FT}

    # --- 1. SETUP ---
    thermo_params = TCP.thermodynamics_params(param_set)
    has_prognostic = !isnan(ql_true) && !isnan(qi_true)
    
    p = TD.air_pressure(thermo_params, ts)
    θli = TD.liquid_ice_pottemp(thermo_params, ts)
    qt = TD.total_specific_humidity(thermo_params, ts)

    # LIQUID BASELINE (PDF Structure)
    ts_pdf = PhaseEquil_pθq_given_liquid_fraction(thermo_params, p, θli, qt, FT(1.0))
    T_pdf = TD.air_temperature(thermo_params, ts_pdf)
    ρ_pdf = TD.air_density(thermo_params, ts_pdf)
    Π_pdf = TD.exner(thermo_params, ts_pdf)
    q_sat = TD.q_vap_saturation_generic(thermo_params, T_pdf, ρ_pdf, TD.Liquid())
    dqsat_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(1.0), T_pdf, q_sat)

    # --- 2. PDF MOMENTS ---
    s_mean = qt - q_sat 
    dsdqt = one(FT)
    dsdθli = -dqsat_dT * Π_pdf
    
    # Variance of saturation deficit (s)
    # Floor to prevent division by zero in standard deviation
    s_var = max(dsdqt^2 * qt_var + dsdθli^2 * h_var + FT(2) * dsdqt * dsdθli * h_qt_cov, FT(1e-12))
    s_std = sqrt(s_var)
    
    α = s_mean / s_std
    inv_sqrt2 = one(FT) / sqrt(FT(2))
    Φ = FT(0.5) * (one(FT) + SpecialFunctions.erf(α * inv_sqrt2))
    
    # PDF Potentials (Mean and Flux)
    # qc_pdf_pot: The mean condensate predicted by the Liquid PDF
    qc_pdf_pot = s_mean * Φ + s_std * (exp(-FT(0.5) * α^2) / sqrt(FT(2) * FT(π)))
    w_qc_pot = (dsdqt * w_qt_cov + dsdθli * w_h_cov) * Φ 
    
    # Scalar Potentials
    cov_qc_qt_pot = (dsdqt * qt_var + dsdθli * h_qt_cov) * Φ
    cov_qc_h_pot  = (dsdqt * h_qt_cov + dsdθli * h_var) * Φ
    
    # Proxy for the Standard Deviation of the "Cloudy" part of the PDF
    # Used to normalize the Transport Velocity.
    sigma_pot = s_std * sqrt(Φ)

    # --- 3. FLUX CALCULATION ---
    ql_out, qi_out = (has_prognostic) ? (ql_true, qi_true) : (qc_pdf_pot, zero(FT))
    w_ql_cov = zero(FT)
    w_qi_cov = zero(FT)
    cov_ql_qt = zero(FT)
    cov_ql_h = zero(FT)

    if has_prognostic
        # =====================================================================
        # PATH A: LIQUID (Strictly Thermodynamic)
        # =====================================================================
        # Liquid flux follows the PDF potential, strictly clamped to the ceiling.
        # This prevents unphysical extrapolation (w'ql > w'qc_pot).
        if qc_pdf_pot > FT(1e-10)
            eff_w = w_qc_pot / qc_pdf_pot
            w_ql_cov = min(eff_w * ql_true, w_qc_pot)
            
            eff_qt = cov_qc_qt_pot / qc_pdf_pot
            cov_ql_qt = min(eff_qt * ql_true, cov_qc_qt_pot)
            
            eff_h = cov_qc_h_pot / qc_pdf_pot
            cov_ql_h = min(eff_h * ql_true, cov_qc_h_pot)
        end

        # =====================================================================
        # PATH B: ICE (Physical Blending)
        # =====================================================================
        
        # 1. Transport Velocity (V_trans)
        # "How fast is the turbulent cloudy air moving?"
        # V_trans = Flux_Potential / Sigma_Potential
        transport_vel = zero(FT)
        if sigma_pot > FT(1e-10)
            transport_vel = w_qc_pot / sigma_pot
        end
        
        # 2. Determine Regimes (Liquid Fraction)
        qc_true_total = ql_true + qi_true
        f_liq_prog = (qc_true_total > FT(1e-10)) ? ql_true / qc_true_total : zero(FT)
        
        # 3. Estimate Ice Variance (Sigma_qi)
        
        #   Regime 1: Passive (Total Water Scaling)
        #   Sigma_Passive = qi * (sigma_qt / qt)
        sigma_ice_passive = zero(FT)
        if qt > FT(1e-6)
            sigma_qt = sqrt(qt_var)
            sigma_ice_passive = qi_true * (sigma_qt / qt)
        end
        
        #   Regime 2: Active (Liquid Potential Scaling)
        #   Sigma_Active = qi * (sigma_pot / qc_pot)
        #   We use the Coefficient of Variation from the PDF.
        #   Supersat contributes to ice supersat but only tangentially... so high qi in dry areas has low variance...
        sigma_ice_active = zero(FT)
        if qc_pdf_pot > FT(1e-10)
            sigma_ice_active = qi_true * (sigma_pot / qc_pdf_pot)
        end
        sigma_ice_active *= Φ # scale by cloud fraction to prevent overestimation in subsaturated conditions

        #   Blend the Standard Deviations
        sigma_ice_final = ((one(FT) - f_liq_prog) * sigma_ice_passive) + (f_liq_prog * sigma_ice_active) # the first term is tiny, almost everything is in the second term...
        
        sigma_ice_final *= (qc_pdf_pot / qt)


        # 4. Final Ice Flux
        # w'qi' = V_trans * Sigma_qi
        w_qi_cov = transport_vel * sigma_ice_final
        # w_qi_cov = sigma_ice_final

    else
        # Diagnostic Case
        w_ql_cov = w_qc_pot
        w_qi_cov = zero(FT)
        cov_ql_qt = cov_qc_qt_pot
        cov_ql_h  = cov_qc_h_pot
    end

    # --- 4. OUTPUT PREP ---
    qt_std = sqrt(qt_var)
    qc_total = ql_out + qi_out
    final_z = (qt_std > eps(FT) && qc_total > FT(1e-8)) ? 
        clamp(cov_ql_qt / (qc_total * qt_std), zero(FT), FT(2.5)) : zero(FT)

    f_liq = (qc_total > eps(FT)) ? ql_out / qc_total : one(FT)
    f_ice = (qc_total > eps(FT)) ? qi_out / qc_total : zero(FT)
    
    ql_pdf = qc_pdf_pot * f_liq
    qi_pdf = qc_pdf_pot * f_ice
    w_ql_pot = w_qc_pot * f_liq
    w_qi_pot = w_qc_pot * f_ice
    
    scale_l = (w_ql_pot > eps(FT)) ? min(w_ql_cov / w_ql_pot, FT(1.0)) : one(FT)
    scale_i = (w_qi_pot > eps(FT)) ? min(w_qi_cov / w_qi_pot, FT(1.0)) : one(FT)

    return ql_out, qi_out, cov_ql_qt, cov_ql_h, w_ql_cov, w_qi_cov, final_z, final_z, qc_pdf_pot, w_qc_pot, ql_pdf, qi_pdf, Φ, α, f_liq, f_ice, scale_l, scale_i, w_ql_pot, w_qi_pot
end

# """
#     get_saturation_adjust_cond_evap(param_set, ts, qt_var, h_var, h_qt_cov, Δt, ql, qi)

# Returns the gross SGS condensation and evaporation rates for latent heating and TKE production.

# Computes the target condensate from the SGS supersaturation integral and compares to the current
# condensate (ql + qi) to determine whether condensation or evaporation occurred:
# - If qc_old > qc_target: evaporation removed the excess
# - If qc_old < qc_target: condensation formed the difference

# Both rates capture the SGS processes driving latent heat release, even when net changes are small.
# For nonequilibrium, preserves the existing liquid fraction when computing the saturation state.
# """
# function get_saturation_adjust_cond_evap(
#     param_set::APS,
#     ts::TD.ThermodynamicState{FT},
#     qt_var::FT,
#     h_var::FT,
#     h_qt_cov::FT,
#     Δt::FT,
#     ql::FT,
#     qi::FT,
# ) where {FT}

#     error("We can't raelly do this without knowing the past qt, and rates of thigns going to precip etc...")

#     thermo_params = TCP.thermodynamics_params(param_set)

#     if Δt <= zero(FT)
#         return zero(FT), zero(FT)
#     end

#     p = TD.air_pressure(thermo_params, ts)
#     θli = TD.liquid_ice_pottemp(thermo_params, ts)
#     qt = TD.total_specific_humidity(thermo_params, ts)

#     # Preserve existing liquid fraction (for noneq) when defining the saturation-adjusted state.
#     denom = ql + qi
#     liq_frac = (denom > eps(FT)) ? clamp(ql / denom, zero(FT), one(FT)) : TD.liquid_fraction(thermo_params, ts)

#     # Linearize saturation deficit around the saturation-adjusted mean state.
#     ts_sat = PhaseEquil_pθq_given_liquid_fraction(thermo_params, p, θli, qt, liq_frac)
#     T_sat = TD.air_temperature(thermo_params, ts_sat)
#     q_sat = TD.vapor_specific_humidity(thermo_params, ts_sat)
#     Π_sat = TD.exner(thermo_params, ts_sat)
#     dqsat_dT = TD.∂q_vap_sat_∂T(thermo_params, liq_frac, T_sat, q_sat)

#     s_mean = qt - q_sat
#     dsdqt = one(FT)
#     dsdθli = -dqsat_dT * Π_sat

#     qt_var = resolve_nan(qt_var, zero(FT))
#     h_var = resolve_nan(h_var, zero(FT))
#     h_qt_cov = resolve_nan(h_qt_cov, zero(FT))

#     # SGS saturation-deficit variance from linearized thermodynamics.
#     s_var = max(dsdqt^2 * qt_var + dsdθli^2 * h_var + FT(2) * dsdqt * dsdθli * h_qt_cov, FT(1e-12))
#     s_std = sqrt(s_var)

#     # Target condensate from SGS supersaturation integral: qc_target = E[max(s, 0)]
#     qc_target = if s_std <= eps(FT)
#         max(s_mean, zero(FT))
#     else
#         α = s_mean / s_std
#         inv_sqrt2 = one(FT) / sqrt(FT(2))
#         Φ = FT(0.5) * (one(FT) + SpecialFunctions.erf(α * inv_sqrt2))
#         φ = exp(-FT(0.5) * α^2) / sqrt(FT(2) * FT(π))
#         max(s_mean * Φ + s_std * φ, zero(FT))
#     end

#     # Current condensate before adjustment.
#     qc_old = max(ql + qi, zero(FT))

#     # Determine gross rates from the change.
#     Δqc = qc_target - qc_old
#     cond_rate = max(Δqc, zero(FT)) / Δt
#     evap_rate = max(-Δqc, zero(FT)) / Δt

#     return cond_rate, evap_rate
# end