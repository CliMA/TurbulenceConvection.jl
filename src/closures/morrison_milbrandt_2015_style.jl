"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function morrison_milbrandt_2015_style(param_set::APS, area::FT, ρ::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q_eq::TD.PhasePartition, Δt::Real) where {FT}
"""
See https://doi.org/10.1175/JAS-D-14-0065.1
we are ignoring the mixing and radiation terms for our short timesteps, as well as rain effects

this *shouldn't* need limiters the way it's defined but be careful I suppose lol...
"""
    thermo_params = TCP.thermodynamics_params(param_set) # currently repeated in both places, pare down later
    # microphys_params = TCP.microphysics_params(param_set)
    S_ql = FT(0)
    S_qi = FT(0)
    if area > 0

        g   = TCP.grav(param_set) # acceleration of gravity
        L_i =  TD.latent_heat_sublim(thermo_params, ts) # Latent heat for ice (L_s)
        L_l =  TD.latent_heat_vapor(thermo_params, ts)  # Latent heat for water

        e_sl = TD.partial_pressure_vapor(thermo_params,p,q_eq.liq) # water vapour pressure or TD.saturation_vapor_pressure(thermo_params, T, ρ, TD.Liquid())

        # analytical forms for change in saturation vapor presure with temperature
        dqsl_dT  = error("not implemented yet")
        dqsi_dT  = error("not implemented yet")

        δ_0 = q_vap - q_eq.liq # supersaturation over liquid

        Γ_l = 1 + L_l/c_p * dqsl_dT # Eqn C3
        Γ_i = 1 + L_i/c_p * dqsi_dT # Eqn C3

        dqdt_mix = FT(0) # ignore for now
        dTdt_rad = FT(0) # ignore for now

        τ = 1/(1/τ_liq + 1/τ_ice) + (1 + (L_s/c_p) * dqsldT) * ((1/τ_ice)/(Γ_i) ) # Eqn C2

        A_c = dqdt_mix - (q_eq.liq * ρ * g * w)/(p-e_sl) - dqsl_dT * (dTdt_rad + dTdt_mix - (w*g)/c_p) - (q_eq.liq - q_eq.ice)/(τ_ice * Γ_l) * (1 + (L_s/c_p) * dqsldT) # Eq C4

        S_ql = A_c * τ/(τ_liq * Γ_l) +  (δ_0 - A_c * τ)* τ/(Δt * τ_liq * Γ_l) * (1-exp(-Δt/τ))   # QCCON EQN C6
        S_qi = A_c * τ/(τ_ice * Γ_i) +  (δ_0 - A_c * τ)* τ/(Δt * τ_ice * Γ_i) * (1-exp(-Δt/τ)) + (q_eq.liq - q_eq.ice)/(τ_ice * Γ_i)  # QICON Eqn C7

    end

    return S_ql, S_qi
end


function morrison_milbrandt_2015_style_exponential_part_only(param_set::APS, area::FT, ρ::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q_eq::TD.PhasePartition, Δt::Real) where {FT}
    """
    This will be a fcn that handles only the exponential decay part for supersaturation w/ liquid and ice, but ignores vertical velocity and everything else (temperature changes and stuff)
    """
    # error("Not implemented yet")

    # w = FT(0) # ignore for now
    # return morrison_milbrandt_2015_style(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap, q_eq, Δt)

    S_ql::FT = FT(0)
    S_qi::FT = FT(0)

    # can we cache repeated terms? compiler probably figures out that anyway...
    if area > 0
        δ_0::FT = q_vap - q_eq.liq # supersaturation over liquid
        τ::FT = 1/(1/τ_liq + 1/τ_ice) 
        A_c::FT =  - (q_eq.liq - q_eq.ice)/(τ_ice )# Eq C4
        S_ql = A_c * τ/(τ_liq) +  (δ_0 - A_c * τ) * τ/(Δt * τ_liq ) * (1-exp(-Δt/τ))   # QCCON EQN C6
        S_qi = A_c * τ/(τ_ice) +  (δ_0 - A_c * τ) * τ/(Δt * τ_ice ) * (1-exp(-Δt/τ)) + (q_eq.liq - q_eq.ice)/(τ_ice)  # QICON Eqn C7
    end

    return S_ql, S_qi

end




# note we hand
# @inline pow_hack(x, y) = exp(y * log(x))

# @inline function saturation_vapor_pressure(
#     param_set::APS,
#     T::FT,
#     LH_0::FT,
#     Δcp::FT,
# ) where {FT <: Real}
#     press_triple = TP.press_triple(param_set)
#     R_v = TP.R_v(param_set)
#     T_triple = TP.T_triple(param_set)
#     T_0 = TP.T_0(param_set)

#     return press_triple *
#            # (T / T_triple)^(Δcp / R_v) *
#            pow_hack(T / T_triple, Δcp / R_v) *
#            exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))

# end


@inline function saturation_vapor_pressure_derivative_Temperature(
    param_set::APS,
    T::FT,
    LH_0::FT,
    Δcp::FT,
) where {FT <: Real}
    press_triple = TP.press_triple(param_set)
    R_v = TP.R_v(param_set)
    T_triple = TP.T_triple(param_set)
    T_0 = TP.T_0(param_set)

    return press_triple *
           # (T / T_triple)^(Δcp / R_v) *
           pow_hack(T / T_triple, Δcp / R_v) *
           exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))

end

# we had 
# @inline function q_vap_saturation_generic(
#     param_set::APS,
#     T::FT,
#     ρ::FT,
#     phase::Phase,
# ) where {FT <: Real}
#     p_v_sat = saturation_vapor_pressure(param_set, T, phase)
#     return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
# end
# q_vap_saturation_generic(param_set::APS, T, ρ, phase::Phase) =
#     q_vap_saturation_generic(param_set, promote(T, ρ)..., phase)

@inline function q_vap_saturation_derivative_Temperature()
end