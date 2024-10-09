"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function morrison_milbrandt_2015_style(
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::Real,
    ts::TD.ThermodynamicState,
) where {FT}
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

        """
        # Milbrandt's equations are in mixing ratio not specific humidity, to convert we need to use

        # q_spe = q_mix / (1 + q_mix)   [ in reality this is more like q_spe_c = q_mix_c / (1 + q_tot_mix) , bc 1+q_tot = total ratio all water / dry air...]
        # q_mix = q_spe / (1 - q_spe)   [ in reality this is more like q_mix_c = q_spe_c / (1 - q_tot_spe) , bc 1-q_tot = dry mixing ratio...]

        # [NOTE: q_tot appears instead of just q because we have more than one water species so we need to be careful...]

        # note then d/dy (q_spe) = d/dy(q_mix / (1 + q_mix)) = d/dy(q_mix) / (1 + q_mix)^2  [ really it's more like d/dy(q_spe_c) = d/dy(q_mix_c / (1 + q_tot_mix)) and q_tot_mix should be a constant under phase transformation... ]
        # note then d/dy (q_mix) = d/dy(q_spe / (1 - q_spe)) = d/dy(q_spe) / (1 - q_spe)^2  [ really it's more like d/dy(q_mix_c) = d/dy(q_spe_c / (1 - q_tot_spe)) and q_tot_spe should be a constant under phase transformation... ]

        NOTE: These mixing ratio, specific humidity, and their derivatives are essentially identical at earthlike conditions...
        """

        # T = TD.air_temperature(thermo_params, ts) # air temperature
        # ρ = TD.air_density(thermo_params, ts) # air density
        # p = TD.air_pressure(thermo_params, ts) # air pressure

        g = TCP.grav(param_set) # acceleration of gravity
        L_i = TD.latent_heat_sublim(thermo_params, ts) # Latent heat for ice (L_s)
        L_l = TD.latent_heat_vapor(thermo_params, ts)  # Latent heat for water
        # c_p = TD.cp_d(param_set) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?
        c_p = TD.cp_m(thermo_params, q) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?

        e_sl = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) # saturation vapor pressure, or TD.partial_pressure_vapor(thermo_params, p, q_eq) 

        # analytical forms for change in saturation vapor presure with temperature
        dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, ts) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix... is it a good enough approx?
        dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, ts) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...

        dqsl_dT /= (1 - q_eq.liq)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.liq = q_sl at T
        dqsi_dT /= (1 - q_eq.ice)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.ice = q_si at T

        # So, we can do everything in mixing ratio and then convert back at the end instead of trying to convert everythig in their derivation to figure out what q_tot_mix
        q_vap = TD.shum_to_mixing_ratio(q_vap, q.tot)
        q_sl = TD.shum_to_mixing_ratio(q_eq.liq, q_eq.tot) # q_eq is the equilibrium mixing ratio, q_sl is the saturation mixing ratio, q_eq we made to contain only saturation vapor pressure values...
        q_si = TD.shum_to_mixing_ratio(q_eq.ice, q_eq.tot)

        # saturation_mixing_ratio_liq = TD.shum_to_mixing_ratio( TD.q_vap_saturation_generic(thermo_params, T, p, TD.Liquid()), q.tot)
        # saturation_mixing_ratio_ice = TD.shum_to_mixing_ratio( TD.q_vap_saturation_generic(thermo_params, T, p, TD.Ice()), q.tot)


        δ_0::FT = q_vap - q_sl # supersaturation over liquid

        Γ_l = 1 + L_l / c_p * dqsl_dT # Eqn C3
        Γ_i = 1 + L_i / c_p * dqsi_dT # Eqn C3

        dTdt_mix = FT(0) # ignore for now
        dTdt_rad = FT(0) # ignore for now

        τ::FT = 1 / (1 / τ_liq + (1 + (L_i / c_p) * dqsl_dT) * ((1 / τ_ice) / (Γ_i))) # Eqn C2

        A_c::FT =
            dTdt_mix - (q_sl * ρ * g * w) / (p - e_sl) - dqsl_dT * (dTdt_rad + dTdt_mix - (w * g) / c_p) -
            (q_sl - q_si) / (τ_ice * Γ_l) * (1 + (L_i / c_p) * dqsl_dT) # Eq C4

        S_ql = A_c * τ / (τ_liq * Γ_l) + (δ_0 - A_c * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ))   # QCCON EQN C6
        S_qi =
            A_c * τ / (τ_ice * Γ_i) +
            (δ_0 - A_c * τ) * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ)) +
            (q_sl - q_si) / (τ_ice * Γ_i)  # QICON Eqn C7


        # Go from mixing ratio world back to specific humidity world
        q_tot_mix = q.tot / (1 - q.tot)
        S_ql /= (1 + q_tot_mix) # bc we showed  d/dt(q_spe) = d/dt(q_mix) / (1+q_tot_mix) under phase transformation, which makes sense, it's smaller
        S_qi /= (1 + q_tot_mix)

    end

    return S_ql, S_qi
end



function morrison_milbrandt_2015_style_exponential_part_only(
    param_set::APS,
    area::FT,
    ρ::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    q_eq::TD.PhasePartition,
    Δt::Real,
) where {FT}
    """
    This will be a fcn that handles only the exponential decay part for supersaturation w/ liquid and ice, but ignores vertical velocity and everything else (temperature changes and stuff)
    """
    # w = FT(0) # ignore for now
    # return morrison_milbrandt_2015_style(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap, q_eq, Δt)

    S_ql::FT = FT(0)
    S_qi::FT = FT(0)

    # can we cache repeated terms? compiler probably figures out that anyway...
    if area > 0
        δ_0::FT = q_vap - q_eq.liq # supersaturation over liquid
        τ::FT = 1 / (1 / τ_liq + 1 / τ_ice)
        A_c::FT = -(q_eq.liq - q_eq.ice) / (τ_ice) # Eq C4
        S_ql = A_c * τ / (τ_liq) + (δ_0 - A_c * τ) * τ / (Δt * τ_liq) * (1 - exp(-Δt / τ))   # QCCON EQN C6
        S_qi =
            A_c * τ / (τ_ice) +
            (δ_0 - A_c * τ) * τ / (Δt * τ_ice) * (1 - exp(-Δt / τ)) +
            (q_eq.liq - q_eq.ice) / (τ_ice)  # QICON Eqn C7
    end

    return S_ql, S_qi

end


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


# @inline function saturation_vapor_pressure_derivative_Temperature(
#     param_set::APS,
#     T::FT,
#     LH_0::FT,
#     Δcp::FT,
# ) where {FT <: Real} # d p_s / dT
#     press_triple = TP.press_triple(param_set)
#     R_v = TP.R_v(param_set)
#     T_triple = TP.T_triple(param_set)
#     T_0 = TP.T_0(param_set)

#     return press_triple *
#            # (T / T_triple)^(Δcp / R_v) *
#            pow_hack(T / T_triple, Δcp / R_v) *
#            exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))

# end

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

# @inline function q_vap_saturation_derivative_Temperature(
#     param_set::APS,
#     ts::TD.ThermodynamicState
# ) where {FT <: Real}


#     return TD.∂q_vap_sat_∂T(param_set, ts)

# end
