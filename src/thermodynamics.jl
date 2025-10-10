function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT) where {FT}
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    thermo_params = TCP.thermodynamics_params(param_set)
    return TD.PhaseEquil_pθq(thermo_params, p, θ_liq_ice, q_tot, config...)
end
function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    config = ()
    # q_tot = max(q_tot, q_liq + q_ice + eps(FT)) # ensure that total specific humidity is at least the sum of liquid and ice
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    thermo_params = TCP.thermodynamics_params(param_set)
    return TD.PhaseNonEquil_pθq(thermo_params, p, θ_liq_ice, q, config...)
end

function geopotential(param_set, z::Real)
    FT = eltype(param_set)
    grav = FT(TCP.grav(param_set))
    return grav * z
end

function enthalpy(mse::FT, e_pot::FT) where {FT}
    return mse - e_pot
end




function relative_humidity_over_ice(thermo_params::TD.APS, T::FT, p::FT, ::Type{phase_type},
    q::TD.PhasePartition{FT} = TD.q_pt_0(FT),
) where {FT <: Real, phase_type <: TD.ThermodynamicState}
    R_v::FT = TD.TP.R_v(thermo_params)
    q_vap = TD.vapor_specific_humidity(q)
    p_vap = q_vap * TD.air_density(thermo_params, T, p, q) * R_v * T
    p_vap_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    return p_vap / p_vap_sat
end
relative_humidity_over_ice(param_set::APS, T::FT, p::FT, ::Type{phase_type}) where {FT <: Real, phase_type <: TD.ThermodynamicState} = 
    relative_humidity_over_ice(TCP.thermodynamics_params(param_set), T, p, phase_type)

"""
The default assumes over liquid, or if you provide a q, it could be over ice but it doesnt pass that q to saturation_vapor_pressure anyway...
"""
relative_humidity_over_ice(thermo_params::TD.APS, ts::TD.ThermodynamicState{FT}) where {FT <: Real} = 
    relative_humidity_over_ice(
        thermo_params,
        TD.air_temperature(thermo_params, ts),
        TD.air_pressure(thermo_params, ts),
        typeof(ts),
        TD.PhasePartition(thermo_params, ts),
    )


relative_humidity_over_ice(param_set::APS, ts::TD.ThermodynamicState{FT}) where {FT <: Real} = relative_humidity_over_ice(
    TCP.thermodynamics_params(param_set),
    ts,
)
