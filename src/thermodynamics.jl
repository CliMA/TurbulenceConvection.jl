function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT) where {FT}
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    return TD.PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
end
function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    config = ()
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    return TD.PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q, config...)
end


function thermo_state_peq(param_set::APS, p::FT, e_int::FT, q_tot::FT) where {FT}
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    return TD.PhaseEquil_peq(param_set, p, e_int, q_tot, config...)
end

function thermo_state_peq(param_set::APS, p::FT, e_int::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    config = ()
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    return TD.PhaseNonEquil_peq(param_set, p, e_int, q, config...)
end

function geopotential(param_set, z::FT) where {FT}
    grav = FT(ICP.grav(param_set))
    return grav * z
end

function anelastic_total_enthalpy(param_set::APS, e_tot::FT, ts) where {FT}
    Rm = TD.gas_constant_air(param_set, ts)
    T = TD.air_temperature(param_set, ts)
    return e_tot + Rm * T
end

function anelastic_total_enthalpy(param_set::APS, e_tot::FT, p::FT, ρ::FT) where {FT}
    return e_tot + p / ρ
end

function kinetic_energy(u::FT, v::FT, w::FT) where {FT}
    return FT(0.5) * (u^2 + v^2 + w^2)
end


function liquid_ice_pottemp_given_pressure(
    param_set::APS,
    T::FT,
    p::FT,
    q::TD.PhasePartition{FT} = TD.q_pt_0(FT),
) where {FT <: Real}
    # liquid-ice potential temperature, approximating latent heats
    # of phase transitions as constants
    return TD.dry_pottemp_given_pressure(param_set, T, p, q) *
           exp(- TD.latent_heat_liq_ice(param_set, q) / (TD.cp_m(param_set, q) * T))
end

function liquid_ice_pottemp(
    param_set::APS,
    T::FT,
    ρ::FT,
    q::TD.PhasePartition{FT} = TD.q_pt_0(FT),
) where {FT <: Real}
    return liquid_ice_pottemp_given_pressure(
        param_set,
        T,
        TD.air_pressure(param_set, T, ρ, q),
        q,
    )
end

liquid_ice_pottemp(param_set::APS, ts) = liquid_ice_pottemp(
    param_set,
    TD.air_temperature(param_set, ts),
    TD.air_density(param_set, ts),
    TD.PhasePartition(param_set, ts),
)
