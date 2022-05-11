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
    grav = CPP.grav(param_set)
    return grav * z
end

function anelastic_total_enthalpy(param_set::APS, e_tot, ts_gm) where {FT}
    Rm_gm = TD.gas_constant_air(param_set, ts_gm)
    T_gm = TD.air_temperature(param_set, ts_gm)
    return e_tot + Rm_gm * T_gm
end

function anelastic_total_enthalpy(param_set::APS, e_tot, p, ρ) where {FT}
    return e_tot + p / ρ
end
