function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT) where {FT}
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    return TD.PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
end
function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    config = ()
    q = TD.PhasePartition(q_tot, q_liq, q_ice)

    # eq_state     = TD.PhaseEquil_pθq(   param_set, p, θ_liq_ice, q_tot,config... )
    # non_eq_state = TD.PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q, config...)

    # @info "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    # @info "eq\\_state: `$(eq_state.ρ)`, `$(eq_state.e_int)`, `$((eq_state.q_tot, TD.liquid_specific_humidity(param_set, eq_state), TD.ice_specific_humidity(param_set, eq_state)))`, `$(TD.air_temperature(param_set, eq_state))`, `$(TD.liquid_ice_pottemp(param_set, eq_state))`"
    # @info "non\\_eq\\_state: `$(non_eq_state.ρ)`, `$(non_eq_state.e_int)`, `$((non_eq_state.q.tot, non_eq_state.q.liq, non_eq_state.q.ice))`, `$(TD.air_temperature(param_set, non_eq_state))`, `$(TD.liquid_ice_pottemp(param_set, non_eq_state))`"
    # @info "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    
    # return non_eq_state

    return TD.PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q, config...)
end
