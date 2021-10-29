function satadjust(gm::GridMeanVariables, grid, state)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    param_set = parameter_set(gm)
    @inbounds for k in real_center_indices(grid)
        θ_liq_ice = prog_gm.θ_liq_ice[k]
        q_tot = prog_gm.q_tot[k]
        ts = thermo_state_pθq(param_set, p0_c[k], θ_liq_ice, q_tot)
        aux_gm.q_liq[k] = TD.liquid_specific_humidity(ts)
        aux_gm.q_ice[k] = TD.ice_specific_humidity(ts)
        aux_gm.T[k] = TD.air_temperature(ts)
        ρ = TD.air_density(ts)
        aux_gm.buoy[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
        aux_gm.RH[k] = TD.relative_humidity(ts)
    end
    return
end
