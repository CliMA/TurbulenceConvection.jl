function initialize_io(precip::PrecipVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "rwp_mean")
    add_ts(Stats, "cutoff_precipitation_rate")
    return
end

function io(
    precip::PrecipVariables,
    grid,
    state,
    Stats::NetCDFIO_Stats,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
    TS::TimeStepping,
)
    precipitation_diagnostics(precip, grid, state, up_thermo, en_thermo)
    write_ts(Stats, "rwp_mean", precip.mean_rwp)

    #TODO - change to rain rate that depends on rain model choice
    write_ts(Stats, "cutoff_precipitation_rate", precip.cutoff_precipitation_rate)
    return
end

function precipitation_diagnostics(
    precip::PrecipVariables,
    grid,
    state,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
)
    ρ0_c = center_ref_state(state).ρ0
    prog_pr = center_prog_precipitation(state)

    precip.mean_rwp = 0.0
    precip.cutoff_precipitation_rate = 0.0

    @inbounds for k in real_center_indices(grid)
        precip.mean_rwp += ρ0_c[k] * prog_pr.qr[k] * grid.Δz

        # precipitation rate from cutoff microphysics scheme defined as a total amount of removed water
        # per timestep per EDMF surface area [mm/h]
        if (precip.precipitation_model == "cutoff")
            precip.cutoff_precipitation_rate -=
                (en_thermo.qt_tendency_rain_formation[k] + up_thermo.qt_tendency_rain_formation_tot[k]) *
                ρ0_c[k] *
                grid.Δz / rho_cloud_liq *
                3.6 *
                1e6
        end
    end
    return
end

"""
Computes the qr advection (down) tendency
"""
function compute_rain_advection_tendencies(precip::PrecipPhysics, grid, state, gm, TS::TimeStepping)
    param_set = parameter_set(gm)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = 0.5

    ρ_0_c = center_ref_state(state).ρ0
    tendencies_pr = center_tendencies_precipitation(state)
    prog_pr = center_prog_precipitation(state)

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel = center_field(grid)
    @inbounds for k in real_center_indices(grid)
        term_vel[k] = CM1.terminal_velocity(param_set, rain_type, ρ_0_c[k], prog_pr.qr[k])
    end

    @inbounds for k in reverse(real_center_indices(grid))
        # check stability criterion
        CFL_out = Δt / Δz * term_vel[k]
        if is_toa_center(grid, k)
            CFL_in = 0.0
        else
            CFL_in = Δt / Δz * term_vel[k + 1]
        end
        if max(CFL_in, CFL_out) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end

        ρ_0_cut = ccut_downwind(ρ_0_c, grid, k)
        QR_cut = ccut_downwind(prog_pr.qr, grid, k)
        w_cut = ccut_downwind(term_vel, grid, k)
        ρQRw_cut = ρ_0_cut .* QR_cut .* w_cut
        ∇ρQRw = c∇_downwind(ρQRw_cut, grid, k; bottom = FreeBoundary(), top = SetValue(0))

        ρ_0_c_k = ρ_0_c[k]

        # TODO - some positivity limiters are needed
        precip.qr_tendency_advection[k] = ∇ρQRw / ρ_0_c_k
        tendencies_pr.qr[k] += precip.qr_tendency_advection[k]
    end
    return
end

"""
Computes the tendencies to θ_liq_ice, qt and qr due to rain evaporation
"""
function compute_rain_evap_tendencies(precip::PrecipPhysics, grid, state, gm, TS::TimeStepping)
    param_set = parameter_set(gm)
    Δt = TS.dt
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    @inbounds for k in real_center_indices(grid)
        q_tot_gm = prog_gm.q_tot[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = TD.PhaseEquil_pTq(param_set, p0_c[k], T_gm, q_tot_gm)
        q = TD.PhasePartition(ts)

        # TODO - move limiters elsewhere
        qt_tendency =
            min(prog_pr.qr[k] / Δt, -CM1.evaporation_sublimation(param_set, rain_type, q, prog_pr.qr[k], ρ0_c[k], T_gm))

        precip.qt_tendency_rain_evap[k] = qt_tendency
        precip.qr_tendency_rain_evap[k] = -qt_tendency
        tendencies_pr.qr[k] += precip.qr_tendency_rain_evap[k]
        precip.θ_liq_ice_tendency_rain_evap[k] = θ_liq_ice_helper(ts, qt_tendency)
    end
    return
end

"""
Updates qr based on all microphysics tendencies
"""
function update_precipitation(
    precip_var::PrecipVariables,
    grid,
    state,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
    precip_phys::PrecipPhysics,
    TS::TimeStepping,
)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)
    @inbounds for k in real_center_indices(grid)
        prog_pr.qr[k] += tendencies_pr.qr[k] * TS.dt
    end
    return
end
