function initialize_io(rain::RainVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "qr_mean")
    add_ts(Stats, "rwp_mean")
    add_ts(Stats, "cutoff_rain_rate")
    return
end

function io(
    rain::RainVariables,
    Stats::NetCDFIO_Stats,
    ref_state::ReferenceState,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
    TS::TimeStepping,
)
    write_profile(Stats, "qr_mean", rain.QR.values)

    rain_diagnostics(rain, ref_state, up_thermo, en_thermo)
    write_ts(Stats, "rwp_mean", rain.mean_rwp)

    #TODO - change to rain rate that depends on rain model choice
    write_ts(Stats, "cutoff_rain_rate", rain.cutoff_rain_rate)
    return
end

function rain_diagnostics(
    rain::RainVariables,
    ref_state::ReferenceState,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
)
    rain.mean_rwp = 0.0
    rain.cutoff_rain_rate = 0.0
    grid = rain.grid

    @inbounds for k in real_center_indices(grid)
        rain.mean_rwp += ref_state.rho0_half[k] * rain.QR.values[k] * grid.Δz

        # rain rate from cutoff microphysics scheme defined as a total amount of removed water
        # per timestep per EDMF surface area [mm/h]
        if (rain.rain_model == "cutoff")
            rain.cutoff_rain_rate -=
                (en_thermo.qt_tendency_rain_formation[k] + up_thermo.qt_tendency_rain_formation_tot[k]) *
                ref_state.rho0_half[k] *
                grid.Δz / rho_cloud_liq *
                3.6 *
                1e6
        end
    end
    return
end

function compute_subdomain_rain_tendency_sum(
    rain::RainVariables,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
    TS::TimeStepping,
)
    # TODO - just move it to the place where all rain tendecies are summed up
    @inbounds for k in real_center_indices(rain.grid)
        rain.QR.values[k] -=
            (en_thermo.qt_tendency_rain_formation[k] + up_thermo.qt_tendency_rain_formation_tot[k]) * TS.dt
    end
    return
end

"""
Computes the qr advection (down) tendency
"""
function compute_rain_advection_tendencies(rain::RainPhysics, gm::GridMeanVariables, TS::TimeStepping, QR::RainVariable)
    param_set = parameter_set(gm)
    grid = get_grid(gm)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = 0.5

    ρ_0_c_field = rain.ref_state.rho0_half

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel = center_field(grid)
    @inbounds for k in real_center_indices(grid)
        term_vel[k] = CM1.terminal_velocity(param_set, rain_type, rain.ref_state.rho0_half[k], QR.values[k])
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

        ρ_0_cut = ccut_downwind(ρ_0_c_field, grid, k)
        QR_cut = ccut_downwind(QR.values, grid, k)
        w_cut = ccut_downwind(term_vel, grid, k)
        ρQRw_cut = ρ_0_cut .* QR_cut .* w_cut
        ∇ρQRw = c∇_downwind(ρQRw_cut, grid, k; bottom = FreeBoundary(), top = SetValue(0))

        ρ_0_c = ρ_0_c_field[k]
        rain.qr_tendency_advection[k] = ∇ρQRw / ρ_0_c

        # TODO - apply the tendecies elsewhere
        # TODO - some positivity limiters are needed
        QR.values[k] += rain.qr_tendency_advection[k] * Δt
    end
    return
end

"""
Computes the tendencies to θ_liq_ice, qt and qr due to rain evaporation
"""
function compute_rain_evap_tendencies(rain::RainPhysics, gm::GridMeanVariables, TS::TimeStepping, QR::RainVariable)
    param_set = parameter_set(gm)
    Δt = TS.dt
    p0_c = rain.ref_state.p0_half
    ρ0_c = rain.ref_state.rho0_half

    @inbounds for k in real_center_indices(rain.grid)
        q_tot_gm = gm.QT.values[k]
        T_gm = gm.T.values[k]
        # When we fuse loops, this should hopefully disappear
        ts = TD.PhaseEquil_pTq(param_set, p0_c[k], T_gm, q_tot_gm)
        q = TD.PhasePartition(ts)

        # TODO - move limiters elsewhere
        qt_tendency =
            min(QR.values[k] / Δt, -CM1.evaporation_sublimation(param_set, rain_type, q, QR.values[k], ρ0_c[k], T_gm))
        rain.qt_tendency_rain_evap[k] = qt_tendency
        rain.qr_tendency_rain_evap[k] = -qt_tendency
        rain.θ_liq_ice_tendency_rain_evap[k] = θ_liq_ice_helper(ts, qt_tendency)

        # TODO - apply tendencies elsewhere
        QR.values[k] += rain.qr_tendency_rain_evap[k] * Δt
    end
    return
end
