
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
                (en_thermo.prec_source_qt[k] + up_thermo.prec_source_qt_tot[k]) * ref_state.rho0_half[k] * grid.Δz /
                rho_cloud_liq *
                3.6 *
                1e6
        end
    end
    return
end

function sum_subdomains_rain(
    rain::RainVariables,
    up_thermo::UpdraftThermodynamics,
    en_thermo::EnvironmentThermodynamics,
    TS::TimeStepping,
)
    @inbounds for k in real_center_indices(rain.grid)
        rain.QR.values[k] -= (en_thermo.prec_source_qt[k] + up_thermo.prec_source_qt_tot[k]) * TS.dt
    end
    return
end

function solve_rain_fall(rain::RainPhysics, gm::GridMeanVariables, TS::TimeStepping, QR::RainVariable)
    param_set = parameter_set(gm)
    grid = get_grid(gm)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = 0.5

    term_vel = center_field(grid)
    term_vel_new = center_field(grid)
    ρ_0_c_field = rain.ref_state.rho0_half

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    @inbounds for k in real_center_indices(grid)
        term_vel[k] = CM1.terminal_velocity(param_set, rain_type, rain.ref_state.rho0_half[k], QR.values[k])
    end

    # rain falling through the domain
    @inbounds for k in reverse(real_center_indices(grid))
        CFL_out = Δt / Δz * term_vel[k]

        if is_toa_center(grid, k)
            CFL_in = 0.0
        else
            CFL_in = Δt / Δz * term_vel[k + 1]
        end

        if max(CFL_in, CFL_out) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end
        ρ_0_c = rain.ref_state.rho0_half[k]

        ρ_0_cut = ccut_downwind(ρ_0_c_field, grid, k)
        QR_cut = ccut_downwind(QR.values, grid, k)
        w_cut = ccut_downwind(term_vel, grid, k)
        ρQRw_cut = ρ_0_cut .* QR_cut .* w_cut
        ∇ρQRw = c∇_downwind(ρQRw_cut, grid, k; bottom = FreeBoundary(), top = SetValue(0))

        # TODO: incorporate area fraction into this equation (right now assume a = 1)
        QR.new[k] = (QR.values[k] * ρ_0_c + Δt * ∇ρQRw) / ρ_0_c
        QR.new[k] = max(QR.new[k], 0)

        term_vel_new[k] = CM1.terminal_velocity(param_set, rain_type, ρ_0_c, QR.new[k])
    end

    QR.values .= QR.new

    term_vel .= term_vel_new
    return
end

function solve_rain_evap(rain::RainPhysics, gm::GridMeanVariables, TS::TimeStepping, QR::RainVariable)
    param_set = parameter_set(gm)
    Δt = TS.dt
    flag_evaporate_all = false
    p0_c = rain.ref_state.p0_half
    ρ0_c = rain.ref_state.rho0_half

    @inbounds for k in real_center_indices(rain.grid)
        flag_evaporate_all = false

        q = TD.PhasePartition(gm.QT.values[k], gm.QL.values[k], 0.0)

        tmp_evap_rate = -CM1.evaporation_sublimation(param_set, rain_type, q, QR.values[k], ρ0_c[k], gm.T.values[k])

        if tmp_evap_rate * Δt > QR.values[k]
            flag_evaporate_all = true
            tmp_evap_rate = QR.values[k] / Δt
        end

        rain.rain_evap_source_qt[k] = tmp_evap_rate

        # TODO add ice
        rain_source_to_thetal(param_set, p0_c[k], gm.T.values[k], gm.QT.values[k], gm.QL.values[k], 0.0, -tmp_evap_rate)

        if flag_evaporate_all
            QR.values[k] = 0.0
        else
            QR.values[k] -= tmp_evap_rate * Δt
        end

    end
    return
end
