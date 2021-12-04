"""
Computes the rain and snow advection (down) tendency
"""
function compute_precipitation_advection_tendencies(edmf, grid, state, gm)
    FT = eltype(grid)

    ρ0_c = center_ref_state(state).ρ0
    tendencies_pr = center_tendencies_precipitation(state)
    prog_pr = center_prog_precipitation(state)
    aux_tc = center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    q_rai = prog_pr.q_rai
    q_sno = prog_pr.q_sno

    If = CCO.DivergenceF2C()
    RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0)))
    ∇ = CCO.DivergenceF2C(; bottom = CCO.Extrapolate())
    wvec = CC.Geometry.WVector

    # TODO - some positivity limiters are needed
    @. aux_tc.qr_tendency_advection = ∇(wvec(RB(ρ0_c * q_rai * term_vel_rain))) / ρ0_c
    @. aux_tc.qs_tendency_advection = ∇(wvec(RB(ρ0_c * q_sno * term_vel_snow))) / ρ0_c

    @. tendencies_pr.q_rai += aux_tc.qr_tendency_advection
    @. tendencies_pr.q_sno += aux_tc.qs_tendency_advection
    return
end

"""
Computes the tendencies to θ_liq_ice, qt and qr due to rain evaporation
"""
function compute_precipitation_sink_tendencies(grid, state, gm, TS::TimeStepping)
    param_set = parameter_set(gm)
    Δt = TS.dt
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
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
        rain_evap =
            -min(
                prog_pr.q_rai[k] / Δt,
                -CM1.evaporation_sublimation(param_set, rain_type, q, prog_pr.q_rai[k], ρ0_c[k], T_gm),
            )

        aux_tc.qr_tendency_evap[k] = rain_evap
        aux_tc.qs_tendency_melt[k] = 0.0
        aux_tc.qs_tendency_dep_sub[k] = 0.0

        aux_tc.qt_tendency_precip_sinks[k] = -rain_evap
        aux_tc.θ_liq_ice_tendency_precip_sinks[k] = θ_liq_ice_helper(ts, -rain_evap)

        tendencies_pr.q_rai[k] += rain_evap
    end
    return
end
