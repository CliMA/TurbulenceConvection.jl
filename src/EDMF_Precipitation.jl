"""
Computes the rain and snow advection (down) tendency
"""
compute_precipitation_advection_tendencies(precip_model, edmf, grid, state, gm) = nothing

function compute_precipitation_advection_tendencies(::Clima1M, edmf, grid, state, gm)
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
Computes the tendencies to θ_liq_ice, q_tot, q_rain and q_snow
due to rain evaporation, snow deposition and sublimation and snow melt
"""
compute_precipitation_sink_tendencies(precip_model, grid, state, gm, Δt) = nothing

function compute_precipitation_sink_tendencies(::Clima1M, grid, state, gm, Δt::Real)
    param_set = parameter_set(gm)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    @inbounds for k in real_center_indices(grid)
        qr = prog_pr.q_rai[k]
        qs = prog_pr.q_sno[k]
        p0 = p0_c[k]
        ρ0 = ρ0_c[k]
        q_tot_gm = prog_gm.q_tot[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = TD.PhaseEquil_pTq(param_set, p0, T_gm, q_tot_gm)
        q = TD.PhasePartition(ts)
        qv = TD.vapor_specific_humidity(ts)

        Π_m = TD.exner(ts)
        c_pm = TD.cp_m(ts)
        c_vm = TD.cv_m(ts)
        R_m = TD.gas_constant_air(ts)
        R_v = CPP.R_v(param_set)
        L_v0 = CPP.LH_v0(param_set)
        L_s0 = CPP.LH_s0(param_set)
        L_v = TD.latent_heat_vapor(ts)
        L_s = TD.latent_heat_sublim(ts)
        L_f = TD.latent_heat_fusion(ts)

        # TODO - move limiters elsewhere
        # TODO - when using adaptive timestepping we are limiting the source terms
        #        with the previous timestep dt
        S_qr_evap = -min(qr / Δt, -CM1.evaporation_sublimation(param_set, rain_type, q, qr, ρ0, T_gm))
        S_qs_melt = -min(qs / Δt, CM1.snow_melt(param_set, qs, ρ0, T_gm))
        tmp = CM1.evaporation_sublimation(param_set, snow_type, q, qs, ρ0, T_gm)
        if tmp > 0.0
            S_qs_sub_dep = min(qv / Δt, tmp)
        else
            S_qs_sub_dep = -min(qs / Δt, -tmp)
        end

        aux_tc.qr_tendency_evap[k] = S_qr_evap
        aux_tc.qs_tendency_melt[k] = S_qs_melt
        aux_tc.qs_tendency_dep_sub[k] = S_qs_sub_dep

        tendencies_pr.q_rai[k] += S_qr_evap - S_qs_melt
        tendencies_pr.q_sno[k] += S_qs_sub_dep + S_qs_melt

        aux_tc.qt_tendency_precip_sinks[k] = -S_qr_evap - S_qs_sub_dep
        aux_tc.θ_liq_ice_tendency_precip_sinks[k] =
            1 / Π_m / c_pm * (
                S_qr_evap * (L_v - R_v * T_gm) * (1 + R_m / c_vm) +
                S_qs_sub_dep * (L_s - R_v * T_gm) * (1 + R_m / c_vm) +
                S_qs_melt * L_f * (1.0 + R_m / c_vm)
            )
    end
    return
end
