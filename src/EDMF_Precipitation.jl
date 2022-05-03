"""
Computes the rain and snow advection (down) tendency
"""
compute_precipitation_advection_tendencies(
    ::AbstractPrecipitationModel,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
) = nothing

function compute_precipitation_advection_tendencies(
    ::Clima1M,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
)
    FT = eltype(grid)

    ρ0_c = center_ref_state(state).ρ0
    tendencies_pr = center_tendencies_precipitation(state)
    prog_pr = center_prog_precipitation(state)
    aux_tc = center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
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
    return nothing
end

"""
Computes the tendencies to θ_liq_ice, q_tot, q_rain and q_snow
due to rain evaporation, snow deposition and sublimation and snow melt
"""
compute_precipitation_sink_tendencies(
    ::AbstractPrecipitationModel,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
    Δt::Real,
) = nothing

function compute_precipitation_sink_tendencies(
    ::Clima1M,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
    Δt::Real,
)
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)
    ts_gm = aux_gm.ts

    @inbounds for k in real_center_indices(grid)
        qr = prog_pr.q_rai[k]
        qs = prog_pr.q_sno[k]
        ρ0 = ρ0_c[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = ts_gm[k]
        q = TD.PhasePartition(param_set, ts)
        qv = TD.vapor_specific_humidity(param_set, ts)

        I_l = TD.internal_energy_liquid(ts)
        I_i = TD.internal_energy_ice(ts)
        #Φ = gravitational_potential(atmos.orientation, aux) #TODO how to use it here?
        g = CPP.grav(param_set)
        Φ = g * grid.zc[k]
        Lf = TD.latent_heat_fusion(param_set, ts)

        α_evp = CPMP.microph_scaling(param_set)
        α_dep_sub = CPMP.microph_scaling_dep_sub(param_set)
        α_melt = CPMP.microph_scaling_melt(param_set)

        # TODO - move limiters elsewhere
        # TODO - when using adaptive timestepping we are limiting the source terms
        #        with the previous timestep dt
        S_qr_evap = -min(qr / Δt, -α_evp * CM1.evaporation_sublimation(param_set, rain_type, q, qr, ρ0, T_gm))
        S_qs_melt = -min(qs / Δt, α_melt * CM1.snow_melt(param_set, qs, ρ0, T_gm))
        tmp = α_dep_sub * CM1.evaporation_sublimation(param_set, snow_type, q, qs, ρ0, T_gm)
        if tmp > 0
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
        aux_tc.e_tot_tendency_precip_sinks[k] = -S_qr_evap * (I_l + Φ) - S_qs_sub_dep * (I_i + Φ) + S_qs_melt * Lf
    end
    return nothing
end
