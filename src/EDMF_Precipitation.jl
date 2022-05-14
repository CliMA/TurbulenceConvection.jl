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

    tendencies_pr = center_tendencies_precipitation(state)
    prog_pr = center_prog_precipitation(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ

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
    @. aux_tc.qr_tendency_advection = ∇(wvec(RB(ρ_c * q_rai * term_vel_rain))) / ρ_c
    @. aux_tc.qs_tendency_advection = ∇(wvec(RB(ρ_c * q_sno * term_vel_snow))) / ρ_c

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
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    ρ_c = prog_gm.ρ
    tendencies_pr = center_tendencies_precipitation(state)
    ts_gm = aux_gm.ts

    @inbounds for k in real_center_indices(grid)
        qr = prog_pr.q_rai[k]
        qs = prog_pr.q_sno[k]
        ρ = ρ_c[k]
        q_tot_gm = aux_gm.q_tot[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = ts_gm[k]
        q = TD.PhasePartition(param_set, ts)
        qv = TD.vapor_specific_humidity(param_set, ts)

        Π_m = TD.exner(param_set, ts)
        c_pm = TD.cp_m(param_set, ts)
        c_vm = TD.cv_m(param_set, ts)
        R_m = TD.gas_constant_air(param_set, ts)
        R_v = ICP.R_v(param_set)
        L_v0 = ICP.LH_v0(param_set)
        L_s0 = ICP.LH_s0(param_set)
        L_v = TD.latent_heat_vapor(param_set, ts)
        L_s = TD.latent_heat_sublim(param_set, ts)
        L_f = TD.latent_heat_fusion(param_set, ts)

        α_evp = ICP.microph_scaling(param_set)
        α_dep_sub = ICP.microph_scaling_dep_sub(param_set)
        α_melt = ICP.microph_scaling_melt(param_set)

        # TODO - move limiters elsewhere
        # TODO - when using adaptive timestepping we are limiting the source terms
        #        with the previous timestep dt
        S_qr_evap = -min(qr / Δt, -α_evp * CM1.evaporation_sublimation(param_set, rain_type, q, qr, ρ, T_gm))
        S_qs_melt = -min(qs / Δt, α_melt * CM1.snow_melt(param_set, qs, ρ, T_gm))
        tmp = α_dep_sub * CM1.evaporation_sublimation(param_set, snow_type, q, qs, ρ, T_gm)
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
        aux_tc.θ_liq_ice_tendency_precip_sinks[k] =
            1 / Π_m / c_pm * (
                S_qr_evap * (L_v - R_v * T_gm) * (1 + R_m / c_vm) +
                S_qs_sub_dep * (L_s - R_v * T_gm) * (1 + R_m / c_vm) +
                S_qs_melt * L_f * (1 + R_m / c_vm)
            )
    end
    return nothing
end
