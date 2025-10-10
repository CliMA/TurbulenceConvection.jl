"""
    compute_precip_fraction

Computes diagnostic precipitation fraction
"""
compute_precip_fraction(edmf::EDMFModel, state::State) = compute_precip_fraction(edmf.precip_fraction_model, state)

compute_precip_fraction(precip_fraction_model::PrescribedPrecipFraction, ::State) =
    precip_fraction_model.prescribed_precip_frac_value

function compute_precip_fraction(precip_fraction_model::DiagnosticPrecipFraction, state::State)
    FT = float_type(state)
    aux_gm = center_aux_grid_mean(state)
    maxcf = maximum(aux_gm.cloud_fraction)
    return max(maxcf, precip_fraction_model.precip_fraction_limiter)
end

"""
Computes the rain and snow advection (down) tendency
"""
compute_precipitation_advection_tendencies(
    ::AbstractPrecipitationModel,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
) = nothing

function compute_precipitation_advection_tendencies(
    ::Clima1M,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
)
    FT = float_type(state)

    N_up = n_updrafts(edmf)

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

    precip_fraction = compute_precip_fraction(edmf, state)

    q_rai = prog_pr.q_rai #./ precip_fraction
    q_sno = prog_pr.q_sno #./ precip_fraction

    If = CCO.DivergenceF2C()
    RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0))) # sedimentation
    LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0))) # updraft
    # ∇ = CCO.DivergenceF2C(; bottom = CCO.Extrapolate())
    ∇ = CCO.DivergenceF2C(;) # extrapolation makes the bottom cell value the same as the one above (preventing buildup at sfc bc 0 gradient...), we don't need that I think? unless it's some stability thing...
    Ic = CCO.InterpolateF2C()

    
    wvec = CC.Geometry.WVector

    # TODO - some positivity limiters are needed
    # @. aux_tc.qr_tendency_advection = ∇(wvec(RB(ρ_c * q_rai * term_vel_rain))) / ρ_c # * precip_fraction
    # @. aux_tc.qs_tendency_advection = ∇(wvec(RB(ρ_c * q_sno * term_vel_snow))) / ρ_c # * precip_fraction



    # Calculate net tendencies allowing for both sedimentation and advection

    # F2Cw::CCO.InterpolateF2C = CCO.InterpolateF2C(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0))) # shouldnt need bcs for interior interp right?
    sedimentation_differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme # reuse from cloud_sedimentation_model since it's all the same user_args option

    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_en_f = face_aux_environment(state)
    aux_bulk_f = face_aux_bulk(state)
    aux_up_f = face_aux_updrafts(state)
    aux_gm = center_aux_grid_mean(state)

    # aux_tc_f = face_aux_turbconv(state) # aux_tc_f.bulk.w should be the same as aux_bulk_f.w, so we don't need this one.


    # w = F2Cw.(aux_bulk_f.w) # this is the bulk wind
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_rai, term_vel_rain, aux_bulk_f.w, aux_bulk.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false)
    mph_sed_q_rai = mph_sed #.* aux_bulk.area # all goes to same bucket
    mph_sed_q_rai .+= mph_sed_other #.* aux_en.area # all goes to same bucket
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_sno, term_vel_snow, aux_bulk_f.w, aux_bulk.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false)
    mph_sed_q_sno = mph_sed #.* aux_bulk.area # all goes to same bucket
    mph_sed_q_sno .+= mph_sed_other #.* aux_en.area # all goes to same bucket


    # [[ ignore this part, I dont think the env response should be so locally constrained or overly important. ]] [[ maybe we should've done the same for liq/ice now that we allow high updraft areas... idk]]
    # @. w = F2Cw(aux_en_f.w) # this is the environment wind 
    # mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_rai, term_vel_rain, w, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = true)
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_rai, term_vel_rain, aux_en_f.w .* 0, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false) # test ignoring the environment downdraft
    mph_sed_q_rai .+= mph_sed #.* aux_en.area
    mph_sed_q_rai .+= mph_sed_other #.* aux_bulk.area # all goes to same bucket [these are fields so idk how this works...]
    # mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_sno, term_vel_snow, w, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = true)
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_sno, term_vel_snow, aux_en_f.w .* 0, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false) # test ignoring the environment downdraft
    mph_sed_q_sno .+= mph_sed #.* aux_en.area # all goes to same bucket
    mph_sed_q_sno .+= mph_sed_other #.* aux_bulk.area # all goes to same bucket

    # these come out of sedimentation sources already area weighted, so we just sum them up.

    # storage
    @. aux_gm.qr_tendency_sedimentation =  ∇(wvec(RB(ρ_c * q_rai * term_vel_rain))) / ρ_c # * precip_fraction
    @. aux_gm.qs_tendency_sedimentation =  ∇(wvec(RB(ρ_c * q_sno * term_vel_snow))) / ρ_c # * precip_fraction

    # @. aux_gm.qr_tendency_vert_adv = mph_sed_q_rai.q_tendency - aux_gm.qr_tendency_sedimentation # take out the sed only part to get the advection part
    # @. aux_gm.qs_tendency_vert_adv = mph_sed_q_sno.q_tendency - aux_gm.qs_tendency_sedimentation # take out the sed only part to get the advection part

    # @. aux_gm.qr_tendency_vert_adv = -∇(wvec(LB(Ic(aux_en_f.w) * ρ_c * aux_en.area * q_rai))) / ρ_c #  [[ i think the env part is fake ... a dryish downdraft shoould take care of it.... ]]
    # @. aux_gm.qs_tendency_vert_adv = -∇(wvec(LB(Ic(aux_en_f.w) * ρ_c * aux_en.area * q_sno))) / ρ_c # [[ i think the env part is fake ... a dryish downdraft shoould take care of it.... ]]
    @. aux_gm.qr_tendency_vert_adv = FT(0)
    @. aux_gm.qs_tendency_vert_adv = FT(0)
    for i in 1:N_up
        @. aux_gm.qr_tendency_vert_adv +=  -∇(wvec(LB(Ic(aux_up_f[i].w) * ρ_c * aux_up[i].area * q_sno))) / ρ_c
        @. aux_gm.qs_tendency_vert_adv +=  -∇(wvec(LB(Ic(aux_up_f[i].w) * ρ_c * aux_up[i].area * q_sno))) / ρ_c
    end

    #

    @. aux_tc.qr_tendency_advection = mph_sed_q_rai.q_tendency + aux_gm.qr_tendency_vert_adv # / precip_fraction
    @. aux_tc.qs_tendency_advection = mph_sed_q_sno.q_tendency + aux_gm.qs_tendency_vert_adv # / precip_fraction


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
    use_fallback_tendency_limiters::Bool,
) = nothing

function compute_precipitation_sink_tendencies(
    ::Clima1M,
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    param_set::APS,
    Δt::Real,
    use_fallback_tendency_limiters::Bool,
)
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params = TCP.microphysics_params(param_set)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    ρ_c = prog_gm.ρ
    tendencies_pr = center_tendencies_precipitation(state)
    ts_gm = aux_gm.ts

    precip_fraction = compute_precip_fraction(edmf, state)

    FT = float_type(state)

    # ptl = edmf.tendency_limiters.precipitation_tendency_limiter
    ptl = get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters)


    @inbounds for k in real_center_indices(grid)
        # TODO - limiters and positivity checks should be done elsewhere

        qr = max(FT(0), prog_pr.q_rai[k]) / precip_fraction
        qs = max(FT(0), prog_pr.q_sno[k]) / precip_fraction
        ρ = ρ_c[k]
        q_tot_gm = aux_gm.q_tot[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = ts_gm[k]
        q = TD.PhasePartition(thermo_params, ts)
        qv = TD.vapor_specific_humidity(thermo_params, ts)

        Π_m = TD.exner(thermo_params, ts)
        c_pm = TD.cp_m(thermo_params, ts)
        c_vm = TD.cv_m(thermo_params, ts)
        R_m = TD.gas_constant_air(thermo_params, ts)
        R_v = TCP.R_v(param_set)
        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)
        L_v = TD.latent_heat_vapor(thermo_params, ts)
        L_s = TD.latent_heat_sublim(thermo_params, ts)
        L_f = TD.latent_heat_fusion(thermo_params, ts)

        I_l = TD.internal_energy_liquid(thermo_params, ts)
        I_i = TD.internal_energy_ice(thermo_params, ts)
        I = TD.internal_energy(thermo_params, ts)
        Φ = geopotential(param_set, grid.zc.z[k])

        α_evp = TCP.microph_scaling(param_set)
        α_dep_sub = TCP.microph_scaling_dep_sub(param_set)
        α_melt = TCP.microph_scaling_melt(param_set)

        # TODO - move limiters elsewhere
        # TODO - when using adaptive timestepping we are limiting the source terms
        #        with the previous timestep dt
        # S_qr_evap = -min(qr / Δt, -α_evp * CM1.evaporation_sublimation(microphys_params, rain_type, q, qr, ρ, T_gm)) * precip_fraction
        S_qr_evap = limit_tendency(ptl, α_evp * CM1.evaporation_sublimation(microphys_params, rain_type, q, qr, ρ, T_gm), qr, Δt) * precip_fraction # i guess rain only evaporates but never condenses?

        # S_qs_melt = -min(qs / Δt, α_melt * CM1.snow_melt(microphys_params, qs, ρ, T_gm)) * precip_fraction
        S_qs_melt = limit_tendency(ptl, -α_melt * CM1.snow_melt(microphys_params, qs, ρ, T_gm), qs, Δt) * precip_fraction

        tmp = α_dep_sub * CM1.evaporation_sublimation(microphys_params, snow_type, q, qs, ρ, T_gm) * precip_fraction
        # Note if T makes a very low excursion (or too high) , while S*G might mostly cancel and give a real result, you can get NaN from 0 * Inf or something, but I think this is just a legitimate model crash bc those temps are super extreme.
        if tmp > 0
            S_qs_sub_dep = -limit_tendency(ptl, -tmp, qv, Δt)
            qvsat_ice = TD.q_vap_saturation_generic(thermo_params, T_gm, ρ, TD.Ice()) # we don't have this stored for grid-mean and we can't calculate form en/up bc it's non-linear...
            δi = qv - qvsat_ice 
            S_qs_sub_dep = -limit_tendency(ptl, -tmp, max(FT(0), δi), Δt) # presumably if tmp > 0 then δi > 0 but can't be too careful
        else
            S_qs_sub_dep = limit_tendency(ptl, tmp, qs, Δt)
        end

        # [[ We could move from precipitation_formation() to here, but the outcome of that is dispatched based on local temperature in the environment and updraft(s) so we leave it there. ]]

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
