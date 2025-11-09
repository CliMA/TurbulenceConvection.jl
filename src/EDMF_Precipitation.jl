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
    state::State,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
) = nothing

function compute_precipitation_advection_tendencies(
    ::Clima1M,
    edmf::EDMFModel,
    state::State,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
)
    FT = float_type(state)
    grid = Grid(state)

    N_up = n_updrafts(edmf)

    tendencies_pr = center_tendencies_precipitation(state)
    prog_pr = center_prog_precipitation(state)
    aux_tc = center_aux_turbconv(state)
    aux_tc_f = face_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)

    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_en_f = face_aux_environment(state)
    aux_bulk_f = face_aux_bulk(state)
    aux_up_f = face_aux_updrafts(state)
    aux_gm = center_aux_grid_mean(state)

    ρ_c = prog_gm.ρ

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    # precip_fraction = compute_precip_fraction(edmf, state)

    q_rai = prog_pr.q_rai #./ precip_fraction
    q_sno = prog_pr.q_sno #./ precip_fraction

    RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0))) # sedimentation
    LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0))) # updraft
    # ∇ = CCO.DivergenceF2C(; bottom = CCO.Extrapolate())
    ∇ = CCO.DivergenceF2C(;) # extrapolation makes the bottom cell value the same as the one above (preventing buildup at sfc bc 0 gradient...), we don't need that I think? unless it's some stability thing...

    
    wvec = CC.Geometry.WVector

    # TODO - some positivity limiters are needed
    # @. aux_tc.qr_tendency_advection = ∇(wvec(RB(ρ_c * q_rai * term_vel_rain))) / ρ_c # * precip_fraction
    # @. aux_tc.qs_tendency_advection = ∇(wvec(RB(ρ_c * q_sno * term_vel_snow))) / ρ_c # * precip_fraction



    # Calculate net tendencies allowing for both sedimentation and advection

    sedimentation_differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme # reuse from cloud_sedimentation_model since it's all the same user_args option



    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_rai, term_vel_rain, aux_bulk_f.w, aux_bulk.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false, scratch1 = aux_tc.temporary_1, scratch2 = aux_tc.temporary_2, scratch3 = aux_tc.temporary_3, scratch4 = aux_tc.temporary_4, scratch1F = aux_tc_f.temporary_f1, scratch2F = aux_tc_f.temporary_f2)
    mph_sed_q_rai = mph_sed #.* aux_bulk.area # all goes to same bucket
    mph_sed_q_rai .+= mph_sed_other #.* aux_en.area # all goes to same bucket
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_sno, term_vel_snow, aux_bulk_f.w, aux_bulk.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false, scratch1 = aux_tc.temporary_1, scratch2 = aux_tc.temporary_2, scratch3 = aux_tc.temporary_3, scratch4 = aux_tc.temporary_4, scratch1F = aux_tc_f.temporary_f1, scratch2F = aux_tc_f.temporary_f2)
    mph_sed_q_sno = mph_sed #.* aux_bulk.area # all goes to same bucket
    mph_sed_q_sno .+= mph_sed_other #.* aux_en.area # all goes to same bucket


    # [[ ignore this part, I dont think the env response should be so locally constrained or overly important. ]] [[ maybe we should've done the same for liq/ice now that we allow high updraft areas... idk]]
    # mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_rai, term_vel_rain, w, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = true)
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_rai, term_vel_rain, aux_en_f.w .* 0, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false, scratch1 = aux_tc.temporary_1, scratch2 = aux_tc.temporary_2, scratch3 = aux_tc.temporary_3, scratch4 = aux_tc.temporary_4, scratch1F = aux_tc_f.temporary_f1, scratch2F = aux_tc_f.temporary_f2) # test ignoring the environment downdraft
    mph_sed_q_rai .+= mph_sed #.* aux_en.area
    mph_sed_q_rai .+= mph_sed_other #.* aux_bulk.area # all goes to same bucket [these are fields so idk how this works...]
    # mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_sno, term_vel_snow, w, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = true)
    mph_sed, mph_sed_other = calculate_sedimentation_sources(param_set, ρ_c, q_sno, term_vel_snow, aux_en_f.w .* 0, aux_en.area, grid; differencing_scheme = sedimentation_differencing_scheme, grid_mean = false, use_relative_w = false, scratch1 = aux_tc.temporary_1, scratch2 = aux_tc.temporary_2, scratch3 = aux_tc.temporary_3, scratch4 = aux_tc.temporary_4, scratch1F = aux_tc_f.temporary_f1, scratch2F = aux_tc_f.temporary_f2) # test ignoring the environment downdraft
    mph_sed_q_sno .+= mph_sed #.* aux_en.area # all goes to same bucket
    mph_sed_q_sno .+= mph_sed_other #.* aux_bulk.area # all goes to same bucket

    # these come out of sedimentation sources already area weighted, so we just sum them up.

    # storage
    @. aux_gm.qr_tendency_sedimentation =  ∇(wvec(RB(ρ_c * q_rai * term_vel_rain))) / ρ_c # * precip_fraction
    @. aux_gm.qs_tendency_sedimentation =  ∇(wvec(RB(ρ_c * q_sno * term_vel_snow))) / ρ_c # * precip_fraction

    # @. aux_gm.qr_tendency_vert_adv = mph_sed_q_rai.q_tendency - aux_gm.qr_tendency_sedimentation # take out the sed only part to get the advection part
    # @. aux_gm.qs_tendency_vert_adv = mph_sed_q_sno.q_tendency - aux_gm.qs_tendency_sedimentation # take out the sed only part to get the advection part


    @. aux_gm.qr_tendency_vert_adv = FT(0)
    @. aux_gm.qs_tendency_vert_adv = FT(0)


    if edmf.area_partition_model isa CoreCloakAreaPartitionModel
        #=
        unlike with up/en, we can't just use f_cm directly to mix an updraft and an environment.
        We'd like to partition q_rai and q_sno into the updraft and environment, e.g. q_rai_up and q_rai_en, such that q_rai_up * a_up_tot + q_rai_en a_en = q_rai  ... However there are infinitely many pairs that solve this so we're kinda stuck
        So instead we pivot. We assume the updraft and updraft cloak will have the same qr, and downdraft cloak will close the budget with remaining env unchanged...
        So for example,  want the value q_rai in a_en_remaining = a_en - a_cloak_up - a_cloak_dn to be unchanged from q_rai, ditto for snow
        So, for example, for liquid need to solve for q_rai_up and q_rai_en, and then q_rai_cloak_up and q_rai_cloak_dn
        We can't really do that easily, so we'll just put all the precip in the updraft and cloak up and leave nothing the downdraft, with the remaining env unchanged.

        So we partition as q_rai - aux_en.a_en_remaining * q_rai, leaving nothing for cloak_dn

        I think in reality we should do what we did from the `precipitation_sink_tendencies()` fcn... otherwise we will risk the same fate of being overly at the whims of updraft area...
        =#
        # f_cm = edmf.area_partition_model.cloak_mix_factor #  Not using this biasing and relying on LB and RB wasn't unstable per se, but it drastically undercounts vert adv

        # Basically in the updraft you'd have q_rai_up * a_up + ((f_cm * q_rai_up) + (1-f_cm) * q_rai_en)

        # ============================================================== #
            # Partition same way as `compute_precipitation_sink_tendencies()`
        # ============================================================== #
        # qs_up :: ifelse(aux_gm.q_ice > FT(0), (aux_bulk.q_ice * aux_bulk.area) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Updraft
        # qs_cup :: ifelse(aux_gm.q_ice > FT(0), (aux_en.q_ice_cloak_up * aux_en.a_cloak_up) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Updraft cloak
        # qs_cdn :: ifelse(aux_gm.q_ice > FT(0), (aux_en.q_ice_cloak_dn * aux_en.a_cloak_dn) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Downdraft cloak
        # qs_enr :: ifelse(aux_gm.q_ice > FT(0), (aux_en.q_ice_en_remaining * aux_en.a_en_remaining) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Env remaining
 
        # qr_up :: ifelse(aux_gm.q_liq > FT(0), (aux_bulk.q_liq * aux_bulk.area) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Updraft
        # qr_cup :: ifelse(aux_gm.q_liq > FT(0), (aux_en.q_liq_cloak_up * aux_en.a_cloak_up) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Updraft cloak
        # qr_cdn :: ifelse(aux_gm.q_liq > FT(0), (aux_en.q_liq_cloak_dn * aux_en.a_cloak_dn) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Downdraft cloak
        # qr_enr :: ifelse(aux_gm.q_liq > FT(0), (aux_en.q_liq_en_remaining * aux_en.a_en_remaining) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Env remaining

        q_rai_here = aux_tc.temporary_1
        q_sno_here = aux_tc.temporary_2

        for i in 1:N_up
            @. q_rai_here = ifelse(aux_gm.q_liq > FT(0), (aux_up[i].q_liq * aux_up[i].area) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Updraft
            @. q_sno_here = ifelse(aux_gm.q_ice > FT(0), (aux_up[i].q_ice * aux_up[i].area) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Updraft
            @. aux_gm.qr_tendency_vert_adv += -∇(wvec(LB(Ic(aux_up_f[i].w) * ρ_c * aux_up[i].area * q_rai_here))) / ρ_c
            @. aux_gm.qs_tendency_vert_adv += -∇(wvec(LB(Ic(aux_up_f[i].w) * ρ_c * aux_up[i].area * q_sno_here))) / ρ_c
        end
        
        if edmf.area_partition_model.confine_all_downdraft_to_cloak # left of right biased based on direction. This does have the risk to be unstable though... hopefully running out of w helps moderate.... otherwise we could try the second order correction lol...
            regions = (
                # (aux_en.a_cloak_up, aux_en_f.w_cloak_up, qr_cup, qs_cup, LB),
                (aux_en.a_cloak_up, aux_en_f.w_cloak_up, LB, CloakUpDomain),
                # (aux_en.a_cloak_dn, aux_en_f.w_cloak_dn, qr_cdn, qs_cdn, RB),
                (aux_en.a_cloak_dn, aux_en_f.w_cloak_dn, RB, CloakDownDomain),
                )
        else
            regions = (
                # (aux_en.a_cloak_up, aux_en_f.w_cloak_up, qr_cup, qs_cup, LB),
                (aux_en.a_cloak_up, aux_en_f.w_cloak_up, LB, CloakUpDomain),
                # (aux_en.a_cloak_dn, aux_en_f.w_cloak_dn, qr_cdn, qs_cdn, RB),
                (aux_en.a_cloak_dn, aux_en_f.w_cloak_dn, RB, CloakDownDomain),
                # (aux_en.a_en_remaining, aux_en_f.w, qr_enr, qs_enr, RB),
                (aux_en.a_en_remaining, aux_en_f.w, RB, EnvRemainingDomain),
                )
        end

        # for (area_region, w_region, q_rai_here, q_sno_here, LRB) in regions
        for (area_region, w_region, LRB, domain_type) in regions
            if domain_type === CloakUpDomain
                @. q_rai_here = ifelse(aux_gm.q_liq > FT(0), aux_en.q_liq_cloak_up * aux_en.a_cloak_up / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Updraft
                @. q_sno_here = ifelse(aux_gm.q_ice > FT(0), aux_en.q_ice_cloak_up * aux_en.a_cloak_up / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Updraft
            elseif domain_type === CloakDownDomain
                @. q_rai_here = ifelse(aux_gm.q_liq > FT(0), (aux_en.q_liq_cloak_dn * aux_en.a_cloak_dn) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Downdraft cloak
                @. q_sno_here = ifelse(aux_gm.q_ice > FT(0), (aux_en.q_ice_cloak_dn * aux_en.a_cloak_dn) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Downdraft cloak
            elseif domain_type === EnvRemainingDomain
                @. q_rai_here = ifelse(aux_gm.q_liq > FT(0), (aux_en.q_liq_en_remaining * aux_en.a_en_remaining) / aux_gm.q_liq , FT(1)) * q_rai # fraction of liquid in Env remaining
                @. q_sno_here = ifelse(aux_gm.q_ice > FT(0), (aux_en.q_ice_en_remaining * aux_en.a_en_remaining) / aux_gm.q_ice , FT(1)) * q_sno # fraction of ice in Env remaining
            end
            @. aux_gm.qr_tendency_vert_adv += -∇(wvec(LRB(Ic(w_region) * ρ_c * area_region * q_rai_here))) / ρ_c
            @. aux_gm.qs_tendency_vert_adv += -∇(wvec(LRB(Ic(w_region) * ρ_c * area_region * q_sno_here))) / ρ_c
        end

    else
        for i in 1:N_up
            @. aux_gm.qr_tendency_vert_adv += -∇(wvec(LB(Ic(aux_up_f[i].w) * ρ_c * aux_up[i].area * q_rai))) / ρ_c
            @. aux_gm.qs_tendency_vert_adv += -∇(wvec(LB(Ic(aux_up_f[i].w) * ρ_c * aux_up[i].area * q_sno))) / ρ_c
        end

        if param_set.user_params.use_convective_tke
            w_conv = aux_tc_f.temporary_f1
            @. w_conv = sqrt(max(Ifx(2 * aux_en.tke_convective), FT(0))) # could add a limiter here if needed
            @. aux_gm.qr_tendency_vert_adv += -∇(wvec(LB(Ic(w_conv) * ρ_c * aux_en.area/2 * q_rai))) / ρ_c
            @. aux_gm.qs_tendency_vert_adv += -∇(wvec(LB(Ic(w_conv) * ρ_c * aux_en.area/2 * q_sno))) / ρ_c    
            # Ignore tke down draft for now...
        else
            # ignore env contribution for now....
            # @. aux_gm.qr_tendency_vert_adv = -∇(wvec(LB(Ic(aux_en_f.w) * ρ_c * aux_en.area * q_rai))) / ρ_c #  [[ i think the env part is fake ... a dryish downdraft shoould take care of it.... ]]
            # @. aux_gm.qs_tendency_vert_adv = -∇(wvec(LB(Ic(aux_en_f.w) * ρ_c * aux_en.area * q_sno))) / ρ_c # [[ i think the env part is fake ... a dryish downdraft shoould take care of it.... ]]
        end

    end


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
    state::State,
    param_set::APS,
    Δt::Real,
    use_fallback_tendency_limiters::Bool,
) = nothing

function compute_precipitation_sink_tendencies(
    ::Clima1M,
    edmf::EDMFModel,
    state::State,
    param_set::APS,
    Δt::Real,
    use_fallback_tendency_limiters::Bool,
)
    grid = Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params = TCP.microphysics_params(param_set)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    ρ_c = prog_gm.ρ
    p_c = aux_gm.p
    tendencies_pr = center_tendencies_precipitation(state)
    ts_gm = aux_gm.ts
    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)

    precip_fraction = compute_precip_fraction(edmf, state)

    FT = float_type(state)

    # ptl = edmf.tendency_limiters.precipitation_tendency_limiter
    ptl = get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters)


    @inbounds for k in real_center_indices(grid)
        # TODO - limiters and positivity checks should be done elsewhere

        qr = max(FT(0), prog_pr.q_rai[k]) / precip_fraction
        qs = max(FT(0), prog_pr.q_sno[k]) / precip_fraction
        # ρ = ρ_c[k]
        q_tot_gm = aux_gm.q_tot[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = ts_gm[k]
        # q = TD.PhasePartition(thermo_params, ts)
        qv = TD.vapor_specific_humidity(thermo_params, ts)
        # qvsat_ice_gm = qv - TD.q_vap_saturation_generic(thermo_params, T_gm, ρ_c[k], TD.Ice())
        Π_m = TD.exner(thermo_params, ts)
        c_pm = TD.cp_m(thermo_params, ts)
        c_vm = TD.cv_m(thermo_params, ts)
        c_p = TD.TP.cp_d(thermo_params)
        R_m = TD.gas_constant_air(thermo_params, ts)
        R_v = TCP.R_v(param_set)
        L_v = TD.latent_heat_vapor(thermo_params, ts)
        L_s = TD.latent_heat_sublim(thermo_params, ts)
        L_f = TD.latent_heat_fusion(thermo_params, ts)

        # I_l = TD.internal_energy_liquid(thermo_params, ts)
        # I_i = TD.internal_energy_ice(thermo_params, ts)
        # I = TD.internal_energy(thermo_params, ts)
        Φ = geopotential(param_set, grid.zc.z[k])

        α_evp = TCP.microph_scaling(param_set)
        α_dep_sub = TCP.microph_scaling_dep_sub(param_set)
        α_melt = TCP.microph_scaling_melt(param_set)

        # TODO - move limiters elsewhere
        # TODO - when using adaptive timestepping we are limiting the source terms
        #        with the previous timestep dt

        # we could try to constrain lambda for ice to respect atleast N_INP_top - N_ice but... 1) too complex 2) we don't currently store N_INP_top 3) it wouldnt be clear how to blend the top for env and updrafts
        # we could at least enforce at least local N_INP - N_i_mean? idk but with boost it's almost certainly above that...


        w0 = FT(0) # use this instead of doing the face to center conversion, it would only be used for getting NINP
        if edmf.area_partition_model isa CoreCloakAreaPartitionModel # Igore w for now
            # As with advection, we will again partition everything outside of en_remaining into the updrafts and cloak up, leaving nothing in cloak_dn, and en_remainig unchanged
            #=
            [[ I think this is too strong. The problem is that variability in updraft area can lead to wild swings in overall sub/dep rate. We had tried e.g. initializing with large updraft area and that was disastrous. at a=0, sure it's just env, but in the middle with large cloaks you bias highly to only growth. It's worse than doubling the rate because you lose the offsetting sum if a_en_remaining is very small]]
            [[ Instead, we'll just keep the ratio 1:1 with the ratios of q_liq, q_ice as a heuristic. It gels w/ local autoconv as well, esp dep driven.. though becomes less accurate with further growth and sedimentation... ]]
            =#
            # combined_up_area = aux_bulk.area[k] + aux_en.a_cloak_up[k]
            # f_q = (combined_up_area > sqrt(eps(FT))) ? ((one(FT) - aux_en.a_en_remaining[k]) / combined_up_area) : one(FT) # how much q is enhanced in updraft and cloak up by squeezing out of downdraft cloak
            # f_q = FT(1) # testing no biasing for now


            f_qi_up = (aux_gm.q_ice[k] > FT(0)) ? (aux_bulk.q_ice[k] * aux_bulk.area[k]) / aux_gm.q_ice[k] : FT(1) # fraction of ice in Updraft
            f_qi_cup = (aux_gm.q_ice[k] > FT(0)) ? (aux_en.q_ice_cloak_up[k] * aux_en.a_cloak_up[k]) / aux_gm.q_ice[k] : FT(1) # fraction of ice in Updraft cloak
            f_qi_cdn = (aux_gm.q_ice[k] > FT(0)) ? (aux_en.q_ice_cloak_dn[k] * aux_en.a_cloak_dn[k]) / aux_gm.q_ice[k] : FT(1) # fraction of ice in Downdraft cloak
            f_qi_enr = (aux_gm.q_ice[k] > FT(0)) ? (aux_en.q_ice_en_remaining[k] * aux_en.a_en_remaining[k]) / aux_gm.q_ice[k] : FT(1) # fraction of ice in Env remaining

            f_ql_up = (aux_gm.q_liq[k] > FT(0)) ? (aux_bulk.q_liq[k] * aux_bulk.area[k]) / aux_gm.q_liq[k] : FT(1) # fraction of liquid in Updraft
            f_ql_cup = (aux_gm.q_liq[k] > FT(0)) ? (aux_en.q_liq_cloak_up[k] * aux_en.a_cloak_up[k]) / aux_gm.q_liq[k] : FT(1) # fraction of liquid in Updraft cloak
            f_ql_cdn = (aux_gm.q_liq[k] > FT(0)) ? (aux_en.q_liq_cloak_dn[k] * aux_en.a_cloak_dn[k]) / aux_gm.q_liq[k] : FT(1) # fraction of liquid in Downdraft cloak
            f_ql_enr = (aux_gm.q_liq[k] > FT(0)) ? (aux_en.q_liq_en_remaining[k] * aux_en.a_en_remaining[k]) / aux_gm.q_liq[k] : FT(1) # fraction of liquid in Env remaining


            if edmf.moisture_model isa NonEquilibriumMoisture
                ts_bulk = thermo_state_pθq(param_set, p_c[k], aux_bulk.θ_liq_ice[k], aux_bulk.q_tot[k], aux_bulk.q_liq[k], aux_bulk.q_ice[k])
            else
                ts_bulk = thermo_state_pθq(param_set, p_c[k], aux_bulk.θ_liq_ice[k], aux_bulk.q_tot[k])
            end
            regions = (
                (aux_bulk.area[k], aux_bulk.T[k], ts_bulk, w0, aux_bulk.N_i[k], TD.vapor_specific_humidity(thermo_params, ts_bulk), f_ql_up*qr, f_qi_up*qs, true), # we don't store bulk ts, so just use grid mean ts (we only really need density)
                (aux_en.a_cloak_up[k], aux_en.T_cloak_up[k], aux_en.ts_cloak_up[k], w0, aux_en.N_i[k], TD.vapor_specific_humidity(thermo_params, aux_en.ts_cloak_up[k]), f_ql_cup*qr, f_qi_cup*qs, true),
                (aux_en.a_cloak_dn[k], aux_en.T_cloak_dn[k], aux_en.ts_cloak_dn[k], w0, aux_en.N_i[k], TD.vapor_specific_humidity(thermo_params, aux_en.ts_cloak_dn[k]), f_ql_cdn*qr, f_qi_cdn*qs, false), # all zeroed out
                (aux_en.a_en_remaining[k], aux_en.T[k], aux_en.ts[k], w0, aux_en.N_i[k], TD.vapor_specific_humidity(thermo_params, aux_en.ts[k]), f_ql_enr*qr, f_qi_enr*qs, false),
                # (FT(1), T_gm, ts, w0, aux_gm.N_i[k], qv, qr, qs, false), # testing
            )
        else
            regions = (
                (FT(1), T_gm, ts, w0, aux_gm.N_i[k], qv, qr, qs, false), 
                ) # we don't store gm w, so just use 0
        end


        S_qr_evap = FT(0)
        S_qs_melt = FT(0)
        S_qs_sub_dep = FT(0)
        for (area_region, T_region, ts_region, w, N_i, qv_region, qr_region, qs_region, is_updraft) in regions
            ρ_region = TD.air_density(thermo_params, ts_region)
            q_region = TD.PhasePartition(thermo_params, ts_region)

            # Note if T makes a very low excursion (or too high) , while S*G might mostly cancel and give a real result, you can get NaN from 0 * Inf or something, but I think this is just a legitimate model crash bc those temps are super extreme.
            qvsat_liq = TD.q_vap_saturation_generic(thermo_params, T_region, ρ_region, TD.Liquid()) # we don't have this stored for grid-mean and we can't calculate form en/up bc it's non-linear...
            dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(1), T_region, qvsat_liq) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...
            Γ_l = one(FT) + (L_v / c_p) * dqsl_dT  # Eqn C3
            δ = (qv_region - qvsat_liq) / Γ_l
            δ = min(δ, FT(0)) # only consider subsaturation for evaporation
            # S_qr_evap = -min(qr / Δt, -α_evp * CM1.evaporation_sublimation(microphys_params, rain_type, q_region, qr_region, ρ_region, T_region)) * precip_fraction
            S_qr_evap_here = limit_tendency(ptl, α_evp * CM1.evaporation_sublimation(microphys_params, rain_type, q_region, qr_region, ρ_region, T_region), min(qr_region, -δ), Δt) * precip_fraction # i guess rain only evaporates but never condenses?

            # S_qs_melt = -min(qs / Δt, α_melt * CM1.snow_melt(microphys_params, qs_region, ρ_region, T_region)) * precip_fraction
            S_qs_melt_here = limit_tendency(ptl, -α_melt * CM1.snow_melt(microphys_params, qs_region, ρ_region, T_region), qs_region, Δt) * precip_fraction

            # -------------------------- #
            # if edmf.moisture_model isa NonEquilibriumMoisture
            #     N_INP = get_INP_concentration(param_set, edmf.moisture_model.scheme, q_region, T_region, ρ_region, w)
            # else
            #     N_INP = get_N_i_Cooper_curve(T_region; clamp_N=true)
            # end

            # N_s_min = max(N_INP - N_i, FT(0)) # the fewest N,
            # if N_s_min > FT(0) && (qs_region > 1e-9) # we could add dry aerosol mass but we don't know N, # we will do this by varying lambda...
            #     _, λ_min = get_n0_lambda(param_set, snow_type, qs_region, ρ_region, N_s_min) # the largest droplets we should allow... [so lambda is smallest]
            # else
            #     λ_min = eps(FT)
            # end
            # λ_s = CM1.lambda(microphys_params, snow_type, qs_region, ρ_region)
            # λ_s = max(λ_s, λ_min) # so that we don't sublimate more than we have INP for
            λ_s = FT(NaN)
            # -------------------------- #

            tmp = α_dep_sub * my_evaporation_sublimation(microphys_params, snow_type, q_region, qs_region, ρ_region, T_region; _λ=λ_s) * precip_fraction

            # Note if T makes a very low excursion (or too high) , while S*G might mostly cancel and give a real result, you can get NaN from 0 * Inf or something, but I think this is just a legitimate model crash bc those temps are super extreme.
            qvsat_ice = TD.q_vap_saturation_generic(thermo_params, T_region, ρ_region, TD.Ice()) # we don't have this stored for grid-mean and we can't calculate form en/up bc it's non-linear...
            dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(0), T_region, qvsat_ice) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...
            Γ_i = one(FT) + (L_s / c_p) * dqsi_dT  # Eqn C3
            δi = (qv_region - qvsat_ice) / Γ_i

            dv = FT(0)
            # dv += Δt * max(aux_tc.dqvdt[k] * aux_tc.area[k], FT(0)) # hopefully this is ok to just add like this...
            if is_updraft
                # dv += Δt * max(( (aux_tc.massflux_tendency_qt[k] + aux_tc.diffusive_tendency_qt[k]) - max(aux_tc.massflux_tendency_ql[k] + aux_tc.diffusive_tendency_ql[k], FT(0)) - max(aux_tc.massflux_tendency_qi[k] + aux_tc.diffusive_tendency_qi[k], FT(0))), FT(0))
                # dv += Δt * max(aux_gm.qt_tendency_ls_vert_adv[k] - max(aux_gm.ql_tendency_ls_vert_adv[k], FT(0)) - max(aux_gm.qi_tendency_ls_vert_adv[k], FT(0)), FT(0))
            end
            
            if tmp > 0
                δi = max(δi, FT(0)) # only consider depositional supersaturation for sublimation/deposition
                # since we don't wanna do the noneq_moisture_sources() style thing, we'll settle for just adding dqvdt, and qt advection here. sed doesn't count.., lsadv does but is slow
                
                # S_qs_sub_dep_here = -limit_tendency(ptl, -tmp, qv_region+dv, Δt) # might be too restrictive at cloud top w/o the dv
                S_qs_sub_dep_here = -limit_tendency(ptl, -tmp, max(FT(0), min(δi+dv/2, qv_region+dv/2, δi*Δt/2)), Δt) # presumably if tmp > 0 then δi > 0 but can't be too careful, δi*Δt would be assuming you're in steady state...
            else
                δi = min(δi, FT(0)) # only consider subsaturational undersaturation for sublimation/deposition
                S_qs_sub_dep_here = limit_tendency(ptl, tmp, max(FT(0), min(qs_region, -δi)), Δt) # Hard to use dv here because gains in vapor are fighting against sublimation.
            end

            # accumulate over regions
            S_qr_evap += area_region * S_qr_evap_here
            S_qs_melt += area_region * S_qs_melt_here
            S_qs_sub_dep += area_region * S_qs_sub_dep_here
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




"""
    A copy just so we can override the default n0, _λ assumptions
"""
function my_evaporation_sublimation(
    prs::ACMP,
    snow_type::CMT.SnowType,
    q::TD.PhasePartition{FT},
    q_sno::FT,
    ρ::FT,
    T::FT;
    _λ::FT = FT(NaN),
) where {FT <: Real}
    evap_subl_rate = FT(0)
    if q_sno > FT(0)
        _ν_air::FT = CMP.ν_air(prs)
        _D_vapor::FT = CMP.D_vapor(prs)

        thermo_params = CMP.thermodynamics_params(prs)
        _S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        _G::FT = CM.Common.G_func(prs, T, TD.Ice())

        _n0::FT = n0(prs, q_sno, ρ, snow_type)
        _r0::FT = r0(prs, snow_type)
        _χv::FT = χv(prs, snow_type)
        _v0::FT = v0(prs, ρ, snow_type)
        _ve::FT = ve(prs, snow_type)
        _Δv::FT = Δv(prs, snow_type)

        _a_vent::FT = a_vent(prs, snow_type)
        _b_vent::FT = b_vent(prs, snow_type)

        _λ::FT = isnan(_λ) ? CM1.lambda(prs, snow_type, q_sno, ρ) : _λ

        evap_subl_rate =
            4 * FT(π) * _n0 / ρ * _S * _G / _λ^FT(2) * (
                _a_vent +
                _b_vent * (_ν_air / _D_vapor)^FT(1 / 3) /
                (_r0 * _λ)^((_ve + _Δv) / FT(2)) *
                (FT(2) * _v0 * _χv / _ν_air / _λ)^FT(1 / 2) *
                CM1.SF.gamma((_ve + _Δv + FT(5)) / FT(2))
            )
    end

    if !isfinite(evap_subl_rate)
        @error "Got non-finite evap_subl_rate $evap_subl_rate at T = $T; q_sno = $q_sno; ρ = $ρ; λ = $_λ; n0 = $(_n0); S = $(_S); G = $(_G); a_vent = $(_a_vent); b_vent = $(_b_vent); ν_air = $(_ν_air); D_vapor = $(_D_vapor); r0 = $(_r0); v0 = $(_v0); χv = $(_χv); ve = $(_ve); Δv = $(_Δv) "
    end

    # return evap_subl_rate
    return resolve_nan(evap_subl_rate, FT(0))
end