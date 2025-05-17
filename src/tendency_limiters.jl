"""
limit tendencies to improve model stability

These are deprecated bc I can't figure out a good way to keep them self-consistent (we already have filtering)

    * Haven't yet come up with self-consistent limiters
    > What we really need is to track all the edges between input and output variables we care about, and use an algorithm to scale down edges so that lead to domain violations
        - that's not easy though, bc scaling down one edge might lead to another edge being violated, so you need some kind of iterative or other clever algorithm, and you need to store all the tendencies which could greatly increase allocations... (though we're storing a lot of them now for output, but not for calibration)


    Note we still have the process limiters we introduced and adaptive timestepping... and the fallback to the old timestepping and limiting after that. If that doesn't work out, idk what to tell ya...
"""



function limit_gm_tendencies!(edmf::EDMFModel, grid::Grid, state::State, Δt::Real, include::Bool)
    # deprecated bc you can't really just change the grid mean without checking if the cause was the updrafts/env tendencies...
    return nothing
end

function limit_gm_tendencies_old!(edmf::EDMFModel, grid::Grid, state::State, param_set::APS, surf::SurfaceBase, Δt::Real)
    # Currently we have no control on gm tendencies. We should. Otherwise, you can get to unphysical states.
    # Still, this isn't foolproof. Some tendencies from noneq/precip etc do map onto the gm tendencies -- these need to be limited in tandem.
    # Still, other tendencies are from external forcings and are not shared by updraft and gm.
    # How to compose the two is a hard problem.

    # However, it should always be the case that prog_gm > prog_up/bulk, so if updraft tendencies are consistent with prog_up, then prog_gm can handle those same tendencies and any limiting is truly likely to have come from external grid mean additions.
    # Note prog_gm > prog_up/bulk is not a guarantee, but I added it as a stabilizing condition in filter_updraft_vars()
    prog_gm = center_prog_grid_mean(state)
    tendencies_gm = center_tendencies_grid_mean(state)

    # force individual gm tendencies to not deplete themselves
    # @. tendencies_gm.ρq_tot = max(tendencies_gm.ρq_tot, -prog_gm.ρq_tot / Δt)
    # @. tendencies_gm.ρθ_liq_ice = max(tendencies_gm.ρθ_liq_ice, -prog_gm.ρθ_liq_ice / Δt)
    # if edmf.moisture_model isa NonEquilibriumMoisture
    #     @. tendencies_gm.q_liq = max(tendencies_gm.q_liq, -prog_gm.q_liq / Δt) # should these happen here or when updraft tendencies are adjusted? we (for now) have no external q_liq forcing
    #     @. tendencies_gm.q_ice = max(tendencies_gm.q_ice, -prog_gm.q_ice / Δt) # should these happen here or when updraft tendencies are adjusted? we (for now) have no external q_ice forcing
    # end

    # dtl = edmf.tendency_limiters.default_tendency_limiter
    dtl = get_tendency_limiter(edmf.tendency_limiters, Val(:default), use_fallback_tendency_limiters)
    @. tendencies_gm.ρq_tot = limit_tendency(dtl, tendencies_gm.ρq_tot, prog_gm.ρq_tot, Δt)
    @. tendencies_gm.ρθ_liq_ice = limit_tendency(dtl, tendencies_gm.ρθ_liq_ice, prog_gm.ρθ_liq_ice, Δt)
    if edmf.moisture_model isa NonEquilibriumMoisture
        @. tendencies_gm.q_liq = limit_tendency(dtl, tendencies_gm.q_liq, prog_gm.q_liq, Δt)
        @. tendencies_gm.q_ice = limit_tendency(dtl, tendencies_gm.q_ice, prog_gm.q_ice, Δt)
    end

    return nothing
end

""" 
Deprecate limiting up tendencies for now -- there is no world in which we can do this successfully w/o a more detailed balance of all contributing tendencies sources and sinks
    e.g. if we are using up too much q_liq, we could scale down all the tendencies that contribute to q_liq depletion, but we would need to follow that graph backwards...
    Also an optimal order of operations etc is not an easy problem... (I'm not even sure you can guarantee a stable solution by going sequentially)

    Maybe there is some graph algorithm for this...

    For now -- either use tendency_limited adaptive_dt with NO tendencies limited or use an implicit solver...
        or
    limit all tendencies the same way and hope your timesteps are small enough for some accuracy... You're in danger of tendencies saturating at different times in the timestep, but I guess that's also true w/ filtering...

    The general problem is that we can make many more guarantees about shrinking than we can about growth... but it's hard to do that w/o tracking everything carefully...
"""
function limit_up_tendencies!(edmf::EDMFModel, grid::Grid, state::State, param_set::APS, surf::SurfaceBase, Δt::Real, is_before_gm_calc::Bool)
    return nothing
end

"""

This function limits tendencies after gm tendencies are calculated to force them to be stable
Do we need to do anything around the surface?

Technically this can be bad because the tendencies are tied together
... if say the change in q requires limiting but the change in θ doesn't, we get leaks in energy/moisture conservation etc. The env should take up some of this but...


 NOTE: We can't determine which adjustments affect the grid-mean and which do not... So this isn't an ideal solution.
    It could be better to limit tendencies at the point of creation, though that introduces biases in composability. Take short steps and pray or use an implicit solver.
    This sol'n however is perhaps closest to the former sol'n by Anna, Charlie et al of stepping and filtering. This is a bit more robust than that but still not perfect.

    Probably the only true sol'n is to limit the sum of the updraft tendencies before the gm tendencies are calculated.
    So at the end of compute_up_tendencies!()

    What do we do about en tendencies? those just go in to the gm() so limiting the gm tendencies should be enough? idk...

    is_before_gm_calc::Bool -- if true, we're limiting the updraft tendencies before the gm tendencies are calculated, if false, we're limiting the gm tendencies after they're calculated
        - ideally it's true, though this is not robust to the gm nudging. However, that should only matter for q_tot and θ_liq_ice, which should be very unlikely to run into an issue. The other gm tendencies come from the underlying up/env temdemcies.
            - if we do limit the gm tendencies, that *should* handle env, but if we dont, env tendencies are never limited. Is that bad? I feel we should not limit up but not env...
        - if it's false, when you update the updrafts, you have to update the gm tendencies again, which becomes an iterative process...


    This is also not a magic pill for subprocesses. E.g. if you limit the updraft ρaq_liq tendency, but don't edit the rain formation tendency, you could have a leak.
        The right way would be to use the precipitation_tendency_limiter or use a dt limited timestep with NO limiting anywhere or implicit solver again with NO limiting anywhere.

"""
function limit_up_tendencies_old!(edmf::EDMFModel, grid::Grid, state::State, param_set::APS, surf::SurfaceBase, Δt::Real, is_before_gm_calc::Bool)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    tendencies_gm = center_tendencies_grid_mean(state)

    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    aux_bulk = center_aux_bulk(state)
    aux_bulk_f = face_aux_bulk(state) # same as using aux_tc_f = face_aux_turbconv(state) and then using aux_tc_f.bulk.w for example

    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ

    N_up = n_updrafts(edmf)
    FT = float_type(state)


    # deprecate, bc of cutoffs, min areas, etc, this can go south and detach from what the tendencies show
    area_bulk = aux_bulk.area
    θ_liq_ice_bulk = aux_bulk.θ_liq_ice
    q_tot_bulk = aux_bulk.q_tot
    if edmf.moisture_model isa NonEquilibriumMoisture
        q_liq_bulk = aux_bulk.q_liq
        q_ice_bulk = aux_bulk.q_ice
    end
    

    prog_bulk = center_prog_bulk(state)
    prog_bulk_f = face_prog_bulk(state)

    ρarea_bulk = prog_bulk.ρarea
    ρaθ_liq_ice_bulk = prog_bulk.ρaθ_liq_ice
    ρaq_tot_bulk = prog_bulk.ρaq_tot
    if edmf.moisture_model isa NonEquilibriumMoisture
        ρaq_liq_bulk = prog_bulk.ρaq_liq
        ρaq_ice_bulk = prog_bulk.ρaq_ice
    end
    ρaw_bulk = prog_bulk_f.ρaw

    area_en = aux_en.area
    q_tot_en = aux_en.q_tot
    θ_liq_ice_en = aux_en.θ_liq_ice
    if edmf.moisture_model isa NonEquilibriumMoisture
        q_liq_en = aux_en.q_liq
        q_ice_en = aux_en.q_ice
    end


    # calculate the bulk updraft tendency
    # check if it violates the limit
    # if it does, scale down all the updrafts equally (not ideal but best for rn)

    # ==== should be stable up to a gm tendency that can't be absorbed by the env... (is that possible?) ==== #

    tendencies_bulk = center_tendencies_bulk(state)
    tendencies_bulk_f = face_tendencies_bulk(state)

    tends_ρarea_bulk = tendencies_bulk.ρarea
    tends_ρaθ_liq_ice_bulk = tendencies_bulk.ρaθ_liq_ice
    tends_ρaq_tot_bulk = tendencies_bulk.ρaq_tot
    @. tends_ρarea_bulk = FT(0) # we're gonna recalculate this in the limiter
    @. tends_ρaθ_liq_ice_bulk = FT(0) # we're gonna recalculate this in the limiter
    @. tends_ρaq_tot_bulk = FT(0) # we're gonna recalculate this in the limiter

    if edmf.moisture_model isa NonEquilibriumMoisture
        tends_ρaq_liq_bulk = tendencies_bulk.ρaq_liq
        tends_ρaq_ice_bulk = tendencies_bulk.ρaq_ice
        @. tends_ρaq_liq_bulk = FT(0) # we're gonna recalculate this in the limiter
        @. tends_ρaq_ice_bulk = FT(0) # we're gonna recalculate this in the limiter
    end
    tends_ρaw_bulk = tendencies_bulk_f.ρaw
    @. tends_ρaw_bulk = FT(0) # we're gonna recalculate this in the limiter

    # force individual updraft tendencies to not deplete themselves and then calculate bulk tendencies
    # these were zeroed out in update_aux so should be fine...

    # uetl = edmf.tendency_limiters.up_en_tendency_limiter
    uetl = get_tendency_limiter(edmf.tendency_limiters, Val(:default), use_fallback_tendency_limiters)
    @inbounds for i in 1:N_up 
        # @. tendencies_up[i].ρarea = max(tendencies_up[i].ρarea, -ρ_c * aux_up[i].area / Δt)
        # @. tendencies_up[i].ρaθ_liq_ice = max(tendencies_up[i].ρaθ_liq_ice, -ρ_c * aux_up[i].area * aux_up[i].θ_liq_ice / Δt)
        # @. tendencies_up[i].ρaq_tot = max(tendencies_up[i].ρaq_tot, -ρ_c * aux_up[i].area * aux_up[i].q_tot / Δt)
        @. tendencies_up[i].ρarea = limit_tendency(uetl, tendencies_up[i].ρarea, ρ_c * aux_up[i].area, Δt)
        @. tendencies_up[i].ρaθ_liq_ice = limit_tendency(uetl, tendencies_up[i].ρaθ_liq_ice, ρ_c * aux_up[i].area * aux_up[i].θ_liq_ice, Δt)
        @. tendencies_up[i].ρaq_tot = limit_tendency(uetl, tendencies_up[i].ρaq_tot, ρ_c * aux_up[i].area * aux_up[i].q_tot, Δt)

        # @. tendencies_up_f[i].ρaw = max(tendencies_up_f[i].ρaw, -ρ_f * ᶠinterp_a(aux_up[i].area) * aux_up_f[i].w / Δt)
        @. tendencies_up_f[i].ρaw = limit_tendency(uetl, tendencies_up_f[i].ρaw, ρ_f * ᶠinterp_a(aux_up[i].area) * aux_up_f[i].w, Δt)

        @. tends_ρarea_bulk += tendencies_up[i].ρarea
        @. tends_ρaθ_liq_ice_bulk += tendencies_up[i].ρaθ_liq_ice
        @. tends_ρaq_tot_bulk += tendencies_up[i].ρaq_tot

        if edmf.moisture_model isa NonEquilibriumMoisture
            # @. tendencies_up[i].ρaq_liq = max(tendencies_up[i].ρaq_liq, -ρ_c * aux_up[i].area * aux_up[i].q_liq / Δt)
            # @. tendencies_up[i].ρaq_ice = max(tendencies_up[i].ρaq_ice, -ρ_c * aux_up[i].area * aux_up[i].q_ice / Δt)
            @. tendencies_up[i].ρaq_liq = limit_tendency(uetl, tendencies_up[i].ρaq_liq, ρ_c * aux_up[i].area * aux_up[i].q_liq, Δt)
            @. tendencies_up[i].ρaq_ice = limit_tendency(uetl, tendencies_up[i].ρaq_ice, ρ_c * aux_up[i].area * aux_up[i].q_ice, Δt)

            @. tends_ρaq_liq_bulk += tendencies_up[i].ρaq_liq
            @. tends_ρaq_ice_bulk += tendencies_up[i].ρaq_ice
        end
        @. tends_ρaw_bulk += tendencies_up_f[i].ρaw
    end


    # calculate scaling to bring into range... you can't lose more than you have, and you can't gain more than the environment can give plus the grid mean tendency. use safe_clamp bc gm_tendency could change make the max value less than the min value and clamp won't work correctly.
    # I guess this could fail if say tends_ρaq_tot_bulk is e.g. 0 but both [[ -ρaq_tot_bulk / Δt ]] and [[ (ρ_c * area_en * q_tot_en) / Δt  + tendencies_gm.ρq_tot ]] are both negative... shouldn't really happen though..? even if it does, we sort of are just enforcing the gm tendency...
    # The only problem is you still can't lose more than you have so you'd need to ensure that e.g. [[ (ρ_c * area_en * q_tot_en) / Δt + tendencies_gm.ρq_tot ]] is still greater than [[ -ρaq_tot_bulk / Δt ]]
    # Really though what that boils down to is ensuring tendencies_gm is not more negative than the updraft plus the environment can provide...
    # maybe additive is more robust for when the original tendency is 0 but gm tendency necessitates otherwise...
    # If the original tendency is 0, we this will give inf, which we should also just rescale to 1 or 0


    # Fix gm tendencies to not be more negative than the existing values can provide
    # - moved to limit_gm_tendencies!()

    #=
    Should gm tendencies even be allowed to affect the updraft though? 
    Should it be 'far removed' enough? also really you only want the LS tendencies to matter, you don't want environment things for example to feed back on the updraft limiting...
 
    If this stability scaling works, we should add permanent variables to the state to store the scaling factors so we can stop allocating
    =#

    tendencies_bulk_adjustments = center_tendencies_bulk_adjustments(state)
    tendencies_bulk_adjustments_f = face_tendencies_bulk_adjustments(state)

    tends_adjustment_ρarea = tendencies_bulk_adjustments.ρarea
    tends_adjustment_ρaθ_liq_ice = tendencies_bulk_adjustments.ρaθ_liq_ice
    tends_adjustment_ρaq_tot = tendencies_bulk_adjustments.ρaq_tot
    tends_adjustment_ρaw = tendencies_bulk_adjustments_f.ρaw

    if edmf.moisture_model isa NonEquilibriumMoisture
        tends_adjustment_ρaq_liq = tendencies_bulk_adjustments.ρaq_liq
        tends_adjustment_ρaq_ice = tendencies_bulk_adjustments.ρaq_ice
    end


    a_min = edmf.minimum_area
    a_max = edmf.max_area


    if is_before_gm_calc # we can't use the gm nudging values bc they haven't been calculated yet. To be fair though, that was the default state before any of my edits.
        # en and gm tendencies wille be calc'd later...
        # This is BAD bc we have no upper bound on growth (the env does not provide this bound bc it can't know e.g. sedimentation and advection that are nonlocal)
        # - maybe removing sedimentation and advection it should work? then limit growth from cond_evap, etc... 
        # - but then incoming tendencies from  

        ∑tendencies_ρaq_tot_local_bulk = @. tends_ρaq_tot_bulk - ρ_c * area_bulk * aux_bulk.qt_tendency_sedimentation - ρ_c * area_bulk * aux_bulk.qt_tendency_precip_formation
        ∑tendencies_ρaθ_liq_ice_local_bulk = @. tends_ρaθ_liq_ice_bulk - ρ_c * area_bulk * aux_bulk.θ_liq_ice_tendency_sedimentation - ρ_c * area_bulk * aux_bulk.θ_liq_ice_tendency_precip_formation

        # Additive scaling would be more robust right, no division problems, and gm tendencies that necessitate a change in the updraft tendency can be handled more easily
        # bc processes can change the env, we'd like to but can't constrain the change in env to be much more than remaining above 0...
        # @. tends_adjustment_ρarea = safe_clamp(tends_ρarea_bulk, - ρ_c * area_bulk / Δt, (ρ_c * a_max - ρarea_bulk) / Δt) - tends_ρarea_bulk # updrfat can only consume up to a_max, not all of env.
        # @. tends_adjustment_ρaθ_liq_ice = safe_clamp(tends_ρaθ_liq_ice_bulk, -ρ_c * area_bulk * θ_liq_ice_bulk / Δt, (ρ_c * area_en * θ_liq_ice_en) / Δt + tendencies_gm.ρθ_liq_ice) - tends_ρaθ_liq_ice_bulk
        # @. tends_adjustment_ρaq_tot = safe_clamp(tends_ρaq_tot_bulk, -ρ_c * area_bulk * q_tot_bulk / Δt, (ρ_c * area_en * q_tot_en) / Δt + tendencies_gm.ρq_tot) - tends_ρaq_tot_bulk
        # @. tends_adjustment_ρaw = safe_clamp(tends_ρaw_bulk, -ρ_f * ᶠinterp_a(area_bulk) * aux_bulk_f.w / Δt, Inf) - tends_ρaw_bulk # no upper bound bc of buoyancy production and nh_pressure

        # No real upper bound bc we don't know the gm tendency...
        @. tends_adjustment_ρarea = limit_tendency(uetl, tends_ρarea_bulk, ρ_c * area_bulk, ρ_c * a_max - ρarea_bulk, Δt) - tends_ρarea_bulk
        @. tends_adjustment_ρaθ_liq_ice = limit_tendency(uetl, tends_ρaθ_liq_ice_bulk, Δt) - tends_ρaθ_liq_ice_bulk
        @. tends_adjustment_ρaq_tot = limit_tendency(uetl, tends_ρaq_tot_bulk, Δt) - tends_ρaq_tot_bulk
        @. tends_adjustment_ρaw = limit_tendency(uetl, tends_ρaw_bulk, ρ_f * ᶠinterp_a(area_bulk) * aux_bulk_f.w, Δt) - tends_ρaw_bulk

        if edmf.moisture_model isa NonEquilibriumMoisture
            # @. tends_adjustment_ρaq_liq = safe_clamp(tends_ρaq_liq_bulk, -ρ_c * area_bulk * q_liq_bulk / Δt, (ρ_c * area_en * q_liq_en) / Δt + ρ_c * tendencies_gm.q_liq) - tends_ρaq_liq_bulk
            # @. tends_adjustment_ρaq_ice = safe_clamp(tends_ρaq_ice_bulk, -ρ_c * area_bulk * q_ice_bulk / Δt, (ρ_c * area_en * q_ice_en) / Δt + ρ_c * tendencies_gm.q_ice) - tends_ρaq_ice_bulk
            @. tends_adjustment_ρaq_liq = limit_tendency(uetl, tends_ρaq_liq_bulk, Δt) - tends_ρaq_liq_bulk
            @. tends_adjustment_ρaq_ice = limit_tendency(uetl, tends_ρaq_ice_bulk, Δt) - tends_ρaq_ice_bulk
        end

    else # after gm is calculated, the nudge values are correct, but the gm tendencies would need adjustment...

        # we have : up_tend + env_tend + nudging = gm_tend = up_tend_no_entr_detr + entr_detr + nudging  = gm_tend
        # we have no way to limit the growth of anything here but entr/detr... [and area in general]
        # we can limit the decrease of anything here but entr/detr... but we would also need to figure out mappings between those tendencies and their sources (e.g. does q_liq  tend come from evaporation [q_vap] or from sedimentation [q_liq, q_tot] or precip formation [q_liq, q_tot] etc.)
        # otherwise we'll have leaks in the system...

        # we can limit updraft - entr/detr to be within the bounds of the env [no we can't, they have their own sources... f]
        tends_gm_ρθ_liq_ice_nudge_only = @. tendencies_gm.ρθ_liq_ice - tends_ρaθ_liq_ice_bulk - tends_ρaθ_liq_ice_en
        tends_gm_ρq_tot_nudge_only = @. tendencies_gm.ρq_tot - tends_ρaq_tot_bulk - tends_ρaq_tot_en
        if edmf.moisture_model isa NonEquilibriumMoisture
            tends_gm_ρq_liq_nudge_only = @. ρ_c * tendencies_gm.q_liq - tends_ρaq_liq_bulk - tends_ρaq_liq_en
            tends_gm_ρq_ice_nudge_only = @. ρ_c * tendencies_gm.q_ice - tends_ρaq_ice_bulk - tends_ρaq_ice_en
        end

        # Additive scaling would be more robust right, no division problems, and gm tendencies that necessitate a change in the updraft tendency can be handled more easily
        # bc processes can change the env, we'd like to but can't constrain the change in env to be much more than remaining above 0...
        # @. tends_adjustment_ρarea = safe_clamp(tends_ρarea_bulk, - ρ_c * area_bulk / Δt, (ρ_c * a_max - ρarea_bulk) / Δt) - tends_ρarea_bulk # updrfat can only consume up to a_max, not all of env.
        # @. tends_adjustment_ρaθ_liq_ice = safe_clamp(tends_ρaθ_liq_ice_bulk, -ρ_c * area_bulk * θ_liq_ice_bulk / Δt, (ρ_c * area_en * θ_liq_ice_en) / Δt + tendencies_gm.ρθ_liq_ice) - tends_ρaθ_liq_ice_bulk
        # @. tends_adjustment_ρaq_tot = safe_clamp(tends_ρaq_tot_bulk, -ρ_c * area_bulk * q_tot_bulk / Δt, (ρ_c * area_en * q_tot_en) / Δt + tendencies_gm.ρq_tot) - tends_ρaq_tot_bulk
        # @. tends_adjustment_ρaw = safe_clamp(tends_ρaw_bulk, -ρ_f * ᶠinterp_a(area_bulk) * aux_bulk_f.w / Δt, Inf) - tends_ρaw_bulk # no upper bound bc of buoyancy production and nh_pressure

        # No real upper bound bc we don't know the gm tendency...
        @. tends_adjustment_ρarea = limit_tendency(uetl, tends_ρarea_bulk, ρ_c * area_bulk, ρ_c * a_max - ρarea_bulk, Δt) - tends_ρarea_bulk
        @. tends_adjustment_ρaθ_liq_ice = limit_tendency(uetl, tends_ρaθ_liq_ice_bulk, Δt) - tends_ρaθ_liq_ice_bulk
        @. tends_adjustment_ρaq_tot = limit_tendency(uetl, tends_ρaq_tot_bulk, Δt) - tends_ρaq_tot_bulk
        @. tends_adjustment_ρaw = limit_tendency(uetl, tends_ρaw_bulk, ρ_f * ᶠinterp_a(area_bulk) * aux_bulk_f.w, Δt) - tends_ρaw_bulk

        if edmf.moisture_model isa NonEquilibriumMoisture
            # @. tends_adjustment_ρaq_liq = safe_clamp(tends_ρaq_liq_bulk, -ρ_c * area_bulk * q_liq_bulk / Δt, (ρ_c * area_en * q_liq_en) / Δt + ρ_c * tendencies_gm.q_liq) - tends_ρaq_liq_bulk
            # @. tends_adjustment_ρaq_ice = safe_clamp(tends_ρaq_ice_bulk, -ρ_c * area_bulk * q_ice_bulk / Δt, (ρ_c * area_en * q_ice_en) / Δt + ρ_c * tendencies_gm.q_ice) - tends_ρaq_ice_bulk
            @. tends_adjustment_ρaq_liq = limit_tendency(uetl, tends_ρaq_liq_bulk, Δt) - tends_ρaq_liq_bulk
            @. tends_adjustment_ρaq_ice = limit_tendency(uetl, tends_ρaq_ice_bulk, Δt) - tends_ρaq_ice_bulk
        end

        # recalc gm tendencies... (this elides the problem of figuring out if )
        @. tendencies_gm.ρq_tot = tends_gm_ρq_tot_nudge_only + tends_ρaq_tot_bulk + tends_adjustment_ρaθ_liq_ice + tends_ρaq_tot_en 
        @. tendencies_gm.ρθ_liq_ice = tends_gm_ρθ_liq_ice_nudge_only + tends_ρaθ_liq_ice_bulk + tends_adjustment_ρaθ_liq_ice + tends_ρaθ_liq_ice_en
        if edmf.moisture_model isa NonEquilibriumMoisture
            @. tendencies_gm.q_liq = (tends_gm_ρq_liq_nudge_only + tends_ρaq_liq_bulk + tends_adjustment_ρaq_liq + tends_ρaq_liq_en) / ρ_c
            @. tendencies_gm.q_ice = (tends_gm_ρq_ice_nudge_only + tends_ρaq_ice_bulk + tends_adjustment_ρaq_ice + tends_ρaq_ice_en) / ρ_c
        end

    end


    # apply scaling to all updrafts (additive version, let the contribution to the tendency be relative to the updraft's contribution to the bulk total -- should ensure stability? It's not perfect for >1 updraft but...

    # For more than 1 updraft, it's unclear which updraft should get the adjustment... (e.g. if some are negative and can't go more negative then only the ones postive can bear the load...)
    # we would need some logic to disentagle this... essentially updraft shrinking is already limited so we'd jus need to constrain growing updrafts together perhaps...
    if N_up == 1
        tends_ρarea = tendencies_up[1].ρarea
        tends_ρaθ_liq_ice = tendencies_up[1].ρaθ_liq_ice
        tends_ρaq_tot = tendencies_up[1].ρaq_tot
        if edmf.moisture_model isa NonEquilibriumMoisture
            tends_ρaq_liq = tendencies_up[1].ρaq_liq
            tends_ρaq_ice = tendencies_up[1].ρaq_ice
        end
        
        @inbounds for k in real_center_indices(grid)
            tends_ρarea[k] += tends_adjustment_ρarea[k] 
            tends_ρaθ_liq_ice[k] += tends_adjustment_ρaθ_liq_ice[k]
            tends_ρaq_tot[k] += tends_adjustment_ρaq_tot[k] 
            # perhaps there would be a conservation term here but idk what it would be... changes in q_tot could come from relaxation but that is not conserved. vapor doesn't matter for liq ice pottemp so maybe it's fine?
            if edmf.moisture_model isa NonEquilibriumMoisture
                thermo_params = TCP.thermodynamics_params(param_set)
                ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ_c[k], aux_bulk.T[k], TD.PhasePartition(aux_bulk.q_tot[k], aux_bulk.q_liq[k], aux_bulk.q_ice[k]))
                L_v0 = TCP.LH_v0(param_set)
                L_s0 = TCP.LH_s0(param_set)
                Π_m = TD.exner(thermo_params, ts)
                c_pm = TD.cp_m(thermo_params, ts)
                tends_ρaq_liq[k] += tends_adjustment_ρaq_liq[k]
                tends_ρaq_ice[k] += tends_adjustment_ρaq_ice[k] 
                tends_ρaq_tot[k] += tends_adjustment_ρaq_liq[k] + tends_adjustment_ρaq_ice[k] # conservation (though not really since relaxation is not conserved)
                tends_ρaθ_liq_ice[k] += tends_adjustment_ρaq_liq[k] / Π_m / c_pm * L_v0  + tends_adjustment_ρaq_ice[k] / Π_m / c_pm * L_s0 # conservation (though not really since relaxation is not conserved)
            end
        end
        # limit ρaw
        tends_ρaw = tendencies_up_f[1].ρaw
        @inbounds for k in real_face_indices(grid)
            tends_ρaw[k] += tends_adjustment_ρaw[k]
        end

        # we really may need to limit the gm tendencies again... the `adjustments` we made could be quite large, especially if the original tendencies were not at all limited...
    else
        error("Updraft tendency limiting for more than 1 updraft not implemented yet because specialized logic is needed to determine which updrafts should be affected by the limiting adjustments.")
        # need to enforce constraints only on growing updrafts bc shrinking updrafts are already limited and shouldn't cause any limiter issues and go to env just fine... 
        # a framework would be needed to determine which updrafts are growing and which are shrinking...
    end

   
    return nothing
end

function limit_en_tendencies_old!(edmf::EDMFModel, grid::Grid, state::State, param_set::APS, surf::SurfaceBase, Δt::Real)
    return nothing
end

"""
There is no option to do this after limiting the gm tendencies here since theyre essentially one and the same...
Deprecate this rn for similar reasons to the updraft limiting...
"""
function limit_en_tendencies_old!(edmf::EDMFModel, grid::Grid, state::State, param_set::APS, surf::SurfaceBase, Δt::Real)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    aux_bulk = center_aux_bulk(state)
    aux_bulk_f = face_aux_bulk(state) # same as using aux_tc_f = face_aux_turbconv(state) and then using aux_tc_f.bulk.w for example

    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ

    N_up = n_updrafts(edmf)
    FT = float_type(state)

    # en has no area change mechanism, only q_tot, q_liq, q_ice, θ_liq_ice

    prog_en = center_prog_environment(state)
    prog_en_f = face_prog_environment(state)

    ρaθ_liq_ice_en = prog_en.ρaθ_liq_ice
    ρaq_tot_en = prog_en.ρaq_tot
    if edmf.moisture_model isa NonEquilibriumMoisture
        ρaq_liq_en = prog_en.ρaq_liq
        ρaq_ice_en = prog_en.ρaq_ice
    end

    ρaθ_liq_ice_en = @. aux_en.θ_liq_ice * prog_en.ρarea
    ρaq_tot_en = @. aux_en.q_tot * prog_en.ρarea
    if edmf.moisture_model isa NonEquilibriumMoisture
        ρaq_liq_en = @. aux_en.q_liq * prog_en.ρarea
        ρaq_ice_en = @. aux_en.q_ice * prog_en.ρarea
    end

    # uetl = edmf.tendency_limiters.up_en_tendency_limiter
    uetl = get_tendency_limiter(edmf.tendency_limiters, Val(:default), use_fallback_tendency_limiters)

    # We don't have a sum of env tendencies... there's _tendency_noneq, _tendency_precip_Formation, _tendency_sedimentation
    # bc env is diagnosed, there's no need to prevent consuming the entire updraft or anything, any growth here is already gonna be reflected in the gm... except entr/detr whose limiting is done in the updraft tendencies
    # for preventing excessive growth in general, use the individual tendency limiters...
    ∑en_tendencies_q_tot = @. aux_en.qt_tendency_precip_formation + aux_en.qt_tendency_sedimentation
    ∑en_tendencies_θ_liq_ice = @.  aux_en.θ_liq_ice_tendency_precip_formation + aux_en.θ_liq_ice_tendency_sedimentation

    limit_tendency(uetl, ∑en_tendencies_q_tot, ρ_c * aux_en.area * aux_en.q_tot, Δt)
    limit_tendency(uetl, ∑en_tendencies_θ_liq_ice, ρ_c * aux_en.area * aux_en.θ_liq_ice, Δt)
    if edmf.moisture_model isa NonEquilibriumMoisture
        ∑en_tendencies_q_liq = @. aux_en.ql_tendency_noneq + aux_en.ql_tendency_precip_formation + aux_en.ql_tendency_sedimentation
        ∑en_tendencies_q_ice = @. aux_en.qi_tendency_noneq + aux_en.qi_tendency_precip_formation + aux_en.qi_tendency_sedimentation
        limit_tendency(uetl, ∑en_tendencies_q_liq, ρ_c * aux_en.area * aux_en.q_liq, Δt)
        limit_tendency(uetl, ∑en_tendencies_q_ice, ρ_c * aux_en.area * aux_en.q_ice, Δt)
    end

end