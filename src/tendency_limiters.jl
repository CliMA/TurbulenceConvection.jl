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
   # deprecated, see an old commit that had it
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