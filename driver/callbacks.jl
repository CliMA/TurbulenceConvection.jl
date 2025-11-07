function condition_io(u, t, integrator)
    UnPack.@unpack TS, Stats = integrator.p
    TS.dt_io += TS.dt
    io_flag = false
    # if TS.dt_io > Stats[CC.Fields.ColumnIndex((1, 1), 1)].frequency
    if TS.dt_io ≥ Stats[CC.Fields.ColumnIndex((1, 1), 1)].frequency # allow exact dt outputs
        TS.dt_io = 0
        io_flag = true
    end
    return io_flag || t ≈ 0 || t ≈ TS.t_max
end

condition_every_iter(u, t, integrator) = true

function affect_io!(integrator)
    integrator.p.skip_io && return nothing

    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, calibrate_io, precip_model, aux, io_nt, diagnostics, surf_params, param_set, Stats) = integrator.p
        t = integrator.t
        prog = integrator.u
        stats = Stats[colidx]
        # TODO: remove `vars` hack that avoids
        # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
        # opening/closing files every step should be okay. #removeVarsHack
        # TurbulenceConvection.io(sim) # #removeVarsHack
        write_simulation_time(stats, t) # #removeVarsHack

        state = TC.column_prog_aux(prog, aux, colidx, calibrate_io)
        grid = TC.Grid(state)
        diag_col = TC.column_diagnostics(diagnostics, colidx)

        # TODO: is this the best location to call diagnostics?
        compute_diagnostics!(edmf, param_set, state, diag_col, stats, surf_params, t, calibrate_io)

        cent = TC.Cent(1)
        diag_svpc = svpc_diagnostics_grid_mean(diag_col)
        diag_tc_svpc = svpc_diagnostics_turbconv(diag_col)
        write_ts(stats, "lwp_mean", diag_svpc.lwp_mean[cent])
        write_ts(stats, "iwp_mean", diag_svpc.iwp_mean[cent])
        write_ts(stats, "rwp_mean", diag_svpc.rwp_mean[cent])
        write_ts(stats, "swp_mean", diag_svpc.swp_mean[cent])
        # write_ts(stats, "ipwp_mean", diag_svpc.swp_mean[cent]) # for all ice precipitation comparison w/ other models

        if !calibrate_io
            write_ts(stats, "integ_total_flux_qt", diag_svpc.integ_total_flux_qt[cent])
            write_ts(stats, "integ_total_flux_s", diag_svpc.integ_total_flux_s[cent])
            write_ts(stats, "updraft_cloud_cover", diag_tc_svpc.updraft_cloud_cover[cent])
            write_ts(stats, "updraft_cloud_base", diag_tc_svpc.updraft_cloud_base[cent])
            write_ts(stats, "updraft_cloud_top", diag_tc_svpc.updraft_cloud_top[cent])
            write_ts(stats, "env_cloud_cover", diag_tc_svpc.env_cloud_cover[cent])
            write_ts(stats, "env_cloud_base", diag_tc_svpc.env_cloud_base[cent])
            write_ts(stats, "env_cloud_top", diag_tc_svpc.env_cloud_top[cent])
            write_ts(stats, "env_lwp", diag_tc_svpc.env_lwp[cent])
            write_ts(stats, "env_iwp", diag_tc_svpc.env_iwp[cent])
            write_ts(stats, "Hd", diag_tc_svpc.Hd[cent])
            write_ts(stats, "updraft_lwp", diag_tc_svpc.updraft_lwp[cent])
            write_ts(stats, "updraft_iwp", diag_tc_svpc.updraft_iwp[cent])

            write_ts(stats, "cutoff_precipitation_rate", diag_svpc.cutoff_precipitation_rate[cent])
            write_ts(stats, "cloud_cover_mean", diag_svpc.cloud_cover_mean[cent])
            write_ts(stats, "cloud_base_mean", diag_svpc.cloud_base_mean[cent])
            write_ts(stats, "cloud_top_mean", diag_svpc.cloud_top_mean[cent])
        end

        io(io_nt.aux, stats, state)
        io(io_nt.diagnostics, stats, diag_col)

        surf = get_surface(surf_params, state, t, param_set) # calls Tsurface, shf, lhf, ustar, wstar automatically
        io(surf, surf_params, stats, t)
        nothing
    end

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end

function affect_filter!(integrator) #  affect_filter!() is also called in the driver directly... is that not redundant? should we also be calling the limiter again in callbacks?
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, param_set, aux, case, surf_params) = integrator.p
        t = integrator.t
        prog = integrator.u

        cfl_limit = integrator.p.TS.cfl_limit
        TS = integrator.p.TS
        Δt = TS.dt_limit_tendencies_factor * (TS.limit_tendencies_by_dt_min ? TS.dt_min : TS.dt)

        state = TC.column_prog_aux(prog, aux, colidx, integrator.p.calibrate_io)
        surf = get_surface(surf_params, state, t, param_set)
        TC.affect_filter!(edmf, state, param_set, surf, cfl_limit, Δt)
        nothing
    end

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end



"""
"""
function get_tendency_Δt(TS)    
    return TS.dt_limit_tendencies
end



"""
    When you are calculating tendencies that are functions of the timestep, you can call this function to reset to the longest possible timestep before recalculating the tendencies.

    NOTE: I'm not sure if this is ideal -- particularly, if the tendencies are scaled to stretch over a timestep, this could cause us to always err towards that long timestep. - however, you would never be able to violate dtmax even if the tendencies would allow it
    # resetting to dt_min would still allow for longer timesteps in the case in which no tendency would saturate, without editing the tendencies towards that longest timestep. - however, you would need to violate dt_min more often. 

    maybe you could err towards the midpoint of dt_min and dt_max or something but i think dt_min is hte cleanest given our N_dt_max_edmf_violate_dt_min_remaining capability...  revisit if problematic [or make it a parameter e.g. dt_guess = dt_min * 2 or something]

    This is potentially not that bad -- before dt_max was just calculated from CFL but note that microphysical tendencies were limited by Δt which probably wasn't ideal...) so ideally dt_min and dt_max would not vary much...

    maybe it should also be max(dt_min, TS.dt)
"""
function reset_dt!(integrator)
    UnPack.@unpack edmf, TS, dt_min = integrator.p
    proposed_dt = SciMLBase.get_proposed_dt(integrator)
    @debug "Resetting dt: TS.dt = $(TS.dt) | resetting to TS.dt = dt_min = $(dt_min)"
    TS.dt = dt_min
    
    SciMLBase.set_proposed_dt!(integrator, TS.dt)
    ODE.u_modified!(integrator, false)
end



function adaptive_dt!(integrator; depth::Int=0)
    UnPack.@unpack edmf, TS, dt_min = integrator.p
    t = integrator.t

    dt_old = TS.dt
    N_dt_max_edmf_violate_dt_min_remaining_old = TS.N_dt_max_edmf_violate_dt_min_remaining
    if TS.adapt_dt # calculate adaptive timestep
        if (TS.dt_max_edmf < dt_min) && (TS.N_dt_max_edmf_violate_dt_min_remaining > 0)
            TS.N_dt_max_edmf_violate_dt_min_remaining -= 1
            @debug " t = $t | dt_max_edmf = $(TS.dt_max_edmf) is less than dt_min = $(dt_min) (probably due to tendency or CFL constraint). Will allow up to $(TS.N_dt_max_edmf_violate_dt_min_remaining) more violations."
            TS.dt = TS.dt_max_edmf 
        else   
            TS.dt = min(TS.dt_max, max(TS.dt_max_edmf, dt_min))
        end
    end
 



    if (TS.spinup_half_t_max > 0) && (t < 2 * TS.spinup_half_t_max)
        if !TS.allow_spinup_adapt_dt # this will be false if adapt_dt is false, otherwise we could have an issue (this blocks here prevents merging of the if condition blocks)
            TS.dt = TS.dt_min # reset dt to dt_min if adapt_dt is false so we don't repeatedly adjust it, if allow_spinup_adapt_dt is true, adapt_dt is true and adapt_dt gets set before
            # should cfl override this if allow_cfl_dt_max_violate_dt_min is false?
        end
        # spinup
        if t < TS.spinup_half_t_max
            TS.dt *= TS.spinup_dt_factor # can we fade this out gradually between t=spinup_half_t_max and t=2*spinup_half_t_max?, and ramp up from Δt*spinup_dt_factor to Δt?
        else # if spinup_half_t_max is 0, we'll bypass this
            TS.dt *= TS.spinup_dt_factor + (1 - TS.spinup_dt_factor) * (t - TS.spinup_half_t_max) / TS.spinup_half_t_max # a smooth ramp from TS.spinup_dt_factor to 1 between t=spinup_half_t_max and t=2*spinup_half_t_max
        end

        if TS.use_fallback_during_spinup
            TS.use_fallback_tendency_limiters = true
        end
    end

    # Enforce this last, allow cfl dt min to override all other settings if applicable..., even if spinup_adapt_dt is not allowed (for now)
    if (TS.cfl_dt_max < TS.dt) && (TS.cfl_dt_max < TS.dt_min) && TS.allow_cfl_dt_max_violate_dt_min # is this right? is cfl_dt_max < dt a problem? bc dt hasn't necessarily been updated yet right? (actually it has if adapt_dt is true)
        @debug " t = $t | cfl_dt_max = $(TS.cfl_dt_max) is less than both dt = $(TS.dt) and dt_min = $(TS.dt_min), setting dt to cfl_dt_max."
        TS.dt = TS.cfl_dt_max
    end

    # actually, any use of basiclimiter should require recalculation of those tendencies if dt changes......

    # # Note if we're using methods that return tendencies as a function of Δt, those may need to be recalculated

    if TS.dt ≠ dt_old # we should recompute or tendencies and be sure
        reset_violate_dt_min_remaining::Bool = false
        recalculate_dt_edmf::Bool = false

        prog = integrator.u
        params = integrator.p

        tendencies = ODE.get_du(integrator) # should be integrator.fsallast for Euler
        if isnothing(tendencies)
            @warn "No tendencies available for timestep calculation, might be during initialization"
            return nothing
        else # This is an attempt to not recalc all the tendencies, if it becomes too hard we can just recalculate all of them with ∑tendencies!()
            # default_tendency_limiter = edmf.tendency_limiters.default_tendency_limiter
            # moisture_sources_limiter = edmf.tendency_limiters.moisture_sources_limiter
            # entr_detr_tendency_limiter = edmf.tendency_limiters.entr_detr_tendency_limiter
            # precipitation_tendency_limiter = edmf.tendency_limiters.precipitation_tendency_limiter

            default_tendency_limiter = TC.get_tendency_limiter(edmf.tendency_limiters, Val(:default), TS.use_fallback_tendency_limiters)
            moisture_sources_limiter = TC.get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources),  TS.use_fallback_tendency_limiters)
            entr_detr_tendency_limiter = TC.get_tendency_limiter(edmf.tendency_limiters, Val(:entr_detr),  TS.use_fallback_tendency_limiters)
            precipitation_tendency_limiter = TC.get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation),  TS.use_fallback_tendency_limiters)

            # if TS.use_fallback_tendency_limiters
            #     @info "fallback:  default_tendency_limiter  = $(default_tendency_limiter) | moisture_sources_limiter = $(moisture_sources_limiter) | entr_detr_tendency_limiter = $(entr_detr_tendency_limiter) | precipitation_tendency_limiter = $(precipitation_tendency_limiter)"
            # else
            #     @info "default_tendency_limiter  = $(default_tendency_limiter) | moisture_sources_limiter = $(moisture_sources_limiter) | entr_detr_tendency_limiter = $(entr_detr_tendency_limiter) | precipitation_tendency_limiter = $(precipitation_tendency_limiter)"
            # end

            if !isa(default_tendency_limiter, TC.NoTendencyLimiter) || # these are used currently for both up and gm. changing up feeds back on gm. changes in gm limit up tendencies. recalculation is best. 
                ((edmf.moisture_model isa TC.NonEquilibriumMoisture) && !isa(moisture_sources_limiter, TC.NoMoistureSourcesLimiter)) || !isa(default_tendency_limiter, TC.NoTendencyLimiter) # while we could try to use update_noneq_moisture_sources_tendencies!(), how that composes w/ gm tendencies is less clear depending on gm limiters. recalculation is best. This could be reduced to just an updated update_noneq_moisture_sources_tendencies!() in the [[ default_tendency_limiter isa NoTendencyLimiter ]] case.
                !isa(entr_detr_tendency_limiter, TC.NoTendencyLimiter) || # these are not stored anywhere, so changing them changes updraft which changes gm. recalculation is best.
                !isa(precipitation_tendency_limiter, TC.NoTendencyLimiter) # these also have an unclear relationship with gm. recalculation is best.

                reset_violate_dt_min_remaining = true
                recalculate_dt_edmf = true
                # update_gm_tendencies!(tendencies, prog, params, t) # update the tendencies again
                # keep track of where else that tendency gets used...

                # this is currently used for limiting the gm and updraft tendencies. Those feed back on each other in a way that would require iteration.
                # easiest thing to do is to just recalculate all the tendencies, but if we ever split out the limiter for gm, up, ∑tendencies, etc. we could just recalculate the gm tendencies maybe?
                @debug "Recalculating tendencies for adaptive timestep calculation iteration... "

                ∑tendencies!(tendencies, prog, params, t) # update the tendencies again, also does limiting etc
            else
                # This block could be valid if [[ default_tendency_limiter is NoTendencyLimiter]] and update_noneq_moisture_sources_tendencies!() is updated to properly update gm tendencies.
                if (edmf.moisture_model isa TC.NonEquilibriumMoisture) && !isa(moisture_sources_limiter, TC.NoMoistureSourcesLimiter) && isa(default_tendency_limiter, TC.NoTendencyLimiter)
                    reset_violate_dt_min_remaining = true
                    recalculate_dt_edmf = true
                    @debug "Recalculating noneq_moisture_sources tendencies for adaptive timestep calculation iteration... "
                    update_noneq_moisture_sources_tendencies!(tendencies, prog, params, t, TS.use_fallback_tendency_limiters) # update the tendencies again
                end
            end

            # recalculate tendencies (bc they change with dt)
            if reset_violate_dt_min_remaining
                if (TS.dt_max_edmf < dt_min) && (TS.N_dt_max_edmf_violate_dt_min_remaining > 0)
                    TS.N_dt_max_edmf_violate_dt_min_remaining += 1  # reset
                end
            end


            if recalculate_dt_edmf
                if depth < TS.adaptive_depth_limit
                    @debug "Recursing to find dt for adaptive timestep calculation w/ updated tendencies... , new depth = $(depth+1)"
                    adaptive_dt!(integrator; depth = depth+1) # recall this function again and pray it doesn't recurse indefinitely...
                    return nothing
                end
            end

        end
    end  


    # if after all the work/recursion etc we get down to the point where dt is fixed.
    # if we don't reach here, it should mean the model has not stabilized, right?
    # the alternative would be to decrement N_dt_max_edmf_violate_dt_min_remaining anbove and set use_fallback_tendency_limiters = true if it reaches 0, but then when recursing, unset it again (you'd still have to track if it was just set? actually no you'd just unset it because N_dt_max_edmf_violate_dt_min_remaining would be > 0, but potentially if N > 1, the call will be redundant... so maybe just do it here)
    if (TS.N_dt_max_edmf_violate_dt_min_remaining == 0) && (N_dt_max_edmf_violate_dt_min_remaining_old == 1)
        TS.N_dt_max_edmf_violate_dt_min_remaining = -1 # we're done with this, so we'll set it to -1 to indicate that we're done with it (so that we can check e.g. w/ use_fallback_during_spinup)
        TS.use_fallback_tendency_limiters = true # switch to fallback limiters if we're out of violations (this is needed so that timesteps are calculated with the fallback limiters from now on and, called here, is only done if the initial N_remaining was > 0)
    end

    if (TS.use_fallback_during_spinup) && (TS.spinup_half_t_max > 0) && (t < 2 * TS.spinup_half_t_max) &&  ( t + TS.dt > 2 * TS.spinup_half_t_max) # if we're about to go over the spinup time, reset dt to dt_min
        if !(TS.N_dt_max_edmf_violate_dt_min_remaining == -1) # if this was -1, it would mean we ran out of violations and we're using the fallback limiters by default now. otherwise, just leave as true
            TS.use_fallback_tendency_limiters = false # we don't want to do this if we've passed N_dt_max_edmf_violate_dt_min_remaining though...
        end
    end

    # we really should limit tendencies according to TS.dt...

    @debug "t = $t | TS.dt = $(TS.dt) (TS.dt_max_edmf = $(TS.dt_max_edmf), TS.cfl_dt_max = $(TS.cfl_dt_max)), TS.use_fallback_tendency_limiters  = $(TS.use_fallback_tendency_limiters)"
    
    SciMLBase.set_proposed_dt!(integrator, TS.dt)
    ODE.u_modified!(integrator, false)
end

# function spinup_dt!(integrator)
#     UnPack.@unpack edmf, TS, dt_min = integrator.p

#     t = integrator.t

#     # Reduce dt during spinup for stability | When adapt_dt is not true, TS.dt never gets updated and is just set to dt_min... so we need to be mindful of that 
#     if t < 2 * TS.spinup_half_t_max
#         if !TS.allow_spinup_adapt_dt # this will be false if adapt_dt is false, otherwise we could have an issue
#             TS.dt = TS.dt_min # reset dt to dt_min if adapt_dt is false so we don't repeatedly adjust it, if allow_spinup_adapt_dt is true, adapt_dt is true and adapt_dt gets set in adaptive_dt!() from driver/callbacks.jl
#         end
#         if t < TS.spinup_half_t_max
#             TS.dt *= TS.spinup_dt_factor # can we fade this out gradually between t=spinup_half_t_max and t=2*spinup_half_t_max?, and ramp up from Δt*spinup_dt_factor to Δt?
#         else # if spinup_half_t_max is 0, we'll bypass this
#             TS.dt *= TS.spinup_dt_factor + (1 - TS.spinup_dt_factor) * (t - TS.spinup_half_t_max) / TS.spinup_half_t_max # a smooth ramp from TS.spinup_dt_factor to 1 between t=spinup_half_t_max and t=2*spinup_half_t_max
#         end
#     end
#     SciMLBase.set_proposed_dt!(integrator, TS.dt)
#     ODE.u_modified!(integrator, false)
# end

"""
Helper function to compute the minimum time step for a given tendency, used in compute_tendency_dt_max()
"""
function Δt_max_helper(Δt_max::FT, up_tendency::FT, gm_tendency::FT, updraft_source::FT, env_source::FT; message::String="") where {FT}

    if up_tendency < zero(FT)
        if updraft_source > zero(FT) # could happen bc of the gm tendency addition
            # If updraft tendency is negative, we need the minimum of up / (-tend_up) and env / (tend_up - tend_gm).
            # If tend_gm is positive enough, (tend_up - tend_gm) becomes negative and the timesecale is Inf, achievable with max(up_tendency - gm_tendency, FT(0))
            if !iszero(env_source)
                env_time = env_source / (max(up_tendency - gm_tendency, FT(0)))
                if env_time < zero(FT)
                    env_time = FT(Inf)
                end
                Δt_max = min(Δt_max, -updraft_source / up_tendency, env_time) 
            else # there's no env_source (should be 0), we just have to hope that gm is providing updraft, don't want to have 0/0
                Δt_max = min(Δt_max, -updraft_source / up_tendency) 
            end
        else
            # the updraft can't consume anything but the tendency is negative, so this is an error
            # @warn("up_tendency $(up_tendency) is negative or 0 but updraft available to consume is zero or neg $(updraft_source) ... weird | $(message) ")
        end

    elseif up_tendency > zero(FT)
        if !iszero(env_source) # if env_source / (max(up_tendency - gm_tendency, FT(0)))) is negative, should mean inf
            env_time = env_source / (max(up_tendency - gm_tendency, FT(0)))
            if env_time < zero(FT)
                env_time = FT(Inf)
            end
            Δt_max = min(Δt_max, env_time)
        else
            # just have to pray and hope that the env has no sinks bc there's no stable step if it does.
            # @warn("up_tendency $(up_tendency) is positive but environment available $(env_source) is negative or 0 (probably because of a grid-mean tendency addition). Thus there is no stable timestep. Will leave Δt_max unchanged and continue...")
        end
    else
        #
    end

    return Δt_max
end

function Δt_max_helper(Δt_max::FT, up_tendency::FT, updraft_source::FT; message::String = "") where {FT}
    
    if up_tendency < zero(FT)
        if updraft_source > zero(FT) # could happen bc of the gm up_tendency addition
            Δt_max = min(Δt_max, -updraft_source / up_tendency)
        else
            # the updraft can't consume anything but the up_tendency is negative, so this is an error
            # @warn("Tendency $(up_tendency) is negative or 0 but updraft $(updraft_source) avalaible to consume is zero or neg ... weird [$(message)]")
        end
    else
        #
    end

    return Δt_max
end

"""
    Compute dt_max using tendencies and available sources/sinks
    TODO: Automate this to just loop over all prognostic variables (center and face) in a type-stable way...
    # -- hard to do bc there's dichotomy between e.g. [θ_liq_ice] and [ρaθ_liq_ice] or [ρ_c * area] and [ρ_c * area * θ_liq_ice]
"""
function compute_tendency_dt_max(state::TC.State, edmf::TC.EDMFModel)

    FT = TC.float_type(state)

    grid = TC.Grid(state)
    N_up = TC.n_updrafts(edmf)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)

    aux_tc = TC.center_aux_turbconv(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_en = TC.center_aux_environment(state)
    aux_en_f = TC.face_aux_environment(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_bulk = TC.center_aux_bulk(state)
    aux_bulk_f = TC.face_aux_bulk(state) # same as using aux_tc_f = face_aux_turbconv(state) and then using aux_tc_f.bulk.w for example
    prog_up = TC.center_prog_updrafts(state)
    prog_up_f = TC.face_prog_updrafts(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)

    tendencies_up = TC.center_tendencies_updrafts(state)
    tendencies_up_f = TC.face_tendencies_updrafts(state)

    # tendencies_bulk = TC.center_tendencies_bulk(state)
    # tendencies_bulk_f = TC.face_tendencies_bulk(state)

    tendencies_gm = TC.center_tendencies_grid_mean(state)

    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    Δt_max = FT(Inf)


    # === Individual Updrafts === #
    @inbounds for i in 1:N_up
        tends_ρarea = tendencies_up[i].ρarea
        tends_ρaθ_liq_ice = tendencies_up[i].ρaθ_liq_ice
        tends_ρaq_tot = tendencies_up[i].ρaq_tot
        # area = aux_up[i].area
        # θ_liq_ice = aux_up[i].θ_liq_ice
        # q_tot = aux_up[i].q_tot
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tends_ρaq_liq = tendencies_up[i].ρaq_liq
            tends_ρaq_ice = tendencies_up[i].ρaq_ice
            # q_liq = aux_up[i].q_liq
            # q_ice = aux_up[i].q_ice
        end
        tends_ρaw = tendencies_up_f[i].ρaw

        ρarea = prog_up[i].ρarea
        @inbounds for k in TC.real_center_indices(grid)
            old_Δt_max = Δt_max
            Δt_max = Δt_max_helper(Δt_max, tends_ρarea[k], prog_up[i].ρarea[k]; message = "ρarea")
            # if Δt_max < old_Δt_max
            #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | tends_ρarea[k] = $(tends_ρarea[k]) | ρarea[k] = $(prog_up[i].ρarea[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
            # end

            old_Δt_max = Δt_max
            Δt_max = Δt_max_helper(Δt_max, tends_ρaθ_liq_ice[k], prog_up[i].ρaθ_liq_ice[k]; message = "ρaθ_liq_ice")
            # if Δt_max < old_Δt_max
            #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | tends_ρaθ_liq_ice[k] = $(tends_ρaθ_liq_ice[k]) | ρaθ_liq_ice[k] = $(prog_up[i].ρaθ_liq_ice[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
            # end

            old_Δt_max = Δt_max
            Δt_max = Δt_max_helper(Δt_max, tends_ρaq_tot[k], prog_up[i].ρaq_tot[k]; message = "ρaq_tot")
            # if Δt_max < old_Δt_max
            #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | tends_ρaq_tot[k] = $(tends_ρaq_tot[k]) | ρaq_tot[k] = $(prog_up[i].ρaq_tot[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
            # end

            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                old_Δt_max = Δt_max
                Δt_max = Δt_max_helper(Δt_max, tends_ρaq_liq[k], prog_up[i].ρaq_liq[k]; message = "ρaq_liq")
                # if Δt_max < old_Δt_max
                #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | tends_ρaq_liq[k] = $(tends_ρaq_liq[k]) | ρaq_liq[k] = $(prog_up[i].ρaq_liq[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
                # end

                old_Δt_max = Δt_max
                Δt_max = Δt_max_helper(Δt_max, tends_ρaq_ice[k], prog_up[i].ρaq_ice[k]; message = "ρaq_ice")
                # if Δt_max < old_Δt_max
                #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | tends_ρaq_ice[k] = $(tends_ρaq_ice[k]) | ρaq_ice[k] = $(prog_up[i].ρaq_ice[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
                # end


                ρaq_vap = prog_up[i].ρaq_tot[k] - prog_up[i].ρaq_liq[k] - prog_up[i].ρaq_ice[k] # should this be prog? or aux?
                ρaq_vap_sat_liq = ρarea[k] * aux_up[i].q_vap_sat_liq[k] 
                ρaq_vap_sat_ice = ρarea[k] * aux_up[i].q_vap_sat_ice[k]
    
                if ρaq_vap > ρaq_vap_sat_liq # if supersaturated
                    ql_tendency_cond_evap = aux_up[i].ql_tendency_cond_evap[k]
                    if ql_tendency_cond_evap > 0 
                        old_Δt_max = Δt_max
                        Δt_max = Δt_max_helper(Δt_max, -ql_tendency_cond_evap, (ρaq_vap - ρaq_vap_sat_liq) / ρarea[k] ; message = "ρaq_vap_sat_liq") # ql tendency could come from vap or ice, but to be safe assume from vap, ql tendency is negative of the vap tendency
                        # if Δt_max < old_Δt_max
                        #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | ql_tendency_cond_evap = $(ql_tendency_cond_evap) | ρaq_vap = $(ρaq_vap) | ρaq_vap_sat_liq = $(ρaq_vap_sat_liq) | ρarea = $(ρarea[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
                        # end
                    end # otherwise, ql is shrinking despite supersaturation
                end
    
                if ρaq_vap > ρaq_vap_sat_ice # if supersaturated
                    qi_tendency_sub_dep = aux_up[i].qi_tendency_sub_dep[k]
                    if qi_tendency_sub_dep > 0
                        old_Δt_max = Δt_max
                        Δt_max = Δt_max_helper(Δt_max, -qi_tendency_sub_dep, (ρaq_vap - ρaq_vap_sat_ice) / ρarea[k] ; message = "ρaq_vap_sat_ice") # qi tendency could come from vap or liq, but to be safe assume from vap, qi tendency is negative of the vap tendency
                        # if Δt_max < old_Δt_max
                        #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | qi_tendency_sub_dep = $(qi_tendency_sub_dep) | ρaq_vap = $(ρaq_vap) | ρaq_vap_sat_ice = $(ρaq_vap_sat_ice) | ρarea = $(ρarea[k]) | z = $(grid.zc.z[k]) | T = $(aux_up[i].T[k])"
                        # end
                    end # otherwise, qi is shrinking despite supersaturation
                end

            end
            
        end

        ρaw = prog_up_f[i].ρaw
        area_f = aux_up_f[i].area # This should have been set in update_aux
        @inbounds for k in TC.real_face_indices(grid)

            old_Δt_max = Δt_max
            Δt_max = Δt_max_helper(Δt_max, tends_ρaw[k] , prog_up_f[i].ρaw[k]; message = "ρaw")
            # if Δt_max < old_Δt_max
            #     @warn "Δt_max decreased from $old_Δt_max to Δt_max = $Δt_max | tends_ρaw[k] = $(tends_ρaw[k]) | ρaw[k] = $(prog_up_f[i].ρaw[k]) | z = $(grid.zf.z[k])"
            # end
        end
    end



    ρaq_tot_en = aux_en.ρaq_tot # this is the prognostic value, not the clipped one

    # area_en = aux_en.area
    # θ_liq_ice_en = aux_en.θ_liq_ice
    # q_tot_en = aux_en.q_tot

    # tends_ρarea_bulk = tendencies_bulk.ρarea
    # tends_ρaθ_liq_ice_bulk = tendencies_bulk.ρaθ_liq_ice
    # tends_ρaq_tot_bulk = tendencies_bulk.ρaq_tot
    if edmf.moisture_model isa TC.NonEquilibriumMoisture
        # tends_ρaq_liq_bulk = tendencies_bulk.ρaq_liq
        # tends_ρaq_ice_bulk = tendencies_bulk.ρaq_ice

        q_liq_bulk = aux_bulk.q_liq
        q_ice_bulk = aux_bulk.q_ice
        # q_liq_en = aux_en.q_liq
        # q_ice_en = aux_en.q_ice
    end

    a_min = edmf.minimum_area # hard to use this for limiting bc prognostic doesn't respectit but aux does...
    a_max = edmf.max_area
    area_bulk = aux_bulk.area


    # - gm tendency added to updraft handles 
    # the env always feels (up tendency + gm tendency)
    # the updraft always feels (up tendency) only

    # if updraft tendency is positive, we need env / (tend_up + tend_gm)
    # if updraft tendency is negative, we need the minimum of up / (-tend_up) and env / (tend_up - tend_gm). if tend_gm is positive enough, (tend_up - tend_gm) becomes negative and the timesecale is Inf

    # the updraft never feels the gm tendency, so we don't need to consider it.

    # so really you need the smaller one...
    @inbounds for k in TC.real_center_indices(grid) 

        ρarea_bulk = sum(i-> prog_up[i].ρarea[k], 1:N_up)
        ρaθ_liq_ice_bulk = sum(i-> prog_up[i].ρaθ_liq_ice[k], 1:N_up)
        ρaq_tot_bulk = sum(i-> prog_up[i].ρaq_tot[k], 1:N_up)
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            ρaq_liq_bulk = sum(i-> prog_up[i].ρaq_liq[k], 1:N_up)
            ρaq_ice_bulk = sum(i-> prog_up[i].ρaq_ice[k], 1:N_up)
        end

        tends_ρarea_bulk = sum(i-> tendencies_up[i].ρarea[k], 1:N_up)
        tends_ρaθ_liq_ice_bulk = sum(i-> tendencies_up[i].ρaθ_liq_ice[k], 1:N_up)
        tends_ρaq_tot_bulk = sum(i-> tendencies_up[i].ρaq_tot[k], 1:N_up)
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tends_ρaq_liq_bulk = sum(i-> tendencies_up[i].ρaq_liq[k], 1:N_up)
            tends_ρaq_ice_bulk = sum(i-> tendencies_up[i].ρaq_ice[k], 1:N_up)
        end

        ρaq_tot_en = max(prog_gm.ρq_tot[k] - ρaq_tot_bulk[k], 0)
        ρaθ_liq_ice_en = max(prog_gm.ρaθ_liq_ice[k] - ρaθ_liq_ice_bulk[k], 0)
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            ρaq_liq_en = max(ρ_c[k] * prog_gm.q_liq[k] - ρaq_liq_bulk[k], 0)
            ρaq_ice_en = max(ρ_c[k] * prog_gm.ρq_ice[k] - ρaq_ice_bulk[k], 0)
        end



        Δt_max_old = Δt_max
        Δt_max = Δt_max_helper(Δt_max, tends_ρarea_bulk, FT(0), ρarea_bulk, ρ_c[k] * a_max - ρarea_bulk; message = "ρarea bulk")

        Δt_max_old = Δt_max
        Δt_max = Δt_max_helper(Δt_max, tends_ρaθ_liq_ice_bulk, tendencies_gm.ρθ_liq_ice[k], ρaθ_liq_ice_bulk, ρaθ_liq_ice_en; message = "ρaθ_liq_ice bulk")

        Δt_max_old = Δt_max
        Δt_max = Δt_max_helper(Δt_max, tends_ρaq_tot_bulk, tendencies_gm.ρq_tot[k], ρaq_tot_bulk, ρaq_tot_en; message = "ρaq_tot bulk")

        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            Δt_max_old = Δt_max
            Δt_max = Δt_max_helper(Δt_max, tends_ρaq_liq_bulk, ρ_c[k] * tendencies_gm.q_liq[k], ρaq_liq_bulk, ρaq_liq_en; message = "ρaq_liq bulk")

            Δt_max_old = Δt_max
            Δt_max = Δt_max_helper(Δt_max, tends_ρaq_ice_bulk, ρ_c[k] * tendencies_gm.q_ice[k], ρaq_ice_bulk, ρaq_ice_en; message = "ρaq_ice bulk")

            # Prevent positive supersaturation from being more than depleted in one timestep (negative supersaturation is limited by the condensate itself... ql_tendency_noneq and qi_tendency_noneq should already be 0 or positive (barring melting etc)
            # -- if we have very low supersaturation for example, we dont need it to take many timesteps to go from liq to ice by homogeneous freezing
            # Δt_max_old = Δt_max

            # bulk liq
            ρaq_vap = ρaq_tot_bulk - ρaq_liq_bulk - ρaq_ice_bulk # should this be prog? or aux?
            ρaq_vap_sat_liq = ρarea_bulk *  aux_bulk.q_vap_sat_liq[k]
            ρaq_vap_sat_ice = ρarea_bulk *  aux_bulk.q_vap_sat_ice[k]

            if ρaq_vap > ρaq_vap_sat_liq # if supersaturated
                ql_tendency_cond_evap = aux_bulk.ql_tendency_cond_evap[k]
                if ql_tendency_cond_evap > 0 
                    old_Δt_max = Δt_max
                    Δt_max = Δt_max_helper(Δt_max, -ql_tendency_cond_evap, (ρaq_vap - ρaq_vap_sat_liq) / ρarea_bulk ; message = "ρaq_vap_sat_liq bulk") # ql tendency could come from vap or ice, but to be safe assume from vap, ql tendency is negative of the vap tendency
                end # otherwise, ql is shrinking despite supersaturation
            end

            # bulk ice
            if ρaq_vap > ρaq_vap_sat_ice # if supersaturated
                qi_tendency_sub_dep = aux_bulk.qi_tendency_sub_dep[k]
                if qi_tendency_sub_dep > 0
                    old_Δt_max = Δt_max
                    Δt_max = Δt_max_helper(Δt_max, -qi_tendency_sub_dep, (ρaq_vap - ρaq_vap_sat_ice) / ρarea_bulk ; message = "ρaq_vap_sat_ice bulk") # qi tendency could come from vap or liq, but to be safe assume from vap, qi tendency is negative of the vap tendency
                end # otherwise, qi is shrinking despite supersaturation
            end


        end
    end

    area_bulk_f = TC.face_aux_turbconv(state).bulk.a_up # this originally was only updated in compute_diagnostics() but we'll do it here
    @. area_bulk_f = TC.ᶠinterp_a(area_bulk)

    @inbounds for k in TC.real_face_indices(grid)
        Δt_max_old = Δt_max
        ρaw_bulk = sum(i-> prog_up_f[i].ρaw[k], 1:N_up)
        tends_ρaw_bulk = sum(i-> tendencies_bulk_f[i].ρaw[k], 1:N_up) # this is the tendency of the bulk updrafts, not the face updrafts
        Δt_max = Δt_max_helper(Δt_max, tends_ρaw_bulk, ρaw_bulk; message = "ρaw bulk") # no env limit, env can go neg... (is that true?)
    end 



    return Δt_max

end

function compute_dt_max(state::TC.State, edmf::TC.EDMFModel, dt_max::FT, CFL_limit::FT, use_tendency_timestep_limiter::Bool) where {FT <: Real}
    grid = TC.Grid(state)

    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    Δzc = TC.get_Δz(prog_gm.ρ) # I'm not entirely sure what this returns, it's not the grid spacing precisely, some sort of weighted jacobian
    Δzf = TC.get_Δz(prog_gm_f.w) # I'm not entirely sure what this returns, it's not the grid spacing precisely, some sort of weighted jacobian

    N_up = TC.n_updrafts(edmf)

    aux_tc = TC.center_aux_turbconv(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_en_f = TC.face_aux_environment(state)

    aux_en = TC.center_aux_environment(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_gm = TC.center_aux_grid_mean(state)

    KM = aux_tc.KM
    KH = aux_tc.KH
    KQ = aux_tc.KQ

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow
    ε = FT(eps(FT)) 


    local cfl_dt_max::FT = FT(Inf)
    
    @inbounds for k in TC.real_face_indices(grid)
        TC.is_surface_face(grid, k) && continue
        @inbounds for i in 1:N_up
            dt_max = min(dt_max, CFL_limit * Δzf[k] / (abs(aux_up_f[i].w[k]) + ε))
        end
        dt_max = min(dt_max, CFL_limit * Δzf[k] / (abs(aux_en_f.w[k]) + ε))
    end
    @inbounds for k in TC.real_center_indices(grid)
        k_in = TC.is_toa_center(grid, k) ? k : (k + 1) # To match the way the terminal velocity is calculated in the monitor_cfl callback... I'm not sure which is right, so maybe you want to do k:(k+1) to be sure?
        k_out = k


        # the way this is implemented in the monitor_cfl callback, the terminal velocity should be shifted to k+1 relative to Δzc, and then use 0 term_vel for toa...
        # I'm not sure which is precisely right, maybe check both?
        # IIUC, the logic goes that 1 is sfc, k_max is toa, so when you do the z diffs, those are going downwards, and each dz diff at k is from k -> k+1 and so affected by the u(z) at k+1.
        # - Except at k_max (toa), there's no k + 1? (also idk how dz has the same length as z but ok...)
        if edmf.cloud_sedimentation_model isa TC.CloudSedimentationModel
            if edmf.cloud_sedimentation_model.grid_mean
                vel_max = max(term_vel_rain[k_in], term_vel_rain[k_out], term_vel_snow[k_in], term_vel_snow[k_out], aux_gm.term_vel_liq[k_in], aux_gm.term_vel_liq[k_out], aux_gm.term_vel_ice[k_in], aux_gm.term_vel_ice[k_out])
            else
                vel_max = max(term_vel_rain[k_in], term_vel_rain[k_out], term_vel_snow[k_in], term_vel_snow[k_out],
                    aux_en.term_vel_liq[k_in], aux_en.term_vel_liq[k_out], (aux_up[i].term_vel_liq[k_in] for i in 1:N_up)...,  (aux_up[i].term_vel_liq[k_out] for i in 1:N_up)...,
                    aux_en.term_vel_ice[k_in], aux_en.term_vel_ice[k_out], (aux_up[i].term_vel_ice[k_in] for i in 1:N_up)...,  (aux_up[i].term_vel_ice[k_out] for i in 1:N_up)...,
                    )
            end
        else
            vel_max = max(term_vel_rain[k_in], term_vel_rain[k_out], term_vel_snow[k_in],  term_vel_snow[k_out])
        end
        # Check terminal rain/snow velocity CFL
        dt_max = min(dt_max, CFL_limit * Δzc[k] / (vel_max + ε))
        cfl_dt_max = min(cfl_dt_max, CFL_limit * Δzc[k] / (vel_max + ε) ) # we have limited updrafts to cfl speeds, but if that is removed this could be a problem if combined with allow_cfl_dt_max_violate_dt_min = true
        # Check diffusion CFL (i.e., Fourier number)
        dt_max = min(dt_max, CFL_limit * Δzc[k]^2 / (max(KH[k], KM[k], KQ[k]) + ε))
    end

    if use_tendency_timestep_limiter
        dt_max_tendency = compute_tendency_dt_max(state, edmf)
        error("Deprecated this, to save space on not storing prog_en... if you ever bring that back, you can use this again")

        # dt_max_tendency = max(dt_max_tendency, ε) # arguably it doesn't matter, if we step tendency times dt and t doesnt change that's probably still ok
        # if dt_max_tendency < dt_max
        #     @warn "Tendency-limited timestep of $dt_max_tendency is smaller than CFL-limited timestep of $dt_max"
        # end
        dt_max = min(dt_max, dt_max_tendency)
    end

    @inbounds (max_term_vel_rain, ind_max) = findmax([aux_tc.term_vel_rain[k] for k in TC.real_center_indices(grid)])
    @inbounds inds = [k for k in TC.real_center_indices(grid)]

    @debug "dt_max = $dt_max | cfl_dt_max = $cfl_dt_max | max_term_vel_rain = $max_term_vel_rain, dz_max_termvel = $(Δzc[inds[ind_max]]))"
    return dt_max, cfl_dt_max
end

dt_max!(integrator) = dt_max!(integrator, integrator.p.TS.dt_max, integrator.p.TS.cfl_dt_max)
function dt_max!(integrator, dt_max, cfl_dt_max)
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, aux, TS) = integrator.p
        prog = integrator.u
        tendencies = ODE.get_du(integrator) # should be integrator.fsallast for Euler
        if isnothing(tendencies)
            @warn "No tendencies available for timestep calculation, might be during initialization"
            state = TC.column_prog_aux(prog, aux, colidx, integrator.p.calibrate_io)
        else
            state = TC.column_state(prog, aux, tendencies, colidx, integrator.p.calibrate_io) # use this so that we can use tendencies to control the timestep
        end
        (dt_max, cfl_dt_max) = compute_dt_max(state, edmf, dt_max, TS.cfl_limit, TS.use_tendency_timestep_limiter)
    end
    to_float(f) = f isa ForwardDiff.Dual ? ForwardDiff.value(f) : f
    (; TS) = integrator.p
    TS.dt_max_edmf = to_float(dt_max)
    TS.cfl_dt_max = to_float(cfl_dt_max)
    ODE.u_modified!(integrator, false)
end

"""
Call the ∑tendencies! function in a callback so you can move its position around
"""
function call_∑tendencies!(integrator)
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, aux, TS) = integrator.p
        tendencies = ODE.get_du(integrator) # should be integrator.fsallast for Euler
        prog = integrator.u
        params = integrator.p
        t = integrator.t
        ∑tendencies!(tendencies, prog, params, t)
    end
    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(integrator)
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, aux, TS) = integrator.p
        prog = integrator.u
        state = TC.column_prog_aux(prog, aux, colidx, integrator.p.calibrate_io)
        monitor_cfl!(state, edmf, TS.dt, TS.cfl_limit)
    end
    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(state, edmf, Δt, CFL_limit)
    grid = TC.Grid(state)
    prog_gm = TC.center_prog_grid_mean(state)
    Δz = TC.get_Δz(prog_gm.ρ) # I'm not entirely sure what this returns, it's not the grid spacing precisely, some sort of weighted jacobian
    FT = TC.float_type(state)
    frac_tol = FT(0.05)

    aux_tc = TC.center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    if edmf.cloud_sedimentation_model isa TC.CloudSedimentationModel
        if edmf.cloud_sedimentation_model.grid_mean
            aux_gm = TC.center_aux_grid_mean(state)
            term_vel_liqs = (aux_gm.term_vel_liq, )
            term_vel_ices = (aux_gm.term_vel_ice, )
            aux_ice = (aux_gm, ) # q_liq, q_ice are stored in prog... the values in aux are different... but prog_gm  also doesn't have q_tot, T etc...
        else
            aux_en = TC.center_aux_environment(state)
            N_up = TC.n_updrafts(edmf)
            aux_up = TC.center_aux_updrafts(state)
            @inbounds term_vel_liqs = (aux_en.term_vel_liq, (aux_up[i].term_vel_liq for i in 1:N_up)...)
            @inbounds term_vel_ices = (aux_en.term_vel_ice, (aux_up[i].term_vel_ice for i in 1:N_up)...)
            aux_ice = (aux_en, (aux_up[i] for i in 1:N_up)...)
        end
    else
        term_vel_liqs = ()
        term_vel_ices = ()
        aux_ice = ()
    end



    @inbounds for k in TC.real_center_indices(grid)
        # check stability criterion
        CFL_out_rain = Δt / Δz[k] * term_vel_rain[k]
        CFL_out_snow = Δt / Δz[k] * term_vel_snow[k]
        CFL_out_liqs = (Δt / Δz[k] * term_vel_liq[k] for term_vel_liq in term_vel_liqs)
        CFL_out_ices = (Δt / Δz[k] * term_vel_ice[k] for term_vel_ice in term_vel_ices)
        if TC.is_toa_center(grid, k)
            CFL_in_rain = FT(0)
            CFL_in_snow = FT(0)
            CFL_in_liqs = (FT(0) for _ in term_vel_liqs)
            CFL_in_ices = (FT(0) for _ in term_vel_ices)
        else
            CFL_in_rain = Δt / Δz[k] * term_vel_rain[k + 1]
            CFL_in_snow = Δt / Δz[k] * term_vel_snow[k + 1]
            CFL_in_liqs = (Δt / Δz[k] * term_vel_liq[k + 1] for term_vel_liq in term_vel_liqs)
            CFL_in_ices = (Δt / Δz[k] * term_vel_ice[k + 1] for term_vel_ice in term_vel_ices)
        end
        CFL_meteor, i_CFL_meteor = findmax((CFL_in_rain, CFL_out_rain, CFL_in_snow, CFL_out_snow, Iterators.flatten(zip(CFL_in_liqs, CFL_out_liqs))..., Iterators.flatten(zip(CFL_in_ices, CFL_out_ices))...))

        i_liq_0 = 4 + 1 # the first liquid index
        i_ice_0 = 4 + 2*length(term_vel_liqs) + 1 # the first ice index

        if CFL_meteor > CFL_limit * (1 + frac_tol)
            prog_pr = TC.center_prog_precipitation(state)
            prog_gm = TC.center_prog_grid_mean(state)
            aux_gm = TC.center_aux_grid_mean(state)
            aux_bulk = TC.center_aux_bulk(state)
            aux_en = TC.center_aux_environment(state)

            if i_CFL_meteor == 1
                term_vel = term_vel_rain[k + 1]
                species = "rain"
                q = prog_pr.q_rai[k + 1]
                q_liq, q_ice = prog_gm.q_liq[k + 1], prog_gm.q_ice[k + 1]
                q_tot = aux_gm.q_tot[k + 1]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k + 1], aux_bulk.q_liq[k + 1], aux_bulk.q_ice[k + 1]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k + 1], aux_en.q_liq[k + 1], aux_en.q_ice[k + 1]
                updraft_area = aux_bulk.area[k + 1]
                T = aux_gm.T[k + 1]
            elseif i_CFL_meteor == 2
                term_vel = term_vel_rain[k]
                species = "rain"
                q = prog_pr.q_rai[k]
                q_liq, q_ice = prog_gm.q_liq[k], prog_gm.q_ice[k]
                q_tot = aux_gm.q_tot[k]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k], aux_bulk.q_liq[k], aux_bulk.q_ice[k]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k]
                updraft_area = aux_bulk.area[k]
                T = aux_gm.T[k]
            elseif i_CFL_meteor == 3
                term_vel = term_vel_snow[k + 1]
                species = "snow"
                q = prog_pr.q_sno[k + 1]
                q_liq, q_ice = prog_gm.q_liq[k + 1], prog_gm.q_ice[k + 1]
                q_tot = aux_gm.q_tot[k + 1]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k + 1], aux_bulk.q_liq[k + 1], aux_bulk.q_ice[k + 1]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k + 1], aux_en.q_liq[k + 1], aux_en.q_ice[k + 1]
                updraft_area = aux_bulk.area[k + 1]
                T = aux_gm.T[k + 1]
            elseif i_CFL_meteor == 4
                term_vel = term_vel_snow[k]
                species = "snow"
                q = prog_pr.q_sno[k]
                q_liq, q_ice = prog_gm.q_liq[k], prog_gm.q_ice[k]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k], aux_bulk.q_liq[k], aux_bulk.q_ice[k]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k]
                updraft_area = aux_bulk.area[k]
                T = aux_gm.T[k]

            elseif 5 ≤ i_CFL_meteor < (5 + 2 * length(term_vel_liqs)) # 5 is the first liquid hydrometeor, and then we have 2 per updraft
                species = "liquid"
                i_liq = i_CFL_meteor - i_liq_0 + 1 # the index we are in the liquid list
                i_liq = (i_liq - 1) ÷ 2 + 1 # the index we are in the aux list
                k_liq = (i_CFL_meteor % 2) == 1 ? (k + 1) : k # continue the alternating in/out k/k+1 pattern

                term_vel = term_vel_liqs[i_liq][k_liq]
                q = aux_ice[i_liq].q_liq[k_liq]
                q_liq, q_ice = aux_ice[i_liq].q_liq[k_liq], aux_ice[i_liq].q_ice[k_liq]
                q_tot = aux_ice[i_liq].q_tot[k_liq]

                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k_liq], aux_bulk.q_liq[k_liq], aux_bulk.q_ice[k_liq]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k_liq], aux_en.q_liq[k_liq], aux_en.q_ice[k_liq]
                updraft_area = aux_bulk.area[k_liq]
                T = aux_gm.T[k_liq]
            else
                species = "ice"
                # i_ice = i_CFL_meteor - 4 # the index we are in the ice list [[ from when we had only ice]]
                i_ice = i_CFL_meteor - i_ice_0 + 1 # the index we are in the ice list
                i_ice = (i_ice - 1) ÷ 2 + 1 # the index we are in the aux list
                k_ice = (i_CFL_meteor % 2) == 1 ? (k + 1) : k # continue the alternating in/out k/k+1 pattern

                term_vel = term_vel_ices[i_ice][k_ice]
                q = aux_ice[i_ice].q_ice
                q_liq, q_ice = aux_ice[i_ice].q_liq[k_ice], aux_ice[i_ice].q_ice[k_ice]
                q_tot = aux_ice[i_ice].q_tot[k_ice]

                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k_ice], aux_bulk.q_liq[k_ice], aux_bulk.q_ice[k_ice]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k_ice], aux_en.q_liq[k_ice], aux_en.q_ice[k_ice]
                updraft_area = aux_bulk.area[k_ice]
                T = aux_gm.T[k_ice]

            end
            error(
                "Time step is too large for hydrometeor fall velocity! Hydrometeor CFL is $(CFL_meteor) > $(CFL_limit), Δt = $(Δt), Δz = $(Δz[k]), term_vel = $(term_vel), species = $(species), q_$(species) = $(q), q_liq = $(q_liq), q_ice = $(q_ice), q_tot = $(q_tot) | q_tot_up = $(q_tot_up), q_liq_up = $(q_liq_up), q_ice_up = $(q_ice_up) | q_tot_env = $(q_tot_env), q_liq_env = $(q_liq_env), q_ice_env = $(q_ice_env) | updraft_area = $(updraft_area) | T = $(T)",
            )
        end
    end
    return nothing
end
