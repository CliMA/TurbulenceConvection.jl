#=

The general idea here is that processes on the grid can los a lot of accuracy as the grid becomes coarse. 
We can do better by accounting for the processes at the edges of the grid.

While ordinarily, this is just increasing the resolution, it can play a bigger role where there are extrema in the grid.

Because it becomes vanishingly unlikely that a an extrema in the grid matches up with an extrema in the process,
we can reweight the processes to account for the possibility that the true extrema is anywhere within (or outside of) the grid cell.

=#

has_extrema(x_lo::FT, x_mid::FT, x_hi::FT) where {FT} = ((x_mid < x_lo) & (x_mid < x_hi)) || ((x_mid > x_lo) & (x_mid > x_hi))
has_extrema(x::Tuple{FT,FT,FT}) where {FT} = has_extrema(x[1], x[2], x[3])
linear_interpolate_extrapolate(z::FT, zs::Tuple{FT,FT}, vs::Tuple{FT,FT}) where {FT} = (z - zs[1]) / (zs[2] - zs[1]) * (vs[2] - vs[1]) + vs[1] # fast linear interpolation/extrapolation


"""


While it might seem atractive to loop over k inside the fcn, it's more versatile to not, and then you can insert it anywhere in the EDMF calculation stack.
This is actually not really optional for us as mph_neq is not stored, so we'd be recalculating it every time, and assigning storage to mph_neq w/ more allocations.
"""
function reweight_noneq_moisture_sources_for_grid(
    k::Cent{Int64}, # really this is kc.
    grid::Grid,
    param_set::APS,
    thermo_params::TDPS, # we already have it available so just pass it in instead of pullling it out of the param_set again
    aux_domain,
    aux_domain_f,
    ts_vec,
    mph_neq::NoneqMoistureSources{FT},
    nonequilibrium_moisture_scheme::AbstractRelaxationTimescaleType,
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    Δt::Real,
    ρ_c::CC.Fields.Field, # this is from gm so just pass In
    p_c::CC.Fields.Field, # this is from gm so just pass In
    w::CC.Fields.Field, # on centers, we've already done the interpolation from aux_<>_f so may as well reuse it
    ;
    reweight_extrema_only::Bool = true # should we pass this in or get it from param_set? maybe make it explicity so you can cache it slash it's more viislbe in the code
    ) where {FT}

    # reweight_extrema_only = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false) # pass in from outside
    ε = eps(FT)


    # check for extrema in prognostic variables [θ_li, qt, ql, qi]
    # Ideally we'd do [T] bc supersat and thus q_vap_sat depend on T, but it's derived from the other variables so we can't do that.
    # we would still need to recalculate q_vap_sat_liq and q_vap_sat_ice, but we can do that in the extrapolation step.

    # we can't be withn one grid point of surf or top of atm
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    

    # reweight [if we reweight everywhere] OR [if we reweight extrema only, we're not at the surface or top of atm, and we have extrema in the prognostic variables]
    if (!reweight_extrema_only) || 
        (reweight_extrema_only &&  !((k == kc_surf) || (k == kc_toa))  && any(has_extrema,
            ((aux_domain.θ_liq_ice[k - 1], aux_domain.θ_liq_ice[k], aux_domain.θ_liq_ice[k + 1]),
            (aux_domain.q_tot[k - 1], aux_domain.q_tot[k], aux_domain.q_tot[k + 1]),
            (aux_domain.q_liq[k - 1], aux_domain.q_liq[k], aux_domain.q_liq[k + 1]), # idk if this really matters, doesn't impact thermo/supersat
            (aux_domain.q_ice[k - 1], aux_domain.q_ice[k], aux_domain.q_ice[k + 1]), # idk if this really matters, doesn't impact thermo/supersat
            (aux_domain.T[k - 1], aux_domain.T[k], aux_domain.T[k + 1])),
            ( TD.vapor_specific_humidity(thermo_params, ts_vec[k]) - aux_domain.q_vap_sat_liq[k - 1],  TD.vapor_specific_humidity(thermo_params, ts_vec[k - 1]) - aux_domain.q_vap_sat_liq[k],  TD.vapor_specific_humidity(thermo_params, ts_vec[k + 1]) - aux_domain.q_vap_sat_liq[k + 1]), # liq supersat (bc both q_vap_sat and T are prolly falling off pretty fast, we can have supersat peaks even w/o peaks in q, T, esp bc theta has even fewer peaks than T)
            ( TD.vapor_specific_humidity(thermo_params, ts_vec[k]) - aux_domain.q_vap_sat_ice[k - 1],  TD.vapor_specific_humidity(thermo_params, ts_vec[k - 1]) - aux_domain.q_vap_sat_ice[k],  TD.vapor_specific_humidity(thermo_params, ts_vec[k + 1]) - aux_domain.q_vap_sat_ice[k + 1]), # ice supersat
            )
        )

        #=
        For anything in the `hi` region, the effect inside our cell is the same. Similarly inside the `lo` region. (the continued extrapolation doesn't impact inside this cell, only the next one.)

        We weight by distances <center to face, face to face, face to center>. While it might seem overly weighted towards the central value (given we assume the cener value for that whole window), sharp gradients are probably less likely so it's probably a fair trade.
        This is a probability weight.

        So we have 3 worlds, with weights prob_weight.
        1. Extreme is in the high region, both extrap from center and bottom center
        2. Extreme is in the center regular extrap. regular extrap from next centers
        3. Extreme is in the low region, both extrap from center and top center

        We assume the area weights are (; lo = 0.25, mid = 0.5, hi = 0.25). This is a `closest neighbor` assumption.
        =#

        kf = CCO.PlusHalf(k.i)

        if !(k == kc_toa)
            zc_hi = FT(grid.zc[k+1].z)
        end

        zf_hi = FT(grid.zf[kf+1].z)
        zc = FT(grid.zc[k].z)
        zf_lo = FT(grid.zf[kf].z)

        if !(k == kc_surf)
            zc_lo = FT(grid.zc[k-1].z)
        end

            
        if k == kc_surf # No extrema, and you can only extrapolate from above... so you have one move
            # these are no longer probabilities of where the extrema are but just 
            scenarios = (
                (; top_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), bottom_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), prob = FT(1) ), # extreme in top
                )
        elseif k == kc_toa # No extrema, and you can only extrapolate from below... so you have one move
            scenarios = (
                (; top_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), bottom_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), prob = FT(1) ), # extreme in bottom
                )
        else
            scenarios = (
                (; top_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), bottom_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), prob = (zc_hi - zf_hi) / (zc_hi - zc_lo) ), # extreme in top / extrap from bottom and center always
                (; top_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), bottom_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), prob = (zf_hi - zf_lo) / (zc_hi - zc_lo)), # extreme in center / extrap from middle out
                (; top_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), bottom_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), prob = (zf_lo - zc_lo) / (zc_hi - zc_lo)), # extreme in bottom / extrap from top and center always
                )
        end

        # (isone(sum(scenario->scenario.prob, scenarios))) || error("Probabilities do not sum to 1") # for debugging but we're good now

        area_weight = (;lo = FT(0.25), mid = FT(0.5), hi = FT(0.25)) # we'll keep this both for extrema and otherwise, as well as for top and bottom value of each face should count as much as each center in the end, really just an implied res doubling.

        qi_tendency = FT(0)
        ql_tendency = FT(0)

        for scenario in scenarios

            θ_liq_ice_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.θ_liq_ice[scenario.bottom_extrap.ks[1]], aux_domain.θ_liq_ice[scenario.bottom_extrap.ks[2]]))
            θ_liq_ice_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.θ_liq_ice[scenario.top_extrap.ks[1]], aux_domain.θ_liq_ice[scenario.top_extrap.ks[2]]))
            q_tot_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.q_tot[scenario.bottom_extrap.ks[1]], aux_domain.q_tot[scenario.bottom_extrap.ks[2]]))
            q_tot_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.q_tot[scenario.top_extrap.ks[1]], aux_domain.q_tot[scenario.top_extrap.ks[2]]))
            q_liq_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.q_liq[scenario.bottom_extrap.ks[1]], aux_domain.q_liq[scenario.bottom_extrap.ks[2]]))
            q_liq_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.q_liq[scenario.top_extrap.ks[1]], aux_domain.q_liq[scenario.top_extrap.ks[2]]))
            q_ice_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.q_ice[scenario.bottom_extrap.ks[1]], aux_domain.q_ice[scenario.bottom_extrap.ks[2]]))
            q_ice_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.q_ice[scenario.top_extrap.ks[1]], aux_domain.q_ice[scenario.top_extrap.ks[2]]))

            p_c_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (p_c[scenario.bottom_extrap.ks[1]], p_c[scenario.bottom_extrap.ks[2]]))
            p_c_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (p_c[scenario.top_extrap.ks[1]], p_c[scenario.top_extrap.ks[2]]))

            ρ_c_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (ρ_c[scenario.bottom_extrap.ks[1]], ρ_c[scenario.bottom_extrap.ks[2]]))
            ρ_c_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (ρ_c[scenario.top_extrap.ks[1]], ρ_c[scenario.top_extrap.ks[2]]))

            area_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.area[scenario.bottom_extrap.ks[1]], aux_domain.area[scenario.bottom_extrap.ks[2]]))
            area_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.area[scenario.top_extrap.ks[1]], aux_domain.area[scenario.top_extrap.ks[2]]))

            # For w we do not extrapolate because we have it on the faces already
            w_lo = aux_domain_f.w[kf]
            w_hi = aux_domain_f.w[kf+1]



            ts_lo = thermo_state_pθq(param_set, p_c_lo, θ_liq_ice_lo, q_tot_lo, q_liq_lo, q_ice_lo)
            ts_hi = thermo_state_pθq(param_set, p_c_hi, θ_liq_ice_hi, q_tot_hi, q_liq_hi, q_ice_hi)
            T_lo = TD.air_temperature(thermo_params, ts_lo)
            T_hi = TD.air_temperature(thermo_params, ts_hi)

        
            # These are extrapolated points so we need to recalculate...
            q_vap_sat_liq_lo = TD.q_vap_saturation_generic(thermo_params, T_lo, ρ_c_lo, TD.Liquid())
            q_vap_sat_ice_lo = TD.q_vap_saturation_generic(thermo_params, T_lo, ρ_c_lo, TD.Ice())
            q_vap_sat_liq_hi = TD.q_vap_saturation_generic(thermo_params, T_hi, ρ_c_hi, TD.Liquid())
            q_vap_sat_ice_hi = TD.q_vap_saturation_generic(thermo_params, T_hi, ρ_c_hi, TD.Ice())

            mph_neq_lo = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, area_lo, ρ_c_lo, Δt + ε, ts_lo, w_lo, zf_lo, q_vap_sat_liq_lo, q_vap_sat_ice_lo) #ts_LCL = nothing)
            mph_neq_hi = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, area_hi, ρ_c_hi, Δt + ε, ts_hi, w_hi, zf_hi, q_vap_sat_liq_hi, q_vap_sat_ice_hi) #ts_LCL = nothing)

            #=
                There's an argument that we should weigh [mph_neq_lo, mph, mph_neq_hi] by the some fraction of [area_lo, area, area_hi]...  so that in the end multiplying by area gives a weighted answer
                I guess you could also argue that they're just local relative rates and do composse normally... but i think that fals apart as the gradients get sharper, e.g. right below an upddraft spike the are areally does matter...
            =#
            if area_lo > FT(0)
                f = (aux_domain.area[k] / area_lo) # then, tendency times aux_domain.area gives the true total from the top edge
                mph_neq_lo = NoneqMoistureSources{FT}(mph_neq_lo.ql_tendency * f, mph_neq_lo.qi_tendency * f) # don't edit this bc it needs to go in if it's to be counted as coming out of qi in precip_formation()
            end
            if area_hi > FT(0)
                f = (aux_domain.area[k] / area_hi) # then, tendency times aux_domain.area gives the true total from the bottom edge
                mph_neq_hi = NoneqMoistureSources{FT}(mph_neq_hi.ql_tendency * f, mph_neq_hi.qi_tendency * f) # don't edit this bc it needs to go in if it's to be counted as coming out of qi in precip_formation()
            end

            ql_tendency += scenario.prob * (mph_neq_lo.ql_tendency * area_weight.lo + mph_neq.ql_tendency * area_weight.mid + mph_neq_hi.ql_tendency * area_weight.hi)
            qi_tendency += scenario.prob * (mph_neq_lo.qi_tendency * area_weight.lo + mph_neq.qi_tendency * area_weight.mid + mph_neq_hi.qi_tendency * area_weight.hi)


            if any(!isfinite, (ql_tendency, qi_tendency))
                @error "Non-equilibrium moisture sources returned non-finite values: $mph_neq; from inputs ts = $ts; w = $w[k]; zc = $zc; q_vap_sat_liq = $q_vap_sat_liq; q_vap_sat_ice = $q_vap_sat_ice; ρ_c = $(ρ_c[k]); aux_domain.area = $(aux_domain.area[k])"
                @info "parts were mph_neq_lo = $mph_neq_lo; mph_neq = $mph_neq; mph_neq_hi = $mph_neq_hi"
                @info "more parts were θ_liq_ice_lo = $θ_liq_ice_lo; θ_liq_ice_hi = $θ_liq_ice_hi; q_tot_lo = $q_tot_lo; q_tot_hi = $q_tot_hi; q_liq_lo = $q_liq_lo; q_liq_hi = $q_liq_hi; q_ice_lo = $q_ice_lo; q_ice_hi = $q_ice_hi; p_c_lo = $p_c_lo; p_c_hi = $p_c_hi; ρ_c_lo = $ρ_c_lo; ρ_c_hi = $ρ_c_hi; area_lo = $area_lo; area_hi = $area_hi; w_lo = $w_lo; w_hi = $w_hi; ts_lo = $ts_lo; ts_hi = $ts_hi; T_lo = $T_lo; T_hi = $T_hi; q_vap_sat_liq_lo = $q_vap_sat_liq_lo; q_vap_sat_ice_lo = $q_vap_sat_ice_lo; q_vap_sat_liq_hi = $q_vap_sat_liq_hi; q_vap_sat_ice_hi = $q_vap_sat_ice_hi"
                error("Non-equilibrium moisture sources returned non-finite values")
            end

        end

        mph_neq = NoneqMoistureSources{FT}(ql_tendency, qi_tendency) # don't edit this bc it needs to go in if it's to be counted as coming out of qi in precip_formation()          

    end

    return mph_neq
end


"""
This one just calculates saturation adjustment in the same way as the noneq but doesn't calculate tendencies, instead just the literal q values
"""
function reweight_equilibrium_saturation_adjustment_for_grid(
    k::Cent{Int64}, # really this is kc.
    grid::Grid,
    param_set::APS,
    thermo_params::TDPS, # we already have it available so just pass it in instead of pullling it out of the param_set again
    aux_domain,
    # ts_vec, # we would use it but it's not stored for the updraft
    q::TD.PhasePartition, # alreayd calcluated, may as well pass it in
    p_c::CC.Fields.Field, # this is from gm so just pass In
    ;
    reweight_extrema_only::Bool = true # should we pass this in or get it from param_set? maybe make it explicity so you can cache it slash it's more viislbe in the code
    )

    FT = eltype(param_set)

    # reweight_extrema_only = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false) # pass in from outside


    # So here, it's mostly the same but we need to do saturation adjustment and weight. These will never feed back to the prognostic state.

    # check for extrema in prognostic variables [θ_li, qt, ql, qi]
    # Ideally we'd do [T] bc supersat and thus q_vap_sat depend on T, but it's derived from the other variables so we can't do that.
    # we would still need to recalculate q_vap_sat_liq and q_vap_sat_ice, but we can do that in the extrapolation step.

    # we can't be withn one grid point of surf or top of atm
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)

    # reweight [if we reweight everywhere] OR [if we reweight extrema only, we're not at the surface or top of atm, and we have extrema in the prognostic variables]
    if (!reweight_extrema_only) || 
        (reweight_extrema_only &&  !((k == kc_surf) || (k == kc_toa))  && any(has_extrema,
            ((aux_domain.θ_liq_ice[k - 1], aux_domain.θ_liq_ice[k], aux_domain.θ_liq_ice[k + 1]),
            (aux_domain.q_tot[k - 1], aux_domain.q_tot[k], aux_domain.q_tot[k + 1]),
            # (aux_domain.q_liq[k - 1], aux_domain.q_liq[k], aux_domain.q_liq[k + 1]), # idk if this really matters, doesn't impact thermo/supersat [ if we want this in sat, we need to move this fully out of the loop  where we update aux_en/aux_up
            # (aux_domain.q_ice[k - 1], aux_domain.q_ice[k], aux_domain.q_ice[k + 1]), # idk if this really matters, doesn't impact thermo/supersat [ if we want this in sat, we need to move this fully out of the loop  where we update aux_en/aux_up]
            (aux_domain.T[k - 1], aux_domain.T[k], aux_domain.T[k + 1])),
            )
        )

        #=
        For anything in the `hi` region, the effect inside our cell is the same. Similarly inside the `lo` region. (the continued extrapolation doesn't impact inside this cell, only the next one.)

        We weight by distances <center to face, face to face, face to center>. While it might seem overly weighted towards the central value (given we assume the cener value for that whole window), sharp gradients are probably less likely so it's probably a fair trade.
        This is a probability weight.

        So we have 3 worlds, with weights prob_weight.
        1. Extreme is in the high region, both extrap from center and bottom center
        2. Extreme is in the center regular extrap. regular extrap from next centers
        3. Extreme is in the low region, both extrap from center and top center

        We assume the area weights are (; lo = 0.25, mid = 0.5, hi = 0.25). This is a `closest neighbor` assumption.
        =#

        kf = CCO.PlusHalf(k.i)

        if !(k == kc_toa)
            zc_hi = FT(grid.zc[k+1].z)
        end

        zf_hi = FT(grid.zf[kf+1].z)
        zc = FT(grid.zc[k].z)
        zf_lo = FT(grid.zf[kf].z)

        if !(k == kc_surf)
            zc_lo = FT(grid.zc[k-1].z)
        end

            
        if k == kc_surf # No extrema, and you can only extrapolate from above... so you have one move
            # these are no longer probabilities of where the extrema are but just 
            scenarios = (
                (; top_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), bottom_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), prob = FT(1) ), # extreme in top
                )
        elseif k == kc_toa # No extrema, and you can only extrapolate from below... so you have one move
            scenarios = (
                (; top_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), bottom_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), prob = FT(1) ), # extreme in bottom
                )
        else
            scenarios = (
                (; top_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), bottom_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), prob = (zc_hi - zf_hi) / (zc_hi - zc_lo) ), # extreme in top / extrap from bottom and center always
                (; top_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), bottom_extrap = (;ks=(k-1, k), zs=(zc_lo, zc)), prob = (zf_hi - zf_lo) / (zc_hi - zc_lo)), # extreme in center / extrap from middle out
                (; top_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), bottom_extrap = (;ks=(k, k+1), zs=(zc, zc_hi)), prob = (zf_lo - zc_lo) / (zc_hi - zc_lo)), # extreme in bottom / extrap from top and center always
                )
        end

        # (isone(sum(scenario->scenario.prob, scenarios))) || error("Probabilities do not sum to 1") # for debugging but we're good now

        area_weight = (;lo = FT(0.25), mid = FT(0.5), hi = FT(0.25)) # we'll keep this both for extrema and otherwise, as well as for top and bottom value of each face should count as much as each center in the end, really just an implied res doubling.

        ql = FT(0)
        qi = FT(0)

        for scenario in scenarios

            θ_liq_ice_lo = linear_interpolate_extrapolate(zf_lo,  scenario.bottom_extrap.zs, (aux_domain.θ_liq_ice[scenario.bottom_extrap.ks[1]], aux_domain.θ_liq_ice[scenario.bottom_extrap.ks[2]]))
            θ_liq_ice_hi = linear_interpolate_extrapolate(zf_hi,  scenario.top_extrap.zs, (aux_domain.θ_liq_ice[scenario.top_extrap.ks[1]], aux_domain.θ_liq_ice[scenario.top_extrap.ks[2]]))
            q_tot_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.q_tot[scenario.bottom_extrap.ks[1]], aux_domain.q_tot[scenario.bottom_extrap.ks[2]]))
            q_tot_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.q_tot[scenario.top_extrap.ks[1]], aux_domain.q_tot[scenario.top_extrap.ks[2]]))

            p_c_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (p_c[scenario.bottom_extrap.ks[1]], p_c[scenario.bottom_extrap.ks[2]]))
            p_c_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (p_c[scenario.top_extrap.ks[1]], p_c[scenario.top_extrap.ks[2]]))


            ts_lo = thermo_state_pθq(param_set, p_c_lo, θ_liq_ice_lo, q_tot_lo)
            ts_hi = thermo_state_pθq(param_set, p_c_hi, θ_liq_ice_hi, q_tot_hi)


            q_lo = TD.PhasePartition(thermo_params, ts_lo) # we have to construct the phase partition based on q_tot in Eq
            q_hi = TD.PhasePartition(thermo_params, ts_hi) # we have to construct the phase partition based on q_tot in Eq

            

            ql += scenario.prob * (q_lo.liq * area_weight.lo + q.liq * area_weight.mid + q_hi.liq * area_weight.hi)
            qi += scenario.prob * (q_lo.ice * area_weight.lo + q.ice * area_weight.mid + q_hi.ice * area_weight.hi)


            if any(!isfinite, (ql, qi))
                @error "Sat adjust moisture sources returned non-finite values: $q; from inputs ts = $ts; w = $w[k]; zc = $zc; aux_domain.area = $(aux_domain.area[k])"
                @info "parts were q_lo = $q_lo; q = $q; q_hi = $q_hi"
                @info "more parts were θ_liq_ice_lo = $θ_liq_ice_lo; θ_liq_ice_hi = $θ_liq_ice_hi; q_tot_lo = $q_tot_lo; q_tot_hi = $q_tot_hi; p_c_lo = $p_c_lo; p_c_hi = $p_c_hi; ts_lo = $ts_lo; ts_hi = $ts_hi"
                error("Equilibrium sat adjust returned non-finite values")
            end

        end
  

    end

    return ql, qi

end
