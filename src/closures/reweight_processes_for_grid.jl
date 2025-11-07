#=

The general idea here is that processes on the grid can los a lot of accuracy as the grid becomes coarse. 
We can do better by accounting for the processes at the edges of the grid.

While ordinarily, this is just increasing the resolution, it can play a bigger role where there are extrema in the grid.

Because it becomes vanishingly unlikely that a an extrema in the grid matches up with an extrema in the process,
we can reweight the processes to account for the possibility that the true extrema is anywhere within (or outside of) the grid cell.

=#

has_extrema(x_lo::FT, x_mid::FT, x_hi::FT) where {FT} = ((x_mid < x_lo) & (x_mid < x_hi)) || ((x_mid > x_lo) & (x_mid > x_hi))
has_extrema(x::Tuple{FT,FT,FT}) where {FT} = has_extrema(x[1], x[2], x[3])

"""


While it might seem atractive to loop over k inside the fcn, it's more versatile to not, and then you can insert it anywhere in the EDMF calculation stack.
This is actually not really optional for us as mph_neq is not stored, so we'd be recalculating it every time, and assigning storage to mph_neq w/ more allocations.
"""
function reweight_noneq_moisture_sources_for_grid(
    k::Cent{Int64}, # really this is kc.
    grid::Grid,
    param_set::APS,
    thermo_params::TDPS, # we already have it available so just pass it in instead of pullling it out of the param_set again
    aux_domain::CC.Fields.Field, # FieldVector was wrong..
    aux_domain_f::CC.Fields.Field, # this is the aux domain on faces, we need it to get the aux variables on faces
    mph_neq::NoneqMoistureSources{FT},
    nonequilibrium_moisture_scheme::AbstractRelaxationTimescaleType,
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    Δt::Real,
    ρ_c::CC.Fields.Field, # this is from gm so just pass In
    p_c::CC.Fields.Field, # this is from gm so just pass In
    w::CC.Fields.Field, # on centers, we've already done the interpolation from aux_<>_f so may as well reuse it
    dqvdt::CC.Fields.Field, # this is the tendency of q_vap, if we have it, we can use it to reweight the sources
    dTdt::CC.Fields.Field, # this is the tendency of T, if we have it, we can use it to reweight the sources
    ;
    reweight_extrema_only::Bool = true, # should we pass this in or get it from param_set? maybe make it explicity so you can cache it slash it's more viislbe in the code
    region::AbstractDomain = Env,
    ) where {FT}

    microphys_params = TCP.microphysics_params(param_set)

    if iszero(aux_domain.area[k]) # bc our output is indexed to this point, aint nun we can do for ya. dw tho, a neighbor nonzero point can still benefit when called from there. unfortunate edge case, this is.
        return mph_neq
    end

    ε = eps(FT)


    # no cloak
    if !((k == kc_surf) || (k == kc_toa))
        if (region isa EnvOrUp) || (region isa EnvRemainingDomain) # no cloak
            q_tot_km1 = aux_domain.q_tot[k - 1]
            q_tot = aux_domain.q_tot[k]
            q_tot_kp1 = aux_domain.q_tot[k + 1]
            #
            q_liq_km1 = aux_domain.q_liq[k - 1]
            q_liq = aux_domain.q_liq[k]
            q_liq_kp1 = aux_domain.q_liq[k + 1]
            #
            q_ice_km1 = aux_domain.q_ice[k - 1]
            q_ice = aux_domain.q_ice[k]
            q_ice_kp1 = aux_domain.q_ice[k + 1]
            #
            θ_liq_ice_km1 = aux_domain.θ_liq_ice[k - 1]
            θ_liq_ice = aux_domain.θ_liq_ice[k]
            θ_liq_ice_kp1 = aux_domain.θ_liq_ice[k + 1]
            #
            T_km1 = aux_domain.T[k - 1]
            T = aux_domain.T[k]
            T_kp1 = aux_domain.T[k + 1]
            #
            ts_km1 = aux_domain.ts[k - 1]
            ts = aux_domain.ts[k]
            ts_kp1 = aux_domain.ts[k + 1]
            #
            q_vap_sat_liq_km1 = aux_domain.q_vap_sat_liq[k - 1]
            q_vap_sat_liq = aux_domain.q_vap_sat_liq[k]
            q_vap_sat_liq_kp1 = aux_domain.q_vap_sat_liq[k + 1]
        elseif (region isa CloakUp) # cloak up
            q_tot_km1 = aux_domain.q_tot_cloak_up[k - 1]
            q_tot = aux_domain.q_tot_cloak_up[k]
            q_tot_kp1 = aux_domain.q_tot_cloak_up[k + 1]
            #
            q_liq_km1 = aux_domain.q_liq_cloak_up[k - 1]
            q_liq = aux_domain.q_liq_cloak_up[k]
            q_liq_kp1 = aux_domain.q_liq_cloak_up[k + 1]
            #
            q_ice_km1 = aux_domain.q_ice_cloak_up[k - 1]
            q_ice = aux_domain.q_ice_cloak_up[k]
            q_ice_kp1 = aux_domain.q_ice_cloak_up[k + 1]
            #
            θ_liq_ice_km1 = aux_domain.θ_liq_ice_cloak_up[k - 1]
            θ_liq_ice = aux_domain.θ_liq_ice_cloak_up[k]
            θ_liq_ice_kp1 = aux_domain.θ_liq_ice_cloak_up[k + 1]
            #
            T_km1 = aux_domain.T_cloak_up[k - 1]
            T = aux_domain.T_cloak_up[k]
            T_kp1 = aux_domain.T_cloak_up[k + 1]
            #
            ts_km1 = aux_domain.ts_cloak_up[k - 1]
            ts = aux_domain.ts_cloak_up[k]
            ts_kp1 = aux_domain.ts_cloak_up[k + 1]
            #
            q_vap_sat_liq_km1 = TD.q_vap_saturation_generic(thermo_params, T_km1, TD.air_density(thermo_params, ts_km1), TD.Liquid())
            q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, TD.air_density(thermo_params, ts), TD.Liquid())
            q_vap_sat_liq_kp1 = TD.q_vap_saturation_generic(thermo_params, T_kp1, TD.air_density(thermo_params, ts_kp1), TD.Liquid())
        elseif (region isa CloakDown) # cloak down
            q_tot_km1 = aux_domain.q_tot_cloak_down[k - 1]
            q_tot = aux_domain.q_tot_cloak_down[k]
            q_tot_kp1 = aux_domain.q_tot_cloak_down[k + 1]
            #
            q_liq_km1 = aux_domain.q_liq_cloak_down[k - 1]
            q_liq = aux_domain.q_liq_cloak_down[k]
            q_liq_kp1 = aux_domain.q_liq_cloak_down[k + 1]
            #
            q_ice_km1 = aux_domain.q_ice_cloak_down[k - 1]
            q_ice = aux_domain.q_ice_cloak_down[k]
            q_ice_kp1 = aux_domain.q_ice_cloak_down[k + 1]
            #
            θ_liq_ice_km1 = aux_domain.θ_liq_ice_cloak_down[k - 1]
            θ_liq_ice = aux_domain.θ_liq_ice_cloak_down[k]
            θ_liq_ice_kp1 = aux_domain.θ_liq_ice_cloak_down[k + 1]
            #
            T_km1 = aux_domain.T_cloak_down[k - 1]
            T = aux_domain.T_cloak_down[k]
            T_kp1 = aux_domain.T_cloak_down[k + 1]
            #
            ts_km1 = aux_domain.ts_cloak_down[k - 1]
            ts = aux_domain.ts_cloak_down[k]
            ts_kp1 = aux_domain.ts_cloak_down[k + 1]
            #
            q_vap_sat_liq_km1 = TD.q_vap_saturation_generic(thermo_params, T_km1, TD.air_density(thermo_params, ts_km1), TD.Liquid())
            q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, TD.air_density(thermo_params, ts), TD.Liquid())
            q_vap_sat_liq_kp1 = TD.q_vap_saturation_generic(thermo_params, T_kp1, TD.air_density(thermo_params, ts_kp1), TD.Liquid())
        end
    end

    # check for extrema in prognostic variables [θ_li, qt, ql, qi]
    # Ideally we'd do [T] bc supersat and thus q_vap_sat depend on T, but it's derived from the other variables so we can't do that.
    # we would still need to recalculate q_vap_sat_liq and q_vap_sat_ice, but we can do that in the extrapolation step.

    # we can't be withn one grid point of surf or top of atm
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)

    # reweight [if we reweight everywhere] OR [if we reweight extrema only, we're not at the surface or top of atm, and we have extrema in the prognostic variables]
    if (!reweight_extrema_only) || 
        (reweight_extrema_only &&  !((k == kc_surf) || (k == kc_toa))  && any(has_extrema, (
            # (aux_domain.θ_liq_ice[k - 1], aux_domain.θ_liq_ice[k], aux_domain.θ_liq_ice[k + 1]),
            # (aux_domain.q_tot[k - 1], aux_domain.q_tot[k], aux_domain.q_tot[k + 1]),
            # (aux_domain.q_liq[k - 1], aux_domain.q_liq[k], aux_domain.q_liq[k + 1]), # idk if this really matters, doesn't impact thermo/supersat
            # (aux_domain.q_ice[k - 1], aux_domain.q_ice[k], aux_domain.q_ice[k + 1]), # idk if this really matters, doesn't impact thermo/supersat
            # (aux_domain.T[k - 1], aux_domain.T[k], aux_domain.T[k + 1]),
            # (TD.vapor_specific_humidity(thermo_params, aux_domain.ts[k-1]) - aux_domain.q_vap_sat_liq[k - 1],  TD.vapor_specific_humidity(thermo_params, aux_domain.ts[k]) - aux_domain.q_vap_sat_liq[k],  TD.vapor_specific_humidity(thermo_params, aux_domain.ts[k+1]) - aux_domain.q_vap_sat_liq[k + 1]), # liq supersat (bc both q_vap_sat and T are prolly falling off pretty fast, we can have supersat peaks even w/o peaks in q, T, esp bc theta has even fewer peaks than T)
            # (TD.vapor_specific_humidity(thermo_params, aux_domain.ts[k-1]) - aux_domain.q_vap_sat_ice[k - 1],  TD.vapor_specific_humidity(thermo_params, aux_domain.ts[k]) - aux_domain.q_vap_sat_ice[k],  TD.vapor_specific_humidity(thermo_params, aux_domain.ts[k+1]) - aux_domain.q_vap_sat_ice[k + 1]), # ice supersat
            #
            (θ_liq_ice_km1, θ_liq_ice, θ_liq_ice_kp1),
            (q_tot_km1, q_tot, q_tot_kp1),
            (q_liq_km1, q_liq, q_liq_kp1),
            (q_ice_km1, q_ice, q_ice_kp1),
            (T_km1, T, T_kp1),
            (TD.vapor_specific_humidity(thermo_params, ts_km1) - q_vap_sat_liq_km1,  TD.vapor_specific_humidity(thermo_params, ts) - q_vap_sat_liq,  TD.vapor_specific_humidity(thermo_params, ts_kp1) - q_vap_sat_liq_kp1),
            (TD.vapor_specific_humidity(thermo_params, ts_km1) - q_vap_sat_ice_km1,  TD.vapor_specific_humidity(thermo_params, ts) - q_vap_sat_ice,  TD.vapor_specific_humidity(thermo_params, ts_kp1) - q_vap_sat_ice_kp1)
            #
            ))
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


            if (region isa EnvOrUp) || (region isa EnvRemainingDomain) # no cloak
                area = (region isa EnvRemainingDomain) ? aux_domain.a_en_remaining : aux_domain.area
                θ_liq_ice = aux_domain.θ_liq_ice
                q_tot = aux_domain.q_tot
                q_liq = aux_domain.q_liq
                q_ice = aux_domain.q_ice
            elseif (region isa CloakUp) # cloak up
                area = aux_domain.area_cloak_up
                θ_liq_ice = aux_domain.θ_liq_ice_cloak_up
                q_tot = aux_domain.q_tot_cloak_up
                q_liq = aux_domain.q_liq_cloak_up
                q_ice = aux_domain.q_ice_cloak_up
            elseif (region isa CloakDown) # cloak down
                area = aux_domain.area_cloak_down
                θ_liq_ice = aux_domain.θ_liq_ice_cloak_down
                q_tot = aux_domain.q_tot_cloak_down
                q_liq = aux_domain.q_liq_cloak_down
                q_ice = aux_domain.q_ice_cloak_down
            end

            θ_liq_ice_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (θ_liq_ice[scenario.bottom_extrap.ks[1]], θ_liq_ice[scenario.bottom_extrap.ks[2]]))
            θ_liq_ice_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (θ_liq_ice[scenario.top_extrap.ks[1]], θ_liq_ice[scenario.top_extrap.ks[2]]))
            q_tot_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (q_tot[scenario.bottom_extrap.ks[1]], q_tot[scenario.bottom_extrap.ks[2]]))
            q_tot_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (q_tot[scenario.top_extrap.ks[1]], q_tot[scenario.top_extrap.ks[2]]))
            q_liq_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (q_liq[scenario.bottom_extrap.ks[1]], q_liq[scenario.bottom_extrap.ks[2]]))
            q_liq_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (q_liq[scenario.top_extrap.ks[1]], q_liq[scenario.top_extrap.ks[2]]))
            q_ice_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (q_ice[scenario.bottom_extrap.ks[1]], q_ice[scenario.bottom_extrap.ks[2]]))
            q_ice_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (q_ice[scenario.top_extrap.ks[1]], q_ice[scenario.top_extrap.ks[2]]))

            p_c_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (p_c[scenario.bottom_extrap.ks[1]], p_c[scenario.bottom_extrap.ks[2]]))
            p_c_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (p_c[scenario.top_extrap.ks[1]], p_c[scenario.top_extrap.ks[2]]))


            p_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (aux_domain.p[scenario.bottom_extrap.ks[1]], aux_domain.p[scenario.bottom_extrap.ks[2]]))
            p_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (aux_domain.p[scenario.top_extrap.ks[1]], aux_domain.p[scenario.top_extrap.ks[2]]))

            # ρ_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (ρ[scenario.bottom_extrap.ks[1]], ρ[scenario.bottom_extrap.ks[2]]))
            # ρ_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (ρ[scenario.top_extrap.ks[1]], ρ[scenario.top_extrap.ks[2]]))

            area_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (area[scenario.bottom_extrap.ks[1]], area[scenario.bottom_extrap.ks[2]]))
            area_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (area[scenario.top_extrap.ks[1]], area[scenario.top_extrap.ks[2]]))

            dqvdt_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (dqvdt[scenario.bottom_extrap.ks[1]], dqvdt[scenario.bottom_extrap.ks[2]]))
            dqvdt_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (dqvdt[scenario.top_extrap.ks[1]], dqvdt[scenario.top_extrap.ks[2]]))

            dTdt_lo = positive_linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (dTdt[scenario.bottom_extrap.ks[1]], dTdt[scenario.bottom_extrap.ks[2]]))
            dTdt_hi = positive_linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (dTdt[scenario.top_extrap.ks[1]], dTdt[scenario.top_extrap.ks[2]]))


            # if θ or qt are 0, set them back to the k value
            if iszero(θ_liq_ice_lo)
                θ_liq_ice_lo = θ_liq_ice[k]
            end
            if iszero(θ_liq_ice_hi)
                θ_liq_ice_hi = θ_liq_ice[k]
            end
            if iszero(q_tot_lo)
                q_tot_lo = q_tot[k]
            end
            if iszero(q_tot_hi)
                q_tot_hi = q_tot[k]
            end


            # For w we do not extrapolate because we have it on the faces already
            w_lo = aux_domain_f.w[kf]
            w_hi = aux_domain_f.w[kf+1]



            ts_lo = thermo_state_pθq(param_set, p_c_lo, θ_liq_ice_lo, q_tot_lo, q_liq_lo, q_ice_lo) # these states actually have their own density, we operate assuming the background 
            ts_hi = thermo_state_pθq(param_set, p_c_hi, θ_liq_ice_hi, q_tot_hi, q_liq_hi, q_ice_hi)
            T_lo = TD.air_temperature(thermo_params, ts_lo)
            T_hi = TD.air_temperature(thermo_params, ts_hi)

            ρ_lo = TD.air_density(thermo_params, ts_lo)
            ρ_hi = TD.air_density(thermo_params, ts_hi)



        
            # These are extrapolated points so we need to recalculate...
            q_vap_sat_liq_lo = TD.q_vap_saturation_generic(thermo_params, T_lo, ρ_lo, TD.Liquid())
            q_vap_sat_ice_lo = TD.q_vap_saturation_generic(thermo_params, T_lo, ρ_lo, TD.Ice())
            q_vap_sat_liq_hi = TD.q_vap_saturation_generic(thermo_params, T_hi, ρ_hi, TD.Liquid())
            q_vap_sat_ice_hi = TD.q_vap_saturation_generic(thermo_params, T_hi, ρ_hi, TD.Ice())

            # for the NN we could consider extrapolating but if it's expnential maybe you want a geometric mean or something? idk. [[ now that we only apply at low resolution now, it's better... ]]
            q_lo = TD.PhasePartition(q_tot_lo, q_liq_lo, q_ice_lo)
            q_hi = TD.PhasePartition(q_tot_hi, q_liq_hi, q_ice_hi)
            τ_liq_lo, τ_ice_lo =  get_τs(param_set, microphys_params, nonequilibrium_moisture_scheme, q_lo, T_lo, p_lo, ρ_lo, w_lo)
            τ_liq_hi, τ_ice_hi =  get_τs(param_set, microphys_params, nonequilibrium_moisture_scheme, q_hi, T_hi, p_hi, ρ_hi, w_hi)
            # end

            mph_neq_lo = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, area_lo, ρ_lo, p_c_lo, T_lo, Δt + ε, ts_lo, w_lo, q_vap_sat_liq_lo, q_vap_sat_ice_lo, dqvdt_lo, dTdt_lo, τ_liq_lo, τ_ice_lo)
            mph_neq_hi = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, area_hi, ρ_hi, p_c_hi, T_hi, Δt + ε, ts_hi, w_hi, q_vap_sat_liq_hi, q_vap_sat_ice_hi, dqvdt_hi, dTdt_hi, τ_liq_hi, τ_ice_hi)

            #=
                There's an argument that we should weigh [mph_neq_lo, mph, mph_neq_hi] by the some fraction of [area_lo, area, area_hi]...  so that in the end multiplying by area gives a weighted answer
                I guess you could also argue that they're just local relative rates and do compose normally... but i think that fals apart as the gradients get sharper, e.g. right below an upddraft spike the are areally does matter...

                This scaling can backfire. if aux_domain.area[k] is near 0 but nonzero (common), and next to a large area, the implied tendency from this rescaling could be massive! We would need a broader framework built into the model where rates act bin-wide to make that work fully
                this would be most common in the updraft. Our best bet, then, is to ignore the part that comes from changes in area, and assume area is constant.

                The questions of if we should still bother with the density weighting still remains... for now we are not, even though the problem there is less severe.
            =#

            # Area and density weight. Essentially it's (ρ_hi*a_hi*q_hi * a_weight_hi + ρaq * a_weigh_mid + ρ_lo*a_lo*q_lo * a_weight_lo) / (ρ*a), but we did the other weights before..
            ql_tendency += scenario.prob * (mph_neq_lo.ql_tendency * area_weight.lo + mph_neq.ql_tendency * area_weight.mid + mph_neq_hi.ql_tendency * area_weight.hi) # area and density weight
            qi_tendency += scenario.prob * (mph_neq_lo.qi_tendency * area_weight.lo + mph_neq.qi_tendency * area_weight.mid + mph_neq_hi.qi_tendency * area_weight.hi) # area and density weight
            

            if any(!isfinite, (ql_tendency, qi_tendency))
                @error "Non-equilibrium moisture sources returned non-finite values: $mph_neq; from inputs; w = $(w[k]); zc = $zc; q_vap_sat_liq = $(aux_domain.q_vap_sat_liq[k]); q_vap_sat_ice = $(aux_domain.q_vap_sat_ice[k]); ρ_c = $(ρ_c[k]); aux_domain.area = $(aux_domain.area[k])"
                @info "parts were mph_neq_lo = $mph_neq_lo; mph_neq = $mph_neq; mph_neq_hi = $mph_neq_hi"
                @info "more parts were θ_liq_ice_lo = $θ_liq_ice_lo; θ_liq_ice_hi = $θ_liq_ice_hi; q_tot_lo = $q_tot_lo; q_tot_hi = $q_tot_hi; q_liq_lo = $q_liq_lo; q_liq_hi = $q_liq_hi; q_ice_lo = $q_ice_lo; q_ice_hi = $q_ice_hi; p_c_lo = $p_c_lo; p_c_hi = $p_c_hi; ρ_lo = $ρ_lo; ρ_hi = $ρ_hi; area_lo = $area_lo; area_hi = $area_hi; w_lo = $w_lo; w_hi = $w_hi; ts_lo = $ts_lo; ts_hi = $ts_hi; T_lo = $T_lo; T_hi = $T_hi; q_vap_sat_liq_lo = $q_vap_sat_liq_lo; q_vap_sat_ice_lo = $q_vap_sat_ice_lo; q_vap_sat_liq_hi = $q_vap_sat_liq_hi; q_vap_sat_ice_hi = $q_vap_sat_ice_hi"
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
    aux_domain::CC.Fields.Field, # FieldVector was wrong
    q::TD.PhasePartition, # alreayd calcluated, may as well pass it in
    p_c::CC.Fields.Field, # this is from gm so just pass In
    ;
    reweight_extrema_only::Bool = true, # should we pass this in or get it from param_set? maybe make it explicity so you can cache it slash it's more viislbe in the code
    region::AbstractDomain = Env,
    )

    FT = eltype(param_set)


    # no cloak
    if reweight_extrema_only && !((k == kc_surf) || (k == kc_toa))
        if (region isa EnvOrUpDomain) || (region isa EnvRemainingDomain) # no cloak
            q_tot_km1 = aux_domain.q_tot[k - 1]
            q_tot = aux_domain.q_tot[k]
            q_tot_kp1 = aux_domain.q_tot[k + 1]
            #
            θ_liq_ice_km1 = aux_domain.θ_liq_ice[k - 1]
            θ_liq_ice = aux_domain.θ_liq_ice[k]
            θ_liq_ice_kp1 = aux_domain.θ_liq_ice[k + 1]
            #
            T_km1 = aux_domain.T[k - 1]
            T = aux_domain.T[k]
            T_kp1 = aux_domain.T[k + 1]
        elseif (region isa CloakUp) # cloak up
            q_tot_km1 = aux_domain.q_tot_cloak_up[k - 1]
            q_tot = aux_domain.q_tot_cloak_up[k]
            q_tot_kp1 = aux_domain.q_tot_cloak_up[k + 1]
            #
            θ_liq_ice_km1 = aux_domain.θ_liq_ice_cloak_up[k - 1]
            θ_liq_ice = aux_domain.θ_liq_ice_cloak_up[k]
            θ_liq_ice_kp1 = aux_domain.θ_liq_ice_cloak_up[k + 1]
            #
            T_km1 = aux_domain.T_cloak_up[k - 1]
            T = aux_domain.T_cloak_up[k]
            T_kp1 = aux_domain.T_cloak_up[k + 1]
        elseif (region isa CloakDown) # cloak down
            q_tot_km1 = aux_domain.q_tot_cloak_down[k - 1]
            q_tot = aux_domain.q_tot_cloak_down[k]
            q_tot_kp1 = aux_domain.q_tot_cloak_down[k + 1]
            #
            θ_liq_ice_km1 = aux_domain.θ_liq_ice_cloak_down[k - 1]
            θ_liq_ice = aux_domain.θ_liq_ice_cloak_down[k]
            θ_liq_ice_kp1 = aux_domain.θ_liq_ice_cloak_down[k + 1]
            #
            T_km1 = aux_domain.T_cloak_down[k - 1]
            T = aux_domain.T_cloak_down[k]
            T_kp1 = aux_domain.T_cloak_down[k + 1]
        end
    end



    # So here, it's mostly the same but we need to do saturation adjustment and weight. These will never feed back to the prognostic state.

    # check for extrema in prognostic variables [θ_li, qt, ql, qi]
    # Ideally we'd do [T] bc supersat and thus q_vap_sat depend on T, but it's derived from the other variables so we can't do that.
    # we would still need to recalculate q_vap_sat_liq and q_vap_sat_ice, but we can do that in the extrapolation step.

    # we can't be withn one grid point of surf or top of atm
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)

    # reweight [if we reweight everywhere] OR [if we reweight extrema only, we're not at the surface or top of atm, and we have extrema in the prognostic variables]
    if (!reweight_extrema_only) || 
        (reweight_extrema_only &&  !((k == kc_surf) || (k == kc_toa))  && any(has_extrema, (
            # (aux_domain.θ_liq_ice[k - 1], aux_domain.θ_liq_ice[k], aux_domain.θ_liq_ice[k + 1]),
            # (aux_domain.q_tot[k - 1], aux_domain.q_tot[k], aux_domain.q_tot[k + 1]),
            # (aux_domain.T[k - 1], aux_domain.T[k], aux_domain.T[k + 1])),
            (θ_liq_ice_km1, θ_liq_ice, θ_liq_ice_kp1),
            (q_tot_km1, q_tot, q_tot_kp1),
            (T_km1, T, T_kp1)
            ),
        ) )

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

            if (region isa EnvOrUpDomain) || (region isa EnvRemainingDomain) # no cloak
                area = (region isa EnvRemainingDomain) ? aux_domain.a_en_remaining : aux_domain.area
                θ_liq_ice = aux_domain.θ_liq_ice
                q_tot = aux_domain.q_tot
            elseif (region isa CloakUpDomain) # cloak up
                θ_liq_ice = aux_domain.θ_liq_ice_cloak_up
                q_tot = aux_domain.q_tot_cloak_up
                area = aux_domain.area_cloak_up
            elseif (region isa CloakDownDomain) # cloak down
                θ_liq_ice = aux_domain.θ_liq_ice_cloak_down
                q_tot = aux_domain.q_tot_cloak_down
                area = aux_domain.area_cloak_down
            end


            θ_liq_ice_lo = linear_interpolate_extrapolate(zf_lo,  scenario.bottom_extrap.zs, (θ_liq_ice[scenario.bottom_extrap.ks[1]], θ_liq_ice[scenario.bottom_extrap.ks[2]]))
            θ_liq_ice_hi = linear_interpolate_extrapolate(zf_hi,  scenario.top_extrap.zs, (θ_liq_ice[scenario.top_extrap.ks[1]], θ_liq_ice[scenario.top_extrap.ks[2]]))
            q_tot_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (q_tot[scenario.bottom_extrap.ks[1]], q_tot[scenario.bottom_extrap.ks[2]]))
            q_tot_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (q_tot[scenario.top_extrap.ks[1]], q_tot[scenario.top_extrap.ks[2]]))

            p_c_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (p_c[scenario.bottom_extrap.ks[1]], p_c[scenario.bottom_extrap.ks[2]]))
            p_c_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (p_c[scenario.top_extrap.ks[1]], p_c[scenario.top_extrap.ks[2]]))



            area_lo = linear_interpolate_extrapolate(zf_lo, scenario.bottom_extrap.zs, (area[scenario.bottom_extrap.ks[1]], area[scenario.bottom_extrap.ks[2]]))
            area_hi = linear_interpolate_extrapolate(zf_hi, scenario.top_extrap.zs, (area[scenario.top_extrap.ks[1]], area[scenario.top_extrap.ks[2]]))


            ts_lo = thermo_state_pθq(param_set, p_c_lo, θ_liq_ice_lo, q_tot_lo)
            ts_hi = thermo_state_pθq(param_set, p_c_hi, θ_liq_ice_hi, q_tot_hi)

            # ρ_lo = TD.air_density(thermo_params, ts_lo)
            # ρ_hi = TD.air_density(thermo_params, ts_hi)


            q_lo = TD.PhasePartition(thermo_params, ts_lo) # we have to construct the phase partition based on q_tot in Eq
            q_hi = TD.PhasePartition(thermo_params, ts_hi) # we have to construct the phase partition based on q_tot in Eq

            
             #=
                There's an argument that we should weigh by the some fraction of [area_lo, area, area_hi]...  so that in the end multiplying by area gives a weighted answer
                I guess you could also argue that they're just local relative rates and do compose normally... but i think that fals apart as the gradients get sharper, e.g. right below an upddraft spike the are areally does matter...

                This scaling can backfire. if aux_domain.area[k] is near 0 but nonzero (common), and next to a large area, the implied tendency from this rescaling could be massive!
                this would be most common in the updraft. Our best bet, then, is to ignore the part that comes from changes in area, and assume area is constant.

                The question of if we should still bother with the density weighting still remains... for now we are not, even though the problem there is less severe.
            =#


            
            # Area and density weight. Essentially it's (ρ_hi*a_hi*q_hi * a_weight_hi + ρaq * a_weigh_mid + ρ_lo*a_lo*q_lo * a_weight_lo) / (ρ*a), but we did the other weights before..
            ql += scenario.prob * (q_lo.liq * area_weight.lo + q.liq * area_weight.mid + q_hi.liq * area_weight.hi)
            qi += scenario.prob * (q_lo.ice * area_weight.lo + q.ice * area_weight.mid + q_hi.ice * area_weight.hi)


            if any(!isfinite, (ql, qi))
                @error "Sat adjust moisture sources returned non-finite values: $q; from inputs; zc = $zc; aux_domain.area = $(aux_domain.area[k])"
                @info "parts were q_lo = $q_lo; q = $q; q_hi = $q_hi"
                @info "more parts were θ_liq_ice_lo = $θ_liq_ice_lo; θ_liq_ice_hi = $θ_liq_ice_hi; q_tot_lo = $q_tot_lo; q_tot_hi = $q_tot_hi; p_c_lo = $p_c_lo; p_c_hi = $p_c_hi; ts_lo = $ts_lo; ts_hi = $ts_hi"
                error("Equilibrium sat adjust returned non-finite values")
            end

        end
  

    else
        ql = q.liq
        qi = q.ice
    end

    return ql, qi

end
