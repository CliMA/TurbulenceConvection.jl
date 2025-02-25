#=
These are methods to take timescales, sources, etc and limit them for our timestep to ensure that we don't violate any physical constraints and crash
=#




"""
No limiter means just return the sources as they are
"""
function calculate_timestep_limited_sources(moisture_sources_limiter::NoMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState) where {FT}

    S_ql = (q_vap - q_eq.liq) / τ_liq  
    S_qi = (q_vap - q_eq.ice) / τ_ice 

    return calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, S_ql, S_qi)
end

function calculate_timestep_limited_sources(moisture_sources_limiter::NoMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT,  q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState, S_ql::FT, S_qi::FT) where {FT}
    
    if iszero(q.liq) && (S_ql < 0) # don't allow evap if there's already no liquid (edge case)
        S_ql = FT(0)
    end

    if iszero(q.ice) && (S_qi < 0)  # don't allow sublimation if there's already no ice (edge case)
        S_qi = FT(0)
    end

    return S_ql, S_qi
end


# ========================================================================================================================================================================================================================================================= #
"""
The basic limiter only does a quick check against the timestep to ensure we don't consume more of our source than we have available.
However, note it's checking q_vap, not supersaturation. So this is remains appropriate for say relax_to_equilibrium but wont't save you from progressing to subsaturation.
"""
function calculate_timestep_limited_sources(moisture_sources_limiter::BasicMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState) where {FT}
    S_ql = (q_vap - q_eq.liq) / τ_liq # | microphys_params.τ_cond_evap | CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, TD.PhasePartition(FT(0),q_vap,FT(0)))  
    S_qi = (q_vap - q_eq.ice) / τ_ice # -(source to vapor) = source to condensate
    return calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, S_ql, S_qi)
end

function calculate_timestep_limited_sources(moisture_sources_limiter::BasicMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT,  q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState, S_ql::FT, S_qi::FT) where {FT}

    # TODO - handle limiters elswhere (note, these break our allowance for liquid to ice compensation, no?)
    if S_ql >= FT(0)
        S_ql = min(S_ql, q_vap / Δt)
    else
        S_ql = -min(-S_ql, q.liq / Δt)
    end
    if S_qi >= FT(0)
        S_qi = min(S_qi, q_vap / Δt)
    else
        S_qi = -min(-S_qi, q.ice / Δt)
    end
    return S_ql, S_qi
end

# ========================================================================================================================================================================================================================================================= #

"""
The morrison milbrandt calculations are a bit more complex... essentially they asymptote to a steady state.
They are not a magic bullet then -- as Δt increases, these can become diluted. If other sources are present, and say if Δt is changing, the ratio of this to other sources will continually change.
This complicates balancing and limiters.

NOTE: We could try to make a version of this that does not cutoff but it would have WBF going off both directions forever... This may not be desirable. Also it's time evolution is dependent on the evolution of q_liq, q_ice, q_vap, so composing with other sources isn't as simple as adding the tendencies together. Just use with caution.
"""
function calculate_timestep_limited_sources(moisture_sources_limiter::MorrisonMilbrandt2015MoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState) where {FT}
    # this *shouldn't* need limiters the way it's defined but... lol we'll see...
    if iszero(Δt)
        # @error "bad things could happen, Δt is zero"
        return FT(0), FT(0)
    end
    S_ql, S_qi = morrison_milbrandt_2015_style(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts; emit_warnings = false)
    if !isfinite(S_ql) || !isfinite(S_qi)
        error("Morrison-Milbrandt 2015 source calculations failed. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts")
    end
    return S_ql, S_qi
end

function calculate_timestep_limited_sources(moisture_sources_limiter::MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState) where {FT}
    if iszero(Δt)
        # @error "bad things could happen, Δt is zero"
        return FT(0), FT(0)
    end
    S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts; emit_warnings = false)
    if !isfinite(S_ql) || !isfinite(S_qi)
        error("Morrison-Milbrandt 2015 (exponential part only) source calculations failed. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts")
    end
    return S_ql, S_qi
end


# ========================================================================================================================================================================================================================================================= #

"""
This is the most basic version of the supersaturation relaxation limiter. It just ensures that we don't consume more of our source than we have available, including not consuming more than the supersaturation we have.
The logic is more compilcated in WBF regimes and when long enough timesteps would cause regime transitions. Consider using the Morrison-Milbrandt source calculations to elide some of these issues.å
"""
function calculate_timestep_limited_sources(moisture_sources_limiter::StandardSupersaturationMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState) where {FT}
    S_ql = (q_vap - q_eq.liq) / τ_liq # | microphys_params.τ_cond_evap | CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, TD.PhasePartition(FT(0),q_vap,FT(0)))  
    S_qi = (q_vap - q_eq.ice) / τ_ice # -(source to vapor) = source to condensate

    # if T > TCP.thermodynamics_params(param_set).T_freeze # handled in other method now...
    #     S_qi = min(S_qi, FT(0)) # no ice growth above freezing.
    # end
    return calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, S_ql, S_qi,)
end

function t_out_of_x(x::FT, S::FT) where {FT} # helper func, how long does source S take to send x to 0
    if iszero(S) || iszero(x)
        return FT(Inf)
    end
    t = x / -S
    return ( t < FT(0) ) ? FT(Inf) : t
end

function calculate_timestep_limited_sources(moisture_sources_limiter::StandardSupersaturationMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT,  q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState, S_ql::FT, S_qi::FT) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)
    # Local sources, note these will always have the same sign as the source terms above since they just have a different Δt instead of τ
 
    # @error "This function is not temperature agnostic either... Above T_freeze = 273 when q_eq flips, you can get wrong results.... look into how to fix.... There's also som other bugs elsewhere... that's why we're using Morrison-Milbrandt type calculations instead, limiters are much easier..."
    # @warn("use Morrison-Milbrand exponential part only instead. In the WBF regime, we have no way of balancing the liquid and ice sources properly or limiting without using an analytic solution or substepping...")
    # e.g. if supersat over liquid is zero, and over ice is >0, what sources should this function return?  it's unclear... we could give the ice growth with no liquid loss and then let the liquid evaporate at the next timestep to make up for it..."

    below_freezing = (T < thermo_params.T_freeze)

    #=
    Eq_point means -S_ql = S_qi.
    This means (eq_point - q_eq.ice)/τ_ice = (q_eq.liq - eq_point)/τ_liq  ==> eq_point = (q_eq.ice * τ_liq + q_eq.liq * τ_ice)/(τ_liq + τ_ice)
    see https://www.wolframalpha.com/input?i2d=true&i=Divide%5B%5C%2840%29x+-+Subscript%5Bq%2Ci%5D%5C%2841%29%2Ck%5D%3DDivide%5B%5C%2840%29Subscript%5Bq%2CL%5D-x%5C%2841%29%2Cc%5D&assumption=%22UnitClash%22+-%3E+%7B%22L%22%2C+%7B%22Liters%22%2C+%22dflt%22%7D%7D&assumption=%7B%22C%22%2C+%22c%22%7D+-%3E+%7B%22Variable%22%2C+%22dflt%22%7D&assumption=%7B%22C%22%2C+%22L%22%7D+-%3E+%7B%22Variable%22%7D&assumption=%7B%22C%22%2C+%22q%22%7D+-%3E+%7B%22Variable%22%2C+%22dflt%22%7D
    =#
    # q_eq_point = (q_eq.ice * τ_liq + q_eq.liq * τ_ice)/(τ_liq + τ_ice) # Below-freezing WBF : at this point we should have S_ql = S_qi = 0, so we can just update the equilibrium point
    
    # safer floating point form
    q_eq_point = q_eq.ice * (τ_liq/(τ_liq + τ_ice)) + q_eq.liq * (τ_ice/(τ_liq + τ_ice))  # [Floating point safer form] Below-freezing WBF : at this point we should have S_ql = S_qi = 0, so we can just update the equilibrium point

    # Fix floating point errors that can give answers out of bounds -- notably e.g. q_eq.ice  * τ_liq / τ_liq cam actually shrink the answer! [using safer floting point form instead ]
    if below_freezing
        if q_eq_point < q_eq.ice
            q_eq_point = q_eq.ice
        elseif q_eq_point > q_eq.liq
            q_eq_point = q_eq.liq
        end
    else
        if q_eq_point < q_eq.liq
            q_eq_point = q_eq.liq
        elseif q_eq_point > q_eq.ice
            q_eq_point = q_eq.ice
        end
    end


    # iterative sol'n should work better and take care of all the WBF drama...
    Δt_left = Δt
    q_liq_now = q.liq
    q_ice_now = q.ice
    q_vap_now = q_vap
    # S_ql_now = S_ql
    # S_qi_now = S_qi
    S_ql_now = iszero(q.liq) ? max(S_ql, FT(0)) : S_ql # no evap if no liquid
    if below_freezing
        S_qi_now = iszero(q.ice) ? max(S_qi, FT(0)) : S_qi # no sublimation if no ice
    else
        S_qi_now = iszero(q.ice) ? FT(0) : min(S_qi, FT(0)) # no sublimation if no ice and no ice growth above freezing
    end

    # @debug "q_liq_now, $q_liq_now, q_ice_now = $q_ice_now, q_vap_now = $q_vap_now, S_ql_now = $S_ql_now, S_qi_now = $S_qi_now, q_eq_point = $q_eq_point"

    while Δt_left > FT(0)

    

        t_out_of_liq = t_out_of_x(q_liq_now, S_ql_now)
        t_out_of_ice = t_out_of_x(q_ice_now, S_qi_now)
        q_vap_avail_liq = q_vap_now - q_eq.liq # by definition these have the same sign as S_ql, S_qi (or both are 0). either way is safe for t_out_of_x
        q_vap_avail_ice = q_vap_now - q_eq.ice # by definition these have the same sign as S_ql, S_qi (or both are 0). either way is safe for t_out_of_x
        t_hit_liq_sat = t_out_of_x(q_vap_avail_liq, -(S_ql_now + S_qi_now))
        t_hit_ice_sat = t_out_of_x(q_vap_avail_ice, -(S_ql_now + S_qi_now))

        q_vap_avail_eq_point = q_vap_now - q_eq_point

        t_hit_eq_point = t_out_of_x(q_vap_avail_eq_point, -(S_ql_now + S_qi_now)) # can happen both above and below freezing, its just the point where gains of one phase are equal to losses of the other

        # @debug "q_vap_avail_liq = $q_vap_avail_liq; q_vap_avail_ice = $q_vap_avail_ice | (S_ql_now + S_qi_now) = $(S_ql_now + S_qi_now)"

        
        # findmin goes left to right, so this is also priority if there's a tie. Eq point gets priority over saturation points bc hitting it stops oscillations
        min_t, i_min_t = findmin([t_out_of_liq, t_out_of_ice, t_hit_eq_point, t_hit_liq_sat, t_hit_ice_sat, ]) # no need for out of vap, you get there you're cooked anyway...


        if min_t < Δt_left
            # take the step
            Δt_left -= min_t

            # This way is robust to floating point error (e.g. if min_t could be very small or 0 due to underflow/overflow which could lead to us being stuck)
            if i_min_t == 1 #  out of liq
                q_vap_now += q_liq_now
                q_liq_now = FT(0)
                q_ice_now += S_qi_now * min_t
                q_vap_now -= S_qi_now * min_t
            elseif i_min_t == 2 # out of ice
                q_liq_now += S_ql_now * min_t
                q_vap_now -= S_ql_now * min_t
                q_vap_now += q_ice_now
                q_ice_now = FT(0)
            elseif i_min_t == 3 # hit eq point
                q_liq_now += S_ql_now * min_t
                q_ice_now += S_qi_now * min_t
                q_vap_now = q_eq_point
            elseif i_min_t == 4 # hit liq sat
                q_liq_now += S_ql_now * min_t
                q_ice_now += S_qi_now * min_t
                q_vap_now = q_eq.liq
            else # i_min_t == 5 # hit ice sat
                q_liq_now += S_ql_now * min_t
                q_ice_now += S_qi_now * min_t
                q_vap_now = q_eq.ice
            end

            # # This way fails on floating point min_t = 0
            # q_liq_now += S_ql_now * min_t
            # q_ice_now += S_qi_now * min_t
            # q_vap_now -= (S_ql_now + S_qi_now) * min_t

            # - we should have no issues, for our t_out_of_x we are basically doing q + (q/-S)*S = q-q = 0 always even in floating point. Tiny negative numbers the model discards anyway
            # - however, for very small numbers, we could have issues (underflow -Inf). Those arent fixable though... setting minimum/maximum in tau should help avoid that ever happening...
            # if perform_checks # add a kw arg for this if you end up needing it... model does the same thing though genearlally and if it's a large error (like >> eps(FT)) we have a bug
            #     q_liq_now = max(q_liq_now, FT(0))
            #     q_ice_now = max(q_ice_now, FT(0))
            # end


            # update rates
            # - (you could do a different eq based on which one you hit and eliminate floating point risks, but you then also need to check for ties...
            # - this is the most versatile way with no tie checks... (easiest for developer maintenance))
            
            S_ql_now = iszero(q_liq_now) ? max((q_vap_now - q_eq.liq) / τ_liq, FT(0)) : ((q_vap_now - q_eq.liq) / τ_liq) # no evaporation if q_liq is 0
            if below_freezing
                S_qi_now = iszero(q_ice_now) ? max((q_vap_now - q_eq.ice) / τ_ice, FT(0))  : ((q_vap_now - q_eq.ice) / τ_ice) # no sublimation if q_ice is 0
            else
                S_qi_now = iszero(q_ice_now) ? FT(0) : min((q_vap_now - q_eq.ice) / τ_ice, FT(0)) # no sublimation if q_ice is 0, no growth if we're above freezing
            end

            if i_min_t == 3 # hit eq point
                if q_eq_point == q_eq.liq # eq point is too close to liq to calculate liq tendency, use -S_qi
                    # @debug "eq point is too close to liq to calculate liq tendency, got S_ql_now = $S_ql_now, S_qi_now = $S_qi_now"
                    S_ql_now = -S_qi_now
                elseif q_eq_point == q_eq.ice # eq point is too close to ice to calculate ice tendency, use -S_ql
                    # @debug "eq point is too close to ice to calculate ice tendency, got S_ql_now = $S_ql_now, S_qi_now = $S_qi_now"
                    S_qi_now = -S_ql_now # eq point is too close to ice to calculate ice tendency, use -S_ql
                else
                    # @debug "eq point is not too close to either phase, got S_ql_now = $S_ql_now, S_qi_now = $S_qi_now"
                    # This just means we're so close that the floating-point math becomes hard.... the real solution is between this q_eq_point and next/prevfloat(q_vap_now) so we can't represent it
                    # In principle, q_eq_point is incorrect, it's not immediately clear which way to move for correctness, but in principle it's towards the faster growing phase, meaning the slower rate is probably better
                    # S_fixed = (abs(S_ql_now) + abs(S_qi_now)) / 2
                    S_fixed = min(abs(S_ql_now), abs(S_qi_now))
                    S_ql_now = sign(S_ql_now) * S_fixed
                    S_qi_now = sign(S_qi_now) * S_fixed
                    # @assert iszero(S_ql_now + S_qi_now) "Didn't hit eq point properly, got S_ql_now = $S_ql_now, S_qi_now = $S_qi_now from q_vap_now = $q_vap_now; q_eq_point = $q_eq_point; q_eq_liq = $(q_eq.liq); q_eq_ice = $(q_eq.ice); τ_liq = $τ_liq; τ_ice = $τ_ice"
                end
            end


        else 
            # finalize
            q_liq_now += S_ql_now * Δt_left
            q_ice_now += S_qi_now * Δt_left
            Δt_left = FT(0)
            # q_vap_now -= (S_ql + S_qi) * Δt_left # we don't need to update this, we don't iuse it anymore
            break
        end
    end

    S_ql = (q_liq_now - q.liq) / Δt
    S_qi = (q_ice_now - q.ice) / Δt

    return  S_ql, S_qi
end


#= Deprecated because:
- It's temperature agnostic, doesn't handle above and below freezing properly.
- For WBF the scaling is hard. We often just are scaling things down uniformly, which while okay-ish in super/subsaturated settings, is bad in WBF.
    Often, we'd just give preference to ice first, but if say liq is faster, we're moving in the wrong direction. scaling things down helps but it's hard to reason about.
- It's complex and still had some bugs I couldn't figure out


=#
# function calculate_timestep_limited_sources(moisture_sources_limiter::StandardSupersaturationMoistureSourcesLimiter, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT,  q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState, S_ql::FT, S_qi::FT) where {FT}
#     thermo_params = TCP.thermodynamics_params(param_set)
#     # Local sources, note these will always have the same sign as the source terms above since they just have a different Δt instead of τ
 

#     # Limiters, where we can rely on the fact that Q_vl,Q_vi are always the same sign as S_ql,S_qi
#     Q_vl = (q_vap - q_eq.liq) / Δt # vapor available for liquid growth in this timestep
#     Q_vi = (q_vap - q_eq.ice) / Δt # vapor available for ice growth in this timestep (if positive, this is larger?)

#     @info( "current", S_ql, S_qi, "--",Q_vl, Q_vi, "---", q.liq / Δt, q.ice / Δt, q_vap/Δt)
#     if (S_ql > 0) && (S_qi > 0) # both growing (this is the hardest with the different saturation vapor pressures) - we'll partition based on respective source magnitudes limited by Q_vi which is larger...
#         S_tot = S_ql + S_qi
#         # create liquid and ice in proportion to their source magnitudes and see if we reach liquid saturation

#         if Q_vl < Q_vi
#             if S_tot > Q_vl # this handles S_ql whether or not it's larger than Q_vl while simultaneously handling the first part of S_qi
#                 S_ql_Qvl = S_ql * Q_vl / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
#                 S_qi_Qvl = S_qi * Q_vl / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

#                 S_ql = S_ql_Qvl
#                 S_qi_addit = min(S_qi - S_qi_Qvl, Q_vi - Q_vl)  # we should still be able to consume some more ice as long as there was more room in the original ice source, (should be up to Q_vi - S_ql_Qvl which is total growth potential minus current accounted for growth)
#                 S_qi = S_qi_Qvl + S_qi_addit # ice growth limited by vapor supersat over liquid and simultaneous liquid growth
#             end # if not > Q_vl, it's def not greater than Q_vi so we don't need to do anything else I think
#         else
#             if S_tot > Q_vi # this handles S_ql whether or not it's larger than Q_vl while simultaneously handling the first part of S_qi
#                 S_ql_Qvi = S_ql * Q_vi / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
#                 S_qi_Qvi = S_qi * Q_vi / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

#                 S_qi = S_qi_Qvi
#                 S_ql_addit = min(S_ql - S_ql_Qvi, Q_vl - Q_vi)  # we should still be able to consume some more ice as long as there was more room in the original ice source, (should be up to Q_vi - S_ql_Qvl which is total growth potential minus current accounted for growth)
#                 S_ql = S_qi_Qvi + S_ql_addit # ice growth limited by vapor supersat over liquid and simultaneous liquid growth
#             end # if not > Q_vl, it's def not greater than Q_vi so we don't need to do anything else I think
#         end

#     elseif (S_ql > 0) && (S_qi < 0) # liquid growth but ice depletion (this can't happen I thnk because we'd be supersaturated over ice if we are over liquid)
#         if T < thermo_params.T_freeze
#             error(
#                 "This should not happen (though I guess it can above freezing but that will be handled in limiters later?)",
#             ) # This is just a check that all is working properly, remove later
#         else
#             # ice sublimates, liquid grows.
#             # hard to do bc the total source cannot pass supersaturation:
#             # - liquid growth limited by vapor supersat over liquid and simultaneous ice sublimation
#             # - ice sublimation limited by vapor supersat over ice and boosted by liquid growth

#             # Note : This isn't quite right, it doesnt correctly account for ice evap being limited by sub...

#             S_qi = -min(-S_qi, q.ice / Δt)  # don't let ice deplete more than exists or can evaporate before reaching saturation
#             S_ql = min(S_ql, Q_vl - S_qi)  # don't let liquid grow more than the vapor can support plus the ice that's sublimating (-S_qi which is negative)
#         end

#     elseif (S_ql < 0) && (S_qi > 0) # ice growth but liquid depletion ( WBF)

#         # we go equally until until either liquid or ice reaches saturation, then the remainder is a balance at that saturation
#         S_ql = -min(-S_ql, q.liq / Δt)  # don't let liquid deplete more than exists (smaller absolute value)
#         S_tot = S_ql + S_qi
#         q_vap_after = q_vap + S_tot * Δt

#         # Something is wrong in this block, we're getting evap sources that exceed Q_vi when there's no liquid...
#         if Q_vl < Q_vi
#             if q_eq.ice <= q_vap_after <= q_eq.liq # we stay between the equilibrium saturation range
#             # all good
#             elseif q_eq.ice > q_vap_after # we can't evaporate enough liquid to grow all the ice it wants to grow, so we scale both down so that we land on the ice saturation line
#                 scaling = -Q_vi / (S_tot * Δt) # should be a positive number < 1 (S_tot should be negative, and  q_eq.ice - q_vap should be negative)
#                 S_qi = S_qi * scaling
#                 S_ql = S_ql * scaling
#             elseif q_vap_after > q_eq.liq # we can't grow enough ice to support the liquid we want to evaporate
#                 scaling = -Q_vl / (S_tot * Δt) # should be a positive number < 1 (S_tot should be positive, and q_eq.liq - q_vap should be positive)
#                 S_qi *= scaling
#                 S_ql *= scaling
#             end
#         else
#             if q_eq.liq <= q_vap_after <= q_eq.ice # we stay between the equilibrium saturation range
#             # all good
#             elseif q_eq.liq > q_vap_after # we can't evaporate enough liquid to grow all the ice it wants to grow, so we scale both down so that we land on the ice saturation line
#                 scaling = (q_eq.liq - q_vap) / (S_tot * Δt) # should be a positive number < 1 (S_tot should be negative, and  q_eq.ice - q_vap should be negative)
#                 S_qi = S_qi * scaling
#                 S_ql = S_ql * scaling
#             elseif q_vap_after > q_eq.ice # we can't grow enough ice to support the liquid we want to evaporate
#                 scaling = (q_eq.ice - q_vap) / (S_tot * Δt) # should be a positive number < 1 (S_tot should be positive, and q_eq.liq - q_vap should be positive)
#                 S_qi *= scaling
#                 S_ql *= scaling
#             end
#         end

#     elseif (S_ql < 0) && (S_qi < 0) # both are depleting (limit by condensate amount and enough to return to saturation)
#         S_ql = -min(-S_ql, q.liq / Δt)  # don't let liquid deplete more than exists (smaller absolute value)
#         S_qi = -min(-S_qi, q.ice / Δt)  # don't let ice deplete more than exists (smaller absolute value)

#         S_ql = -min(-S_ql, -Q_vl)  # don't let liquid deplete more than the vapor can support (smaller absolute value)
#         S_qi = -min(-S_qi, -Q_vi)  # don't let ice deplete more than the vapor can support (smaller absolute value)

#         S_tot = S_ql + S_qi
#         # consume liquid and ice in proportion to their source magnitudes and see if we reach ice saturation

#         if Q_vl < Q_vi # should mean we're below freezing and can evaporate more liquid than we can ice 
#             if S_tot < Q_vi # we would create more vapor than the air wants
#                 # go combined until reaching ice saturation, then remainder is liquid evaporation only... (with no re-condensation of ice)
#                 S_ql_Qvi = S_ql * Q_vi / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
#                 S_qi_Qvi = S_qi * Q_vi / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

#                 S_qi = S_qi_Qvi # ice until ice reaching ice saturation vapor pressure... (ice T > freezing, ice is actually more )
#                 S_ql_addit = max(S_ql - S_ql_Qvi, Q_vl - Q_vi)  # we should be able to consume some more liquid up until we exhaust the original liquid source (S_ql - S_ql_Qvi), or reach the liquid saturation vapor pressure (Q_vl - Q_vi), whichever is smaller (max of 2 negative numbers)
#                 S_ql = S_ql_Qvi + S_ql_addit # liquid growth limited by vapor supersat over liquid and simultaneous ice growth

#                 if S_ql_addit > 0
#                     error("this shouldnt happen lol, should be max of two negative numbers")
#                 end

#             end # if not < Q_vi, we dont even reach ice saturation so we don't need to do anything else
#         else # this is technically weird bc all the ice should evaporate anyway when we're above freezing, but we'll handle anyway... (actually it's slightly different due to difference between thermo_params.T_triple and thermo_params.T_freeze which is 0.01 K but good to handle properly...)
#             if S_tot < Q_vl # we would create more vapor than the air wants
#                 # go combined until reaching ice saturation, then remainder is liquid evaporation only... (with no re-condensation of ice)
#                 S_ql_Qvl = S_ql * Q_vl / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
#                 S_qi_Qvl = S_qi * Q_vl / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

#                 S_ql = S_ql_Qvl # ice until ice reaching ice saturation vapor pressure... (ice T > freezing, ice is actually more )
#                 S_qi_addit = max(S_qi - S_qi_Qvl, Q_vi - Q_vl)  # we should be able to consume some more liquid up until we exhaust the original liquid source (S_ql - S_ql_Qvi), or reach the liquid saturation vapor pressure (Q_vl - Q_vi), whichever is smaller (max of 2 negative numbers)
#                 S_qi = S_qi_Qvl + S_qi_addit # liquid growth limited by vapor supersat over liquid and simultaneous ice growth

#                 if S_qi_addit > 0
#                     error("this shouldnt happen lol, should be max of two negative numbers")
#                 end

#             end # 

#         end

#     end # otherwise we have the negative limiters below for sublimation, evaporation... if both are neg that's sufficient....

#     @info( "output", S_ql, S_qi, Δt)

#     return S_ql, S_qi
# end


