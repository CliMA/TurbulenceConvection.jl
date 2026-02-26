"""
# Milbrandt's equations are in mixing ratio not specific humidity, to convert we need to use

# q_spe = q_mix / (1 + q_mix)   [ in reality this is more like q_spe_c = q_mix_c / (1 + q_tot_mix) , bc 1+q_tot = total ratio all water / dry air...]
# q_mix = q_spe / (1 - q_spe)   [ in reality this is more like q_mix_c = q_spe_c / (1 - q_tot_spe) , bc 1-q_tot = dry mixing ratio...]

# [NOTE: q_tot appears instead of just q because we have more than one water species so we need to be careful...]

# note then d/dy (q_spe) = d/dy(q_mix / (1 + q_mix)) = d/dy(q_mix) / (1 + q_mix)^2  [ really it's more like d/dy(q_spe_c) = d/dy(q_mix_c / (1 + q_tot_mix)) and q_tot_mix should be a constant under phase transformation... ]
# note then d/dy (q_mix) = d/dy(q_spe / (1 - q_spe)) = d/dy(q_spe) / (1 - q_spe)^2  [ really it's more like d/dy(q_mix_c) = d/dy(q_spe_c / (1 - q_tot_spe)) and q_tot_spe should be a constant under phase transformation... ]

NOTE: These mixing ratio, specific humidity, and their derivatives are essentially identical at earthlike conditions...


Because dqvdt can be non-zero, we are never done until we reach Δt, so we remove the early exits.

Note this negelects any change in q_sl - q_si over the timestep. it accounts for those hanges through Γ, but we only calculate one dδ, and assume the liquid and ice change the same amount.
This is in line w/ neglecting T changes, unlike the other regular MM2015 that finds the new T at each step.

"""

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #



function do_standard_fallback(milestone_t::FT, milestone::MilestoneType, time_tolerance::FT, S_ql::FT, S_qi::FT, q_liq::FT, q_ice::FT, δ_eq::FT, δi_eq::FT, dδdt_no_S::FT, Γ_l::FT, Γ_i::FT,
    regime::AbstractSaturationRegime, param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState; opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}()
    )::Tuple{FT, FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    # @debug "do_standard_fallback: milestone_t = $milestone_t; milestone = $milestone; time_tolerance = $time_tolerance; S_ql = $S_ql; S_qi = $S_qi; q_liq = $q_liq; q_ice = $q_ice; δ_eq = $δ_eq; δi_eq = $δi_eq; dδdt_no_S = $dδdt_no_S; Γ_l = $Γ_l; Γ_i = $Γ_i; regime = $regime; area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q_eq = $q_eq; Δt = $Δt"


    # old_δ_0 = δ_0
    # old_δ_0i = δ_0i

    # do step
    dt = min(milestone_t, Δt)
    q_liq, q_ice, δ_0, δ_0i, new_regime_enum_type = step(regime, StandardSupersaturationMoistureSourcesLimiter(), dt, q_liq, q_ice, δ_0, δ_0i, δ_eq, δi_eq, q_eq, S_ql, S_qi, ((milestone_t < Δt) ? milestone : NotAtSupersaturationMilestone); dδdt_no_S=dδdt_no_S, Γ_l=Γ_l, Γ_i=Γ_i) # if milestone_t < Δt then we do the step, otherwise we don't (we just return the current state)
    new_regime_type = get_saturation_regime_type(Val(new_regime_enum_type)) # 1 allocation

    new_regime = new_regime_type{q_liq>0, q_ice>0, is_below_freezing(regime)}() # reconstruct regime with new q_liq and q_ice
    Δt_left = Δt - dt


    # @debug("status: regime = $regime; milestone = $milestone; is_below_freezing(regime) = $(is_below_freezing(regime)); q_liq = $q_liq; q_ice = $q_ice; q_eq = $q_eq; S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; Δt_left = $Δt_left; dt = $dt; Γ_l = $Γ_l; Γ_i = $Γ_i; dqvdt = $dqvdt; dTdt = $dTdt; T = $T; p = $p; area = $area; ρ = $ρ; τ_liq = $τ_liq; τ_ice = $τ_ice; dδdt_no_S = $dδdt_no_S")

    if (milestone == AtSupersaturationStationaryPointMilestone) && (Δt_left > 0) # hit eq point, not a real milestone for morrison milbrandt (only for standard) so we still want to continue on to the next milestone
        # do this again
        milestone_t, milestone, S_ql_addit, S_qi_addit, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, is_below_freezing(regime), τ_liq, τ_ice; dδdt_no_S = dδdt_no_S, Γ_l=Γ_l, Γ_i=Γ_i, at_δ_eq_point = true) # can we use the sources we already have here?
        # @debug "milestone_t = $milestone_t; milestone = $milestone; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; δ_eq=$δ_eq; δi_eq=$δi_eq"
        
        # if (milestone_t < time_tolerance) 
        #=
            If we do not do this step (regardless of how long it is), and, to floating point error, we are just off the saturation point, we cannot guarantee that the next call will not immediately send us back here.
            Alternatively, our calls to calculate_next_standard_milestone_time() outside of this function, should not permit equilibrium points.

            We have opted for the latter sol'n, meaning we can't actually enter this block at all!

            But, if you ever allow falling back to i_min_t == 3 (at_δ_eq_point = true), then you will need to do this step (regardless of how long it is!) to ensure you don't return conrol to MM2015_EPA at the equilibrium point, but slightly off due to numerical error.

            UPDATE :: Simply ignoring the eq point was a mistake -- you can easily reach a situation where two MM2015_EPA calls bounce back and forth between liq and ice sat, losing some liq, forming some ice, losing some liq, forming some ice, and this can go on forever...
                        We used to even disallow such transitions, becauase we assumed that transitions are unidirectional in Standard()! (they're allowed in MM2015_EPA() under the idea that they shouldn't happen frequently, maybe one foray into supersat that is eventually reversed by dTdt or something, but not repeated oscillation)
                        So, instead, when faced with Standard() fallback hitting the eq point, no matter how long the step, we will take it.

                        If this is not satisfactory to you, you could implement always falling back to BigFloat when the milestone_t is too short. This would very likely work, but would be more expensive (BigFloats seem to be like 100-200x more expensive in our usage)
                        It is also not bulletproof if the time is very short, whereas Standard() is bulletproof down to floatmin() [[ I have seen BigFloat fail as high as 1e-20 ]]
        =#
        
        # do step
        dt_here = min(milestone_t - dt, Δt_left) # go up to the milestone time, but not past the remaining time
        Δt_left -= dt_here
        q_liq, q_ice, δ_0, δ_0i, new_regime_enum_type = step(new_regime, StandardSupersaturationMoistureSourcesLimiter(), dt_here, q_liq, q_ice, δ_0, δ_0i, δ_eq, δi_eq, q_eq, S_ql_addit, S_qi_addit, ((milestone_t < Δt_left) ? milestone : NotAtSupersaturationMilestone); at_δ_eq_point = true, dδdt_no_S=dδdt_no_S, Γ_l=Γ_l, Γ_i=Γ_i) # if milestone_t < Δt_left then we do the step, otherwise we don't (we just return the current state)
        new_regime_type = get_saturation_regime_type(Val(new_regime_enum_type))
        new_regime = new_regime_type{q_liq>0, q_ice>0, is_below_freezing(regime)}() # reconstruct regime with new q_liq and q_ice
        
        S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, dt, S_ql_addit, S_qi_addit, dt_here, dt + dt_here) # rescale to the timestep
        dt += dt_here

        # regime = get_saturation_regime(δ_0, δ_0i, q_liq, q_ice, is_below_freezing(regime)) #
        # @debug("status: regime = $regime; milestone = $milestone; is_below_freezing(regime) = $(is_below_freezing(regime)); q_liq = $q_liq; q_ice = $q_ice; q_eq = $q_eq; S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; Δt_left = $Δt_left; dt = $dt; Γ_l = $Γ_l; Γ_i = $Γ_i; dqvdt = $dqvdt; dTdt = $dTdt; T = $T; p = $p; area = $area; ρ = $ρ; τ_liq = $τ_liq; τ_ice = $τ_ice; dδdt_no_S = $dδdt_no_S")
        # regime = add_regime_parameters(get_new_saturation_regime_type_from_milestone(milestone, regime, δ_0, δ_0i ), q_liq, q_ice, is_below_freezing(regime)) # add the parameters to the regime [[ i think this way is safer as it assures transitions ]]
        # end

        !(milestone == AtSupersaturationStationaryPointMilestone) || error("milestone should not be AtSupersaturationStationaryPoint() here, i dont think you should be able to get AtSupersaturationStationaryPoint() twice in a row")
    end


    if iszero(Δt_left)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    else



        # @debug("status: regime = $regime; milestone = $milestone; is_below_freezing(regime) = $(is_below_freezing(regime)); q_liq = $q_liq; q_ice = $q_ice; q_eq = $q_eq; S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; Δt_left = $Δt_left; dt = $dt; Γ_l = $Γ_l; Γ_i = $Γ_i; dqvdt = $dqvdt; dTdt = $dTdt; T = $T; p = $p; area = $area; ρ = $ρ; τ_liq = $τ_liq; τ_ice = $τ_ice; dδdt_no_S = $dδdt_no_S")
        # regime_type = get_new_saturation_regime_type_from_milestone(milestone, regime, old_δ_0, old_δ_0i) # get the new regime type from the milestone


        new_q = TD.PhasePartition(q.tot, q_liq, q_ice) # don't actually let q_vap change...
        if TD.vapor_specific_humidity(new_q) < FT(0)
            # we shouldn't be able to get below 0 unless the saturatoin point is wrong I believe... it's possible q_vap_sat is at 0 I guess if the model is crashing.
            error("Got <0 vapor specific humidity, new_q_vap = $(TD.vapor_specific_humidity(new_q)); new_q = $new_q; old q=$q; q_eq = $q_eq; q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i; Δt_left = $Δt_left; milestone = $milestone; milestone_t = $milestone_t")
        end


        # @debug "new_regime = $new_regime; milestone_t = $milestone_t; Δt_left = $Δt_left; milestone = $milestone; S_ql = $S_ql; S_qi = $S_qi; q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i"
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))::Tuple{FT, FT}
        S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, dt, S_ql_addit, S_qi_addit, Δt_left, Δt)


        # @debug "S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; S_ql = $S_ql; S_qi = $S_qi; milestone_t = $milestone_t; Δt_left = $Δt_left; Δt = $Δt; q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i"
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ==================================================================================================================================================================================================================================================================================== #





function get_params_and_go_to_mixing_ratio_exponential_part_only(
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    ) where {FT}
    (; use_fix) = opts

    thermo_params = TCP.thermodynamics_params(param_set) # currently repeated in both places, pare down later
    g = TCP.grav(param_set) # acceleration of gravity
    L_i = TD.latent_heat_sublim(thermo_params, ts) # Latent heat for ice (L_i)
    L_l = TD.latent_heat_vapor(thermo_params, ts)  # Latent heat for water
    c_p = TD.TP.cp_d(thermo_params) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?

    # We have chosen not to bother with mixing ratios in this func... It shouldn't be neceessary.. it's just a scaling and qt is constant... when used as a fallback in morrison_milbrandt_2015_style() we just adjust accordingly....
    
    q_vap = TD.vapor_specific_humidity(q) # ensure it's not in mixing ratio units... [ note this  can cause a problem if q_vap doesn't equal the q_vap passed in due to floating point errors... ]
    q_liq = q.liq
    q_ice = q.ice
    q_sl = q_eq.liq
    q_si = q_eq.ice

    # CAN WE GET AWAY WITH NOT USING MIXING RATIO ANYWHERE? WE HAVE return_mixing_ratio = true EVERYWHERE! It's all just a scaling factor of 1/(1+q_tot) or 1/(1-q_tot)... don't think it changes the equations at all...
    # q_vap = TD.shum_to_mixing_ratio(TD.vapor_specific_humidity(q), q.tot) # convert to mixing ratio (use this version instead of q_vap bc q_vap passed from other functions into this fcn can be in mix form already than )
    # q_liq = TD.shum_to_mixing_ratio(q.liq, q.tot)
    # q_ice = TD.shum_to_mixing_ratio(q.ice, q.tot)
    # q_sl = TD.shum_to_mixing_ratio(q_eq.liq, q.tot) # q_eq is the equilibrium mixing ratio, q_sl is the saturation mixing ratio, q_eq we made to contain only saturation vapor pressure values...
    # q_si = TD.shum_to_mixing_ratio(q_eq.ice, q.tot)


    # T_freeze = TCP.T_freeze(param_set)
    T_freeze = TCP.T_triple(param_set) # In the code this is actually where the saturation vapor pressures are equal... I think it's wrong though



    # ==  We pass δ_0, δ_0i directly through each function now. == #
    # δ_0 = q_vap - q_sl # supersaturation over liquid
    # # δ_0i = δ_0 + (q_sl - q_si) # supersaturation over ice
    # δ_0i = q_vap - q_si # =  δ_0 + (q_sl - q_si) but fewer operations, less floating point precise ? # supersaturation over ice
    # δ_0 = limit_δ(δ_0, q_vap, q_liq, q_ice) # because we're always heading towards WBF, if we're `close enough` we'll just call it WBF to avoid misses.
    # δ_0i = limit_δ(δ_0i, q_vap, q_liq, q_ice) # because we're always heading towards WBF, if we're `close enough` we'll just call it WBF to avoid misses.
    # δ_0 and δ_0i will be wrong unless Γ_l, Γ_i are 1 as they change the conversion rate between cloud liquid and supersaturation and we do not track supersaturation

    # δ_0 = FT(NaN)
    # δ_0i = FT(NaN)


    # ======= Get Gammas ======= #
    # Γ_l = FT(1)
    # Γ_i = FT(1) # Eqn C3
    # dqsl_dT = FT(0)
    # dqsi_dT = FT(0) # Eqn C3

    # we don't let T update, so in the substeps applying the Gammas would change the asymptotes, but basically we put a ratio on the conversion of δ to q...
    # So if we can't figure out the timing part, we can at least keep the Γ_l and Γ_i just for the ratio in vapor to liquid and vapor to ice conversions...


    dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(1), T, q_eq.liq) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix... is it a good enough approx?
    dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(0), T, q_eq.ice) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...
    # dqsl_dT /= (1 - q_eq.liq)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.liq = q_sl at T
    # dqsi_dT /= (1 - q_eq.ice)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.ice = q_si at T
    # dqsl_dT = TD.shum_to_mixing_ratio(dqsl_dT, q.tot) # convert to mixing ratio (testinggg)
    # dqsi_dT = TD.shum_to_mixing_ratio(dqsi_dT, q.tot) # convert to mixing ratio (testinggg)

    # Balanced at large timestep, but needs limiters (otherwise, may consume more liq/ice than exists... but handles WBF directions fine and most stuff for short timesteps
    if use_fix # setting these makes the asymptote ok, not sure if it breaks other things....
        # i think their equation balances liq evap <-> ice growth, but doesn't limit it properly over the timestep, technically q_sl-q_si addition/subtraction should decay as you approach saturation...
        L_i = L_l
        dqsi_dT = dqsl_dT
        # Does evaporating water releasing more latent heat than forming ice takeup risk evaporating a cloud completely?
    end

    Γ_l = 1 + (L_l / c_p) * dqsl_dT  # Eqn C3
    Γ_i = 1 + (L_i / c_p) * dqsi_dT  # Eqn C3
    # ============================ #

    e_sl = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())


    return (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i, dqvdt, dTdt)
end

# τ_func_exponential_part_only(τ_liq::FT, τ_ice::FT) where {FT} = 1 / (1 / τ_liq +  1 / τ_ice) # Eqn C2
function τ_func_exponential_part_only(τ_liq::FT, τ_ice::FT, L_i::FT, c_p::FT, dqsl_dT::FT, Γ_i::FT ) where {FT}
    # τ = min( 1 / (1 / τ_liq +  1 / τ_ice), τ_liq, τ_ice) # Eqn C2 [ Better for Floating point safety... ]
    τ = min(τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i), τ_liq, τ_ice) # Eqn C2 [ min better for floating point safety....]
    return iszero(τ) ? min(τ_liq, τ_ice) : τ # if τ is 0, we'll just use the minimum of the two... Probably we faced an underflow...
end

const τ_func_EPA = τ_func_exponential_part_only

# THESE ARE WRITTEN ASSUMING WE AREN'T USING MIXING RATIO ANYWHERE, IS THAT BAD? OTHERWISE WE CAN JUST USE THE REGULAR Q_NEW FUNCTIONS FROM morrison_milbrandt_2015_style.jl
# Eventually we need to deprecate these, we really only care about ql, qi, and δ_0,δ_0i. Now we track δ_0, δ_0i directly, so we really don't need to create a phase partition or calcluate a q_vap at all...
q_new_exponential_part_only(q::TD.PhasePartition, S_ql_Δt::FT, S_qi_Δt::FT, dqvdt_Δt::FT) where {FT} = TD.PhasePartition(q.tot + dqvdt_Δt, q.liq + S_ql_Δt, q.ice + S_qi_Δt) # Eqn C3
q_new_exponential_part_only(q::TD.PhasePartition, S_ql::FT, S_qi::FT, Δt::FT, dqvdt::FT) where {FT} = TD.PhasePartition(q.tot + dqvdt * Δt, q.liq + S_ql * Δt, q.ice + S_qi * Δt) # Eqn C3
q_new_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql_Δt::FT, S_qi_Δt::FT, dqvdt_Δt::FT) where {FT} = TD.PhasePartition(q.tot + dqvdt_Δt, q_liq + S_ql_Δt, q_ice + S_qi_Δt) # Eqn C3
q_new_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql::FT, S_qi::FT, Δt::FT, dqvdt::FT) where {FT} = TD.PhasePartition(q.tot + dqvdt * Δt, q_liq + S_ql * Δt, q_ice + S_qi * Δt) # Eqn C3
const q_new_EPA = q_new_exponential_part_only


function morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q::TD.PhasePartition)
    return new_q
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, S_ql_Δt::FT, S_qi_Δt::FT, dqvdt_Δt::FT) where {FT}
    new_q = q_new_exponential_part_only(q, S_ql_Δt, S_qi_Δt, dqvdt_Δt)
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, S_ql::FT, S_qi::FT, Δt::FT, dqvdt::FT) where {FT}
    new_q = q_new_exponential_part_only(q, S_ql, S_qi, Δt, dqvdt)
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql_Δt::FT, S_qi_Δt::FT, dqvdt_Δt::FT) where {FT}
    # This version should be more floating point stable -- whatever you calculate in mixing ratio you should get exactly out.
    new_q = q_new_exponential_part_only(q, q_liq, q_ice, S_ql_Δt, S_qi_Δt, dqvdt_Δt) # use safe version
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql::FT, S_qi::FT, Δt::FT, dqvdt::FT) where {FT}
    # This version should be more floating point stable -- whatever you calculate in mixing ratio you should get exactly out.
    new_q = q_new_exponential_part_only(q, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt) # use safe version
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end
const morrison_milbrandt_2015_get_new_status_helper_EPA = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only


# A_c_func_no_WBF_exponential_part_only(FT) = FT(0) # don't use this just use FT(0) where this comes up...
# A_c_func_exponential_part_only(τ_ice::FT, q_sl::FT, q_si::FT, Γ_l::FT, L_i::FT, dqsl_dT::FT, c_p::FT) where {FT} =  - (q_sl - q_si) / (τ_ice * Γ_l) * (1+(L_i/c_p * dqsl_dT)) # Eq C4
const A_c_func_no_WBF_exponential_part_only = A_c_func_no_WBF
const A_c_func_exponential_part_only = A_c_func

const A_c_func_no_WBF_EPA = A_c_func_no_WBF_exponential_part_only
const A_c_func_EPA = A_c_func_exponential_part_only

function δ_func_exponential_part_only(A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT}
    term_1 = (δ_0 - A_c * τ)
    term_2 = exp(-Δt / τ)
    prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
    return  A_c * τ  + prod # EQN C5
end
const δ_func_EPA = δ_func_exponential_part_only

δi_func_exponential_part_only(δ_0::FT, A_c::FT, τ::FT, Δt::FT, q_sl::FT, q_si::FT) where {FT} = δ_func_exponential_part_only(A_c, τ, δ_0, Δt) + (q_sl - q_si) # EQN C5
const δi_func_EPA = δi_func_exponential_part_only
# indiv for if we need this for just ice with δi and no liquid change
const δi_func_indiv_exponential_part_only = δ_func_exponential_part_only
const δi_func_indiv_EPA = δi_func_indiv_exponential_part_only


# Could be useful for floating point accuracy when dδ < nextfloat(δ) maybe -- seems not, prod_1 and prod_2 differ by neg sign and small amount can lead to prod = nextfloat(prod_1)-prod_1 rather than true value
function dδ_func_exponential_part_only(A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT} 
    # term = 1 - exp(-Δt / τ)
    term = -expm1(-Δt / τ) # this is more precise
    prod_1 = iszero(term) ? FT(0) : ((A_c * τ) * term) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
    prod_2 = iszero(term) ? FT(0) : (-δ_0 * term) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
    prod = prod_1 + prod_2
    return iszero(term) ? FT(0) : prod # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
end
const dδ_func_EPA = dδ_func_exponential_part_only

""" One species τ_species = τ, no WBF """
function S_func_indiv_no_WBF_exponential_part_only( A_c::FT, τ::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    if isfinite(τ) # avoid NaN problems at τ = Inf
        term_1 = (δ_0 - A_c * τ) / (Δt) # can't be nan bc it's finite * finite, neither of which should be 0
        # term_2 = (1 - exp(-Δt / τ))
        term_2 = -expm1(-Δt / τ) # this is more precise for small exp(-Δt / τ) leading to floating point problems... [expm1 is more precise than 1-exp]
        prod_ = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem

        # if !isfinite(term_1) || !isfinite(term_2) || !isfinite(prod) # if any of these are NaN, we'll just return 0
        # end

        return (A_c + prod_) / Γ # QCCON EQN C6 # Use this form bc we can keep the first term if prod was gonna be NaN
        # return A_c + (δ_0 - A_c * τ)  / (Δt) * (1 - exp(-Δt / τ)) # QCCON EQN C6 [ did I deprecate this to create the lower structure? need to document better ] This way is ad bc if prod is NaN you also lose the first term.
    else
        return FT(0) # I know it looks like it should be A_c bc τ and τ_c canceled, but infinite τ means no change... Really what this is saying is τ_c == τ == τ_<other> = ∞ , the limit at τ -> ∞ of τ/τ_c goes to 0 not 1. It's a simplification of S_func_no_WBF_exponential_part_only().
    end
end
# const S_func_indiv_no_WBF_EPA = S_func_indiv_no_WBF_exponential_part_only
#
const S_ql_func_indiv_exponential_part_only = S_func_indiv_no_WBF_exponential_part_only
const S_qi_func_indiv_exponential_part_only = S_func_indiv_no_WBF_exponential_part_only
const S_ql_func_indiv_EPA = S_ql_func_indiv_exponential_part_only
const S_qi_func_indiv_EPA = S_qi_func_indiv_exponential_part_only



""" One of two species i.e. τ_species ≠ τ, no WBF """
function S_func_no_WBF_exponential_part_only(A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    # if τ is Inf, then τ_c and the other τ must have both been Inf
        # if (!isfinite(τ) && isfinite(τ_c)); error("τ is Inf so τ_c and the other τ must have both been Inf but got τ = Inf and τ_c = $τ_c"); end

    # if τ_c is Inf but τ is not, then τ other must not have been Inf. Either way, S is 0.
    if !isfinite(τ_c) # avoid NaN problems at τ = Inf
        return FT(0)
    else
        # # @debug "S_func_no_WBF_exponential_part_only: A_c = $A_c; τ = $τ; τ_c = $τ_c; δ_0 = $δ_0; Δt = $Δt"
        # term_1 = (δ_0 - A_c * τ) *  τ / (Δt * τ_c) # can't be nan bc it's finite * finite, neither of which should be 0
        term_1 = (δ_0 - A_c * τ) *  (τ / τ_c)/(Δt) # can't be nan bc it's finite * finite, neither of which should be 0
        if isinf(term_1) #
            # We'll just return the largest number we can... it's still smaller than the correct answer so hopefully we don't end up in any infinite loops... [the tests should hopefully catch problems?]
            # return floatmax(FT) * sign(term_1) # this is a problem, we need to fix it. 
            alt = SΔt_func_no_WBF_exponential_part_only(A_c, τ, τ_c, δ_0, Δt, Γ)  / Δt # this is a problem, we need to fix it.
            return isinf(alt) ? floatmax(FT) * sign(alt) : (alt)
        end

        # term_2 = (1 - exp(-Δt / τ)) # this term can be imprecise for small exp(-Δt / τ) leading to floating point problems...
        term_2 = -expm1(-Δt / τ) # this is more precise for small exp(-Δt / τ) leading to floating point problems... [expm1 is more precise than 1-exp]
        prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
        S =  A_c * τ / τ_c + prod # QCCON EQN C6 # Use this form bc we can keep the first term if prod was gonna be NaN
        # return (sign(S) == sign(δ_0)) ? S : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
        return S / Γ
    end
end

const S_func_no_WBF_EPA = S_func_no_WBF_exponential_part_only
#
const S_ql_func_no_WBF_exponential_part_only = S_func_no_WBF_exponential_part_only
const S_ql_func_EPA = S_ql_func_no_WBF_exponential_part_only
const S_qi_func_no_WBF_exponential_part_only = S_func_no_WBF_exponential_part_only
const S_qi_func_no_WBF_EPA = S_qi_func_no_WBF_exponential_part_only

function SΔt_func_no_WBF_exponential_part_only(A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    # if τ_c is Inf but τ is not, then τ other must not have been Inf. Either way, S is 0.
    if !isfinite(τ_c) # avoid NaN problems at τ = Inf
        return FT(0)
    else
        # term_1 = (δ_0 - A_c * τ) *  τ / (Δt * τ_c) # can't be nan bc it's finite * finite, neither of which should be 0
        term_1 = (δ_0 - A_c * τ) *  (τ / τ_c) # can't be nan bc it's finite * finite, neither of which should be 0
        if isinf(term_1) #
            return floatmax(FT) * sign(term_1) # this is a problem, we need to fix it. We'll just return the largest number we can... it's still smaller than the correct answer so hopefully we don't end up in any infinite loops... [the tests should hopefully catch problems?]
        end

        # term_2 = (1 - exp(-Δt / τ))
        term_2 = -expm1(-Δt / τ) # this is more precise for small exp(-Δt / τ) leading to floating point problems... [expm1 is more precise than 1-exp]
        prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
        SΔt =  (A_c * τ * Δt) / τ_c + prod # QCCON EQN C6 # Use this form bc we can keep the first term if prod was gonna be NaN

        # I think this is bad bc WBF is in A so you can get a different sign than δ_0 from WBF
        # return (sign(SΔt) == sign(δ_0)) ? SΔt : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
        return SΔt / Γ
    end
end

const SΔt_func_no_WBF_EPA = SΔt_func_no_WBF_exponential_part_only
const SΔt_ql_func_no_WBF_exponential_part_only = SΔt_func_no_WBF_exponential_part_only
const SΔt_ql_func_EPA = SΔt_ql_func_no_WBF_exponential_part_only
const SΔt_qi_func_no_WBF_exponential_part_only = SΔt_func_no_WBF_exponential_part_only
const SΔt_qi_func_no_WBF_EPA = SΔt_qi_func_no_WBF_exponential_part_only


""" Both species - default is WBF on ice, but can call the other way with appropriate replacements"""
function S_func_exponential_part_only(A_c::FT, τ::FT, τ_ice::FT, δ_0::FT, Δt::FT, Γ::FT, q_sl::FT, q_si::FT) where {FT}
    S = S_func_no_WBF_exponential_part_only( A_c, τ, τ_ice, δ_0, Δt, Γ) + (q_sl - q_si) / (τ_ice * Γ) #  # QICON Eqn C7 # if τ_ice is Inf this works fine,  no issues...
    # @error("I think this line is bad bc if you're passing in δ_0 and operating on WBF sign can flip")
    # return (sign(S) == sign(δ_0)) ? S : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
    return S
end

const S_func_EPA = S_func_exponential_part_only
const S_qi_func_exponential_part_only = S_func_exponential_part_only
const S_qi_func_EPA = S_qi_func_exponential_part_only

""" Both species - default is WBF on ice, but can call the other way with appropriate replacements"""
function SΔt_func_exponential_part_only(A_c::FT, τ::FT, τ_ice::FT, δ_0::FT, Δt::FT, Γ::FT, q_sl::FT, q_si::FT) where {FT}
    SΔt = SΔt_func_no_WBF_exponential_part_only( A_c, τ, τ_ice, δ_0, Δt, Γ) + (q_sl - q_si) * Δt / (τ_ice * Γ) #  # QICON Eqn C7 # if τ_ice is Inf this works fine,  no issues...
    # @error("I think this line is bad bc if you're passing in δ_0 and operating on WBF sign can flip")
    # return (sign(SΔt) == sign(δ_0)) ? SΔt : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
    return SΔt
end
const SΔt_func_EPA = SΔt_func_exponential_part_only
const SΔt_qi_func_exponential_part_only = SΔt_func_exponential_part_only
const SΔt_qi_func_EPA = SΔt_qi_func_exponential_part_only



get_t_out_of_q_no_WBF_EPA(δ_0::FT, A_c::FT, τ::FT, τ_c::FT, q_c::FT, Γ::FT, exit_if_fail::Bool = true) where{FT} = get_t_out_of_q_no_WBF(δ_0, A_c, τ, τ_c, q_c, Γ, exit_if_fail) # easier to use this than rewrite, just fill in Γ = 1
const get_t_out_of_q_liq_EPA = get_t_out_of_q_no_WBF_EPA
const get_t_out_of_q_ice_no_WBF_EPA = get_t_out_of_q_no_WBF_EPA

get_t_out_of_q_WBF_EPA(δ_0::FT, A_c::FT, τ::FT, τ_c::FT, q_ice::FT, Γ::FT, q_sl::FT, q_si::FT, exit_if_fail::Bool = true) where {FT} = get_t_out_of_q_WBF(δ_0, A_c, τ, τ_c, q_ice, Γ, q_sl, q_si, exit_if_fail) # easier to use this than rewrite, just fill in Γ = 1
const get_t_out_of_q_ice_EPA = get_t_out_of_q_WBF_EPA

# Deprecated for now, I don't think you can actually reach the exact equilibrium point in q_new_exponential_part_only bc we're ignoring external forcings...
# see https://www.wolframalpha.com/input?i2d=true&i=solve+-Divide%5Bx%2Cc%5D%3DDivide%5B%5C%2840%29x%2Ba%5C%2841%29%2Ck%5D+for+x
# function calc_WBF_δ_equilibrium(q_sl::FT, q_si::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT) where {FT}
#     # see writeup
#     return -(q_sl - q_si) * (Γ_l * τ_liq) / ( (Γ_i * τ_ice + Γ_l * τ_liq))
# end

# # see https://www.wolframalpha.com/input?i2d=true&i=solve+-Divide%5Bx%2Cc%5D%3DDivide%5B%5C%2840%29x%2Ba%5C%2841%29%2Ck%5D+for+x
# function calc_WBF_δ_equilibrium_exponential_part_only(q_sl::FT, q_si::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT) where {FT}
#     # see writeup
#     return -(q_sl - q_si) * (τ_liq) /  (τ_ice + τ_liq)
# end
# calc_WBF_δ_equilibrium_EPA = calc_WBF_δ_equilibrium_exponential_part_only


# S_func_exponential_part_only() can have nextfloat problem for large numbers [e.g. nextfloat(S_func_no_WBF_exponential_part_only( A_c, τ, τ_ice, δ_0, Δt)) can be arbitrarily large as S_func_no_WBF_exponential_part_only grows, sometimes leading to S_ql, S_qi that are too big...]
# we could switch to BigFloats which would solve thigns but they're oh so slow and allocate like crazy. Just use this stopgap for now...

# [[ these don't properly account for w or dTdt so deprecate... ]]
# """ 
#     Unlike S_above_floor(), this does rely on the fact that δ_0, δ_0i are not affected by external forcings and so are only consumed by cond/evap sub/dep.
#     In the face of w, mixing, and other forcings, maybe you would jus want to fall back to big float. [or maybe something like DoubleFloats.jl for speed]

#     `q_other` is the other species... q_liq if q_ice is being updated and vice versa.

#     We also need to decrease gains but not increase losses, so we use switch state
# """
# function S_below_ceiling_EPA(S_qc::FT, q_other::FT, δ_0::FT, Δt::FT, dqvdt::FT=FT(0)) where {FT} # Don't exceed consuming all supersat + q_other in a timestep
#     δ_0 += max(dqvdt, 0) * Δt # add in the change from dqvdt, as that's also available for use
#     (S_qc > FT(0)) ? min(S_qc, (δ_0 + q_other) / Δt, floatmax(FT)) : max(S_qc, -floatmax(FT)) # Assuming gains are limited by first removing any subsaturation then either consuming other species or consuming supersaturation. if not gaining, then don't use limit.
   
#     # temp form for debugging problems. The line above should work fine... once all the kinks (elsewhere not in this func) are worked out
#     if S_qc > FT(0)
#         if (δ_0 + q_other) / Δt < FT(0)
#             error("δ_0 + q_other < 0 while S_qc > 0")
#         end
#         return min(S_qc, (δ_0 + q_other) / Δt, floatmax(FT)) # Assuming gains are limited by first removing any subsaturation then either consuming other species or consuming supersaturation.
#     else
#         return max(S_qc, -floatmax(FT)) # if not gaining, then don't use limit.
#     end
# end
# @inline S_below_ceiling_EPA(S_ql::FT, S_qi::FT, q_liq::FT, q_ice::FT, δ_0::FT, δ_0i::FT, Δt::FT, dqvdt::FT=FT(0)) where {FT} = S_below_ceiling_EPA(S_ql, q_ice, δ_0, Δt, dqvdt), S_below_ceiling_EPA(S_qi, q_liq, δ_0i, Δt, dqvdt) 


"""
From a list of numbers, return the number with the smallest magnitude
"""
@inline function smallest_magnitude(numbers::Vararg{FT}) where {FT}
    abs_values = abs.(collect(numbers)) # Collect to Vector for argmin
    idx = argmin(abs_values)
    return numbers[idx]
end

# # only returns times in the future. x can be of any sign. so can dxdt
# @inline function get_t_hit_zero(x::FT, dxdt::FT) where {FT}
#     if iszero(dxdt)
#         return FT(Inf) # if dxdt is zero, then it will never hit zero, so we return Inf
#     end

#     # now we need to calculate the time it takes to hit zero. if -x/dxdt is negative, then it will never hit zero.
#     t = -x / dxdt # if dxdt is zero, then t will be Inf, which is what we want
#     if t < FT(0)
#         return FT(Inf)
#     else
#         return t
#     end
# end



function clamp_δ(δ_0::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, q_sl::FT, q_si::FT) where {FT}
    if below_freezing
        if regime_type <: Subsaturated
            return min(δ_0, -(q_sl - q_si))
        elseif regime_type <: WBF
            return clamp(δ_0, -(q_sl - q_si), FT(0))
        elseif regime_type <: Supersaturated
            return max(δ_0, FT(0))
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if regime_type <: Subsaturated
            return min(δ_0, FT(0))
        elseif regime_type <: WBF
            return clamp(δ_0, FT(0), q_si - q_sl)
        elseif regime_type <: Supersaturated
            return max(δ_0, q_si - q_sl)
        else
            error("Unknown regime type: $regime_type")
        end
    end
end         
clamp_δ(δ_0::FT, regime_type::Type{<:AbstractSaturationRegime}, q_sl::FT, q_si::FT) where {FT} = clamp_δ(δ_0, regime_type, regime_type.parameters[3], q_sl, q_si) # don't use this one
clamp_δ(δ_0::FT, regime::AbstractSaturationRegime, q_sl::FT, q_si::FT) where {FT} = clamp_δ(δ_0, typeof(regime), is_below_freezing(regime), q_sl, q_si)


function clamp_δi(δ_0i::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, q_sl::FT, q_si::FT) where {FT}
    if below_freezing
        if regime_type <: Subsaturated
            return min(δ_0i, FT(0))
        elseif regime_type <: WBF
            return clamp(δ_0i, FT(0), q_sl - q_si)
        elseif regime_type <: Supersaturated
            return max(δ_0i, q_sl - q_si)
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if regime_type <: Subsaturated
            return min(δ_0i, -(q_si - q_sl))
        elseif regime_type <: WBF
            return clamp(δ_0i, -(q_si - q_sl), FT(0))
        elseif regime_type <: Supersaturated
            return max(δ_0i, FT(0))
        else
            error("Unknown regime type: $regime_type")
        end
    end
end
clamp_δi(δ_0i::FT, regime_type::Type{<:AbstractSaturationRegime}, q_sl::FT, q_si::FT) where {FT} = clamp_δi(δ_0i, regime_type, regime_type.parameters[3], q_sl, q_si) # don't use this one
clamp_δi(δ_0i::FT, regime::AbstractSaturationRegime, q_sl::FT, q_si::FT) where {FT} = clamp_δi(δ_0i, typeof(regime), is_below_freezing(regime), q_sl, q_si)
    

#= [[ I think the clamps are fine bc on the decline side, you still can't lose more than you have, and on the growth side, we always have Γ ≥ 1, so if anything the bound is too loose... ]]
    However, we're not getting the portion that comes from dTdt making more vapor available...


    This doesn't account for `w`, we really should be using A_c_no_WBF = dδdt_no_S

=#
function clamp_S_ql(S_ql::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, δ_0::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT}

    #  if dTdt < FT(0)
    #     return S_qi # it's unclear how much vapor this will make avilable but it certainly will bring some
    # end

    # dδdt_no_S = dqvdt + (-dqsl_dT * dTdt) # this is the change in δ_0 due to external forcings, not due to phase changes. This is the change in supersaturation due to external forcings. Note the change in dqsl_dT is a negative since δ = q_v - q_sl

    if below_freezing
        if (regime_type <: Subsaturated) || (regime_type <: WBF)
            return safe_clamp(S_ql, max(-q_liq / Δt, -floatmax(FT)), FT(0))
        elseif regime_type <: Supersaturated
            # @debug "clamping S_ql = $S_ql to between 0 and $(min((δ_0 + max(dδdt_no_S*Δt, FT(0)) + q_ice) / Δt, floatmax(FT))); δ_0=$δ_0; q_ice=$q_ice; Δt=$Δt; dδdt_no_S=$dδdt_no_S; out = $(clamp(S_ql, FT(0), min((δ_0 + max(dδdt_no_S*Δt, FT(0)) + q_ice) / Δt, floatmax(FT))))" # this is just growth limting so i think this is fine.
            # return safe_clamp(S_ql, FT(0), min((δ_0 + q_ice) / Δt, floatmax(FT)))
            # return clamp(S_ql, extrema((FT(0), min((δ_0 + q_ice) / Δt, floatmax(FT)), min((δ_0_new + q_ice) / Δt, floatmax(FT))))) # use extrema to avoid problems with sign flipping [should work]
            return clamp(S_ql, FT(0), min((δ_0 + max(dδdt_no_S*Δt, FT(0)) + q_ice) / Δt, floatmax(FT))) # this is just growth limting so i think this is fine. clamp() instead of safe_clamp() should work, if you get an error you're probably in the wrong regime.
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if (regime_type <: Subsaturated)
            return clamp(S_ql, max(-q_liq / Δt, -floatmax(FT)), FT(0))
        elseif (regime_type <: WBF) || (regime_type <: Supersaturated) # growth
            # return clamp(S_ql, FT(0), min((δ_0 + q_ice) / Δt, floatmax(FT)))
            # return clamp(S_ql, extrema((FT(0), min((δ_0 + q_ice) / Δt, floatmax(FT)), min((δ_0_new + q_ice) / Δt, floatmax(FT))))) # use extrema to avoid problems with sign flipping [ should work ]
            return  clamp(S_ql, FT(0), min((δ_0 + max(dδdt_no_S*Δt, FT(0)) + q_ice) / Δt, floatmax(FT))) # this is just growth limting so i think this is fine. clamp() instead of safe_clamp() should work, if you get an error you're probably in the wrong regime.
        else
            error("Unknown regime type: $regime_type")
        end
    end
end
clamp_S_ql(S_ql::FT, regime_type::Type{<:AbstractSaturationRegime}, δ_0::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S_ql(S_ql, regime_type, regime_type.parameters[3], δ_0, q_liq, q_ice, Δt, dδdt_no_S) # don't use this one
clamp_S_ql(S_ql::FT, regime::AbstractSaturationRegime, δ_0::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S_ql(S_ql, typeof(regime), is_below_freezing(regime), δ_0, q_liq, q_ice, Δt, dδdt_no_S)

function clamp_S_qi(S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT}

    # if dTdt < FT(0)
    #     return S_qi # it's unclear how much vapor this will make avilable but it certainly will bring some
    # end

    # dδdt_no_S = dqvdt + (-dqsi_dT * dTdt) # this is the change in δ_0i due to external forcings, not due to phase changes. This is the change in supersaturation due to external forcings. Note the change in dqsi_dT is a negative since δ = q_v - q_si

    if below_freezing
        if (regime_type <: Subsaturated)
            return safe_clamp(S_qi, max(-q_ice / Δt, -floatmax(FT)), FT(0))
        elseif (regime_type <: WBF) || (regime_type <: Supersaturated)
            # @debug "clamping S_qi = $S_qi to between 0 and $(min((δ_0i + max(dδdt_no_S*Δt, FT(0)) + q_liq) / Δt, floatmax(FT))); δ_0i=$δ_0i; q_liq=$q_liq; Δt=$Δt; dδdt_no_S=$dδdt_no_S; out = $(clamp(S_qi, FT(0), min((δ_0i + max(dδdt_no_S*Δt, FT(0)) + q_liq) / Δt, floatmax(FT))))" # this is just growth limting so i think this is fine.
            # return safe_clamp(S_qi, FT(0), min((δ_0i + q_liq) / Δt, floatmax(FT)))
            # return clamp(S_qi, extrema((FT(0), min((δ_0i + q_liq) / Δt, floatmax(FT)), min((δ_0i_new + q_liq) / Δt, floatmax(FT))))) # use extrema to avoid problems with sign flipping
            return clamp(S_qi, FT(0), min((δ_0i + max(dδdt_no_S*Δt, FT(0)) + q_liq) / Δt, floatmax(FT))) #  this is just growth limting so i think this is fine, clamp() instead of safe_clamp() should work, if you get an error you're probably in the wrong regime.
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if (regime_type <: Subsaturated) || (regime_type <: WBF)
            return safe_clamp(S_qi, max(-q_ice / Δt, -floatmax(FT)), FT(0))
        elseif (regime_type <: Supersaturated)
            return FT(0) # clamp(S_qi, FT(0), FT(0)) # no ice growth below freezing, no decay while supersaturated. w/o crossing over to WBF
        else
            error("Unknown regime type: $regime_type")
        end
    end
end
clamp_S_qi(S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S_qi(S_qi, regime_type, regime_type.parameters[3], δ_0i, q_liq, q_ice, Δt, dδdt_no_S)
clamp_S_qi(S_qi::FT, regime::AbstractSaturationRegime, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S_qi(S_qi, typeof(regime), is_below_freezing(regime), δ_0i, q_liq, q_ice, Δt, dδdt_no_S)

clamp_S(S_ql::FT, S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S_ql(S_ql, regime_type, below_freezing, δ_0, q_liq, q_ice, Δt, dδdt_no_S), clamp_S_qi(S_qi, regime_type, below_freezing, δ_0i, q_liq, q_ice, Δt, dδdt_no_S)
clamp_S(S_ql::FT, S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S(S_ql, S_qi, regime_type, regime_type.parameters[3], δ_0, δ_0i, q_liq, q_ice, Δt, dδdt_no_S) # don't use this one
clamp_S(S_ql::FT, S_qi::FT, regime::AbstractSaturationRegime, δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT, dδdt_no_S::FT) where {FT} = clamp_S(S_ql, S_qi, typeof(regime), is_below_freezing(regime), δ_0, δ_0i, q_liq, q_ice, Δt, dδdt_no_S)
            




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

function morrison_milbrandt_2015_style_exponential_part_only(
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
) where {FT}
    (; return_mixing_ratio, fallback_to_standard_supersaturation_limiter, emit_warnings, time_tolerance) = opts

    # TODO: Track supersaturation directly... should fix floating point problems...

    if area > FT(0)

        if emit_warnings && (Δt < eps(FT))
            # @debug "Timestep $(Δt) is very small (smaller than eps(FT) = $(eps(FT))), may cause numerical issues..."
        end

        # if iszero(q_vap) || iszero(TD.vapor_specific_humidity(q))
        #     error("q_vap = 0; q = $q; q_eq = $q_eq; T = $T; p = $p; ρ = $ρ; area = $area; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; ts = $ts; dqvdt = $dqvdt; dTdt = $dTdt")
        # end

        δ_0 = q_vap - q_eq.liq # supersaturation over liquid
        δ_0i = δ_0 + (q_eq.liq - q_eq.ice) # supersaturation over ice



        # below_freezing::Bool = T < TCP.T_freeze(param_set)
        below_freezing::Bool = T < TCP.T_triple(param_set) # In the code this is actually where the saturation vapor pressures are equal... I think it's wrong though
        # regime = get_saturation_regime(q_vap, q, q_eq, below_freezing)

        if iszero(δ_0) || iszero(δ_0i) # possible but unlikely
            (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)
            dδdt_no_S = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)  # Eq C4 no WBF
            dδdt_0 = get_dδdt_0(δ_0, δ_0i, q.liq, q.ice, τ_liq, τ_ice, dδdt_no_S, below_freezing)
            regime = get_saturation_regime(δ_0, δ_0i, q.liq, q.ice, below_freezing; dδdt = dδdt_0) # use this version to break ties and make sure we start going the right direction.
        else
            regime = get_saturation_regime(δ_0, δ_0i, q.liq, q.ice, below_freezing)
        end
        return morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true))::Tuple{FT,FT} # right now we aren't using mixing ratio, so setting return_mixing_ratio to true means that no conversions happen at all. if you ever bring it back, you'd want to combine false on the output here, with true on all nested calls for summations, but we'd need to ensure everything passed to the next layer is not in mixing ratio which we haven't done
    else
        return FT(0), FT(0)
    end
end

# ====================================================================================================================================================================================================== #
# ====================================================================================================================================================================================================== #
# ====================================================================================================================================================================================================== #





# -=-=-=-=-=-=-= Main Dispatcher =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

#=
    Instead of parametric types which complicate stability, we just have concrete types with fields and use a dispatcher function to call the right one.
=#
# function morrison_milbrandt_2015_style_exponential_part_only_dispatcher(
#     regime::Union{Supersaturated{true, true, true}, Supersaturated{true, false, true}, Supersaturated{false, true, true}, Supersaturated{false, false, true}}, # these should all work the same, right? you'll end up with some liq/ice at the end no matter what
#     param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
#     use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
#     return_mixing_ratio::Bool = false,
#     depth::Int = 0,
#     dqvdt::FT = FT(0),
#     dTdt::FT = FT(0),
#     fallback_to_standard_supersaturation_limiter::Bool = false,
#     time_tolerance::FT = FT(1e-8),
#     ) where {FT}

#     error("not implemented yet")
# end

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



"""
    Supersaturated (below freezing)
    - Both can grow!
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Supersaturated{true, true, true}, Supersaturated{true, false, true}, Supersaturated{false, true, true}, Supersaturated{false, false, true}}, # these should all work the same, right? you'll end up with some liq/ice at the end no matter what
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts
    
    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end

    

    # @debug "Calling Supersaturated{$(q.liq > FT(0)), $(q.ice > FT(0)), true}..."


    # --- Thermo  constants ------------------------------------------------------------------------------------ #
   (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)


    # Both are growing..., no explicit need for WBF right
    # A_c = A_c_func_EPA(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    # A_c_no_WBF = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)  # for clamping and fallback
    (; A_c, A_c_no_WBF) = A_c_func_with_and_without_WBF(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    
    τ = τ_func_EPA(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)


    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c_no_WBF, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #

    t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ) # Eq C5

    # @debug "t_hit_liq_sat = $t_hit_liq_sat; A_c = $A_c; A_c_no_WBF = $A_c_no_WBF; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q_sl = $q_sl; q_si = $q_si; q_liq = $q_liq; q_ice = $q_ice; Δt = $Δt; dqvdt = $dqvdt; dTdt = $dTdt; Γ_l = $Γ_l; Γ_i = $Γ_i"

    min_t, i_min_t = find_min_t(t_hit_liq_sat)  # find_min_t helps resolve if min_t is 0 for example, don't skip this call
    if min_t < Δt
        # @debug "will hit liq sat before timestep is over... will transition to wbf at t = $(t_hit_liq_sat)..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # This includes the wbf part though...
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si) #

        if isinf((S_ql+S_qi)*min_t) # scale down just to hit WBF. This can happen when the timescale is too short to calculate. This will crash the model when you go to calculate new_q for example
            liq_frac = τ/τ_liq
            S_ql = δ_0 * liq_frac
            S_qi = δ_0 * (1-liq_frac)
        end

        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF)

        Δt_left = Δt - min_t

        new_δ_0 = FT(0) # hit liq sat from above, so δ_0 = 0
        new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si


        new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt)
        new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, true) # not sure if we can have underflow problems here... don't think so
        
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
        S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
        # @debug "q_ice = $q_ice; S_qi = $S_qi; Δt = $Δt; (q_ice +  S_qi * Δt) = $(q_ice + S_qi * Δt);"
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si)
        # @debug "Before clamp: S_ql = $S_ql; S_qi = $S_qi; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q_sl = $q_sl; q_si = $q_si; q_liq = $q_liq; q_ice = $q_ice; Δt = $Δt; dqvdt = $dqvdt; dTdt = $dTdt; Γ_l = $Γ_l; Γ_i = $Γ_i; dqsl_dT = $dqsl_dT; dqsi_dT = $dqsi_dT;"
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF)
        # @debug "After clamp: S_ql = $S_ql; S_qi = $S_qi;"

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Supersaturated (above freezing)
    - No ice growth, liquid only.
    - Existing ice remains, do not enforce melting here.
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Supersaturated{true, true, false}, Supersaturated{true, false, false}, Supersaturated{false, true, false}, Supersaturated{false, false, false}}, # these should all work the same, right? you'll end up with some liq/ice at the end no matter what
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end
    
    # @debug "Calling Supersaturated{$(q.liq > FT(0)), $(q.ice > FT(0)), false}..."


   (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)

    A_c = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_c_no_WBF = A_c # for clamping and fallback
    τ = τ_liq # supersaturated, so no ice decay, but above freezing, so no ice growth

    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #


    t_hit_ice_sat = t_δ_hit_value(q_si - q_sl, δ_0, A_c, τ_liq) # Eq C5 (can only happen bc of A_c) [ we're above freezing so we hit ice sat first?]
    min_t, i_min_t = find_min_t(t_hit_ice_sat)  # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    if min_t < Δt # bc of A_c

        # @debug "will hit ice sat before timestep is over... will transition to wbf at t = $(min_t)..."
        S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, min_t, Γ_l) # This includes the wbf part though...
        S_qi = FT(0)

        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_ql to not exceed the amount of liquid available

        new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt)

        new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
        new_δ_0i = FT(0) # hit ice sat
        
        Δt_left = Δt - min_t
        # maybe we should check that new_q actually has values, bc if τ or δ or Δt is really small it might not...
        new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, false)
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))

        S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    else
        # exponential decay, so will never hit liq sat
        S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, Δt, Γ_l)
        S_qi = FT(0)
        # @debug "Before clamp: S_ql = $S_ql; S_qi = $S_qi; A_c = $A_c; τ_liq = $τ_liq; δ_0 = $δ_0; q_sl = $q_sl; q_si = $q_si; q_liq = $q_liq; q_ice = $q_ice; Δt = $Δt; dqvdt = $dqvdt; dTdt = $dTdt; Γ_l = $Γ_l; dqsl_dT = $dqsl_dT;"
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt, A_c_no_WBF) 
        # @debug "After clamp: S_ql = $S_ql; S_qi = $S_qi;"
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF has liquid (below freezing)
    - liq can shrink, ice can grow

    - w/ no dqvdt, should asymptote to an equilibrium but not reach it,    
    - you shouldn't be able to hit sat unless dqvdt pulls you there.
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{true, true, true}, WBF{true, false, true}}, # has liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end


    # @debug "Calling WBF{true, $(q.ice > FT(0)), true}..."


    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)

    # One is growing, one is shrinking...
    # A_c = A_c_func_EPA(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    # A_c_no_WBF = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)  # for clamping and fallback
    (; A_c, A_c_no_WBF) = A_c_func_with_and_without_WBF(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4

    τ = τ_func_EPA(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)


    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c_no_WBF, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #


    (t_out_of_liq, t_out_of_liq_valid) = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # Eq C6
    t_hit_ice_sat = t_δ_hit_value(q_si - q_sl, δ_0, A_c, τ) # Eq C5 (can only happen bc of A_c) [ we're above freezing so we hit ice sat first?]
    t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ)

    if !t_out_of_liq_valid # upgrade to BigFloat Call
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            q_vap = TD.vapor_specific_humidity(q)
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) 
        end

        A_c_big = A_c_func_EPA(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τ_big = τ_func_EPA(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
        t_out_of_liq = FT(get_t_out_of_q_liq_EPA(big(δ_0), A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't exit again if fail...

        # @debug "After upgrading to BigFloat, t_out_of_liq = $(t_out_of_liq)"
    end



    # δ_equil = calc_WBF_δ_equilibrium_EPA(q_sl, q_si, τ_liq, τ_ice)
    # t_hit_δ_equil = t_δ_hit_value(δ_equil, δ_0, A_c, τ) # I think this should always be Inf


    # can run out of liq or reach equilibrium (but you'll never reach equilibrium right?)
    # can't reach ice sat bc liq evap would pull you up, can't hit liq sat bc ice would pull you down. equil is guaranteed to be between 0 and (q_sl - q_si) [see writeup]


    min_t, i_min_t = find_min_t(t_out_of_liq, t_hit_liq_sat, t_hit_ice_sat )  # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    if min_t < Δt
        if i_min_t == 1
            # @debug "liq will run out first before timestep is over... will transition at t = $(min_t) to just ice growth"
            S_ql = -q_liq / min_t
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)  # however much happens in that time, rescaled to the timestep [ if this underflows, then what? ]
            Δt_left = Δt - min_t


            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_qi to the timestep, so that it doesn't exceed the amount of ice that can grow in that time


            # dδ = smallest_magnitude(((dqvdt - S_qi) * min_t + q_liq), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc. [deprecate bc it doesn't handle Γ_l properly]
            dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)

            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)
            
            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, -q_liq, S_qi*min_t, dqvdt*min_t) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, true)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1)) 


            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt) # rescale to the timestep
            # @debug "q_ice = $q_ice; S_qi = $S_qi; Δt = $Δt; (q_ice +  S_qi * Δt) = $(q_ice + S_qi * Δt);"
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        elseif i_min_t == 2 # if you hit liq sat
            # @debug "liq will hit sat first before timestep is over... will transition at t = $(min_t) to Supersaturated"
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep [ if this underflows, then what? ]
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)  # however much happens in that time, rescaled to the timestep [ if this underflows, then what? ]
            Δt_left = Δt - min_t

            # @debug "S_ql = $S_ql; S_qi = $S_qi; min_t = $min_t; dqvdt = $dqvdt; dTdt = $dTdt; dqsl_dT = $dqsl_dT; dqsi_dT = $dqsi_dT; δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq; q_ice = $q_ice;"
            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_ql to the timestep, so that it doesn't exceed the amount of liq that can shrink in that time
            # @debug "S_ql after clamp = $S_ql; min_t = $min_t; dqvdt = $dqvdt; dTdt = $dTdt; dqsl_dT = $dqsl_dT;"
            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_qi to the timestep, so that it doesn't exceed the amount of ice that can grow in that time

            # dδ = smallest_magnitude((dqvdt - S_ql - S_qi) * min_t, dδ_func_EPA(A_c, τ, δ_0, min_t)) # [deprecate bc it doesn't handle Γ properly]
            # dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
            

            new_δ_0 = FT(0)
            new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si
            

            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(Supersaturated, new_q.liq, new_q.ice, true)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            # @debug "S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; Δt_left = $Δt_left; Δt = $Δt; min_t = $min_t; dqvdt = $dqvdt; dTdt = $dTdt; dqsl_dT = $dqsl_dT; dqsi_dT = $dqsi_dT;"
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt) # rescale to the timestep
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        else # i_min_t == 3 # if you hit ice sat
            # @debug "ice will hit sat first before timestep is over... will transition at t = $(min_t) to Subsaturated"
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep [ if this underflows, then what? ]
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)  # however much happens in that time, rescaled to the timestep [ if this underflows, then what? ]
            Δt_left = Δt - min_t

            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_ql to the timestep, so that it doesn't exceed the amount of liq that can shrink in that time
            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_qi to the timestep, so that it doesn't exceed the amount of ice that can grow in that time
            # dδ = smallest_magnitude(((dqvdt - S_qi) * min_t + q_liq), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc. [deprecate bc it doesn't handle Γ properly]
            dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
            new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
            new_δ_0i = FT(0) # hit ice sat
            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, true)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt) # rescale to the timestep
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        end
        
    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si)
        # @debug "S_ql = $S_ql; S_qi = $S_qi; Δt = $Δt; (q_liq +  S_ql * Δt) = $(q_liq + S_ql * Δt); (q_ice +  S_qi * Δt) = $(q_ice + S_qi * Δt); δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq; q_ice = $q_ice; dqvdt = $dqvdt"
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF) # clamp S_ql and S_qi to the timestep, so that it doesn't exceed the amount of liq/ice that can shrink/grow in that time
        # @debug "q_ice = $q_ice; S_qi = $S_qi; Δt = $Δt; (q_ice +  S_qi * Δt) = $(q_ice + S_qi * Δt);"
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF no liquid (below freezing)
    - no liq, just ice growth

    -- you should never hit saturation right? it's just exponential decay. though now w/ dqvdt we can...
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{false, true, true}, WBF{false, false, true}}, # no liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end

    # @debug "Calling WBF{false, $(q.ice > FT(0)), true}..."

   (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)


    # no liq so this is just ice growth.
    A_c = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_c_no_WBF = A_c

    τ = τ_ice # we have no liq
     # We only have ice growth so it's like QCCON but with just ice 

    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #

    # can reach ice sat
    t_hit_ice_sat = t_δ_hit_value(FT(0), δ_0i, A_c, τ) # switch to using δ_0i and let that hit 0
    t_hit_liq_sat = t_δ_hit_value(q_sl-q_si, δ_0i, A_c, τ) # (only possible if A_c does it...)

    min_t, i_min_t = find_min_t(t_hit_ice_sat, t_hit_liq_sat)
    # @debug "min_t = $min_t, i_min_t = $i_min_t | t_hit_ice_sat = $t_hit_ice_sat, t_hit_liq_sat = $t_hit_liq_sat"

    # # @debug "A_c = $A_c; τ = $τ; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; Δt = $Δt; q_sl = $q_sl; q_si = $q_si; q_vap = $q_vap; q_liq = $q_liq; q_ice = $q_ice; T=$T"

    if min_t < Δt

        if iszero(dqvdt + dTdt + w)
            error("should be unreachable bc we're doing this w/ no external forcing")
        else
            if i_min_t == 1
                # @debug "will hit ice sat before timestep is over... will transition to subsaturated"
                S_ql = FT(0)
                S_qi = S_qi_func_indiv_EPA(A_c, τ_ice, δ_0i, min_t, Γ_i)
                Δt_left = Δt - min_t

                new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
                new_δ_0i = FT(0) # we're at ice sat so δ_0i = 0
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt)

                new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, true)
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1)) # should have no liq still, no ice above freezing
                # @debug "S_ql = $S_ql; S_qi = $S_qi; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; Δt_left = $Δt_left; min_t = $min_t; q_ice = $q_ice; new_q = $new_q"
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                # @debug "q_ice = $q_ice; S_qi = $S_qi; Δt = $Δt; (q_ice +  S_qi * Δt) = $(q_ice + S_qi * Δt);"
                return  return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
               
            else # i_min_t == 2 # A_c has raised us up to liq_sat, transition to supersaturated? or stay on WBF line?
                # @debug "Reached liq sat, will terminate unless T > T_freeze"
                S_ql = FT(0)
                S_qi = S_qi_func_indiv_EPA(A_c, τ_ice, δ_0i, min_t, Γ_i) #   Δδ_0i = (q_sl - q_si) - δ_0i --> Δq_ice = -Δδ_0i
                Δt_left = Δt - min_t

                new_δ_0 = FT(0) # we're at liq sat so δ_0 = 0
                new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si

                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt)
                new_regime = add_regime_parameters(Supersaturated, new_q.liq, new_q.ice, true) # if T < T_freeze then we have subsat ice, otherwise we have supersat liq
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1)) # should have no liq still
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))

            end
        end
    else

        # @debug "nothing of note through end of timestep..."
        S_ql = FT(0)
        S_qi = S_qi_func_indiv_EPA(A_c, τ, δ_0i, Δt, Γ_i) # use the no WBF version w/ just exponential decay of δ_0i
        S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF) 

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

    
"""
    WBF above freezing yes ice (above freezing) 
    
    can lose ice, hit liq sat or ice sat
    liq can only grow
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{false, true, false},  WBF{true, true, false}},
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(), # if time is shorter than the tolerance, we do a StandardSupersaturation step first, since we can't guarantee success w/ the lambert W methods...
    )::Tuple{FT,FT} where {FT} 
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end
    # @debug "Calling WBF{false, $(q.ice > FT(0)), false}..."

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)
   
    # A_c = A_c_func_EPA(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    # A_c_no_WBF = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)  # for clamping and fallback
    (; A_c, A_c_no_WBF) = A_c_func_with_and_without_WBF(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4

    τ = τ_func_EPA(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)
    

    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c_no_WBF, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #



    (t_out_of_ice, t_out_of_ice_valid) = has_ice(regime) ? get_t_out_of_q_ice_EPA(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si) : (FT(Inf), true) # Eq C6
    if !t_out_of_ice_valid # upgrade to BigFloat Call
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            q_vap = TD.vapor_specific_humidity(q)
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts)
        end
        A_c_big = A_c_func_EPA(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τ_big = τ_func_EPA(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
        t_out_of_ice = FT(get_t_out_of_q_ice_EPA(big(δ_0i), A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_ice = $(t_out_of_ice)"
    end
    t_hit_liq_sat = t_δ_hit_value(FT(0.), δ_0, A_c, τ_liq) # Eq C5
    t_hit_ice_sat = t_δ_hit_value(q_si - q_sl, δ_0, A_c, τ_liq) # Eq C5 (can only happen bc of A_c) [ we're above freezing so we hit ice sat first?]

    min_t, i_min_t = find_min_t(t_out_of_ice, t_hit_liq_sat, t_hit_ice_sat) # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    # @debug "min_t = $min_t, i_min_t = $i_min_t | t_out_of_ice = $t_out_of_ice, t_hit_liq_sat = $t_hit_liq_sat, t_hit_ice_sat = $t_hit_ice_sat"

    if min_t < Δt
        if i_min_t == 1 # out of ice, continue with WBF (we're above freezing so we can still have been making liquid)
            # @debug "Out of ice, continuing with WBF..."
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = -q_ice / min_t

            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t, A_c_no_WBF) # clamp S_ql to the timestep, so that it doesn't exceed the amount of liq that can shrink in that time
            
            Δt_left = Δt - min_t
            # dδ = smallest_magnitude(((dqvdt-S_ql) * min_t + q_ice), dδ_func_EPA(A_c, τ, δ_0, min_t))
            dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ # Eq C8
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)

            # @debug "S_ql = $S_ql; S_qi = $S_qi; dδ = $dδ; new_δ_0 = $new_δ_0; new_δ_0i = $new_δ_0i; A_c = $A_c; τ = $τ; δ_0 = $δ_0; min_t = $min_t; Γ_l = $Γ_l; dqvdt = $dqvdt; dTdt = $dTdt"
            # @debug "new_δ_0 = $new_δ_0; new_δ_0i = $new_δ_0i; dδ = $dδ;"

            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql*min_t, -q_ice, dqvdt*min_t) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(WBF, new_q.liq, false, false)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))

        elseif i_min_t == 2 # hit liq sat, cross into subsaturated
            # @debug "Hit liq sat, crossing into subsaturated regime..."
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)
            # S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt, A_c_no_WBF) 
            # S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF) 
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF) # clamp S_ql and S_qi to the timestep, so that it doesn't exceed the amount of liq/ice that can shrink/grow in that time

            if iszero(dqvdt + dTdt + w)
                error("we shouldn't be able to hit liquid sat from WBF{false, $(has_ice(regime)), false} and dqvdt=dTdt=w=0, it should just asymptote...")
            else
                Δt_left = Δt - min_t
                new_δ_0 = FT(0)
                new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si
                
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
                new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, false) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            end
        else # i_min_t == 3, hit ice sat, go to Supersaturated but above freezing.
            # @debug "Hit ice saturation, transitioning to Supersaturated but above freezing regime"
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si) # However much happens in this amount of time
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF)

            # if iszero(dqvdt + dTdt + w)
            #     error("we shouldn't be able to get up to ice sat from WBF{false, $(has_ice(regime)), false} and dqvdt=dTdt=w=0")  # this is alse, evap can hit ice sat easily...
            # else
            Δt_left = Δt - min_t
            new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
            new_δ_0i = FT(0)
            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(Supersaturated, new_q.liq, new_q.ice, false) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            # end
        end

    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) # Eq C2
        # @debug "Before clamping: S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq; q_ice = $q_ice; Δt = $Δt; dqvdt = $dqvdt; dTdt = $dTdt"
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF)
        # @debug "After clamping: S_ql = $S_ql; S_qi = $S_qi;"

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

    

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
    
"""
    WBF above freezing, no ice.
    
    can have liq growth or hit liq or ice sat.
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{false, false, false}, WBF{true, false, false}},
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(), # if time is shorter than the tolerance, we do a StandardSupersaturation step first, since we can't guarantee success w/ the lambert W methods...
    )::Tuple{FT,FT} where {FT} 
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end
    # @debug "Calling WBF{false, $(q.ice > FT(0)), false}..."

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)

    A_c = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    A_c_no_WBF = A_c # no WBF so A_c = A_c_no_WBF
    τ = τ_liq # we have no ice so τ = τ_liq


    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #

    # @debug "A_c = $A_c; τ = $τ; τ_ice=$τ_ice; τ_liq = $τ_liq; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; Δt = $Δt; q_sl = $q_sl; q_si = $q_si; Γ_l = $Γ_l; Γ_i = $Γ_i; L_l = $L_l; L_i = $L_i; c_p = $c_p; e_sl = $e_sl; e_si = $e_si; dqsl_dT=$dqsl_dT; dqsi_dT=$dqsi_dT; w=$w; g=$g; p=$p; ρ=$ρ; q_vap = $(TD.vapor_specific_humidity(q)); q_liq = $q_liq; q_ice = $q_ice; T=$T; T_freeze=$T_freeze; q_sl = $q_sl; q_si = $q_si; q_liq = $q_liq; q_ice = $q_ice; dqvdt = $dqvdt; dTdt = $dTdt"


    t_hit_liq_sat = t_δ_hit_value(FT(0.), δ_0, A_c, τ_liq) # Eq C5
    t_hit_ice_sat = t_δ_hit_value(q_si - q_sl, δ_0, A_c, τ_liq) # Eq C5 (can only happen bc of A_c) [ we're above freezing so we hit ice sat first?]

    min_t, i_min_t = find_min_t(t_hit_liq_sat, t_hit_ice_sat) # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    # @debug "min_t = $min_t, i_min_t = $i_min_t | t_hit_liq_sat = $t_hit_liq_sat, t_hit_ice_sat = $t_hit_ice_sat"

    if min_t < Δt
        if i_min_t == 1 # hit liq sat, cross into subsaturated
            # @debug "Hit liq sat, crossing into subsaturated regime..."
            S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, min_t, Γ_l)
            S_qi = FT(0)
            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt, A_c_no_WBF) 

            if iszero(dqvdt + dTdt + w)
                error("we shouldn't be able to run out of liquid sat from WBF{false, $(has_ice(regime)), false} and dqvdt=dTdt=w0, it should just asymptote...")
            else
                Δt_left = Δt - min_t
                new_δ_0 = FT(0)
                new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
                new_regime = add_regime_parameters(Subsaturated, new_q.liq, false, false) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            end
        else # i_min_t == 2, hit ice sat, go to Supersaturated but above freezing.
            # @debug "Hit ice saturation, transitioning to Supersaturated but above freezing regime"
            S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, min_t, Γ_l)
            S_qi = FT(0)
            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt, A_c_no_WBF) 

            if iszero(dqvdt + dTdt + w)
                error("we shouldn't be able to get up to ice sat from WBF{false, $(has_ice(regime)), false} and dqvdt=dTdt=w=0")
            else
                Δt_left = Δt - min_t
                new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
                new_δ_0i = FT(0)
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
                new_regime = add_regime_parameters(Supersaturated, new_q.liq, false, false) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            end
        end

    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, Δt, Γ_l)
        S_qi = FT(0)
        # @debug "Before clamping: S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq; q_ice = $q_ice; Δt = $Δt; dqvdt = $dqvdt; dTdt = $dTdt"
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt, A_c_no_WBF) 
        # @debug "After clamping: S_ql = $S_ql;"

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated with liq and ice
    - both shrink
    - if below freezing, would hit ice_sat first and stop
    - if above freezing, would hit liq_sat first and stop
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{true, true, true}, Subsaturated{true, true, false}}, # can run out of either first
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts


    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end


   (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)

    BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{true, true, $BF}..."

    # A_c = A_c_func_EPA(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    # A_c_no_WBF = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)  # for clamping and fallback
    (; A_c, A_c_no_WBF) = A_c_func_with_and_without_WBF(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4

    τ = τ_func_EPA(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)



    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"

        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c_no_WBF, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #

    t_hit_sat = BF ? t_δ_hit_value(q_si-q_sl, δ_0, A_c, τ) : t_δ_hit_value(FT(0), δ_0, A_c, τ) # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower

    (t_out_of_liq, t_out_of_liq_valid) = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq, Γ_l)
    (t_out_of_ice, t_out_of_ice_valid) = get_t_out_of_q_ice_EPA(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si)

    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            q_vap = TD.vapor_specific_humidity(q)
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts)
        end
        # upgrade to BigFloat Call
        A_c_big = A_c_func_EPA(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τ_big = τ_func_EPA(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
        t_out_of_liq = FT(get_t_out_of_q_liq_EPA(big(δ_0), A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_liq = $(t_out_of_liq)"
    end
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            q_vap = TD.vapor_specific_humidity(q)
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts)
        end
        # upgrade to BigFloat Call
        A_c_big = (@isdefined A_c_big) ? A_c_big :  A_c_func_EPA(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        A_c_big += big(dqvdt)
        τ_big = (@isdefined τ_big) ? τ_big : τ_func_EPA(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
        t_out_of_ice = FT(get_t_out_of_q_ice_EPA(big(δ_0), A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_ice = $(t_out_of_ice)"
    end

    t_hit_sat = BF ? t_δ_hit_value(q_si-q_sl, δ_0, A_c, τ) : t_δ_hit_value(FT(0), δ_0, A_c, τ) # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower
    
    min_t, i_min_t = find_min_t(t_out_of_liq, t_out_of_ice, t_hit_sat)

    # # @debug "min_t = $min_t; i_min_t = $i_min_t; Δt = $Δt | t_out_of_liq = $t_out_of_liq; t_out_of_ice = $t_out_of_ice; t_hit_sat = $t_hit_sat"

    if min_t < Δt
        if i_min_t == 1
            # @debug "liq will run out first before timestep is over... will transition at t = $(min_t) to just ice decay if ice is present"
            S_ql = -q_liq / min_t
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)
            Δt_left = Δt - min_t

            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF)

            # dδ = smallest_magnitude(((dqvdt-S_qi) * min_t + q_liq), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc. [deprecate bc it doesn't handle Γ properly]
            dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)

            
            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, -q_liq, S_qi * min_t, dqvdt*min_t) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        elseif i_min_t == 2
            # @debug "ice will run out first before timestep is over... will transition at t = $(min_t) to just liq decay if liq is present"
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi =  -q_ice / min_t
            Δt_left = Δt - min_t

            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF)
 
            # dδ = smallest_magnitude(((dqvdt-S_ql) * min_t + q_ice), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc. [deprecate bc it doesn't handle Γ properly]
            dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)


            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql * min_t, -q_ice, dqvdt*min_t) # use multiplied form for floating point accuracy
            new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        else # i_min_t == 3
            # @debug "will hit sat before timestep is over... will transition to wbf at t = $(min_t) if liq is present otherwise we're stuck..."
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)
            Δt_left = Δt - min_t

            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF)

            if BF
                new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
                new_δ_0i = FT(0) # hit ice sat from below so δ_0i = 0
            else
                new_δ_0 = FT(0) # hit liq sat from below so δ_0 = 0
                new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si
            end

            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt)
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            # @debug "S_ql = $S_ql; S_qi = $S_qi; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit"
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
            # @debug "after resolve_S_S_addit: S_ql = $S_ql; S_qi = $S_qi"
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        end
    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si)
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end


end


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated liq but no ice
    - no ice so liq shrink
    - if below freezing, would hit ice_sat first and stop
    - if above freezing, would hit liq_sat first and stop
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{true, false, true}, Subsaturated{true, false, false}}, # can run out of liq first. if not and we make it to ice sat, transition to WBF
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end

   (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)

    BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{true, false, $BF}..."


    A_c = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_c_no_WBF = A_c # for clamping and fallback
    τ = τ_liq # we have no ice so ice can't grow or shrink while subsaturated

    
    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #

    (t_out_of_liq, t_out_of_liq_valid) = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # Eq C6
    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            q_vap = TD.vapor_specific_humidity(q)
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts)
        end
        # upgrade to BigFloat Call
        A_c_big = A_c_func_no_WBF_EPA(big(q_sl), big(g), big(w), big(c_p), big(e_sl), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τ_big = big(τ_liq)
        t_out_of_liq = FT(get_t_out_of_q_liq_EPA(big(δ_0), A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_liq = $(t_out_of_liq)"
    end

    t_hit_sat = BF ? t_δ_hit_value(q_si-q_sl, δ_0, A_c, τ) : t_δ_hit_value(FT(0), δ_0, A_c, τ) # below freezing, stop at ice_sat which is lower, above freezing, stop at liq_sat which is lower

    min_t, i_min_t = find_min_t(t_out_of_liq, t_hit_sat)

    if min_t < Δt
        if i_min_t == 1
            # @debug "liq will run out first at t = $(min_t) before timestep is over..."
            S_ql = -q_liq / min_t # can't continue once have nothing
            S_qi = FT(0) # no ice so no ice growth

            if iszero(dqvdt + dTdt + w) # we're done bc we're not subsaturated with nothing happening
                S_ql = -q_liq / Δt # rescale to the timestep
                return return_mixing_ratio ? (S_ql, FT(0)) : (S_mixing_ratio_to_shum(S_ql, q.tot), FT(0))
            else
                Δt_left = Δt - min_t
                # dδ = smallest_magnitude(((dqvdt) * min_t + q_liq), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc. 
                dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
                new_δ_0 = δ_0 + dδ
                new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
                new_δ_0i = δ_0i + dδ
                new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, -q_liq, S_qi*min_t, dqvdt*min_t) # use multiplied form for floating point accuracy
                new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, BF)
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
                # @debug "S_ql = $S_ql; S_qi = $S_qi; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit"
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                # @debug "after resolve_S_S_addit: S_ql = $S_ql; S_qi = $S_qi"
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            end


        else # i_min_t == 2
            # @debug "will hit sat before timestep is over... will transition at t = $(min_t) to wbf"
            S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, min_t, Γ_l)
            S_qi = FT(0)
            Δt_left = Δt - min_t

            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t, A_c_no_WBF) 

            if BF
                new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl 
                new_δ_0i = FT(0) # hit ice sat from below so δ_0i = 0
            else
                new_δ_0 = FT(0) # hit liq sat from below so δ_0 = 0
                new_δ_0i = q_sl - q_si # diff bewteen where we're at (q_sl) and q_si
            end

            new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) 
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
            # @debug "S_ql = $S_ql; S_qi = $S_qi; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit"
            S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
            # @debug "after resolve_S_S_addit: S_ql = $S_ql; S_qi = $S_qi"
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        end
    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_indiv_EPA(A_c, τ_liq, δ_0, Δt, Γ_l)
        S_qi = FT(0)
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt, A_c_no_WBF) 

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated ice but no liq 
    - ice can shrink
    - below freezing, would hit ice_sat first, above freezing, would hit liq_sat first
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{false, true, true}, Subsaturated{false, true, false}},
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts

    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end


    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)

    BF::Bool = T < T_freeze
    # ================================================================================= #

    # @debug "Calling Subsaturated{true, false, $BF}..."

    A_c = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_c_no_WBF = A_c # for clamping and fallback
    τ = τ_ice # we have no liq

    # ===== Fallback Block ===== #
    standard_milestone_t, standard_milestone, S_ql, S_qi, δ_eq, δi_eq = calculate_next_standard_milestone_time(regime, q_eq, q_liq, q_ice, δ_0, δ_0i, T<T_freeze, τ_liq, τ_ice; dδdt_no_S = A_c_no_WBF, Γ_l=Γ_l, Γ_i=Γ_i, allow_δ_eq_point = true) # we need to allow the eq point because otherwise we risk WBF oscillations, see note in do_standard_fallback()
    # @debug "standard_milestone_t = $standard_milestone_t; standard_milestone = $standard_milestone; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dqvdt = $dqvdt; dTdt = $dTdt"
    if (standard_milestone_t < time_tolerance) && !(standard_milestone == NotAtSupersaturationMilestone) && !(standard_milestone == AtSupersaturationStationaryPointMilestone) # 0 means never hitting a milestone again, 3 means eq point which we don't recognize in this framwork.
        # @debug "falling bacc"
        return do_standard_fallback(
            standard_milestone_t, standard_milestone, time_tolerance, S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, A_c, Γ_l, Γ_i,
            regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; opts = update(opts; return_mixing_ratio = true)
            )
    end
    # =========================== #


    
    (t_out_of_ice, t_out_of_ice_valid) = get_t_out_of_q_ice_no_WBF_EPA(δ_0i, A_c, τ, τ_ice, q_ice, Γ_i) # just like QCCON but w/ only ice
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            q_vap = TD.vapor_specific_humidity(q)
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts)
        end
        # upgrade to BigFloat Call
        A_c_big = A_c_func_no_WBF_EPA(big(q_sl), big(g), big(w), big(c_p), big(e_sl), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τ_big = big(τ_ice)
        t_out_of_ice = FT(get_t_out_of_q_ice_no_WBF_EPA(big(δ_0i), A_c_big, τ_big, τ_big, big(q_ice), big(Γ_i), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_ice = $(t_out_of_ice)"
    end

    t_hit_sat = BF ? t_δ_hit_value(FT(0), δ_0i, A_c, τ) : t_δ_hit_value(q_sl-q_si, δ_0i, A_c, τ) # below freezing, stop at ice_sat which is lower, above freezing, stop at liq_sat which is lower [ in δ_0i terms bc that's what _indiv methods are accurate for w/ no WBF ?]

    # ================================================================================= #

    min_t, i_min_t = find_min_t(t_out_of_ice, t_hit_sat)

    if min_t < Δt # either way it's the same up to 
        if i_min_t == 1
            # @debug "ice will run out first before timestep is over at t = $(min_t)..."
            S_ql = FT(0)
            S_qi = -q_ice / min_t # we're done so just scale to the entire timestep

            if iszero(dqvdt + dTdt + w)
                S_qi = -q_ice / Δt # if dqvdt is 0, then we just scale to the entire timestep
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            else
                Δt_left = Δt - min_t
                dδ = dδ_func_EPA(A_c, τ, δ_0, min_t)
                new_δ_0 = δ_0 + dδ
                new_δ_0i = δ_0i + dδ
                new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
                new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql*min_t, -q_ice, dqvdt*min_t)
                new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, BF)
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
                # @debug "S_ql = $S_ql; S_qi = $S_qi; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit"
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                # @debug "after resolve_S_S_addit: S_ql = $S_ql; S_qi = $S_qi"
                # @debug "q_ice = $q_ice; S_qi = $S_qi; Δt = $Δt; (q_ice +  S_qi * Δt) = $(q_ice + S_qi * Δt);"
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            end

        else # i_min_t == 2
            # @debug "will hit sat before timestep is over... Then we're done bc we don't have liq in BF and if above freezing, we can't form ice."
            S_ql = FT(0)
            S_qi = S_qi_func_indiv_EPA(A_c, τ_ice, δ_0i, min_t, Γ_i)
            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t, A_c_no_WBF) # clamp to the regime


            if iszero(dqvdt + dTdt + w)
                S_qi *= min_t / Δt
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            else
                Δt_left = Δt - min_t
                new_δ_0 = q_si - q_sl # diff bewteen where we're at (q_si) and q_sl
                new_δ_0i = FT(0) # we're at ice sat so δ_0i = 0
                new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt)
                new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, BF)
                S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
                # @debug "S_ql = $S_ql; S_qi = $S_qi; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit"
                S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
                # @debug "after resolve_S_S_addit: S_ql = $S_ql; S_qi = $S_qi"
                return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
            end
            
        end

    else
        # @debug "nothing of note through end of timestep..."
        S_ql = FT(0)
        S_qi = S_qi_func_indiv_EPA(A_c, τ_ice, δ_0i, Δt, Γ_i)
        S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, Δt, A_c_no_WBF) 
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
"""
    Subsaturated with nothing
    - return 0, 0
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{false, false, true}, Subsaturated{false, false, false}}, # kind of null, you have nothing and can't make anything
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, δ_0::FT, δ_0i::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015EPAOpts{FT} = MM2015EPAOpts{FT}(),
    )::Tuple{FT,FT} where {FT}
    (; return_mixing_ratio, depth, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, time_tolerance) = opts



    if depth ≥ 10
        @error "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; δ_0 = $δ_0; δ_0i = $δ_0i; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dTdt = $dTdt; dqvdt = $dqvdt;"
        error("Failed to converge after 10 iterations")
    end


    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; opts = opts)
    BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{false, false, $(BF)}..."


    
    A_c = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) # do the full thing bc we dont know when dqvdt, dTdt, w will make you hit sat... but no WBF bc no condensate...
    # A_c_no_WBF = A_c # for clamping and fallback

    # no tau bc no condensate...., essentially τ is Inf, so dδ/dt = A_c...

    # I believe we don't need a fallback block since we're just using t_out_of_x... 

    if BF
        t_hit_sat = t_out_of_x(δ_0, A_c)  # if we use this we don't have to use anything else
    else
        t_hit_sat = t_out_of_x(δ_0i, A_c)
    end

    min_t = t_hit_sat # no other time to hit anything, so just use this

    if min_t < Δt
        S_ql = FT(0)
        S_qi = FT(0)
        new_δ_0 = BF ? (q_sl-q_si) : FT(0) # if below freezing, we hit ice sat first so δ_0 = δ_0 + t_hit_sat * dqvdt, if above freezing, we hit liq sat first so δ_0 = 0, which should just be (q_sl - q_si)
        new_δ_0i = BF ? FT(0) : (q_si - q_sl) # if below freezing, we hit ice sat first so δ_0i = 0, if above freezing, we hit liq sat first so δ_0i = δ_0i + t_hit_sat * dqvdt, which should just be (q_si - q_sl)

        Δt_left = Δt - min_t

        new_regime = add_regime_parameters(WBF, FT(0), FT(0), BF)
        new_q = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt) # use multiplied form for floating point accuracy
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_δ_0, new_δ_0i, new_q, q_eq, Δt_left, ts; opts = update(opts; return_mixing_ratio = true, depth = depth+1))
        S_ql, S_qi = resolve_S_S_addit(FT(0), FT(0), min_t, S_ql_addit, S_qi_addit, Δt_left, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    else
        return FT(0), FT(0) 
    end
end




