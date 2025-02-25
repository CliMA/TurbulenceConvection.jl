"""
These rates are inherently timestep dependent, because they average the source over the timestep...
If you want the WBF continuous calculations without the timestep averageing, you can only pray.

-- this means that if your dt is highly adaptive, the relative importance of the tendencies from this source relative to other sources will change significantly, potentially causing discontinuties in the solution
    Perhaps only use when you can be reasonably confident that dt is not changing too much, e.g. within a factor of 2 or so or smoothly changing from a spinup to a steady state dt
    This is not as bad as for some other tendencies though, because the source should scale with the timestep for reasonably small timesteps.
    If τ_liq/τ_ice are very fast relative to your timestep (i.e. cond/evap sub/dep are generally expected to saturate before the timestep is over), then this becomes more of a concern.

    It would be nice if there were a way to get that point out for the timestepper, but it's not so trivial.


    - compute_tendency_dt_max() in callbacks now does respect supersaturation depletion, so it should be more okay now to not use this method for the timestep limiter

    We don't have a good way of making this timestep agnostic but still getting the limiter-free WBF logic
    - this is because cond/evap sub/dep are naturally limited by the exponential decay
    - WBF can be left unbounded to compose w/ other tendencies 
        > i.e WBF from liq to ice supported by advective convergence or sgs flux rather than local formation... However that does mix up euleriean and lagrangian views
        > This also still doesn't solve anything for cond/evap sub/dep, which perhaps would also want to compose in the same way...
        > When using adaptive timestepping based on tendency limiters, we calculate the tendency and then find dt assuming those tendencies are constant over the timestep.
        > For this method, this is not true, so if Δt is changed, in principle the tendency should be recalculated...
            >> So really, this would need an adaptive timestepper with steps proposed then tendency calculated then accepted or rejected based on stability, but we don't have that rn.
        
        - compute_tendency_dt_max() in callbacks now does respect supersaturation depletion, so it should be more okay now rely on that timestep limiter


Note also I've noted some numerical issues below ~eps(FT)... maybe a warning should be issued for that
"""
function morrison_milbrandt_2015_style(
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::Real,
    ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    emit_warnings::Bool = true,
) where {FT}
    """
    See https://doi.org/10.1175/JAS-D-14-0065.1
    we are ignoring the mixing and radiation terms for our short timesteps, as well as rain effects

    this *shouldn't* need limiters the way it's defined but be careful I suppose lol...
    """
    thermo_params = TCP.thermodynamics_params(param_set) # currently repeated in both places, pare down later
    # microphys_params = TCP.microphysics_params(param_set)
    S_ql = FT(0)
    S_qi = FT(0)
    if area > 0

        if emit_warnings && (Δt < eps(FT))
            @warn "Timestep $(Δt) is very small (smaller than eps(FT) = $(eps(FT))), may cause numerical issues..."
        end

        """
        # Milbrandt's equations are in mixing ratio not specific humidity, to convert we need to use

        # q_spe = q_mix / (1 + q_mix)   [ in reality this is more like q_spe_c = q_mix_c / (1 + q_tot_mix) , bc 1+q_tot = total ratio all water / dry air...]
        # q_mix = q_spe / (1 - q_spe)   [ in reality this is more like q_mix_c = q_spe_c / (1 - q_tot_spe) , bc 1-q_tot = dry mixing ratio...]

        # [NOTE: q_tot appears instead of just q because we have more than one water species so we need to be careful...]

        # note then d/dy (q_spe) = d/dy(q_mix / (1 + q_mix)) = d/dy(q_mix) / (1 + q_mix)^2  [ really it's more like d/dy(q_spe_c) = d/dy(q_mix_c / (1 + q_tot_mix)) and q_tot_mix should be a constant under phase transformation... ]
        # note then d/dy (q_mix) = d/dy(q_spe / (1 - q_spe)) = d/dy(q_spe) / (1 - q_spe)^2  [ really it's more like d/dy(q_mix_c) = d/dy(q_spe_c / (1 - q_tot_spe)) and q_tot_spe should be a constant under phase transformation... ]

        NOTE: These mixing ratio, specific humidity, and their derivatives are essentially identical at earthlike conditions...
        """

        local A_c::FT
        local τ::FT
        local δ_0::FT


        # --- Thermo ------------------------------------------------------------------------------------ #
        g = TCP.grav(param_set) # acceleration of gravity
        L_i = TD.latent_heat_sublim(thermo_params, ts) # Latent heat for ice (L_s)
        L_l = TD.latent_heat_vapor(thermo_params, ts)  # Latent heat for water
        # c_p = TD.TP.cp_d(thermo_params) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?
        c_p = TD.cp_m(thermo_params, q) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?

        e_sl = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) # saturation vapor pressure, or TD.partial_pressure_vapor(thermo_params, p, q_eq) 

        # analytical forms for change in saturation vapor presure with temperature
        dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, ts) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix... is it a good enough approx?
        dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, ts) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...

        q_vap_sat = TD.vapor_specific_humidity(thermo_params, ts)
        dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(1), T, q_vap_sat) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix... is it a good enough approx?
        dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(0), T, q_vap_sat) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...

        dqsl_dT /= (1 - q_eq.liq)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.liq = q_sl at T
        dqsi_dT /= (1 - q_eq.ice)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.ice = q_si at T

        # So, we can do everything in mixing ratio and then convert back at the end instead of trying to convert everythig in their derivation to figure out what q_tot_mix
        q_vap = TD.shum_to_mixing_ratio(q_vap, q.tot)
        q_liq = TD.shum_to_mixing_ratio(q.liq, q.tot)
        q_ice = TD.shum_to_mixing_ratio(q.ice, q.tot)
        q_sl = TD.shum_to_mixing_ratio(q_eq.liq, q.tot) # q_eq is the equilibrium mixing ratio, q_sl is the saturation mixing ratio, q_eq we made to contain only saturation vapor pressure values...
        q_si = TD.shum_to_mixing_ratio(q_eq.ice, q.tot)

        T_freeze = TCP.T_freeze(param_set)

        # --- Constants ------------------------------------------------------------------------------------ #

        δ_0 = q_vap - q_sl # supersaturation over liquid
        δ_0i = δ_0 + (q_sl - q_si) # supersaturation over ice

        # Balanced at large timestep, but needs limiters (otherwise, may consume more liq/ice than exists... but handles WBF directions fine and most stuff for short timesteps
        if use_fix # setting these makes the asymptote ok, not sure if it breaks other things....
            # i think their equation balances liq evap <-> ice growth, but doesn't limit it properly over the timestep, technically q_sl-q_si addition/subtraction should decay as you approach saturation...
            L_i = L_l
            dqsi_dT = dqsl_dT
            # Does evaporating water releasing more latent heat than forming ice takeup risk evaporating a cloud completely?
        end

        Γ_l = 1 + L_l / c_p * dqsl_dT  # Eqn C3
        Γ_i = 1 + L_i / c_p * dqsi_dT  # Eqn C3
        @info "" L_i, L_l, c_p, dqsl_dT, dqsi_dT
        @info "Γ_l, Γ_i" Γ_l, Γ_i
        # Γ_l = 1 + L_l / c_p * dqsi_dT 
        # Γ_i = 1 + L_l / c_p * dqsi_dT 


        dTdt_mix = FT(0) # ignore for now
        dTdt_rad = FT(0) # ignore for now




        if δ_0i < 0 # always evap for both... (only lasts until run out of one...)
            if iszero(q_liq)
                τ = τ_ice
            elseif iszero(q_ice)
                τ = τ_liq
            else
                τ = 1 / (1 / τ_liq + (1 + (L_i / c_p) * dqsl_dT) * ((1 / τ_ice) / (Γ_i))) # Eqn C2
            end

            A_c =
                dTdt_mix - (q_sl * ρ * g * w) / (p - e_sl) - dqsl_dT * (dTdt_rad + dTdt_mix - (w * g) / c_p) -
                (q_sl - q_si) / (τ_ice * Γ_l) * (1 + (L_i / c_p) * dqsl_dT) # Eq C4

            # t_out_of_liq [δ - δ_0 = -q_liq/Δt]
            print(( A_c * τ - δ_0 ) / ( A_c * τ - -q_liq  - δ_0 ))
            t_out_of_liq = τ * log( ( A_c * τ - δ_0 ) / ( A_c * τ - q_liq  - δ_0 )) # solution

            # t_out_of_ice [δ - δ_0 = -q_ice/Δt]
            t_out_of_ice = τ * log( ( A_c * τ - δ_0 ) / ( A_c * τ - q_ice  - δ_0 )) # solution

            # t_hit_ice_sat [δ = -(q_sl - q_si)]
            @info "" (δ_0 - A_c * τ ) (  -(q_sl - q_si)  - A_c * τ )
            t_hit_ice_sat = τ * log( ( δ_0 - A_c * τ ) / (  -(q_sl - q_si)  - A_c * τ )) # solution

            @info "t_out_of_liq, t_out_of_ice, t_hit_ice_sat" t_out_of_liq, t_out_of_ice, t_hit_ice_sat

            min_t, i_min_t = findmin([t_out_of_liq, t_out_of_ice, t_hit_ice_sat])
            if min_t < Δt
                if i_min_t == 1
                    @warn "liq will run out before timestep is over... will transition to just ice decay if ice is present"
                    if !iszero(q_ice) # if liq ran out first and ice is 0, then we just do nothing
                        S_ql_here = -q_liq / Δt
                        S_qi_here = error("not implemented")
                        q_new = TD.PhasePartition(q.tot, q.liq - S_ql_here * Δt, q.ice - S_qi_here * Δt)
                        Δt_new = Δt - t_out_of_liq
                        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q_new, q_eq, Δt - t_out_of_liq, ts)
                    else
                        S_ql_here = -q_liq / Δt
                        S_qi
                    end

                elseif i_min_t == 2
                    @warn "ice will run out before timestep is over... will transition to just liq decay if liq is present"
                else
                    @warn "will hit ice sat before timestep is over... will transition to wbf if liq is present"
                end
            end

            # make a recurisve call starting from that point?
            # end
            @warn("doesn't cover what happens when one becomes zero first...")
        elseif (δ_0 < 0 < δ_0i) # WBF (evap and sublimation)
            τ = 1 / (1 / τ_liq + (1 + (L_i / c_p) * dqsl_dT) * ((1 / τ_ice) / (Γ_i))) # Eqn C2
            # t_out_of_liq

            # t_hit_ice_sat
            @warn("Not implemented yet")
        else # growth for both (only lasts up to liq sat then WBF begins)
            τ = 1 / (1 / τ_liq + (1 + (L_i / c_p) * dqsl_dT) * ((1 / τ_ice) / (Γ_i))) # Eqn C2
            # t_hit_liq_sat
        end

        @info "τ" τ τ_liq τ_ice (τ / τ_liq) (τ / τ_ice)

        # --- Calculations ------------------------------------------------------------------------------------ #

        if T < T_freeze # WBF flips at T = T_freeze (not that it's well defined) We are not robust do dT/dt causing a flip aross freezing.

            A_c =
                dTdt_mix - (q_sl * ρ * g * w) / (p - e_sl) - dqsl_dT * (dTdt_rad + dTdt_mix - (w * g) / c_p) -
                (q_sl - q_si) / (τ_ice * Γ_l) * (1 + (L_i / c_p) * dqsl_dT) # Eq C4


            # S_ql = A_c * τ / (τ_liq * Γ_l) + (δ_0 - A_c * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ))   # QCCON EQN C6 [ did I deprecate this to create the lower structure? need to document better ]

            # -- components of QCCON Eq C6 -- #
            β = -(q_sl - q_si) / (τ_ice * Γ_l) * (1 + (L_i / c_p) * dqsl_dT) # WBF part of A_c
            WBF_cont = β * τ / (τ_liq * Γ_l) + (-β * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ)) # WBF part contribution to QCCON Eq C6 from all terms [ goes from liq → ice] # While the ice side is ice growth ltd, why is the loss of liq at the same rate? Are we just assuming ice is faster?
            α = (A_c - β) # non-WBF part of A_c
            α_cont = α * τ / (τ_liq * Γ_l) + (-α * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ)) # non-WBF part contribution to QCCON Eq C6 from all terms
            δ_0_cont = δ_0 * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ)) # supersaturation contirbution to QCCON Eq C6 from all terms

            @info "liq" WBF_cont α_cont δ_0_cont


            #= WBF is parameterized as a constant flux from liquid to ice for all time, thus some care is required to limits its effect =#
                # we can move neg cond/evap to the WBF bucket since WBF is just evap
                # cond increases the WBF contribution limit
                # WBF plus cond/evap cannot exceed what liq can provide
                # Thus, cond_evap + WBF ≥ -q_liq/ Δt --> cond_evap ≥ -q_liq / Δt - WBF

                # WBF ≤ 0 always, and if WBF + cond_evap < -q_liq / Δt, shrink cond_evap not WBF

            cond_evap_cont = δ_0_cont + α_cont # This is equivalent to QCCON - WBF_cont. Some of this could go to ice by WBF if WBF is not 0... let's just say if it would have been WBF'd remove it

            # WBF_cont_orig = WBF_cont # save so can use in ice
            if cond_evap_cont ≥ 0 # condensation
                WBF_cont = max(WBF_cont, -q_liq / Δt - cond_evap_cont) # WBF cannot exceed what liq can provide from existing plus cond. evap doesn't count against WBF, we just assume it's part of it.
                @info "" WBF_cont cond_evap_cont  Δt
                @info "" δ_0 cond_evap_cont*Δt
            else
                # cond_evap_cont = max(cond_evap_cont, -q_liq / Δt) # no more evap than liq can provides
                WBF_cont = max(WBF_cont, -q_liq / Δt) # no more WBF than liq can provide
                cond_evap_cont = max(cond_evap_cont, -q_liq / Δt - WBF_cont) # WBF has priority, cond_evap_cont gets the leftovers if those are smaller than cond_evap_cont was originally
            end

            # @info "WBF_cont" WBF_cont Δt

            @info "liq after " WBF_cont cond_evap_cont cond_evap_cont*Δt (-q_liq / Δt ) δ_0
            S_ql = cond_evap_cont + WBF_cont # This is the total contribution to liquid tendency from cond_Evap and WBF
            @info "S_ql" S_ql
            # WBF_limit = -max(WBF_cont, -q_liq / Δt) # WBF cannot exceed what what liq can provide by evaporation  [ why is this now in S_i units? also why do we not just call this WBF_cont?]
            # @info "WBF_limit" WBF_limit

            # evap_cont = max(evap_cont - WBF_cont,) # no negative evap
            # S_ql = max(WBF_cont + evap_cont, -q_liq / Δt) # This seems bad, we just re-added the WBF back in? should we not be using the WBF limit?
            # S_ql = max(WBF_limit + cond_evap_cont, -q_liq / Δt) # Ok so cond_evap can give or take, WBF can only take, sum of the two can't take more than liq can provide. Bc WBF is really just evap, we'll assume that WBF is not shrinking.

            # -- components of QICON Eq C7 -- #
            # We find it hard to use the components bc they disagree in Γ and some other such constants between the liq and ice eqn. so we just take the ice one and modify it...
            # Let us instead presume that WBF liq losses do not need to be balanced by WBF ice gains, some can be lost in vapor.

        
            # S_qi =
            #     A_c * τ / (τ_ice * Γ_i) +
            #     (δ_0 - A_c * τ) * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ)) +
            #     (q_sl - q_si) / (τ_ice * Γ_i) #+  # QICON Eqn C7

            
            # Then, the WBF term is

            WBF_cont_ice = β * τ / (τ_ice * Γ_i) + (-β * τ) * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ)) +
                ((q_sl - q_si) / (τ_ice * Γ_i)  / Δt) * (1 - exp(-Δt / τ))
                # (q_sl - q_si) / (τ_ice * Γ_i) / Δt   # supersaturation contribution to QICON Eq C7 from all terms ( I think we add this term here and not in WBF)

            WBF_cont_ice *= 0 # fix for now
            @warn("tesing w/ no wbf_cont_ice")
            @info "" WBF_cont_ice
            # @info "(q_sl - q_si) / (τ_ice * Γ_i)" (q_sl - q_si) / (τ_ice * Γ_i) Δt

        
            α_cont_ice = α * τ / (τ_ice * Γ_i) + (-α * τ) * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ_ice)) # non-WBF part contribution to QICON Eq C7 from all terms
            δ_0_cont_ice = (δ_0 * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ_ice))) + # supersaturation contribution to QICON Eq C7 from all terms ( I think we add this term here and not in WBF)
            ((q_sl - q_si) * τ / (τ_ice * Γ_i *  Δt)) * (1 - exp(-Δt / τ_ice))
    
            # my fix... really though it's not perfect when not in a WBF regime...


            @info "ice before " WBF_cont_ice α_cont_ice δ_0_cont_ice
            sub_dep_cont = δ_0_cont_ice + α_cont_ice # This is equivalent to QICON - WBF_cont. Some of this could go to ice by WBF if WBF is not 0... let's just say if it would have been WBF'd remove it
            # WBF_cont_ice_orig = WBF_cont_ice # save so can use in ice
            if sub_dep_cont ≥ 0 # deposition
                WBF_cont_ice = min(WBF_cont_ice, -WBF_cont) # Don't let it get too positive
                @info "" WBF_cont_ice WBF_cont
                # WBF_cont_ice = max(WBF_cont_ice, -q_ice / Δt - sub_dep_cont) # WBF cannot exceed what ice can provide from existing plus sub/dep. sub/dep doesn't count against WBF, we just assume it's part of it.
            else # sublimation
                # @info "before sub_dep_cont" sub_dep_cont
                # @info "WBF_cont" WBF_cont
                WBF_cont_ice = min(WBF_cont_ice, -WBF_cont) # Don't let it get too positive
                # @info "WBF_cont_ice" WBF_cont_ice
                sub_dep_cont = max(sub_dep_cont, (-q_ice / Δt - WBF_cont_ice)) # WBF has priority, sub_dep_cont gets the leftovers if those are smaller than sub_dep_cont was originally
                # @info "after sub_dep_cont" sub_dep_cont
                # @info " -q_ice / Δt"  -q_ice / Δt

            end

            @info "ice after " WBF_cont_ice sub_dep_cont (-q_ice / Δt) δ_0i sub_dep_cont*Δt
            S_qi = sub_dep_cont + WBF_cont_ice # This is the total contribution to ice tendency from sub/dep and WBF


            # δ_q_cont = (q_sl - q_si) / (τ_ice * Γ_i)
            # β_ice_cont = β * τ / (τ_ice * Γ_i) +  (-β * τ) * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ)) 
            # δ_0_ice_cont = δ_0 *  τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ))

            # WBF_cont_ice = -max(-WBF_cont, -q_liq / Δt) # not more than liq can provide (is this not redundant?)
            # @debug "this is just min(WBF_cont, liq_provide)? and isn't this redundant? why do we even need a separate wbf_cont_ice?"

            # if ice growing
            # if S_qi > 0
            #     # liq_provide = min(max(-S_ql, 0), q_liq / Δt) # [is this redundant?] How much can liq provide
            #     # @debug "is this redundant?"
            #     # WBF_cont_ice = -max(-WBF_cont, -liq_provide) # not more than liq can provide [this is just min(WBF_cont, liq_provide)?] how is this different than the WBF_limit calc?
            #     # WBF_cont_ice = min(WBF_cont, liq_provide) # smaller of the two (both should be positive of 0 I believe)
            #     # WBF_cont = min(S_qi, WBF_limit) # ? This is nevver used again...
            #     # @debug "why is this never used again? should we not be taking WBF out of the equation?"
            #     S_qi = min(S_qi, (δ_0 + q_sl - q_si) / Δt - S_ql) # ice growth cannot exceed ice supersaturation (δ_0 + q_sl - q_si) minus liquid growth taken out (why does liq get priority?)
            #     @debug "why does liq get priority here"
            # # if ice shrinking
            # else # if ice is shrinking,
            #     # Any WBF going to ice is less than what's being evaporated...
            #     # @info "WBF_cont_ice" WBF_cont_ice

            #     # @warn("something here is wrong...")
            #     # WBF_cont = min(S_qi, WBF_limit) # Whwat does this do? The WBF contributoin is the min of the WBF_limit and S_qi? why?
            #     # S_qi += (WBF_cont_ice - WBF_cont) # switch to our limited WBF contribution from liquid

            #     # We'd like to change from using WBF_cont_ice. 
            #     # @info "WBF_limit" WBF_limit
            #     # @info "S_qi" S_qi
            #     # @info "WBF_cont_ice" WBF_cont_ice
            #     # @info "WBF_cont" WBF_cont
            #     # This is pointless bc WBF_limit is positive and S_qi is negative
            #     # WBF_cont = max(S_qi, WBF_limit) # Whwat does this do? The WBF contributoin is the min of the WBF_limit and S_qi? why? changed to max bc it's neg but i think it wants the smaller number...
            #     # WBF_cont = WBF_limit
            #     S_qi += (WBF_cont - WBF_cont_orig) # switch out the unlimited WBF_cont_ice for the limited WBF_limit
            # end
            # S_qi = max(S_qi, -q_ice / Δt) # don't let it take more ice than exists

        # T > T_freeze (no WBF, no ice formation, etc...)    
        else
            # WBF applies backwards above T_freeze (in principle from saturation vapor pressures, this is not a well defined state)
            # while it may be prudent to stil allow evap, refreezing is not something that should be happening.
            # so we just do nothing
            # it's debatable if we should reformulate all the equations... right now we just do nothing.

            A_c = dTdt_mix - (q_sl * ρ * g * w) / (p - e_sl) - dqsl_dT * (dTdt_rad + dTdt_mix - (w * g) / c_p) # Eq C4 w/ WBF q_si - q_sl term removed

            # S_ql = A_c * τ / (τ_liq * Γ_l) + (δ_0 - A_c * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ))   # QCCON EQN C6


            β = FT(0)
            β_cont = β * τ / (τ_liq * Γ_l) + (-β * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ))
            α = (A_c - β)
            α_cont = α * τ / (τ_liq * Γ_l) + (-α * τ) * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ))
            δ_0_cont = δ_0 * τ / (Δt * τ_liq * Γ_l) * (1 - exp(-Δt / τ))

            evap_cont = δ_0_cont + α_cont # none of this should go to ice... there's no WBF above freezing...
            S_ql = max(evap_cont, -q_liq / Δt) #  what happend to eq 6?

            # WBF cannot happen above freezing

            S_qi =
                A_c * τ / (τ_ice * Γ_i) +
                (δ_0 - A_c * τ) * τ / (Δt * τ_ice * Γ_i) * (1 - exp(-Δt / τ))  # QICON Eqn C w/ WBF q_si - q_sl term removed

            S_qi = max(S_qi, -q_ice / Δt) # don't let it take more ice than exists

            S_qi = clamp(S_qi, -q_ice / Δt, 0) # don't let it take more ice than exists, don't let ice grow

        end


        # Go from mixing ratio world back to specific humidity world
        q_tot_mix = q.tot / (1 - q.tot) # there's no such thing as q_tot_mix... right?
        S_ql /= (1 + q_tot_mix) # bc we showed  d/dt(q_spe) = d/dt(q_mix) / (1+q_tot_mix) under phase transformation, which makes sense, it's smaller
        S_qi /= (1 + q_tot_mix)

    end

    # fix NaNs from infinite timescales
    if isinf(τ_liq) 
        S_ql = FT(0)
    end
    if isinf(τ_ice) 
        S_qi = FT(0)
    end

    # S_ql *= Γ_l * (τ_liq / τ) # why are these factors left over?
    # S_qi *= Γ_i * (τ_ice / τ) # why are these factors left over?

    @info "S_ql, S_qi" S_ql, S_qi
    # S_ql *= Γ_l # why are these factors left over?
    # S_qi *= Γ_i # why are these factors left over?


    if emit_warnings && (Δt < eps(FT))
        @info "T = $(T), q_tot = $(q.tot), q_vap = $(q_vap), q_liq = $(q.liq), q_ice = $(q.ice), S_l = $(q_vap - q_eq.liq), S_i = $(q_vap - q_eq.ice) | τ_liq = $(τ_liq), τ_ice = $(τ_ice) | S_ql = $(S_ql), S_qi = $(S_qi), Δt = $(Δt)"
    end

    return S_ql, S_qi #, (δ_Δt - δ_0)/Δt
end



function morrison_milbrandt_2015_style_exponential_part_only(
    param_set::APS,
    area::FT,
    ρ::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    q_eq::TD.PhasePartition,
    Δt::Real,
) where {FT}
    """
    This will be a fcn that handles only the exponential decay part for supersaturation w/ liquid and ice, but ignores vertical velocity and everything else (temperature changes and stuff)
    """
    # w = FT(0) # ignore for now
    # return morrison_milbrandt_2015_style(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap, q_eq, Δt)

    S_ql::FT = FT(0)
    S_qi::FT = FT(0)

    local δ_0::FT
    local τ::FT
    local A_c::FT

    # can we cache repeated terms? compiler probably figures out that anyway...
    if area > 0

        T_freeze = TCP.T_freeze(param_set)
        if T < T_freeze # WBF flips at T = T_freeze (not that it's well defined bc freezing should not be allowed)

            δ_0 = q_vap - q_eq.liq # supersaturation over liquid
            τ = 1 / (1 / τ_liq + 1 / τ_ice)
            A_c = -(q_eq.liq - q_eq.ice) / (τ_ice) # Eq C4
            S_ql = A_c * τ / (τ_liq) + (δ_0 - A_c * τ) * τ / (Δt * τ_liq) * (1 - exp(-Δt / τ))   # QCCON EQN C6

            S_qi =
                A_c * τ / (τ_ice) +
                (δ_0 - A_c * τ) * τ / (Δt * τ_ice) * (1 - exp(-Δt / τ)) +
                (q_eq.liq - q_eq.ice) / (τ_ice)  # QICON Eqn C7
            
        else # we cannot allow ice formation (there's also no WBF)

            # NOTE: These equations probably should be changed somehow above freezing. just clamping shouldn't be enough.

            δ_0 = q_vap - q_eq.liq # supersaturation over liquid
            τ = 1 / (1 / τ_liq + 1 / τ_ice)
            A_c = FT(0) # Eq C4 w/ extraneous and WBF terms removed
            S_ql = A_c * τ / (τ_liq) + (δ_0 - A_c * τ) * τ / (Δt * τ_liq) * (1 - exp(-Δt / τ))   # QCCON EQN C6

            S_qi =
                A_c * τ / (τ_ice) +
                (δ_0 - A_c * τ) * τ / (Δt * τ_ice) * (1 - exp(-Δt / τ))  # QICON Eqn C7 w/ extraneous and wbf terms removed

            S_qi = clamp(S_qi, -q.ice / Δt, 0) # don't let it take more ice than exists, don't let ice grow

        end
    end


    # # liquid consumption limiting
    # if S_ql < -q.liq/Δt
    #     # S_qi -= -S_ql - q.liq/Δt
    #     S_extra = S_ql + q.liq/Δt
    #     S_ql -= S_extra  # Discard the extra part.
    #     # What part of this is from vapor and what part went to ice?


    #     @info "before" S_qi, S_extra, S_qi - S_extra
    #     # @info δ_0i, S_extra, S_qi - δ_0i, S_qi
    #     if S_qi > q.ice/Δt + S_ql
    #         # S_qi += min(S_extra, max(S_qi - δ_0i, 0), 0) # only if ice source is more than initial ice supersat, remove extra liquid contribution, otherwise leave.
    #         I_extra = 0
    #         # @info "I_extra" I_extra
    #         S_qi +=  max(I_extra - S_extra, 0) # only if ice source is more than initial ice supersat, remove extra liquid contribution, otherwise leave.
    #         # S_qi += min(S_extra, S_qi - δ_0i, 0) # neg + pos, remove contribution up to initial ice supersat

    #         # @info "after" S_qi
    #     end
    #     # S_extra = 0
    #     # I_extra = max( δ_0i/Δt - S_qi, 0)
    #     # @info "I_extra" I_extra, S_qi ,  δ_0i
    #     # S_qi +=  min(I_extra + S_extra, 0) 

    #     # S_qi = min(S_qi + I_extra, S_qi + S_extra)
    #     # S_qi += min( I_extra, S_extra)
    #     # I_extra = min(S_qi - δ_0i, 0)
    #     # S_qi += min(I_extra, 0) # only if ice source is more than initial ice supersat, remove extra liquid contribution, otherwise leave.  
    #     # @info "after" S_qi, δ_0i/Δt, S_extra
    # end


    if isinf(τ_liq) 
        S_ql = FT(0)
    end
    if isinf(τ_ice) 
        S_qi = FT(0)
    end

    return S_ql, S_qi

end


# @inline pow_hack(x, y) = exp(y * log(x))

# @inline function saturation_vapor_pressure(
#     param_set::APS,
#     T::FT,
#     LH_0::FT,
#     Δcp::FT,
# ) where {FT <: Real}
#     press_triple = TP.press_triple(param_set)
#     R_v = TP.R_v(param_set)
#     T_triple = TP.T_triple(param_set)
#     T_0 = TP.T_0(param_set)

#     return press_triple *
#            # (T / T_triple)^(Δcp / R_v) *
#            pow_hack(T / T_triple, Δcp / R_v) *
#            exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))

# end


# @inline function saturation_vapor_pressure_derivative_Temperature(
#     param_set::APS,
#     T::FT,
#     LH_0::FT,
#     Δcp::FT,
# ) where {FT <: Real} # d p_s / dT
#     press_triple = TP.press_triple(param_set)
#     R_v = TP.R_v(param_set)
#     T_triple = TP.T_triple(param_set)
#     T_0 = TP.T_0(param_set)

#     return press_triple *
#            # (T / T_triple)^(Δcp / R_v) *
#            pow_hack(T / T_triple, Δcp / R_v) *
#            exp((LH_0 - Δcp * T_0) / R_v * (1 / T_triple - 1 / T))

# end

# @inline function q_vap_saturation_generic(
#     param_set::APS,
#     T::FT,
#     ρ::FT,
#     phase::Phase,
# ) where {FT <: Real}
#     p_v_sat = saturation_vapor_pressure(param_set, T, phase)
#     return q_vap_saturation_from_density(param_set, T, ρ, p_v_sat)
# end
# q_vap_saturation_generic(param_set::APS, T, ρ, phase::Phase) =
#     q_vap_saturation_generic(param_set, promote(T, ρ)..., phase)

# @inline function q_vap_saturation_derivative_Temperature(
#     param_set::APS,
#     ts::TD.ThermodynamicState
# ) where {FT <: Real}


#     return TD.∂q_vap_sat_∂T(param_set, ts)

# end
