"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""

function handle_expr(expr::String; kwargs...)
    # see https://stackoverflow.com/a/57749395, takes a function that accepts kwargs
    # create tuple from args and create func

    if ~contains(expr, "->")
        expr = "(" * join([string(x) for x in keys(kwargs)], ',') * ")" * " -> " * expr
    end

    expr = eval(Meta.parse(expr))
    return Base.invokelatest(expr, values(NamedTuple(kwargs))...) # works but very slow to use invokelatest, maybe try https://github.com/SciML/RuntimeGeneratedFunctions.jl
end


function noneq_moisture_sources(param_set::APS, area::FT, ρ::FT, Δt::Real, ts::TD.ThermodynamicState, w::FT, z::FT; ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    if area > 0
        use_supersat = get_isbits_nt(param_set.user_args, :use_supersat, false) # (:use_supersat in keys(param_set.user_args)) ? param_set.user_args.use_supersat : false # so we dont have to set everything we dont know is in user_args in the defaults...
        use_supersat = isa(use_supersat, Val) ? typeof(use_supersat).parameters[1] : use_supersat # extract the value from the Val type (we put it in there in parameter_set.jl to make it isbits)
        use_korolev_mazin = get_isbits_nt(param_set.user_args, :use_korolev_mazin, false) # (:use_supersat in keys(param_set.user_args)) ? param_set.user_args.use_supersat : false # so we dont have to set everything we dont know is in user_args in the defaults...
        use_korolev_mazin = isa(use_korolev_mazin, Val) ? typeof(use_korolev_mazin).parameters[1] : use_korolev_mazin # extract the value from the Val type (we put it in there in parameter_set.jl to make it isbits)

        τ_use = get_isbits_nt(param_set.user_args, :τ_use, :standard)  # this has w built in though, no way around it, maybe we should write one that's just the exponential decay part w/o anything else... ( # change the default to at least :morrison_milbrandt_2015_style_exponential_part_only soon )
        τ_use = isa(τ_use, Val) ? typeof(τ_use).parameters[1] : τ_use

        q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
        T::FT = TD.air_temperature(thermo_params, ts)
        p::FT = TD.air_pressure(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

        if !(use_supersat == false) # it's either true or a specified value
            # use phase partition in case we wanna use the conv_q_vap fcn but maybe not best for supersat since it's not really a phase partition (all 3 are vapor amounts)

            τ_liq, τ_ice = get_τ(param_set, microphys_params, use_supersat, q, T, p, w, z)

            q_eq = TD.PhasePartition(
                q.tot,
                TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
                TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
            ) # all 3 are vapor amounts
            # @info ("τliq =  $τ_liq,  τice =  $τ_ice, T =  $T, w = $w  qliq =  $(q.liq), qice = $(q.ice) qvap = $q_vap, qveq liq = $(q_eq.liq), qveq ice = $(q_eq.ice) ")

        elseif use_korolev_mazin
            # need to get w into here somewhere...
            S_ql, S_qi = korolev_mazin_2007(param_set, area, ρ, Δt, ts, w)
            q_eq = TD.PhasePartition(
                q.tot,
                TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
                TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
            ) # all 3 are vapor amounts

            # korolev_mazin fix effective tau > 1 (or whatever timestep was realy but we run w/ 1 mostly rn)
            τ_liq_eff = (q_vap - q_eq.liq) / S_ql
            τ_ice_eff = (q_vap - q_eq.ice) / S_qi
            # println("effective τ_liq = ",τ_liq_eff, " effective τ_ice = ",τ_ice_eff, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )
            if (0 < τ_liq_eff) && (τ_liq_eff < 1) # fast source
                # @show("rate limiting cond: ", S_ql, (q_vap - q_eq.liq) / 1)
                S_ql = (q_vap - q_eq.liq) / 1
            elseif (0 < τ_ice_eff) && (τ_ice_eff < 1) # fast source
                # @show("rate limiting dep: ", S_qi, (q_vap - q_eq.ice) / 1)
                S_qi = (q_vap - q_eq.ice) / 1
            elseif (-1 < τ_liq_eff) && (τ_liq_eff < 0) # fast sink
                # @show("rate limiting evap: ", S_ql, (q_vap - q_eq.liq) / 1)
                S_ql = (q_vap - q_eq.liq) / -1
            elseif (-1 < τ_ice_eff) && (τ_ice_eff < 0) # fast sink
                # @show("rate limiting sub: ", S_qi, (q_vap - q_eq.ice) / 1)
                S_qi = (q_vap - q_eq.ice) / -1
            else
            end

        else # basic noneq (no supersat formulation so not likely to be right)
            # TODO - is that the state we want to be relaxing to?
            ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q.tot)
            q_eq = TD.PhasePartition(thermo_params, ts_eq)

            S_ql = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, q)
            S_qi = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, ice_type, q_eq, q)
        end



        # LIMITER (and calculate source)
        if !(use_supersat == false) # might need to do these first bc ql,qc tendencies maybe are applied individually and can still crash the code if one is too large...

            if τ_use == :standard
                S_ql = (q_vap - q_eq.liq) / τ_liq # | microphys_params.τ_cond_evap | CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, TD.PhasePartition(FT(0),q_vap,FT(0)))  
                S_qi = (q_vap - q_eq.ice) / τ_ice # -(source to vapor) = source to condensate

                # Local sources, note these will always have the same sign as the source terms above since they just have a different Δt instead of τ
                Q_vl = (q_vap - q_eq.liq) / Δt # vapor available for liquid growth in this timestep
                Q_vi = (q_vap - q_eq.ice) / Δt # vapor available for ice growth in this timestep (if positive, this is larger?)

                # Limiters, where we can rely on the fact that Q_vl,Q_vi are always the same sign as S_ql,S_qi

                # debug check            
                # if Q_vl < Q_vi
                #     if T > thermo_params.T_triple
                #         @info(
                #             "status",
                #             T,
                #             thermo_params.T_freeze,
                #             thermo_params.T_triple,
                #             Q_vl,
                #             Q_vi,
                #             q_vap,
                #             q_eq.liq,
                #             q_eq.ice,
                #             S_ql,
                #             S_qi
                #         )
                #         error("This should not happen (we should be below triple point?)") # This is just a check that all is working properly, remove later
                #     end
                # elseif Q_vl > Q_vi
                #     if T < thermo_params.T_triple
                #         @info("status", T, Q_vl, Q_vi, q_vap, q_eq.liq, q_eq.ice, S_ql, S_qi)
                #         error("This should not happen (we should be above triple point?)") # This is just a check that all is working properly, remove later
                #     end
                # else
                #     # all good, we're at the triple point or indistinguashable within numerical precision (e.g. at very cold temperatures)
                # end

                # @info( "current", S_ql, S_qi, "--",Q_vl, Q_vi, "---", q.liq / Δt, q.ice / Δt, q_vap/Δt)
                if (S_ql > 0) && (S_qi > 0) # both growing (this is the hardest with the different saturation vapor pressures) - we'll partition based on respective source magnitudes limited by Q_vi which is larger...
                    S_tot = S_ql + S_qi
                    # create liquid and ice in proportion to their source magnitudes and see if we reach liquid saturation

                    if Q_vl < Q_vi
                        if S_tot > Q_vl # this handles S_ql whether or not it's larger than Q_vl while simultaneously handling the first part of S_qi
                            S_ql_Qvl = S_ql * Q_vl / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
                            S_qi_Qvl = S_qi * Q_vl / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

                            S_ql = S_ql_Qvl
                            S_qi_addit = min(S_qi - S_qi_Qvl, Q_vi - Q_vl)  # we should still be able to consume some more ice as long as there was more room in the original ice source, (should be up to Q_vi - S_ql_Qvl which is total growth potential minus current accounted for growth)
                            S_qi = S_qi_Qvl + S_qi_addit # ice growth limited by vapor supersat over liquid and simultaneous liquid growth
                        end # if not > Q_vl, it's def not greater than Q_vi so we don't need to do anything else I think
                    else
                        if S_tot > Q_vi # this handles S_ql whether or not it's larger than Q_vl while simultaneously handling the first part of S_qi
                            S_ql_Qvi = S_ql * Q_vi / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
                            S_qi_Qvi = S_qi * Q_vi / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

                            S_qi = S_qi_Qvi
                            S_ql_addit = min(S_ql - S_ql_Qvi, Q_vl - Q_vi)  # we should still be able to consume some more ice as long as there was more room in the original ice source, (should be up to Q_vi - S_ql_Qvl which is total growth potential minus current accounted for growth)
                            S_ql = S_qi_Qvi + S_ql_addit # ice growth limited by vapor supersat over liquid and simultaneous liquid growth
                        end # if not > Q_vl, it's def not greater than Q_vi so we don't need to do anything else I think
                    end

                elseif (S_ql > 0) && (S_qi < 0) # liquid growth but ice depletion (this can't happen I thnk because we'd be supersaturated over ice if we are over liquid)
                    if T < thermo_params.T_freeze
                        error(
                            "This should not happen (though I guess it can above freezing but that will be handled in limiters later?)",
                        ) # This is just a check that all is working properly, remove later
                    else
                        # the ice is sublimating anyway lol idk...
                    end

                elseif (S_ql < 0) && (S_qi > 0) # ice growth but liquid depletion ( WBF)

                    # we go equally until until either liquid or ice reaches saturation, then the remainder is a balance at that saturation
                    S_ql = -min(-S_ql, q.liq / Δt)  # don't let liquid deplete more than exists (smaller absolute value)
                    S_tot = S_ql + S_qi
                    q_vap_after = q_vap + S_tot * Δt

                    if Q_vl < Q_vi
                        if q_eq.ice <= q_vap_after <= q_eq.liq # we stay between the equilibrium saturation range
                        # all good
                        elseif q_eq.ice > q_vap_after # we can't evaporate enough liquid to grow all the ice it wants to grow, so we scale both down so that we land on the ice saturation line
                            scaling = (q_eq.ice - q_vap) / (S_tot * Δt) # should be a positive number < 1 (S_tot should be negative, and  q_eq.ice - q_vap should be negative)
                            S_qi = S_qi * scaling
                            S_ql = S_ql * scaling
                        elseif q_vap_after > q_eq.liq # we can't grow enough ice to support the liquid we want to evaporate
                            scaling = (q_eq.liq - q_vap) / (S_tot * Δt) # should be a positive number < 1 (S_tot should be positive, and q_eq.liq - q_vap should be positive)
                            S_qi *= scaling
                            S_ql *= scaling
                        end
                    else
                        if q_eq.liq <= q_vap_after <= q_eq.ice # we stay between the equilibrium saturation range
                        # all good
                        elseif q_eq.liq > q_vap_after # we can't evaporate enough liquid to grow all the ice it wants to grow, so we scale both down so that we land on the ice saturation line
                            scaling = (q_eq.liq - q_vap) / (S_tot * Δt) # should be a positive number < 1 (S_tot should be negative, and  q_eq.ice - q_vap should be negative)
                            S_qi = S_qi * scaling
                            S_ql = S_ql * scaling
                        elseif q_vap_after > q_eq.ice # we can't grow enough ice to support the liquid we want to evaporate
                            scaling = (q_eq.ice - q_vap) / (S_tot * Δt) # should be a positive number < 1 (S_tot should be positive, and q_eq.liq - q_vap should be positive)
                            S_qi *= scaling
                            S_ql *= scaling
                        end
                    end

                elseif (S_ql < 0) && (S_qi < 0) # both are depleting (limit by condensate amount and enough to return to saturation)
                    S_ql = -min(-S_ql, q.liq / Δt)  # don't let liquid deplete more than exists (smaller absolute value)
                    S_qi = -min(-S_qi, q.ice / Δt)  # don't let ice deplete more than exists (smaller absolute value)

                    S_ql = -min(-S_ql, -Q_vl)  # don't let liquid deplete more than the vapor can support (smaller absolute value)
                    S_qi = -min(-S_qi, -Q_vi)  # don't let ice deplete more than the vapor can support (smaller absolute value)

                    S_tot = S_ql + S_qi
                    # consume liquid and ice in proportion to their source magnitudes and see if we reach ice saturation

                    if Q_vl < Q_vi # should mean we're below freezing and can evaporate more liquid than we can ice 
                        if S_tot < Q_vi # we would create more vapor than the air wants
                            # go combined until reaching ice saturation, then remainder is liquid evaporation only... (with no re-condensation of ice)
                            S_ql_Qvi = S_ql * Q_vi / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
                            S_qi_Qvi = S_qi * Q_vi / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

                            S_qi = S_qi_Qvi # ice until ice reaching ice saturation vapor pressure... (ice T > freezing, ice is actually more )
                            S_ql_addit = max(S_ql - S_ql_Qvi, Q_vl - Q_vi)  # we should be able to consume some more liquid up until we exhaust the original liquid source (S_ql - S_ql_Qvi), or reach the liquid saturation vapor pressure (Q_vl - Q_vi), whichever is smaller (max of 2 negative numbers)
                            S_ql = S_ql_Qvi + S_ql_addit # liquid growth limited by vapor supersat over liquid and simultaneous ice growth

                            if S_ql_addit > 0
                                error("this shouldnt happen lol, should be max of two negative numbers")
                            end

                        end # if not < Q_vi, we dont even reach ice saturation so we don't need to do anything else
                    else # this is technically weird bc all the ice should evaporate anyway when we're above freezing, but we'll handle anyway... (actually it's slightly different due to difference between thermo_params.T_triple and thermo_params.T_freeze which is 0.01 K but good to handle properly...)
                        if S_tot < Q_vl # we would create more vapor than the air wants
                            # go combined until reaching ice saturation, then remainder is liquid evaporation only... (with no re-condensation of ice)
                            S_ql_Qvl = S_ql * Q_vl / S_tot # liquid growth limited by vapor supersat over liquid and simultaneous ice growth (S_ql_Qvl + S_qi_Qvl = Q_vl)
                            S_qi_Qvl = S_qi * Q_vl / S_tot # ice growth limited by vapor supersat over liquid and simultaneous liquid growth

                            S_ql = S_ql_Qvl # ice until ice reaching ice saturation vapor pressure... (ice T > freezing, ice is actually more )
                            S_qi_addit = max(S_qi - S_qi_Qvl, Q_vi - Q_vl)  # we should be able to consume some more liquid up until we exhaust the original liquid source (S_ql - S_ql_Qvi), or reach the liquid saturation vapor pressure (Q_vl - Q_vi), whichever is smaller (max of 2 negative numbers)
                            S_qi = S_qi_Qvl + S_qi_addit # liquid growth limited by vapor supersat over liquid and simultaneous ice growth

                            if S_qi_addit > 0
                                error("this shouldnt happen lol, should be max of two negative numbers")
                            end

                        end # 

                    end

                end # otherwise we have the negative limiters below for sublimation, evaporation... if both are neg that's sufficient....

            elseif τ_use == :morrison_milbrandt_2015_style
                # this *shouldn't* need limiters the way it's defined but... lol we'll see...
                S_ql, S_qi =
                    morrison_milbrandt_2015_style(
                        param_set,
                        area,
                        ρ,
                        p,
                        T,
                        w,
                        τ_liq,
                        τ_ice,
                        q_vap,
                        q,
                        q_eq,
                        Δt,
                        ts)

                    # Still needs limiters (don't consume liquid that isn't there mainly...)
            elseif τ_use == :morrison_milbrandt_2015_style_exponential_part_only
                # this *shouldn't* need limiters the way it's defined but... lol we'll see...
                S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(
                    param_set,
                    area,
                    ρ,
                    T,
                    w,
                    τ_liq,
                    τ_ice,
                    q_vap,
                    q_eq,
                    Δt,
                )
            else
                error("Unsupported τ_use: $(τ_use)")
            end



            #### REALLY WE SHOULD ONLY MELT ENOUGH TO KNOCK THE EMPERATURE BACK DOWN TO freezing
            #### REALLY WE SHOULD ONLY FREEZE ENOUGH TO KNOCK THE TEMPERATURE BACK UP ABOVE HOMOGENOUS 

            # These would mess with the limiters (because the limiters are saturation bound but these are not) so we do them after... we assume they're agnostic to the antecedent saturation vapor pressures (but no hypsometric correction so can overshoot... but at least they work in one timestep...)
            # Also we only add the remainder after the sources that already exist (debatable the order of operations with the new saturation vapor pressures after the timestep but...)
            if T >= thermo_params.T_freeze # could go insisde supersat bc base eq is already 0 above freezng
                S_qi = min(0, S_qi) # allow sublimation but not new ice formation ( or should we just melt and then let evaporation work? ), i think this is better than redirecting bc that could go on forever with a different timescale, and ice doesn't form above freezing and then melt
                # send existing ice to liquid
                S_ql += (q.ice / Δt + S_qi) # send any existing ice that is not sublimating to liquid (maybe make this a rate later)
                S_qi = -q.ice / Δt # destroy all existing ice

            # homogeneous_freezing
            elseif (T < thermo_params.T_icenuc)
                # Here we allow the V -> L --> I pathway to happen, since it's probably smoother idk and mayb still be physical (form liquid then instafreee unlike ice forming above freezing and then melting...)
                if (S_ql > 0)
                    S_qi += S_ql # any vapor coming to liquid goes to ice instead (smoother in total condensate than setting it to 0 suddenly?)
                    S_ql = 0
                end
                if q.liq > 0
                    S_qi += (q.liq / Δt + S_ql) # send any existing liquid to ice that isn't already evaporating (if S_ql was positive, it's already been sent to ice and set to 0)
                    S_ql = -q.liq / Δt # remove all liquid
                end
            # Heterogeneous ice nucleation
            else
                if get_isbits_nt(param_set.user_args, :use_heterogeneous_ice_nucleation, false)
                    heterogeneous_ice_nuclation_coefficient = get_isbits_nt(param_set.user_aux, :heterogeneous_ice_nuclation_coefficient, FT(1))
                    heterogeneous_ice_nuclation_exponent = get_isbits_nt(param_set.user_aux, :heterogeneous_ice_nuclation_exponent, FT(1))

                    c_1 = ρ * heterogeneous_ice_nuclation_coefficient * exp(heterogeneous_ice_nuclation_exponent * (T - thermo_params.T_freeze) - 1)
                    heterogeneous_ice_nuclation =  q.liq * ( 1 - exp(- c_1 * Δt)) # positive number, how much liquid is losing


                    # we have two exponentials so if we'd deplete all our liquid we just scale down.
                    if S_ql < 0
                        heterogeneous_ice_nuclation = max( -q.liq/Δt - S_ql, -heterogeneous_ice_nuclation ) # S_ql should be smaller than q.liq/Δt after doing limiters above, but maybe move them down here?
                    end
                    S_ql -= heterogeneous_ice_nuclation
                    S_qi += heterogeneous_ice_nuclation
                end
            end





            # # This limiter I believe is wrong, we should be limited by supersaturation, not by the vapor amount and this lets fast timescales do bizarre things and jump to states that condense or create a ton of papper that lead to unrealizeable temperatures...
            # # Too much vapor consumption (note in a supersaturation context, that's a little bit ridiculous to worry about since supersaturation is usually single digit percent at most and timesteps are short-ish)
            # if S > ( Qv - min(0,S_ql) - min(0,S_qi)  ) # only add positive sources to vapor (i.e. subtract if negative S_q)
            #     if (S_qi > 0) && (S_ql > 0)
            #             S_ql *= Qv/S
            #             S_qi *= Qv/S
            #     elseif (S_qi > 0) && (S_ql < 0)
            #         S_qi *= (Qv - S_ql)/S # source to ice not to exceed vapor plus addition from liquid... (S=S_qi here) # (are these stable?, theyre potentially big leaps)
            #     elseif (S_qi < 0) && (S_ql > 0)
            #         S_ql *= (Qv - S_qi)/S # source to liquid not to exceed vapor plus addition from ice... (S=S_ql here) # (are these stable?, theyre potentially big leaps)
            #     end # otherwise we have the negative limiters below for sublimation, evaporation... if both are neg that's sufficient....
            # end



        else # let other limiters do their thing... (this is only acting on sat adjust)
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
        end



        # debug
        # S_ql *= 0. 
        # S_qi *= 0.

        ql_tendency += S_ql
        qi_tendency += S_qi
        # if area < 1e-2
        # @info( "current", Qv, S_ql, S_qi, q.liq / Δt, q.ice / Δt, T, area)
        # end
    end
    return NoneqMoistureSources{FT}(ql_tendency, qi_tendency)
end

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function precipitation_formation(
    param_set::APS,
    precip_model::AbstractPrecipitationModel,
    rain_formation_model::AbstractRainFormationModel,
    qr::FT,
    qs::FT,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts,
    precip_fraction,
) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)

    microphys_params = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    qt_tendency = FT(0)
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    qr_tendency = FT(0)
    qs_tendency = FT(0)
    θ_liq_ice_tendency = FT(0)

    α_acnv = TCP.microph_scaling_acnv(param_set)
    α_accr = TCP.microph_scaling_accr(param_set)

    if area > 0

        q = TD.PhasePartition(thermo_params, ts)

        Π_m = TD.exner(thermo_params, ts)
        c_pm = TD.cp_m(thermo_params, ts)
        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)
        I_i = TD.internal_energy_ice(thermo_params, ts)
        I = TD.internal_energy(thermo_params, ts)

        if precip_model isa Clima0M
            qsat = TD.q_vap_saturation(thermo_params, ts)
            λ = TD.liquid_fraction(thermo_params, ts)

            S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(microphys_params, q, qsat))

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1 - λ)
            qt_tendency += S_qt
            ql_tendency += S_qt * λ
            qi_tendency += S_qt * (1 - λ)
            θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))
        end

        if precip_model isa Clima1M
            T = TD.air_temperature(thermo_params, ts)
            T_fr = TCP.T_freeze(param_set)
            c_vl = TCP.cv_l(param_set)
            c_vm = TD.cv_m(thermo_params, ts)
            Rm = TD.gas_constant_air(thermo_params, ts)
            Lf = TD.latent_heat_fusion(thermo_params, ts)

            # TODO - limiters and positivity checks should be done elsewhere
            qr = max(qr, FT(0)) / precip_fraction
            qs = max(qs, FT(0)) / precip_fraction

            # Autoconversion of cloud ice to snow is done with a simplified rate.
            # The saturation adjustment scheme prevents using the
            # 1-moment snow autoconversion rate that assumes
            # that the supersaturation is present in the domain.
            if rain_formation_model isa Clima1M_default
                S_qt_rain = -min(q.liq / Δt, α_acnv * CM1.conv_q_liq_to_q_rai(microphys_params, q.liq))
            elseif rain_formation_model isa Clima2M
                S_qt_rain =
                    -min(
                        q.liq / Δt,
                        α_acnv * CM2.conv_q_liq_to_q_rai(
                            microphys_params,
                            rain_formation_model.type,
                            q.liq,
                            ρ,
                            N_d = rain_formation_model.prescribed_Nd,
                        ),
                    )
            else
                error("Unrecognized rain formation model")
            end

            use_supersat = get_isbits_nt(param_set.user_args, :use_supersat, false)
            if use_supersat == false # it's either true or a specified value
                S_qt_snow = -min(q.ice / Δt, α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice))
            else
                p::FT = TD.air_pressure(thermo_params, ts)
                S_qt_snow = -min(q.ice / Δt, α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p, use_supersat))
            end
            
            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            ql_tendency += S_qt_rain
            qi_tendency += S_qt_snow

            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

            # accretion cloud water + rain
            if rain_formation_model isa Clima1M_default
                S_qr =
                    min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ)) *
                    precip_fraction
            elseif rain_formation_model isa Clima2M
                if rain_formation_model.type isa CMT.LD2004Type
                    S_qr =
                        min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ)) *
                        precip_fraction
                elseif rain_formation_model.type isa CMT.KK2000Type || rain_formation_model.type isa CMT.B1994Type
                    S_qr =
                        min(
                            q.liq / Δt,
                            α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr, ρ),
                        ) * precip_fraction
                elseif rain_formation_model.type isa CMT.TC1980Type
                    S_qr =
                        min(
                            q.liq / Δt,
                            α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr),
                        ) * precip_fraction
                else
                    error("Unrecognized 2-moment rain formation model type")
                end
            else
                error("Unrecognized rain formation model")
            end
            qr_tendency += S_qr
            qt_tendency -= S_qr
            ql_tendency -= S_qr
            θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

            # accretion cloud ice + snow
            S_qs =
                min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, snow_type, q.ice, qs, ρ)) *
                precip_fraction
            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

            # sink of cloud water via accretion cloud water + snow
            S_qt =
                -min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, snow_type, q.liq, qs, ρ)) *
                precip_fraction
            if T < T_fr # cloud droplets freeze to become snow)
                qs_tendency -= S_qt
                qt_tendency += S_qt
                ql_tendency += S_qt
                θ_liq_ice_tendency -= S_qt / Π_m / c_pm * Lf * (1 + Rm / c_vm)
            else # snow melts, both cloud water and snow become rain
                α::FT = c_vl / Lf * (T - T_fr)
                qt_tendency += S_qt
                ql_tendency += S_qt
                qs_tendency += S_qt * α
                qr_tendency -= S_qt * (1 + α)
                θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1 + Rm / c_vm) * α - L_v0)
            end

            # sink of cloud ice via accretion cloud ice - rain
            S_qt =
                -min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, rain_type, q.ice, qr, ρ)) *
                precip_fraction
            # sink of rain via accretion cloud ice - rain
            S_qr = -min(qr / Δt, α_accr * CM1.accretion_rain_sink(microphys_params, q.ice, qr, ρ)) * precip_fraction
            qt_tendency += S_qt
            qi_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)

            # accretion rain - snow
            if T < T_fr
                S_qs =
                    min(qr / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, snow_type, rain_type, qs, qr, ρ)) *
                    precip_fraction *
                    precip_fraction
            else
                S_qs =
                    -min(qs / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, rain_type, snow_type, qr, qs, ρ)) *
                    precip_fraction *
                    precip_fraction
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
        end
    end
    return PrecipFormation{FT}(θ_liq_ice_tendency, qt_tendency, ql_tendency, qi_tendency, qr_tendency, qs_tendency)
end
