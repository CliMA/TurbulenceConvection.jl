"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function noneq_moisture_sources(param_set::APS, area::FT, ρ::FT, Δt::Real, ts, p::FT) where {FT}

    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    if area > 0

        q = TD.PhasePartition(param_set, ts)
        T = TD.air_temperature(param_set, ts)
        q_vap = TD.vapor_specific_humidity(param_set, ts) 

        supersat_formulation = true
        # supersat_formulation = false

        if ! supersat_formulation
            # This seems to match the sat adjust one the best (ρTq not as good)
            ts_eq = TD.PhaseEquil_pTq(param_set, p, T, q.tot) # pressure not density
            q_eq  = TD.PhasePartition(param_set, ts_eq)

            # calculate the saturation adustment state here... and print both out (and maybe diff?)
            # θ_liq_ice = TD.liquid_ice_pottemp(param_set, ts) 
            # ts_eq_sat_adj = TD.PhaseEquil_pθq(param_set, p, θ_liq_ice, q.tot) # true eq is sat adjustmen phaseequil p theta q (rho, theta_liq_ice, rho q_tot/rho)
            # q_eq_sat_adj  = TD.PhasePartition(param_set, ts_eq_sat_adj)
            # ts_eq = ts_eq_sat_adj # test if this is the issue between them... ( matches very well w/ this one from testing so far )
            # q_eq  = q_eq_sat_adj

            # This is from CloudMicrophysics.jl::MicrophysicsNonEq.jl
            S_ql = CMNe.conv_q_vap_to_q_liq_ice(param_set, liq_type, q_eq, q)
            S_qi = CMNe.conv_q_vap_to_q_liq_ice(param_set, ice_type, q_eq, q)

        else #if supersat_formulation:
            _T_freeze =    ICP.T_freeze(param_set) # this one defined locally, do we need some sort of lower bound different rate once we get down to ice_nuc temperature? Forming a lot of liquid instead of ice up high rn...
            _T_icenuc = TD.ICP.T_icenuc(param_set) # this one currently only exists in TD, should be 233K
            # test here a native supersat formulation? then eq would be... somehow need to partition the source additions based on the relative supersats and taus? kind of a roundabout way to do that... but timestepping is elsewhere so..        
            q_vap_sat_liq = TD.q_vap_saturation_liquid( param_set, ts )
            q_vap_sat_ice = TD.q_vap_saturation_ice(param_set, ts )

            τ_liq = CMNe.τ_relax(param_set, liq_type) #  τ_cond_evap(params_set)
            τ_ice = CMNe.τ_relax(param_set, ice_type) #  we also would have our turbulence adjusted tau, q_liq from field, like zhang 2019, but also are there vertical velocity distributions at which sublimation vs evap happens as in Storelvmo 2010, re Korolev 2003?
            # I think this is susceptible to, because we're not using a state w/ consistent energy/ theta_liq_ice, problems w/ stability...
            S_ql = (q_vap - q_vap_sat_liq) / τ_liq # should we add the 'psychrometric correction to account for the release of latent heat'?
            S_qi = (q_vap - q_vap_sat_ice) / τ_ice
            if (T > _T_freeze) && (S_qi > 0) # not sure if this is necessary but i saw some small amounts of ice above 0 and ion think this function is range limited so i'ma test that jawn out
                # @info "limiting the ice source: `S\_i: $(S_qi)` | T: `$(T)`"
                S_qi = 0 # no gaining of ice, allows loss though I feel like there should be melting above freezing not just sub... should the timescale also be T dependent? idk.
            end
            if T < _T_icenuc
                if S_ql > 0 # T < T_icenuc and T > T_freeze should be mutually exclusive but break up if loop just incase rather than use elseif
                    @info "limiting liquid formation: S\\_l: `$(S_ql)` | T: `$(T)`"
                    # stop liquid formation
                    S_qi += S_ql # not sure if just to get rid of it or add it to ice and presumably this shouldn't be an absolute heaviside since it's not a phase transition value? idk... seems to end up in liquid still idk why
                    S_ql = 0 # do existing liquid drops autoconvert to ice as well? WBF here should sweep them up but idk if they should autoconvert faster if they're mixed/transported above that level....
                end
                if q.liq > 0
                    # do we need homogenous freezing? this works but then the rates look all weird in the cloud w/ jumps in ice...
                    S_ql = -q.liq/Δt - S_ql # limit it so if we already have some negative liquid they add up 
                    S_qi =  q.liq/Δt
                end
            end

            # liquid being transported by the updraft seems to be the main problem, w/ no freezing mech besides WBF, idk why env liquid has a spike tho...


            # proportionally, these can't add up to more than q_vap, we should probably limit one but for now let's just rescale
            S = S_ql + S_qi
            S_max = q_vap / Δt 
            # If S > 0, this will have the effect of limiting vapor depletion
            # If S < 0, i.e. we are getting vapor from our condensates, the limiters below should still hold fine and this won't get triggered because the ratio will be negative.
            if S > S_max
                S_scale = S_max/S # If this works we'll move the scaling below up to here and clean this jawn up, still never gets triggered but at least it's not setting things to 1 that should be neg etc...
            else
                S_scale = 1
            end
            
            S_ql = S_ql * S_scale
            S_qi = S_qi * S_scale

        end

        # TODO - handle limiters elswhere
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

        ql_tendency += S_ql
        qi_tendency += S_qi

        # ## add liq <--> ice transition (could also be a partition on the sources above based on the relative amounts of q.liq, q.ice, q_vap)
        # # τ_frz_mlt = ICP.τ_frz_mlt(param_set)
        # # τ_frz_mlt = param_set.τ_frz_mlt
        # τ_frz_mlt(params_set) # right way to call this.
        # τ_frz_mlt = 1e-5

        # # split into liq/ice separately, keep only liq to ice at first and use 2 timescales, ensure conservation
        # # @info("limiting") # doesnt work in repl? idk works in jobs
        # S_ql_to_qi = (q.liq - q_eq.liq)/τ_frz_mlt - (q.ice - q_eq.ice)/τ_frz_mlt # remove excess liquid and gain liq from excess ice
    
        # # TODO - handle limiters elswhere
        # if S_ql_to_qi >= FT(0)
        #     S_ql_to_qi =  min( S_ql_to_qi, q.liq / Δt + ql_tendency) # a positive flux from liq to ice is limited by liq, and we gotta include the tendency we already utilized
        # else
        #     S_ql_to_qi = -min(-S_ql_to_qi, q.ice / Δt + qi_tendency) # a negative flux to liq is limitied by ice
        # end

        # ql_tendency -= S_ql_to_qi # these must be equal and opposite for conservation
        # qi_tendency += S_ql_to_qi

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
    qr::FT,
    qs::FT,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts,
    precip_fraction,
) where {FT}

    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    qt_tendency = FT(0)
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    qr_tendency = FT(0)
    qs_tendency = FT(0)
    θ_liq_ice_tendency = FT(0)

    if area > 0

        q = TD.PhasePartition(param_set, ts)

        Π_m = TD.exner(param_set, ts)
        c_pm = TD.cp_m(param_set, ts)
        L_v0 = ICP.LH_v0(param_set)
        L_s0 = ICP.LH_s0(param_set)

        if precip_model isa Clima0M
            qsat = TD.q_vap_saturation(param_set, ts)
            λ = TD.liquid_fraction(param_set, ts)

            S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(param_set, q, qsat))

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1 - λ)
            qt_tendency += S_qt
            ql_tendency += S_qt * λ
            qi_tendency += S_qt * (1 - λ)
            θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))
        end

        if precip_model isa Clima1M
            T = TD.air_temperature(param_set, ts)
            T_fr = ICP.T_freeze(param_set)
            c_vl = ICP.cv_l(param_set)
            c_vm = TD.cv_m(param_set, ts)
            Rm = TD.gas_constant_air(param_set, ts)
            Lf = TD.latent_heat_fusion(param_set, ts)

            qr = qr / precip_fraction
            qs = qs / precip_fraction

            # Autoconversion of cloud ice to snow is done with a simplified rate.
            # The saturation adjustment scheme prevents using the
            # 1-moment snow autoconversion rate that assumes
            # that the supersaturation is present in the domain.
            S_qt_rain = -min(q.liq / Δt, CM1.conv_q_liq_to_q_rai(param_set, q.liq))
            S_qt_snow = -min(q.ice / Δt, CM1.conv_q_ice_to_q_sno_no_supersat(param_set, q.ice))
            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            ql_tendency += S_qt_rain
            qi_tendency += S_qt_snow
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

            # accretion cloud water + rain
            S_qr = min(q.liq / Δt, CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ)) * precip_fraction
            qr_tendency += S_qr
            qt_tendency -= S_qr
            ql_tendency -= S_qr
            θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

            # accretion cloud ice + snow
            S_qs = min(q.ice / Δt, CM1.accretion(param_set, ice_type, snow_type, q.ice, qs, ρ)) * precip_fraction
            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

            # sink of cloud water via accretion cloud water + snow
            S_qt = -min(q.liq / Δt, CM1.accretion(param_set, liq_type, snow_type, q.liq, qs, ρ)) * precip_fraction
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
            S_qt = -min(q.ice / Δt, CM1.accretion(param_set, ice_type, rain_type, q.ice, qr, ρ)) * precip_fraction
            # sink of rain via accretion cloud ice - rain
            S_qr = -min(qr / Δt, CM1.accretion_rain_sink(param_set, q.ice, qr, ρ)) * precip_fraction
            qt_tendency += S_qt
            qi_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)

            # accretion rain - snow
            if T < T_fr
                S_qs =
                    min(qr / Δt, CM1.accretion_snow_rain(param_set, snow_type, rain_type, qs, qr, ρ)) *
                    precip_fraction *
                    precip_fraction
            else
                S_qs =
                    -min(qs / Δt, CM1.accretion_snow_rain(param_set, rain_type, snow_type, qr, qs, ρ)) *
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
