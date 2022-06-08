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

        # TODO - is that the state we want to be relaxing to?
        ts_eq = TD.PhaseEquil_ρTq(param_set, ρ, T, q.tot)
        q_eq = TD.PhasePartition(param_set, ts_eq)
        #alt
        ts_eq_alt = TD.PhaseEquil_pTq(param_set, p, T, q.tot)
        q_eq_alt = TD.PhasePartition(param_set, ts_eq_alt)
        # calculate the saturation adustment state here... and print both out (and maybe diff?)
        θ_liq_ice = TD.liquid_ice_pottemp(param_set, ts) 
        ts_eq_sat_adj = TD.PhaseEquil_ρθq(param_set, ρ, θ_liq_ice, q.tot) # true eq is sat adjustmen phaseequil rho theta q (rho, theta_liq_ice, rhoq_ot/rho)
        q_eq_sat_adj = TD.PhasePartition(param_set, ts_eq_sat_adj)

        # println("==================================================================")
        # @info "q_eq noneq `$(q_eq)`"
        # @info "q_eq noneq alt `$(q_eq_alt)`"
        # @info "q_eq sat adjust `$(q_eq_sat_adj)`"
        # @info "q existing `$(q)`"
        # println("-----------------")
        # @info "T noneq `$(TD.air_temperature(param_set, ts_eq))` | T noneq alt p `$(TD.air_temperature(param_set, ts_eq_alt))` | T sat adjust `$(TD.air_temperature(param_set, ts_eq_sat_adj))` | T existing `$(T)`"
        # println("-----------------")
        # @info "θ\\_liq\\_ice noneq `$(TD.liquid_ice_pottemp(param_set, ts_eq))` | θ\\_liq\\_ice noneq alt p `$(TD.liquid_ice_pottemp(param_set, ts_eq_alt))` | θ\\_liq\\_ice sat adjust `$(TD.liquid_ice_pottemp(param_set, ts_eq_sat_adj))` | θ\\_liq\\_ice existing `$(TD.liquid_ice_pottemp(param_set, ts))`"

        ts_eq = ts_eq_alt # use the alt formluation to run
        q_eq  = q_eq_alt

        S_ql = CMNe.conv_q_vap_to_q_liq_ice(param_set, liq_type, q_eq, q)
        S_qi = CMNe.conv_q_vap_to_q_liq_ice(param_set, ice_type, q_eq, q)

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
