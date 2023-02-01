"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function noneq_moisture_sources(param_set::APS, area::FT, ρ::FT, Δt::Real, ts, w; ts_LCL=nothing) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep Δt
    ql_tendency = FT(0)
    qi_tendency = FT(0)
    if area > 0
        use_supersat      = get(param_set.user_args, :use_supersat, false) # (:use_supersat in keys(param_set.user_args)) ? param_set.user_args.use_supersat : false # so we dont have to set everything we dont know is in user_args in the defaults...
        use_korolev_mazin = get(param_set.user_args, :use_korolev_mazin, false) # (:use_supersat in keys(param_set.user_args)) ? param_set.user_args.use_supersat : false # so we dont have to set everything we dont know is in user_args in the defaults...
        raymond_ice_test  = get(param_set.user_args, :raymond_ice_test, false) # the timescale parameterization from raymond...
        N_r_closure       = get(param_set.user_args, :N_r_closure, :monodisperse)

        q = TD.PhasePartition(thermo_params, ts)
        T = TD.air_temperature(thermo_params, ts)
        q_vap = TD.vapor_specific_humidity(thermo_params, ts)

        if use_supersat # use phase partition in case we wanna use the conv_q_vap fcn but maybe not best for supersat since it's not really a phase partition (all 3 are vapor amounts)
            # nonmutable param_set so would have to edit tau usage everywhere -- let's just do in supersat
            τ_liq = CMNe.τ_relax(microphys_params, liq_type)
            τ_ice = CMNe.τ_relax(microphys_params, ice_type)
            if raymond_ice_test
                # problems: 
                # - POTENTIAL CCN SHOULD GO UP W/ SUPERSATURATION (for a given q this also means up with decreasing temp?)
                # - INP CONC HERE COULD GROW W/O BOUND, WE HAVE NO METRIC FOR GROWING CCN  conc
                # - if you trust something like https://doi.org/10.5194/acp-14-81-2014 then INP/CDNC_i concentrations are way lower than liq but droplet sizes way higher... how to handle?
                # either way our CCN concentratinos seems to be far too large for our INP numbers to give relevant taus until far too cold...
                # gets even worse cause CDNC,CCN etc are advected in updraft so... it's not just local...
                # if we assume CDCN/CCN/INP is constant ish and it's just activation fraction of CCN/INP that changes with T and/or SS, then what...

                N_0 = get(param_set.user_args, :N_0, FT(100*10e6))  # per m^3 both ccn and N_d in this range https://doi.org/10.1029/2020MS002205 (or  FT(1) # params tuned to fit this better... )
                N_m = get(param_set.user_args, :N_m, FT(-0.2)) # log slope https://doi.org/10.1073/pnas.1514034112
                N_b = get(param_set.user_args, :N_b, FT(-5 - 273.15 * N_m)) # -5 - 273.15 * N_m

                N_INP = 10^(N_m*T + N_b) * 10^3 # per liter to per m^3
                f_INP = N_INP/(N_0+N_INP) # clamp( N_INP/(N_0+N_INP), 0, 1)  # let N_0 and N_INP vary freely...

                # -- hopefully this gets around the problem of draining water vapor at initialization of clouds but also allows speedup as droplets grow (assuming fixed drop concenctration)
                D, ρ_w = FT(0.0000226), FT(1) 
                # R    =  max(((q.liq + q.ice)/(4/3*π*ρ_w*N_0))^(1/3), FT(0.2*10^-6)) # bound to be at least ~micron size...something like kohler crit radius

                # @show(ts_LCL)
                if N_r_closure == :inhomogeneous
                    _,R_liq = NR_inhomogeneous_mixing_liquid(thermo_params, N_0, TD.air_pressure(thermo_params,ts), q.liq, ts_LCL) # testing inhomogeneous mixing would have r fixed and then let N vary... set r based on adiabatic extrapolation from cloud base 
                    _,R_ice = NR_monodisperse(N_INP,q.ice)
                elseif N_r_closure == :monodisperse # uniform size for all droplets, liquid and Ice I guess
                    _,R_liq = NR_monodisperse(N_0  ,q.liq)
                    _,R_ice = NR_monodisperse(N_INP,q.ice)
                else
                    error("Unsupported size distribution closure (N_r_closure): $(N_r_closure)")
                end

                base = 1/(4*π*D) # as q goes up, R goes 
                τ_liq = base / (N_0 * R_liq)
                τ_ice = base / (N_INP * R_ice)
                # @show((T, q, τ_liq, τ_ice))

                # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster...
                # println("effective τ_liq = ",τ_liq_eff, " effective τ_ice = ",τ_ice_eff, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )
                if     (0 < τ_liq) && (τ_liq < 1) # fast source
                    # @show("rate limiting cond: ", S_ql, (q_vap - q_eq.liq) / 1)
                    τ_liq = FT(1)
                elseif (0 < τ_ice) && (τ_ice < 1) # fast source
                    # @show("rate limiting dep: ", S_qi, (q_vap - q_eq.ice) / 1)
                    τ_ice = FT(1)
                elseif (-1 < τ_liq) && (τ_liq < 0) # fast sink
                    # @show("rate limiting evap: ", S_ql, (q_vap - q_eq.liq) / 1)
                    τ_liq = FT(-1)
                elseif (-1 < τ_ice) && (τ_ice < 0) # fast sink
                    # @show("rate limiting sub: ", S_qi, (q_vap - q_eq.ice) / 1)
                    τ_ice = FT(-1)
                else
                end

            end

            q_eq = TD.PhasePartition(q.tot, TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()), TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())) # all 3 are vapor amounts
            S_ql = (q_vap - q_eq.liq) / τ_liq # | microphys_params.τ_cond_evap | CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, TD.PhasePartition(FT(0),q_vap,FT(0)))  
            S_qi = (q_vap - q_eq.ice) / τ_ice # -(source to vapor) = source to condensate
        elseif use_korolev_mazin
            # need to get w into here somewhere...
            S_ql,S_qi = korolev_mazin_2007(param_set, area, ρ, Δt, ts, w)
            q_eq = TD.PhasePartition(q.tot, TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()), TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())) # all 3 are vapor amounts
            # println("effective τ_liq = ",(q_vap - q_eq.liq)/S_ql, " effective τ_ice = ",(q_vap - q_eq.ice)/S_qi, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )

            # korolev_mazin fix effective tau > 1 (or whatever timestep was realy but we run w/ 1 mostly rn)
            τ_liq_eff = (q_vap - q_eq.liq)/S_ql 
            τ_ice_eff = (q_vap - q_eq.ice)/S_qi
            # println("effective τ_liq = ",τ_liq_eff, " effective τ_ice = ",τ_ice_eff, " T = ", T, " w = ",w, " q = ",q, " q_eq = ", q_eq )
            if     (0 < τ_liq_eff) && (τ_liq_eff < 1) # fast source
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
            
        else
            # TODO - is that the state we want to be relaxing to?
            ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q.tot)
            q_eq = TD.PhasePartition(thermo_params, ts_eq)

            S_ql = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, q) 
            S_qi = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, ice_type, q_eq, q)
        end

        if T >= thermo_params.T_freeze # could go insisde supersat bc base eq is already 0 above freezng
            S_qi = min(0,S_qi) # allow sublimation but not new ice formation ( or should we just melt and then let evaporation work? )
            
            # send existing ice to liquid
            S_qi -= q.ice/Δt # send any existing ice to liquid
            S_ql += q.ice/Δt # send any existing ice to liquid
        end

        if (use_supersat || use_korolev_mazin) # might need to do these first bc ql,qc tendencies maybe are applied individually and can still crash the code if one is too large...
            # print("limiting")
            # add homogenous freezing (first so can limit later)s
            if (T < thermo_params.T_icenuc)
                if (S_ql > 0)
                    S_qi += S_ql # any vapor coming to liquid goes to ice instead (smoother than setting it to 0 suddenly?)
                    S_ql = 0
                end
                if q.liq > 0
                    S_ql -= q.liq/Δt # send any existing liquid to ice (maybe make this a rate later)
                    S_qi += q.liq/Δt # send any existing liquid to ice (maybe make this a rate later)
                end
            end
        
            # limiter (maybe we dont need this cause of our later limiters?) was useful when we did this first....
            # could be bad in wbf if ice drain on vapor is too large but liquid is supposed to be going to vapor but less so 
            # S = S_ql + S_qi
            S = max(0, S_ql) + max(0,S_qi) # only add if positive sources to condensate
            Qv = q_vap / Δt
            # if S > Qv # not enough vapor to create condensate (reverse are both handled already)
            if S > ( Qv - min(0,S_ql) - min(0,S_qi)  ) # only add if positive source to vapor (subtracat the negative)
                if (S_qi > 0) && (S_ql > 0)
                        S_ql *= Qv/S
                        S_qi *= Qv/S
                elseif (S_qi > 0) && (S_ql < 0)
                    S_qi *= (Qv + S_ql)/S # source to ice not to exceed vapor plus addition from liquid... (S=S_qi here)
                elseif (S_qi < 0) && (S_ql > 0)
                    S_ql *= (Qv + S_qi)/S
                end # otherwise we have the negative limiters below for sublimation, evaporation... if both are neg that's sufficient....
                # S_ql *= Qv/S
                # S_qi *= Qv/S
            end
        end

        # TODO - handle limiters elswhere
        if S_ql >= FT(0)
            S_ql = min(S_ql, q_vap / Δt) # source to liquid not more than current vapor
        else
            S_ql = -min(-S_ql, q.liq / Δt) # loss of liquid not to exceed current liquid
        end
        if S_qi >= FT(0)
            S_qi = min(S_qi, q_vap / Δt)
        else
            S_qi = -min(-S_qi, q.ice / Δt)
        end

        # numerical_mismatch = q_vap - S_ql*Δt - S_qi*Δt
        # if numerical_mismatch < 0 # there would be negative vapor (fix doesnt fix numeral instabilities tho so maybe exclude, is an FT error)
        #     if S_ql  > 0
        #         S_ql += numerical_mismatch
        #     else
        #         S_qi += numerical_mismatch
        #     end
        # end

        ql_tendency += S_ql
        qi_tendency += S_qi
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
            S_qt_snow = -min(q.ice / Δt, α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice))
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
