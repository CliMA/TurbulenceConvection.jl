"""
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

include("microphysics_coupling_limiters.jl")

"""
for dispatch backwards compat. we've moved to caching saturation sepcific humidities, so in general those methods can just accept an argument.
"""
function noneq_moisture_sources(param_set::APS,
    noneq_sources_type::AbstractNonEquillibriumSourcesType,
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    w::FT,
    z::FT
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
)  where {FT}

    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    
    T::FT = TD.air_temperature(thermo_params, ts)

    q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
    q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())

    return noneq_moisture_sources(param_set, noneq_sources_type, moisture_sources_limiter, area, ρ, Δt, ts, w, z, q_vap_sat_liq, q_vap_sat_ice)
end

# -------------------------------------------------------------------------------------------------------------------------- #

function noneq_moisture_sources(
    param_set::APS,
    ::RelaxToEquilibrium,
    moisture_sources_limiter::Union{NoMoistureSourcesLimiter, BasicMoistureSourcesLimiter},
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    w::FT,
    z::FT,
    q_vap_sat_liq::FT,
    q_vap_sat_ice::FT,
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
)  where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    if area > 0
        # basic noneq (no supersat formulation so not likely to be right) should be a RelaxToEquilibrium()
        q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
        T::FT = TD.air_temperature(thermo_params, ts)
        p::FT = TD.air_pressure(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

        # TODO - is that the state we want to be relaxing to?
        ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ, T, q.tot)
        q_eq = TD.PhasePartition(thermo_params, ts_eq)

        S_ql = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, q)
        S_qi = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, ice_type, q_eq, q)
        S_ql, S_qi = calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, FT(0), FT(0), FT(0), FT(0),FT(0), q_vap, q, q_eq, Δt, ts, S_ql, S_qi)
        ql_tendency += S_ql
        qi_tendency += S_qi
    end
    return NoneqMoistureSources{FT}(ql_tendency, qi_tendency)
end

function noneq_moisture_sources(
    param_set::APS,
    noneq_moisture_scheme::AbstractRelaxationTimescaleType,
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    w::FT,
    z::FT,
    q_vap_sat_liq::FT,
    q_vap_sat_ice::FT,
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    if area > 0

        q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
        T::FT = TD.air_temperature(thermo_params, ts)
        p::FT = TD.air_pressure(thermo_params, ts) # really this should come from aux_gm.p_c but I don't want to change the form of non_eq_moisture_sources for backwards compat...
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

        # use phase partition in case we wanna use the conv_q_vap fcn but maybe not best for supersat since it's not really a phase partition (all 3 are vapor amounts)

        τ_liq, τ_ice = get_τ(param_set, microphys_params, noneq_moisture_scheme, q, T, p, ρ, w, z)

        q_eq = TD.PhasePartition(
            q.tot,
            q_vap_sat_liq,
            q_vap_sat_ice,
            # TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
            # TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
        ) # all 3 are vapor amounts
        S_ql, S_qi = calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
        ql_tendency += S_ql
        qi_tendency += S_qi
    end
    return NoneqMoistureSources{FT}(ql_tendency, qi_tendency)
end


function noneq_moisture_sources(
    param_set::APS,
    ::KorolevMazin2007,
    moisture_sources_limiter::Union{BasicMoistureSourcesLimiter, StandardSupersaturationMoistureSourcesLimiter},
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    w::FT,
    z::FT,
    q_vap_sat_liq::FT,
    q_vap_sat_ice::FT,
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
)  where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    if area > 0
        q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
        T::FT = TD.air_temperature(thermo_params, ts)
        p::FT = TD.air_pressure(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)
        
        S_ql, S_qi = korolev_mazin_2007(param_set, area, ρ, Δt, ts, w)
        q_eq = TD.PhasePartition(
            q.tot,
            q_vap_sat_liq,
            q_vap_sat_ice,
            # TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
            # TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
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

    S_ql, S_qi = calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, FT(0), FT(0), FT(0), FT(0),FT(0), q_vap, q, q_eq, Δt, ts,  S_ql, S_qi)
    ql_tendency += S_ql
    qi_tendency += S_qi    
    end
    return NoneqMoistureSources{FT}(ql_tendency, qi_tendency)
end

# -------------------------------------------------------------------------------------------------------------------------- #

"""
Allow for other microphysics processes to be handled outside of the main microphysics noneq calls
This is probably good bc we had it in noneq_moisture_sources but it's conflated w/ cond/evap sub/dep then and the limiters we had didn't account for precip tendencies anyway...

For consistency I've kept passing in S_ql, S_qi here but maybe we should just deprecate that and just have one limiter at the end? idk...
"""
function other_microphysics_processes(
    param_set::APS,
    heterogeneous_ice_nucleation::Tuple{Bool, FT, FT},
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    w::FT,
    z::FT,
    S_ql::FT,
    S_qi::FT
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
) where {FT}

    S_ql_init::FT = S_ql
    S_qi_init::FT = S_qi

    # S_ql = FT(0) # for backwards compat, but in the end should just remove this altogether... this might be unstable though bc it has no relaxation so might run up against limiters? doesn't account for incoming noneq fluxes... (old way didn't account for precip but that's within one condensate species at least)
    # S_qi = FT(0) # for backwards compat, but in the end should just remove this altogether... this might be unstable though bc it has no relaxation so might run up against limiters? doesn't account for incoming noneq fluxes... (old way didn't account for precip but that's within one condensate species at least)
    # if it becomes a problem we could try also passing state and pulling out the sources from there... but we'll have to know if we're in env or up etc...

    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    qi_tendency_homogeneous_freezing::FT = FT(0)
    qi_tendency_heterogeneous_icenuc::FT = FT(0)
    qi_tendency_heterogeneous_freezing::FT = FT(0)
    qi_tendency_icemult::FT = FT(0)
    qi_tendency_melting::FT = FT(0)
    if area > 0
        q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
        T::FT = TD.air_temperature(thermo_params, ts)
        p::FT = TD.air_pressure(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

        if T >= thermo_params.T_freeze # could go insisde supersat bc base eq is already 0 above freezng
            S_qi = min(0, S_qi) # allow sublimation but not new ice formation ( or should we just melt and then let evaporation work? ), i think this is better than redirecting bc that could go on forever with a different timescale, and ice doesn't form above freezing and then melt
            # send existing ice to liquid
            qi_tendency_melting = (q.ice / Δt + S_qi) # send any existing ice that is not sublimating to liquid (maybe make this a rate later)
            S_ql += qi_tendency_melting # send any existing ice that is not sublimating to liquid (maybe make this a rate later)
            S_qi = -q.ice / Δt # destroy all existing ice  # limiters will maybe fight against this... could also make some very large number? idk.. if this fails in one timestep, the next should take care of it... as long as those sources are aware of temp regimes?)

        # homogeneous_freezing
        elseif (T < thermo_params.T_icenuc)
            # Here we allow the V -> L --> I pathway to happen, since it's probably smoother idk and mayb still be physical (form liquid then instafreee unlike ice forming above freezing and then melting...)
            if (S_ql > 0)
                qi_tendency_homogeneous_freezing = S_ql # any vapor coming to liquid goes to ice instead (smoother in total condensate than setting it to 0 suddenly?)
                S_qi += qi_tendency_homogeneous_freezing # any vapor coming to liquid goes to ice instead (smoother in total condensate than setting it to 0 suddenly?)
                S_ql = 0
            end
            if q.liq > 0
                qi_tendency_homogeneous_freezing = q.liq / Δt + S_ql
                S_qi += qi_tendency_homogeneous_freezing # send any existing liquid to ice that isn't already evaporating (if S_ql was positive, it's already been sent to ice and set to 0)
                S_ql = -q.liq / Δt # remove all liquid (is this the best way to do this? limiters will maybe fight against this... could also make some very large number? idk.. if this fails in one timestep, the next should take care of it... maybe?)
            end
            # Heterogeneous ice nucleation
        else
            if heterogeneous_ice_nucleation[1]
                heterogeneous_ice_nuclation_coefficient, heterogeneous_ice_nuclation_exponent = heterogeneous_ice_nucleation[2:3]
                c_1 =
                    ρ *
                    heterogeneous_ice_nuclation_coefficient *
                    exp(heterogeneous_ice_nuclation_exponent * (T - thermo_params.T_freeze) - 1)
                qi_tendency_heterogeneous_icenuc = q.liq * (1 - exp(-c_1 * Δt)) # positive number, how much liquid is losing and ice is gaining


                # we have two exponentials so if we'd deplete all our liquid we just scale down.
                if S_ql < 0
                    # qi_tendency_heterogeneous_icenuc = max(-q.liq / Δt - S_ql, -qi_tendency_heterogeneous_icenuc) # S_ql should be smaller than q.liq/Δt after doing limiters above, but maybe move them down here? [[ i think this is wrong..., should be -max(...) ]]
                    # qi_tendency_heterogeneous_icenuc = min(q.liq / Δt + S_ql, qi_tendency_heterogeneous_icenuc) # S_ql should be smaller than q.liq/Δt after doing limiters above, but maybe move them down here?
                    qi_tendency_heterogeneous_icenuc = -limit_tendency(moisture_sources_limiter, -qi_tendency_heterogeneous_icenuc, q.liq + S_ql * Δt, Δt)
                end
                S_ql -= qi_tendency_heterogeneous_icenuc
                S_qi += qi_tendency_heterogeneous_icenuc
            end
        end
    end

    return OtherMicrophysicsSources{FT}(
        S_ql - S_ql_init,
        S_qi - S_qi_init,
        # my additions for storage
        qi_tendency_homogeneous_freezing,
        qi_tendency_heterogeneous_freezing,
        qi_tendency_heterogeneous_icenuc,
        qi_tendency_melting,
    )
end


# -------------------------------------------------------------------------------------------------------------------------- #

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function precipitation_formation(
    param_set::APS,
    moisture_model::AbstractMoistureModel,
    precip_model::AbstractPrecipitationModel,
    cloud_sedimentation_model::AbstractCloudSedimentationModel,
    rain_formation_model::AbstractRainFormationModel,
    snow_formation_model::AbstractSnowFormationModel,
    qr::FT,
    qs::FT,
    ql::FT,
    qi::FT,
    Ni::FT,
    term_vel_ice::FT, # maybe this isn't needed lol
    term_vel_rain::FT,
    term_vel_snow::FT,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    w::FT,
    z::FT,
    qi_tendency_sub_dep::FT,
    qi_tendency_sed::FT,
    precip_fraction::FT,
    precipitation_tendency_limiter::AbstractTendencyLimiter,
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


    ql_tendency_acnv = FT(0)
    qi_tendency_acnv = FT(0)
    qi_tendency_acnv_dep = FT(0) # is an input bc it's calculated right after noneq_moisture_sources
    qi_tendency_acnv_agg = FT(0)
    ql_tendency_accr_liq_ice = FT(0)
    ql_tendency_accr_liq_rai = FT(0)
    ql_tendency_accr_liq_sno = FT(0)
    qi_tendency_accr_ice_liq = FT(0)
    qi_tendency_accr_ice_rai = FT(0)
    qi_tendency_accr_ice_sno = FT(0)

    α_acnv = TCP.microph_scaling_acnv(param_set)
    α_accr = TCP.microph_scaling_accr(param_set)


    

    if area > 0
        q = TD.PhasePartition(thermo_params, ts)
        # ql = q.liq
        # qi = q.ice

        Π_m = TD.exner(thermo_params, ts)
        c_pm = TD.cp_m(thermo_params, ts)
        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)
        I_i = TD.internal_energy_ice(thermo_params, ts)
        I = TD.internal_energy(thermo_params, ts)

        if precip_model isa Clima0M
            qsat = TD.q_vap_saturation(thermo_params, ts)
            λ = TD.liquid_fraction(thermo_params, ts)

            # S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(microphys_params, q, qsat))
            S_qt = limit_tendency(precipitation_tendency_limiter, CM0.remove_precipitation(microphys_params, q, qsat), (q.liq + q.ice), Δt)

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
                # S_qt_rain = -min(q.liq / Δt, α_acnv * CM1.conv_q_liq_to_q_rai(microphys_params, q.liq))
                
                S_qt_rain = limit_tendency(precipitation_tendency_limiter, -α_acnv * CM1.conv_q_liq_to_q_rai(microphys_params, q.liq), max(q.liq - get_q_threshold(param_set, CMT.LiquidType()), 0), Δt) # don't allow reducing below threshold
                # should we have a noneq version here? that uses my_conv_q_liq_to_q_rai()
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

            ql_tendency_acnv = S_qt_rain # for storage

            if snow_formation_model isa NonEquilibriumSnowFormationModel # should also be a NonEquilibriumMoisture model
                # if moisture_model isa NonEquilibriumMoisture # should be true
                q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
                T::FT = TD.air_temperature(thermo_params, ts)
                p::FT = TD.air_pressure(thermo_params, ts)
                # q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)
                # τ_liq, τ_ice = get_τ(param_set, microphys_params, moisture_model.scheme, q, T, p, ρ, w, z)

                # S_qt_snow_ice_dep = limit_tendency(precipitation_tendency_limiter, -α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p; τ = τ_ice), q.ice, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_dep = qi_tendency_acnv_dep # the qi tendencency is some fraction of positive deposition, and then this is a loss from that to snow

                # [ a lot of `agg` is also just threshold driven convergence from things like sed, so also don't count sed against this let it go straight to snow if needed and can evap there ]

                # not sure if we should only use positive sedimentation convergence... i think so bc losing ice doesnt make things smaller
                # qi_tendency_sed =  max(FT(0), qi_tendency_sed)
                qi_tendency_sed =  resolve_nan(max(FT(0), qi_tendency_sed), FT(0)) # Only accept deposition, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))
                # qi_tendency_sed = FT(0)

                # qi_tendency_sub_dep = max(FT(0), qi_tendency_sub_dep) # this is the part that is new ice formation, so we can use that to limit the snow formation
                qi_tendency_sub_dep = resolve_nan(max(FT(0), qi_tendency_sub_dep), FT(0)) # Only accept deposition, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))

                # techincally acnv really can deplete existing ice too.. so add taht to limiter. include sed just bc
                S_qt_snow_ice_dep = limit_tendency(precipitation_tendency_limiter, -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor * my_conv_q_ice_to_q_sno(param_set, q, T, p; N = Ni), q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_dep = -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor * qi_tendency_sub_dep # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]


                # This check may not be true (i.e. in principle you can really acnv more than the new ice formation in thier (not morrison 05) framework -- this check also encourages boosting acnv scaling w/o penalty which may lead to unrealistic outcomes from balanced parameters.
                # S_qt_snow_ice_dep = max(-max(0, qi_tendency_sub_dep), S_qt_snow_ice_dep) # Now that they're calculated differently, ensure S_qt_snow_ice_dep takes lss than deposition from qi_tendency_noneq
                


                # the limiter should probably pick one to give priority, dep i would assume

                # my_conv_q_ice_to_q_sno_thresh() is really just dispatching to my_conv_q_ice_to_q_sno_no_supersat attm, but multiplying by (r/r_ice_snow)^3


                #=
                get minimum q threshold for snow formation
                    path of calls is my_conv_q_ice_to_q_sno_thresh -> my_conv_q_ice_to_q_sno_no_supersat -> get_q_threshold_acnv(param_set, CMT.IceType(), q, T, p; N=N, assume_N_is = assume_N_is)
                    this is kinda fragile though... the target threshold is kinda buried in there...
                    and if it's a q based threshold, does it even matter? maybe we just fall back to the default min threshold regardless? either way it would be nice to get out the threshold the fcn is using...
                =#
                    
                S_qt_snow_ice_agg, q_ice_thresh = get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set, q, T, p; N = Ni)

                S_qt_snow_ice_agg = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_agg , max(q.ice - q_ice_thresh, 0), Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_agg = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_agg , max(q.ice - q_ice_thresh, 0) + qi_tendency_sed * Δt, Δt)  # I think allowing sed to join the party is a mistake bc we do not calculate acnv including that so you still get jumps. i.e. you hit the tresh and then all of a sudden sed ca take you wherever you want, you'd need to add it when calculating the thresh, as in the input being q.ice + sed tendency etc... but that get's more and more complicated, you could make the argument with more and more tendencies... etc... this is simple and in line wit everything else...


                # S_qt_snow_ice_agg = limit_tendency(precipitation_tendency_limiter, -α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p; N = Ni, assume_N_is=true), q.ice + qi_tendency_sed * Δt, Δt) # based on aggregation but threshold grows with N
                # S_qt_snow_ice_agg = limit_tendency(precipitation_tendency_limiter, -α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p;), q.ice, Δt) # based on aggregation [ trialing thresh not being N depdendent q_0 small anyway, use same threshold as growth then.] otherwise we need another parameter for the stable q,N based on r_is

                # S_qt_snow_ice_agg = FT(0) # nothing works and arguably i think it's pretty small? the RF09 is prolly an illusion... the other spikes seem fake too. maybe there's room for one that scales w/ q but i think it's just super small? Most of it i think is very large r driven and we just don't have that...





                
                # not sure if we should only use positive sedimentation convergence...
                S_qt_snow = limit_tendency(precipitation_tendency_limiter, S_qt_snow_ice_dep + S_qt_snow_ice_agg, q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth and aggregation # technically the part from sub_dep is new ice so add that in to the limiter (net bc the acnv tendency is neg of the sub tendency)
                    
                # if S_qt_snow > (S_qt_snow_ice_dep + S_qt_snow_ice_agg)
                #     error("S_qt_snow < S_qt_snow_ice_dep + S_qt_snow_ice_agg, this should not happen; S_qt_snow = $(S_qt_snow); S_qt_snow_ice_dep = $(S_qt_snow_ice_dep); S_qt_snow_ice_agg = $(S_qt_snow_ice_agg); qi_tendency_sed = $(qi_tendency_sed); q_i = $(q.ice); Δt = $(Δt);")
                # end


                if any(!isfinite, (S_qt_snow, S_qt_snow_ice_dep, S_qt_snow_ice_agg))
                    error("S_qt_snow or S_qt_snow_ice_dep or S_qt_snow_ice_agg is not finite; S_qt_snow = $(S_qt_snow); S_qt_snow_ice_dep = $(S_qt_snow_ice_dep); S_qt_snow_ice_agg = $(S_qt_snow_ice_agg); qi_tendency_sub_dep = $(qi_tendency_sub_dep); qi_tendency_sed = $(qi_tendency_sed); q_i = $(q.ice); Δt = $(Δt); q = $(q); T = $(T); p = $(p);")
                end

                # S_qt_snow = -min(q.ice / Δt, α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p))
                # S_qt_snow = limit_tendency(precipitation_tendency_limiter, -α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p; N0 = Ni), q.ice, Δt)
            else
                #=
                NOTE: NonEquilibriumSnowFormationModel was maybe a bad naming convention -- it's not 1:1 with NonEquilibriumMoistureModel... even EquilibriumMoistureModel can use this...
                So this is really just a fallback method...
                =#
                # S_qt_snow = -min(q.ice / Δt, α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice))
                S_qt_snow = limit_tendency(precipitation_tendency_limiter, -α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice), q.ice, Δt)
                S_qt_snow_ice_dep = FT(NaN) # no easy breakdown
                S_qt_snow_ice_agg = FT(NaN) # no easy breakdown
            end

            qi_tendency_acnv = S_qt_snow # for storage
            qi_tendency_acnv_agg = S_qt_snow_ice_agg # for storage
            qi_tendency_acnv_dep = S_qt_snow_ice_dep # for storage [ was already an input but may be removing ]

            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            ql_tendency += S_qt_rain
            qi_tendency += S_qt_snow

            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

            # accretion cloud water + rain
            if rain_formation_model isa Clima1M_default
                # S_qr = min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ)) * precip_fraction
                S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ), q.liq, Δt) * precip_fraction
            elseif rain_formation_model isa Clima2M
                if rain_formation_model.type isa CMT.LD2004Type
                    # S_qr = min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ)) * precip_fraction
                    S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ), q.liq, Δt) * precip_fraction
                elseif rain_formation_model.type isa CMT.KK2000Type || rain_formation_model.type isa CMT.B1994Type
                    # S_qr = min(q.liq / Δt, α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr, ρ), ) * precip_fraction
                    S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr, ρ), q.liq, Δt) * precip_fraction
                elseif rain_formation_model.type isa CMT.TC1980Type
                    # S_qr = min(q.liq / Δt, α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr),) * precip_fraction
                    S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr), q.liq, Δt) * precip_fraction
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
            ql_tendency_accr_liq_rai = -S_qr # for storage

            # accretion cloud ice + snow
            # S_qs = min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, snow_type, q.ice, qs, ρ)) * precip_fraction
            S_qs = -limit_tendency(precipitation_tendency_limiter,  -α_accr * CM1.accretion(microphys_params, ice_type, snow_type, q.ice, qs, ρ), q.ice, Δt) * precip_fraction

            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0
            qi_tendency_accr_ice_sno = -S_qs # for storage

            # sink of cloud water via accretion cloud water + snow
            # S_qt = -min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, snow_type, q.liq, qs, ρ)) * precip_fraction
            S_qt = limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, liq_type, snow_type, q.liq, qs, ρ), q.liq, Δt) * precip_fraction
            if T < T_fr # cloud droplets freeze to become snow)
                qs_tendency -= S_qt
                qt_tendency += S_qt
                ql_tendency += S_qt
                θ_liq_ice_tendency -= S_qt / Π_m / c_pm * Lf * (1 + Rm / c_vm)
                ql_tendency_accr_liq_sno = S_qt # for storage
            else # snow melts, both cloud water and snow become rain
                α = c_vl / Lf * (T - T_fr)
                qt_tendency += S_qt
                ql_tendency += S_qt
                qs_tendency += S_qt * α
                qr_tendency -= S_qt * (1 + α)
                θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1 + Rm / c_vm) * α - L_v0)
                ql_tendency_accr_liq_sno = S_qt # for storage
            end


            # sink of cloud ice via accretion cloud ice - rain
            # S_qt = -min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, rain_type, q.ice, qr, ρ)) * precip_fraction
            S_qt = limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, ice_type, rain_type, q.ice, qr, ρ), q.ice, Δt) * precip_fraction
            # sink of rain via accretion cloud ice - rain
            # S_qr = -min(qr / Δt, α_accr * CM1.accretion_rain_sink(microphys_params, q.ice, qr, ρ)) * precip_fraction
            S_qr = limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion_rain_sink(microphys_params, q.ice, qr, ρ), qr, Δt) * precip_fraction
            qt_tendency += S_qt
            qi_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s0)
            qi_tendency_accr_ice_rai = S_qt # for storage

            # accretion rain - snow
            if T < T_fr
                # S_qs = min(qr / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, snow_type, rain_type, qs, qr, ρ)) * precip_fraction * precip_fraction
                S_qs = -limit_tendency(precipitation_tendency_limiter, -α_accr * my_accretion_snow_rain(microphys_params, snow_type, rain_type, qs, qr, ρ, term_vel_snow, term_vel_rain), qr, Δt) * precip_fraction * precip_fraction
            else
                # S_qs = -min(qs / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, rain_type, snow_type, qr, qs, ρ)) * precip_fraction * precip_fraction
                S_qs = limit_tendency(precipitation_tendency_limiter, -α_accr * my_accretion_snow_rain(microphys_params, rain_type, snow_type, qr, qs, ρ, term_vel_rain, term_vel_snow), qs, Δt) * precip_fraction * precip_fraction
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm

            # liq - ice accretion

            if T < T_fr # losing liq
                S_qi = accretion_liq_ice(microphys_params, ql, qi, ρ, term_vel_ice, cloud_sedimentation_model.ice_terminal_velocity_scheme, Ni, cloud_sedimentation_model.ice_Dmax, cloud_sedimentation_model.liq_ice_collision_efficiency; liq_ice_accr_scaling_factor = cloud_sedimentation_model.liq_ice_collision_scaling_factor ) # should be a positive number
                # Tendencies get multiplied by area outside of this fcn and cloud fraction is 1 if has condensate and 0 if not so no need here
                S_qi = max(S_qi, FT(0)) # don't allow negative values (depending on ai, bi, ci in Chen that can just barely happen)
                S_qi = -limit_tendency(precipitation_tendency_limiter, -α_accr * S_qi , ql, Δt) # signs should match accretion_snow_rain except snow ⟹ ice, rain ⟹ liq
            else # losing ice
                # S_qi = accretion_liq_ice(microphys_params, ql, qi, ρ, term_vel_ice, cloud_sedimentation_model.ice_terminal_velocity_scheme, Ni, cloud_sedimentation_model.ice_Dmax, cloud_sedimentation_model.liq_ice_collision_efficiency; liq_ice_accr_scaling_factor = cloud_sedimentation_model.liq_ice_collision_scaling_factor ) # should be a positive number
                # S_qi = max(S_qi, FT(0)) # don't allow negative values (depending on ai, bi, ci  Chen that can just barely happen)
                # S_qi = limit_tendency(precipitation_tendency_limiter, -α_accr * S_qi , qi, Δt) # signs should match accretion_snow_rain except snow ⟹ ice, rain ⟹ liq
                S_qi = FT(0) # PSACWI doesn't cover this, either way we have instant melting so...
            end
            qi_tendency += S_qi
            ql_tendency -= S_qi   
            # Because both ql and qi are part of the working fluid and θ_liq_ice_tendency is 0, we really should have  no change in θ_liq_ice_tendency
            # θ_liq_ice_tendency += FT(0)
            # θ_liq_ice_tendency += S_qi * Lf / Π_m / c_vm # fusion liq ⟹ ice
            # if !iszero(S_qi)
            #     @info "S_qi = $S_qi; qi = $qi; ql = $ql; Δt = $Δt; ρ = $ρ; term_vel_ice = $term_vel_ice; Ni = $Ni; θ_liq_ice_tendency_addit = $θ_liq_ice_tendency_addit"
            # end
            qi_tendency_accr_ice_liq = S_qi # for storage
            ql_tendency_accr_liq_ice = -S_qi # for storage
        end
    end
    return PrecipFormation{FT}(
        θ_liq_ice_tendency,
        qt_tendency,
        ql_tendency,
        qi_tendency,
        qr_tendency,
        qs_tendency,
        # my additions
        ql_tendency_acnv,
        qi_tendency_acnv,
        qi_tendency_acnv_dep,
        qi_tendency_acnv_agg,
        ql_tendency_accr_liq_rai,
        # FT(0), # ql_tendency_accr_liq_ice
        ql_tendency_accr_liq_ice,
        ql_tendency_accr_liq_sno,
        # FT(0), # qi_tendency_accr_ice_liq,
        qi_tendency_accr_ice_liq,
        qi_tendency_accr_ice_rai,
        qi_tendency_accr_ice_sno,
    )
end

 precipitation_formation(
    param_set::APS,
    moisture_model::AbstractMoistureModel,
    precip_model::AbstractPrecipitationModel,
    cloud_sedimentation_model::AbstractCloudSedimentationModel,
    rain_formation_model::AbstractRainFormationModel,
    snow_formation_model::AbstractSnowFormationModel,
    qr::FT,
    qs::FT,
    ql::FT,
    qi::FT,
    Ni::FT,
    term_vel_ice::FT, # maybe this isn't needed lol
    term_vel_rain::FT,
    term_vel_snow::FT,
    area::FT,
    ρ::FT,
    Δt::Real,
    ts::TD.ThermodynamicState,
    precip_fraction::FT,
    precipitation_tendency_limiter::AbstractTendencyLimiter
) where {FT} = 
    precipitation_formation(
        param_set,
        moisture_model,
        precip_model,
        cloud_sedimentation_model,
        rain_formation_model,
        snow_formation_model,
        qr,
        qs,
        ql,
        qi,
        Ni,
        term_vel_ice,
        term_vel_rain,
        term_vel_snow,
        area,
        ρ,
        Δt,
        ts,
        FT(NaN), # placeholder for backwards compat w
        FT(NaN), # placeholder for backwards compat z
        FT(NaN), # placeholder for backwards compat qi_tendency_sub_dep
        FT(NaN), # placeholder for backwards compat qi_tendency_sed
        precip_fraction,
        precipitation_tendency_limiter,
    )


"""
    We'd copy from accretion_snow_rain() in Microphysics1M.jl bc it has framework for dueling terminal velocities we don't have all parameters defined there.

    Liq also doesn't sediment rn so the implementation is more like accretion(). Then, the source is to the moving one (qi) by default.
    But, like accretion_rain_snow(), we change that based on T < T_fr.

    accretion_liq_sno(prs, type_i, type_j, q_i, q_j, ρ)

    - `i` - snow for temperatures below freezing
            or rain for temperatures above freezing
    - `j` - rain for temperatures below freezing
            or snow for temperatures above freezing
    - `prs` - abstract set with Earth parameters
    - `type_i`, `type_j` - a type for snow or rain
    - `q_` - specific humidity of snow or rain
    - `ρ` - air density

    Returns the accretion rate between rain and snow.
    Collisions between liq and ice result in ice at temperatures below freezing and in liq at temperatures above freezing.

    As written, using already calculated terminal velocities, nothing here depends on the terminal velocity scheme for Blk1MVel or Chen2022, though if you didn't have that precomputed you would need it.
    This is just because they both assume the same form for m(r)
        
"""
function accretion_liq_ice(
    prs::ACMP,
    ql::FT, # stationary
    qi::FT, # moving
    ρ::FT,
    term_vel_ice::FT,
    velo_scheme::Union{CMT.Chen2022Type, CMT.Blk1MVelType},
    N_i::FT,
    ice_Dmax::FT,
    E_liq_ice::FT;
    liq_ice_accr_scaling_factor::FT = FT(1),
    # D_transition::FT = FT(0.625e-3), # .625mm is the transition from small to large ice crystals in Chen paper, leave unchanged
) where {FT <: Real}

    accr_rate = FT(0)
    if (ql > FT(0) && qi > FT(0))
        
        integral_nav = int_nav_dr(prs, ice_type, velo_scheme, qi, ρ, FT(0); Nt=N_i, Dmax=ice_Dmax) # assume μ in n(r) is 0 (exponential distribution)

        # what do we do if nm_integral_value is inf?

        # if isnan(integral_nav) || !isfinite(integral_nav)
        #     error("Got nm_integral_value = $integral_nav from nm_integral() in accretion_liq_ice() for inputs qi = $qi; ρ = $ρ; N_i = $N_i; ice_Dmax = $ice_Dmax")
        # end

        # if isnan(integral_nav) || !isfinite(integral_nav)
        #     error("Got integral_nav = $integral_nav from term_vel_ice * integral_nav in accretion_liq_ice() for inputs qi = $qi; ρ = $ρ; N_i = $N_i; ice_Dmax = $ice_Dmax; term_vel_ice = $term_vel_ice; nm_integral_value = $nm_integral_value")
        # end
        
        accr_rate = E_liq_ice * ql * integral_nav * liq_ice_accr_scaling_factor  # Eq 16 https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/#Accretion

        # if isnan(accr_rate) || !isfinite(accr_rate)
        #     error("Got accr_rate = $accr_rate from E_liq_ice * ql * int in accretion_liq_ice() for inputs qi = $qi; ρ = $ρ; N_i = $N_i; ice_Dmax = $ice_Dmax; term_vel_ice = $term_vel_ice; nm_integral_value = $nm_integral_value; E_liq_ice = $E_liq_ice; ql = $ql")
        # end

    end
    return accr_rate 
end


"""
My version of this function that uses the true terminal_velocity() [stored or recalculated]
"""
function my_accretion_snow_rain(
    prs::ACMP,
    type_i::CMT.AbstractPrecipType,
    type_j::CMT.AbstractPrecipType,
    q_i::FT,
    q_j::FT,
    ρ::FT,
    term_vel_i::FT, # ordinarily CloudMicrophysics.jl would recalculate this but we already have it stored in TC.jl
    term_vel_j::FT, # ordinarily CloudMicrophysics.jl would recalculate this but we already have it stored in TC.jl
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        _n0_i::FT = CM1.n0(prs, q_i, ρ, type_i)
        _n0_j::FT = CM1.n0(prs, q_j, ρ, type_j)

        _r0_j::FT = CM1.r0(prs, type_j)
        _m0_j::FT = CM1.m0(prs, type_j)
        _me_j::FT = CM1.me(prs, type_j)
        _Δm_j::FT = CM1.Δm(prs, type_j)
        _χm_j::FT = CM1.χm(prs, type_j)

        _E_ij::FT = CM1.E(prs, type_i, type_j)

        _λ_i::FT = CM1.lambda(prs, type_i, q_i, ρ)
        _λ_j::FT = CM1.lambda(prs, type_j, q_j, ρ)

        # _v_ti = terminal_velocity(prs, type_i, CT.Blk1MVelType(), ρ, q_i)
        # _v_tj = terminal_velocity(prs, type_j, CT.Blk1MVelType(), ρ, q_j)
        _v_ti = term_vel_i
        _v_tj = term_vel_j

        accr_rate =
            FT(π) / ρ *
            _n0_i *
            _n0_j *
            _m0_j *
            _χm_j *
            _E_ij *
            abs(_v_ti - _v_tj) / _r0_j^(_me_j + _Δm_j) * (
                FT(2) * CM1.SF.gamma(_me_j + _Δm_j + FT(1)) / _λ_i^FT(3) /
                _λ_j^(_me_j + _Δm_j + FT(1)) +
                FT(2) * CM1.SF.gamma(_me_j + _Δm_j + FT(2)) / _λ_i^FT(2) /
                _λ_j^(_me_j + _Δm_j + FT(2)) +
                CM1.SF.gamma(_me_j + _Δm_j + FT(3)) / _λ_i /
                _λ_j^(_me_j + _Δm_j + FT(3))
            )
    end
    # return accr_rate # we shouldnt need resolve nan bc we're using the cm_1 values which shouldn't blow up? but just to be safe
    return resolve_nan(accr_rate, FT(0))
end