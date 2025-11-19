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

include("microphysics_coupling_helpers.jl") # some definitions
include("microphysics_coupling_limiters.jl")

"""
for dispatch backwards compat. we've moved to caching saturation sepcific humidities, so in general those methods can just accept an argument.
"""
function noneq_moisture_sources(param_set::APS,
    noneq_sources_type::AbstractNonEquillibriumSourcesType,
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    area::FT,
    ρ_c::FT,
    p::FT,
    T::FT,
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    w::FT,
    # z::FT,
    dqvdt::FT = FT(0),
    dTdt::FT = FT(0),
    τ_liq::FT = FT(NaN),
    τ_ice::FT = FT(NaN),
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
)  where {FT}

    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    ρ = TD.air_density(thermo_params, ts)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    
    q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
    q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())

    return noneq_moisture_sources(param_set, noneq_sources_type, moisture_sources_limiter, area, ρ, p, T, Δt, ts, w, q_vap_sat_liq, q_vap_sat_ice, dqvdt, dTdt, τ_liq, τ_ice)
end

# -------------------------------------------------------------------------------------------------------------------------- #

function noneq_moisture_sources(
    param_set::APS,
    ::RelaxToEquilibrium,
    moisture_sources_limiter::Union{NoMoistureSourcesLimiter, BasicMoistureSourcesLimiter},
    area::FT,
    ρ_c::FT,
    p::FT,
    T::FT,
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    w::FT,
    # z::FT,
    q_vap_sat_liq::FT,
    q_vap_sat_ice::FT,
    dqvdt::FT,
    dTdt::FT,
    τ_liq::FT = FT(NaN),
    τ_ice::FT = FT(NaN),
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
        q::TD.PhasePartition{FT} = TD.PhasePartition(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)
        ρ::FT = TD.air_density(thermo_params, ts)

        # TODO - is that the state we want to be relaxing to?
        ts_eq = TD.PhaseEquil_ρTq(thermo_params, ρ_c, T, q.tot)
        q_eq = TD.PhasePartition(thermo_params, ts_eq)

        S_ql = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, q)
        S_qi = CMNe.conv_q_vap_to_q_liq_ice(microphys_params, ice_type, q_eq, q)
        S_ql, S_qi = calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, FT(0), FT(0), FT(0), FT(0), FT(0), q_vap, dqvdt, dTdt, q, q_eq, Δt, ts, S_ql, S_qi)
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
    ρ_c::FT,
    p::FT,
    T::FT,
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    w::FT,
    # z::FT, # deprecate, we only had it for raymond ice and that is no longer used
    q_vap_sat_liq::FT,
    q_vap_sat_ice::FT,
    dqvdt::FT,
    dTdt::FT,
    τ_liq::FT, # not optional here bc it's used
    τ_ice::FT,
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    if area > 0

        q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
        ρ::FT = TD.air_density(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

        # use phase partition in case we wanna use the conv_q_vap fcn but maybe not best for supersat since it's not really a phase partition (all 3 are vapor amounts)
        # τ_liq, τ_ice = get_τs(param_set, microphys_params, noneq_moisture_scheme, q, T, p, ρ, w) # no more z, there's just not valid world in which z matters sine we're not doing raymond ice anymore.

        q_eq = TD.PhasePartition(
            q.tot,
            q_vap_sat_liq,
            q_vap_sat_ice,
        ) # all 3 are vapor amounts
        S_ql, S_qi = calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, q, q_eq, Δt, ts)
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
    ρ_c::FT,
    p::FT, # pass this bc it's not a no-op, unlike ρ and q which are stored in nonequil
    T::FT, # pass this bc it's not a no-op, unlike ρ and q which are stored in nonequil. note pressure is not stored in nonequil
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    w::FT,
    # z::FT,
    q_vap_sat_liq::FT,
    q_vap_sat_ice::FT,
    dqvdt::FT,
    dTdt::FT,
    τ_liq::FT = FT(NaN),
    τ_ice::FT = FT(NaN),
    ;
    # ts_LCL::Union{Nothing, TD.ThermodynamicState} = nothing,
)  where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    # TODO - when using adaptive timestepping we are limiting the source terms with the previous timestep Δt
    ql_tendency::FT = FT(0)
    qi_tendency::FT = FT(0)
    if area > 0
        q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
        ρ::FT = TD.air_density(thermo_params, ts)
        q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)
        
        q_eq = TD.PhasePartition(
            q.tot,
            q_vap_sat_liq,
            q_vap_sat_ice,
        ) # all 3 are vapor amounts

    S_ql, S_qi = calculate_timestep_limited_sources(moisture_sources_limiter, param_set, area, ρ, FT(0), FT(0), FT(0), FT(0),FT(0), q_vap, dqvdt, dTdt, q, q_eq, Δt, ts)
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
    moisture_model::AbstractMoistureModel,
    noneq_moisture_scheme::AbstractRelaxationTimescaleType,
    moisture_sources_limiter::AbstractMoistureSourcesLimiter,
    area::FT,
    ρ_c::FT,
    p::FT,
    T::FT, # pass this bc it's not a no-op, unlike ρ and q which are stored in nonequil
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    w::FT,
    # z::FT, # idek why we had this, i think it just got copied from noneq_moisture_sources() but we don't need it here
    S_ql::FT,
    S_qi::FT;
    N_l::FT = FT(NaN),
    N_i::FT = FT(NaN),
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
        if moisture_model isa NonEquilibriumMoisture
            q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
            q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)
            ρ = TD.air_density(thermo_params, ts)

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
                # liq-ice collisions [[ nucleation only, separate from accretion liq_ice ]]
                heterogeneous_ice_nucleation = moisture_model.heterogeneous_ice_nucleation
                # if heterogeneous_ice_nucleation[1]
                if heterogeneous_ice_nucleation.use_heterogeneous_ice_nucleation
                    # heterogeneous_ice_nucleation_coefficient, heterogeneous_ice_nucleation_exponent = heterogeneous_ice_nucleation[2:3]
                    heterogeneous_ice_nucleation_coefficient = heterogeneous_ice_nucleation.heterogeneous_ice_nucleation_coefficient
                    heterogeneous_ice_nucleation_exponent = heterogeneous_ice_nucleation.heterogeneous_ice_nucleation_exponent
                    c_1 =
                        ρ *
                        heterogeneous_ice_nucleation_coefficient *
                        exp(heterogeneous_ice_nucleation_exponent * (T - thermo_params.T_freeze) - 1)
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

                # direct activation
                S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())

                # if (S_i >= 1.08) || (S_i >= 0.999 && T <= 265.15)
                if (S_i > FT(0)) && (T < TCP.T_freeze(param_set))

                    q_activation = FT(0)
                    N_i_activation = FT(0)
                    N_INP_max = FT(0)

                    r0_activation = FT(10e-6) # 10 μm particle post activation, matches https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L1089
                    m0_activation = particle_mass(param_set, ice_type, r0_activation)


                    # N scales w/ T, do we wanna add an S scaling here? Maybe that's overcomplex...


                    # Cooper curve: max number of INPs per L, convert to per kg and clamp
                    if !(noneq_moisture_scheme isa INP_Aware_Timescale)
                        N_INP_max = get_N_i_Cooper_curve(T; clamp_N=true)
                        # N_INP_max = clamp(N_INP_max / ρ_c, FT(0), FT(500e3) / ρ_c)
                    else
                        N_INP_max = get_INP_concentration(param_set, noneq_moisture_scheme, q, T, ρ_c, w) # get the number of ice nuclei per m³
                    end
                    if isnan(N_i) # just assume q / MI0 (actually that's too small, we'll assume r_is...)
                        # N_i = q.ice / m0_activation
                        r_is = CMP.r_ice_snow(microphys_params)
                        N_i = N_from_qr(param_set, ice_type, q.ice, r_is; monodisperse=false, μ=FT(0), ρ=ρ_c)
                    end

                    # This means within 10 seconds, we will be at 10 micron average, which is about 100x smaller than 30 micron average, regardless of what the later growth rate we calibrated is...
                    # I guess in principle we could have just kept N = N_i always when supersaturated, the problem becomes then you need to set <r>.  (10/0.2)^3 is about 1.25e5, so we would have had to have tau be dependent on <r> or something to facilitate activation... basically assume <r> is never less than 10 micrometer when supersaturated or something like that. but we can just park that responsibility here.
                    if N_INP_max > N_i
                        N_i_activation = (N_INP_max - N_i) / FT(10) # just go w/ 10 seconds, it should be pretty fast... [[ note this means that in practice, tau is very fast at the start to get you to 10 micrometer droplets, hence prolly the jump we see in our LES data plots from the low tau,N at cold to suddenly hella ice. ]]
                        q_activation = N_i_activation * m0_activation / ρ_c
                    end
            
                    qi_tendency_heterogeneous_freezing += q_activation
                    qi_tendency_heterogeneous_freezing = -limit_tendency(moisture_sources_limiter, -qi_tendency_heterogeneous_freezing, max(N_INP_max - N_i, FT(0)) * m0_activation, Δt)
                    S_qi += qi_tendency_heterogeneous_freezing
                end
            end
        else
            # Equilibrium just pass, do nothing...
        end
    end

    return OtherMicrophysicsSources{FT}(
        S_ql - S_ql_init,
        S_qi - S_qi_init,
        # my additions for storage
        qi_tendency_homogeneous_freezing, # homogenous liquid to ice feeezing
        qi_tendency_heterogeneous_freezing, # primary ice nucleation (aerosol freezing)
        qi_tendency_heterogeneous_icenuc, # contact freezing liq_ice
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
    area_partition_model::AbstractAreaPartitionModel,
    qr::FT,
    qs::FT,
    ql::FT,
    qi::FT,
    Ni::FT,
    term_vel_ice::FT, # maybe this isn't needed lol
    term_vel_rain::FT,
    term_vel_snow::FT,
    area::FT,
    ρ_c::FT, # for some reason they always used ρ_c here, i'm not sure if that's good or bad.
    p::FT,
    T::FT,
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    ql_tendency_cond_evap::FT,
    qi_tendency_sub_dep::FT,
    qi_tendency_sed::FT,
    precip_fraction::FT,
    precipitation_tendency_limiter::AbstractTendencyLimiter,
    tke_var::FT,
    N_i_no_boost::FT,
    S_i::FT,
    τ_sub_dep::FT,
    dN_i_dz::FT,
    dqidz::FT,
    N_INP::FT,
    massflux::FT,
    domain::AbstractDomain,
    qt_var::FT,
    dqtdz::FT,
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
    qi_tendency_acnv_dep_is = FT(0)
    qi_tendency_acnv_dep_above = FT(0)
    qi_tendency_acnv_agg_mix = FT(0)
    qi_tendency_acnv_thresh = FT(0)
    ql_tendency_accr_liq_ice = FT(0)
    ql_tendency_accr_liq_rai = FT(0)
    ql_tendency_accr_liq_sno = FT(0)
    qi_tendency_accr_ice_liq = FT(0)
    qi_tendency_accr_ice_rai = FT(0)
    qi_tendency_accr_ice_sno = FT(0)
    #
    qs_tendency_accr_rai_sno = FT(0)

    α_acnv = TCP.microph_scaling_acnv(param_set)
    α_accr = TCP.microph_scaling_accr(param_set)


    

    if area > 0
        q = TD.PhasePartition(thermo_params, ts)
        # ql = q.liq
        # qi = q.ice
        ρ = TD.air_density(thermo_params, ts)

        Π_m = TD.exner(thermo_params, ts)
        c_pm = TD.cp_m(thermo_params, ts)

        I_i = TD.internal_energy_ice(thermo_params, ts)
        I = TD.internal_energy(thermo_params, ts)

        if precip_model isa Clima0M
            qsat = TD.q_vap_saturation(thermo_params, ts)
            λ = TD.liquid_fraction(thermo_params, ts)

            L_v = TD.latent_heat_vapor(thermo_params, ts)
            L_s = TD.latent_heat_sublim(thermo_params, ts)

            # S_qt = -min((q.liq + q.ice) / Δt, -CM0.remove_precipitation(microphys_params, q, qsat))
            S_qt = limit_tendency(precipitation_tendency_limiter, CM0.remove_precipitation(microphys_params, q, qsat), (q.liq + q.ice), Δt)

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1 - λ)
            qt_tendency += S_qt
            ql_tendency += S_qt * λ
            qi_tendency += S_qt * (1 - λ)
            θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v * λ + L_s * (1 - λ))
        end

        if precip_model isa Clima1M
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
                               
                # S_qt_rain = limit_tendency(precipitation_tendency_limiter, -α_acnv * CM1.conv_q_liq_to_q_rai(microphys_params, q.liq), q.liq, Δt)
                ql_tendency_cond_evap = resolve_nan(max(FT(0), ql_tendency_cond_evap), FT(0)) # Only accept evaporation, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))
                S_qt_rain = limit_tendency(precipitation_tendency_limiter, -α_acnv * CM1.conv_q_liq_to_q_rai(microphys_params, q.liq), max(q.liq + ql_tendency_cond_evap*Δt - get_q_threshold(param_set, CMT.LiquidType()), 0), Δt) # don't allow reducing below threshold
                # should we have a noneq version here? that uses my_conv_q_liq_to_q_rai()
            elseif rain_formation_model isa Clima2M
                S_qt_rain =
                    -min(
                        q.liq / Δt,
                        α_acnv * CM2.conv_q_liq_to_q_rai(
                            microphys_params,
                            rain_formation_model.type,
                            q.liq,
                            ρ_c,
                            N_d = rain_formation_model.prescribed_Nd,
                        ),
                    )
            else
                error("Unrecognized rain formation model")
            end

            ql_tendency_acnv = S_qt_rain # for storage

            if (snow_formation_model isa NonEquilibriumSnowFormationModel)
                q::TD.PhasePartition{FT} =TD.PhasePartition(thermo_params, ts)
                # q_vap::FT = TD.vapor_specific_humidity(thermo_params, ts)

                # S_qt_snow_ice_dep = limit_tendency(precipitation_tendency_limiter, -α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q; N = Ni, τ = τ_ice, r_acnv_scaling_factor=snow_formation_model.r_ice_acnv_scaling_factor), q.ice, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_dep = qi_tendency_acnv_dep # the qi tendencency is some fraction of positive deposition, and then this is a loss from that to snow

                # [ a lot of `agg` is also just threshold driven convergence from things like sed, so also don't count sed against this let it go straight to snow if needed and can evap there ]

                # not sure if we should only use positive sedimentation convergence... i think so bc losing ice doesnt make things smaller
                # qi_tendency_sed =  max(FT(0), qi_tendency_sed)
                qi_tendency_sed =  resolve_nan(max(FT(0), qi_tendency_sed), FT(0)) # Only accept deposition, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))
                # qi_tendency_sed = FT(0)

                # qi_tendency_sub_dep = max(FT(0), qi_tendency_sub_dep) # this is the part that is new ice formation, so we can use that to limit the snow formation
                qi_tendency_sub_dep = resolve_nan(max(FT(0), qi_tendency_sub_dep), FT(0)) # Only accept deposition, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))

                # techincally acnv really can deplete existing ice too.. so add taht to limiter. include sed just bc

                # this version technically doesn't respect whatever the limiter did in creating qi_tendency_sub_dep, it's a static calculation from the timestep beginning.
                # S_qt_snow_ice_dep = limit_tendency(precipitation_tendency_limiter, -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor * my_conv_q_ice_to_q_sno(param_set, q, T, p; N = Ni), q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_dep = -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor * qi_tendency_sub_dep # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_dep_is = limit_tendency(precipitation_tendency_limiter, -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor * my_conv_q_ice_to_q_sno(param_set, q, T, p; N = Ni, thresh_only=true), q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]

                # this version respects whatever the limiter did in creating qi_tendency_sub_dep, even for the MM2015 variants
                S_qt_snow_ice_dep_is = -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor * my_conv_q_ice_to_q_sno_at_r_is_given_qi_tendency_sub_dep(param_set, qi_tendency_sub_dep, q.ice, Ni, ρ) # This is better i think  bc it respects the real sub_dep tendency
                S_qt_snow_ice_dep_above = -α_acnv * snow_formation_model.ice_dep_acnv_scaling_factor_above * my_conv_q_ice_to_q_sno_by_fraction(param_set, qi_tendency_sub_dep, q.ice, Ni, ρ) # I feel like this shouldn't really be rescalable idk...
                S_qt_snow_ice_dep = limit_tendency(precipitation_tendency_limiter, S_qt_snow_ice_dep_is + S_qt_snow_ice_dep_above, q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # technically this is the way to do it that respects the sub_dep from the limiter, even MM2015

                # This check may not be true (i.e. in principle you can really acnv more than the new ice formation in thier (not morrison 05) framework -- this check also encourages boosting acnv scaling w/o penalty which may lead to unrealistic outcomes from balanced parameters.
                # S_qt_snow_ice_dep = max(-max(0, qi_tendency_sub_dep), S_qt_snow_ice_dep) # Now that they're calculated differently, ensure S_qt_snow_ice_dep takes lss than deposition from qi_tendency_noneq
                


                # the limiter should probably pick one to give priority, dep i would assume

                # my_conv_q_ice_to_q_sno_thresh() is really just dispatching to my_conv_q_ice_to_q_sno_no_supersat attm, but multiplying by (r/r_ice_snow)^3


                #=
                get minimum q threshold for snow formation
                    path of calls is my_conv_q_ice_to_q_sno_thresh -> my_conv_q_ice_to_q_sno_no_supersat -> get_q_threshold_acnv(param_set, CMT.IceType(), q; N=N, assume_N_is = assume_N_is, r_acnv_scaling_factor=snow_formation_model.r_ice_acnv_scaling_factor)
                    this is kinda fragile though... the target threshold is kinda buried in there...
                    and if it's a q based threshold, does it even matter? maybe we just fall back to the default min threshold regardless? either way it would be nice to get out the threshold the fcn is using...
                =#
                    
                #=
                    This is a challenging implementation 
                    We would like, in steady state, for acnv loss to be supported by incoming vapor. Thus, we set the limit to allow for contributions from qi_tendency_sub_dep and qi_tendency_sed.
                    However, this is sub-ideal. When we pass qₜₕ based on current q, we gain sudden access to all the other tendencies that we didn't have before, which can lead to large jumps in the tendency.

                    So either we provide none, or always provide all. Right now we are opting for always providing all tendencies but this has the downside that acnv is no longer based on current q, but rather on the projected q.
                    Perhaps that is better? But you run the risk of e.g. having wildly misleading q relative to the acnv rates and production rates. e.g. if acnv and dep valance, q could be ~0., even if in steady state q should be much higher.

                    We also take the positive part only, so say dep/sed are nearly in balance, you gradually pass the threshold, acnv comes in and balances all of dep, then you add sed and crater suddenly. Not using a forward projection helps fix this.
                    Furthermore we should assume we're nearly in steady state, in which case the forward projection is just plain wrong. at very long timesteps it also becomes problematic.

                    A potential middle ground, since dep acnv is handled separately, would be to restrict the contribution to agg-acnv to the deposition part.


                =#

                # S_qt_snow_ice_agg, q_ice_thresh = get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set, q.ice; N = Ni, r_acnv_scaling_factor=snow_formation_model.r_ice_acnv_scaling_factor, ice_acnv_power = snow_formation_model.ice_acnv_power, w = w, dwdz = dwdz, use_sed=true) #  The positive source to snow is what this returns, so we do need to flip the sign for S_qs in the limit_tendency() call just as before
                # S_qt_snow_ice_agg, q_ice_thresh = get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set, q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt; N = Ni, r_acnv_scaling_factor=snow_formation_model.r_ice_acnv_scaling_factor, ice_acnv_power = snow_formation_model.ice_acnv_power) # For balance, allow incoming to contribute to acnv. Deposition acnv is already accounted for, so maybe don't double count? sedimentation should lead to acnv though.
                # S_qt_snow_ice_agg = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_agg , max(q.ice + (qi_tendency_sub_dep + S_qt_snow_ice_dep)*Δt + qi_tendency_sed*Δt - q_ice_thresh, 0), Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                # S_qt_snow_ice_agg = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_agg , max(q.ice - q_ice_thresh, 0) + qi_tendency_sed * Δt, Δt)  # I think allowing sed to join the party is a mistake bc we do not calculate acnv including that so you still get jumps. i.e. you hit the tresh and then all of a sudden sed ca take you wherever you want, you'd need to add it when calculating the thresh, as in the input being q.ice + sed tendency etc... but that get's more and more complicated, you could make the argument with more and more tendencies... etc... this is simple and in line wit everything else...
                
                # Nearly impossible to work with, almost always far too large in the wrong places...
                # S_qt_snow_ice_agg_mix = tke_eddy_diffusivity_driven_acnv(param_set, ice_type, qi, Ni, CMP.ρ_cloud_ice(microphys_params), ρ_c, tke_var; r_acnv_scaling_factor = snow_formation_model.r_ice_acnv_scaling_factor, C = cloud_sedimentation_model.E_ice_ice_mix, Dmax = cloud_sedimentation_model.ice_Dmax)
                # S_qt_snow_ice_agg_mix = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_agg_mix, max(q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, 0), Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                S_qt_snow_ice_agg_mix = FT(0)

                # [[ since we assume evap is constantly raising <r>, there is actually no need to enforce q_thresh... just proceed at τ or even faster -- you should never get <r> below your target value... honestly the longer time timestep the worse it should get lol. No threshold kind of means assuming <r> constant and so the rate is fixed.. maybe not a bad guess if you assume steady state idk... ]]
                use_cloak = (area_partition_model isa CoreCloakAreaPartitionModel)
                qi_tendency_sub_dep_for_thresh = max(qi_tendency_sub_dep + S_qt_snow_ice_dep_above + S_qt_snow_ice_dep_is, FT(0)) # remove the part already going to snow from sub dep
                S_qt_snow_ice_thresh, q_thresh_acnv = threshold_driven_acnv(param_set, ice_type, qi, Ni, ρ; Δt = Δt, Dmax=cloud_sedimentation_model.ice_Dmax, r_threshold_scaling_factor=snow_formation_model.r_ice_snow_threshold_scaling_factor, r_acnv_scaling_factor=snow_formation_model.r_ice_acnv_scaling_factor, add_dry_aerosol_mass=true, N_i_no_boost=N_i_no_boost, S_i=S_i, τ_sub_dep=τ_sub_dep, dqdt_sed=qi_tendency_sed, dqdt_dep=qi_tendency_sub_dep_for_thresh, dN_i_dz=dN_i_dz, dqidz=dqidz, w=term_vel_ice, N_INP=N_INP, massflux=massflux, domain=domain, use_cloak=use_cloak, tke=tke_var, qt_var=qt_var, qt = q.tot, dqtdz=dqtdz, qisedflux = term_vel_ice*qi) # pass in sed tendency so it can be accounted for in the threshold calculation
                S_qt_snow_ice_thresh = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_thresh, max(q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt - q_thresh_acnv, 0), Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]

                # TEST!
                # S_qt_snow_ice_thresh = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_thresh, max(q.ice + qi_tendency_sub_dep*Δt - q_thresh_acnv, 0), Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]

                # S_qt_snow_ice_thresh -= S_qt_snow_ice_dep # don't double dip [ i still don't know if this is right.... physics is iffy imo...]
                # S_qt_snow_ice_thresh = min(S_qt_snow_ice_thresh, FT(0)) # don't allow growth.

                # not sure if we should only use positive sedimentation convergence...
                # S_qt_snow = limit_tendency(precipitation_tendency_limiter, S_qt_snow_ice_dep + S_qt_snow_ice_agg_mix, q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth and aggregation # technically the part from sub_dep is new ice so add that in to the limiter (net bc the acnv tendency is neg of the sub tendency)
                

                # agg no longer just comes from us, so remove here...
                S_qt_snow = limit_tendency(precipitation_tendency_limiter, S_qt_snow_ice_dep + S_qt_snow_ice_agg_mix + S_qt_snow_ice_thresh, q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth and aggregation # technically the part from sub_dep is new ice so add that in to the limiter (net bc the acnv tendency is neg of the sub tendency)


                # S_qt_snow = -min(q.ice / Δt, α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q))
                # S_qt_snow = limit_tendency(precipitation_tendency_limiter, -α_acnv * my_conv_q_ice_to_q_sno_no_supersat(param_set, q; N0 = Ni, r_acnv_scaling_factor=snow_formation_model.r_ice_acnv_scaling_factor), q.ice, Δt)

            elseif snow_formation_model isa DefaultSnowFormationModel
                #=
                NOTE: NonEquilibriumSnowFormationModel was maybe a bad naming convention -- it's not 1:1 with NonEquilibriumMoistureModel... even EquilibriumMoistureModel can use this...
                So this is really just a fallback method...
                =#
                # S_qt_snow = -min(q.ice / Δt, α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice))

                qi_tendency_sed =  resolve_nan(max(FT(0), qi_tendency_sed), FT(0)) # Only accept deposition, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))
                qi_tendency_sub_dep = resolve_nan(max(FT(0), qi_tendency_sub_dep), FT(0)) # Only accept deposition, and if we passed in FT(NaN) as a fallback (e.g. in Equilibrium, just use FT(0))
                
                S_qt_snow = limit_tendency(precipitation_tendency_limiter, -α_acnv * CM1.conv_q_ice_to_q_sno_no_supersat(microphys_params, q.ice), max(q.ice + qi_tendency_sub_dep*Δt - CMP.q_ice_threshold(TCP.microphysics_params(param_set)), 0), Δt) # don't allow reducing below threshold

                if !isfinite(S_qt_snow)
                    @error("Got nonfinite S_qt_snow = $(S_qt_snow) from inputs qi = $(qi); Ni = $(Ni); ρ_c = $(ρ_c); p = $(p)")
                end

                S_qt_snow_ice_dep = FT(0) # no easy breakdown [ don't actually put NaN because it messes with the sums below and with the regridder afterwards], just put 0
                S_qt_snow_ice_dep_is = FT(0)
                S_qt_snow_ice_dep_above = FT(0)
                S_qt_snow_ice_agg_mix = FT(0) # no easy breakdown
                use_cloak = (area_partition_model isa CoreCloakAreaPartitionModel)
                qi_tendency_sub_dep_for_thresh = max(qi_tendency_sub_dep + S_qt_snow_ice_dep_above + S_qt_snow_ice_dep_is, FT(0)) # remove the part already going to snow from sub dep
                S_qt_snow_ice_thresh, q_thresh_acnv = threshold_driven_acnv(param_set, ice_type, qi, Ni, ρ_c; Δt=Δt, Dmax=cloud_sedimentation_model.ice_Dmax, r_threshold_scaling_factor=param_set.user_params.r_ice_snow_threshold_scaling_factor, r_acnv_scaling_factor=param_set.user_params.r_ice_acnv_scaling_factor, add_dry_aerosol_mass=true, N_i_no_boost=N_i_no_boost, S_i=S_i, τ_sub_dep=τ_sub_dep, dqdt_sed=qi_tendency_sed, dqdt_dep=qi_tendency_sub_dep_for_thresh, dN_i_dz=dN_i_dz, dqidz=dqidz, w=term_vel_ice, N_INP=N_INP, massflux=massflux, domain=domain, use_cloak=use_cloak, tke=tke_var, qt=q.tot, qt_var=qt_var, dqtdz=dqtdz, qisedflux=term_vel_ice*qi) # pass in sed tendency so it can be accounted for in the threshold calculation

                if !isfinite(S_qt_snow_ice_thresh) || ~isfinite(q_thresh_acnv)
                    @error("Got nonfinite S_qt_snow_ice_thresh = $(S_qt_snow_ice_thresh) or q_thresh_acnv = $(q_thresh_acnv) from inputs qi = $(qi); Ni = $(Ni); ρ_c = $(ρ_c)")
                end

                S_qt_snow_ice_thresh = limit_tendency(precipitation_tendency_limiter, -α_acnv * S_qt_snow_ice_thresh, max(q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt - q_thresh_acnv, 0), Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
                
                if !isfinite(S_qt_snow_ice_thresh) || ~isfinite(q_thresh_acnv)
                    @error("After limiting, got nonfinite S_qt_snow_ice_thresh = $(S_qt_snow_ice_thresh) or q_thresh_acnv = $(q_thresh_acnv) from inputs qi = $(qi); Ni = $(Ni); ρ_c = $(ρ_c)")
                end


                S_qt_snow = limit_tendency(precipitation_tendency_limiter, S_qt_snow + S_qt_snow_ice_thresh, q.ice + qi_tendency_sub_dep*Δt + qi_tendency_sed*Δt, Δt) # based on growth and aggregation # technically the part from sub_dep is new ice so add that in to the limiter (net bc the acnv tendency is neg of the sub tendency)
            else
                error("Unrecognized snow formation model $(snow_formation_model)")
            end

            
            qi_tendency_acnv = S_qt_snow # for storage
            qi_tendency_acnv_dep = S_qt_snow_ice_dep # for storage [ was already an input but may be removing ]
            qi_tendency_acnv_dep_is = S_qt_snow_ice_dep_is
            qi_tendency_acnv_dep_above = S_qt_snow_ice_dep_above
            qi_tendency_acnv_agg_mix = S_qt_snow_ice_agg_mix # for storage
            qi_tendency_acnv_thresh = S_qt_snow_ice_thresh

            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            ql_tendency += S_qt_rain
            qi_tendency += S_qt_snow

            L_v = TD.latent_heat_vapor(thermo_params, ts)
            L_s = TD.latent_heat_sublim(thermo_params, ts)

            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (L_v * S_qt_rain + L_s * S_qt_snow)

            # accretion cloud water + rain
            if rain_formation_model isa Clima1M_default
                # S_qr = min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ_c)) * precip_fraction
                S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ_c), q.liq, Δt) * precip_fraction
            elseif rain_formation_model isa Clima2M
                if rain_formation_model.type isa CMT.LD2004Type
                    # S_qr = min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ_c)) * precip_fraction
                    S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, liq_type, rain_type, q.liq, qr, ρ_c), q.liq, Δt) * precip_fraction
                elseif rain_formation_model.type isa CMT.KK2000Type || rain_formation_model.type isa CMT.B1994Type
                    # S_qr = min(q.liq / Δt, α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr, ρ_c), ) * precip_fraction
                    S_qr = -limit_tendency(precipitation_tendency_limiter, -α_accr * CM2.accretion(microphys_params, rain_formation_model.type, q.liq, qr, ρ_c), q.liq, Δt) * precip_fraction
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
            θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v
            ql_tendency_accr_liq_rai = -S_qr # for storage

            # accretion cloud ice + snow
            # S_qs = min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, snow_type, q.ice, qs, ρ_c)) * precip_fraction
            S_qs = -limit_tendency(precipitation_tendency_limiter,  -α_accr * CM1.accretion(microphys_params, ice_type, snow_type, q.ice, qs, ρ_c), q.ice, Δt) * precip_fraction

            qs_tendency += S_qs
            qt_tendency -= S_qs
            qi_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s
            qi_tendency_accr_ice_sno = -S_qs # for storage

            # sink of cloud water via accretion cloud water + snow
            # S_qt = -min(q.liq / Δt, α_accr * CM1.accretion(microphys_params, liq_type, snow_type, q.liq, qs, ρ_c)) * precip_fraction
            S_qt = limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, liq_type, snow_type, q.liq, qs, ρ_c), q.liq, Δt) * precip_fraction # microphysics params should contain E_liq_sno which can be a tuning knob here.
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
                θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1 + Rm / c_vm) * α - L_v)
                ql_tendency_accr_liq_sno = S_qt # for storage
            end


            # sink of cloud ice via accretion cloud ice - rain
            # S_qt = -min(q.ice / Δt, α_accr * CM1.accretion(microphys_params, ice_type, rain_type, q.ice, qr, ρ_c)) * precip_fraction
            S_qt = limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion(microphys_params, ice_type, rain_type, q.ice, qr, ρ_c), q.ice, Δt) * precip_fraction
            # sink of rain via accretion cloud ice - rain
            # S_qr = -min(qr / Δt, α_accr * CM1.accretion_rain_sink(microphys_params, q.ice, qr, ρ_c)) * precip_fraction
            S_qr = limit_tendency(precipitation_tendency_limiter, -α_accr * CM1.accretion_rain_sink(microphys_params, q.ice, qr, ρ_c), qr, Δt) * precip_fraction
            qt_tendency += S_qt
            qi_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1 / Π_m / c_pm * (S_qr * Lf * (1 + Rm / c_vm) + S_qt * L_s)
            qi_tendency_accr_ice_rai = S_qt # for storage

            # accretion rain - snow
            if T < T_fr
                # S_qs = min(qr / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, snow_type, rain_type, qs, qr, ρ_c)) * precip_fraction * precip_fraction
                S_qs = -limit_tendency(precipitation_tendency_limiter, -α_accr * my_accretion_snow_rain(microphys_params, snow_type, rain_type, qs, qr, ρ_c, term_vel_snow, term_vel_rain), qr, Δt) * precip_fraction * precip_fraction
            else
                # S_qs = -min(qs / Δt, α_accr * CM1.accretion_snow_rain(microphys_params, rain_type, snow_type, qr, qs, ρ_c)) * precip_fraction * precip_fraction
                S_qs = limit_tendency(precipitation_tendency_limiter, -α_accr * my_accretion_snow_rain(microphys_params, rain_type, snow_type, qr, qs, ρ_c, term_vel_rain, term_vel_snow), qs, Δt) * precip_fraction * precip_fraction
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
            qs_tendency_accr_rai_sno = S_qs # for storage [we calculate this here for the temperature dependence of the outcome]

            # liq - ice accretion [[ should this still be happening in equilibrium thermo? i think so but idk... ]]

            if moisture_model isa NonEquilibriumMoisture
                if T < T_fr # losing liq
                    S_qi = accretion_liq_ice(param_set, ql, qi, ρ_c, term_vel_ice, cloud_sedimentation_model.ice_terminal_velocity_scheme, Ni, cloud_sedimentation_model.ice_Dmax, cloud_sedimentation_model.E_liq_ice) #; liq_ice_accr_scaling_factor = cloud_sedimentation_model.liq_ice_collision_scaling_factor ) # should be a positive number
                    # Tendencies get multiplied by area outside of this fcn and cloud fraction is 1 if has condensate and 0 if not so no need here
                    S_qi = max(S_qi, FT(0)) # don't allow negative values (depending on ai, bi, ci in Chen that can just barely happen)
                    S_qi = -limit_tendency(precipitation_tendency_limiter, -α_accr * S_qi , ql, Δt) # signs should match accretion_snow_rain except snow ⟹ ice, rain ⟹ liq
                else # losing ice
                    # We don't have accretion liq-ice because we only send the result to liq or ice so it doesn't matter -- if it ever gets sent to precip we can bring this back...
                    # S_qi = accretion_liq_ice(microphys_params, ql, qi, ρ_c, term_vel_ice, cloud_sedimentation_model.ice_terminal_velocity_scheme, Ni, cloud_sedimentation_model.ice_Dmax, cloud_sedimentation_model.E_liq_ice; liq_ice_accr_scaling_factor = cloud_sedimentation_model.liq_ice_collision_scaling_factor ) # should be a positive number
                    # S_qi = max(S_qi, FT(0)) # don't allow negative values (depending on ai, bi, ci  Chen that can just barely happen)
                    # S_qi = limit_tendency(precipitation_tendency_limiter, -α_accr * S_qi , qi, Δt) # signs should match accretion_snow_rain except snow ⟹ ice, rain ⟹ liq
                    S_qi = FT(0) # PSACWI doesn't cover this, either way we have instant melting so...
                end
                qi_tendency += S_qi
                ql_tendency -= S_qi   
                # Because both ql and qi are part of the working fluid and θ_liq_ice_tendency is 0, we really should have  no change in θ_liq_ice_tendency
                # θ_liq_ice_tendency += FT(0)
                # θ_liq_ice_tendency += S_qi * Lf / Π_m / c_vm # fusion liq ⟹ ice
                qi_tendency_accr_ice_liq = S_qi # for storage
                ql_tendency_accr_liq_ice = -S_qi # for storage
            else
                # do nothing, bc the tendency isn't real bc it doesn't change qt so it doesn't matter...
            end
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
        qi_tendency_acnv_dep_is, # from crossing ice_snow boundary.
        qi_tendency_acnv_dep_above, # from already above threshold
        qi_tendency_acnv_agg_mix, # this is now from (us->other) and (other) so really, it shouldn't come at all from us! it's to other, on other's area!
        qi_tendency_acnv_thresh, # this is now from (us->other) and (other) so really, it shouldn't come at all from us! it's to other, on other's area!
        #
        ql_tendency_accr_liq_rai,
        ql_tendency_accr_liq_ice,
        ql_tendency_accr_liq_sno,
        #
        qi_tendency_accr_ice_liq,
        qi_tendency_accr_ice_rai,
        qi_tendency_accr_ice_sno,
        #
        qs_tendency_accr_rai_sno,
    )
end

 precipitation_formation(
    param_set::APS,
    moisture_model::AbstractMoistureModel,
    precip_model::AbstractPrecipitationModel,
    cloud_sedimentation_model::AbstractCloudSedimentationModel,
    rain_formation_model::AbstractRainFormationModel,
    snow_formation_model::AbstractSnowFormationModel,
    area_partition_model::AbstractAreaPartitionModel,
    qr::FT,
    qs::FT,
    ql::FT,
    qi::FT,
    Ni::FT,
    term_vel_ice::FT, # maybe this isn't needed lol
    term_vel_rain::FT,
    term_vel_snow::FT,
    area::FT,
    ρ_c::FT,
    p::FT,
    T::FT,
    Δt::FT,
    ts::TD.ThermodynamicState{FT},
    precip_fraction::FT,
    precipitation_tendency_limiter::AbstractTendencyLimiter,
    tke_var::FT,
    ;
    N_i_no_boost::FT = Ni, # for threshold acnv calculation
    qi_tendency_sed::FT = FT(NaN),
    S_i::FT = FT(NaN),
    τ_sub_dep::FT = FT(NaN),
    dN_i_dz::FT = FT(0),
    dqidz::FT = FT(0),
    N_INP::FT = FT(NaN),
    massflux::FT = FT(0),
    domain::AbstractDomain = Env,
    qt_var::FT = FT(0),
    dqtdz::FT = FT(0),
) where {FT} = 
    precipitation_formation(
        param_set,
        moisture_model,
        precip_model,
        cloud_sedimentation_model,
        rain_formation_model,
        snow_formation_model,
        area_partition_model,
        qr,
        qs,
        ql,
        qi,
        Ni,
        term_vel_ice,
        term_vel_rain,
        term_vel_snow,
        area,
        ρ_c,
        p,
        T,
        Δt,
        ts,
        # FT(NaN), # placeholder for backwards compat w [removed]
        # FT(NaN), # placeholder for backwards compat z [removed]
        FT(NaN), # placeholder for backwards compat ql_tendency_cond_evap
        FT(NaN), # placeholder for backwards compat qi_tendency_sub_dep
        qi_tendency_sed, # placeholder for backwards compat qi_tendency_sed
        precip_fraction,
        precipitation_tendency_limiter,
        tke_var,
        N_i_no_boost,
        S_i,
        τ_sub_dep,
        dN_i_dz,
        dqidz,
        N_INP,
        massflux,
        domain,
        qt_var,
        dqtdz,
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
    # prs::ACMP,
    param_set::APS,
    ql::FT, # stationary [[ note this actually does move now... ]]
    qi::FT, # moving
    ρ::FT,
    term_vel_ice::FT,
    velo_scheme::Union{CMT.Chen2022Type, CMT.Blk1MVelType},
    N_i::FT,
    ice_Dmax::FT,
    E_liq_ice::FT;
    # liq_ice_accr_scaling_factor::FT = FT(1),
    # D_transition::FT = FT(0.625e-3), # .625mm is the transition from small to large ice crystals in Chen paper, leave unchanged
) where {FT <: Real}

    accr_rate = FT(0)
    if (ql > FT(0) && qi > FT(0))

        # I think the chen velocity goes off the rails here, Chen doesnt stop it well enough.. so maybe we truncate? idk.
        # The plots at https://clima.github.io/CloudMicrophysics.jl/dev/TerminalVelocity/ assume N = 500e6 which is wildly unrealistic.
        # At <r> = 100 μm, it does suggest single droplet velocities around 0.2 m/s, but that grows quickly, so if N is very small it can get out of hand.


        integral_nav = int_nav_dr(param_set, ice_type, velo_scheme, qi, ρ, FT(NaN); Nt=N_i, Dmax=ice_Dmax) # assume μ in n(r) is 0 (exponential distribution)
        
        accr_rate = E_liq_ice * ql * integral_nav  # Eq 16 https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/#Accretion, ql * ρ goes from /kgair to /m^3 like nav, but to go back then we need to divide by ρ so we just ignore ρ completely


    end
    return accr_rate 
end

"""
A version w/ both moving. Note some of these would ned to use rain parameters, idk if they're all defined for liquid, but i havent made that update yet.
"""
function my_collision_coalescence(
    param_set::APS,
    type_i::CMTWaterTypes,
    type_j::CMTWaterTypes,
    qi::FT, # stationary [[ note this actually does move now... ]]
    qj::FT, # moving
    ρ::FT,
    term_vel_i::FT,
    term_vel_j::FT,
    # velo_scheme::Union{CMT.Chen2022Type, CMT.Chen2022Type},
    N_i::FT,
    N_j::FT,
    Dmax::FT,
    E::FT,
    ;
    # D_transition::FT = FT(0.625e-3), # .625mm is the transition from small to large ice crystals in Chen paper, leave unchanged
    liq_type::CMT.AbstractCloudType = liq_type,
    ice_type::CMT.AbstractCloudType = ice_type,
    rain_type::CMT.AbstractPrecipType = rain_type,
    snow_type::CMT.AbstractPrecipType = snow_type,
    μ_i::FT = FT(NaN),
    μ_j::FT = FT(NaN),
    add_dry_aerosol_mass::Bool = false,
) where {FT <: Real}

    prs = TCP.microphysics_params(param_set)



    Dmin = FT(0)

    accr_rate = FT(0)
    if (qi > FT(0) && qj > FT(0))

        if add_dry_aerosol_mass
            qi += mass(param_set, type_i, param_set.user_params.particle_min_radius, N_i; monodisperse=true) / ρ # add the dry aerosol mass to q when calculating λ and n0
            qj += mass(param_set, type_j, param_set.user_params.particle_min_radius, N_j; monodisperse=true) / ρ # add the dry aerosol mass to q when calculating λ and n0
        end

        # _n0_i::FT = CM1.n0(prs, qi, ρ, type_i)
        # _n0_j::FT = CM1.n0(prs, qj, ρ, type_j)
        _n0_i::FT = isnan(N_i) ? n0(prs, qi, ρ, type_i) : n0(param_set, qi, ρ, type_i, N_i; μ=μ_i) # use wanted N_i If given
        _n0_j::FT = isnan(N_j) ? n0(prs, qj, ρ, type_j) : n0(param_set, qj, ρ, type_j, N_j; μ=μ_j) # use wanted N_j If given

        _r0_j::FT = r0(prs, type_j)
        _m0_j::FT = m0(prs, type_j)
        _me_j::FT = me(prs, type_j)
        _Δm_j::FT = Δm(prs, type_j)
        _χm_j::FT = χm(prs, type_j)

        # _λ_i::FT = CM1.lambda(prs, type_i, qi, ρ)
        # _λ_j::FT = CM1.lambda(prs, type_j, qj, ρ)
        _λ_i::FT = lambda(param_set, type_i, qi, ρ, N_i, Dmin, Dmax; μ=μ_i) # use our lambda from terminal_velocity.jl
        _λ_j::FT = lambda(param_set, type_j, qj, ρ, N_j, Dmin, Dmax; μ=μ_j) # use our lambda from terminal_velocity.jl



        accr_rate =
            FT(π) / ρ *
            _n0_i *
            _n0_j *
            _m0_j *
            _χm_j *
            E *
            abs(term_vel_i - term_vel_j) / _r0_j^(_me_j + _Δm_j) * (
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

        _λ_i::FT = CM1.lambda(prs, type_i, q_i, ρ) # we have no N for precip so these are fine...
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


"""
This calculates self-collection rates, via collisions across domains.
For the condensate species (liquid and ice) crossing from env to updraft and vice versa was ordinarilly treated by just merging the cross flux in gently and without consequence.
In reality though, this should be a source of collision and self collection, which can lead to autoconversion..

In principle this should also lead to cross-species and other autoconversions, but we'll make the argument here that this is most critically relevant for ice-ice collisions creating snow, due to the large ice sizes in sedimentation

The other terms defintely do exist though and could be added to liq acnv and liq-ice accretion

The key point is that even if we assumed homogeneous development, this cross env flow will induce collisions


    #= note, the original flux `to_other` already could define the area swept out, the relation between that and this shouldn't bee too extravagant i'd imagine.
    The complex part is that there's a sedimentation rate and a rate of flux over.
    But note the flux over is similar to just taking the sedimentation flux of itself or another species., in reality it is just some fraction of that.

    So really, S_other / term_vel_other ≈ q_other, and we should just work in that framework...
    =#

"""

function cross_domain_self_collection_or_autoconversion(
    param_set::APS,
    type_s::CMTWaterTypes, # rn we assume precip is well mixed but we recently made it advection aware so perhaps this function could help precip too
    type_o::CMTWaterTypes,
    S_o::FT, # flux from self to other
    q_s::FT,
    q_o::FT, # other
    q_thresh_s::FT,
    q_thresh_o::FT,
    T::FT,
    ρ::FT,
    term_vel_s::FT, # should be term_vel + w really
    term_vel_o::FT, # should be term_vel + w really
    N_s::FT,
    N_o::FT,
    Dmin_s::FT,
    Dmax_s::FT,
    Dmin_o::FT,
    Dmax_o::FT,
    E::FT,
    acnv_power::FT,
    μ_s::FT = FT(NaN), # μ for self, if not given, use n0
    μ_o::FT = FT(NaN), # μ for other, if not given, use n0  
) where {FT <: Real}

    accr_rate = FT(0)

    prs = TCP.microphysics_params(param_set)

    if iszero(S_o) || iszero(q_o) || (iszero(term_vel_s) && iszero(term_vel_o))
        return accr_rate # no flux or no terminal velocity, no self collection
    end


    # is this really necessary? in principle, it's still just coming over as q, just with a different velocity no? I guess the area weighting is the hard part? idk i feel like it should still just be 1...
    q_s = S_o / abs(term_vel_s) # this is the flux over, so we can use it to define the self collection rate. This absolves us of the burden of using da/dz.

    # I think this would be close to q_s * a/(1-a) * w_s but no way to be sure, but the slope of da/dt mediates how much comes over, so it's not really that.

    # Literally only E's for collision are defined for liquid so use rain types instead



    # _n0_i::FT = CM1.n0(prs, q_s, ρ, type_s)
    # _n0_o::FT = CM1.n0(prs, q_o, ρ, type_o)
    _n0_s::FT = isnan(N_s) ? n0(prs, q_s, ρ, type_s) : n0(param_set, q_s, ρ, type_s, N_s; μ=μ_s) # use wanted N_s If given
    _n0_o::FT = isnan(N_o) ? n0(prs, q_o, ρ, type_o) : n0(param_set, q_o, ρ, type_o, N_o; μ=μ_o) # use wanted N_o If given
    _r0_o::FT = r0(prs, type_o)
    _m0_o::FT = m0(prs, type_o)
    _me_o::FT = me(prs, type_o)
    _Δm_o::FT = Δm(prs, type_o)
    _χm_o::FT = χm(prs, type_o)


    # _λ_s::FT = CM1.lambda(prs, type_s, q_s, ρ)
    # _λ_o::FT = CM1.lambda(prs, type_o, q_o, ρ)
    _λ_s::FT = lambda(param_set, type_s, q_s, ρ, N_s, Dmin_s, Dmax_s; μ=μ_s) # use our lambda from terminal_velocity.jl
    _λ_o::FT = lambda(param_set, type_o, q_o, ρ, N_o, Dmin_o, Dmax_o; μ=μ_o) # use our lambda from terminal_velocity.jl

    # unlike normal collisions, these things don't really coexist so we'll only have collisions `self` terminal velocity is faster than `other`
    # dv = if term_vel_s > 0 # we are falling, so the other one should be falling slower or rising
    #     (term_vel_o > -term_vel_s) ? abs(term_vel_s - term_vel_o) : FT(0) # going down (term_vel_s > 0), other has to be greater than -term_vel_s
    # else # we are rising, so the other one should be rising slower or falling
    #     (term_vel_o > term_vel_s) ? abs(term_vel_s - term_vel_o) : FT(0) # going up (term_vel_s < 0), other has to be greater than term_vel_s
    # end

    # i think the above is wrong... just because you're now crossing into the other domain doesn't mean the other domain doesn't have particles that can hit you from above.
    dv = abs(term_vel_s - term_vel_o) # this is the relative velocity, so we can use it to define the self collection rate

    accr_rate =
        FT(π) / ρ *
        _n0_s *
        _n0_o *
        _m0_o *
        _χm_o *
        E *
        dv / _r0_o^(_me_o + _Δm_o) * (
            FT(2) * CM1.SF.gamma(_me_o + _Δm_o + FT(1)) / _λ_s^FT(3) /
            _λ_o^(_me_o + _Δm_o + FT(1)) +
            FT(2) * CM1.SF.gamma(_me_o + _Δm_o + FT(2)) / _λ_s^FT(2) /
            _λ_o^(_me_o + _Δm_o + FT(2)) +
            CM1.SF.gamma(_me_o + _Δm_o + FT(3)) / _λ_s /
            _λ_o^(_me_o + _Δm_o + FT(3))
        )

    # return accr_rate # we shouldnt need resolve nan bc we're using the cm_1 values which shouldn't blow up? but just to be safe
    accr_rate = resolve_nan(accr_rate, FT(0))

    if !isequal(type_s, type_o) || iszero(accr_rate) # different types so not acnv, one will just take the result, so just return the accr_rate
        return accr_rate
    else # we only care if we acnv, self collection does nothing given we don't have prognostic N

        # estimate acnv probability based on how close each is to its own threshold (a measure of existing r...)
        
        # this really ought to asymptote to 1 or something, no?
        # also  for small acnv powers, this makes them bigger, not smaller...
        # if ((q_s / q_thresh_s) = .001, and acnv_power = 1/3, we get (.001 * .001)^(1/3 / 6) = 0.46, which is perhaps less than ideal...
        # also if 1 is close, the other shouldnt necessarily stop it, even if the rate is smaller, so maybe it should sum the qs? idk.
        # prob_acnv = ((q_s / q_thresh_s) * (q_o / q_thresh_o))^(acnv_power) # this is the probability of autoconversion, so we can use it to scale the accr_rate
        prob_acnv = ((q_s / q_thresh_s) + (q_o / q_thresh_o))^(acnv_power) # this is the probability of autoconversion, so we can use it to scale the accr_rate

        # acnv_power assumed to be in q/q_thresh units

        # prob_acnv = min(prob_acnv, FT(1)) # don't allow it to be greater than 1 (idk it's prolly fine... idk)

        # now we can scale the accr_rate by the probability of acnv

        return accr_rate * prob_acnv
    end
end


"""
Environment only acnv driven by collisions, driven by sgs tke or eddy diffusivity (not sure which yet)...

# For now we'll go w/ raw TKE
"""
function tke_eddy_diffusivity_driven_acnv(
    param_set::APS,
    q_type::CMTWaterTypes, # type of q, should be liquid or ice
    q::FT,
    N::FT,
    ρ_q::FT,
    ρ_a::FT, # air density
    # w::FT, # idk if we need this anymore... tke should do most of the talking... feels like terminal velocity should matter though idk.
    # eddy_diffusivity::FT, # eddy diffusivity [ i think tke is more approprate here]
    # tke::FT, # TKE
    tke_var::FT; # TKE or TKE dissipation rate [saffman-turner uses dissipation, but bc of our sed, i think raw TKE is more appropriate]
    r_acnv_scaling_factor::FT = FT(1.0), # scaling factor for the r threshold, used to adjust the threshold for acnv
    # acnv_power::FT = FT(1.0), # power for the acnv probability, used to adjust the threshold for acnv
    C::FT = FT(1.3), # empirical, we already have scaling parameters so maybe don't need this? or it should be calibrateable -- we need something to go from tke units to something more meaningful
    Dmax = FT(Inf),
    rain_type::CMT.RainType = rain_type,
) where {FT}


    microphys_params = TCP.microphysics_params(param_set)

    acnv_rate = FT(0)

    if q > FT(0) && (tke_var > FT(0))
        # w_tke = sqrt(tke)

        if isnan(N)
            if q_type isa CMT.IceType
                q_is_here::FT = q_is(param_set, q_type)
                q_threshold = CMP.q_ice_threshold(microphys_params)
                N = q_threshold / q_is_here # the implied number concentration at threshold for the largest possible droplets... so N and threshold are linked... ( i guess that makes sense? idk...)
            elseif q_type isa CMT.LiquidType
                q_lr_here = q_lr(param_set, q_type)
                q_threshold = CMP.q_liq_threshold(microphys_params)
                N = q_threshold / q_lr_here # the implied number concentration at threshold for the largest possible droplets... so N and threshold are linked... ( i guess that makes sense? idk...)
            else
                error("tke acnv not supported for q_type $q_type")
            end
        end

        # proportional to area which goes as r^2 * w, w here comes from (sqrt TKE)
        # r = r_from_qN(param_set, q_type, q, N; monodisperse = true)
        r = r_from_qN(param_set, q_type, q, N; monodisperse = false, ρ=ρ_a, Dmax=Dmax) # ρ=ρ_a, Dmax=Dmax)  # rn monodisperse is fine i think...? hopefully it's even a little pessimistic bc i think it's smaller than the real one? 

        # Note Saffman-Turner should only work at small scales, but we're at massive scales
        ϵ = tke_var # saffman-turner uses TKE dissipation rate
        w_tke = sqrt(2tke_var/ρ_a) # assuming tke is 1/2 ρ_a * w^2 per cubic meter, this solves the vertical velocity.
        # K_st = C * w_tke * (2r)^2 # saffman would be (2r^3) but we only care about area, mostly and vertical coherent motions, not isotropic volume like they care about for turbulence
        # coll_rate = FT(0.5) * K_st * N^2

        # Now we need to assess the probability this flux goes to actually acnv...
        # r_thresh = CMP.r_ice_snow(microphys_params) * r_acnv_scaling_factor
        # r_thresh = r_acnv(param_set, q_type, r_acnv_scaling_factor = r_acnv_scaling_factor) # this is the threshold for acnv, which is a function of the type and scaling factor
        # P_acnv = (r/r_thresh)^(3*acnv_power) # should we have a separate p here? [ if not, we keep units of q/q_thresh]

        #=
        Estimate the probability and mass fraction that two self-colliding particles 
        from a Marshall–Palmer distribution (n(r) ∝ exp(-λ * r)) coalesce into a 
        particle with radius exceeding a threshold `r_thresh`.

        Assumptions:
        - Collisions are volume additive: r_out = (r₁³ + r₂³)^(1/3)
        - Radii are drawn independently from an exponential (marshall-palmer) distribution with slope λ
        - The resulting r_out distribution is approximated as also being exponential with mean:
            r_mean_out = (2 * ⟨r³⟩)^(1/3) = (12 / λ⁴)^(1/3),  since ⟨r³⟩ = 6 / λ⁴
            i.e. two particles of average mass colliding (working in mass space is fairly weighted for this question)
            
        - No collision kernel or sticking efficiency

        Then, we can observe these results:
        1. Probability that r_out > r_thresh:
            P ≈ exp( - r_thresh / r_mean_out )

        2. Fraction of total mass contributed by such collisions:
            f_mass ≈ Γ(4, r_thresh / r_mean_out) / 6

        where Γ(4, x) is the upper incomplete gamma function:
            Γ(4, x) = ∫ₓ^∞ t³ e^{-t} dt

            This is because:
            The total mass of all coalesced particles is proportional to:
                ∫₀^∞ r³ · p(r) dr = ∫₀^∞ r³ · (1 / r_mean_out) · e^{-r / r_mean_out} dr
                = r_mean_out⁴ · Γ(4) = 6 · r_mean_out⁴

            The mass of coalesced particles with r > r_thresh is:
                ∫_{r_thresh}^∞ r³ · (1 / r_mean_out) · e^{-r / r_mean_out} dr
                = r_mean_out⁴ · Γ(4, r_thresh / r_mean_out)


            The total mass of all coalesced particles is proportional to:
            ∫₀^∞ r³ · p(r) dr
            = ∫₀^∞ r³ · (1 / r_mean_out) · e^{-r / r_mean_out} dr

            Let x = r / r_mean_out ⇒ dr = r_mean_out · dx, r³ = x³ · r_mean_out³:

            = r_mean_out⁴ · ∫₀^∞ x³ e^{-x} dx
            = r_mean_out⁴ · Γ(4)

            The mass of coalesced particles with r > r_thresh is:
                ∫_{r_thresh}^∞ r³ · p(r) dr
            = r_mean_out⁴ · ∫_{x = r_thresh / r_mean_out}^∞ x³ e^{-x} dx
            = r_mean_out⁴ · Γ(4, r_thresh / r_mean_out)

            Therefore:
                f_mass = Γ(4, r_thresh / r_mean_out) / Γ(4)
                    = Γ(4, r_thresh / r_mean_out) / 6

            So:
                f_mass ≈ Γ(4, r_thresh / r_mean_out) / 6
                where x = r_thresh / r_mean_out
            =#

        r_is = CMP.r_ice_snow(microphys_params)
        μ = μ_from_qN(param_set, q_type, q, N; ρ=ρ_a)
        # _n0::FT = isnan(N) ? n0(microphys_params, q, ρ_a, q_type) : n0(param_set, q, ρ_a, q_type, N; μ=μ, Dmax=Dmax) # use wanted Nt If given
        # λ = lambda(microphys_params, q_type, q, ρ_a, N, μ; Dmin=FT(0), Dmax=Dmax, _n0=_n0)
        # mean_r_out = cbrt(FT(12)/ λ^4) # this is the mean radius of the output particle, assuming self-collisions combine volumes and follow MP distribution
        # P_acnv = exp(-r_is / mean_r_out) # gpt derivation...
        F_acnv = acnv_mass_fraction(r, r_is; μ=μ)
        F_acnv = resolve_nan(F_acnv, FT(0)) # just to be safe, though it should be fine (the gammma can underfly at vry large values) Γ(4, 1e70) = 0 and Γ(4, 1e80) = NaN for rexample

        # now we need to go from coll_rate to acnv_rate
        # coll_rate *= P_acnv
        # acnv_rate = coll_rate * (q/N)
        # acnv_rate = coll_rate * q_from_r(param_set, r) * F_acnv # total mass in collisions, times the fraction that should autoconvert...

        #=
        Assume particle velocities are distributed as v ~ N(0, σ²),
        with mean zero but known mean absolute velocity: ⟨|v|⟩ = w.

        For a zero-mean normal distribution:
            ⟨|v|⟩ = sqrt(2/π) * σ  ⇒  σ = w * sqrt(π/2)

        The expected absolute velocity difference between two independent particles is:
            ⟨|v₁ - v₂|⟩ = sqrt(4/π) * σ = sqrt(2) * w

        This gives a simple estimate for relative motion in a normally
        distributed ensemble of particles with known |v|.
        =#

        # keep q so the distribution is the same, technically the collision rate will be about a factor of 4 too high, so /4
        acnv_rate = my_collision_coalescence(param_set, q_type, q_type, q, q, ρ_a, FT(sqrt(2)/2) * w_tke, -FT(sqrt(2)/2)*w_tke, N, N, Dmax, C; add_dry_aerosol_mass=true) * F_acnv # assume these are all tke-fueled (bc they're falling at the same speed) and one has -w_tke and the other w_tke speed
    end

    return resolve_nan(acnv_rate, FT(0))
end




"""
These use domain interactions so they need both sides to have been calculated already, and thus can't be folded into env and up specific functions.
"""
function compute_domain_interaction_microphysics_tendencies!(
    state::State,
    edmf::EDMFModel,
    Δt::Real,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
) 
    grid = Grid(state)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params = TCP.microphysics_params(param_set)

    tendencies_pr = center_tendencies_precipitation(state)
    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_en_f = face_aux_environment(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_bulk_f = face_aux_bulk(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    ρ_c = prog_gm.ρ
    # precip_fraction = compute_precip_fraction(edmf, state)
    N_up = n_updrafts(edmf)

    cloud_sedimentation_model = edmf.cloud_sedimentation_model
    # rain_formation_model = edmf.rain_formation_model
    snow_formation_model = edmf.snow_formation_model    
    precipitation_tendency_limiter = get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters)



    # acnv from environment to updraft and vice versa


    w_en::CC.Fields.Field = Ic.(aux_en_f.w)
    w_up::CC.Fields.Field = Ic.(aux_bulk_f.w)



    Dmin = FT(0)
    Dmax = cloud_sedimentation_model.ice_Dmax
    ice_acnv_power = snow_formation_model.ice_acnv_power
    r_ice_acnv_scaling_factor = snow_formation_model.r_ice_acnv_scaling_factor
    E_ii = cloud_sedimentation_model.E_ice_ice


    #=
        This should be a flux from other and us, to snow. We will keep the instantaneous part and assume that the rest of the ice mixes in and assumes the other terminal velocity quickly.
        Thus the limit must include the sedimentation tendency from `other`
    =#

    kc_toa = kc_top_of_atmos(grid)
    kc_surf = kc_surface(grid)

    area = aux_en.area # this is the area of the environment, which is what we use to define the fluxes
    _area_top = isa(area, FT) ? area : area[kc_toa]
    _area_bottom = isa(area, FT) ? area : area[kc_surf]
    C2Fa = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(_area_bottom)), top = CCO.SetValue(FT(_area_top)))
    wvec = CC.Geometry.WVector
    # ∇c = CCO.DivergenceF2C() # F2C to come back from C2F
    ∇a = @. ∇c(wvec(C2Fa(area))) # seems stable w/ either sign w...? (Left bias this?)



    @inbounds for k in real_center_indices(grid)
        kf = CCO.PlusHalf(k.i)
        dz = grid.zf.z[kf+1] - grid.zf.z[kf]  # this is how it'z defined



        if (aux_bulk.area[k] > FT(0)) && (aux_bulk.qi_tendency_sedimentation_other[k] > FT(0))
            
            
            # env -> bulk means env area increases with z, da/dz > 0 
            interface_area = ∇a[k] * (∇a[k] > 0) * dz[k] # the real area where these collisions can happen...


            # we do env -> bulk, the bulk part can be divided among each updraft (bc to_other is only calculated for bulk)
            qi_tendency_sedimentation_other_blk = aux_bulk.qi_tendency_sedimentation_other[k] / aux_bulk.area[k] # de-area weight
            S_qt_snow_ice_agg_other_env_blk = cross_domain_self_collection_or_autoconversion(microphys_params,
                ice_type, ice_type,
                qi_tendency_sedimentation_other_blk, # flux from env to bulk
                aux_en.q_ice[k], aux_bulk.q_ice[k],
                get_q_threshold_acnv(param_set, ice_type; N = aux_en.N_i[k], assume_N_is = false, r_acnv_scaling_factor = r_ice_acnv_scaling_factor),
                get_q_threshold_acnv(param_set, ice_type; N = aux_bulk.N_i[k], assume_N_is = false, r_acnv_scaling_factor = r_ice_acnv_scaling_factor),
                aux_en.T[k], ρ_c[k], # air density
                aux_en.term_vel_ice[k] - w_en[k], aux_bulk.term_vel_ice[k] - w_up[k],
                aux_en.N_i[k], aux_bulk.N_i[k],
                Dmin, Dmax, Dmin, Dmax,
                E_ii, ice_acnv_power,
                )
            S_qt_snow_ice_agg_other_env_blk = limit_tendency(precipitation_tendency_limiter, -S_qt_snow_ice_agg_other_env_blk, aux_bulk.q_ice[k] + qi_tendency_sedimentation_other_blk*Δt, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
        
            # Agg cross-domain acnv || Env --> Bulk
            qi_tendency_agg_acnv_other_env_blk = S_qt_snow_ice_agg_other_env_blk 
            @inbounds for i in 1:N_up # split the env->bulk over the updrafts
                qi_tendency_agg_acnv_other_env_blk_i = S_qt_snow_ice_agg_other_env_blk # already should be area weighted
                Π_m = TD.exner(thermo_params, aux_up[i].ts[k])
                c_pm = TD.cp_m(thermo_params, aux_up[i].ts[k])
                L_s = TD.latent_heat_sublim(thermo_params, aux_up[i].ts[k])
                θ_liq_ice_tendency_agg_acnv_other_env_blk = S_qt_snow_ice_agg_other_env_blk / Π_m / c_pm * L_s

                aux_up[i].qi_tendency_acnv_agg_other[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_up[i].area[k] # add back area weight
                aux_up[i].qi_tendency_acnv[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_up[i].area[k] # add back area weight
                aux_up[i].qt_tendency_precip_formation[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_up[i].area[k] # add back area weight
                aux_up[i].θ_liq_ice_tendency_precip_formation[k] += θ_liq_ice_tendency_agg_acnv_other_env_blk * interface_area / aux_up[i].area[k] # add back area weight
                aux_bulk.θ_liq_ice_tendency_precip_formation[k]  += θ_liq_ice_tendency_agg_acnv_other_env_blk * interface_area / aux_up[i].area[k]  # this one is T depdendent so add here
                if edmf.moisture_model isa NonEquilibriumMoisture
                    aux_up[i].qi_tendency_precip_formation[k] += qi_tendency_agg_acnv_other_env_blk * aux_up[i].area[k] # this is the same as qt_tendency_precip_formation
                end
            end
            aux_bulk.qi_tendency_acnv_agg_other[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_bulk.area[k] # this is the total tendency for the bulk
            aux_bulk.qi_tendency_acnv[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_bulk.area[k] # this is the total tendency for the bulk
            aux_bulk.qt_tendency_precip_formation[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_bulk.area[k] # this is the total tendency for the bulk
            if edmf.moisture_model isa NonEquilibriumMoisture
                aux_bulk.qi_tendency_precip_formation[k] += qi_tendency_agg_acnv_other_env_blk * interface_area / aux_bulk.area[k] # this is the same as qt_tendency_precip_formation
            end

            tendencies_pr.q_sno[k] -= qi_tendency_agg_acnv_other_env_blk * interface_area / aux_bulk.area[k] # this is the total tendency for the environment and bulk, so we need to subtract it from the total

        end

        ## -------------------------- ##

        # do the reverse, from bulk to env (really only one of these should be nonzero at a time I guess...)
        if (aux_en.area[k] > FT(0)) && (aux_en.qi_tendency_sedimentation_other[k] > FT(0))

            interface_area = -∇a[k] * (∇a[k] < 0) * dz[k] # the real area where these collisions can happen...


            qi_tendency_sedimentation_other_en = aux_en.qi_tendency_sedimentation_other[k] / aux_en.area[k] # de-area weight
            S_qt_snow_ice_agg_other_blk_env = cross_domain_self_collection_or_autoconversion(microphys_params,
                ice_type, ice_type,
                qi_tendency_sedimentation_other_en,
                aux_bulk.q_ice[k], aux_en.q_ice[k],
                get_q_threshold_acnv(param_set, ice_type; N = aux_bulk.N_i[k], assume_N_is = false, r_acnv_scaling_factor = r_ice_acnv_scaling_factor),
                get_q_threshold_acnv(param_set, ice_type; N = aux_en.N_i[k], assume_N_is = false, r_acnv_scaling_factor = r_ice_acnv_scaling_factor),
                aux_bulk.T[k], ρ_c[k], # air density
                aux_bulk.term_vel_ice[k] - w_up[k], aux_en.term_vel_ice[k] - w_en[k],
                aux_bulk.N_i[k], aux_en.N_i[k],
                Dmin, Dmax, Dmin, Dmax,
                E_ii, ice_acnv_power,
                )
            S_qt_snow_ice_agg_other_blk_env = limit_tendency(precipitation_tendency_limiter, -S_qt_snow_ice_agg_other_blk_env, aux_en.q_ice[k] + qi_tendency_sedimentation_other_en*Δt, Δt) # based on growth but threshold is fixed [ needs a scaling factor for how much is over the thresh]
        
        
            # Agg cross-domain acnv || Bulk --> Env
            qi_tendency_agg_acnv_other_blk_env = S_qt_snow_ice_agg_other_blk_env # already should be area weighted
            Π_m = TD.exner(thermo_params, aux_en.ts[k])
            c_pm = TD.cp_m(thermo_params, aux_en.ts[k])
            L_s = TD.latent_heat_sublim(thermo_params, aux_en.ts[k])
            θ_liq_ice_tendency_agg_acnv_other_blk_env = S_qt_snow_ice_agg_other_blk_env / Π_m / c_pm * L_s

            aux_en.qi_tendency_acnv_agg_other[k] += qi_tendency_agg_acnv_other_blk_env * interface_area/ aux_en.area[k] # this is the total tendency for the environment
            aux_en.qi_tendency_acnv[k] += qi_tendency_agg_acnv_other_blk_env * interface_area / aux_en.area[k]
            aux_en.qt_tendency_precip_formation[k] += qi_tendency_agg_acnv_other_blk_env * interface_area/ aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] += θ_liq_ice_tendency_agg_acnv_other_blk_env * interface_area / aux_en.area[k]
            if edmf.moisture_model isa NonEquilibriumMoisture
                aux_en.qi_tendency_precip_formation[k] += qi_tendency_agg_acnv_other_blk_env * interface_area / aux_en.area[k] # this is the same as qt_tendency_precip_formation
            end

            tendencies_pr.q_sno[k] -= qi_tendency_agg_acnv_other_blk_env * interface_area / aux_en.area[k] # this is the total tendency for the environment and bulk, so we need to subtract it from the total

        end


        # now we can add these to the tendencies





        # sum total
    end
    
    return nothing
end



"""
Because we diagnose n_0, μ, λ from q and N, in principle there's a radius assumption.
Inherently, some of the moisture exists above the threshold between ice and snow.

The greater the fraction that exists above the threshold, the less consistent we are being about dividing between ice and snow (which is arbitrary to begin with)

To mitigate, this function allows relaxation above the threshold. 

To avoid continually relaxing away your ice, rather than attempting to autoconvert all the ice above the threshold (when we know a good fraction must exist),
we instead try to keep <r> less than 1/4 the threshold.

When `N` is not provided, this function cannot be assumed to work correctly...

The default q, N relations for liquid and ice assumes a fixed n_0 (though snow does not)
This means that N = n0/λ

recall λ = (_χm * _m0 * _n0 * CM1.SF.gamma(_me + _Δm + μ + FT(1)) / ρ / q / _r0^(_me + _Δm))^FT(1 / (_me + _Δm + μ + 1))

This means in general λ decreases with q, and N increases, meaning the result will always scale w/ q.

Note we have <r> = 1/λ so we know <r> decreases as lambda increases or as q shrinks... Thus we cannot have larger radii at smaller q.
This is perhaps bad for ice down low, albeit being correct for ice up high....

Ideally then, when N is not provided, we just assume that the background profile is scaled such that no acnv happens...

When we do have N, we assume acnv ramps up past r_acnv which is as low as .2 r_is.

Notably thse also rn share autoconversion_timescale, so if you have N it is good if this turns off... perhaps... idk.


This is similar to Ferrier et al 1994, acnv. 


This may very well still grow w/ q...
Some strategies would be to:
    1. scale strongly by r/r_thresh, so if we have a lot of large particles, even if not over the threshold, we go faster than if we had smaller particles...
    2. Subtract out deposition under the guise that this is part of the same process, and that part was already taken care of, so we don't need to worry about it. This is just a backup process.



This should be slowed in superaturated environments since you are constantly losing the largest particles to dep acnv and sed and generating new ones so threshold acnv is less likely.
In subsaturated environments it can be come more likely, since we both lose the small particles to evap and only have large particles since they had to have sedimented in.
So τ should be mediated to be slower in supersat and faster in subsat.
At the same time, even the subsat rate should be mediated by ∂qi/∂z, since it's the size sorting that matters.

At the same time, we know we have PITOSN at supersat due to unresolved supersaturation variability in updrafts/downdrafts, so maybe massflux should mediate in suprsat?

# =================================================================================== #
Threshold autoconversion for a gamma PSD under the <r> fixed assumption.

The standard relaxation toward a critical radius rₜₕ is:

    dq/dt = -(q - q_th)/τ
    q_th = q * (rₜₕ / <r>)³

Enforcing <r> fixed implies that the number concentration N decreases proportionally to q:

    dN/dt = (N / q) * dq/dt   ⇒   <r> = constant

The effective relaxation timescale can be written as:

    τ_eff = τ / [1 - (rₜₕ / <r>)³] > τ

so the relaxation toward the threshold is slightly slower, since we would start to more gently asymptote. No hard threshold enforcement is required; q naturally relaxes toward q_th while <r> remains constant.

Semi-analytic timestep update:

    q_new = q_old * exp(-Δt / τ_eff)
    N_new = N_old * (q_new / q_old)

You can choose to use the semi-analytic or not. I think we should not, so that we can assume steady state since this will balance sedimentation.

When supersaturated, we do not move towards 0, but keep the move only towards q_th. (though if it's supersat fluctuation based who knows lol)

If we're going up (domain = Up() or CloakUp() or Bulk()), we'll forgo supersaturated pitosn...

Note while, when subsaturated, this is driven by evaporation raising <r>, there's no clear exchange rate. small losses in q to sublimation could mean larger increases in <r> that require removing large particles (lots of mass) to rectify.
So it's sort of hard to turn this on and off... you basically get up to <r> = r_th and then lose everything suddenly...
We will try to ease into it as we approach r = r_th if subsat to represent underlying variability, and changes in <r> over a timestep even if the initial <r> is below r_th.


# ================================================================================== #

When S_i < 0 , since evaporation keeps increasing <r>, we go all the way to q = 0 eventually. because of this strong bound, we can afford to not really turn on until we get close to r_th.
Still, we're not as into instantaneous acnv as Morrison, so we want to ramp up as we approach r_th.
By r = r_th, we should still have fairly quick growth since small changes in N from the low end evaporating can be canceled by huge changes in q from the large end of N autoconverting in response....
Should converge at S_i = 0 to 0 though (up to variance)

However, when S_i > 0, we only go to q_th since new small particles are being formed all the time... Technically this wouldn't be until the 2q_is limit.
For continuity, however, and because we don't believe in instantaneous acnv in general like Morrison, we'd like the target to converge to q_th as S_i -> 0+. We don't want to cross and go from a target of 2*r_is to r_thresh instantly. Arguably it would be safer to just use r_thresh for everything, and we might... idk.
The rate also converges to lower number then too though since sub_dep goes down.


Basically the LES believes in letting everything grow to r_max = 2r_is and then transferring everything at once, we take a more gradual stance where starting sooner, we start bringing things over such that we never actually get up to r_thresh.
In that way our `r_thresh` is a unique entity -- it both represents the radius at which mixing with dry air ensures acnv, but on the growth side is a contract for when to limit growth. if you ever see <r> = 2 r_is, we've gone too far since to match the LES we would need sudden catastrophic acnv at S_i > 0...
You could argue that this evap r_thresh, and growth r_thresh could be different, for simplicity we use the same number...

Maybe one could look at the fraction of mass above 2r_is and try to keep it down?


A complicating factor is in the environment where we have σ(q) that can drive thresh acnv even at S_i > 0... Our heuristic for this is :


Finally, since we do have this heuristic of having some acnv before getting to r_thresh, q_target just snapping to 0 at S_i = 0 is a little strong... since the decay rate slows over the timestep, we can get an effective rate...
then, q_target would be q at the end of the timestep = q * exp(-Δt / τ)

When we lose some small particles to evap and then autoconvert some large ones to balance it out, we in this framework keep the q in the middle... so we should start to slow down as <r> comes down as we lose more of the largest particles...
losing one 10 micron particle balanced by losing 1 100 micron particle is a huge mass loss yes, but there are so many more 1 micron particles that soon you'll run out of large particles and stop...
However, even if we restricted the loss to particles above 2*r_is, that is like 95% of the mass still...

"""
function threshold_driven_acnv(
    param_set::APS,
    q_type::Union{CMT.LiquidType, CMT.IceType}, # type of q, should be liquid or ice
    q::FT,
    N::FT,
    ρ_a::FT; # air density
    Δt::FT = FT(Inf),
    Dmax = FT(Inf), # maximum diameter for the particles, used to limit the
    μ::FT = FT(NaN),
    r_threshold_scaling_factor::FT = FT(1),
    r_acnv_scaling_factor::FT = FT(1), # scaling factor for the r threshold, used to adjust the threshold for acnv
    add_dry_aerosol_mass::Bool = false,
    N_i_no_boost::FT = FT(NaN),
    S_i::FT = FT(0),
    ∂qi∂z::FT = FT(0), # vertical gradient of qi, used to adjust the timescale for acnv
    drdz::FT = FT(0),
    τ_sub_dep::FT = FT(NaN),
    dqdt_sed::FT = FT(0),
    dqdt_dep::FT = FT(0),
    dN_i_dz::FT = FT(0),
    dqidz::FT = FT(0),
    w::FT = FT(0),
    N_INP::FT = FT(NaN),
    massflux::FT = FT(0),
    domain = Env,
    use_cloak::Bool = false,
    tke::FT = FT(0), # not sure which vars i want
    qt::FT = FT(NaN), # not sure which vars i want
    qt_var::FT = FT(0), # not sure which vars i want
    dqtdz::FT = FT(0), # not sure which vars i want
    qisedflux::FT = FT(0),
) where {FT <: Real}


    acnv_rate = FT(0)
    q_target = q # default to doing nothing

    if (q > 0) && !isnan(N)
        if iszero(N)
            N = FT(NaN) # have an invalid prediction
            error("Got positive q = $q but N = 0, cannot calculate threshold_driven_acnv")
        end

        microphys_params = TCP.microphysics_params(param_set)


        # do_print = (rand() < 1e-3) && (domain isa EnvDomain) # for debugging

        is_growing = S_i > FT(0) # need to use this because we use S_i for exponentiation and such
        # is_growing = iszero(dqdt_dep) ? (S_i > FT(0)) : (dqdt_dep > FT(0)) # if dep is zero, use S_i, else use dep rate. dep rate might be better bc better captures the entire timestep...
        # r_i_acnv = r_ice_acnv(param_set, r_acnv_scaling_factor)

        r_is = get_r_cond_precip(param_set, q_type)
        r_thresh = r_is * r_threshold_scaling_factor

        if is_growing
            r_max = 2*r_is # the MG max value... maybe we should ramp up from something small at r = r_is to max power at 2*r_is

            # we start acnv long before r_max...
            # some factor on r_thresh but yeah instead of letting this go all the way to 2 r_is and then crashing down, we just start siphoning off mass earlier so it's a continous process. Staying near to r_thresh for sublimation acnv is good cause it means things are more continous.
            # in highly supersaturated environments, maybe we tolerate slightly higher <r>... but not that much. at 20% supersat maybe we go halfway to r_max (cubic mean of r^3 and r_thresh^3?, e.g. blend Ns, assuming m_e = 3). at 0% supersat we go to r_thresh Only
            r_max = ((1 - min(S_i / 0.2, 1)) * r_thresh^3 + min(S_i / 0.2, 1) * ((r_thresh^3 + (2*r_is)^3) / 2))^(1/3)


            
            if ((domain isa UpDomain) || (domain isa BulkDomain)) #  || (domain isa CloakUpDomain)
                # let's have this converge to r_thresh at S_i = 0+, and to full r_max at S_i = 0.05 (5% supersat)
                r_target = r_thresh + (r_max - r_thresh) * clamp(S_i / FT(0.05), FT(0), FT(1))
            else # Growing, but in Env or Cloak. Because of assumed variance, we will stay close to r_thresh for longer, not reaching r_max until 20% supersaturation
                r_target = r_thresh + (r_max - r_thresh) * clamp(S_i / FT(0.2), FT(0), FT(1))
            end
        else # Shrinking. Here we want to reach r_thresh quickly since evaporation is pushing <r> up.
            r_target = r_thresh
        end

        if isnan(μ)
            μ = μ_from_qN(param_set, q_type, q, N; ρ=ρ_a)
        end
        if add_dry_aerosol_mass
            q += mass(param_set, q_type, param_set.user_params.particle_min_radius, N; monodisperse=true) / ρ_a # add the dry aerosol mass to q when calculating λ and n0
        end
        r = r_from_qN(param_set, q_type, q, N; monodisperse = false, ρ=ρ_a, Dmax=Dmax, μ=μ)


        # == Undo Boost == # [[ Prolly shouldn't do it for the cloak... ]] The idea is that the boost is small particles coming down in dryer air and maybe aren't the bulk of the threshold acnv, but .... debatable lol. you could also argue theyre the ones drying out more and autoconverting... maybe we should just deprecate this.
        if !isnan(N_i_no_boost) && !use_cloak
            # at 5 percent subsat reach r_thresh, at 5% supersat reach r_i_acnv?
            if !is_growing
                r_boost = r
                N_i_boost = max(N - N_i_no_boost, FT(0)) # if subsat, any imported N_i should be large...
                q_boost = (N_i_boost/N) * q
                q_no_boost = (N_i_no_boost/N) * q
                r_no_boost = r
            else
                # in supersaturated regimers, it's harder to know what r is. we'll assume it's smaller than the current r, catching up to r_i_acnv at 5% supersat
                
                if !use_cloak
                    # Basically the downwelling particles are dryer, so maybe don't grow as much (there's some supersaturation deficit)
                    # We observe PITOSN at almost any RH, I think any real MF should engender r_boost to hit r_thresh regardless of supersat
                    # Say that deficit is 5-10%. Then r_boost hits r at 15% supersat the maximum of that and 1.1 r_thresh at 5% supersat

                    # start at original r, reach 1.2r_th at MF = 0.02 [[[ this is too strong... you can go from <r> = 20 micron all the way to way above threshold instantly... even at super high supersaturation.... and it also puts most mass in the downcloak which probably isn't true.... ]]]

                    # We're growing so supersat, but we need to converge to meet the !isgrowing branch at r_boost = r_no_boost at S_i = 0. When more supersaturated we should probably assume even less in the downdraft? but it's not clear... We know $τ is slower than N, r suggest so maybe that's less evap because more is in the non cloak... So we'll go 0.9 r at 10% supersat, 1.0 r at 0% supersat
                    r_boost = r - (r - 0.9r) * clamp(S_i / FT(0.1), FT(0), FT(1))

                    # if subsat, any imported N_i should be large...
                    N_i_boost = max(N - N_i_no_boost, FT(0)) # might need to be some fraction of this idk... going from r of say 50 microns to r_boost of 70 microns is a 5x boost in q, might be too extreme.

                    
                    q_boost = (N_i_boost/N) * q
                    N_i_boost = N_from_qr(param_set, q_type, q_boost, r_boost; ρ=ρ_a, monodisperse=false) # change (mostly reduce) N_i_boost to match r_boost
                    q_no_boost = (N_i_no_boost/N) * q
                    r_no_boost = r
                    # If the above is too strong, maybe try just relaxing the part above r_thresh but keeping N_boost or something idk... or changing the timescale.
                else # Cloak the boost is inherent, so ignore
                    r_boost = r
                    N_i_boost = N
                    q_boost = q
                    r_no_boost = r
                    q_no_boost = q
                end
            end
        end
        

        if q_type isa CMT.IceType
            τ = FT(param_set.user_params.τ_acnv_sno_threshold) # use separate timescale bc this should be order timestep scale (in M2005), while real acnv should be like 10^4 seconds.
        elseif q_type isa CMT.LiquidType
            # τ = FT(param_set.user_params.τ_acnv_liq_threshold)
            error("threshold_driven_acnv not implemented for liquid")
        else
            error("threshold_driven_acnv not supported for q_type $q_type")
        end
        τ_orig = τ


        # TODO: Turn this into a smooth transition, allow more masflux effect at supersat as proxy for qt_var.
        acnv_rate_other = FT(0)
        acnv_rate_fluc = FT(0)
        if is_growing # growing, so we only go to q_th, not 0

            if ((domain isa UpDomain) || (domain isa BulkDomain))
                # we're at max r so autoconvert quickly, tapering off as we approach r_target since there we are essentially balancing depositional growth once you reach the target.. so use the same timescale as a rough heuristic, but ttight taper so small γ
                f = clamp((r - r_target) / (FT(0.05) * r_target), FT(0), FT(1))  # normalized ramp
                γ = FT(0.05)  # adjust to control how fast τ drops
                τ = τ_sub_dep^(1 - f^γ) * τ^(f^γ)

                # redirect deposition: by r = r_max we want to redirect all of it to acnv_thresh, at r=r_th we want none of it redirected
                acnv_rate_other += dqdt_dep * clamp((r - r_target) / (r_max - r_target), FT(0), FT(1)) # re-direct deposition -- by r = r_max, redirecting all of it to acnv_thresh
            else # Env or Cloak growing.
                #=
                    Here it's slightly more complex. We of course don't want to hit r_target. But we also have some variance like metrics for staying close to r_thresh since σ(q) is real in the environment
                    The growth part right now is tied to 2*r_is, but the variance part should be tied to r_thresh... To compromise, we used a lower r_target
                =#

                # Don't let growth get out of hand. At r = r_thresh, redirect none, at r = r_max redirect all.
                # acnv_rate_other += dqdt_dep * clamp((r_no_boost - r_thresh) / (0.2r_thresh), FT(0), FT(1)) # re-direct deposition -- by r = 1.2 r_th, redirecting all of it to acnv_thresh
                acnv_rate_other += dqdt_dep * clamp((r - r_target) / (r_max - r_target), FT(0), FT(1)) # re-direct deposition -- by r = r_max, redirecting all of it to acnv_thresh

                # # redirect sedimentation
                # acnv_rate_other += dqdt_sed * clamp((r - r_target) / (r_max - r_target), FT(0), FT(1)) # re-direct sedimentation -- by r = r_max, redirecting all of it to acnv_thresh

                # growth regime so probably no excursions that force acnv... can still tie to MF
                # [[ unlike subsat where the default acnv rate is fast, but is slowed by MF bringing down particles and moisture and spreading the size dist]], we're supersaturated here, so the default rate is 0, not fast.
                # so here, MF can speed things up by making qt variability. have supersaturation fluctuations etc. So this will depend on MF and S_i
                # τ = τ_sub_dep + (τ_orig - τ_sub_dep) * clamp((r - r_thresh) / (FT(0.05) * r_thresh), FT(0), FT(1)) # tau gets fast as we approach r_thresh, peaks at 1.05 r_thresh...

                # Adjust to asymptote to target
                f = clamp((r - r_target) / (FT(0.05) * r_target), FT(0), FT(1))  # normalized ramp
                γ = FT(0.05)  # adjust to control how fast τ drops
                τ = τ_sub_dep^(1 - f^γ) * τ^(f^γ)

                # Adjustment towards slower τ as S_i → 0+ not needed because growth acnv isn't sub/dep dependent and we have a defiend target that also shrinks as S_i → 0+
                
                
                # Mass flux (MF) generates fluctuations that *enable* autoconversion. [[ this probably matters more when not cloaking...]] This pushes us towards faster acnv from runaway losses on the dry side... this keeps τ closer to τ_orig, and also might push the target q down a bit...
                
               
                # Fluctuations [[ fluctuations also do have a target limit, since not all of the particles are susceptible to fluctuations...]]
                if use_cloak # smaller but still existent fluctuations in env
                    M_norm = clamp(massflux / FT(0.02), FT(0), FT(1)) # At MF=0 → τ=∞ (no acnv), At MF=0.02 → τ=τ (normal, fast), Clamp ensures saturation at MF=0.02
                    S_i_max = FT(0.05) * M_norm  # vanish at 5% supersat for MF=0.02, linearly less for smaller MF
                    S_norm = clamp(one(FT) - S_i / max(S_i_max, eps(FT)), FT(0), FT(1)) # at S_i = S_i_max, S_norm = 0, so τ → ∞, at S_i = 0, S_norm = 1, so τ = τ
                    τ_fluc = τ / max(M_norm * S_norm, eps(FT)) # as MF→0, τ→∞; as MF→0.02, τ→τ]
                    
                else # Entertain more MF based fluctuations when not cloaking. Here we update acnv_rate_fluc
                    M_norm = clamp(massflux / FT(0.02), FT(0), FT(1)) # At MF=0 → τ=∞ (no acnv), At MF=0.02 → τ=τ (normal, fast), Clamp ensures saturation at MF=0.02
                    S_i_max = FT(0.05) * M_norm  # vanish at 5% supersat for MF=0.02, linearly less for smaller MF
                    S_norm = clamp(one(FT) - S_i / max(S_i_max, eps(FT)), FT(0), FT(1)) # at S_i = S_i_max, S_norm = 0, so τ → ∞, at S_i = 0, S_norm = 1, so τ = τ
                    τ_fluc = τ / max(M_norm * S_norm, eps(FT)) # as MF→0, τ→∞; as MF→0.02, τ→τ]
                end



                if !use_cloak # Boost for downwelling particles
                    # The boost part (downwelling particles) are assumed to be dryer, so more likely to autoconvert. We calculate their part separately
                    # τ_boost = τ_sub_dep + (τ_orig - τ_sub_dep) * clamp((r_boost - r_thresh) / (FT(0.05) * r_thresh), FT(0), FT(1)) # tau gets fast as we approach r_thresh, peaks at 1.05 r_thresh...
                    f_boost = clamp((r_boost - r_thresh) / (FT(0.05) * r_thresh), FT(0), FT(1))  # normalized ramp
                    γ_boost = FT(0.05)  # adjust to control how fast τ drops
                    τ_boost = τ_sub_dep^(1 - f_boost^γ_boost)
                    R = FT(10) # how much to slowdown
                    M_norm = clamp(massflux / FT(0.02), FT(0), FT(1)) # saturate slowdown at MF = 0.02
                    τ_boost = min(τ_sub_dep, τ_boost * (one(FT) + (R - one(FT)) * M_norm)) # no S_N here since N is being converted out anyway
                end
            end

        else # subsaturated... We are losing q. We show in the docstring that evaporation raises <r> while acnv combats this by shrinking <r>, all the way till everything evaporates and autoconverts. This should go faster at higher subsat, and slower the closer to the target you are... though here the target r=r_thresh but q is going to 0 eventually.
            # τ = τ_sub_dep + (τ - τ_sub_dep) * clamp((r - r_thresh) / (FT(0.05) * r_thresh), FT(0), FT(1))
            f = clamp((r - r_target) / (FT(0.05) * r_target), FT(0), FT(1))  # normalized ramp
            γ = FT(0.05)  # adjust to control how fast τ drops
            τ = τ_sub_dep^(1 - f^γ) * τ^(f^γ) # As we get very close to the target, we do slow down as we rely on evaporation to raise <r> and drive acnv, though the ratio is between how much q needs to evaporate vs how much q needs to autoconvert are massively different, asymptote to the same timescale as a heuristic.

            # at S_i = 0, we'd like to slow down a lot since evaporation is not helping raise <r>, at S_i = -0.05 (5% subsat) we'd like to be at full speed since evaporation is strongly helping raise <r>. By S_i = -0.05, we want full speed
            # τ /= clamp(S_i / FT(0.05), FT(0), FT(1))
            γ_S = FT(0.25)
            fS = clamp((-S_i / FT(0.05))^γ_S, FT(0), one(FT))  # use -S_i so base is positive
            τ = τ_sub_dep^(one(FT) - fS) * τ^(fS)  # geometric blend, slow near S_i=0, ramps near max subsaturation

            # # redirect sedimentation
            # acnv_rate_other += dqdt_sed * clamp((r - r_target) / (1.2*r_thresh - r_target), FT(0), FT(1)) # re-direct sedimentation -- by r = r_max, redirecting all of it to acnv_thresh

        end


        λ_thresh = (μ + FT(1)) / r_thresh # this is the λ for the threshold radius, so we can use it to scale the acnv rate
        λ_target = (μ + FT(1)) / r_target
        _, λ = get_n0_lambda(microphys_params, q_type, q, ρ_a, N, μ; Dmax=Dmax)



        if !is_growing
            if λ < λ_target
                # No boost/no-boost separation here...

                # TEST
                # acnv_rate = (1 - (λ/λ_target)^3) * (q / τ) # q goes as λ^-4, this is variable N, fixed n_0 intercept
                # q_target = (λ/λ_target)^3 * q # q goes as λ^-4, this is variable N, fixed n_0 intercept

                # τ_eff = τ / (1 - (λ/λ_target)^3) # effective timescale [[ see docstring,     τ_eff = τ / [1 - (rₜₕ / <r>)³] > τ  ]] However, this goes away as we approach the threshold despite continued sublimation driving <r> up. We know  (1 - (λ/λ_target)^3) * (q / τ)  is the instantaneous rate at t=0... it just isn't representative of the entire timeste as evap keeps driving <r> up and massive amounts of q need to autoconvert to keep up.
                acnv_rate = q / τ
                # q_target = FT(0) # we go all the way to 0 in subsaturated regions

                q_below_r_max = get_fraction_below_or_above_thresh_from_qN(param_set, ice_type, q, N, ρ_a, 2*r_is; return_N = false, below_thresh = true, return_fraction = false, μ=μ, Dmax=FT(Inf)) # we'll only allow losing the particles above 2*r_max, at some point we'll run out of large particles to acnv and then <r> won't go up anymore...
                q_target = q_below_r_max * exp(-Δt / τ) # target at end of timestep
            else
                # We're not big enough but could get there...
                # Since our formulation has no lower bound, we don't want it to turn off suddenly... so we want to ramp into things... at λ = λ_target, we'd have acnv_rate = 0...
                # It's very hard to know the exchange rate between evaporation driving <r> up and acnv driving <r> down... so we just assume a smooth ramps. The λ at which we hit 0 gets lower with increasing subsat since evaporation is helping more...
                # At S_i = 0, we'd like to be off by r = 0.75 r_target, at S_i = -0.25 (5% subsat) we'd like to be off by r = 0.95 r_target. At 'off' τ → ∞ [or here just τ_sub_dep]
                λ_off = λ_target / (FT(0.75) + (FT(0.95) - FT(0.75)) * clamp(-S_i / FT(0.05), FT(0), one(FT))) # at S_i = 0, λ_off = λ_target / 0.75; at S_i = -0.05, λ_off = λ_target / 0.95
                γ_λ = FT(0.25)
                f = clamp((λ - λ_off) / (λ_target - λ_off), FT(0), one(FT))
                τ_eff = τ_sub_dep^(1 - f^γ_λ) * τ^(f^γ_λ)
                acnv_rate = q / τ_eff
                # q_target = FT(0) # we go all the way to 0 in subsaturated regions
                q_below_r_max = get_fraction_below_or_above_thresh_from_qN(param_set, ice_type, q, N, ρ_a, 2*r_is; return_N = false, below_thresh = true, return_fraction = false, μ=μ, Dmax=FT(Inf)) # we'll only allow losing the particles above 2*r_max, at some point we'll run out of large particles to acnv and then <r> won't go up anymore...
                q_target = q_below_r_max * exp(-Δt / τ) # target at end of timestep

            end
        else  # supersaturated S_i > 0
            # we separate boost and no boost regions, and treat them differently.


            if ((domain isa UpDomain) || (domain isa BulkDomain)) # just direct acnv... since we're supersaturated though, no runaway, just go to target
                if λ < λ_target
                    acnv_rate = (1 - (λ/λ_target)^3) * (q / τ) # q goes as λ^-4, this is variable N, fixed n_0 intercept
                    q_target = (λ/λ_target)^3 * q # q goes as λ^-4, this is variable N, fixed n_0 intercept
                end
            else
            
                if !use_cloak # Entertain more MF based fluctuations when not cloaking. Here we update acnv_rate_fluc
                    # Boost
                    _, λ_boost = get_n0_lambda(microphys_params, q_type, q_boost, ρ_a, N_i_boost, μ; Dmax=Dmax)
                    _, λ_no_boost = get_n0_lambda(microphys_params, q_type, q_no_boost, ρ_a, N_i_no_boost, μ; Dmax=Dmax)

                    # if do_print
                    #     @warn("Boost region: r_boost = $r_boost; N_i_boost = $N_i_boost; q_boost = $q_boost; λ_boost = $λ_boost; τ_boost = $τ_boost;")
                    #     @warn("No-boost region: r_no_boost = $r_no_boost; N_i_no_boost = $N_i_no_boost; q_no_boost = $q_no_boost; λ_no_boost = $λ_no_boost; τ = $τ;")
                    # end

                    if (λ_boost < λ_target) && (q_boost > 0)
                        acnv_rate += (1 - (λ_boost/λ_target)^3) * (q_boost / τ_boost) # q goes as λ^-4, this is variable N, fixed n_0 intercept (we're losing particles without new generation, so not fixed N)
                        q_target_boost = (λ_boost/λ_target)^3 * q_boost # q goes as λ^-4, this is variable N, fixed n_0 intercept

                        # if do_print
                        #     @warn("Boost region: λ_boost = $λ_boost; λ_target = $λ_target; q_boost = $q_boost; N_i_boost = $N_i_boost; acnv_rate contribution = $((1 - (λ_boost/λ_target)^3) * (q_boost / τ_boost)); q_target_boost = $q_target_boost; S_i = $S_i; τ_boost = $τ_boost; r_boost = $r_boost;")
                        # end

                    else
                        q_target_boost = FT(q_boost) # no acnv in boost region
                    end


                    if (λ_no_boost < λ_target) && (q_no_boost > 0)
                        acnv_rate += (1 - (λ_no_boost/λ_target)^3) * (q_no_boost / τ) # q goes as λ^-4, this is variable N, fixed n_0 intercept (we're losing particles without new generation, so not fixed N)
                        q_target_no_boost = (λ_no_boost/λ_target)^3 * q_no_boost # q goes as λ^-4, this is variable N, fixed n_0 intercept

                        # if do_print
                        #     @warn("No-boost region: λ_no_boost = $λ_no_boost; λ_target = $λ_target; q_no_boost = $q_no_boost; N_i_no_boost = $N_i_no_boost; acnv_rate contribution = $((1 - (λ_no_boost/λ_target)^3) * (q_no_boost / τ)); q_target_no_boost = $q_target_no_boost; S_i = $S_i; τ = $τ; r_no_boost = $r_no_boost;")
                        # end

                    else
                        q_target_no_boost = FT(q_no_boost) # no acnv in no boost region
                    end

                    q_target = q_target_boost + q_target_no_boost

                else # regular growth, stop at threshold
                   
                    if λ < λ_target
                        acnv_rate = ( 1 - (λ/λ_target)^3) * (q / τ) # q goes as λ^-4, this is variable N, fixed n_0 intercept
                        q_target = (λ/λ_target)^3 * q # q goes as λ^-4, this is variable N, fixed n_0 intercept
                    end
                end

                 # Fluctuations
                # acnv_rate_fluc = q/τ_fluc
                # This is hard because yes the target should go down but it's not clear how much... We need to estimate how much of the original mass i guess was susceptible to fluctuations into subsat...
                # Maybe one could use `get_qs_from_saturation_excesses()` and a closure but that's a bit complex...
                # We'll aim for q_target = q at S_i = 10%, and 0 at S_i = 0% 
                # q_target_fluc = q * clamp(S_i / FT(0.1), FT(0), one(FT))
                q_target_fluc = q * clamp(1 + S_i / FT(0.1), FT(0), one(FT))

                acnv_rate_fluc = (q - q_target_fluc) / τ_fluc # just setting acnv_rate to this (0) seemed to work ok lol
                acnv_rate_fluc = min(acnv_rate_fluc, (q - q_target)/Δt) # ensure we don't go slower than the main acnv rate

                # if do_print
                #     @warn("S_i = $S_i")
                #     f = 1 + S_i / FT(0.1)
                #     f_clamped = clamp(f, FT(0), one(FT))
                #     @warn("f = $f; f_clamped = $f_clamped;")
                #     @warn("q = $q;  q_target_fluc = $(q * f_clamped); diff = $(q - q * f_clamped)")
                #     println("--------------------------------")
                # end

                # acnv_rate = acnv_rate_other
                acnv_rate = dqdt_sed


                q_target = min(q_target, q_target_fluc) # take the lower of the two targets

            end
        end


        acnv_rate += acnv_rate_other # add any other acnv rate from redirected deposition
        # if rand() < 5e-6
        #     @error("We need to take out sub_dep_Above and crossing r_is from being redirected since it already is... maybe just subtract them from sub_Dep before passing in?")
        # end
        # acnv_rate += acnv_rate_fluc
        # acnv_rate = dqdt_sed * clamp((r - r_thresh) / (1.2*r_thresh - r_thresh), FT(0), FT(1))  # / (N/100) # re-direct sedimentation -- by r = r_max, redirecting all of it to acnv_thresh


        # Redirect sedimentation :: on when both (r > r_thresh) and subsaturated (S_i < 0)
        # Off at S_i = 0 and r = r_thresh, fully on by S_i = -0.05 and r = 1.05*r_thresh



        #=
            We should have a slowdown factor tied to TKE/QTvar/S_i_var.... LES seems to go over at -10% subsat automatically...

            They do : FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON :: see [ https://github.com/DOI-USGS/COAWST/blob/6419fc46d737b9703f31206112ff5fba65be400d/WRF/phys/module_mp_morr_two_moment.F#L1361C14-L1366C22 ]

            IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF

            This is probably why it was so hard for us to capture, they have both DCS = 2 r_is acnv happening and this, but this only comes on at q_i = 1e-8 which is borderline on RF09, and then the variability in QT on it as well.
        =#

        # scale down based on how much is actually below -10% subsat
        # calculate q_vap_sat_ice from qt and S_i
        q_vap_sat = qt / (one(FT) + S_i)
        q_subsat = q_vap_sat * 0.9  # 10% subsat
        qt_std = sqrt(qt_var)
        frac_sub10 = (qt_std > 0) ? Distributions.cdf(Distributions.Normal(qt, qt_std), q_subsat) : ((qt < q_subsat) ? one(FT) : zero(FT))
        
        qt_std_min = FT(.01) * qt
        frac_0_to_10 = Distributions.cdf(Distributions.Normal(qt, max(qt_std, qt_std_min)), q_vap_sat) - frac_sub10



        if (domain isa EnvDomain) || (domain isa CloakDownDomain) # only do this in env/cloak env
            f_r  = clamp((r - r_thresh) / (1.05*r_thresh - r_thresh), FT(0), FT(1))
            f_Si = clamp(-S_i / FT(0.1), FT(0), one(FT)) # -10% full force to match LES behavior
            acnv_rate = FT(0)
            # acnv_rate_addit = dqdt_sed * f_r * f_Si # / (N/100)
        else # 
            f_r  = clamp((r - r_thresh) / (1.05*r_thresh - r_thresh), FT(0), FT(1))
            f_Si = clamp(-S_i / FT(0.1), FT(0), one(FT))
            # acnv_rate_addit = dqdt_sed * f_r * f_Si
        end

        # acnv_rate = qisedflux * f_r * f_Si * frac_sub10 / τ # direct sedimentation redirection scaled by subsat fraction [[ pure qi_sed_flux is too susceptible to being 0 and have having no acnv despite being deeply subsat....  we want flux, but maybe enhanced if convergent since that implies more size sorting. ]]
        acnv_rate_addit = (dqdt_sed + q) * f_r * f_Si * frac_sub10 / τ # direct sedimentation redirection scaled by subsat fraction [[ pure qi_sed_flux is too susceptible to being 0 and have having no acnv despite being deeply subsat....  we want flux, but maybe enhanced if convergent since that implies more size sorting. ]]
        acnv_rate += acnv_rate_addit


        # This should probably be calibrated somehow?  or at least taper the timescale from S = 0 to S = -10? we also have the S scaling though so maybe it's fine... idk...
        acnv_rate += (dqdt_sed + q) * f_r * f_Si * frac_0_to_10 / (FT(5) * τ) # weaker acnv for 0-10% subsat region
        q_target = (1-frac_sub10) * q 




        # # TEST :: Scale down based on dNdz [[ 0 at dNdz = 0, 1 at dNdz = +2 ]], 
        # scale_factor = clamp(dN_i_dz / FT(2), FT(0), one(FT))
        # acnv_rate *= scale_factor
        # q_target += (q - q_target) * scale_factor

        # w_0 = 0.5
        # drdt = (r/3)*(-(1/q)*dqdt_dep + 1/N* dN_i_dz * w_0)
        # acnv_rate = max(drdt / 500, 0) * max(((r - r_thresh)/r_thresh)^3, 0) * N/.005 / (N/50)
        # # q_target = q - drdt * (3 * q / r) * Δt

        # # acnv_rate /= (N/100)


    
    end
    return resolve_nan(acnv_rate, FT(0)), q_target # also return q_target for limiters
end





# ====================================================================================== #




"""
    acnv_mass_fraction(mu, mean_r, r_is; verbose=false)

Estimate the fraction of collision mass producing particles with radius above a threshold `r_is`,
assuming the initial particle size distribution follows a gamma distribution:

    n(r) ∝ r^μ * exp(-λ r)

with shape parameter μ and scale parameter λ = (μ + 1) / mean_r.

# Mathematical derivation:

1. The radius distribution  n(r)  is gamma with parameters  α = μ + 1 , scale  β = \frac{\text{mean}_r}{α} :

   n(r) = \frac{r^{α - 1} e^{-r/β}}{β^α Γ(α)}
   

2. Collision output radius distribution is assumed to be governed by the volume  v = r^3 .

3. Calculate the first and second moments of volume  v :

   E[v] = {E}[r^3] = β^3 \frac{Γ(α + 3)}{Γ(α)}
   E[v^2] = {E}[r^6] = β^6 \frac{Γ(α + 6)}{Γ(α)}

   with variance:

   
   Var[v] = E[v^2] - (E[v])^2

4. Fit a gamma distribution for  v  with shape  k_v = \frac{E[v]^2}{Var[v]}  and scale  θ_v = \frac{Var[v]}{E[v]} .

5. The fraction of collision mass with output radius greater than threshold  r_{is}  is the tail probability of volume  v > v_{th} = r_{is}^3 :

   
   F_{acnv} = \frac{Γ(2 k_v + 1, v_{th} / θ_v)}{2 k_v Γ(2 k_v)}
   

where  Γ(s,x)  is the upper incomplete gamma function.

# Arguments
- `mu`: shape parameter of the initial radius gamma distribution.
- `mean_r`: mean radius of the initial distribution (meters).
- `r_is`: threshold radius (meters).
- `verbose`: print intermediate variables if true.

# Returns
- Fraction of collision mass above threshold radius.

# Notes
- Approximates collision output volume distribution as gamma.
- Accuracy decreases for very peaked distributions (large μ).
"""


function acnv_mass_fraction(
    mean_r::FT,
    r_is::FT;
    μ::FT = FT(0)    
) where {FT <: Real}
    if isnan(μ)
        μ = FT(0)
    end
    α = μ + 1
    β = mean_r / α  # scale parameter for Gamma distribution of r

    # Moments of r^3:
    E_r3 = β^3 * CM1.SF.gamma(α + 3) / CM1.SF.gamma(α)
    E_r6 = β^6 * CM1.SF.gamma(α + 6) / CM1.SF.gamma(α)
    Var_r3 = E_r6 - E_r3^2

    # Fit Gamma distribution to v = r^3
    k_v = E_r3^2 / Var_r3
    θ_v = Var_r3 / E_r3

    v_thresh = r_is^3

    # Numerator: upper incomplete gamma Γ(s,x) = CM1.SF.gamma(s,x)
    numer = CM1.SF.gamma(2*k_v + 1, v_thresh / θ_v)
    denom = 2 * k_v * CM1.SF.gamma(2*k_v)

    F_acnv = numer / denom

    # if verbose
    #     println("Input mean radius = $(mean_r) m")
    #     println("Threshold radius = $(r_is) m")
    #     println("Fitted Gamma params for v=r^3: k=$(k_v), θ=$(θ_v)")
    #     println("Estimated mass fraction above threshold: $(F_acnv)")
    # end

    return F_acnv
end

# ====================================================================================== #