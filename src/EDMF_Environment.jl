function microphysics!(
    en_thermo::Union{SGSMean, SGSMeanWQuadratureAdjustedNoneqMoistureSources},
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    moisture_model::AbstractMoistureModel,
    precip_model::AbstractPrecipitationModel,
    cloud_sedimentation_model::AbstractCloudSedimentationModel,
    rain_formation_model::AbstractRainFormationModel,
    snow_formation_model::AbstractSnowFormationModel,
    Δt::Real,
    param_set::APS,
    use_fallback_tendency_limiters::Bool
)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    tendencies_pr = center_tendencies_precipitation(state)
    aux_en = center_aux_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_en_f = face_aux_environment(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    aux_tc_f = face_aux_turbconv(state)
    prog_pr = center_prog_precipitation(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    ts_env = center_aux_environment(state).ts
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_en_sat = aux_en.sat
    aux_en_unsat = aux_en.unsat
    precip_fraction = compute_precip_fraction(edmf, state)

    if (moisture_model isa NonEquilibriumMoisture) || (edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean) || ((edmf.area_partition_model isa CoreCloakAreaPartitionModel) && edmf.area_partition_model.apply_cloak_to_condensate_formation)
        if (edmf.area_partition_model isa CoreCloakAreaPartitionModel)
            if edmf.area_partition_model.confine_all_downdraft_to_cloak
                w::CC.Fields.Field = similar(aux_en.T)
                @. w = FT(0) # we confine the downdraft (which will be drier) to the cloak region, so we can set w = 0 in the environment
            else # not confined to cloak, but also we do need some to be converted to updraft, maybe still go with 0. That biases us towards respecting supersat generation in the cloak updraft...
                w = similar(aux_en.T)
                @. w = FT(0)
            end
            w_cloak_up = Ic.(aux_en_f.w_cloak_up)
            w_cloak_dn = Ic.(aux_en_f.w_cloak_dn)
        else
            w = Ic.(aux_en_f.w)
        end
    else
        # don't think we need w for EquilibriumMoisture... we're using w on faces rn...
    end

    if moisture_model isa NonEquilibriumMoisture
        nonequilibrium_moisture_scheme = moisture_model.scheme # my new addition
        # moisture_sources_limiter = edmf.tendency_limiters.moisture_sources_limiter
        moisture_sources_limiter = get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources), use_fallback_tendency_limiters)
    else
        nonequilibrium_moisture_scheme = moisture_model.scheme
    end

    # should we consider using  Δt + ε for calculating microphysical process rates?
    # - pro: could decrease floating point errors in tendency limited dt calculation, favoring slower rates to ensure rates calculated from Δt do not deplete faster than Δt when backed out through division later (maybe important for our N_remaining_violation counter...)
    # - con: could maybe lead to leftovers from inexact removal of tendencies in the specirfied timestep (though the leftover should maybe be small enough to be elided by our eps(FT) size filter...)
    # - con: for very short Δt, could significantly change the rates (counterpoint -- eps(FT) is a ridiculously fast rate for anything...)
    ε = eps(FT)

    # ======================================================================== # To Do: Separate sedimentation out so it can be called w/ either sgs or mean...
    # sedimentation (should this maybe be a grid mean tendency?)
    # [ do sed first bc this should happen before precip_formation so we can know sed-driven aggregatory autoconversion, and it doesnt depend on anything below. ]
    if edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean  
        
        #= 
            TEST : The argument `doing advection and sedimentation separately bc the dispersion is good for size distribution` argument is kinda weak, in the limit we could double our condensate production
            Instead, we'll trial here calculating total advection via relative w (w_sed - w) and then subtracting a version based on only advection

            The advective part is usually calculated separately. For the updraft, it's in EDMF_functions.jl, but is denominated in `relative to grid mean` terms by default
            There we only store the grid mean total flux, and then take the gradient to get a tendecy for storage, but the raw sedimentation tendency then is biased and dispersive, so this should perhaps be better without relying code changes everywhere.
            
            Note though that we cannot in princile guarantee that the up/env tendencies this generate sum to exactly the same as the tendency from the gradient of the grid mean fluxes (especially with our transfer heuristic)

            Also for the `advective only` parts, those calculate `to_other` ignoring the background w, so `advection only` should cause no `to_other` atall.

        =#

        mph_sed_ice, mph_sed_ice_other = calculate_sedimentation_sources(
            param_set,
            ρ_c,
            aux_en.q_ice,
            aux_en.term_vel_ice,
            aux_en_f.w, # could just be 0 bc we don't use it w/ use_relative_w = false, bc advection is handled elsewhere and we argue that diffusion is maybe ok since we don't model the size distribution
            aux_en.area,
            grid,
            differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
            grid_mean = false,
            use_relative_w = false,
        ) # should this be a grid mean tendency? 

        mph_sed_liq, mph_sed_liq_other = calculate_sedimentation_sources(
            param_set,
            ρ_c,
            aux_en.q_liq,
            aux_en.term_vel_liq,
            aux_en_f.w, # could just be 0 bc we don't use it w/ use_relative_w = false, bc advection is handled elsewhere and we argue that diffusion is maybe ok since we don't model the size distribution
            aux_en.area,
            grid,
            differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
            grid_mean = false,
            use_relative_w = false,
        ) # should this be a grid mean tendency? 


        #=
            We divided out ρ so, instead of ρaq, these are tendencies in just aq... we did this to avoid dividing by then multiplying by a... but for storage we really should right?
        =#

        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)

        # Store sedimentation
        @inbounds for k in real_center_indices(grid)
            # === sed ======================== [ do before precip_formation so we can know sed-driven aggregatory autoconversion, and it doesnt depend on anything below. ]
            if edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean        
                ql_tendency_sedimentation = mph_sed_liq[k].q_tendency
                qi_tendency_sedimentation = mph_sed_ice[k].q_tendency
                qt_tendency_sedimentation = ql_tendency_sedimentation + qi_tendency_sedimentation
                aux_en.ql_tendency_sedimentation[k] += ql_tendency_sedimentation # these get added to gm in compute_gm
                aux_en.qi_tendency_sedimentation[k] += qi_tendency_sedimentation
                aux_en.qt_tendency_sedimentation[k] += qt_tendency_sedimentation # used in dycore.jl (= not += , cause this doesnt seem to get reset every iteration?) (fixed in update_aux)
    
                if !isfinite(qi_tendency_sedimentation)
                    @error "qi_tendency_sedimentation is NaN at k = $k; aux_en.qi_tendency_sedimentation[k] = $(aux_en.qi_tendency_sedimentation[k]); aux_en.qi_tendency_sedimentation[k] = $(aux_en.qi_tendency_sedimentation[k])"
                end
    
                Π_m = TD.exner(thermo_params, ts_env[k])
                c_pm = TD.cp_m(thermo_params, ts_env[k])
                θ_liq_ice_tendency_sedimentation =
                    1 / Π_m / c_pm * (L_v0 * ql_tendency_sedimentation + L_s0 * qi_tendency_sedimentation)
                aux_en.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # adapted from microphysics_coupling.jl: precipitation_formation() | (= not += caue these don't seem to get reset every iteration?)
    
                # sedimentation loss into updraft (should this be allowed lol)
                # How do we know when these get set to 0 and can be used again?
                if aux_en.area[k] < 1 # we have some updrafts
                    N_up = n_updrafts(edmf)
                    @inbounds for i in 1:N_up
                        ql_tendency_sedimentation_other = mph_sed_liq_other[k].q_tendency .* (aux_up[i].area[k] ./ aux_bulk.area[k])
                        qi_tendency_sedimentation_other = mph_sed_ice_other[k].q_tendency .* (aux_up[i].area[k] ./ aux_bulk.area[k])
                        qt_tendency_sedimentation_other = ql_tendency_sedimentation_other + qi_tendency_sedimentation_other
                        θ_liq_ice_tendency_sedimentation_other =
                            1 / Π_m / c_pm *
                            (L_v0 * ql_tendency_sedimentation_other + L_s0 * qi_tendency_sedimentation_other)
                        aux_up[i].ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                        aux_up[i].qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                        aux_up[i].qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                        aux_up[i].θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other
                        aux_bulk.ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                        aux_bulk.qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                        aux_bulk.qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                        aux_bulk.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other

                        # store contributions to other
                        aux_up[i].ql_tendency_sedimentation_other[k] += ql_tendency_sedimentation_other # for storage
                        aux_up[i].qi_tendency_sedimentation_other[k] += qi_tendency_sedimentation_other # for storage
                        aux_up[i].qt_tendency_sedimentation_other[k] += qt_tendency_sedimentation_other # for storage
                        aux_up[i].θ_liq_ice_tendency_sedimentation_other[k] += θ_liq_ice_tendency_sedimentation_other # for storage

                        aux_bulk.ql_tendency_sedimentation_other[k] += ql_tendency_sedimentation_other # for storage
                        aux_bulk.qi_tendency_sedimentation_other[k] += qi_tendency_sedimentation_other # for storage
                        aux_bulk.qt_tendency_sedimentation_other[k] += qt_tendency_sedimentation_other # for storage
                        aux_bulk.θ_liq_ice_tendency_sedimentation_other[k] += θ_liq_ice_tendency_sedimentation_other # for storage
                    end
                end
            end
            # ============================================== #
        end
        #
    end
    # ======================================================================== #

    reweight_processes_for_grid::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_processes_for_grid, false)
    reweight_extrema_only::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false)


    @inbounds for k in real_center_indices(grid)
        # condensation
        ts = ts_env[k]
        ρ = TD.air_density(thermo_params, ts)
        if moisture_model isa NonEquilibriumMoisture

            # Apply massflux correction to qt...
            # - not sure how to do, we would still ned multiple calls to noneq_moisture_sources() right?


            if en_thermo isa SGSMean

                do_cloak = (edmf.area_partition_model isa CoreCloakAreaPartitionModel) && (edmf.area_partition_model.apply_cloak_to_condensate_formation) && (aux_bulk.area[k] > edmf.minimum_area)
                # (region_area, w_region, T, ts, cloak_up) 
                if !do_cloak
                    regions = ((aux_en.area[k], w[k], aux_en.T[k], ts_env[k], Env),)
                else
                    w_en_remaining = (edmf.area_partition_model.confine_all_downdraft_to_cloak) ? FT(0) : w[k] # if we confine all downdraft to cloak, then the remaining env has no vertical velocity
                    regions = (
                        (aux_en.a_cloak_up[k], w_cloak_up[k], aux_en.T_cloak_up[k], aux_en.ts_cloak_up[k], CloakUp),
                        (aux_en.a_cloak_dn[k], w_cloak_dn[k], aux_en.T_cloak_dn[k], aux_en.ts_cloak_dn[k], CloakDown),
                        (aux_en.a_en_remaining[k], w_en_remaining, aux_en.T_en_remaining[k], aux_en.ts_en_remaining[k], EnvRemaining),
                        )
                end


                mph_neq = null_NoneqMoistureSources(FT; fill_value = FT(0)) # initialize...
                mph_neq_other = null_OtherMicrophysicsSources(FT; fill_value = FT(0)) # initialize...
                mph_precip = null_PrecipitationSources(FT; fill_value = FT(0)) # initialize...
                for (region_area, w_region, T, ts, region) in regions

                    ρ = TD.air_density(thermo_params, ts)
                    if (region isa EnvDomain) || (region isa EnvRemainingDomain)
                        q_vap_sat_liq = aux_en.q_vap_sat_liq[k]
                        q_vap_sat_ice = aux_en.q_vap_sat_ice[k]
                        p = aux_en.p[k]
                    else
                        q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
                        q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
                        p = TD.air_pressure(thermo_params, ts)
                    end

                    q_here = TD.PhasePartition(thermo_params, ts) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work

                    # calculate local timescales by scaling the env timescale. We know τ = 1/(4πDNr) ~ 1/(4πDqr) assuming fixed <r>. so τ_here = τ_env * (q_env/q_here) but that blows up when q_here -> 0. so we do the following scaling instead:
                    τ_liq_here = aux_en.τ_liq[k] * (aux_en.q_liq[k] + eps(FT)) / (q_here.liq + eps(FT))
                    τ_ice_here = aux_en.τ_ice[k] * (aux_en.q_ice[k] + eps(FT)) / (q_here.ice + eps(FT))


                    # limit exterior tendencies to give microphysics first swing
                    dqvdt = aux_en.dqvdt[k]
                    dTdt = aux_en.dTdt[k]
                    dqvdt = max(dqvdt, FT(0)) # give microphysics first right of refusal on vapor usage
                    dTdt = min(dTdt, FT(0)) # give microphysics first right of refusal on cooling usage
                    # dqvdt = FT(0)
                    # # add advection, diffusion
                    # dqvdt += (aux_tc.diffusive_tendency_qt[k] - aux_tc.diffusive_tendency_ql[k] - aux_tc.diffusive_tendency_qi[k]) / ρ # diffusion of qt into/outof env goes into vapor
                    # dqvdt += (aux_tc.massflux_tendency_qt[k] - aux_tc.massflux_tendency_ql[k] - aux_tc.massflux_tendency_qi[k]) / ρ # massflux of qt into/outof env goes into vapor
                    # dTdt = FT(0)

                    if ((param_set.user_params.use_convective_tke) || (param_set.user_params.use_convective_tke_production_only)) && (region isa EnvDomain)
                        # assume the condensate is in the tke part going up KE = 1/2 ρ w^2. critical to generating condensate
                        w_noneq = sqrt(2 * aux_en.tke[k] / ρ)
                        
                        # Now, we assume liq is biased with covariance, when we add this eddy mixing so we go up and down 1SD, but assume liq is tied to +1 SD
                    else
                        w_noneq = w[k]
                    end

                    if !iszero(param_set.user_params.condensate_qt_SD) # maybe you'd want different values for liq and ice but theyre irreparably tied together in MM2015 microphysics, etc.
                        # TODO: Consider turning this off for when we are using cloaks.... also set an error in the quadrature branches
                        q_tot_sd = sqrt(aux_en.QTvar[k])
                        qt = q_here.tot + param_set.user_params.condensate_qt_SD * q_tot_sd
                        θ = TD.liquid_ice_pottemp(thermo_params, ts)
                        ts = (edmf.moisture_model isa NonEquilibriumMoisture) ? thermo_state_pθq(param_set, p_c[k], θ, qt, q_here.liq, q_here.ice) : thermo_state_pθq(param_set, p_c[k], θ, qt)
                    end


                    mph_neq_here = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, region_area, ρ_c[k], aux_en.p[k], T, Δt + ε, ts, w_noneq, q_vap_sat_liq, q_vap_sat_ice, dqvdt, dTdt, τ_liq_here, τ_ice_here)

                    if region isa EnvRemainingDomain
                        aux_en.qi_tendency_sub_dep_en_remaining[k] = mph_neq_here.qi_tendency
                    elseif region isa CloakUpDomain
                        aux_en.qi_tendency_sub_dep_cloak_up[k] = mph_neq_here.qi_tendency
                    elseif region isa CloakDownDomain
                        aux_en.qi_tendency_sub_dep_cloak_dn[k] = mph_neq_here.qi_tendency
                    end
                    


                    # if at an extrema, we don't know where the peak is so we'll reweight based on probability of where it could be in between.
                    if reweight_processes_for_grid
                        w_region_full = if region isa EnvOrUp
                            w # is already on centers
                        elseif region isa EnvRemainingDomain
                            edmf.area_partition_model.confine_all_downdraft_to_cloak ? (w .* FT(0)) : w
                        elseif region isa CloakUp
                            w_cloak_up
                        elseif region isa CloakDown
                            w_cloak_dn
                        else
                            error("Unknown region type")
                        end
                            
                        mph_neq_here = reweight_noneq_moisture_sources_for_grid(k, grid, param_set, thermo_params, aux_en, aux_en_f, mph_neq, nonequilibrium_moisture_scheme, moisture_sources_limiter, Δt, ρ_c, p_c, w_region_full, aux_en.dqvdt, aux_en.dTdt; reweight_extrema_only = reweight_extrema_only, region = region)
                    end
                    mph_neq += mph_neq_here * (region_area / aux_en.area[k]) # weight by area fraction
                    
                    mph_neq_other_here = other_microphysics_processes(
                        param_set,
                        moisture_model,
                        # moisture_model.heterogeneous_ice_nucleation,
                        nonequilibrium_moisture_scheme,
                        moisture_sources_limiter,
                        region_area,
                        ρ_c[k],
                        p,
                        T,
                        Δt,
                        ts,
                        w_region,
                        mph_neq.ql_tendency,
                        mph_neq.qi_tendency;
                        N_l = aux_en.N_l[k],
                        N_i = aux_en.N_i[k],
                    )
                    mph_neq_other += mph_neq_other_here * (region_area / aux_en.area[k]) # weight by area fraction

                    q = TD.PhasePartition(thermo_params, ts) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
                    S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
                    N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, T, ρ, w_region)

                    mph_precip_here = precipitation_formation(
                        param_set,
                        moisture_model,
                        precip_model,
                        cloud_sedimentation_model,
                        rain_formation_model,
                        snow_formation_model,
                        edmf.area_partition_model,
                        prog_pr.q_rai[k],
                        prog_pr.q_sno[k],
                        aux_en.q_liq[k],
                        aux_en.q_ice[k],
                        aux_en.N_i[k],
                        aux_en.term_vel_ice[k],
                        aux_tc.term_vel_rain[k],
                        aux_tc.term_vel_snow[k],
                        region_area,
                        ρ_c[k],
                        p,
                        T,
                        Δt,
                        ts,
                        mph_neq_here.ql_tendency, # this is the sub-deposition contribution to precip formation
                        mph_neq_here.qi_tendency, # this is the sub-deposition contribution to precip formation
                        (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation[k] / aux_en.area[k] : FT(0), # this is the sedimentation contribution to precip formation
                        precip_fraction,
                        # edmf.tendency_limiters.precipitation_tendency_limiter,
                        get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                        # aux_en_2m.tke.dissipation[k]
                        aux_en.tke[k],
                        aux_en.N_i_no_boost[k],
                        S_i,
                        # aux_en.τ_ice[k], # keep same in all regions for simplicity
                        τ_ice_here,
                        aux_en.dN_i_dz[k],
                        aux_en.dqidz[k],
                        N_INP,
                        aux_tc.massflux[k],
                        region,
                    )
                    mph_precip += mph_precip_here * (region_area / aux_en.area[k]) # weight by area fraction
                end

            elseif en_thermo isa SGSMeanWQuadratureAdjustedNoneqMoistureSources # Doesn't seem to do much in the end, variance is pretty small...

                if edmf.area_partition_model isa CoreCloakAreaPartitionModel && (edmf.area_partition_model.apply_cloak_to_condensate_formation) && (aux_bulk.area[k] > edmf.minimum_area)
                    error("SGSMeanWQuadratureAdjustedNoneqMoistureSources with CoreCloakAreaPartitionModel not implemented yet")
                else
                    # use quadrature to adjust noneq moisture sources based on subgrid variability in qt and θl
                    region = Env
                end

                if ((param_set.user_params.use_convective_tke) || (param_set.user_params.use_convective_tke_production_only)) && (region isa EnvDomain)
                    # assume the condensate is in the tke part going up KE = 1/2 ρ w^2. critical to generating condensate
                    w_noneq = sqrt(2 * aux_en.tke[k] / ρ)
                    
                    # Now, we assume liq is biased with covariance, when we add this eddy mixing so we go up and down 1SD, but assume liq is tied to +1 SD
                else
                    w_noneq = w[k]
                end

                # limit exterior tendencies to give microphysics first swing
                dqvdt = aux_en.dqvdt[k]
                dTdt = aux_en.dTdt[k]
                dqvdt = max(dqvdt, FT(0))
                dTdt = min(dTdt, FT(0))

                quadrature_type = en_thermo.quadrature_type
                quad_order = quadrature_order(en_thermo)
                χ = en_thermo.a
                weights = en_thermo.w

                # We just are predicting mph_neq which is 2 things
                sqpi_inv = FT(1 / sqrt(π))
                sqrt2 = FT(sqrt(2))
                # Epsilon defined per typical variable fluctuation
                qt′qt′ = aux_en.QTvar[k]
                qt_mean = aux_en.q_tot[k]
                θl′θl′ = aux_en.Hvar[k]
                θl_mean = aux_en.θ_liq_ice[k]
                θl′qt′ = aux_en.HQTcov[k]
                q_liq = aux_en.q_liq[k]
                q_ice = aux_en.q_ice[k]

                eps_q = (qt_mean ≈ FT(0)) ? eps(FT) : (eps(FT) * qt_mean)
                eps_θ = eps(FT)
                if quadrature_type isa LogNormalQuad
                    # Lognormal parameters (ν, s) from mean and variance
                    ν_q = log(qt_mean^2 / max(sqrt(qt_mean^2 + qt′qt′), eps_q))
                    ν_θ = log(θl_mean^2 / sqrt(θl_mean^2 + θl′θl′))
                    s_q = sqrt(log(qt′qt′ / max(qt_mean, eps_q)^2 + 1))
                    s_θ = sqrt(log(θl′θl′ / θl_mean^2 + 1))

                    # Enforce Cauchy-Schwarz inequality, numerically stable compute
                    corr = θl′qt′ / max(sqrt(qt′qt′), eps_q)
                    corr = max(min(corr / max(sqrt(θl′θl′), eps_θ), 1), -1)

                    # Conditionals
                    s2_θq = log(corr * sqrt(θl′θl′ * qt′qt′) / θl_mean / max(qt_mean, eps_q) + 1)
                    s_c = sqrt(max(s_θ^2 - s2_θq^2 / max(s_q, eps_q)^2, 0))

                elseif quadrature_type isa GaussianQuad
                    # limit σ_q to prevent negative qt_hat
                    σ_q_lim = -qt_mean / (sqrt2 * χ[1])
                    σ_q = min(sqrt(qt′qt′), σ_q_lim)
                    σ_θ = sqrt(θl′θl′)

                    # Enforce Cauchy-Schwarz inequality, numerically stable compute
                    corr = θl′qt′ / max(σ_q, eps_q)
                    corr = max(min(corr / max(σ_θ, eps_θ), 1), -1)

                    # Conditionals
                    σ_c = sqrt(max(1 - corr * corr, 0)) * σ_θ
                end

                mph_neq = null_NoneqMoistureSources(FT; fill_value = FT(0)) # initialize...
                mph_neq_other = null_OtherMicrophysicsSources(FT; fill_value = FT(0)) # initialize...
                @inbounds for m_q in 1:quad_order
                    if quadrature_type isa LogNormalQuad
                        qt_hat = exp(ν_q + sqrt2 * s_q * χ[m_q])
                        ν_c = ν_θ + s2_θq / max(s_q, eps_q)^2 * (log(qt_hat) - ν_q)
                    elseif quadrature_type isa GaussianQuad
                        qt_hat = qt_mean + sqrt2 * σ_q * χ[m_q]
                        μ_c = θl_mean + sqrt2 * corr * σ_θ * χ[m_q]
                    end

                    @inbounds for m_h in 1:quad_order
                        if quadrature_type isa LogNormalQuad
                            h_hat = exp(ν_c + sqrt2 * s_c * χ[m_h])
                        elseif quadrature_type isa GaussianQuad
                            h_hat = μ_c + sqrt2 * σ_c * χ[m_h]
                        end
                        ts_hat = thermo_state_pθq(param_set, p_c[k], h_hat, qt_hat, aux_en.q_liq[k], aux_en.q_ice[k])
                        T = TD.air_temperature(thermo_params, ts_hat)

                        q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
                        q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
                        mph_neq_quad = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, aux_en.area[k], ρ_c[k], aux_en.p[k], T, Δt + ε, ts_hat, w_noneq, q_vap_sat_liq, q_vap_sat_ice, dqvdt, dTdt, aux_en.τ_liq[k], aux_en.τ_ice[k])

                        # if at an extrema, we don't know where the peak is so we'll reweight based on probability of where it could be in between.
                        if reweight_processes_for_grid
                            mph_neq_quad = reweight_noneq_moisture_sources_for_grid(k, grid, param_set, thermo_params, aux_en, aux_en_f, mph_neq, nonequilibrium_moisture_scheme, moisture_sources_limiter, Δt, ρ_c, p_c, w, aux_en.dqvdt, aux_en.dTdt; reweight_extrema_only = reweight_extrema_only)
                        end

                        mph_neq += mph_neq_quad * weights[m_h] * sqpi_inv * weights[m_q] * sqpi_inv
                        

                        #                         
                        mph_neq_other_quad = other_microphysics_processes(
                            param_set,
                            moisture_model,
                            # moisture_model.heterogeneous_ice_nucleation,
                            nonequilibrium_moisture_scheme,
                            moisture_sources_limiter,
                            aux_en.area[k],
                            ρ_c[k],
                            aux_en.p[k],
                            T,
                            Δt,
                            ts_hat,
                            w[k],
                            mph_neq_quad.ql_tendency,
                            mph_neq_quad.qi_tendency;
                            N_l = aux_en.N_l[k],
                            N_i = aux_en.N_i[k],
                        )
                        mph_neq_other += mph_neq_other_quad * weights[m_h] * sqpi_inv * weights[m_q] * sqpi_inv
                    end
                end

                q = TD.PhasePartition(thermo_params, ts_env[k]) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
                S_i = TD.supersaturation(thermo_params, q, ρ_c[k], aux_en.T[k], TD.Ice())
                N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, aux_en.T[k], ρ_c[k], w[k])
                mph_precip = precipitation_formation( # we leave this out of quadrature for now...
                    param_set,
                    moisture_model,
                    precip_model,
                    cloud_sedimentation_model,
                    rain_formation_model,
                    snow_formation_model,
                    edmf.area_partition_model,
                    prog_pr.q_rai[k],
                    prog_pr.q_sno[k],
                    aux_en.q_liq[k],
                    aux_en.q_ice[k],
                    aux_en.N_i[k],
                    aux_en.term_vel_ice[k],
                    aux_tc.term_vel_rain[k],
                    aux_tc.term_vel_snow[k],
                    aux_en.area[k],
                    ρ_c[k],
                    aux_en.p[k],
                    aux_en.T[k],
                    Δt,
                    ts,
                    mph_neq.ql_tendency, # this is the sub-deposition contribution to precip formation
                    mph_neq.qi_tendency, # this is the sub-deposition contribution to precip formation
                    (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation[k] / aux_en.area[k] : FT(0), # this is the sedimentation contribution to precip formation
                    precip_fraction,
                    # edmf.tendency_limiters.precipitation_tendency_limiter,
                    get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                    # aux_en_2m.tke.dissipation[k]
                    aux_en.tke[k],
                    aux_en.N_i_no_boost[k],
                    S_i,
                    aux_en.τ_ice[k],
                    aux_en.dN_i_dz[k],
                    aux_en.dqidz[k],
                    N_INP,
                    aux_tc.massflux[k],
                    Env,
                )

            else
                error("en_thermo of type $(typeof(en_thermo)) not supported in this function")
            end

            # update the tendencies
            aux_en.ql_tendency_noneq[k] = mph_neq.ql_tendency * aux_en.area[k]
            aux_en.qi_tendency_noneq[k] = mph_neq.qi_tendency * aux_en.area[k] # add in the sub-deposition conribution to precip formation
            aux_en.ql_tendency_cond_evap[k] = aux_en.ql_tendency_noneq[k] # for storage
            aux_en.qi_tendency_sub_dep[k] = aux_en.qi_tendency_noneq[k] # for storage

            if any(!isfinite, (mph_neq_other.ql_tendency, mph_neq_other.qi_tendency))
                @error "Other microphysics processes returned non-finite values: $mph_neq_other; from inputs ts = $ts; w = $w[k]; mph_neq.ql_tendency = $(mph_neq.ql_tendency); mph_neq.qi_tendency = $(mph_neq.qi_tendency); ρ_c = $(ρ_c[k]); aux_en.area = $(aux_en.area[k])"
                error("Other microphysics processes returned non-finite values")
            end

            aux_en.ql_tendency_noneq[k] += mph_neq_other.ql_tendency * aux_en.area[k] # is storing them w/ noneq best? should i add another one...? would need to tie it into everywhere else as well in dycore.jl etc... but wouldn't be too bad
            aux_en.qi_tendency_noneq[k] += mph_neq_other.qi_tendency * aux_en.area[k]

            aux_en.qi_tendency_hom_frz[k] = mph_neq_other.qi_tendency_homogeneous_freezing * aux_en.area[k] # for storage
            aux_en.qi_tendency_het_frz[k] = mph_neq_other.qi_tendency_heterogeneous_freezing * aux_en.area[k] # for storage
            aux_en.qi_tendency_het_nuc[k] = mph_neq_other.qi_tendency_heterogeneous_icenuc * aux_en.area[k] # for storage
            aux_en.qi_tendency_melt[k] = mph_neq_other.qi_tendency_melting * aux_en.area[k] # for storage
            
        else # EquilibriumMoisture


            do_cloak = (edmf.area_partition_model isa CoreCloakAreaPartitionModel) && (edmf.area_partition_model.apply_cloak_to_condensate_formation) && (aux_bulk.area[k] > edmf.minimum_area)
            # (region_area, T, ts, cloak_up) 
            if !do_cloak
                regions = ((aux_en.area[k], aux_en.T[k], ts_env[k]),)
            else
                regions = (
                    (aux_en.a_cloak_up[k], aux_en.T_cloak_up[k], aux_en.ts_cloak_up[k]),
                    (aux_en.a_cloak_dn[k], aux_en.T_cloak_dn[k], aux_en.ts_cloak_dn[k]),
                    (aux_en.a_en_remaining[k], aux_en.T[k], ts_env[k]),
                    )
            end


            mph_precip = null_PrecipitationSources(FT; fill_value = FT(0)) # initialize...
            for (region_area, T, ts) in regions
                # Compute microphysics tendencies for each region

                q = TD.PhasePartition(thermo_params, ts) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
                ρ = TD.air_density(thermo_params, ts)
                p = TD.air_pressure(thermo_params, ts)
                S_i = TD.supersaturation(thermo_params, q, ρ_c[k], ρ, TD.Ice())
                N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, T, ρ, FT(0)) # Don't need w since it'll just call Cooper)

                mph_precip_here = precipitation_formation(
                    param_set,
                    moisture_model,
                    precip_model,
                    cloud_sedimentation_model,
                    rain_formation_model,
                    snow_formation_model,
                    edmf.area_partition_model,
                    prog_pr.q_rai[k],
                    prog_pr.q_sno[k],
                    aux_en.q_liq[k],
                    aux_en.q_ice[k],
                    aux_en.N_i[k],
                    aux_en.term_vel_ice[k],
                    aux_tc.term_vel_rain[k],
                    aux_tc.term_vel_snow[k],
                    aux_en.area[k],
                    ρ_c[k],
                    p,
                    T,
                    Δt,
                    ts,
                    precip_fraction,
                    # edmf.tendency_limiters.precipitation_tendency_limiter,
                    get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                    # aux_en_2m.tke.dissipation[k],
                    aux_en.tke[k]
                    ;
                    N_i_no_boost = aux_en.N_i_no_boost[k],
                    qi_tendency_sed = (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation[k] / aux_en.area[k] : FT(0), # this is the sedimentation contribution to precip formation
                    S_i = S_i,
                    τ_sub_dep = aux_en.τ_ice[k],
                    dN_i_dz = aux_en.dN_i_dz[k],
                    dqidz = aux_en.dqidz[k],
                    N_INP = N_INP,
                    massflux = aux_tc.massflux[k],
                    domain = Env,
                )
                mph_precip += mph_precip_here * (region_area / aux_en.area[k]) # weight by area fraction
            end
        end


        # Update cloud fraction, q_liq, q_ice
        if moisture_model isa EquilibriumMoisture
            if (edmf.area_partition_model isa CoreCloakAreaPartitionModel)
                for (region_area,  ts, region) in (
                    (aux_en.a_cloak_up[k], aux_en.ts_cloak_up[k], CloakUp),
                    (aux_en.a_cloak_dn[k], aux_en.ts_cloak_dn[k], CloakDown),
                    (aux_en.a_en_remaining[k], aux_en.ts[k], EnvRemaining)
                )

                    aux_en.q_liq[k] = FT(0)
                    aux_en.q_ice[k] = FT(0)
                    if region_area > FT(0) # can be zero for a_en_remaining
                        if reweight_processes_for_grid
                            q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_en, q, p_c; reweight_extrema_only = reweight_extrema_only, region = region)
                            aux_en.q_liq[k] += q_liq * (region_area / aux_en.area[k])
                            aux_en.q_ice[k] += q_ice * (region_area / aux_en.area[k])
                            aux_en.cloud_fraction[k] += ((q_liq > 0) || (q_ice > 0)) ? (region_area / aux_en.area[k]) : 0 
                        else
                            if TD.has_condensate(thermo_params, ts)
                                q_liq = TD.liquid_specific_humidity(thermo_params, ts)
                                q_ice = TD.ice_specific_humidity(thermo_params, ts)
                                aux_en.q_liq[k] += q_liq * (region_area / aux_en.area[k])
                                aux_en.q_ice[k] += q_ice * (region_area / aux_en.area[k])
                                aux_en.cloud_fraction[k] += (region_area / aux_en.area[k])
                            end
                        end
                    end
                end
            else
                ts = ts_env[k]
                if reweight_processes_for_grid
                    q = TD.PhasePartition(thermo_params, ts) # we don't have q cached here but we do in updraft so calculating outside the fcn saves work
                    q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_en, q, p_c; reweight_extrema_only = reweight_extrema_only)
                    aux_en.q_liq[k] = q_liq
                    aux_en.q_ice[k] = q_ice
                    aux_en.cloud_fraction[k] = ((q_liq > 0) || (q_ice > 0)) ? 1 : 0 
                else
                    if TD.has_condensate(thermo_params, ts)
                        aux_en.cloud_fraction[k] = 1
                        aux_en.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
                        aux_en.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)
                    else
                        aux_en.cloud_fraction[k] = 0
                        aux_en.q_liq[k] = 0
                        aux_en.q_ice[k] = 0
                    end
                end
            end
        elseif moisture_model isa NonEquilibriumMoisture
            # aux_en.q_liq[k] and aux_en.q_ice[k] should already be set in update_aux()
            aux_en.cloud_fraction[k] = ((aux_en.q_liq[k] > 0) || (aux_en.q_ice[k] > 0)) ? 1 : 0
        end

        # update_sat_unsat
        # if TD.has_condensate(thermo_params, ts)
        if (aux_en.q_liq[k] > 0) || (aux_en.q_ice[k] > 0) # since we moved hte update above, don't need this
            aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
            if moisture_model isa NonEquilibriumMoisture
                # aux_en_sat.q_vap[k] = min(TD.q_vap_saturation(thermo_params, ts), TD.vapor_specific_humidity(thermo_params, ts)) # having condensate is no guarantee of saturation in noneq
                # aux_en_sat.q_vap[k] = TD.q_vap_saturation(thermo_params, ts) # having condensate is no guarantee of saturation in noneq
                aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts) # i think it is just the local qv, not sure. note we'll be stuck in the sat state until  all condensate is gone
                aux_en_sat.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts) # added
                aux_en_sat.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts) # added
            else
                aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
            end
            aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts)
            aux_en_sat.T[k] = TD.air_temperature(thermo_params, ts)
            aux_en_sat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)
        else
            aux_en.cloud_fraction[k] = 0
            aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
            aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts)
            aux_en_unsat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)

            aux_en_sat.T[k] = aux_en.T[k]
            if moisture_model isa NonEquilibriumMoisture # In noneq, not having condensate doesn't guarantee unsat
                # aux_en_sat.q_vap[k] = min(TD.q_vap_saturation(thermo_params, ts), TD.vapor_specific_humidity(thermo_params, ts)) # having condensate is no guarantee of saturation in noneq
                aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts) # i think it is just the local qv, not sure. note we'll be stuck in the sat state until  all condensate is gone
                aux_en_sat.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts) # added
                aux_en_sat.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts) # added
            else
                aux_en_sat.q_vap[k] = FT(0) #  i think this is wrong...
                # aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
            end
            aux_en_sat.q_tot[k] = aux_en.q_tot[k]
            aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
            aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]
        end




        # Storage  [final]
        # update_env_precip_tendencies
        # TODO: move ..._tendency_precip_formation to diagnostics
        aux_en.qt_tendency_precip_formation[k] = mph_precip.qt_tendency * aux_en.area[k]
        aux_en.θ_liq_ice_tendency_precip_formation[k] = mph_precip.θ_liq_ice_tendency * aux_en.area[k]
        if moisture_model isa NonEquilibriumMoisture
            aux_en.ql_tendency_precip_formation[k] = mph_precip.ql_tendency * aux_en.area[k]
            aux_en.qi_tendency_precip_formation[k] = mph_precip.qi_tendency * aux_en.area[k]
        end
        tendencies_pr.q_rai[k] += mph_precip.qr_tendency * aux_en.area[k]
        tendencies_pr.q_sno[k] += mph_precip.qs_tendency * aux_en.area[k]

        # store autoconversion and accretion for diagnostics (doens't mean much for Eq since liq/ice get set by sat adjust and T...)
        aux_en.ql_tendency_acnv[k] = mph_precip.ql_tendency_acnv * aux_en.area[k]
        aux_en.qi_tendency_acnv[k] = mph_precip.qi_tendency_acnv * aux_en.area[k]
        aux_en.qi_tendency_acnv_dep[k] = mph_precip.qi_tendency_acnv_dep * aux_en.area[k]
        aux_en.qi_tendency_acnv_dep_is[k] = mph_precip.qi_tendency_acnv_dep_is * aux_en.area[k]
        aux_en.qi_tendency_acnv_dep_above[k] = mph_precip.qi_tendency_acnv_dep_above * aux_en.area[k]
        aux_en.qi_tendency_acnv_agg_mix[k] = mph_precip.qi_tendency_acnv_agg_mix * aux_en.area[k] # this shoud have already been set but
        aux_en.qi_tendency_acnv_thresh[k] = mph_precip.qi_tendency_acnv_thresh * aux_en.area[k] # this shoud have already been set but
        aux_en.ql_tendency_accr_liq_rai[k] = mph_precip.ql_tendency_accr_liq_rai * aux_en.area[k]
        aux_en.ql_tendency_accr_liq_ice[k] = mph_precip.ql_tendency_accr_liq_ice * aux_en.area[k]
        aux_en.ql_tendency_accr_liq_sno[k] = mph_precip.ql_tendency_accr_liq_sno * aux_en.area[k]
        aux_en.qi_tendency_accr_ice_liq[k] = mph_precip.qi_tendency_accr_ice_liq * aux_en.area[k]
        aux_en.qi_tendency_accr_ice_rai[k] = mph_precip.qi_tendency_accr_ice_rai * aux_en.area[k]
        aux_en.qi_tendency_accr_ice_sno[k] = mph_precip.qi_tendency_accr_ice_sno * aux_en.area[k]
        #
        aux_tc.qs_tendency_accr_rai_sno[k] += mph_precip.qs_tendency_accr_rai_sno * aux_en.area[k] # we calculate in aux/en for the temperature dependence but store a combined output
    end



    return nothing
end


const QuadratureVarsType{FT} = @NamedTuple{
    qt′qt′::FT,
    qt_mean::FT,
    θl′θl′::FT,
    θl_mean::FT,
    θl′qt′::FT,
    subdomain_area::FT,
    q_rai::FT,
    q_sno::FT,
    ρ_c::FT,
    p_c::FT,
    precip_frac::FT,
    use_fallback_tendency_limiters::Bool,
    aux_en::CC.Fields.Field,
    aux_en_f::CC.Fields.Field,
    k::Cent{Int64},
    q_liq::FT,
    q_ice::FT,
    N_l::FT,
    N_i::FT,
    N_i_no_boost::FT,
    term_vel_liq::FT,
    term_vel_ice::FT,
    term_vel_rain::FT,
    term_vel_snow::FT,
    w::FT,
    tke::FT,
    ql_tendency_sedimentation::FT,
    ql_tendency_sedimentation_other::FT,
    qi_tendency_sedimentation::FT,
    qi_tendency_sedimentation_other::FT,
    dqvdt::FT,
    dTdt::FT,
    S_i::FT,
    N_INP::FT,
    τ_liq::FT,
    τ_ice::FT,
    dN_i_dz::FT,
    dqidz::FT,
    massflux::FT,
}

# doesnt really support NonEq yet
# function quad_loop(en_thermo::SGSQuadrature, grid, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, tendency_limiters, vars, param_set, Δt::Real, use_fallback_tendency_limiters::Bool)
function quad_loop(en_thermo::SGSQuadrature, grid::Grid, edmf::EDMFModel, moisture_model::AbstractMoistureModel, precip_model::AbstractPrecipitationModel, cloud_sedimentation_model::AbstractCloudSedimentationModel, rain_formation_model::AbstractRainFormationModel, snow_formation_model::AbstractSnowFormationModel, Δt::FT, param_set::APS, use_fallback_tendency_limiters::Bool, vars::QuadratureVarsType{FT}) where {FT} # vars is a named_tuple

    env_len = 8
    src_len = 8
    i_ql, i_qi, i_T, i_cf, i_qt_sat, i_qt_unsat, i_T_sat, i_T_unsat = 1:env_len
    i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH, i_Sqr, i_Sqs = 1:src_len


    raw_src_len = 8
    raw_sourcenames = (:SH_qt, :Sqt_H, :SH_H, :Sqt_qt, :Sqt, :SH, :Sqr, :Sqs)
    # i_raw_src_len = (;name => i for (i, name) in enumerate(raw_sourcenames))
    src_mph_neq_len = fieldcount(NoneqMoistureSources)  # [2] these would be extra additions to tendencies... Smph_neq `NoneqMoistureSources`
    # i_NoneqMoistureSources = (;name => i for (i, name) in enumerate(fieldnames(NoneqMoistureSources)))
    src_mph_neq_other_len = fieldcount(OtherMicrophysicsSources)  # [6] these would be extra additions to tendencies... Smph_neq_other `NoneqMoistureSources`
    # i_OtherMicrophysicsSources = (;name => i for (i, name) in enumerate(fieldnames(OtherMicrophysicsSources)))
    src_mph_precip_len = fieldcount(PrecipFormation) # [20] these would be extra additions to tendencies... Smph `PrecipitationFormation`
    # i_PrecipFormation = (;name => i for (i, name) in enumerate(fieldnames(PrecipFormation)))


    all_src_len = raw_src_len + src_mph_neq_len + src_mph_neq_other_len + src_mph_precip_len  # [36] these would be extra
    # combine to one long named tuple in order raw, Noneq, noneq_other, precip
    # all_names = (raw_sourcenames..., fieldnames(NoneqMoistureSources)..., fieldnames(OtherMicrophysicsSources)..., fieldnames(PrecipFormation)...) # This can have clobbers...
    all_names = (
        raw_sourcenames...,
        (Symbol(:NoneqMoistureSources_, name) for name in fieldnames(NoneqMoistureSources))...,
        (Symbol(:OtherMicrophysicsSources_, name) for name in fieldnames(OtherMicrophysicsSources))...,
        (Symbol(:PrecipFormation_, name) for name in fieldnames(PrecipFormation))...
    )
    i_src = NamedTuple{all_names, NTuple{all_src_len, Int}}(Tuple(1:length(all_names)))


    i_NoneqMoistureSources = NamedTuple{fieldnames(NoneqMoistureSources), NTuple{src_mph_neq_len, Int}}(Tuple(i_src[Symbol(:NoneqMoistureSources_, name)] for name in fieldnames(NoneqMoistureSources)))
    i_OtherMicrophysicsSources = NamedTuple{fieldnames(OtherMicrophysicsSources), NTuple{src_mph_neq_other_len, Int}}(Tuple(i_src[Symbol(:OtherMicrophysicsSources_, name)] for name in fieldnames(OtherMicrophysicsSources)))
    i_PrecipFormation = NamedTuple{fieldnames(PrecipFormation), NTuple{src_mph_precip_len, Int}}(Tuple(i_src[Symbol(:PrecipFormation_, name)] for name in fieldnames(PrecipFormation)))


    # additions to tendencies... Smph

    thermo_params = TCP.thermodynamics_params(param_set)
    quadrature_type = en_thermo.quadrature_type
    quad_order = quadrature_order(en_thermo)
    χ = en_thermo.a
    weights = en_thermo.w

    # qt - total water specific humidity
    # θl - liquid ice potential temperature
    # _mean and ′ - subdomain mean and (co)variances
    # q_rai, q_sno - grid mean precipitation
    UnPack.@unpack qt′qt′, qt_mean, θl′θl′, θl_mean, θl′qt′, subdomain_area, q_rai, q_sno, ρ_c, p_c, precip_frac = vars
    # (; q_liq, q_ice, N_i, term_vel_ice, term_vel_rain, term_vel_snow) = vars # property desturcturing exists in julia now...
    # FT = eltype(ρ_c)

    inner_env = SA.MVector{env_len, FT}(undef)
    outer_env = SA.MVector{env_len, FT}(undef)
    # inner_src = SA.MVector{src_len, FT}(undef)
    # outer_src = SA.MVector{src_len, FT}(undef)
    inner_src = SA.MVector{all_src_len, FT}(undef)
    outer_src = SA.MVector{all_src_len, FT}(undef)

    sqpi_inv = FT(1 / sqrt(π))
    sqrt2 = FT(sqrt(2))

    # Epsilon defined per typical variable fluctuation
    eps_q = qt_mean ≈ FT(0) ? eps(FT) : eps(FT) * qt_mean
    eps_θ = eps(FT)

    if quadrature_type isa LogNormalQuad
        # Lognormal parameters (ν, s) from mean and variance
        ν_q = log(qt_mean^2 / max(sqrt(qt_mean^2 + qt′qt′), eps_q))
        ν_θ = log(θl_mean^2 / sqrt(θl_mean^2 + θl′θl′))
        s_q = sqrt(log(qt′qt′ / max(qt_mean, eps_q)^2 + 1))
        s_θ = sqrt(log(θl′θl′ / θl_mean^2 + 1))

        # Enforce Cauchy-Schwarz inequality, numerically stable compute
        corr = θl′qt′ / max(sqrt(qt′qt′), eps_q)
        corr = max(min(corr / max(sqrt(θl′θl′), eps_θ), 1), -1)

        # Conditionals
        s2_θq = log(corr * sqrt(θl′θl′ * qt′qt′) / θl_mean / max(qt_mean, eps_q) + 1)
        s_c = sqrt(max(s_θ^2 - s2_θq^2 / max(s_q, eps_q)^2, 0))

    elseif quadrature_type isa GaussianQuad
        # limit σ_q to prevent negative qt_hat
        σ_q_lim = -qt_mean / (sqrt2 * χ[1])
        σ_q = min(sqrt(qt′qt′), σ_q_lim)
        σ_θ = sqrt(θl′θl′)

        # Enforce Cauchy-Schwarz inequality, numerically stable compute
        corr = θl′qt′ / max(σ_q, eps_q)
        corr = max(min(corr / max(σ_θ, eps_θ), 1), -1)

        # Conditionals
        σ_c = sqrt(max(1 - corr * corr, 0)) * σ_θ
    end

    # zero outer quadrature points
    @inbounds for idx in 1:env_len
        outer_env[idx] = 0
    end
    # @inbounds for idx in 1:src_len
    for idx in 1:all_src_len
        outer_src[idx] = 0
    end
    

    # calculate saturation exceses for partitinoning ql, qi
    if moisture_model isa NonEquilibriumMoisture
        (;q_liq, q_ice) = vars # unpack
        saturation_excesses_liq = SA.MArray{Tuple{quad_order, quad_order}, FT}(undef)
        saturation_excesses_ice = SA.MArray{Tuple{quad_order, quad_order}, FT}(undef)

        q_liq_ens = SA.MArray{Tuple{quad_order, quad_order}, FT}(undef)
        q_ice_ens = SA.MArray{Tuple{quad_order, quad_order}, FT}(undef)

        @inbounds for m_q in 1:quad_order
            if quadrature_type isa LogNormalQuad
                qt_hat = exp(ν_q + sqrt2 * s_q * χ[m_q])
                ν_c = ν_θ + s2_θq / max(s_q, eps_q)^2 * (log(qt_hat) - ν_q)
            elseif quadrature_type isa GaussianQuad
                qt_hat = qt_mean + sqrt2 * σ_q * χ[m_q]
                μ_c = θl_mean + sqrt2 * corr * σ_θ * χ[m_q]
            end

            # zero inner quadrature points
            inner_env .= 0
            inner_src .= 0

            # min_q_vap_sat_liq = FT(Inf)
            # min_q_vap_sat_ice = FT(Inf)
            for m_h in 1:quad_order
                if quadrature_type isa LogNormalQuad
                    h_hat = exp(ν_c + sqrt2 * s_c * χ[m_h])
                elseif quadrature_type isa GaussianQuad
                    h_hat = μ_c + sqrt2 * σ_c * χ[m_h]
                end      

                ts = thermo_state_pθq(param_set, p_c, h_hat, qt_hat, q_liq, q_ice)
                T = TD.air_temperature(thermo_params, ts)
                p = TD.air_pressure(thermo_params, ts)  
                ρ = TD.air_density(thermo_params, ts)
                q = TD.PhasePartition(thermo_params, ts)
                saturation_excesses_liq[m_q, m_h] = TD.vapor_specific_humidity(q) -  TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
                saturation_excesses_ice[m_q, m_h] = TD.vapor_specific_humidity(q) -  TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())

                # saturation_excesses_liq[m_q, m_h] = TD.vapor_specific_humidity(thermo_params, ts)
                # saturation_excesses_ice[m_q, m_h] = TD.vapor_specific_humidity(thermo_params, ts)
                # min_q_vap_sat_liq = min(min_q_vap_sat_liq, TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()))
                # min_q_vap_sat_ice = min(min_q_vap_sat_ice, TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
            end
        end

        # the mean value is 
        # ts = thermo_state_pθq(param_set, p_c, θl_mean, qt_mean, q_liq, q_ice)
        # T = TD.air_temperature(thermo_params, ts)
        # p = TD.air_pressure(thermo_params, ts)  
        # q = TD.PhasePartition(thermo_params, ts)
        # ρ = TD.air_density(thermo_params, ts)
        # mean_q_vap_sat_liq = TD.saturation_excess(thermo_params, T, ρ, TD.Liquid(), q)
        # mean_q_vap_sat_ice = TD.saturation_excess(thermo_params, T, ρ, TD.Ice(), q)

        # we want to reduce if below the average saturation excess, increase if above, and conserve the total, with the weights...


        #=
            We estimate var(q)/q ≈ var(S)/S. but S can be positive or negative and even have mean 0 so you'd need to shift S so all values are positive
            We can for simplicity assume
                va(q) = c var(S)
            cov(q,S) = c √(var(q)*var(S)) = √(c*var(S)*var(S)) = c*var(S)
            As for what c is, maybe go with 1 for liq since it's fast and 1e-3 for ice since it's slow, but these could be calibratable.
        =#
        W = (weights .* sqpi_inv) .* (weights .* sqpi_inv)'
        Wtot = sum(W)

        S̄_liq = sum(W .* saturation_excesses_liq) / Wtot      # Mean of S (μ_S)
        δS_liq = saturation_excesses_liq .- S̄_liq
        var_S_liq = sum(W .* δS_liq .^ 2) / Wtot # Variance of S (σ_S²)

        S̄_ice = sum(W .* saturation_excesses_ice) / Wtot      # Mean of S (μ_S)
        δS_ice = saturation_excesses_ice .- S̄_ice
        var_S_ice = sum(W .* δS_ice .^ 2) / Wtot # Variance of S (σ_S²)

        if (q_liq > 2eps(FT)) && (var_S_liq > eps(FT))
            q′S′ = FT(1) * var_S_liq
            q_liq_ens .= get_qs_from_saturation_excesses(saturation_excesses_liq, weights .* sqpi_inv, q_liq; q′S′ = q′S′, q′q′ = nothing) # for liquid we'll go 1:1 cause it's fast, though supersat struggles to exceed 1. so it might be skewed idk...
            # q_liq_ens .= q_liq
        else
            q_liq_ens .= q_liq
        end
        if (q_ice > 2eps(FT)) && (var_S_ice > eps(FT))
            q′S′ = FT(1e-3) * var_S_ice
            q_ice_ens .= get_qs_from_saturation_excesses(saturation_excesses_ice, weights .* sqpi_inv, q_ice; q′S′ = q′S′, q′q′ = nothing) # for ice we'll go 1:10 because it's slow but vapor excess can be is much higher relative to the ice amount..., and then subimation is slow....
            # q_ice_ens .= q_ice
        else
            q_ice_ens .= q_ice
        end
       

        # total_saturation_excess_liq = sum(saturation_excesses_liq)
        # total_saturation_excess_ice = sum(saturation_excesses_ice)
    end 



    @inbounds for m_q in 1:quad_order
        if quadrature_type isa LogNormalQuad
            qt_hat = exp(ν_q + sqrt2 * s_q * χ[m_q])
            ν_c = ν_θ + s2_θq / max(s_q, eps_q)^2 * (log(qt_hat) - ν_q)
        elseif quadrature_type isa GaussianQuad
            qt_hat = qt_mean + sqrt2 * σ_q * χ[m_q]
            μ_c = θl_mean + sqrt2 * corr * σ_θ * χ[m_q]
        end

        # zero inner quadrature points
        inner_env .= 0
        inner_src .= 0

        for m_h in 1:quad_order

            if quadrature_type isa LogNormalQuad
                h_hat = exp(ν_c + sqrt2 * s_c * χ[m_h])
            elseif quadrature_type isa GaussianQuad
                h_hat = μ_c + sqrt2 * σ_c * χ[m_h]
            end



            

            # We want to try assigning q_l and q_i based on the fraction of total vapor above supersaturation we have. `TD.saturation_excess()`


            if moisture_model isa NonEquilibriumMoisture # Here we do noneqmoisture sources. We assume the same ql, qi everywhere w/ no quadrature though...

                # assume q_liq and q_en do not change and have no variance... ( do not participate in quadrature )
                (; use_fallback_tendency_limiters, aux_en, aux_en_f, k, q_liq, q_ice, N_l, N_i, N_i_no_boost, term_vel_liq, term_vel_ice, term_vel_rain, term_vel_snow, w, tke, ql_tendency_sedimentation, ql_tendency_sedimentation_other, qi_tendency_sedimentation, qi_tendency_sedimentation_other, dqvdt, dTdt, S_i, N_INP, τ_liq, τ_ice, dN_i_dz, dqidz, massflux) = vars # unpack

                ts = thermo_state_pθq(param_set, p_c, h_hat, qt_hat, q_liq, q_ice)
                T = TD.air_temperature(thermo_params, ts)
                p = TD.air_pressure(thermo_params, ts)
                ρ = TD.air_density(thermo_params, ts)
                q_vap_sat_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
                q_vap_sat_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
                # q_liq_en = q_liq
                # q_ice_en = q_ice

                # distribute by saturation excess over minimum
                q_liq_en = q_liq_ens[m_q, m_h]
                q_ice_en = q_ice_ens[m_q, m_h]

                (; mph_neq, mph_neq_other, mph_precip) = microphysics_helper(grid, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, Δt, param_set, use_fallback_tendency_limiters, aux_en, aux_en_f, k,
                    subdomain_area, h_hat, qt_hat, q_liq_en, q_ice_en, q_rai, q_sno, N_l, N_i, N_i_no_boost, term_vel_liq, term_vel_ice, term_vel_rain, term_vel_snow, ρ_c, p_c, T, ts, precip_frac,
                    ql_tendency_sedimentation, ql_tendency_sedimentation_other, qi_tendency_sedimentation, qi_tendency_sedimentation_other;
                    q_vap_sat_liq = q_vap_sat_liq, q_vap_sat_ice = q_vap_sat_ice, dqvdt = dqvdt, dTdt = dTdt, w=w, tke = tke, S_i = S_i, τ_liq = τ_liq, τ_ice = τ_ice, dN_i_dz = dN_i_dz, dqidz = dqidz, N_INP = N_INP, massflux = massflux)
            else
                (; use_fallback_tendency_limiters, aux_en, aux_en_f, k, N_l, N_i, N_i_no_boost, term_vel_liq, term_vel_ice, term_vel_rain, term_vel_snow, w, tke, ql_tendency_sedimentation, ql_tendency_sedimentation_other, qi_tendency_sedimentation, qi_tendency_sedimentation_other) = vars # unpack

                ts = thermo_state_pθq(param_set, p_c, h_hat, qt_hat)
                q_liq_en = TD.liquid_specific_humidity(thermo_params, ts)
                q_ice_en = TD.ice_specific_humidity(thermo_params, ts)
                T = TD.air_temperature(thermo_params, ts)
                p = TD.air_pressure(thermo_params, ts)

                (; mph_neq, mph_neq_other, mph_precip) = microphysics_helper(grid, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, Δt, param_set, use_fallback_tendency_limiters, aux_en, aux_en_f, k,
                    subdomain_area, h_hat, qt_hat, q_liq_en, q_ice_en, q_rai, q_sno, N_l, N_i, N_i_no_boost, term_vel_liq, term_vel_ice, term_vel_rain, term_vel_snow, ρ_c, p_c, T, ts, precip_frac,
                    ql_tendency_sedimentation, ql_tendency_sedimentation_other, qi_tendency_sedimentation, qi_tendency_sedimentation_other)
            end

            # environmental variables
            inner_env[i_ql] += q_liq_en * weights[m_h] * sqpi_inv
            inner_env[i_qi] += q_ice_en * weights[m_h] * sqpi_inv
            inner_env[i_T] += T * weights[m_h] * sqpi_inv
            # cloudy/dry categories for buoyancy in TKE
            if TD.has_condensate(q_liq_en + q_ice_en)
                if moisture_model isa EquilibriumMoisture
                    inner_env[i_cf] += weights[m_h] * sqpi_inv
                    inner_env[i_qt_sat] += qt_hat * weights[m_h] * sqpi_inv
                    inner_env[i_T_sat] += T * weights[m_h] * sqpi_inv
                elseif moisture_model isa NonEquilibriumMoisture # for NonEq, having condensate or not doesn't map to sat/unsat
                    inner_env[i_cf] += weights[m_h] * sqpi_inv
                    inner_env[i_qt_sat] += qt_hat * weights[m_h] * sqpi_inv
                    inner_env[i_T_sat] += T * weights[m_h] * sqpi_inv   
                else
                    error("moisture_model not recognized")
                end
            else
                if moisture_model isa EquilibriumMoisture
                    inner_env[i_qt_unsat] += qt_hat * weights[m_h] * sqpi_inv
                    inner_env[i_T_unsat] += T * weights[m_h] * sqpi_inv
                elseif moisture_model isa NonEquilibriumMoisture # for NonEq, having condensate or not doesn't map to sat/unsat
                    inner_env[i_qt_unsat] += qt_hat * weights[m_h] * sqpi_inv
                    inner_env[i_T_unsat] += T * weights[m_h] * sqpi_inv
                else
                    error("moisture_model not recognized")
                end
            end



            # My additions [[ none of these should change qt or θ_liq_ice_tendency below ]]
            for (_, name) in pairs(fieldnames(NoneqMoistureSources))
                # inner_src[i_src[Symbol(:NoneqMoistureSources_, name)]] += getfield(mph_neq, name) * weights[m_h] * sqpi_inv
                inner_src[i_NoneqMoistureSources[name]] += getfield(mph_neq, name) * weights[m_h] * sqpi_inv
            end
            for (_, name) in pairs(fieldnames(OtherMicrophysicsSources))
                # inner_src[i_src[Symbol(:OtherMicrophysicsSources_, name)]] += getfield(mph_neq_other, name) * weights[m_h] * sqpi_inv
                inner_src[i_OtherMicrophysicsSources[name]] += getfield(mph_neq_other, name) * weights[m_h] * sqpi_inv
            end
            for (_, name) in pairs(fieldnames(PrecipFormation))
                # inner_src[i_src[Symbol(:PrecipFormation_, name)]] += getfield(mph_precip, name) * weights[m_h] * sqpi_inv
                inner_src[i_PrecipFormation[name]] += getfield(mph_precip, name) * weights[m_h] * sqpi_inv
            end

            # products for variance and covariance source terms [ only mph_precip (and sedimentation) impact qt and θ_li]
            inner_src[i_Sqt] += mph_precip.qt_tendency * weights[m_h] * sqpi_inv
            inner_src[i_Sqr] += mph_precip.qr_tendency * weights[m_h] * sqpi_inv
            inner_src[i_Sqs] += mph_precip.qs_tendency * weights[m_h] * sqpi_inv
            inner_src[i_SH] += mph_precip.θ_liq_ice_tendency * weights[m_h] * sqpi_inv
            # covariance terms
            inner_src[i_Sqt_H] += mph_precip.qt_tendency * h_hat * weights[m_h] * sqpi_inv
            inner_src[i_Sqt_qt] += mph_precip.qt_tendency * qt_hat * weights[m_h] * sqpi_inv
            inner_src[i_SH_H] += mph_precip.θ_liq_ice_tendency * h_hat * weights[m_h] * sqpi_inv
            inner_src[i_SH_qt] += mph_precip.θ_liq_ice_tendency * qt_hat * weights[m_h] * sqpi_inv
        end

        for idx in 1:env_len
            outer_env[idx] += inner_env[idx] * weights[m_q] * sqpi_inv
        end
        # for idx in 1:src_len
        for idx in 1:all_src_len
            outer_src[idx] += inner_src[idx] * weights[m_q] * sqpi_inv
        end
    end

    outer_src_nt_Noneq = NamedTuple{keys(i_NoneqMoistureSources), NTuple{src_mph_neq_len, FT}}(map(i -> outer_src[i], values(i_NoneqMoistureSources)))
    outer_src_nt_Noneq_other = NamedTuple{keys(i_OtherMicrophysicsSources), NTuple{src_mph_neq_other_len, FT}}(map(i -> outer_src[i], values(i_OtherMicrophysicsSources)))
    outer_src_nt_Precip = NamedTuple{keys(i_PrecipFormation), NTuple{src_mph_precip_len, FT}}(map(i -> outer_src[i], values(i_PrecipFormation)))

    outer_src_nt = (;
        SH_qt = outer_src[i_SH_qt],
        Sqt_H = outer_src[i_Sqt_H],
        SH_H = outer_src[i_SH_H],
        Sqt_qt = outer_src[i_Sqt_qt],
        Sqt = outer_src[i_Sqt],
        SH = outer_src[i_SH],
        Sqr = outer_src[i_Sqr],
        Sqs = outer_src[i_Sqs],
        #
        NoneqMoistureSources = outer_src_nt_Noneq,
        OtherMicrophysicsSources = outer_src_nt_Noneq_other,
        PrecipFormation = outer_src_nt_Precip,
    )

    outer_env_nt = (;
        ql = outer_env[i_ql],
        qi = outer_env[i_qi],
        T = outer_env[i_T],
        cf = outer_env[i_cf],
        qt_sat = outer_env[i_qt_sat],
        qt_unsat = outer_env[i_qt_unsat],
        T_sat = outer_env[i_T_sat],
        T_unsat = outer_env[i_T_unsat],
    )
    return outer_env_nt, outer_src_nt
end




function microphysics!(
    en_thermo::SGSQuadrature,
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    moisture_model::AbstractMoistureModel,
    precip_model::AbstractPrecipitationModel,
    cloud_sedimentation_model::AbstractCloudSedimentationModel,
    rain_formation_model::AbstractRainFormationModel,
    snow_formation_model::AbstractSnowFormationModel,
    Δt::Real,
    param_set::APS,
    use_fallback_tendency_limiters::Bool
)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    tendencies_pr = center_tendencies_precipitation(state)
    aux_en = center_aux_environment(state)
    aux_en_2m = center_aux_environment_2m(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    aux_en_f = face_aux_environment(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
    prog_pr = center_prog_precipitation(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    ts_env = center_aux_environment(state).ts
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_en_sat = aux_en.sat
    aux_en_unsat = aux_en.unsat
    precip_fraction = compute_precip_fraction(edmf, state)

    # TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)

    epsilon = FT(1e-14) # eps(float)

     if (moisture_model isa NonEquilibriumMoisture) || (edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean) || ((edmf.area_partition_model isa CoreCloakAreaPartitionModel) && edmf.area_partition_model.apply_cloak_to_condensate_formation)
        if (edmf.area_partition_model isa CoreCloakAreaPartitionModel)
            if edmf.area_partition_model.confine_all_downdraft_to_cloak
                w::CC.Fields.Field = similar(aux_en.T)
                @. w = FT(0) # we confine the downdraft (which will be drier) to the cloak region, so we can set w = 0 in the environment
            else # not confined to cloak, but also we do need some to be converted to updraft, maybe still go with 0. That biases us towards respecting supersat generation in the cloak updraft...
                w = similar(aux_en.T)
                @. w = FT(0)
            end
        else
            w = Ic.(aux_en_f.w)
            # w = aux_en.w # gonna deprecate this I think for RAM savings
        end
    else
        # don't think we need w for EquilibriumMoisture... we're using w on faces rn...
    end

    # ======================================================================================================================================================================================== #
    # Calculate Sedimentation
    if edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean  

        mph_sed_ice, mph_sed_ice_other = calculate_sedimentation_sources(
            param_set,
            ρ_c,
            aux_en.q_ice,
            aux_en.term_vel_ice,
            aux_en_f.w, # could just be 0 bc we don't use it w/ use_relative_w = false, bc advection is handled elsewhere and we argue that diffusion is maybe ok since we don't model the size distribution
            aux_en.area,
            grid,
            differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
            grid_mean = false,
            use_relative_w = false,
        ) # should this be a grid mean tendency? 

        mph_sed_liq, mph_sed_liq_other = calculate_sedimentation_sources(
            param_set,
            ρ_c,
            aux_en.q_liq,
            aux_en.term_vel_liq,
            aux_en_f.w, # could just be 0 bc we don't use it w/ use_relative_w = false, bc advection is handled elsewhere and we argue that diffusion is maybe ok since we don't model the size distribution
            aux_en.area,
            grid,
            differencing_scheme = edmf.cloud_sedimentation_model.sedimentation_differencing_scheme,
            grid_mean = false,
            use_relative_w = false,
        ) # should this be a grid mean tendency? 

        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)
    end

    # ======================================================================================================================================================================================== #
    reweight_processes_for_grid::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_processes_for_grid, false)
    reweight_extrema_only::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false)
    # ======================================================================================================================================================================================== #


    # initialize the quadrature points and their labels

    @inbounds for k in real_center_indices(grid)
         # === sed ======================== [ do before precip_formation so we can know sed-driven aggregatory autoconversion, and it doesnt depend on anything below. ]
        if edmf.cloud_sedimentation_model isa CloudSedimentationModel && !edmf.cloud_sedimentation_model.grid_mean        
            ql_tendency_sedimentation = mph_sed_liq[k].q_tendency
            qi_tendency_sedimentation = mph_sed_ice[k].q_tendency
            qt_tendency_sedimentation = ql_tendency_sedimentation + qi_tendency_sedimentation
            aux_en.ql_tendency_sedimentation[k] += ql_tendency_sedimentation # these get added to gm in compute_gm
            aux_en.qi_tendency_sedimentation[k] += qi_tendency_sedimentation
            aux_en.qt_tendency_sedimentation[k] += qt_tendency_sedimentation # used in dycore.jl (= not += , cause this doesnt seem to get reset every iteration?) (fixed in update_aux)

            if !isfinite(qi_tendency_sedimentation)
                @error "qi_tendency_sedimentation is NaN at k = $k; aux_en.qi_tendency_sedimentation[k] = $(aux_en.qi_tendency_sedimentation[k]); aux_en.qi_tendency_sedimentation[k] = $(aux_en.qi_tendency_sedimentation[k])"
            end

            Π_m = TD.exner(thermo_params, ts_env[k])
            c_pm = TD.cp_m(thermo_params, ts_env[k])
            θ_liq_ice_tendency_sedimentation =
                1 / Π_m / c_pm * (L_v0 * ql_tendency_sedimentation + L_s0 * qi_tendency_sedimentation)
            aux_en.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # adapted from microphysics_coupling.jl: precipitation_formation() | (= not += caue these don't seem to get reset every iteration?)

            # sedimentation loss into updraft (should this be allowed lol)
            # How do we know when these get set to 0 and can be used again?
            if aux_en.area[k] < 1 # we have some updrafts
                N_up = n_updrafts(edmf)
                @inbounds for i in 1:N_up
                    ql_tendency_sedimentation_other = mph_sed_liq_other[k].q_tendency .* (aux_up[i].area[k] ./ aux_bulk.area[k])
                    qi_tendency_sedimentation_other = mph_sed_ice_other[k].q_tendency .* (aux_up[i].area[k] ./ aux_bulk.area[k])
                    qt_tendency_sedimentation_other = ql_tendency_sedimentation_other + qi_tendency_sedimentation_other
                    θ_liq_ice_tendency_sedimentation_other =
                        1 / Π_m / c_pm *
                        (L_v0 * ql_tendency_sedimentation_other + L_s0 * qi_tendency_sedimentation_other)
                    aux_up[i].ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                    aux_up[i].qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                    aux_up[i].qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                    aux_up[i].θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other
                    aux_bulk.ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                    aux_bulk.qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                    aux_bulk.qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                    aux_bulk.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other

                    # store contributions to other
                    aux_up[i].ql_tendency_sedimentation_other[k] += ql_tendency_sedimentation_other # for storage
                    aux_up[i].qi_tendency_sedimentation_other[k] += qi_tendency_sedimentation_other # for storage
                    aux_up[i].qt_tendency_sedimentation_other[k] += qt_tendency_sedimentation_other # for storage
                    aux_up[i].θ_liq_ice_tendency_sedimentation_other[k] += θ_liq_ice_tendency_sedimentation_other # for storage

                    aux_bulk.ql_tendency_sedimentation_other[k] += ql_tendency_sedimentation_other # for storage
                    aux_bulk.qi_tendency_sedimentation_other[k] += qi_tendency_sedimentation_other # for storage
                    aux_bulk.qt_tendency_sedimentation_other[k] += qt_tendency_sedimentation_other # for storage
                    aux_bulk.θ_liq_ice_tendency_sedimentation_other[k] += θ_liq_ice_tendency_sedimentation_other # for storage
                end
            end
        end
        # ============================================== #
        ts = ts_env[k]


        if (
            aux_en.QTvar[k] > epsilon &&
            aux_en.Hvar[k] > epsilon &&
            abs(aux_en.HQTcov[k]) > epsilon &&
            aux_en.q_tot[k] > epsilon &&
            sqrt(aux_en.QTvar[k]) < aux_en.q_tot[k]
        )
            #= 
            =#
            q = TD.PhasePartition(thermo_params, ts)
            ρ = TD.air_density(thermo_params, ts)

            vars::QuadratureVarsType{FT} = (;
                qt′qt′ = aux_en.QTvar[k],
                qt_mean = aux_en.q_tot[k],
                θl′θl′ = aux_en.Hvar[k],
                θl_mean = aux_en.θ_liq_ice[k],
                θl′qt′ = aux_en.HQTcov[k],
                subdomain_area = aux_en.area[k],
                q_rai = prog_pr.q_rai[k],
                q_sno = prog_pr.q_sno[k],
                ρ_c = ρ_c[k],
                p_c = p_c[k],
                precip_frac = precip_fraction,
                #
                # My Additions
                #
                use_fallback_tendency_limiters = use_fallback_tendency_limiters, 
                #
                aux_en = aux_en, # pass aux_en for access to ts etc
                aux_en_f = aux_en_f,
                k = k,
                #
                q_liq = aux_en.q_liq[k],
                q_ice = aux_en.q_ice[k],
                #
                N_l = aux_en.N_l[k], 
                N_i = aux_en.N_i[k],
                N_i_no_boost = aux_en.N_i_no_boost[k],
                term_vel_liq = aux_en.term_vel_liq[k],
                term_vel_ice = aux_en.term_vel_ice[k],
                term_vel_rain = aux_tc.term_vel_rain[k],
                term_vel_snow = aux_tc.term_vel_snow[k],
                #
                w = (moisture_model isa NonEquilibriumMoisture) ? w[k] : FT(NaN), # only needed for NonEq
                #
                tke = aux_en.tke[k],
                ql_tendency_sedimentation = (aux_en.area[k] > FT(0)) ? aux_en.ql_tendency_sedimentation[k] / aux_en.area[k] : FT(0),
                ql_tendency_sedimentation_other = (aux_en.area[k] > FT(0)) ? aux_en.ql_tendency_sedimentation_other[k] / aux_en.area[k] : FT(0),
                qi_tendency_sedimentation = (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation[k] / aux_en.area[k] : FT(0),
                qi_tendency_sedimentation_other = (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation_other[k] / aux_en.area[k] : FT(0),
                dqvdt = aux_en.dqvdt[k],
                dTdt = aux_en.dTdt[k],
                #
                S_i = (moisture_model isa NonEquilibriumMoisture) ? TD.supersaturation(thermo_params, q, ρ, aux_en.T[k], TD.Ice()) : FT(NaN), # only needed for NonEq
                N_INP = (moisture_model isa NonEquilibriumMoisture) ? get_INP_concentration(param_set, moisture_model.scheme, q, aux_en.T[k], ρ, w[k]) : FT(NaN), # only needed for NonEq
                τ_liq = aux_en.τ_liq[k], # only needed for NonEq
                τ_ice = aux_en.τ_ice[k], # only needed for NonEq
                dN_i_dz = aux_en.dN_i_dz[k], # only needed for NonEq
                dqidz = aux_en.dqidz[k], # only needed for NonEq
                massflux = aux_tc.massflux[k],
            )
            
            # outer_env, outer_src = quad_loop(en_thermo, grid, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, edmf.tendency_limiters, vars, param_set, Δt, use_fallback_tendency_limiters)
            # function quad_loop(en_thermo::SGSQuadrature, grid::Grid, edmf::EDMFModel, moisture_model::AbstractMoistureModel, precip_model::AbstractPrecipitationModel, cloud_sedimentation_model::AbstractCloudSedimentationModel, rain_formation_model::AbstractRainFormationModel, snow_formation_model::AbstractSnowFormationModel, Δt::FT, param_set::APS, use_fallback_tendency_limiters::Bool, aux_en::CC.Fields.Field, aux_en_f::CC.Fields.Field, vars) where {FT} # vars is a named_tuple
            outer_env, outer_src = quad_loop(en_thermo, grid, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, Δt, param_set, use_fallback_tendency_limiters, vars)
            
            # update environmental cloudy/dry variables for buoyancy in TKE
            # update_env_precip_tendencies
            qt_tendency = outer_src.Sqt # should be the same as outer_src.PrecipFormation.qt_tendency
            θ_liq_ice_tendency = outer_src.SH # should be the same as outer_src.PrecipFormation.θ_liq_ice_tendency
            qr_tendency = outer_src.Sqr # should be the same as outer_src.PrecipFormation.qr_tendency
            qs_tendency = outer_src.Sqs # should be the same as outer_src.PrecipFormation.qs_tendency
            # TODO: move ..._tendency_precip_formation to diagnostics
            aux_en.qt_tendency_precip_formation[k] = qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = θ_liq_ice_tendency * aux_en.area[k]

            tendencies_pr.q_rai[k] += qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += qs_tendency * aux_en.area[k]


            if TD.has_condensate(outer_env.ql + outer_env.qi)
                aux_en.cloud_fraction[k] = outer_env.cf
                if moisture_model isa EquilibriumMoisture
                    aux_en.q_liq[k] = outer_env.ql
                    aux_en.q_ice[k] = outer_env.qi
                end # otherwise, let update_aux() take care of it... We will have set tendencies but not modified q_liq and q_ice

            else
                aux_en.cloud_fraction[k] = 0.0
                if moisture_model isa EquilibriumMoisture
                    aux_en.q_liq[k] = 0.0
                    aux_en.q_ice[k] = 0.0

                    # ---  if noneq, things could be happening, e.g. dep and snow formation, or accretion rain-snow  --- #
                    # aux_en.qt_tendency_precip_formation[k] = 0.0
                    # aux_en.θ_liq_ice_tendency_precip_formation[k] = 0.0
                    # tendencies_pr.q_rai[k] = 0.0
                    # tendencies_pr.q_sno[k] = 0.0
                    # Even in eq, why should these be 0? We could still have e.g. rain snow accretion
                end

                aux_en.Hvar_rain_dt[k] = 0.0
                aux_en.QTvar_rain_dt[k] = 0.0
                aux_en.HQTcov_rain_dt[k] = 0.0

            end

            if aux_en.cloud_fraction[k] < 1
                aux_en_unsat.q_tot[k] = outer_env.qt_unsat / (1 - aux_en.cloud_fraction[k])
                T_unsat = outer_env.T_unsat / (1 - aux_en.cloud_fraction[k])
                ts_unsat = TD.PhaseEquil_pTq(thermo_params, p_c[k], T_unsat, aux_en_unsat.q_tot[k])
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts_unsat)
            else
                aux_en_unsat.q_tot[k] = 0
                aux_en_unsat.θ_dry[k] = 0
            end

            if aux_en.cloud_fraction[k] > 0
                aux_en_sat.T[k] = outer_env.T_sat / aux_en.cloud_fraction[k]
                aux_en_sat.q_tot[k] = outer_env.qt_sat / aux_en.cloud_fraction[k]
                aux_en_sat.q_vap[k] = (outer_env.qt_sat - outer_env.ql - outer_env.qi) / aux_en.cloud_fraction[k]
                ts_sat = TD.PhaseEquil_pTq(thermo_params, p_c[k], aux_en_sat.T[k], aux_en_sat.q_tot[k])
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts_sat)
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts_sat)
            else
                if moisture_model isa EquilibriumMoisture
                    aux_en_sat.T[k] = aux_en.T[k]
                    aux_en_sat.q_vap[k] = FT(0)
                    aux_en_sat.q_tot[k] = aux_en.q_tot[k]
                    aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
                    aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]
                else
                    ts_sat = TD.PhaseEquil_pTq(thermo_params, p_c[k], aux_en.T[k], aux_en.q_tot[k]) # recompute a saturated state at the grid mean qt
                    aux_en_sat.T[k] = TD.air_temperature(thermo_params, ts_sat)
                    aux_en_sat.q_vap[k] = FT(0)
                    aux_en_sat.q_tot[k] = aux_en.q_tot[k]
                    aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts_sat)
                    aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts_sat)
                    aux_en_sat.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts_sat)
                    aux_en_sat.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts_sat)
                end
            end

            # update var/covar rain sources
            aux_en.Hvar_rain_dt[k] = outer_src.SH_H - outer_src.SH * aux_en.θ_liq_ice[k]
            aux_en.QTvar_rain_dt[k] = outer_src.Sqt_qt - outer_src.Sqt * aux_en.q_tot[k]
            aux_en.HQTcov_rain_dt[k] = outer_src.SH_qt - outer_src.SH * aux_en.q_tot[k] + outer_src.Sqt_H - outer_src.Sqt * aux_en.θ_liq_ice[k]

            # ---------------------------------------------------------------------------------- #
            # Storage [w/ quadrature]
            # derived sources
            if moisture_model isa NonEquilibriumMoisture
                # NoneqMoistureSources
                aux_en.ql_tendency_cond_evap[k] = outer_src.NoneqMoistureSources.ql_tendency * aux_en.area[k]
                aux_en.qi_tendency_sub_dep[k] = outer_src.NoneqMoistureSources.qi_tendency * aux_en.area[k]
                # OtherMicrophysicsSources
                aux_en.qi_tendency_hom_frz[k] = outer_src.OtherMicrophysicsSources.qi_tendency_homogeneous_freezing * aux_en.area[k]
                aux_en.qi_tendency_het_frz[k] = outer_src.OtherMicrophysicsSources.qi_tendency_heterogeneous_freezing * aux_en.area[k]
                aux_en.qi_tendency_het_nuc[k] = outer_src.OtherMicrophysicsSources.qi_tendency_heterogeneous_icenuc * aux_en.area[k]
                aux_en.qi_tendency_melt[k] = outer_src.OtherMicrophysicsSources.qi_tendency_melting * aux_en.area[k]
                # NoneqMoistureSources + OtherMicrophysicsSources
                aux_en.ql_tendency_noneq[k] = outer_src.NoneqMoistureSources.ql_tendency * aux_en.area[k] + outer_src.OtherMicrophysicsSources.ql_tendency * aux_en.area[k]
                aux_en.qi_tendency_noneq[k] = outer_src.NoneqMoistureSources.qi_tendency * aux_en.area[k] + outer_src.OtherMicrophysicsSources.qi_tendency * aux_en.area[k]
                # PrecipFormation
                aux_en.ql_tendency_precip_formation[k] = outer_src.PrecipFormation.ql_tendency * aux_en.area[k]
                aux_en.qi_tendency_precip_formation[k] = outer_src.PrecipFormation.qi_tendency * aux_en.area[k]
            end
            # PrecipFormation
            # aux_en.qt_tendency_precip_formation[k] :: Already Set Above
            # tendencies_pr.q_rai[k] :: Already Set Above
            # tendencies_pr.q_sno[k] :: Already Set Above
            aux_en.ql_tendency_acnv[k] = outer_src.PrecipFormation.ql_tendency_acnv * aux_en.area[k]
            aux_en.qi_tendency_acnv[k] = outer_src.PrecipFormation.qi_tendency_acnv * aux_en.area[k]
            aux_en.qi_tendency_acnv_dep[k] = outer_src.PrecipFormation.qi_tendency_acnv_dep * aux_en.area[k]
            aux_en.qi_tendency_acnv_dep_is[k] = outer_src.PrecipFormation.qi_tendency_acnv_dep_is * aux_en.area[k]
            aux_en.qi_tendency_acnv_dep_above[k] = outer_src.PrecipFormation.qi_tendency_acnv_dep_above * aux_en.area[k]
            aux_en.qi_tendency_acnv_agg_mix[k] = outer_src.PrecipFormation.qi_tendency_acnv_agg_mix * aux_en.area[k] # this shoud have already been set but
            aux_en.qi_tendency_acnv_thresh[k] = outer_src.PrecipFormation.qi_tendency_acnv_thresh * aux_en.area[k] # this shoud have already been set but
            aux_en.ql_tendency_accr_liq_rai[k] = outer_src.PrecipFormation.ql_tendency_accr_liq_rai * aux_en.area[k]
            aux_en.ql_tendency_accr_liq_ice[k] = outer_src.PrecipFormation.ql_tendency_accr_liq_ice * aux_en.area[k]
            aux_en.ql_tendency_accr_liq_sno[k] = outer_src.PrecipFormation.ql_tendency_accr_liq_sno * aux_en.area[k]
            aux_en.qi_tendency_accr_ice_liq[k] = outer_src.PrecipFormation.qi_tendency_accr_ice_liq * aux_en.area[k]
            aux_en.qi_tendency_accr_ice_rai[k] = outer_src.PrecipFormation.qi_tendency_accr_ice_rai * aux_en.area[k]
            aux_en.qi_tendency_accr_ice_sno[k] = outer_src.PrecipFormation.qi_tendency_accr_ice_sno * aux_en.area[k]
            #
            aux_tc.qs_tendency_accr_rai_sno[k] += outer_src.PrecipFormation.qs_tendency_accr_rai_sno * aux_en.area[k] # we calculate in aux/en for the temperature dependence but store a combined output
            # ---------------------------------------------------------------------------------- #

        else 
            # if variance and covariance are zero do the same as in SA_mean
            # microphysics!(SGSMean(), grid, state, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, Δt, param_set, use_fallback_tendency_limiters) # this won't work bc microphysics() goes over all k... rippity lipstick. # # whoops, this is just 1000% wrong... we need to do what microphysics() does, but only for this k...  maybe if we broke out the part that writes as a helper fcn of k we could just call that and then have the sedimentation outside... idk.

            if moisture_model isa NonEquilibriumMoisture # Here we do noneqmoisture sources. We assume the same ql, qi everywhere w/ no quadrature though...
                # assume q_liq and q_en do not change and have no variance... ( do not participate in quadrature )
                ρ = TD.air_density(thermo_params, ts)
                q = TD.PhasePartition(thermo_params, ts_env[k]) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
                S_i = TD.supersaturation(thermo_params, q, ρ, aux_en.T[k], TD.Ice())
                N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, aux_en.T[k], ρ, w[k])

                (; mph_neq, mph_neq_other, mph_precip) = microphysics_helper(grid, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, Δt, param_set, use_fallback_tendency_limiters, aux_en, aux_en_f, k,
                    aux_en.area[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k], prog_pr.q_rai[k], prog_pr.q_sno[k], aux_en.N_l[k], aux_en.N_i[k], aux_en.N_i_no_boost[k], aux_en.term_vel_liq[k], aux_en.term_vel_ice[k], aux_tc.term_vel_rain[k], aux_tc.term_vel_snow[k], ρ_c[k], p_c[k], aux_en.T[k], ts, precip_fraction,
                    (aux_en.area[k] > FT(0)) ? aux_en.ql_tendency_sedimentation[k] / aux_en.area[k] : FT(0),
                    (aux_en.area[k] > FT(0)) ? aux_en.ql_tendency_sedimentation_other[k] / aux_en.area[k] : FT(0),
                    (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation[k] / aux_en.area[k] : FT(0),
                    (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation_other[k] / aux_en.area[k] : FT(0)
                    ;
                    q_vap_sat_liq = aux_en.q_vap_sat_liq[k], q_vap_sat_ice = aux_en.q_vap_sat_ice[k], dqvdt = aux_en.dqvdt[k], dTdt = aux_en.dTdt[k], w=w[k], tke = aux_en.tke[k], S_i = S_i, τ_liq = aux_en.τ_liq[k], τ_ice = aux_en.τ_ice[k], dN_i_dz = aux_en.dN_i_dz[k], dqidz = aux_en.dqidz[k], N_INP = N_INP, massflux = aux_tc.massflux[k])
            else                
                (; mph_neq, mph_neq_other, mph_precip) = microphysics_helper(grid, edmf, moisture_model, precip_model, cloud_sedimentation_model, rain_formation_model, snow_formation_model, Δt, param_set, use_fallback_tendency_limiters, aux_en, aux_en_f, k,
                    aux_en.area[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k], prog_pr.q_rai[k], prog_pr.q_sno[k], aux_en.N_l[k], aux_en.N_i[k], aux_en.N_i_no_boost[k], aux_en.term_vel_liq[k], aux_en.term_vel_ice[k], aux_tc.term_vel_rain[k], aux_tc.term_vel_snow[k], ρ_c[k], p_c[k], aux_en.T[k], ts, precip_fraction,
                    (aux_en.area[k] > FT(0)) ? aux_en.ql_tendency_sedimentation[k] / aux_en.area[k] : FT(0),
                    (aux_en.area[k] > FT(0)) ? aux_en.ql_tendency_sedimentation_other[k] / aux_en.area[k] : FT(0),
                    (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation[k] / aux_en.area[k] : FT(0),
                    (aux_en.area[k] > FT(0)) ? aux_en.qi_tendency_sedimentation_other[k] / aux_en.area[k] : FT(0)
                )
            end


            # Storage
            aux_en.ql_tendency_noneq[k] = mph_neq.ql_tendency * aux_en.area[k]
            aux_en.qi_tendency_noneq[k] = mph_neq.qi_tendency * aux_en.area[k] # add in the sub-deposition conribution to precip formation

            aux_en.ql_tendency_cond_evap[k] = aux_en.ql_tendency_noneq[k] # for storage
            aux_en.qi_tendency_sub_dep[k] = aux_en.qi_tendency_noneq[k] # for storage

            # Storage
            aux_en.ql_tendency_noneq[k] += mph_neq_other.ql_tendency * aux_en.area[k] # is storing them w/ noneq best? should i add another one...? would need to tie it into everywhere else as well in dycore.jl etc... but wouldn't be too bad
            aux_en.qi_tendency_noneq[k] += mph_neq_other.qi_tendency * aux_en.area[k]

            aux_en.qi_tendency_hom_frz[k] = mph_neq_other.qi_tendency_homogeneous_freezing * aux_en.area[k] # for storage
            aux_en.qi_tendency_het_frz[k] = mph_neq_other.qi_tendency_heterogeneous_freezing * aux_en.area[k] # for storage
            aux_en.qi_tendency_het_nuc[k] = mph_neq_other.qi_tendency_heterogeneous_icenuc * aux_en.area[k] # for storage
            aux_en.qi_tendency_melt[k] = mph_neq_other.qi_tendency_melting * aux_en.area[k] # for storage

            # update_sat_unsat
            if TD.has_condensate(thermo_params, ts)
                aux_en.cloud_fraction[k] = 1
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
                if moisture_model isa NonEquilibriumMoisture
                    # aux_en_sat.q_vap[k] = min(TD.q_vap_saturation(thermo_params, ts), TD.vapor_specific_humidity(thermo_params, ts)) # having condensate is no guarantee of saturation in noneq
                    aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
                    aux_en_sat.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
                    aux_en_sat.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)
                else
                    aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
                end
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts)
                aux_en_sat.T[k] = TD.air_temperature(thermo_params, ts)
                aux_en_sat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)
                if moisture_model isa EquilibriumMoisture
                    aux_en.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
                    aux_en.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)

                    if reweight_processes_for_grid  # Equivalent call for updraft is in update_aux.jl. This probably shouldn't rely on  
                        q = TD.PhasePartition(thermo_params, ts_env[k]) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
                        q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_en, q, p_c; reweight_extrema_only = reweight_extrema_only)
                        aux_en.q_liq[k] = q_liq
                        aux_en.q_ice[k] = q_ice
                    end

                end # if NonEquilibriumMoisture, these should be already set in update_aux(), so calling doesn't do much (doesn't hurt either I suppose)
            else
                aux_en.cloud_fraction[k] = 0
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
                aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts)
                aux_en_unsat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)

                # technically we should still check the edges of the domain since they could be sat even if center is unsat
                if reweight_processes_for_grid && (moisture_model isa EquilibriumMoisture)
                    q = TD.PhasePartition(thermo_params, ts_env[k]) # we don't have q cached here but we do in updraft so calculating outside the fcn saves work
                    q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_en, q, p_c; reweight_extrema_only = reweight_extrema_only)
                    aux_en.q_liq[k] = q_liq
                    aux_en.q_ice[k] = q_ice
                end

                aux_en_sat.T[k] = aux_en.T[k]
                if moisture_model isa NonEquilibriumMoisture # In noneq, not having condensate doesn't guarantee unsat
                    aux_en_sat.q_vap[k] = min(TD.q_vap_saturation(thermo_params, ts), TD.vapor_specific_humidity(thermo_params, ts)) # having condensate is no guarantee of saturation in noneq
                else
                    aux_en_sat.q_vap[k] = FT(0) #  i think this is wrong...
                    # aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
                end
                aux_en_sat.q_tot[k] = aux_en.q_tot[k]
                aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
                aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]
                if moisture_model isa EquilibriumMoisture
                    aux_en.q_liq[k] = 0
                    aux_en.q_ice[k] = 0
                    aux_en_sat.q_liq[k] = 0
                    aux_en_sat.q_ice[k] = 0
                end # if NonEquilibriumMoisture, these should be already set in update_aux(), so calling doesn't do much (doesn't hurt either I suppose)
            end


            # Storage  [final]
            # update_env_precip_tendencies
            # TODO: move ..._tendency_precip_formation to diagnostics
            aux_en.qt_tendency_precip_formation[k] = mph_precip.qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = mph_precip.θ_liq_ice_tendency * aux_en.area[k]
            if moisture_model isa NonEquilibriumMoisture
                aux_en.ql_tendency_precip_formation[k] = mph_precip.ql_tendency * aux_en.area[k]
                aux_en.qi_tendency_precip_formation[k] = mph_precip.qi_tendency * aux_en.area[k]
            end
            tendencies_pr.q_rai[k] += mph_precip.qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += mph_precip.qs_tendency * aux_en.area[k]

            # store autoconversion and accretion for diagnostics (doens't mean much for Eq since liq/ice get set by sat adjust and T...)
            aux_en.ql_tendency_acnv[k] = mph_precip.ql_tendency_acnv * aux_en.area[k]
            aux_en.qi_tendency_acnv[k] = mph_precip.qi_tendency_acnv * aux_en.area[k]
            aux_en.qi_tendency_acnv_dep[k] = mph_precip.qi_tendency_acnv_dep * aux_en.area[k]
            aux_en.qi_tendency_acnv_dep_is[k] = mph_precip.qi_tendency_acnv_dep_is * aux_en.area[k]
            aux_en.qi_tendency_acnv_dep_above[k] = mph_precip.qi_tendency_acnv_dep_above * aux_en.area[k]
            aux_en.qi_tendency_acnv_agg_mix[k] = mph_precip.qi_tendency_acnv_agg_mix * aux_en.area[k] # this shoud have already been set but
            aux_en.qi_tendency_acnv_thresh[k] = mph_precip.qi_tendency_acnv_thresh * aux_en.area[k] # this shoud have already been set but
            aux_en.ql_tendency_accr_liq_rai[k] = mph_precip.ql_tendency_accr_liq_rai * aux_en.area[k]
            aux_en.ql_tendency_accr_liq_ice[k] = mph_precip.ql_tendency_accr_liq_ice * aux_en.area[k]
            aux_en.ql_tendency_accr_liq_sno[k] = mph_precip.ql_tendency_accr_liq_sno * aux_en.area[k]
            aux_en.qi_tendency_accr_ice_liq[k] = mph_precip.qi_tendency_accr_ice_liq * aux_en.area[k]
            aux_en.qi_tendency_accr_ice_rai[k] = mph_precip.qi_tendency_accr_ice_rai * aux_en.area[k]
            aux_en.qi_tendency_accr_ice_sno[k] = mph_precip.qi_tendency_accr_ice_sno * aux_en.area[k]
            #
            aux_tc.qs_tendency_accr_rai_sno[k] += mph_precip.qs_tendency_accr_rai_sno * aux_en.area[k] # we calculate in aux/en for the temperature dependence but store a combined output
    

            aux_en.Hvar_rain_dt[k] = 0
            aux_en.QTvar_rain_dt[k] = 0
            aux_en.HQTcov_rain_dt[k] = 0
        end
    end

    return nothing
end


"""
Helper function to abstract out sat_adjust and noneq_moisture_sources(), and precipitation_formation() to simplify implementation.
     calculate_sedimentation_sources() has to be precomputed because it's not clear how it works otherwise, because it depends on vertical gradients which are weird with quadrature...
"""
function microphysics_helper(grid::Grid, edmf::EDMFModel, moisture_model::AbstractMoistureModel, precip_model::AbstractPrecipitationModel, cloud_sedimentation_model::AbstractCloudSedimentationModel, rain_formation_model::AbstractRainFormationModel, snow_formation_model::AbstractSnowFormationModel, Δt::FT, param_set::APS, use_fallback_tendency_limiters::Bool,
    aux_en::CC.Fields.Field, aux_en_f::CC.Fields.Field, k::Cent{Int64},
    # vars::NamedTuple{(:qt′qt′, :qt_mean, :θl′θl′, :θl_mean, :θl′qt′, :subdomain_area, :q_rai, :q_sno, :q_liq, :q_ice, :N_i, :term_vel_ice, :term_vel_rain, :term_vel_snow, :ρ_c, :p_c, :precip_frac, :tke), T} where T,
    area::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    q_rai::FT,
    q_sno::FT,
    N_l::FT,
    N_i::FT,
    N_i_no_boost::FT,
    term_vel_liq::FT,
    term_vel_ice::FT,
    term_vel_rain::FT,
    term_vel_snow::FT,
    ρ_c::FT,
    p::FT,
    T::FT,
    ts::TD.ThermodynamicState,
    precip_fraction::FT,
    ql_tendency_sedimentation::FT,
    ql_tendency_sedimentation_other::FT,
    qi_tendency_sedimentation::FT,
    qi_tendency_sedimentation_other::FT,
    ; # Things below aren't needed in Eq
    #
    q_vap_sat_liq::FT = FT(NaN), # only needed for NonEq
    q_vap_sat_ice::FT = FT(NaN), # only needed for
    dqvdt::FT = FT(0), # only needed for NonEq
    dTdt::FT = FT(0),
    w::FT = FT(0),
    #
    tke::FT = FT(NaN),
    S_i::FT = FT(NaN),
    τ_liq::FT = FT(NaN),
    τ_ice::FT = FT(NaN),
    dN_i_dz::FT = FT(0),
    dqidz::FT = FT(0),
    N_INP::FT = FT(NaN),
    massflux::FT = FT(0),
    domain = Env,
) where {FT}


    ε = eps(FT)
    thermo_params = TCP.thermodynamics_params(param_set)
    ρ = TD.air_density(thermo_params, ts)

    reweight_processes_for_grid::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_processes_for_grid, false)
    reweight_extrema_only::Bool = TCP.get_isbits_nt(param_set.user_params, :reweight_extrema_only, false)

    if moisture_model isa NonEquilibriumMoisture
        nonequilibrium_moisture_scheme = moisture_model.scheme # my new addition
        # moisture_sources_limiter = edmf.tendency_limiters.moisture_sources_limiter
        moisture_sources_limiter = get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources), use_fallback_tendency_limiters)
    else
        moisture_sources_limiter = get_tendency_limiter(edmf.tendency_limiters, Val(:moisture_sources), use_fallback_tendency_limiters)
    end




    if moisture_model isa NonEquilibriumMoisture
        mph_neq = noneq_moisture_sources(param_set, nonequilibrium_moisture_scheme, moisture_sources_limiter, area, ρ_c, p, T, Δt + ε, ts, w, q_vap_sat_liq, q_vap_sat_ice, dqvdt, dTdt, τ_liq, τ_ice)

        # if at an extrema, we don't know where the peak is so we'll reweight based on probability of where it could be in between.
        if reweight_processes_for_grid
            error("This won't work because it takes full vector inputs rn... we're only passing in the local scalar")
            mph_neq = reweight_noneq_moisture_sources_for_grid(k, grid, param_set, thermo_params, aux_en, aux_en_f, mph_neq, nonequilibrium_moisture_scheme, moisture_sources_limiter, Δt, ρ_c, p, w, dqvdt, dTdt; reweight_extrema_only = reweight_extrema_only)
        end
        #


        # Other microphysics processes (autoconversion, accretion, aggregation, nucleation, etc)
        mph_neq_other = other_microphysics_processes(
                param_set,
                moisture_model,
                # moisture_model.heterogeneous_ice_nucleation,
                nonequilibrium_moisture_scheme,
                moisture_sources_limiter,
                area,
                ρ,
                p,
                T,
                Δt,
                ts,
                w,
                mph_neq.ql_tendency,
                mph_neq.qi_tendency;
                N_l = N_l,
                N_i = N_i
            )

            q = TD.PhasePartition(thermo_params, ts) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
            S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
            N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, T, ρ, w)
            
            # autoconversion and accretion
            mph_precip = precipitation_formation(
                param_set,
                moisture_model,
                precip_model,
                cloud_sedimentation_model,
                rain_formation_model,
                snow_formation_model,
                edmf.area_partition_model,
                q_rai,
                q_sno,
                q_liq,
                q_ice,
                N_i,
                term_vel_ice,
                term_vel_rain,
                term_vel_snow,
                area,
                ρ,
                p,
                T,
                Δt,
                ts,
                mph_neq.ql_tendency, # this is the sub-deposition contribution to precip formation
                mph_neq.qi_tendency, # this is the sub-deposition contribution to precip formation
                (area > FT(0)) ? qi_tendency_sedimentation : FT(0), # this is the sedimentation contribution to precip formation
                precip_fraction,
                get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                tke,
                N_i_no_boost,
                S_i,
                τ_ice,
                dN_i_dz,
                dqidz,
                N_INP,
                massflux,
                domain,
            )
    else

            mph_neq = null_NoneqMoistureSources(FT; fill_value = FT(NaN))
            mph_neq_other = null_OtherMicrophysicsSources(FT; fill_value = FT(NaN))

            q = TD.PhasePartition(thermo_params, ts) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
            S_i = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
            N_INP = get_INP_concentration(param_set, moisture_model.scheme, q, T, ρ, w)



            mph_precip = precipitation_formation(
                param_set,
                moisture_model,
                precip_model,
                cloud_sedimentation_model,
                rain_formation_model,
                snow_formation_model,
                edmf.area_partition_model,
                q_rai,
                q_sno,
                q_liq,
                q_ice,
                N_i,
                term_vel_ice,
                term_vel_rain,
                term_vel_snow,
                area,
                ρ,
                p,
                T,
                Δt,
                ts,
                precip_fraction,
                get_tendency_limiter(edmf.tendency_limiters, Val(:precipitation), use_fallback_tendency_limiters),
                aux_en.tke
                ;
                N_i_no_boost = N_i_no_boost,
                qi_tendency_sed = qi_tendency_sedimentation,
                S_i = S_i,
                τ_sub_dep = τ_ice,
                dN_i_dz = dN_i_dz,
                dqidz = dqidz,
                N_INP = N_INP,
                massflux = massflux,
                domain = domain,
            )

    end

    return (; mph_neq, mph_neq_other, mph_precip)


end


"""
    For future-proofing and ease of maintenance it might be good to be able to take mph_neq, mph_neq_other, mph_precip and just update tendency fields in a separate eqn.
    Note, this does not update sedimentation tendencies (I guess it could though? idk)
"""
function sgs_mean_update_tendencies!(aux_en::CC.Fields.Field, aux_tc::CC.Fields.Field, k::Cent{Int64}, param_set::APS,
    mph_neq::NoneqMoistureSources, mph_neq_other::OtherMicrophysicsSources, mph_precip::PrecipFormation,
    ts::TD.ThermodynamicState,
    )

    error("Not yet implemented")
    
    # Storage
    aux_en.ql_tendency_noneq[k] = mph_neq.ql_tendency * aux_en.area[k]
    aux_en.qi_tendency_noneq[k] = mph_neq.qi_tendency * aux_en.area[k] # add in the sub-deposition conribution to precip formation

    aux_en.ql_tendency_cond_evap[k] = aux_en.ql_tendency_noneq[k] # for storage
    aux_en.qi_tendency_sub_dep[k] = aux_en.qi_tendency_noneq[k] # for storage

    # Storage
    aux_en.ql_tendency_noneq[k] += mph_neq_other.ql_tendency * aux_en.area[k] # is storing them w/ noneq best? should i add another one...? would need to tie it into everywhere else as well in dycore.jl etc... but wouldn't be too bad
    aux_en.qi_tendency_noneq[k] += mph_neq_other.qi_tendency * aux_en.area[k]

    aux_en.qi_tendency_hom_frz[k] = mph_neq_other.qi_tendency_homogeneous_freezing * aux_en.area[k] # for storage
    aux_en.qi_tendency_het_frz[k] = mph_neq_other.qi_tendency_heterogeneous_freezing * aux_en.area[k] # for storage
    aux_en.qi_tendency_het_nuc[k] = mph_neq_other.qi_tendency_heterogeneous_icenuc * aux_en.area[k] # for storage
    aux_en.qi_tendency_melt[k] = mph_neq_other.qi_tendency_melting * aux_en.area[k] # for storage


    # # update_sat_unsat [[ TODO:  I need to figure out what we'd need to pass in here... ]]
    # if TD.has_condensate(thermo_params, ts)
    #     aux_en.cloud_fraction[k] = 1
    #     aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
    #     if moisture_model isa NonEquilibriumMoisture
    #         aux_en_sat.q_vap[k] = min(TD.q_vap_saturation(thermo_params, ts), TD.vapor_specific_humidity(thermo_params, ts)) # having condensate is no guarantee of saturation in noneq
    #     else
    #         aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
    #     end
    #     aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts)
    #     aux_en_sat.T[k] = TD.air_temperature(thermo_params, ts)
    #     aux_en_sat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)
    #     if moisture_model isa EquilibriumMoisture
    #         aux_en.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
    #         aux_en.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)


    #         if reweight_processes_for_grid  # Equivalent call for updraft is in update_aux.jl. This probably shouldn't rely on  
    #             q = TD.PhasePartition(thermo_params, ts_env[k]) # we don't have it cached here but we do in updraft so calculating outside the fcn saves work
    #             q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_en, q, p_c; reweight_extrema_only = reweight_extrema_only)
    #             aux_en.q_liq[k] = q_liq
    #             aux_en.q_ice[k] = q_ice
    #         end

    #     end # if NonEquilibriumMoisture, these should be already set in update_aux(), so calling doesn't do much (doesn't hurt either I suppose)
    # else
    #     aux_en.cloud_fraction[k] = 0
    #     aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
    #     aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts)
    #     aux_en_unsat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)

    #     # technically we should still check the edges of the domain since they could be sat even if center is unsat
    #     if reweight_processes_for_grid && (moisture_model isa EquilibriumMoisture)
    #         q = TD.PhasePartition(thermo_params, ts_env[k]) # we don't have q cached here but we do in updraft so calculating outside the fcn saves work
    #         q_liq, q_ice = reweight_equilibrium_saturation_adjustment_for_grid(k, grid, param_set, thermo_params, aux_en, q, p_c; reweight_extrema_only = reweight_extrema_only)
    #         aux_en.q_liq[k] = q_liq
    #         aux_en.q_ice[k] = q_ice
    #     end

    #     aux_en_sat.T[k] = aux_en.T[k]
    #     if moisture_model isa NonEquilibriumMoisture # In noneq, not having condensate doesn't guarantee unsat
    #         aux_en_sat.q_vap[k] = min(TD.q_vap_saturation(thermo_params, ts), TD.vapor_specific_humidity(thermo_params, ts)) # having condensate is no guarantee of saturation in noneq
    #     else
    #         aux_en_sat.q_vap[k] = FT(0) #  i think this is wrong...
    #         # aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
    #     end
    #     aux_en_sat.q_tot[k] = aux_en.q_tot[k]
    #     aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
    #     aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]
    #     if moisture_model isa EquilibriumMoisture
    #         aux_en.q_liq[k] = 0
    #         aux_en.q_ice[k] = 0
    #     end # if NonEquilibriumMoisture, these should be already set in update_aux(), so calling doesn't do much (doesn't hurt either I suppose)
    # end


    # Storage  [final]
    # update_env_precip_tendencies
    # TODO: move ..._tendency_precip_formation to diagnostics
    aux_en.qt_tendency_precip_formation[k] = mph_precip.qt_tendency * aux_en.area[k]
    aux_en.θ_liq_ice_tendency_precip_formation[k] = mph_precip.θ_liq_ice_tendency * aux_en.area[k]
    if moisture_model isa NonEquilibriumMoisture
        aux_en.ql_tendency_precip_formation[k] = mph_precip.ql_tendency * aux_en.area[k]
        aux_en.qi_tendency_precip_formation[k] = mph_precip.qi_tendency * aux_en.area[k]
    end
    tendencies_pr.q_rai[k] += mph_precip.qr_tendency * aux_en.area[k]
    tendencies_pr.q_sno[k] += mph_precip.qs_tendency * aux_en.area[k]

    # store autoconversion and accretion for diagnostics (doens't mean much for Eq since liq/ice get set by sat adjust and T...)
    aux_en.ql_tendency_acnv[k] = mph_precip.ql_tendency_acnv * aux_en.area[k]
    aux_en.qi_tendency_acnv[k] = mph_precip.qi_tendency_acnv * aux_en.area[k]
    aux_en.qi_tendency_acnv_dep[k] = mph_precip.qi_tendency_acnv_dep * aux_en.area[k]
    aux_en.qi_tendency_acnv_dep_is[k] = mph_precip.qi_tendency_acnv_dep_is * aux_en.area[k]
    aux_en.qi_tendency_acnv_dep_above[k] = mph_precip.qi_tendency_acnv_dep_above * aux_en.area[k]
    aux_en.qi_tendency_acnv_agg_mix[k] = mph_precip.qi_tendency_acnv_agg_mix * aux_en.area[k] # this shoud have already been set but
    aux_en.qi_tendency_acnv_thresh[k] = mph_precip.qi_tendency_acnv_thresh * aux_en.area[k] # this shoud have already been set but
    aux_en.ql_tendency_accr_liq_rai[k] = mph_precip.ql_tendency_accr_liq_rai * aux_en.area[k]
    aux_en.ql_tendency_accr_liq_ice[k] = mph_precip.ql_tendency_accr_liq_ice * aux_en.area[k]
    aux_en.ql_tendency_accr_liq_sno[k] = mph_precip.ql_tendency_accr_liq_sno * aux_en.area[k]
    aux_en.qi_tendency_accr_ice_liq[k] = mph_precip.qi_tendency_accr_ice_liq * aux_en.area[k]
    aux_en.qi_tendency_accr_ice_rai[k] = mph_precip.qi_tendency_accr_ice_rai * aux_en.area[k]
    aux_en.qi_tendency_accr_ice_sno[k] = mph_precip.qi_tendency_accr_ice_sno * aux_en.area[k]
    #
    aux_tc.qs_tendency_accr_rai_sno[k] += mph_precip.qs_tendency_accr_rai_sno * aux_en.area[k] # we calculate in aux/en for the temperature dependence but store a combined output
    
    return nothing
end