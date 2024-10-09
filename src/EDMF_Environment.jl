function microphysics(
    ::SGSMean,
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    precip_model::AbstractPrecipitationModel,
    rain_formation_model::AbstractRainFormationModel,
    Δt::Real,
    param_set::APS,
)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    tendencies_pr = center_tendencies_precipitation(state)
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    ts_env = center_aux_environment(state).ts
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_en_sat = aux_en.sat
    aux_en_unsat = aux_en.unsat
    precip_fraction = compute_precip_fraction(edmf, state)
    ts_LCL = cloud_base(aux_en, grid, ts_env, :env)[:cloud_base_ts] # cloud base, only keep the thermodynamic state part



    local liq_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type}
    local ice_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} 
    local liq_Dmax::FT
    local ice_Dmax::FT
    local scaling_factor::FT
    local c_1::FT
    local c_2::FT


    @inbounds for k in real_center_indices(grid)
        # condensation
        ts = ts_env[k]
        if edmf.moisture_model isa NonEquilibriumMoisture
            aux_en_f = face_aux_environment(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
            w = CCO.InterpolateF2C(aux_en_f.w)
            k_int = getfield(k, Base.propertynames(k)[1]) # get the value k.i out
            w = Base.getindex(CC.Fields.field_values(getfield(w, Base.propertynames(w)[1])), k_int) # w first field is w.bcs = getfield(w,Base.propertynames(w)[1]) , then get the index from the field value
            zc = FT(grid.zc[k].z)

            mph_neq = noneq_moisture_sources(param_set, aux_en.area[k], ρ_c[k], Δt, ts, w, zc; ts_LCL = ts_LCL)

            aux_en.ql_tendency_noneq[k] = mph_neq.ql_tendency * aux_en.area[k]
            aux_en.qi_tendency_noneq[k] = mph_neq.qi_tendency * aux_en.area[k]

            aux_en.ql_tendency_cond_evap[k] = aux_en.ql_tendency_noneq[k] # for storage
            aux_en.qi_tendency_sub_dep[k] = aux_en.qi_tendency_noneq[k] # for storage


            # MOVE HETEROGENOUS FREEZ/ICENUC, MELTING, HOMOGENOUS FREEZING/ICENUC HERE -- the limiters on those processes in noneq_moisture_sources() aren't even that good because they don't account for precip anyway
            mph_other = other_microphysics_processes(param_set, aux_en.area[k], ρ_c[k], Δt, ts, w, zc, mph_neq.ql_tendency, mph_neq.qi_tendency; ts_LCL = ts_LCL) 

            aux_en.ql_tendency_noneq[k] += mph_other.ql_tendency * aux_en.area[k] # is storing them w/ noneq best? should i add another one...? would need to tie it into everywhere else as well in dycore.jl etc... but wouldn't be too bad
            aux_en.qi_tendency_noneq[k] += mph_other.qi_tendency * aux_en.area[k]


            aux_en.qi_tendency_het_nuc[k] = mph_other.qi_tendency_heterogeneous_icenuc * aux_en.area[k] # for storage



        end

        # autoconversion and accretion
        mph = precipitation_formation(
            param_set,
            precip_model,
            rain_formation_model,
            prog_pr.q_rai[k],
            prog_pr.q_sno[k],
            aux_en.area[k],
            ρ_c[k],
            Δt,
            ts,
            precip_fraction,
        )

        # update_sat_unsat
        if TD.has_condensate(thermo_params, ts)
            aux_en.cloud_fraction[k] = 1
            aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
            aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts)
            aux_en_sat.T[k] = TD.air_temperature(thermo_params, ts)
            aux_en_sat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)
            aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
            aux_en.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
            aux_en.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)
        else
            aux_en.cloud_fraction[k] = 0
            aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
            aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts)
            aux_en_unsat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)
            aux_en_sat.T[k] = aux_en.T[k]
            aux_en_sat.q_vap[k] = 0
            aux_en_sat.q_tot[k] = aux_en.q_tot[k]
            aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
            aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]
            aux_en.q_liq[k] = 0
            aux_en.q_ice[k] = 0
        end

        # update_env_precip_tendencies
        # TODO: move ..._tendency_precip_formation to diagnostics
        aux_en.qt_tendency_precip_formation[k] = mph.qt_tendency * aux_en.area[k]
        aux_en.θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_en.area[k]
        if edmf.moisture_model isa NonEquilibriumMoisture
            aux_en.ql_tendency_precip_formation[k] = mph.ql_tendency * aux_en.area[k]
            aux_en.qi_tendency_precip_formation[k] = mph.qi_tendency * aux_en.area[k]

        end
        tendencies_pr.q_rai[k] += mph.qr_tendency * aux_en.area[k]
        tendencies_pr.q_sno[k] += mph.qs_tendency * aux_en.area[k]

        # store autoconversion and accretion for diagnostics (doens't mean much for Eq since liq/ice get set by sat adjust and T...)
        aux_en.ql_tendency_acnv[k] = mph.ql_tendency_acnv * aux_en.area[k]
        aux_en.qi_tendency_acnv[k] = mph.qi_tendency_acnv * aux_en.area[k]
        aux_en.ql_tendency_accr_liq_rai[k] = mph.ql_tendency_accr_liq_rai * aux_en.area[k]
        aux_en.ql_tendency_accr_liq_ice[k] = mph.ql_tendency_accr_liq_ice * aux_en.area[k]
        aux_en.ql_tendency_accr_liq_sno[k] = mph.ql_tendency_accr_liq_sno * aux_en.area[k]
        aux_en.qi_tendency_accr_ice_liq[k] = mph.qi_tendency_accr_ice_liq * aux_en.area[k]
        aux_en.qi_tendency_accr_ice_rai[k] = mph.qi_tendency_accr_ice_rai * aux_en.area[k]
        aux_en.qi_tendency_accr_ice_sno[k] = mph.qi_tendency_accr_ice_sno * aux_en.area[k]

    end

    # ======================================================================== # To Do: Separate sedimentation out so it can be called w/ either sgs or mean...
    # sedimentation (should this maybe be a grid mean tendency?)
    if  get_isbits_nt(param_set.user_args, :use_sedimentation, false) && !get_isbits_nt(param_set.user_args, :grid_mean_sedimentation, false)# drop this eventually?

        aux_en_f = face_aux_environment(state) # state or does ts include this? I guess you'd want to move this to the calling place... to choose updraft or environment
        w = CCO.InterpolateF2C(aux_en_f.w)
        # @info w
        # @info  Base.propertynames(w)
        # k_int = getfield(k, Base.propertynames(k)[1]) # get the value k.i out
        w = [ Base.getindex(CC.Fields.field_values(getfield(w, Base.propertynames(w)[1])), getfield(k, Base.propertynames(k)[1]) )  for k in real_center_indices(grid)] # w first field is w.bcs = getfield(w,Base.propertynames(w)[1]) , then get the index from the field value (why is it in BCs after F2C? doesn't that mean it's still on F?)
        # @info w

        # sedimentation_liq_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_liq_number_concentration, nothing)
        # sedimentation_ice_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_ice_number_concentration, nothing)

        sedimentation_liq_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_liq_number_concentration, FT(NaN)) # testing NaN over nothing for type stability
        sedimentation_ice_number_concentration = get_isbits_nt(param_set.user_args, :sedimentation_ice_number_concentration, FT(NaN)) # testing NaN over nothing for type stability

        # get liquid number concentration
        if isa(sedimentation_liq_number_concentration, Number)
            N_l = FT(sedimentation_liq_number_concentration)
        elseif isa(sedimentation_liq_number_concentration, Symbol)
            # N_l = get_N_l.(param_set, sedimentation_liq_number_concentration, ts_env, w)
            N_l = get_N_l.(param_set, sedimentation_liq_number_concentration, [ts_env[k] for k in real_center_indices(grid)], w) # broadcasting not working for some reason
            # N_l = FT(NaN)
        # elseif isnothing(sedimentation_liq_number_concentration)
            # N_l = sedimentation_liq_number_concentration
        else
            error("Unsupported liquid number concentration")
        end

        # get ice number concentration
        if isa(sedimentation_ice_number_concentration, Number)
            N_i = FT(sedimentation_ice_number_concentration)
        elseif isa(sedimentation_ice_number_concentration, Symbol)
            # N_i = get_N_i.(param_set, sedimentation_ice_number_concentration, ts_env, w)
            N_i = get_N_i.(param_set, sedimentation_ice_number_concentration, [ts_env[k] for k in real_center_indices(grid)], w) # broadcasting not working for some reason
            # N_i = FT(NaN)
        # elseif isnothing(sedimentation_ice_number_concentration)
            # N_i = sedimentation_ice_number_concentration
        else
            error("Unsupported ice number concentration") 
        end

        liq_velo_scheme = get_termvel_type(get_isbits_nt(param_set.user_args, :liq_velo_scheme, :Blk1MVel)) # we dont have this anyway
        ice_velo_scheme = get_termvel_type(get_isbits_nt(param_set.user_args, :ice_velo_scheme, :Chen2022Vel)) # Blk1MVel was too fast I believe

        sedimentation_integration_method = get_isbits_nt(param_set.user_args, :sedimentation_integration_method, :upwinding)
        liq_Dmax = get_isbits_nt(param_set.user_aux, :liq_sedimentation_Dmax, FT(Inf))
        ice_Dmax = get_isbits_nt(param_set.user_aux, :ice_sedimentation_Dmax, FT(62.5e-6)) # should this default to r_ice_snow?
        liq_sedimentation_scaling_factor = get_isbits_nt(param_set.user_aux, :liq_sedimentation_scaling_factor, FT(1.0))
        ice_sedimentation_scaling_factor = get_isbits_nt(param_set.user_aux, :ice_sedimentation_scaling_factor, FT(1.0))
        mph, mph_other = calculate_sedimentation_sources(param_set, grid, ρ_c, ts_env;
            w = w, 
            area = aux_en.area, 
            grid_mean = false, 
            integration_method = sedimentation_integration_method, 
            liq_velo_scheme = liq_velo_scheme, # defined in update_aux
            ice_velo_scheme = ice_velo_scheme, # defined in update_aux
            liq_Dmax = liq_Dmax,
            ice_Dmax = ice_Dmax, 
            liq_scaling_factor = liq_sedimentation_scaling_factor,
            ice_scaling_factor = ice_sedimentation_scaling_factor,
            Nl = N_l,
            Ni = N_i,
            ) # should this be a grid mean tendency?
        # tendencies_en = TC.center_tendencies_environment(state) 
        # tendencies_gm = TC.center_tendencies_grid_mean(state)

        L_v0 = TCP.LH_v0(param_set)
        L_s0 = TCP.LH_s0(param_set)
        @inbounds for k in real_center_indices(grid)
            ql_tendency_sedimentation = mph[k].ql_tendency
            qi_tendency_sedimentation = mph[k].qi_tendency
            qt_tendency_sedimentation = ql_tendency_sedimentation + qi_tendency_sedimentation
            aux_en.ql_tendency_sedimentation[k] += ql_tendency_sedimentation # these get added to gm in compute_gm
            aux_en.qi_tendency_sedimentation[k] += qi_tendency_sedimentation
            aux_en.qt_tendency_sedimentation[k] += qt_tendency_sedimentation # used in dycore.jl (= not += , cause this doesnt seem to get reset every iteration?) (fixed in update_aux)


            Π_m = TD.exner(thermo_params, ts_env[k])
            c_pm = TD.cp_m(thermo_params, ts_env[k])
            θ_liq_ice_tendency_sedimentation = 1 / Π_m / c_pm * (L_v0 * ql_tendency_sedimentation + L_s0 * qi_tendency_sedimentation)
            aux_en.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation # adapted from microphysics_coupling.jl | precipitation_formation() | (= not += caue these don't seem to get reset every iteration?)

            # sedimentation loss into updraft (should this be allowed lol)
            # How do we know when these get set to 0 and can be used again?
            if aux_en.area[k] <  1 # we have some updrafts
                N_up = n_updrafts(edmf)
                aux_up = center_aux_updrafts(state)
                aux_bulk = center_aux_bulk(state)
                @inbounds for i in 1:N_up
                    ql_tendency_sedimentation_other = mph_other[k].ql_tendency .* (aux_up[i].area[k] ./ aux_bulk.area[k]) # was made with environment area, so aux_en[i].area
                    qi_tendency_sedimentation_other = mph_other[k].qi_tendency .* (aux_up[i].area[k] ./ aux_bulk.area[k])
                    qt_tendency_sedimentation_other = ql_tendency_sedimentation_other + qi_tendency_sedimentation_other
                    θ_liq_ice_tendency_sedimentation_other = 1 / Π_m / c_pm * (L_v0 * ql_tendency_sedimentation_other + L_s0 * qi_tendency_sedimentation_other)
                    aux_up[i].ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                    aux_up[i].qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                    aux_up[i].qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                    aux_up[i].θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other
                    aux_bulk.ql_tendency_sedimentation[k] += ql_tendency_sedimentation_other
                    aux_bulk.qi_tendency_sedimentation[k] += qi_tendency_sedimentation_other
                    aux_bulk.qt_tendency_sedimentation[k] += qt_tendency_sedimentation_other
                    aux_bulk.θ_liq_ice_tendency_sedimentation[k] += θ_liq_ice_tendency_sedimentation_other
                end
            end

        end
    end
    # ======================================================================== #

    return nothing
end

function quad_loop(en_thermo::SGSQuadrature, precip_model, rain_formation_model, vars, param_set, Δt::Real)

    env_len = 8
    src_len = 8
    i_ql, i_qi, i_T, i_cf, i_qt_sat, i_qt_unsat, i_T_sat, i_T_unsat = 1:env_len
    i_SH_qt, i_Sqt_H, i_SH_H, i_Sqt_qt, i_Sqt, i_SH, i_Sqr, i_Sqs = 1:src_len

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

    FT = eltype(ρ_c)

    inner_env = SA.MVector{env_len, FT}(undef)
    outer_env = SA.MVector{env_len, FT}(undef)
    inner_src = SA.MVector{src_len, FT}(undef)
    outer_src = SA.MVector{src_len, FT}(undef)

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
    @inbounds for idx in 1:src_len
        outer_src[idx] = 0
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

            # condensation
            ts = thermo_state_pθq(param_set, p_c, h_hat, qt_hat)
            q_liq_en = TD.liquid_specific_humidity(thermo_params, ts)
            q_ice_en = TD.ice_specific_humidity(thermo_params, ts)
            T = TD.air_temperature(thermo_params, ts)
            # autoconversion and accretion
            mph = precipitation_formation(
                param_set,
                precip_model,
                rain_formation_model,
                q_rai,
                q_sno,
                subdomain_area,
                ρ_c,
                Δt,
                ts,
                precip_frac,
            )

            # environmental variables
            inner_env[i_ql] += q_liq_en * weights[m_h] * sqpi_inv
            inner_env[i_qi] += q_ice_en * weights[m_h] * sqpi_inv
            inner_env[i_T] += T * weights[m_h] * sqpi_inv
            # cloudy/dry categories for buoyancy in TKE
            if TD.has_condensate(q_liq_en + q_ice_en)
                inner_env[i_cf] += weights[m_h] * sqpi_inv
                inner_env[i_qt_sat] += qt_hat * weights[m_h] * sqpi_inv
                inner_env[i_T_sat] += T * weights[m_h] * sqpi_inv
            else
                inner_env[i_qt_unsat] += qt_hat * weights[m_h] * sqpi_inv
                inner_env[i_T_unsat] += T * weights[m_h] * sqpi_inv
            end
            # products for variance and covariance source terms
            inner_src[i_Sqt] += mph.qt_tendency * weights[m_h] * sqpi_inv
            inner_src[i_Sqr] += mph.qr_tendency * weights[m_h] * sqpi_inv
            inner_src[i_Sqs] += mph.qs_tendency * weights[m_h] * sqpi_inv
            inner_src[i_SH] += mph.θ_liq_ice_tendency * weights[m_h] * sqpi_inv
            inner_src[i_Sqt_H] += mph.qt_tendency * h_hat * weights[m_h] * sqpi_inv
            inner_src[i_Sqt_qt] += mph.qt_tendency * qt_hat * weights[m_h] * sqpi_inv
            inner_src[i_SH_H] += mph.θ_liq_ice_tendency * h_hat * weights[m_h] * sqpi_inv
            inner_src[i_SH_qt] += mph.θ_liq_ice_tendency * qt_hat * weights[m_h] * sqpi_inv
        end

        for idx in 1:env_len
            outer_env[idx] += inner_env[idx] * weights[m_q] * sqpi_inv
        end
        for idx in 1:src_len
            outer_src[idx] += inner_src[idx] * weights[m_q] * sqpi_inv
        end
    end


    outer_src_nt = (;
        SH_qt = outer_src[i_SH_qt],
        Sqt_H = outer_src[i_Sqt_H],
        SH_H = outer_src[i_SH_H],
        Sqt_qt = outer_src[i_Sqt_qt],
        Sqt = outer_src[i_Sqt],
        SH = outer_src[i_SH],
        Sqr = outer_src[i_Sqr],
        Sqs = outer_src[i_Sqs],
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

function microphysics(
    en_thermo::SGSQuadrature,
    grid::Grid,
    state::State,
    edmf::EDMFModel,
    precip_model,
    rain_formation_model,
    Δt::Real,
    param_set::APS,
)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en_unsat = aux_en.unsat
    aux_en_sat = aux_en.sat
    tendencies_pr = center_tendencies_precipitation(state)
    ts_env = center_aux_environment(state).ts
    precip_fraction = compute_precip_fraction(edmf, state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    #TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)

    epsilon = 10e-14 # eps(float)

    if edmf.moisture_model isa NonEquilibriumMoisture
        error("The SGS quadrature microphysics is not compatible with non-equilibrium moisture")
    end

    # initialize the quadrature points and their labels

    @inbounds for k in real_center_indices(grid)
        if (
            aux_en.QTvar[k] > epsilon &&
            aux_en.Hvar[k] > epsilon &&
            abs(aux_en.HQTcov[k]) > epsilon &&
            aux_en.q_tot[k] > epsilon &&
            sqrt(aux_en.QTvar[k]) < aux_en.q_tot[k]
        )
            vars = (;
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
                zc = FT(grid.zc[k].z),
            )
            outer_env, outer_src = quad_loop(en_thermo, precip_model, rain_formation_model, vars, param_set, Δt)

            # update environmental cloudy/dry variables for buoyancy in TKE
            # update_env_precip_tendencies
            qt_tendency = outer_src.Sqt
            θ_liq_ice_tendency = outer_src.SH
            qr_tendency = outer_src.Sqr
            qs_tendency = outer_src.Sqs
            # TODO: move ..._tendency_precip_formation to diagnostics
            aux_en.qt_tendency_precip_formation[k] = qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = θ_liq_ice_tendency * aux_en.area[k]

            tendencies_pr.q_rai[k] += qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += qs_tendency * aux_en.area[k]


            if TD.has_condensate(outer_env.ql + outer_env.qi)
                aux_en.cloud_fraction[k] = outer_env.cf
                aux_en.q_liq[k] = outer_env.ql
                aux_en.q_ice[k] = outer_env.qi

            else
                aux_en.cloud_fraction[k] = 0.0
                aux_en.q_liq[k] = 0.0
                aux_en.q_ice[k] = 0.0

                aux_en.qt_tendency_precip_formation[k] = 0.0
                aux_en.θ_liq_ice_tendency_precip_formation[k] = 0.0
                tendencies_pr.q_rai[k] = 0.0
                tendencies_pr.q_sno[k] = 0.0

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
                aux_en_sat.T[k] = aux_en.T[k]
                aux_en_sat.q_vap[k] = 0
                aux_en_sat.q_tot[k] = aux_en.q_tot[k]
                aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
                aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]

            end

            # update var/covar rain sources
            aux_en.Hvar_rain_dt[k] = outer_src.SH_H - outer_src.SH * aux_en.θ_liq_ice[k]
            aux_en.QTvar_rain_dt[k] = outer_src.Sqt_qt - outer_src.Sqt * aux_en.q_tot[k]
            aux_en.HQTcov_rain_dt[k] =
                outer_src.SH_qt - outer_src.SH * aux_en.q_tot[k] + outer_src.Sqt_H - outer_src.Sqt * aux_en.θ_liq_ice[k]


        else
            # if variance and covariance are zero do the same as in SA_mean
            ts = ts_env[k]
            mph = precipitation_formation(
                param_set,
                precip_model,
                rain_formation_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_en.area[k],
                ρ_c[k],
                Δt,
                ts,
                precip_fraction,
            )

            # update_env_precip_tendencies
            # TODO: move ..._tendency_precip_formation to diagnostics
            aux_en.qt_tendency_precip_formation[k] = mph.qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_en.area[k]
            tendencies_pr.q_rai[k] += mph.qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += mph.qs_tendency * aux_en.area[k]

            # update_sat_unsat
            if TD.has_condensate(thermo_params, ts)
                aux_en.cloud_fraction[k] = 1
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(thermo_params, ts)
                aux_en_sat.T[k] = TD.air_temperature(thermo_params, ts)
                aux_en_sat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)
                aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(thermo_params, ts)
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
                aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts)
                aux_en_unsat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)

            else
                aux_en.cloud_fraction[k] = 0
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(thermo_params, ts)
                aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(thermo_params, ts)
                aux_en_unsat.q_tot[k] = TD.total_specific_humidity(thermo_params, ts)

                aux_en_sat.T[k] = aux_en.T[k]
                aux_en_sat.q_vap[k] = 0
                aux_en_sat.q_tot[k] = aux_en.q_tot[k]
                aux_en_sat.θ_dry[k] = aux_en.θ_dry[k]
                aux_en_sat.θ_liq_ice[k] = aux_en.θ_liq_ice[k]

            end

            aux_en.Hvar_rain_dt[k] = 0
            aux_en.QTvar_rain_dt[k] = 0
            aux_en.HQTcov_rain_dt[k] = 0
        end
    end

    return nothing
end
