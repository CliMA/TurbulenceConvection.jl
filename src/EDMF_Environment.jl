function microphysics(
    ::SGSMean,
    grid::Grid,
    state::State,
    precip_model::AbstractPrecipitationModel,
    Δt::Real,
    param_set::APS,
)

    tendencies_pr = center_tendencies_precipitation(state)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    aux_en_sat = aux_en.sat
    aux_en_unsat = aux_en.unsat

    @inbounds for k in real_center_indices(grid)
        # condensation
        q_tot_en = aux_en.q_tot[k]
        ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], q_tot_en)
        # autoconversion and accretion
        mph = precipitation_formation(
            param_set,
            precip_model,
            prog_pr.q_rai[k],
            prog_pr.q_sno[k],
            aux_en.area[k],
            ρ0_c[k],
            Δt,
            ts,
        )

        # update_sat_unsat
        if TD.has_condensate(ts)
            aux_en.cloud_fraction[k] = 1
            aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts)
            aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts)
            aux_en_sat.T[k] = TD.air_temperature(ts)
            aux_en_sat.q_tot[k] = TD.total_specific_humidity(ts)
            aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(ts)
        else
            aux_en.cloud_fraction[k] = 0
            aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts)
            aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(ts)
            aux_en_unsat.q_tot[k] = TD.total_specific_humidity(ts)
        end

        # update_env_precip_tendencies
        # TODO: move qt_tendency_precip_formation and θ_liq_ice_tendency_precip_formation
        # to diagnostics
        aux_en.qt_tendency_precip_formation[k] = mph.qt_tendency * aux_en.area[k]
        aux_en.θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_en.area[k]
        tendencies_pr.q_rai[k] += mph.qr_tendency * aux_en.area[k]
        tendencies_pr.q_sno[k] += mph.qs_tendency * aux_en.area[k]
    end
    return nothing
end

include("quad_loop.jl")

function microphysics(en_thermo::SGSQuadrature, grid::Grid, state::State, precip_model, Δt::Real, param_set::APS)
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    aux_en_unsat = aux_en.unsat
    aux_en_sat = aux_en.sat
    tendencies_pr = center_tendencies_precipitation(state)

    #TODO - remember you output source terms multiplied by Δt (bec. of instantaneous autoconv)
    #TODO - add tendencies for gm H, QT and QR due to rain
    #TODO - if we start using eos_smpl for the updrafts calculations
    #       we can get rid of the two categories for outer and inner quad. points

    # arrays for storing quadarature points and ints for labeling items in the arrays
    # a python dict would be nicer, but its 30% slower than this (for python 2.7. It might not be the case for python 3)

    epsilon = 10e-14 # eps(float)

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
                QTvar_en = aux_en.QTvar[k],
                q_tot_en = aux_en.q_tot[k],
                Hvar_en = aux_en.Hvar[k],
                θ_liq_ice_en = aux_en.θ_liq_ice[k],
                HQTcov_en = aux_en.HQTcov[k],
                area_en = aux_en.area[k],
                q_rai = prog_pr.q_rai[k],
                q_sno = prog_pr.q_sno[k],
                ρ0_c = ρ0_c[k],
                p0_c = p0_c[k],
            )
            outer_env, outer_src = quad_loop(en_thermo, precip_model, vars, param_set, Δt)

            # update environmental variables

            # update_env_precip_tendencies
            qt_tendency = outer_src[i_Sqt]
            θ_liq_ice_tendency = outer_src[i_SH]
            qr_tendency = outer_src[i_Sqr]
            qs_tendency = outer_src[i_Sqs]
            # TODO: move qt_tendency_precip_formation and θ_liq_ice_tendency_precip_formation
            # to diagnostics
            aux_en.qt_tendency_precip_formation[k] = qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = θ_liq_ice_tendency * aux_en.area[k]

            tendencies_pr.q_rai[k] += qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += qs_tendency * aux_en.area[k]

            # update cloudy/dry variables for buoyancy in TKE
            aux_en.cloud_fraction[k] = outer_env[i_cf]
            if aux_en.cloud_fraction[k] < 1
                aux_en_unsat.q_tot[k] = outer_env[i_qt_unsat] / (1 - aux_en.cloud_fraction[k])
                T_unsat = outer_env[i_T_unsat] / (1 - aux_en.cloud_fraction[k])
                ts_unsat = TD.PhaseEquil_pTq(param_set, p0_c[k], T_unsat, aux_en_unsat.q_tot[k])
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts_unsat)
            else
                aux_en_unsat.q_tot[k] = 0
                aux_en_unsat.θ_dry[k] = 0
            end

            if aux_en.cloud_fraction[k] > 0
                aux_en_sat.T[k] = outer_env[i_T_sat] / aux_en.cloud_fraction[k]
                aux_en_sat.q_tot[k] = outer_env[i_qt_sat] / aux_en.cloud_fraction[k]
                aux_en_sat.q_vap[k] =
                    (outer_env[i_qt_sat] - outer_env[i_ql] - outer_env[i_qi]) / aux_en.cloud_fraction[k]
                ts_sat = TD.PhaseEquil_pTq(param_set, p0_c[k], aux_en_sat.T[k], aux_en_sat.q_tot[k])
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts_sat)
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts_sat)
            else
                aux_en_sat.T[k] = 0
                aux_en_sat.q_vap[k] = 0
                aux_en_sat.q_tot[k] = 0
                aux_en_sat.θ_dry[k] = 0
                aux_en_sat.θ_liq_ice[k] = 0
            end

            # update var/covar rain sources
            aux_en.Hvar_rain_dt[k] = outer_src[i_SH_H] - outer_src[i_SH] * aux_en.θ_liq_ice[k]
            aux_en.QTvar_rain_dt[k] = outer_src[i_Sqt_qt] - outer_src[i_Sqt] * aux_en.q_tot[k]
            aux_en.HQTcov_rain_dt[k] =
                outer_src[i_SH_qt] - outer_src[i_SH] * aux_en.q_tot[k] + outer_src[i_Sqt_H] -
                outer_src[i_Sqt] * aux_en.θ_liq_ice[k]

        else
            # if variance and covariance are zero do the same as in SA_mean
            ts = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])
            mph = precipitation_formation(
                param_set,
                precip_model,
                prog_pr.q_rai[k],
                prog_pr.q_sno[k],
                aux_en.area[k],
                ρ0_c[k],
                Δt,
                ts,
            )

            # update_env_precip_tendencies
            # TODO: move qt_tendency_precip_formation and θ_liq_ice_tendency_precip_formation
            # to diagnostics
            aux_en.qt_tendency_precip_formation[k] = mph.qt_tendency * aux_en.area[k]
            aux_en.θ_liq_ice_tendency_precip_formation[k] = mph.θ_liq_ice_tendency * aux_en.area[k]
            tendencies_pr.q_rai[k] += mph.qr_tendency * aux_en.area[k]
            tendencies_pr.q_sno[k] += mph.qs_tendency * aux_en.area[k]

            # update_sat_unsat
            if TD.has_condensate(ts)
                aux_en.cloud_fraction[k] = 1
                aux_en_sat.θ_dry[k] = TD.dry_pottemp(ts)
                aux_en_sat.θ_liq_ice[k] = TD.liquid_ice_pottemp(ts)
                aux_en_sat.T[k] = TD.air_temperature(ts)
                aux_en_sat.q_tot[k] = TD.total_specific_humidity(ts)
                aux_en_sat.q_vap[k] = TD.vapor_specific_humidity(ts)
            else
                aux_en.cloud_fraction[k] = 0
                aux_en_unsat.θ_dry[k] = TD.dry_pottemp(ts)
                aux_en_unsat.θ_virt[k] = TD.virtual_pottemp(ts)
                aux_en_unsat.q_tot[k] = TD.total_specific_humidity(ts)
            end

            aux_en.Hvar_rain_dt[k] = 0
            aux_en.QTvar_rain_dt[k] = 0
            aux_en.HQTcov_rain_dt[k] = 0
        end
    end

    return nothing
end
