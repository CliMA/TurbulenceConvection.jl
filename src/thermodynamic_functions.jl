
function exner_c(p0; kappa = kappa)
    return (p0 / p_tilde)^kappa
end

function theta_c(p0, T)
    return T / exner_c(p0)
end

function thetali_c(param_set, p0, T, qt, ql, qi)
    # Liquid ice potential temperature consistent with Triopoli and Cotton (1981)
    L = TD.latent_heat_vapor(param_set, T)
    return theta_c(p0, T) * exp(-L * (ql / (1.0 - qt) + qi / (1.0 - qt)) / (T * cpd))
end
function theta_virt_c(p0, T, qt, ql)
    # qd = 1 - qt
    # qt = qv + ql + qi
    # qr = mr/md+mv+ml+mi
    # Ignacio: This formulation holds when qt = qv + ql (negligible/no ice)
    return theta_c(p0, T) * (1.0 + 0.61 * (qt - ql) - ql)
end

function density_temperature_c(T, qt, qv)
    return T * (1.0 - qt + eps_vi * qv)
end
function theta_rho_c(p0, T, qt, qv)
    return density_temperature_c(T, qt, qv) / exner_c(p0)
end

function buoyancy_c(param_set, rho0, rho)
    g = CPP.grav(param_set)
    return g * (rho0 - rho) / rho0
end

function alpha_c(p0, T, qt, qv)
    return (Rd * T) / p0 * (1.0 - qt + eps_vi * qv)
end
function rho_c(p0, T, qt, qv)
    return p0 / ((Rd * T) * (1.0 - qt + eps_vi * qv))
end

function t_to_thetali_c(param_set, p0, T, qt, ql, qi)
    return thetali_c(param_set, p0, T, qt, ql, qi)
end

function eos(param_set, p0, qt, prog)
    qv = qt
    ql = 0.0

    _ret = eos_struct()

    pv_1 = p0 * eps_vi * qt / (1.0 - qt + eps_vi * qt)
    pd_1 = p0 - pv_1
    # In the original version a completly dry (qt=0) is set - changing this effect the results!
    pp = TD.PhasePartition(0.0, 0.0, 0.0) 
    T_1 = prog * TD.exner_given_pressure(param_set, pd_1 + pv_1, pp)
    pv_star_1 = 6.1094 * exp((17.625 * (T_1 - 273.15)) / float((T_1 - 273.15) + 243.04)) * 100
    qv_star_1 = eps_v * (1.0 - qt) * pv_star_1 / (p0 - pv_star_1)

    ql_2 = 0.0
    # If not saturated
    if (qt <= qv_star_1)
        _ret.T = T_1
        _ret.ql = 0.0

    else
        ql_1 = qt - qv_star_1
        prog_1 = t_to_thetali_c(param_set, p0, T_1, qt, ql_1, 0.0)
        f_1 = prog - prog_1
        L = TD.latent_heat_vapor(param_set, T_1)
        T_2 = T_1 + ql_1 * L / ((1.0 - qt) * cpd + qv_star_1 * cpv)
        delta_T = abs(T_2 - T_1)
        qv_star_2 = 0
        while delta_T > 1.0e-3 || ql_2 < 0.0
            pv_star_2 = 6.1094 * exp((17.625 * (T_2 - 273.15)) / float((T_2 - 273.15) + 243.04)) * 100
            qv_star_2 = eps_v * (1.0 - qt) * pv_star_2 / (p0 - pv_star_2)
            pv_2 = p0 * eps_vi * qv_star_2 / (1.0 - qt + eps_vi * qv_star_2)
            pd_2 = p0 - pv_2
            ql_2 = qt - qv_star_2
            prog_2 = t_to_thetali_c(param_set, p0, T_2, qt, ql_2, 0.0)
            f_2 = prog - prog_2
            T_n = T_2 - f_2 * (T_2 - T_1) / (f_2 - f_1)
            T_1 = T_2
            T_2 = T_n
            f_1 = f_2
            delta_T = abs(T_2 - T_1)
        end

        _ret.T = T_2
        qv = qv_star_2
        _ret.ql = ql_2
    end

    return _ret
end
