
function buoyancy_c(param_set, rho0, rho)
    g = CPP.grav(param_set)
    return g * (rho0 - rho) / rho0
end

function eos(param_set, p0, qt, prog)
    qv = qt
    ql = 0.0

    _ret = eos_struct()

    pv_1 = p0 * eps_vi * qt / (1.0 - qt + eps_vi * qt)
    pd_1 = p0 - pv_1
    phase_part = TD.PhasePartition(qt, ql, 0.0)
    Π = TD.exner_given_pressure(param_set, pd_1 + pv_1, phase_part)
    T_1 = prog * Π
    pv_star_1 = 6.1094 * exp((17.625 * (T_1 - 273.15)) / float((T_1 - 273.15) + 243.04)) * 100
    qv_star_1 = eps_v * (1.0 - qt) * pv_star_1 / (p0 - pv_star_1)

    ql_2 = 0.0
    # If not saturated
    if (qt <= qv_star_1)
        _ret.T = T_1
        _ret.ql = 0.0

    else
        ql_1 = qt - qv_star_1
        prog_1 = TD.liquid_ice_pottemp_given_pressure(param_set, T_1, p0, TD.PhasePartition(qt, ql_1, 0.0))

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

            prog_2 = TD.liquid_ice_pottemp_given_pressure(param_set, T_2, p0, TD.PhasePartition(qt, ql_2, 0.0))
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
