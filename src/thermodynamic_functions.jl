function sd_c(pd, T)
    return sd_tilde + cpd * log(T / T_tilde) - Rd * log(pd / p_tilde)
end

function sv_c(pv, T)
    return sv_tilde + cpv * log(T / T_tilde) - Rv * log(pv / p_tilde)
end

function theta_virt_c(p0, T, qt, ql)
    # qd = 1 - qt
    # qt = qv + ql + qi
    # qr = mr/md+mv+ml+mi
    # Ignacio: This formulation holds when qt = qv + ql (negligible/no ice)

    exner = (p0 / p_tilde)^kappa

    return T / exner * (1.0 + 0.61 * (qt - ql) - ql)
end

function pd_c(p0, qt, qv)
    return p0 * (1.0 - qt) / (1.0 - qt + eps_vi * qv)
end

function pv_c(p0, qt, qv)
    return p0 * eps_vi * qv / (1.0 - qt + eps_vi * qv)
end

function density_temperature_c(T, qt, qv)
    return T * (1.0 - qt + eps_vi * qv)
end

function theta_rho_c(p0, T, qt, qv)

    exner = (p0 / p_tilde)^kappa

    return density_temperature_c(T, qt, qv) / exner
end

function cpm_c(qt)
    return (1.0 - qt) * cpd + qt * cpv
end

function buoyancy_c(param_set, rho0, rho)
    g = CPP.grav(param_set)
    return g * (rho0 - rho) / rho0
end
