function sd_c(pd, T)
    return sd_tilde + cpd * log(T / T_tilde) - Rd * log(pd / p_tilde)
end

function sv_c(pv, T)
    return sv_tilde + cpv * log(T / T_tilde) - Rv * log(pv / p_tilde)
end

function buoyancy_c(param_set, rho0, rho)
    g = CPP.grav(param_set)
    return g * (rho0 - rho) / rho0
end
