
function buoyancy_c(param_set, rho0, rho)
    g = CPP.grav(param_set)
    return g * (rho0 - rho) / rho0
end

function rho_c(p0, T, qt, qv)
    return p0 / ((Rd * T) * (1.0 - qt + eps_vi * qv))
end