import Thermodynamics
const TD = Thermodynamics

function buoyancy_flux(param_set, shf::FT, lhf, T_b, qt_b, α0_b, ts) where {FT}
    g = FT(CPP.grav(param_set))
    molmass_ratio = FT(TD.molmass_ratio(param_set))
    lv = TD.latent_heat_vapor(param_set, T_b)
    cp_m = TD.cp_m(ts)
    return (g * α0_b / cp_m / T_b * (shf + (molmass_ratio - 1) * cp_m * T_b * lhf / lv))
end

function psi_m_unstable(zeta, zeta0)
    x = (1.0 - gamma_m * zeta)^0.25
    x0 = (1.0 - gamma_m * zeta0)^0.25
    psi_m = (2.0 * log((1.0 + x) / (1.0 + x0)) + log((1.0 + x * x) / (1.0 + x0 * x0)) - 2.0 * atan(x) + 2.0 * atan(x0))
    return psi_m
end

function psi_h_unstable(zeta, zeta0)
    y = sqrt(1.0 - gamma_h * zeta)
    y0 = sqrt(1.0 - gamma_h * zeta0)
    psi_h = 2.0 * log((1.0 + y) / (1.0 + y0))
    return psi_h
end

function psi_m_stable(zeta, zeta0)
    psi_m = -beta_m * (zeta - zeta0)
    return psi_m
end

function psi_h_stable(zeta, zeta0)
    psi_h = -beta_h * (zeta - zeta0)
    return psi_h
end

function compute_ustar(param_set, windspeed, buoyancy_flux, z0, z1)
    vkb = CPSGS.von_karman_const(param_set)
    logz = log(z1 / z0)
    # use neutral condition as first guess
    ustar0 = windspeed * vkb / logz
    monin_obukov(ustar) = -ustar^3 / (buoyancy_flux * vkb)
    if (abs(buoyancy_flux) > 1.0e-20)
        ζ_1 = z1 / monin_obukov(ustar0)
        psi_m = ζ_1 >= 0.0 ? psi_m_stable : psi_m_unstable
        function roots(ustar)
            lmo = monin_obukov(ustar)
            zeta = z1 / lmo
            zeta0 = z0 / lmo
            return windspeed - ustar / vkb * (logz - psi_m(zeta, zeta0))
        end
        sol = RS.find_zero(roots, RS.NewtonsMethodAD(ustar0), RS.CompactSolution(), RS.SolutionTolerance(1e-3), 100)
        sol.converged || error("Unconverged root solve in compute_ustar")
        return sol.root
    else
        return ustar0
    end
end

"""
    exchange_coefficients_byun(param_set, Ri, zb, z0, cm, ch, lmo)

Ref: [Byun1990](@cite)
"""
function exchange_coefficients_byun(param_set, Ri, zb, z0, cm, ch, lmo)
    von_karman_const = CPSGS.von_karman_const(param_set)

    logz = log(zb / z0)
    zfactor = zb / (zb - z0) * logz
    sb = Ri / Pr0

    if Ri > 0.0
        zeta =
            zfactor / (2.0 * beta_h * (beta_m * Ri - 1.0)) *
            ((1.0 - 2.0 * beta_h * Ri) - sqrt(1.0 + 4.0 * (beta_h - beta_m) * sb))
        lmo = zb / zeta
        zeta0 = z0 / lmo
        psi_m = psi_m_stable(zeta, zeta0)
        psi_h = psi_h_stable(zeta, zeta0)
    else
        qb = 1.0 / 9.0 * (1.0 / (gamma_m * gamma_m) + 3.0 * gamma_h / gamma_m * sb * sb)
        pb = 1.0 / 54.0 * (-2.0 / (gamma_m * gamma_m * gamma_m) + 9.0 / gamma_m * (-gamma_h / gamma_m + 3.0) * sb * sb)
        crit = qb * qb * qb - pb * pb
        if crit < 0.0
            tb = cbrt(sqrt(-crit) + abs(pb))
            zeta = zfactor * (1.0 / (3.0 * gamma_m) - (tb + qb / tb))
        else
            angle = acos(pb / sqrt(qb * qb * qb))
            zeta = zfactor * (-2.0 * sqrt(qb) * cos(angle / 3.0) + 1.0 / (3.0 * gamma_m))
        end
        lmo = zb / zeta
        zeta0 = z0 / lmo
        psi_m = psi_m_unstable(zeta, zeta0)
        psi_h = psi_h_unstable(zeta, zeta0)
    end

    cu = von_karman_const / (logz - psi_m)
    cth = von_karman_const / (logz - psi_h) / Pr0
    cm = cu * cu
    ch = cu * cth
    return cm, ch, lmo
end
