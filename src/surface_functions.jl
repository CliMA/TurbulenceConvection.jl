
function buoyancy_flux(param_set, shf::FT, lhf, T_b, qt_b, α0_b, ts) where {FT}
    g = FT(CPP.grav(param_set))
    molmass_ratio = FT(CPP.molmass_ratio(param_set))
    lv = TD.latent_heat_vapor(param_set, T_b)
    cp_m = TD.cp_m(ts)
    return (g * α0_b / cp_m / T_b * (shf + (molmass_ratio - 1) * cp_m * T_b * lhf / lv))
end

function compute_friction_velocity(param_set, windspeed, buoy_flux, z0, z1)
    vkb = CPSGS.von_karman_const(param_set)
    logz = log(z1 / z0)
    # use neutral condition as first guess
    ustar0 = windspeed * vkb / logz
    if (abs(buoy_flux) > 1.0e-20)
        function roots(ustar)
            lmo = obukhov_length(param_set, ustar, buoy_flux)
            uf = UF.Businger(param_set, lmo)
            ζ = z1 / lmo
            ζ_0 = z0 / lmo
            Ψ_m_1 = UF.psi(uf, ζ, UF.MomentumTransport())
            Ψ_m_0 = UF.psi(uf, ζ_0, UF.MomentumTransport())
            return windspeed - ustar / vkb * (logz - Ψ_m_1 + Ψ_m_0)
        end
        sol = RS.find_zero(roots, RS.NewtonsMethodAD(ustar0), RS.CompactSolution(), RS.SolutionTolerance(1e-3), 100)
        sol.converged || error("Unconverged root solve in compute_friction_velocity")
        return sol.root
    else
        return ustar0
    end
end

"""
    exchange_coefficients_byun(param_set, Ri, zb, z0)

Ref: [Byun1990](@cite)
"""
function exchange_coefficients_byun(param_set, Ri, zb, z0)
    von_karman_const = CPSGS.von_karman_const(param_set)

    zb = zb.z
    logz = log(zb / z0)
    zfactor = zb / (zb - z0) * logz
    sb = Ri / Pr0

    if Ri > 0.0
        ζ =
            zfactor / (2.0 * beta_h * (beta_m * Ri - 1.0)) *
            ((1.0 - 2.0 * beta_h * Ri) - sqrt(1.0 + 4.0 * (beta_h - beta_m) * sb))
    else
        qb = 1.0 / 9.0 * (1.0 / (gamma_m * gamma_m) + 3.0 * gamma_h / gamma_m * sb * sb)
        pb = 1.0 / 54.0 * (-2.0 / (gamma_m * gamma_m * gamma_m) + 9.0 / gamma_m * (-gamma_h / gamma_m + 3.0) * sb * sb)
        crit = qb * qb * qb - pb * pb
        if crit < 0.0
            tb = cbrt(sqrt(-crit) + abs(pb))
            ζ = zfactor * (1.0 / (3.0 * gamma_m) - (tb + qb / tb))
        else
            angle = acos(pb / sqrt(qb * qb * qb))
            ζ = zfactor * (-2.0 * sqrt(qb) * cos(angle / 3.0) + 1.0 / (3.0 * gamma_m))
        end
    end
    lmo = zb / ζ
    ζ_0 = z0 / lmo
    uf = SF.Businger(param_set, lmo)
    Ψ_m_1 = UF.psi(uf, ζ, SF.MomentumTransport())
    Ψ_m_0 = UF.psi(uf, ζ_0, SF.MomentumTransport())
    Ψ_h_1 = UF.psi(uf, ζ, SF.HeatTransport())
    Ψ_h_0 = UF.psi(uf, ζ_0, SF.HeatTransport())

    cu = von_karman_const / (logz - Ψ_m_1 + Ψ_m_0)
    cth = von_karman_const / (logz - Ψ_h_1 + Ψ_h_0) / Pr0
    cm = cu * cu
    ch = cu * cth
    return cm, ch, lmo
end

function obukhov_length(param_set, ustar, buoy_flux)
    von_karman_const = CPSGS.von_karman_const(param_set)
    if abs(buoy_flux) > 1e-10
        lmo = -ustar^3 / buoy_flux / von_karman_const
    else
        lmo = 0.0
    end
    return lmo
end

function friction_velocity_given_windspeed(cm, windspeed)
    return sqrt(cm) * windspeed
end
