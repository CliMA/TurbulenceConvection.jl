
function buoyancy_c(param_set, rho0, rho)
    g = CPP.grav(param_set)
    return g * (rho0 - rho) / rho0
end

function rho_c(p0, T, qt, qv)
    return p0 / ((Rd * T) * (1.0 - qt + eps_vi * qv))
end

function eos(param_set, p0, qt, prog)
    ql = 0.0

    _ret = eos_struct()

    pv_1 = p0 * eps_vi * qt / (1.0 - qt + eps_vi * qt)
    pd_1 = p0 - pv_1
    phase_part = TD.PhasePartition(qt, ql, 0.0)
    Π = TD.exner_given_pressure(param_set, pd_1 + pv_1, phase_part)
    T_1 = prog * Π
    pv_star_1 = 6.1094 * exp((17.625 * (T_1 - 273.15)) / float((T_1 - 273.15) + 243.04)) * 100
    qv_star_1 = eps_v * (1.0 - qt) * pv_star_1 / (p0 - pv_star_1)

    # If not saturated
    if (qt <= qv_star_1)
        _ret.T = T_1
        _ret.ql = 0.0

    else

        function compute_q_liq(T)
            pv_star = 6.1094 * exp((17.625 * (T - 273.15)) / float((T - 273.15) + 243.04)) * 100
            ql = qt - (eps_v * (1.0 - qt) * pv_star / (p0 - pv_star))
            return ql
        end

        function roots(T)
            ql = compute_q_liq(T)
            return oftype(T, prog) - TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                oftype(T, p0),
                TD.PhasePartition(oftype(T, qt), oftype(T, ql), oftype(T, 0.0)),
            )
        end

        ql_1 = qt - qv_star_1
        L = TD.latent_heat_vapor(param_set, T_1)
        T_init = T_1 + ql_1 * L / ((1.0 - qt) * cpd + qv_star_1 * cpv)

        maxiter = 50
        tol = RS.SolutionTolerance(1.0e-3)
        sol = RS.find_zero(roots, RS.NewtonsMethodAD(T_init), RS.CompactSolution(), tol, maxiter)
        if !sol.converged
            println("-----------------------------------------\n")
            println("maxiter reached in eos:\n")
            println("    Method=NewtonsMethodAD")
            println(", θ_liq_ice=", prog)
            println(", p_0=", p0)
            println(", q_tot=", qt)
            println(", T=", sol.root)
            println(", maxiter=", maxiter)
            println(", tol=", tol.tol, "\n")
        end
        _ret.T = sol.root
        _ret.ql = compute_q_liq(_ret.T)
    end

    return _ret
end
