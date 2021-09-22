module TCThermodynamics

import Thermodynamics
const TD = Thermodynamics

import RootSolvers
const RS = RootSolvers

import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet

struct eos_struct{PS, FT}
    param_set::PS
    p0::FT
    qt::FT
    T::FT
    ql::FT
end

air_temperature(sa::eos_struct) = sa.T
liquid_specific_humidity(sa::eos_struct) = sa.ql
air_density(sa::eos_struct) = rho_c(sa.param_set, sa.p0, sa.T, sa.qt, sa.qt - sa.ql)

function rho_c(param_set, p0, T, qt, qv::FT) where {FT}
    molmass_ratio = FT(CPP.molmass_ratio(param_set))
    R_d = FT(CPP.R_d(param_set))
    return p0 / ((R_d * T) * (1 - qt + molmass_ratio * qv))
end

function eos(param_set, p0, prog::FT, qt) where {FT}
    ql = 0.0

    molmass_ratio = FT(CPP.molmass_ratio(param_set))
    cp_d = FT(CPP.cp_d(param_set))
    cp_v = FT(CPP.cp_v(param_set))

    pv_1 = p0 * molmass_ratio * qt / (1 - qt + molmass_ratio * qt)
    pd_1 = p0 - pv_1
    phase_part = TD.PhasePartition(qt, ql, 0.0)
    Π = TD.exner_given_pressure(param_set, pd_1 + pv_1, phase_part)
    T_1 = prog * Π
    pv_star_1 = TD.saturation_vapor_pressure(param_set, T_1, TD.Liquid())
    qv_star_1 = (1 / molmass_ratio) * (1 - qt) * pv_star_1 / (p0 - pv_star_1)

    # If not saturated
    if (qt <= qv_star_1)
        T_sol = T_1
        ql_sol = 0.0

    else

        function compute_q_liq(T)
            pv_star = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
            ql = qt - ((1 / molmass_ratio) * (1 - qt) * pv_star / (p0 - pv_star))
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
        T_init = T_1 + ql_1 * L / ((1 - qt) * cp_d + qv_star_1 * cp_v)

        maxiter = 50
        tol = RS.SolutionTolerance(1.0e-3)
        # TODO: remove these hard-coded bounds
        sol = RS.find_zero(roots, RS.SecantMethod(T_1, T_init + FT(10)), RS.CompactSolution(), tol, maxiter)
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
            error("Halting execution")
        end
        T_sol = sol.root
        ql_sol = compute_q_liq(T_sol)
    end

    return eos_struct(param_set, p0, qt, T_sol, ql_sol)
end

end
