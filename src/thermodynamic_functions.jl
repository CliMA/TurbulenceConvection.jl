module TCThermodynamics

# TODO: use CLIMAParameters instead!
using ..TurbulenceConvection: eps_v
using ..TurbulenceConvection: cpd
using ..TurbulenceConvection: cpv
using ..TurbulenceConvection: eps_vi
using ..TurbulenceConvection: Rd

import Thermodynamics
const TD = Thermodynamics

import RootSolvers
const RS = RootSolvers

import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet

Base.@kwdef mutable struct eos_struct
    p0::Float64 = 0
    qt::Float64 = 0
    T::Float64 = 0
    ql::Float64 = 0
end

air_temperature(sa::eos_struct) = sa.T
liquid_specific_humidity(sa::eos_struct) = sa.ql
air_density(sa::eos_struct) = rho_c(sa.p0, sa.T, sa.qt, sa.qt - sa.ql)

function rho_c(p0, T, qt, qv)
    return p0 / ((Rd * T) * (1.0 - qt + eps_vi * qv))
end

function eos(param_set, p0, prog::FT, qt) where {FT}
    ql = 0.0

    _ret = eos_struct(; p0 = p0, qt = qt)

    pv_1 = p0 * eps_vi * qt / (1.0 - qt + eps_vi * qt)
    pd_1 = p0 - pv_1
    phase_part = TD.PhasePartition(qt, ql, 0.0)
    Π = TD.exner_given_pressure(param_set, pd_1 + pv_1, phase_part)
    T_1 = prog * Π
    pv_star_1 = TD.saturation_vapor_pressure(param_set, T_1, TD.Liquid())
    qv_star_1 = eps_v * (1.0 - qt) * pv_star_1 / (p0 - pv_star_1)

    # If not saturated
    if (qt <= qv_star_1)
        _ret.T = T_1
        _ret.ql = 0.0

    else

        function compute_q_liq(T)
            pv_star = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
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
        _ret.T = sol.root
        _ret.ql = compute_q_liq(_ret.T)
    end

    return _ret
end

end
