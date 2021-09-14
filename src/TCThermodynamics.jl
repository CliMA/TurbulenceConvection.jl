"""
    TCThermodynamics

WARNING: This is a temporary module that we're using to replace
the old `eos` with thermodynamic states in order to allow for
refactoring without needing to maintain/rebase this currently
working solution.

When we work out the changes in Thermodynamics, we will delete
this entire file.
"""
module TCThermodynamics

import CLIMAParameters
import RootSolvers
import Thermodynamics
const TD = Thermodynamics
const RS = RootSolvers
const CPP = CLIMAParameters.Planet
const APS = CLIMAParameters.AbstractEarthParameterSet


"""
    AnelasticPhaseEquil{FT} <: AbstractPhaseEquil{FT}

A thermodynamic state assuming thermodynamic equilibrium (therefore, saturation adjustment
may be needed).

# Constructors

    PhaseEquil(param_set, ρ, p, e_int, q_tot, T)

# Fields
"""
struct AnelasticPhaseEquil{FT, PS} <: TD.AbstractPhaseEquil{FT}
    "parameter set, used to dispatch planet parameter function calls"
    param_set::PS
    "full density of air (potentially moist)"
    ρ::FT
    "reference air pressure"
    p::FT
    "internal energy"
    e_int::FT
    "total specific humidity"
    q_tot::FT
    "temperature: computed via [`saturation_adjustment`](@ref)"
    T::FT
end

const ITERTYPE = Union{Int, Nothing}
TOLTYPE(FT) = Union{FT, Nothing}

"""
    saturation_adjustment_given_pθq_anelastic(
        param_set,
        p,
        θ_liq_ice,
        q_tot,
        phase_type,
        maxiter,
        tol
    )

Compute the temperature `T` that is consistent with

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `tol` absolute tolerance for saturation adjustment iterations. Can be one of:
    - `SolutionTolerance()` to stop when `|x_2 - x_1| < tol`
    - `ResidualTolerance()` to stop when `|f(x)| < tol`
 - `maxiter` maximum iterations for non-linear equation solve

by finding the root of

`θ_{liq_ice} - liquid_ice_pottemp_given_pressure(param_set, T, p, phase_type, q_tot) = 0`

This saturation adjustment is consistent with the anelastic approximation by avoiding computing
air_density in the process and receiving the pressure as an input

See also [`saturation_adjustment`](@ref).
"""
function saturation_adjustment_given_pθq_anelastic(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    phase_type::Type{<:AnelasticPhaseEquil},
    maxiter::Int,
    tol::RS.AbstractTolerance,
) where {FT <: Real}
    _T_min::FT = CPP.T_min(param_set)
    _cp_d::FT = CPP.cp_d(param_set)
    _cp_v::FT = CPP.cp_v(param_set)
    T_1 = max(_T_min, TD.air_temperature_given_pθq(param_set, p, θ_liq_ice, TD.PhasePartition(q_tot))) # Assume all vapor
    q_v_sat_1 = q_vap_saturation_from_partial_pressures(param_set, q_tot, p, T_1, TD.Liquid())
    unsaturated = q_tot <= q_v_sat_1
    if unsaturated && T_1 > _T_min
        return T_1
    else
        T_2 = TD.air_temperature_given_pθq(param_set, p, θ_liq_ice, TD.PhasePartition(q_tot, FT(0), q_tot)) # Assume all ice
        T_2 = TD.bound_upper_temperature(T_1, T_2)
        function roots(T)
            q_v_sat = q_vap_saturation_from_partial_pressures(param_set, oftype(T, q_tot), oftype(T, p), T, TD.Liquid())
            q_liq = q_tot - q_v_sat
            return oftype(T, θ_liq_ice) - TD.liquid_ice_pottemp_given_pressure(
                param_set,
                T,
                oftype(T, p),
                TD.PhasePartition(oftype(T, q_tot), oftype(T, q_liq), oftype(T, 0.0)),
            )
        end
        sol = RS.find_zero(roots, RS.SecantMethod(T_1, T_2), RS.CompactSolution(), tol, maxiter)
        if !sol.converged
            print("-----------------------------------------\n")
            print("maxiter reached in saturation_adjustment_given_pθq_anelastic:\n")
            print("    Method=SecantMethod")
            print(", p=", p)
            print(", θ_liq_ice=", θ_liq_ice)
            print(", q_tot=", q_tot)
            print(", T=", sol.root)
            print(", maxiter=", maxiter)
            print(", tol=", tol.tol, "\n")
            error("Failed to converge with printed set of inputs.")
        end
        return sol.root
    end
end


"""
    AnelasticPhaseEquil_pθq(param_set, θ_liq_ice, q_tot)

Constructs a [`AnelasticPhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` anelastic reference air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `temperature_tol` temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
"""
function AnelasticPhaseEquil_pθq(
    param_set::APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
) where {FT <: Real, IT <: ITERTYPE, FTT <: TOLTYPE(FT)}
    maxiter === nothing && (maxiter = 50)
    temperature_tol === nothing && (temperature_tol = FT(1e-3))
    phase_type = AnelasticPhaseEquil
    tol = RS.ResidualTolerance(temperature_tol)
    T = saturation_adjustment_given_pθq_anelastic(param_set, p, θ_liq_ice, q_tot, phase_type, maxiter, tol)
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)
    ρ = TD.air_density(param_set, T, p, q_pt)
    e_int = TD.internal_energy(param_set, T, q_pt)
    return AnelasticPhaseEquil{FT, typeof(param_set)}(param_set, ρ, p, e_int, q_tot, T)
end


TD.air_pressure(ts::AnelasticPhaseEquil) = ts.p

function TD.PhasePartition(ts::AnelasticPhaseEquil)
    return PhasePartition_equil_given_p(
        ts.param_set,
        TD.air_temperature(ts),
        TD.air_pressure(ts),
        TD.total_specific_humidity(ts),
        AnelasticPhaseEquil,
    )
end


"""
    PhasePartition_equil_given_p(param_set, T, p, q_tot, phase_type)

Partition the phases in equilibrium, returning a [`PhasePartition`](@ref) object using the
[`liquid_fraction`](@ref) function where
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `T` temperature
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
The residual `q.tot - q.liq - q.ice` is the vapor specific humidity.
"""
function PhasePartition_equil_given_p(
    param_set::APS,
    T::FT,
    p::FT,
    q_tot::FT,
    phase_type::Type{<:TD.ThermodynamicState},
) where {FT <: Real}

    q_v_sat = q_vap_saturation_from_partial_pressures(param_set, q_tot, p, T, TD.Liquid())
    _liquid_frac = TD.liquid_fraction(param_set, T, phase_type) # fraction of condensate that is liquid
    q_c = q_tot - q_v_sat
    q_liq = _liquid_frac * q_c
    q_ice = (1 - _liquid_frac) * q_c
    return TD.PhasePartition(q_tot, q_liq, q_ice)
end

"""
    q_vap_saturation_from_partial_pressures(param_set, q_tot, p, T, ::Phase)

Alternative method to compute the saturation specific humidity without using density, given
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `q_tot` total water specific humidity,
 - `p` air pressure,
 - `T` air temperature
-  `phase` Phase

This is useful for anelastic model where air_density should be avoided
"""
function q_vap_saturation_from_partial_pressures(
    param_set::APS,
    q_tot::FT,
    p::FT,
    T::FT,
    phase::TD.Phase,
) where {FT <: Real}
    _R_v::FT = CPP.R_v(param_set)
    _R_d::FT = CPP.R_d(param_set)
    p_v_sat = TD.saturation_vapor_pressure(param_set, T, phase)
    return _R_d / _R_v * (1 - q_tot) * p_v_sat / (p - p_v_sat)
end


end
