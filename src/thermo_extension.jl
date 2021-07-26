# TODO: move to Thermodynamics.jl and CLIMAParameters.jl
using Thermodynamics
using KernelAbstractions: @print
using RootSolvers:
    find_zero,
    CompactSolution,
    VerboseSolution,
    ResidualTolerance,
    SolutionTolerance,
    NewtonsMethod,
    NewtonsMethodAD,
    SecantMethod,
    RegulaFalsiMethod

const TD = Thermodynamics
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters.Planet:
    MSLP,
    R_d,
    R_v,
    cp_d,
    cp_v,
    cv_d,
    T_min,
    T_max,
    grav,
    T_freeze,
    entropy_reference_temperature,
    entropy_dry_air,
    entropy_water_vapor

const APS = AbstractEarthParameterSet

# The standard entropy value for dry air is computed based
# on the reference data given in Lemmon et al. [2000], and
# the standard entropy value for water vapor is based on
# the reference data given in Chase [1998].
function entropy(param_set::APS, p::FT, T::FT, q::TD.PhasePartition{FT}) where {FT}
    ps = param_set
    s̃_d = FT(entropy_dry_air(ps))
    s̃_v = FT(entropy_water_vapor(ps))
    T̃ = FT(entropy_reference_temperature(ps))
    p̃ = FT(MSLP(ps))
    q_tot = q.tot
    q_liq = q.liq
    q_ice = q.ice
    εᵥ⁻¹ = FT(R_v(ps)) / FT(R_d(ps))
    q_vap = q_tot - (q_liq + q_ice)
    p_d = p * (1 - q_tot) / (1 - q_tot + εᵥ⁻¹ * q_vap)
    p_v = p * (εᵥ⁻¹ * q_tot) / (1 - q_tot + εᵥ⁻¹ * q_vap)
    s_d = s̃_d + FT(cp_d(ps)) * log(T / T̃) - FT(R_d(ps)) * log(p_d / p̃)
    if p_v / p̃ < eps(FT) # avoid log(0)
        s_v = FT(0)
    else
        s_v = s̃_v + FT(cp_v(ps)) * log(T / T̃) - FT(R_v(ps)) * log(p_v / p̃)
    end
    L_v = TD.latent_heat_vapor(ps, T)
    L_s = TD.latent_heat_sublim(ps, T)
    return (1 - q_tot) * s_d + q_tot * s_v - (q_liq * L_v + q_ice * L_s) / T
end
entropy(param_set::APS, p::FT, T::FT) where {FT} = entropy(param_set, p, T, TD.PhasePartition{FT}(0, 0, 0))

function entropy_pTq(param_set::APS, p::FT, T::FT, q_tot::FT, phase_type::Type{<:PhaseEquil}) where {FT}
    ρ = air_density_equil(param_set, p, T, q_tot)
    q_eq = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    return entropy(param_set, p, T, q_eq)
end


entropy(ts::ThermodynamicState) =
    entropy(ts.param_set, TD.air_pressure(ts), TD.air_temperature(ts), TD.PhasePartition(ts))

function sa_numerical_method_psq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    s::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
) where {FT, NM <: RegulaFalsiMethod}
    T_1 = FT(T_min(param_set))
    T_2 = FT(T_max(param_set))
    return RegulaFalsiMethod(T_1, T_2)
end

"""
    saturation_adjustment_psq(
        sat_adjust_method,
        param_set,
        s,
        p,
        q_tot,
        phase_type,
        maxiter,
        temperature_tol
    )

Compute the temperature that is consistent with

 - `sat_adjust_method` the numerical method to use.
    See the [`Thermodynamics`](@ref) for options.
 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `s` entropy
 - `p` air pressure
 - `q_tot` total specific humidity
 - `phase_type` a thermodynamic state type
 - `maxiter` maximum iterations for non-linear equation solve
 - `temperature_tol` temperature tolerance

by finding the root of

`s - entropy(param_set, p, T, q) = 0`

using the given numerical method `sat_adjust_method`.

See also [`saturation_adjustment`](@ref).
"""
function saturation_adjustment_psq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    s::FT,
    q_tot::FT,
    phase_type::Type{<:PhaseEquil},
    maxiter::Int,
    temperature_tol::FT,
) where {FT <: Real, sat_adjust_method}
    _T_min::FT = T_min(param_set)
    _T_max::FT = T_max(param_set)
    tol = SolutionTolerance(temperature_tol)

    # TODO: use better T_1 estimate
    T_1 = _T_min # Assume all vapor
    ρ_T(T) = air_density(param_set, T, p, PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 > _T_min
        return T_1
    else

        _T_freeze::FT = T_freeze(param_set)
        s_upper = entropy_pTq(
            param_set,
            p,
            _T_freeze + temperature_tol / 2, # /2 => resulting interval is `temperature_tol` wide
            q_tot,
            phase_type,
        )
        s_lower = entropy_pTq(
            param_set,
            p,
            _T_freeze - temperature_tol / 2, # /2 => resulting interval is `temperature_tol` wide
            q_tot,
            phase_type,
        )
        if s_lower < s < s_upper
            return _T_freeze
        end
        sol = find_zero(
            T -> begin
                _q_tot = oftype(T, q_tot)
                _p = oftype(T, p)

                _T = TD.heavisided(T)
                _T = max(_T, _T_min)
                _T = min(_T, _T_max)

                entropy_pTq(param_set, _p, _T, _q_tot, phase_type) - s
            end,
            sa_numerical_method_psq(sat_adjust_method, param_set, p, s, q_tot, phase_type),
            CompactSolution(),
            tol,
            maxiter,
        )
        if !sol.converged
            if TD.print_warning()
                @print("-----------------------------------------\n")
                @print("maxiter reached in saturation_adjustment_psq:\n")
                TD.print_numerical_method(sat_adjust_method)
                @print(", p=", p)
                @print(", s=", s)
                @print(", q_tot=", q_tot)
                @print(", T=", sol.root)
                @print(", maxiter=", maxiter)
                @print(", tol=", tol.tol, "\n")
            end
            if TD.error_on_non_convergence()
                error("Failed to converge with printed set of inputs.")
            end
        end
        return sol.root
    end
end

"""
    air_density_equil(
        param_set,
        p,
        T,
        q_tot,
    )

Compute the (equilibrium-)air density `ρ` that is consistent with

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `T` temperature
 - `q_tot` total specific humidity

by finding the root of

```
T - TD.air_temperature_from_ideal_gas_law(
    param_set,
    p,
    ρ,
    PhasePartition_equil(
                    param_set,
                    T,
                    ρ,
                    q_tot)
)
```

This air density is in sync with that computed from `air_density`,
which relies on the full `PhasePartition`.

See also [`saturation_adjustment`](@ref).

!!! warn
    This is an expensive function, and has internalized
    iterative solver parameters to simplify the interface.
"""
function air_density_equil(param_set::APS, p::FT, T::FT, q_tot::FT) where {FT <: Real}
    # Assume all liquid
    ρ_init = air_density(param_set, T, p, PhasePartition(q_tot))
    phase_type = PhaseEquil
    maxiter = 50
    tol = SolutionTolerance(sqrt(eps(FT)))
    sol = find_zero(
        ρ -> begin
            _q_tot = oftype(ρ, q_tot)
            _p = oftype(ρ, p)

            q_pt = PhasePartition_equil(param_set, oftype(ρ, T), ρ, oftype(ρ, q_tot), phase_type)
            # TODO: is one of these equations preferred over the other?
            # (e.g., physically / numerically)
            # T - TD.air_temperature_from_ideal_gas_law(param_set, _p, ρ, q_pt)
            ρ - TD.air_density(param_set, oftype(ρ, T), _p, q_pt)
        end,
        NewtonsMethodAD(ρ_init),
        CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if TD.print_warning()
            @print("-----------------------------------------\n")
            @print("maxiter reached in air_density_equil:\n")
            @print("    Method=NewtonsMethodAD")
            @print(", p=", p)
            @print(", T=", T)
            @print(", q_tot=", q_tot)
            @print(", ρ=", sol.root)
            @print(", maxiter=", maxiter)
            @print(", tol=", tol.tol, "\n")
        end
        if TD.error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end

"""
    PhaseEquil_psq(param_set, p, s, q_tot)

Constructs a [`PhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` pressure
 - `s` entropy
 - `q_tot` total specific humidity
 - `temperature_tol` temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
"""
function PhaseEquil_psq(
    param_set::APS,
    p::FT,
    s::FT,
    q_tot::FT,
    maxiter::Int = 50,
    temperature_tol::FT = FT(1e-8),
    ::Type{sat_adjust_method} = RegulaFalsiMethod,
) where {FT <: Real, sat_adjust_method}
    phase_type = PhaseEquil
    T = saturation_adjustment_psq(sat_adjust_method, param_set, p, s, q_tot, phase_type, maxiter, temperature_tol)
    ρ = air_density_equil(param_set, p, T, q_tot)
    q = PhasePartition_equil(param_set, T, ρ, q_tot, phase_type)
    e_int = internal_energy(param_set, T, q)
    return PhaseEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q.tot, T)
end
