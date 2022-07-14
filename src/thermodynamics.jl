import RootSolvers
RS = RootSolvers

import KernelAbstractions
KA = KernelAbstractions

CPP = CLIMAParameters.Planet
const ITERTYPE = Union{Int, Nothing}
TOLTYPE(FT) = Union{FT, Nothing}

function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT) where {FT}
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    return TD.PhaseEquil_pθq(param_set, p, θ_liq_ice, q_tot, config...)
end
function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    config = ()
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    return TD.PhaseNonEquil_pθq(param_set, p, θ_liq_ice, q, config...)
end

function geopotential(param_set, z::FT) where {FT}
    grav = FT(CPP.grav(param_set))
    return grav * z
end

function enthalpy(mse::FT, e_pot::FT) where {FT}
    return mse - e_pot
end

#####
##### Thermodynamic variable inputs: p, h, q_tot
#####

function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT, NM <: RS.NewtonsMethodAD, phase_type <: TD.PhaseEquil}
    T_min::FT = CPP.T_min(param_set)
    T_init = max(T_min, air_temperature_from_enthalpy(param_set, h, TD.PhasePartition(q_tot))) # Assume all vapor
    return RS.NewtonsMethodAD(T_init)
end

function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT, NM <: RS.SecantMethod, phase_type <: TD.PhaseEquil}
    T_min::FT = CPP.T_min(param_set)
    q_pt = TD.PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_enthalpy(param_set, h, q_pt)
    T_1 = max(T_min, air_temperature_from_enthalpy(param_set, h, TD.PhasePartition(q_tot))) # Assume all vapor
    T_2 = TD.bound_upper_temperature(T_1, T_2)
    return RS.SecantMethod(T_1, T_2)
end

function sa_numerical_method_phq(
    ::Type{NM},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT, NM <: RS.RegulaFalsiMethod, phase_type <: TD.PhaseEquil}
    T_min::FT = CPP.T_min(param_set)
    q_pt = TD.PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_enthalpy(param_set, h, q_pt)
    T_1 = max(T_min, air_temperature_from_enthalpy(param_set, h, TD.PhasePartition(q_tot))) # Assume all vapor
    T_2 = TD.bound_upper_temperature(T_1, T_2)
    return RS.RegulaFalsiMethod(T_1, T_2)
end

"""
    air_temperature_from_enthalpy(param_set, h, q::PhasePartition)
"""
function air_temperature_from_enthalpy(param_set::APS, h::FT, q::TD.PhasePartition{FT} = q_pt_0(FT)) where {FT <: Real}
    cp_m_ = TD.cp_m(param_set, q)
    T_0::FT = CPP.T_0(param_set)
    R_d::FT = CPP.R_d(param_set)
    LH_v0::FT = CPP.LH_v0(param_set)
    LH_f0::FT = CPP.LH_f0(param_set)
    e_int_i0::FT = CPP.e_int_i0(param_set)
    q_vap::FT = TD.vapor_specific_humidity(q)
    return (h + cp_m_ * T_0 - q_vap * LH_v0 + q.ice * LH_f0 - (1 - q.tot) * R_d * T_0) / cp_m_
end

function saturation_adjustment_given_phq(
    ::Type{sat_adjust_method},
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    temperature_tol::FT,
) where {FT <: Real, sat_adjust_method, phase_type <: TD.PhaseEquil}
    _T_min::FT = CPP.T_min(param_set)
    cp_d::FT = CPP.cp_d(param_set)
    # Convert temperature tolerance to a convergence criterion on internal energy residuals
    tol = RS.ResidualTolerance(temperature_tol * cp_d)

    T_1 = max(_T_min, air_temperature_from_enthalpy(param_set, h, TD.PhasePartition(q_tot))) # Assume all vapor
    ρ_T(T) = TD.air_density(param_set, T, p, TD.PhasePartition(q_tot))
    ρ_1 = ρ_T(T_1)
    q_v_sat = TD.q_vap_saturation(param_set, T_1, ρ_1, phase_type)
    unsaturated = q_tot <= q_v_sat
    if unsaturated && T_1 > _T_min
        return T_1
    end
    _T_freeze::FT = ICP.T_freeze(param_set)
    h_sat(T) = specific_enthalpy_sat(param_set, TD.heavisided(T), ρ_T(T), q_tot, phase_type)

    h_upper = h_sat(_T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    h_lower = h_sat(_T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if h_lower < h < h_upper
        return _T_freeze
    end
    sol = RS.find_zero(
        T -> h_sat(T) - h,
        sa_numerical_method_phq(sat_adjust_method, param_set, p, h, q_tot, phase_type),
        RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if print_warning()
            KA.@print("-----------------------------------------\n")
            KA.@print("maxiter reached in saturation_adjustment_phq:\n")
            TD.print_numerical_method(sat_adjust_method)
            KA.@print(", h=", h)
            KA.@print(", p=", p)
            KA.@print(", q_tot=", q_tot)
            KA.@print(", T=", sol.root)
            KA.@print(", maxiter=", maxiter)
            KA.@print(", tol=", tol.tol, "\n")
        end
        if TD.error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end

function specific_enthalpy(param_set::APS, T::FT, q::TD.PhasePartition{FT} = TD.q_pt_0(FT)) where {FT <: Real}
    R_m = TD.gas_constant_air(param_set, q)
    e_int = TD.internal_energy(param_set, T, q)
    return TD.specific_enthalpy(e_int, R_m, T)
end

function specific_enthalpy_sat(
    param_set::APS,
    T::FT,
    ρ::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT <: Real, phase_type <: TD.ThermodynamicState}
    return specific_enthalpy(param_set, T, TD.PhasePartition_equil(param_set, T, ρ, q_tot, phase_type))
end

function PhaseDry_ph(param_set::APS, p::FT, h::FT) where {FT <: Real}
    T = air_temperature_from_enthalpy(param_set, h)
    ρ = TD.air_density(param_set, T, p)
    e_int = TD.internal_energy(param_set, T)
    return TD.PhaseDry{FT}(e_int, ρ)
end

function PhaseEquil_phq(
    param_set::APS,
    p::FT,
    h::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    temperature_tol::FTT = nothing,
    ::Type{sat_adjust_method} = RS.SecantMethod,
) where {FT <: Real, sat_adjust_method, IT <: ITERTYPE, FTT <: TOLTYPE(FT)}
    maxiter === nothing && (maxiter = 40)
    temperature_tol === nothing && (temperature_tol = FT(1e-2))
    phase_type = TD.PhaseEquil{FT}
    q_tot_safe = TD.clamp(q_tot, FT(0), FT(1))
    T = saturation_adjustment_given_phq(
        sat_adjust_method,
        param_set,
        p,
        h,
        q_tot_safe,
        phase_type,
        maxiter,
        temperature_tol,
    )
    q_pt = TD.PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type)
    ρ = TD.air_density(param_set, T, p, q_pt)
    e_int = TD.internal_energy(param_set, T, q_pt)
    return TD.PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
end

function PhaseNonEquil_phq(param_set::APS, p::FT, h::FT, q_pt::TD.PhasePartition{FT}) where {FT <: Real}
    T = air_temperature_from_enthalpy(param_set, h, q_pt)
    ρ = TD.air_density(param_set, T, p, q_pt)
    e_int = TD.internal_energy(param_set, T, q_pt)
    return TD.PhaseNonEquil{FT}(e_int, ρ, q_pt)
end
