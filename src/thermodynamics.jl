function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT) where {FT}
    # config = (50, 1e-3, RootSolvers.RegulaFalsiMethod)
    # config = (50, 1e-3, RootSolvers.NewtonMethodAD)
    config = ()
    thermo_params = TCP.thermodynamics_params(param_set)
    return TD.PhaseEquil_pθq(thermo_params, p, θ_liq_ice, q_tot, config...)
end
function thermo_state_pθq(param_set::APS, p::FT, θ_liq_ice::FT, q_tot::FT, q_liq::FT, q_ice::FT) where {FT}
    config = ()
    # q_tot = max(q_tot, q_liq + q_ice + eps(FT)) # ensure that total specific humidity is at least the sum of liquid and ice
    q = TD.PhasePartition(q_tot, q_liq, q_ice)
    thermo_params = TCP.thermodynamics_params(param_set)
    return TD.PhaseNonEquil_pθq(thermo_params, p, θ_liq_ice, q, config...)
end

function geopotential(param_set, z::Real)
    FT = eltype(param_set)
    grav = FT(TCP.grav(param_set))
    return grav * z
end

function enthalpy(mse::FT, e_pot::FT) where {FT}
    return mse - e_pot
end



# --------------------------------------------------------------- #
function relative_humidity_over_phase(thermo_params::TD.APS, T::FT, p::FT, phase::TD.Phase, q::TD.PhasePartition{FT} = TD.q_pt_0(FT),
) where {FT <: Real}
    R_v::FT = TD.TP.R_v(thermo_params)
    q_vap = TD.vapor_specific_humidity(q)
    p_vap = q_vap * TD.air_density(thermo_params, T, p, q) * R_v * T
    p_vap_sat = TD.saturation_vapor_pressure(thermo_params, T, phase)
    return p_vap / p_vap_sat
end
relative_humidity_over_phase(param_set::APS, T::FT, p::FT, phase::TD.Phase) where {FT <: Real} = relative_humidity_over_phase(TCP.thermodynamics_params(param_set), T, p, phase)

function relative_humidity_over_phase(thermo_params::TD.APS, ts::TD.ThermodynamicState{FT}, phase::TD.Phase) where {FT <: Real}
    return relative_humidity_over_phase(thermo_params, TD.air_temperature(thermo_params, ts), TD.air_pressure(thermo_params, ts), phase, TD.PhasePartition(thermo_params, ts))
end
# --------------------------------------------------------------- #
"""
The default assumes over equilibrium phase, or if you provide a q, it could be over that q but it doesnt pass that q to saturation_vapor_pressure anyway...
"""
relative_humidity_over_ice(param_set::APS, ts::TD.ThermodynamicState{FT}) where {FT <: Real} = relative_humidity_over_ice(TCP.thermodynamics_params(param_set), ts)
relative_humidity_over_ice(thermo_params::TD.APS, ts::TD.ThermodynamicState{FT}) where {FT <: Real} = relative_humidity_over_phase(thermo_params, ts, TD.Ice())
relative_humidity_over_ice(param_set::APS, T::FT, p::FT, q::TD.PhasePartition{FT} = TD.q_pt_0(FT)) where {FT <: Real} = relative_humidity_over_phase(TCP.thermodynamics_params(param_set), T, p, TD.Ice(), q)

relative_humidity_over_liquid(param_set::APS, ts::TD.ThermodynamicState{FT}) where {FT <: Real} = relative_humidity_over_liquid(TCP.thermodynamics_params(param_set), ts)
relative_humidity_over_liquid(thermo_params::TD.APS, ts::TD.ThermodynamicState{FT}) where {FT <: Real} = relative_humidity_over_phase(thermo_params, ts, TD.Liquid())
relative_humidity_over_liquid(param_set::APS, T::FT, p::FT, q::TD.PhasePartition{FT} = TD.q_pt_0(FT)) where {FT <: Real} = relative_humidity_over_phase(TCP.thermodynamics_params(param_set), T, p, TD.Liquid(), q)
# --------------------------------------------------------------- #




"""
    If we want to set a liquid fraction but still do saturation adjustment,
    there's no trivial way to do it so we add our own

    E.g. for initialization with fixed liquid fraction
"""
function saturation_adjustment_given_pθq_and_liquid_fraction(
    ::Type{sat_adjust_method},
    param_set::TD.APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::FT,
    λ = one(FT), # fixed instead of default
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: TD.PhaseEquil}
    tol = TD.RS.RelativeSolutionTolerance(relative_temperature_tol)
    T_min::FT = TD.TP.T_min(param_set)
    T_freeze::FT = TD.TP.T_freeze(param_set)
    cp_d::FT = TD.TP.cp_d(param_set)
    cp_v::FT = TD.TP.cp_v(param_set)
    air_temp(q) = TD.air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    function θ_liq_ice_closure(T)
        q = TD.PhasePartition(oftype(T, 0))
        # λ = liquid_fraction(param_set, T, phase_type, q)
        q_pt = TD.PhasePartition_equil_given_p(
            param_set,
            T,
            oftype(T, p),
            oftype(T, q_tot),
            phase_type,
            λ, # use fixed liquid fraction
        )
        return TD.liquid_ice_pottemp_given_pressure(
            param_set,
            T,
            oftype(T, p),
            q_pt,
        )
    end
    q_vap_sat(T) =
        TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
    T_1 = max(T_min, air_temp(TD.PhasePartition(q_tot))) # Assume all vapor
    q_v_sat_1 = q_vap_sat(T_1)
    unsaturated = q_tot <= q_v_sat_1
    if unsaturated && T_1 > T_min
        return T_1
    end
    temperature_tol = T_freeze * relative_temperature_tol
    θ_liq_ice_upper = θ_liq_ice_closure(T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    θ_liq_ice_lower = θ_liq_ice_closure(T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if θ_liq_ice_lower < θ_liq_ice < θ_liq_ice_upper
        return T_freeze
    end
    roots(T) = oftype(T, θ_liq_ice) - θ_liq_ice_closure(T)
    sol = TD.RS.find_zero(
        roots,
        TD.sa_numerical_method_pθq(
            sat_adjust_method,
            param_set,
            p,
            θ_liq_ice,
            q_tot,
            phase_type,
            T_guess,
        ),
        TD.RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if TD.print_warning()
            TD.KA.@print("-----------------------------------------\n")
            TD.KA.@print("maxiter reached in saturation_adjustment_given_pθq:\n")
            TD.print_numerical_method(sat_adjust_method)
            TD.print_T_guess(sat_adjust_method, T_guess)
            TD.KA.@print(", p=", p)
            TD.KA.@print(", θ_liq_ice=", θ_liq_ice)
            TD.KA.@print(", q_tot=", q_tot)
            TD.KA.@print(", T=", sol.root)
            TD.KA.@print(", maxiter=", maxiter)
            TD.KA.@print(", tol=", tol.tol, "\n")
        end
        if TD.error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end


"""
    PhaseEquil_pθq(param_set, θ_liq_ice, q_tot[, maxiter, relative_temperature_tol, sat_adjust_method, T_guess])

Constructs a [`PhaseEquil`](@ref) thermodynamic state from:

 - `param_set` an `AbstractParameterSet`, see the [`Thermodynamics`](@ref) for more details
 - `p` air pressure
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `λ` liquid fraction
 - `relative_temperature_tol` relative temperature tolerance for saturation adjustment
 - `maxiter` maximum iterations for saturation adjustment
 - `sat_adjust_method` the numerical method to use.
 - `T_guess` initial guess for temperature in saturation adjustment
"""
PhaseEquil_pθq_given_liquid_fraction(
    param_set::TD.APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    λ::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::FTT = nothing,
    ::Type{sat_adjust_method} = TD.RS.SecantMethod,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, IT <: TD.ITERTYPE, FTT <: TD.TOLTYPE(FT), sat_adjust_method} =
    begin
        maxiter === nothing && (maxiter = 50)
        relative_temperature_tol === nothing &&
            (relative_temperature_tol = FT(1e-4))
        phase_type = TD.PhaseEquil{FT}
        q_tot_safe = clamp(q_tot, FT(0), FT(1))
        T = saturation_adjustment_given_pθq_and_liquid_fraction(
            sat_adjust_method,
            param_set,
            p,
            θ_liq_ice,
            q_tot_safe,
            phase_type,
            maxiter,
            relative_temperature_tol,
            λ,
            T_guess,
        )
        q_pt = TD.PhasePartition_equil_given_p(
            param_set,
            T,
            p,
            q_tot_safe,
            phase_type,
            λ,
        )
        ρ = TD.air_density(param_set, T, p, q_pt)
        e_int = TD.internal_energy(param_set, T, q_pt)
        return TD.PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T)
    end

function PhasePartition_given_liquid_fraction(param_set::TD.APS, ts::TD.AbstractPhaseEquil{FT}, λ::FT) where {FT <: Real}
    T = TD.air_temperature(param_set, ts)
    ρ = TD.air_density(param_set, ts)
    q_tot = TD.total_specific_humidity(param_set, ts)
    phase_type = typeof(ts)
    # λ = liquid_fraction(param_set, T, phase_type) # fraction of condensate that is liquid
    qpt0 = TD.PhasePartition(zero(FT))
    p_vap_sat = TD.saturation_vapor_pressure(param_set, phase_type, T, qpt0, λ)

    return TD.PhasePartition_equil(param_set, T, ρ, q_tot, p_vap_sat, λ)
end

liquid_specific_humidity_given_liquid_fraction(q::TD.PhasePartition, λ::FT) where {FT} = q.liq
liquid_specific_humidity_given_liquid_fraction(param_set::TD.APS, ts::TD.ThermodynamicState{FT}, λ::FT = one(FT)) where{FT} = PhasePartition_given_liquid_fraction(param_set, ts, λ).liq