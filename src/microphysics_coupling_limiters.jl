#=
These are methods to take timescales, sources, etc and limit them for our timestep to ensure that we don't violate any physical constraints and crash
=#

# ======================================================================================================================================================================================================================================================== #

@inline function resolve_S_S_addit(
    S_ql::FT,
    S_qi::FT,
    Δt_S::FT,
    S_ql_addit::FT,
    S_qi_addit::FT,
    Δt_S_addit::FT,
    Δt::FT,
) where {FT}
    # return S_ql * (Δt_S / Δt) + S_ql_addit * (Δt_S_addit / Δt), S_qi * (Δt_S / Δt) + S_qi_addit * (Δt_S_addit / Δt)
    return (S_ql * Δt_S) / Δt + (S_ql_addit * Δt_S_addit) / Δt, (S_qi * Δt_S) / Δt + (S_qi_addit * Δt_S_addit) / Δt # group this way to preserve the total tendency over the full Δt somce `/ Δt * Δt` will cancel, and `S_ql * Δt_S` is the orignal tendency over the substep Δt_S
end

@inline function resolve_S_S_addit(S_q::FT, Δt_S::FT, S_q_addit::FT, Δt_S_addit::FT, Δt::FT) where {FT}
    # return S_q * (Δt_S / Δt) + S_q_addit * (Δt_S_addit / Δt)
    return (S_q * Δt_S) / Δt + (S_q_addit * Δt_S_addit) / Δt # group this way to preserve the total tendency over the full Δt somce `/ Δt * Δt` will cancel, and `S_q * Δt_S` is the orignal tendency over the substep Δt_S
end

@inline function t_out_of_x(x::FT, dxdt::FT) where {FT} # helper func, how long does source S take to send x to 0
    # shortcircuit bc these will never change anything and cause division by 0 problems [ technically we could return to 0 but we use a fixed dxdt so we disallow that]
    if iszero(dxdt) || iszero(x)
        return FT(Inf)
    end
    # now we calcualte the time
    t = x / -dxdt
    return (t < zero(FT)) ? FT(Inf) : t
end






# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #


"""
No limiter means just return the sources as they are
"""
function calculate_timestep_limited_sources(
    moisture_sources_limiter::NoMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState,
) where {FT}

    S_ql = (q_vap - q_eq.liq) / τ_liq
    S_qi = (q_vap - q_eq.ice) / τ_ice

    return calculate_timestep_limited_sources(
        moisture_sources_limiter,
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        q_vap,
        dqvdt,
        dTdt,
        q,
        q_eq,
        Δt,
        ts,
        S_ql,
        S_qi,
    )
end

function calculate_timestep_limited_sources(
    moisture_sources_limiter::NoMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState,
    S_ql::FT,
    S_qi::FT,
) where {FT}

    if iszero(q.liq) && (S_ql < 0) # don't allow evap if there's already no liquid (edge case)
        S_ql = FT(0)
    end

    if iszero(q.ice) && (S_qi < 0)  # don't allow sublimation if there's already no ice (edge case)
        S_qi = FT(0)
    end

    return S_ql, S_qi
end

# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================= #

"""
The basic limiter only does a quick check against the timestep to ensure we don't consume more of our source than we have available.
However, note it's checking q_vap, not supersaturation. So this is remains appropriate for say relax_to_equilibrium but wont't save you from progressing to subsaturation.
"""
function calculate_timestep_limited_sources(
    moisture_sources_limiter::BasicMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState,
) where {FT}
    S_ql = (q_vap - q_eq.liq) / τ_liq # | microphys_params.τ_cond_evap | CMNe.conv_q_vap_to_q_liq_ice(microphys_params, liq_type, q_eq, TD.PhasePartition(FT(0),q_vap,FT(0)))  
    S_qi = (q_vap - q_eq.ice) / τ_ice # -(source to vapor) = source to condensate
    return calculate_timestep_limited_sources(
        moisture_sources_limiter,
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        q_vap,
        dqvdt,
        dTdt,
        q,
        q_eq,
        Δt,
        ts,
        S_ql,
        S_qi,
    )
end

function calculate_timestep_limited_sources(
    moisture_sources_limiter::BasicMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState,
    S_ql::FT,
    S_qi::FT,
) where {FT}

    # TODO - handle limiters elswhere (note, these break our allowance for liquid to ice compensation, no?)
    if S_ql >= FT(0)
        S_ql = min(S_ql, q_vap / Δt)
    else
        S_ql = -min(-S_ql, q.liq / Δt)
    end
    if S_qi >= FT(0)
        S_qi = min(S_qi, q_vap / Δt)
    else
        S_qi = -min(-S_qi, q.ice / Δt)
    end
    return S_ql, S_qi
end

# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================= #

"""
The morrison milbrandt calculations are a bit more complex... essentially they asymptote to a steady state.
They are not a magic bullet then -- as Δt increases, these can become diluted. If other sources are present, and say if Δt is changing, the ratio of this to other sources will continually change.
This complicates balancing and limiters.

NOTE: We could try to make a version of this that does not cutoff but it would have WBF going off both directions forever... This may not be desirable. Also it's time evolution is dependent on the evolution of q_liq, q_ice, q_vap, so composing with other sources isn't as simple as adding the tendencies together. Just use with caution.
"""
function calculate_timestep_limited_sources(
    moisture_sources_limiter::MorrisonMilbrandt2015MoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState;
    liq_fraction::FT = FT(1),
    ice_fraction::FT = FT(1),
    cld_fraction::FT = FT(1),
) where {FT}
    # this *shouldn't* need limiters the way it's defined but... lol we'll see...
    if iszero(Δt)
        return FT(0), FT(0)
    end
    S_ql, S_qi = morrison_milbrandt_2015_style(
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        q_vap,
        dqvdt,
        dTdt,
        q,
        q_eq,
        Δt,
        ts;
        opts = MM2015Opts{FT}(
            emit_warnings = false,
            liq_fraction = liq_fraction,
            ice_fraction = ice_fraction,
            cld_fraction = cld_fraction,
        ),
    )

    # if isnothing(S_ql)  # An underflow fallback, go to :standard
    #     return calculate_timestep_limited_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    # end

    if !isfinite(S_ql) || !isfinite(S_qi)
        error(
            "Morrison-Milbrandt 2015 source calculations failed. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts",
        )
    end
    return S_ql, S_qi
end

# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================= #


function calculate_timestep_limited_sources(
    moisture_sources_limiter::MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState;
    liq_fraction::FT = FT(1),
    ice_fraction::FT = FT(1),
    cld_fraction::FT = FT(1),
) where {FT}
    if iszero(Δt)
        return FT(0), FT(0)
    end

    if q.tot ≥ 1 # should never happen but MM2015 one can't handle it.
        return calculate_timestep_limited_sources(
            StandardSupersaturationMoistureSourcesLimiter(),
            param_set,
            area,
            ρ,
            p,
            T,
            w,
            τ_liq,
            τ_ice,
            dqvdt,
            dTdt,
            q_vap,
            q,
            q_eq,
            Δt,
            ts,
        )
    end

    S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        q_vap,
        dqvdt,
        dTdt,
        q,
        q_eq,
        Δt,
        ts;
        opts = MM2015EPAOpts{FT}(
            emit_warnings = false,
            fallback_to_standard_supersaturation_limiter = moisture_sources_limiter.fallback_to_standard_supersaturation_limiter,
            liq_fraction = liq_fraction,
            ice_fraction = ice_fraction,
            cld_fraction = cld_fraction,
        ),
    )

    if isnothing(S_ql)  # An underflow fallback, go to :standard [not implemented yet]
        return calculate_timestep_limited_sources(
            StandardSupersaturationMoistureSourcesLimiter(),
            param_set,
            area,
            ρ,
            p,
            T,
            w,
            τ_liq,
            τ_ice,
            q_vap,
            dqvdt,
            dTdt,
            q,
            q_eq,
            Δt,
            ts,
        )
    end

    if !isfinite(S_ql) || !isfinite(S_qi)
        error(
            "Morrison-Milbrandt 2015 (exponential part only) source calculations failed. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts",
        )
    end

    thermo_params = TCP.thermodynamics_params(param_set)
    dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(1), T, q_eq.liq)
    dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(0), T, q_eq.ice)

    if (
        (S_ql + S_qi) * Δt ≥ (q_vap + eps(FT) + max(dqvdt, FT(0)) * Δt + max(-dTdt * dqsl_dT, -dTdt * dqsi_dT, 0) * Δt)
    ) && (!iszero(q_vap)) # consume all vapor (should never happen, really we should allow for liq and ice to contribute but all vapor gone is ridiculous, qv is like 2 orders of magnitude larger than ql, qi)
        # if ((S_ql + S_qi)*Δt ≥ (q_vap+eps(FT))) && (!iszero(q_vap)) # consume all vapor (should never happen, really we should allow for liq and ice to contribute but all vapor gone is ridiculous, qv is like 2 orders of magnitude larger than ql, qi)
        @error(
            "Morrison-Milbrandt 2015 (exponential part only) source calculations returned a total source greater than the available vapor. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dqvdt = $dqvdt; dTdt = $dTdt; dqsl_dT = $dqsl_dT; dqsi_dT = $dqsi_dT"
        )
        # fallback to StandardSupersaturationMoistureSourcesLimiter?
        return calculate_timestep_limited_sources(
            StandardSupersaturationMoistureSourcesLimiter(),
            param_set,
            area,
            ρ,
            p,
            T,
            w,
            τ_liq,
            τ_ice,
            q_vap,
            dqvdt,
            dTdt,
            q,
            q_eq,
            Δt,
            ts,
        )
    elseif (S_ql + S_qi) * Δt < -(q.liq + q.ice + eps(FT)) # lose more than all liquid and ice to sublimation...
        error(
            "Morrison-Milbrandt 2015 (exponential part only) source calculations returned a total source less than the available liquid and ice. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dqvdt = $dqvdt; dTdt = $dTdt; dqsl_dT = $dqsl_dT; dqsi_dT = $dqsi_dT",
        )
    end

    # restrict (S_ql + S_qi) * Δt to be less than q_vap / 2 since it seems the nudging tendencies can still cause problems esp in the updraft
    if (S_tot_Δt = ((S_ql + S_qi) * Δt)) > (q_vap / FT(2))
        S_ql *= (q_vap / FT(2)) / S_tot_Δt # scale down to q_vap / 2
        S_qi *= (q_vap / FT(2)) / S_tot_Δt # scale down to q_vap / 2
    end

    return S_ql, S_qi
end


# ======================================================================================================================================================================================================================================================== #
# ======================================================================================================================================================================================================================================================== #
# ========================================================================================================================================================================================================================================================= #

"""
This is the most basic version of the supersaturation relaxation limiter. It just ensures that we don't consume more of our source than we have available, including not consuming more than the supersaturation we have.
The logic is more compilcated in WBF regimes and when long enough timesteps would cause regime transitions. Consider using the Morrison-Milbrandt source calculations to elide some of these issues.å
"""



"""
Instead of S_ql + S_qi = 0, here we look for balance S_ql + S_qi = dδdt_no_S, where dδdt_no_S is the tendency of δ without any sources, from e.g. radiation, etc
We ignore the dTdt and w terms here for simplicity, since they make the system super complicated and in the standard form the rates don't respond to changes.

Note that now, this eq point may not exist! (i.e. the outcome could be ngative for example, or undefined)
"""
@inline function get_qv_eq_point(
    q_eq::TD.PhasePartition,
    τ_liq::FT,
    τ_ice::FT;
    dδdt_no_S::FT = FT(0),
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
) where {FT}
    # Psychrometric corrections [[ I dont think we need these here bc we're in all vapor...]]
    # τ_liq *= Γ_l # convert to 
    # τ_ice *= Γ_i # convert to time

    if isinf(τ_liq)
        if isinf(τ_ice) # no sources, so we just return the equilibrium point
            return FT(Inf)
        else # now it's just for S_qi to balance dδdt_no_S --> [[ δ / τ_ice = -dδdt_no_S  ]] --> [[ δ = - τ_ice * dδdt_no_S ]]
            return τ_ice * dδdt_no_S + (q_eq.ice - q_eq.liq) + q_eq.liq
        end
    elseif isinf(τ_ice) # now it's just for S_ql to balance dδdt_no_S --> [[ δ / τ_liq = - dδdt_no_S  ]] --> [[ δ = - τ_liq * dδdt_no_S ]]
        return τ_liq * dδdt_no_S + q_eq.liq
    end

    # return (q_eq.liq * τ_ice + q_eq.ice * τ_liq + dδdt_no_S * τ_liq * τ_ice) / (τ_liq + τ_ice)

    #=
    Eq_point means -S_ql = S_qi.
    This means (eq_point - q_eq.ice)/τ_ice = (q_eq.liq - eq_point)/τ_liq  ==> eq_point = (q_eq.ice * τ_liq + q_eq.liq * τ_ice)/(τ_liq + τ_ice)
    see https://www.wolframalpha.com/input?i2d=true&i=Divide%5B%5C%2840%29x+-+Subscript%5Bq%2Ci%5D%5C%2841%29%2Ck%5D%3DDivide%5B%5C%2840%29Subscript%5Bq%2CL%5D-x%5C%2841%29%2Cc%5D&assumption=%22UnitClash%22+-%3E+%7B%22L%22%2C+%7B%22Liters%22%2C+%22dflt%22%7D%7D&assumption=%7B%22C%22%2C+%22c%22%7D+-%3E+%7B%22Variable%22%2C+%22dflt%22%7D&assumption=%7B%22C%22%2C+%22L%22%7D+-%3E+%7B%22Variable%22%7D&assumption=%7B%22C%22%2C+%22q%22%7D+-%3E+%7B%22Variable%22%2C+%22dflt%22%7D
    =#

    qv_eq =
        q_eq.liq * (τ_ice / (τ_liq + τ_ice)) +
        q_eq.ice * (τ_liq / (τ_liq + τ_ice)) +
        dδdt_no_S * (τ_liq * τ_ice) / (τ_liq + τ_ice) # Floating point safer form... [[ we don't combine bc first is x τ_ice, second is x τ_liq ]]


    if qv_eq == q_eq.liq # These are indistinguishable so we cannot have equilibrium
        error(
            "Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast, qv_eq = q_eq.liq = $qv_eq (they are indistinguishable at this precision). Consider limiting the timescales or upgrading to BigFloats.",
        )
    elseif qv_eq == q_eq.ice
        error(
            "Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast, qv_eq = q_eq.ice = $qv_eq (they are indistinguishable at this precision). Consider limiting the timescales or upgrading to BigFloats.",
        )
    end

    return qv_eq

end

# reimplement for floating point precision...
@inline function get_δ_eq_point(
    q_eq::TD.PhasePartition,
    τ_liq::FT,
    τ_ice::FT;
    dδdt_no_S::FT = FT(0),
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
) where {FT} # = get_qv_eq_point(q_eq, τ_liq, τ_ice; dδdt_no_S=dδdt_no_S, Γ_l=Γ_l, Γ_i=Γ_i) - q_eq.liq # this is the point where the vapor sources balance out, i.e. the supersaturation is zero
    # @debug "calling get_δ_eq_point with q_eq = $q_eq; τ_liq = $τ_liq; τ_ice = $τ_ice; dδdt_no_S = $dδdt_no_S; Γ_l = $Γ_l; Γ_i = $Γ_i"
    if isinf(τ_liq)
        if isinf(τ_ice) # no sources, so we just return the equilibrium point
            # @debug "here 1"
            return FT(Inf), FT(Inf)

            # note δ0 = δ0i + (q_eq.ice - q_eq.liq) --> δ0i = δ0 - (q_eq.ice - q_eq.liq) --> δ0
        else
            #= 
                We have [[ dδ/dt = dδ/dt_external - δ_0/τ_liq - δ_0i/τ_ice ]], so what we need is [[ dδ/dt_external = δ_0/τ_liq + δ_0i/τ_ice ]].
            =#
            # now it's just for S_qi to balance dδdt_no_S --> [[ δ0i / τ_ice = dδdt_no_S  ]] --> [[ δ0i =  τ_ice * dδdt_no_S --> δ0 = (- τ_ice * dδdt_no_S) + (δ0 - δ0i). ]] Additionally, (δ0 - δ0i) = (q_eq.ice - q_eq.liq)
            # @debug "here 2"
            δ_eq = τ_ice * dδdt_no_S + (q_eq.ice - q_eq.liq)
            δi_eq = τ_ice * dδdt_no_S
            return δ_eq, δi_eq
        end
    elseif isinf(τ_ice) # now it's just for S_ql to balance dδdt_no_S --> [[ δ / τ_liq = - dδdt_no_S  ]] --> [[ δ = - τ_liq * dδdt_no_S ]]
        # @debug "here 3"
        δ_eq = τ_liq * dδdt_no_S
        δi_eq = τ_liq * dδdt_no_S + (q_eq.liq - q_eq.ice) # δi_eq = δ_eq + (q_eq.liq - q_eq.ice)
        return δ_eq, δi_eq
    end

    # @debug "here 4"

    # return (q_eq.liq * τ_ice + q_eq.ice * τ_liq + dδdt_no_S * τ_liq * τ_ice) / (τ_liq + τ_ice) - q_eq.liq
    # return (q_eq.liq * (τ_ice/(τ_liq + τ_ice)) + q_eq.ice * (τ_liq/(τ_liq + τ_ice)) + dδdt_no_S * (τ_liq * τ_ice)/(τ_liq + τ_ice)) - q_eq.liq # Floating point safer form...

    # Simplify the above by taking [[ take (q_eq.liq * τ_ice)/(τ_liq + τ_ice)  - q_eq.liq (τ_liq + τ_ice)/(τ_liq + τ_ice)  = (-q_eq.liq * τ_liq)/(τ_liq + τ_ice)  ]] , same as https://www.wolframalpha.com/input?i2d=true&i=solve+R-Divide%5B%5C%2840%29x%5C%2841%29%2Ck%5D-Divide%5B%5C%2840%29x%2Bd%5C%2841%29%2Cc%5D+%3D+0+for+x
    # return -(q_eq.liq * τ_liq)/(τ_liq + τ_ice) + (q_eq.ice * τ_liq)/(τ_liq + τ_ice) + (τ_liq * τ_ice) * dδdt_no_S / (τ_liq + τ_ice) # Floating point safer form...
    δ_eq = (q_eq.ice - q_eq.liq) * (τ_liq / (τ_liq + τ_ice)) + dδdt_no_S * (τ_liq * τ_ice) / (τ_liq + τ_ice) # Floating point safer form... [[can combine first two terms since they're both on (τ_liq/(τ_liq + τ_ice)) now ]]

    # δi_eq = δ_eq + (q_eq.liq - q_eq.ice)
    # calculate ice equilibrium point separately since, if it's small to floating point, it may be indistinguishable from zero from liquid supersaturation's perspective once it's smaller than nextfloat(δ_eq) - δ_eq
    # we can try the other form, we add back the q_eq.liq and then instead subtract q_eq.ice, so now [[ take (q_eq.ice * τ_liq)/(τ_liq + τ_ice) - q_eq.ice (τ_liq + τ_ice)/(τ_liq + τ_ice)  = (-q_eq.ice * τ_ice)/(τ_liq + τ_ice)  ]]
    δi_eq = (q_eq.liq - q_eq.ice) * (τ_ice / (τ_liq + τ_ice)) + dδdt_no_S * (τ_liq * τ_ice) / (τ_liq + τ_ice) # Floating point safer form... [[can combine first two terms since they're both on (τ_ice/(τ_liq + τ_ice)) now ]] 


    # if δ_eq == (q_eq.ice - q_eq.liq) # These are indistinguishable so we cannot have equilibrium, δi_eq will be 0 
    if iszero(δi_eq)
        # error("Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast, δ_eq = (q_eq.ice - q_eq.liq) = $δ_eq (they are indistinguishable at this precision). Consider limiting the timescales or upgrading to BigFloats.")
        error(
            "Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast, δi_eq = 0 (they are indistinguishable at this precision). Consider limiting the timescales or upgrading to BigFloats.",
        )
    end
    if iszero(δ_eq)
        error(
            "Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast, δ_eq = 0 (they are indistinguishable at this precision). Consider limiting the timescales or upgrading to BigFloats.",
        )
    end

    return δ_eq, δi_eq

end


# """
# Better bc saves cpu cycles.
# """
# @inline function get_qv_and_δ_eq_point(q_eq::TD.PhasePartition, τ_liq::FT, τ_ice::FT; dδdt_no_S::FT = FT(0), Γ_l::FT = FT(1), Γ_i::FT = FT(1)) where {FT}
#     qv_eq_point = get_qv_eq_point(q_eq, τ_liq, τ_ice; dδdt_no_S=dδdt_no_S, Γ_l=Γ_l, Γ_i=Γ_i)
#     return qv_eq_point, qv_eq_point - q_eq.liq # this is the point where the vapor sources balance out, i.e. the supersaturation is zero
# end



"""
    If dδdt_is_full_tendency = true, that means that S_ql and S_qi are already contained in dqvdt, so we don't need to subtract them.
"""
function calculate_next_standard_milestone_time_given_milestones(
    S_ql::FT,
    S_qi::FT,
    dδdt::FT,
    δ_0::FT,
    δ_0i::FT,
    δ_eq::FT,
    δi_eq::FT,
    q_liq::FT,
    q_ice::FT;
    dδdt_is_full_tendency::Bool = false,
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
) where {FT}
    # by definition, δ_0, δ_0i have the same sign as S_ql, S_qi (or both are 0). either way is safe for t_out_of_x
    t_out_of_liq = t_out_of_x(q_liq, S_ql)
    t_out_of_ice = t_out_of_x(q_ice, S_qi)


    if !at_δ_eq_point
        if dδdt_is_full_tendency

            t_hit_liq_sat = t_out_of_x(δ_0, dδdt) # this is the time it takes to hit liquid saturation, i.e. when we have no more liquid left
            t_hit_ice_sat = t_out_of_x(δ_0i, dδdt) # this is the time it takes to hit ice saturation, i.e. when we have no more ice left

            if allow_δ_eq_point
                if is_same_supersaturation_regime(δ_eq, δ_0, δ_0i) # check if the candidate is in the same supersaturation regime as the initial conditions
                    t_hit_eq_point_liq = t_out_of_x(δ_0 - δ_eq, dδdt) # this is the time it takes to hit the equilibrium point, i.e. when we have no more vapor left
                    t_hit_eq_point_ice = t_out_of_x(δ_0i - δi_eq, dδdt) # this is the time it takes to hit the equilibrium point, i.e. when we have no more vapor left
                    t_hit_eq_point = min(t_hit_eq_point_liq, t_hit_eq_point_ice) # we take the minimum of the two, since this is mostly an underflow check
                else
                    t_hit_eq_point = FT(Inf) # we can't hit the equilibrium point if we're not on the same side of it
                end
            else
                t_hit_eq_point = FT(Inf) # we can't hit the equilibrium point if we're not allowed to
            end

        else
            # dδdt_no_S = dδdt
            t_hit_liq_sat = t_out_of_x(δ_0, dδdt - (S_ql * Γ_l + S_qi * Γ_i))
            t_hit_ice_sat = t_out_of_x(δ_0i, dδdt - (S_ql * Γ_l + S_qi * Γ_i))

            if allow_δ_eq_point
                if is_same_supersaturation_regime(δ_eq, δ_0, δ_0i) # same sign means same regime, if not in the same regime we can't hit the equilibrium point first
                    t_hit_eq_point_liq = t_out_of_x(δ_0 - δ_eq, dδdt - (S_ql * Γ_l + S_qi * Γ_i)) # can happen both above and below freezing, its just the point where gains of one phase are equal to losses of the other
                    t_hit_eq_point_ice = t_out_of_x(δ_0i - δi_eq, dδdt - (S_ql * Γ_l + S_qi * Γ_i)) # can happen both above and below freezing, its just the point where gains of one phase are equal to losses of the other
                    t_hit_eq_point = min(t_hit_eq_point_liq, t_hit_eq_point_ice) # we take the minimum of the two, since this is mostly an underflow check
                else
                    t_hit_eq_point = FT(Inf) # we can't hit the equilibrium point if we're not on the same side of it
                end
            else
                t_hit_eq_point = FT(Inf) # we can't hit the equilibrium point if we're not allowed to
            end
        end
    else
        # if we're at the equilibrium point, we don't need to hit it, so we just set the time to 0
        t_hit_eq_point = FT(Inf) # we can't hit the equilibrium point, we're already there
        t_hit_liq_sat = FT(Inf) # we can't hit liquid saturation if we're at the equilibrium point
        t_hit_ice_sat = FT(Inf) # we can't hit ice saturation if we're at the equilibrium point
    end

    # findmin goes left to right, so this is also priority if there's a tie. Eq point gets priority over saturation points bc hitting it stops oscillations
    min_t, i_min_t = findmin((t_out_of_liq, t_out_of_ice, t_hit_eq_point, t_hit_liq_sat, t_hit_ice_sat)) # no need for out of vap, you get there you're cooked anyway... # Type Stable ((should we use find_min_t?) or does not have the MM2015 problems)


    local milestone::MilestoneType
    if isinf(min_t)
        milestone = NotAtSupersaturationMilestone
    else
        milestone = if i_min_t == 1
            OutOfLiquidMilestone
        elseif i_min_t == 2
            OutOfIceMilestone
        elseif i_min_t == 3
            AtSupersaturationStationaryPointMilestone
        elseif i_min_t == 4
            AtSaturationOverLiquidMilestone
        elseif i_min_t == 5
            AtSaturationOverIceMilestone
        end
    end

    # @debug "min_t = $min_t; milestone = $milestone; i_min_t = $i_min_t; t_out_of_liq = $t_out_of_liq; t_out_of_ice = $t_out_of_ice; t_hit_eq_point = $t_hit_eq_point; t_hit_liq_sat = $t_hit_liq_sat; t_hit_ice_sat = $t_hit_ice_sat; S_ql = $S_ql; S_qi = $S_qi; dδdt = $dδdt; δ_0 = $δ_0; δ_0i = $δ_0i; δ_eq = $δ_eq; δi_eq = $δi_eq; q_liq = $q_liq; q_ice = $q_ice; dδdt_is_full_tendency = $dδdt_is_full_tendency; Γ_l = $Γ_l; Γ_i = $Γ_i; allow_δ_eq_point = $allow_δ_eq_point; at_δ_eq_point = $at_δ_eq_point"
    if milestone == AtSupersaturationStationaryPointMilestone
        # @debug "milestone = $milestone from inputs: S_ql = $S_ql; S_qi = $S_qi; dδdt = $dδdt; δ_0 = $δ_0; δ_0i = $δ_0i; δ_eq = $δ_eq; δi_eq = $δi_eq; q_liq = $q_liq; q_ice = $q_ice; dδdt_is_full_tendency = $dδdt_is_full_tendency; Γ_l = $Γ_l; Γ_i = $Γ_i"
    end

    return min_t, milestone::MilestoneType
end
calculate_next_standard_milestone_time_given_milestones(
    S_ql::FT,
    S_qi::FT,
    dδdt::FT,
    δ_0::FT,
    δ_0i::FT,
    δ_eq::FT,
    δi_eq::FT,
    q::TD.PhasePartition;
    dδdt_is_full_tendency::Bool = false,
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
) where {FT} = calculate_next_standard_milestone_time_given_milestones(
    S_ql,
    S_qi,
    dδdt,
    δ_0,
    δ_0i,
    δ_eq,
    δi_eq,
    q.liq,
    q.ice;
    dδdt_is_full_tendency = dδdt_is_full_tendency,
    Γ_l = Γ_l,
    Γ_i = Γ_i,
    at_δ_eq_point = at_δ_eq_point,
    allow_δ_eq_point = allow_δ_eq_point,
)
calculate_next_standard_milestone_time_given_milestones(
    S_ql::FT,
    S_qi::FT,
    dδdt::FT,
    q_vap::FT,
    δ_eq::FT,
    δi_eq::FT,
    q_liq::FT,
    q_ice::FT,
    q_eq::TD.PhasePartition;
    dδdt_is_full_tendency::Bool = false,
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
) where {FT} = calculate_next_standard_milestone_time_given_milestones(
    S_ql,
    S_qi,
    dδdt,
    q_vap - q_eq.liq,
    q_vap - q_eq.ice,
    δ_eq,
    δi_eq,
    q_liq,
    q_ice;
    dδdt_is_full_tendency = dδdt_is_full_tendency,
    Γ_l = Γ_l,
    Γ_i = Γ_i,
    at_δ_eq_point = at_δ_eq_point,
    allow_δ_eq_point = allow_δ_eq_point,
)
calculate_next_standard_milestone_time_given_milestones(
    S_ql::FT,
    S_qi::FT,
    dδdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    δ_eq::FT,
    δi_eq::FT;
    dδdt_is_full_tendency::Bool = false,
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
) where {FT} = calculate_next_standard_milestone_time_given_milestones(
    S_ql,
    S_qi,
    dδdt,
    TD.vapor_specific_humidity(q),
    δ_eq,
    δi_eq,
    q.liq,
    q.ice,
    q_eq;
    dδdt_is_full_tendency = dδdt_is_full_tendency,
    Γ_l = Γ_l,
    Γ_i = Γ_i,
    at_δ_eq_point = at_δ_eq_point,
    allow_δ_eq_point = allow_δ_eq_point,
)



function calculate_next_standard_milestone_time(
    q_eq::TD.PhasePartition,
    q_liq::FT,
    q_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    δ_eq::FT,
    δi_eq::FT,
    S_ql::FT,
    S_qi::FT;
    dδdt::FT = FT(0),
    dδdt_is_full_tendency::Bool = false,
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
) where {FT}
    # @debug "calculate_next_standard_milestone_time: δ_eq = $δ_eq; δi_eq = $δi_eq; q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i; dδdt = $dδdt; Γ_l = $Γ_l; Γ_i = $Γ_i; dδdt_is_full_tendency = $dδdt_is_full_tendency; S_ql = $S_ql; S_qi = $S_qi; at_δ_eq_point = $at_δ_eq_point; allow_δ_eq_point = $allow_δ_eq_point"
    min_t, milestone::MilestoneType = calculate_next_standard_milestone_time_given_milestones(
        S_ql,
        S_qi,
        dδdt,
        δ_0,
        δ_0i,
        δ_eq,
        δi_eq,
        q_liq,
        q_ice;
        Γ_l = Γ_l,
        Γ_i = Γ_i,
        dδdt_is_full_tendency = dδdt_is_full_tendency,
        at_δ_eq_point = at_δ_eq_point,
        allow_δ_eq_point = allow_δ_eq_point,
    ) # dδdt_is_full_tendency = true

    # TODO: We really should just keep track of dδdt everywhere, rather than not  returning it. that would make it easier to handle `at_δ_eq_point` and other cases with better continuity

    return min_t, milestone, S_ql, S_qi, δ_eq, δi_eq
end


"""

"""
function calculate_next_standard_milestone_time(
    regime::AbstractSaturationRegime,
    q_eq::TD.PhasePartition,
    q_liq::FT,
    q_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    below_freezing::Bool,
    τ_liq::FT,
    τ_ice::FT;
    dδdt_no_S::FT = FT(0),
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
) where {FT}

    # using regimes is better becaues it breaks the uncertainty at the saturation points, since we calculate δ_eq assuming that q_liq and q_ice contribute only if the exist. consider being at liquid saturation, if we're going down and have no liquid, then liquid does not contribute. but if we're going up, then liquid does contribute. so we need to know which regime we are going to be in. MM2015_EPA handled this with clear transitions, we will do something similar
    # otherwise, we should use `≥` to defer to growth... really it should depend on which way we're going though... that's why the `regime` model is so crucial

    if below_freezing
        if (regime isa Supersaturated) # δ_0 ≥ FT(0) # Supersaturated
            S_ql = δ_0 / (Γ_l * τ_liq) # Eq C1
            S_qi = δ_0i / (Γ_i * τ_ice) # Eq C2
        elseif (regime isa WBF) #  (δ_0i ≥ FT(0)) # WBF
            S_ql = (q_liq > FT(0)) ? δ_0 / (Γ_l * τ_liq) : FT(0) # Eq C1
            τ_liq = (q_liq > FT(0)) ? τ_liq : FT(Inf) # liq only contributes if it exists
            S_qi = δ_0i / (Γ_i * τ_ice) #
        elseif (regime isa Subsaturated) # subsaturated
            S_ql = (q_liq > FT(0)) ? δ_0 / (Γ_l * τ_liq) : FT(0) # Eq C1
            τ_liq = (q_liq > FT(0)) ? τ_liq : FT(Inf) # liq only contributes if it exists
            S_qi = (q_ice > FT(0)) ? δ_0i / (Γ_i * τ_ice) : FT(0) # Eq C2
            τ_ice = (q_ice > FT(0)) ? τ_ice : FT(Inf) # ice only contributes if it exists
        else
            error("Unknown regime: $regime. Expected one of Supersaturated, WBF, or Subsaturated.")
        end
    else # above freezing
        if (regime isa Supersaturated) # δ_0i ≥ FT(0) # Supersaturated
            S_ql = δ_0 / (Γ_l * τ_liq) # Eq C1
            S_qi = FT(0)
            τ_ice = FT(Inf) # ice can't contribute bc it can't grow above freezing
        elseif (regime isa WBF) # δ_0 ≥ FT(0) # WBF
            S_ql = δ_0 / (Γ_l * τ_liq) # Eq C1
            S_qi = (q_ice > FT(0)) ? δ_0i / (Γ_i * τ_ice) : FT(0) # Eq C2
            τ_ice = (q_ice > FT(0)) ? τ_ice : FT(Inf) # ice only contributes if it exists
        elseif (regime isa Subsaturated) # subsaturated
            S_ql = (q_liq > FT(0)) ? δ_0 / (Γ_l * τ_liq) : FT(0) # Eq C1
            τ_liq = (q_liq > FT(0)) ? τ_liq : FT(Inf) # liq only contributes if it exists
            S_qi = (q_ice > FT(0)) ? δ_0i / (Γ_i * τ_ice) : FT(0) # Eq C2
            τ_ice = (q_ice > FT(0)) ? τ_ice : FT(Inf) # ice only contributes if it exists
        else
            error("Unknown regime: $regime. Expected one of Supersaturated, WBF, or Subsaturated.")
        end
    end

    δ_eq, δi_eq = get_δ_eq_point(q_eq, τ_liq, τ_ice; dδdt_no_S = dδdt_no_S, Γ_l = Γ_l, Γ_i = Γ_i) # this is the point where the vapor sources balance out, i.e. the supersaturation is zero


    dδdt = dδdt_no_S - (S_ql * Γ_l + S_qi * Γ_i) # this is the tendency of the supersaturation, i.e. how fast it is changing without the sources

    # TODO: We really should just keep track of dδdt everywhere, rather than not  returning it. that would make it easier to handle `at_δ_eq_point` and other cases with better continuity

    return calculate_next_standard_milestone_time(
        q_eq,
        q_liq,
        q_ice,
        δ_0,
        δ_0i,
        δ_eq,
        δi_eq,
        S_ql,
        S_qi;
        dδdt = dδdt,
        dδdt_is_full_tendency = true,
        Γ_l = Γ_l,
        Γ_i = Γ_i,
        at_δ_eq_point = at_δ_eq_point,
        allow_δ_eq_point = allow_δ_eq_point,
    ) # dδdt_is_full_tendency is a full tendency here as constructed by definition since it started w/o S
end



"""
A clean fallback to :StandardSupersaturation when the projected timescale is very short.
"""

function get_new_saturation_regime_type_from_milestone(
    milestone::MilestoneType,
    regime::AbstractSaturationRegime,
    old_δ_0::FT,
    old_δ_0i::FT,
) where {FT}
    regime_type =
        if milestone ∈ (
            NotAtSupersaturationMilestone,
            OutOfLiquidMilestone,
            OutOfIceMilestone,
            AtSupersaturationStationaryPointMilestone,
        ) # For these, we should still be in the same regime...
            if regime isa Supersaturated
                Supersaturated
            elseif regime isa WBF
                WBF
            elseif regime isa Subsaturated
                Subsaturated
            else
                error("invalid regime type $regime") # for type stability, we need to remove `nothing` from outputs
            end
        elseif milestone == AtSaturationOverLiquidMilestone # liq sat, above freezing so that's going to Subsat
            if is_below_freezing(regime)
                if regime isa Supersaturated
                    WBF
                elseif regime isa WBF
                    Supersaturated
                elseif (regime isa Subsaturated)
                    if !iszero(old_δ_0i)
                        error(
                            "Below freezing, we shouldn't be able to hit liq sat from Subsaturated unless δ_0i was zero",
                        )
                    else
                        Supersaturated
                    end
                else
                    error("invalid regime type $regime")
                end
            else # above freezing
                if (regime isa Supersaturated)
                    if !iszero(old_δ_0i)
                        error(
                            "Above freezing, we shouldn't be able to hit liq sat from Supersaturated unless δ_0i was zero",
                        )
                    else
                        Subsaturated
                    end
                elseif regime isa WBF
                    Subsaturated
                elseif regime isa Subsaturated
                    WBF
                else
                    error("invalid regime type $regime")
                end
            end

        elseif milestone == AtSaturationOverIceMilestone # ice sat, above freezing so that's going to Supersat
            if is_below_freezing(regime)
                if (regime isa Supersaturated)
                    if !iszero(old_δ_0)
                        error(
                            "Above freezing, we shouldn't be able to hit ice sat from Supersaturated unless δ_0 was zero",
                        )
                    else
                        Subsaturated
                    end
                elseif regime isa WBF
                    Subsaturated
                elseif regime isa Subsaturated
                    WBF
                else
                    error("invalid regime type $regime")
                end
            else # above freezing
                if regime isa Supersaturated
                    WBF
                elseif regime isa WBF
                    Supersaturated
                elseif regime isa Subsaturated
                    if !iszero(old_δ_0)
                        error(
                            "Above freezing, we shouldn't be able to hit ice sat from Subsaturated unless δ_0 was zero",
                        )
                    else
                        Supersaturated
                    end
                else
                    error("invalid regime type $regime")
                end
            end
        else
            error(
                "milestone should be NotAtSupersaturationMilestone, OutOfLiquidMilestone, OutOfIceMilestone, AtSaturationOverLiquidMilestone, or AtSaturationOverIceMilestone, but got $milestone",
            ) # branch to remove `nothing` from outputs and keep type stability
        end

    return regime_type
end
get_new_saturation_regime_enum_type_from_milestone(
    milestone::MilestoneType,
    regime::AbstractSaturationRegime,
    old_δ_0::FT,
    old_δ_0i::FT,
) where {FT} =
    get_saturation_regime_enum_type(get_new_saturation_regime_type_from_milestone(milestone, regime, old_δ_0, old_δ_0i))


function step(
    regime::AbstractSaturationRegime,
    ::StandardSupersaturationMoistureSourcesLimiter,
    Δt::FT,
    q_liq::FT,
    q_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    δ_eq::FT,
    δi_eq::FT,
    q_eq::TD.PhasePartition,
    S_ql::FT,
    S_qi::FT,
    milestone::MilestoneType = NotAtSupersaturationMilestone;
    at_δ_eq_point::Bool = false,
    dδdt_no_S::FT = FT(0),
    Γ_l = FT(1),
    Γ_i = FT(1),
) where {FT}
    isapprox(δ_0 - δ_0i, q_eq.ice - q_eq.liq, atol = 1e-6) || error(
        "δ_0 - δ_0i is not approximately equal to q_eq.ice - q_eq.liq, got δ_0 = $δ_0, δ_0i = $δ_0i, q_eq.ice = $(q_eq.ice), q_eq.liq = $(q_eq.liq)",
    )
    # @debug("starting with: q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i; δ_eq = $δ_eq; δi_eq = $δi_eq; S_ql = $S_ql; S_qi = $S_qi; milestone = $milestone; dδdt_no_S = $dδdt_no_S; Γ_l = $Γ_l; Γ_i = $Γ_i")
    # @debug("starting with (δ_0 - δ_0i) = $(δ_0 - δ_0i); (q_eq.ice - q_eq.liq) = $(q_eq.ice - q_eq.liq)")


    # We could add a check here that milestone != last_milestone, but it shouldn't be necessary if the code has no errors. this is an internal method so it shouldn't really be exposed to end users.

    old_δ_0 = δ_0
    old_δ_0i = δ_0i

    if milestone == NotAtSupersaturationMilestone # not at a milestone or not provided, just step as far as Δt says. (not floating point safe)
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        if (!at_δ_eq_point)
            dδ = -(S_ql * Γ_l * Δt) - (S_qi * Γ_i * Δt)
            δ_0 += dδ + (dδdt_no_S * Δt)
            δ_0i += dδ + (dδdt_no_S * Δt)
        end
    elseif milestone == OutOfLiquidMilestone #  out of liq
        if !at_δ_eq_point
            dδ = +q_liq * Γ_l
            δ_0 += dδ
            δ_0i += dδ
        end
        q_liq = FT(0)
        q_ice += S_qi * Δt

        if !at_δ_eq_point # floating point safe, this sum should be zero but it might not be (to floating point accuracy) since we reconstitute [[ δ_0 =  τ_ice * dδdt_no_S + (q_eq.ice - q_eq.liq) ]], and if [[ τ_ice * dδdt_no_S  << (q_eq.ice - q_eq.liq) ]], τ_ice * dδdt_no_S can be unrecoverable given it is not stored, meaning the new dδdt can differ from 0 by the same floating point error
            dδ = -(S_qi * Γ_i * Δt)
            δ_0 += dδ + (dδdt_no_S * Δt)
            δ_0i += dδ + (dδdt_no_S * Δt)
        end
    elseif milestone == OutOfIceMilestone # out of ice
        q_liq += S_ql * Δt
        if !at_δ_eq_point
            dδ = -(S_ql * Γ_l * Δt)
            δ_0 += dδ + (dδdt_no_S * Δt)
            δ_0i += dδ + (dδdt_no_S * Δt)
            dδ = +q_ice * Γ_i
            δ_0 += dδ
            δ_0i += dδ
        end
        q_ice = FT(0)
    elseif milestone == AtSupersaturationStationaryPointMilestone # hit eq point
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        δ_0i = δi_eq # pass in both explicitly to avoid underflow
        δ_0 = δ_eq
    elseif milestone == AtSaturationOverLiquidMilestone # hit liq sat
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        δ_0i = (δ_0i - δ_0) #+ q_eq.liq
        δ_0 = FT(0)
    elseif milestone == AtSaturationOverIceMilestone # hit ice sat
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        δ_0 = (δ_0 - δ_0i) #+ q_eq.ice
        δ_0i = FT(0)
    else
        error(
            "Unknown milestone type: $milestone. Expected one of NotAtSupersaturationMilestone, OutOfLiquidMilestone, OutOfIceMilestone, AtSaturationOverLiquidMilestone, or AtSaturationOverIceMilestone.",
        )
    end



    # @debug("dδdt = $(dδdt_no_S - S_ql * Γ_l - S_qi * Γ_i); dδdt*Δt = $((dδdt_no_S - S_ql * Γ_l - S_qi * Γ_i) * Δt); at_δ_eq_point = $at_δ_eq_point")

    isapprox(δ_0 - δ_0i, q_eq.ice - q_eq.liq, atol = 1e-6) || error(
        "δ_0 - δ_0i is not approximately equal to q_eq.ice - q_eq.liq; got δ_0 = $δ_0; δ_0i = $δ_0i; q_eq.ice = $(q_eq.ice); q_eq.liq = $(q_eq.liq)",
    )
    # @debug("ending with: q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i; δ_eq = $δ_eq; δi_eq = $δi_eq; S_ql = $S_ql; S_qi = $S_qi; milestone = $milestone; dδdt_no_S = $dδdt_no_S; Γ_l = $Γ_l; Γ_i = $Γ_i")
    # @debug("ending with (δ_0 - δ_0i) = $(δ_0 - δ_0i); (q_eq.ice - q_eq.liq) = $(q_eq.ice - q_eq.liq)")

    # new_regime_type = get_new_saturation_regime_type_from_milestone(milestone, regime, old_δ_0, old_δ_0i) # old_δ_0 and old_δ_0i are just for checks... could prolly get rid of them...
    # new_regime::new_regime_type = add_regime_parameters(new_regime_type, q_liq, q_ice, is_below_freezing(regime)) # This takes 2.436 μs, just simply doing new_regime_type{q_liq>0, q_ice>0, is_below_freezing(regime)}() would have taken 193 ns... so this is 10x slower

    # new_regime::SaturationRegimeEnumTypes = get_saturation_regime_enum_type(new_regime_type)
    new_regime_enum_type = get_new_saturation_regime_enum_type_from_milestone(milestone, regime, old_δ_0, old_δ_0i)
    return q_liq, q_ice, δ_0, δ_0i, new_regime_enum_type # new_regime_type{q_liq>0, q_ice>0, is_below_freezing(regime)}()
end


function calculate_timestep_limited_sources(
    moisture_sources_limiter::StandardSupersaturationMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState,
) where {FT}
    δ_0 = q_vap - q_eq.liq
    δ_0i = q_vap - q_eq.ice

    use_fix = true # same as MM2015_EPA default, but this might be bad at high T, maybe set to false later

    return standard_supersaturation_sources(
        moisture_sources_limiter,
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        δ_0,
        δ_0i,
        dqvdt,
        dTdt,
        q,
        q_eq,
        Δt,
        ts;
        use_fix = use_fix,
    )
end

function standard_supersaturation_sources(
    moisture_sources_limiter::StandardSupersaturationMoistureSourcesLimiter,
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    dqvdt::FT,
    dTdt::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState;
    use_fix::Bool = true,
) where {FT}


    (;
        g,
        L_i,
        L_l,
        c_p,
        e_sl,
        e_si,
        dqsl_dT,
        dqsi_dT,
        q_sl,
        q_si,
        q_liq,
        q_ice,
        T_freeze,
        δ_0,
        δ_0i,
        Γ_l,
        Γ_i,
        dqvdt,
        dTdt,
    ) = get_params_and_go_to_mixing_ratio_exponential_part_only(
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        δ_0,
        δ_0i,
        dqvdt,
        dTdt,
        q,
        q_eq,
        Δt,
        ts;
        use_fix = use_fix,
    )
    dδdt_no_S = A_c_func_no_WBF_EPA(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)  # Eq C4 no WBF

    S_ql = FT(0)
    S_qi = FT(0)

    Δt_left = Δt
    dt = FT(0)
    last_milestone = NotAtSupersaturationMilestone # this is the last milestone we hit, so we can use it to determine the next one
    depth = 0

    #= since the while loop handles the transitions, we just calculate regime internally instead of dispatching
        if we have problems with choosing regime because we start at supersaturation either way, then idk...
        i feel like it's very unlikelyy to happen, and even less likely to cause a problem unless timescales are unreasonably fast.
        when we use the fallback in MM2015_EPA, we pass the regime and only call `step()` so I think it should be fine

        If push came to shove we could calculate the initial dδdt and then use that to break the tie and determine the regime.
    =#



    below_freezing::Bool = T < TCP.T_triple(param_set) # In the code this is actually where the saturation vapor pressures are equal... I think it's wrong though

    # Get our initial regime, being careful to ensure we're going the right direction initially.
    if iszero(δ_0) || iszero(δ_0i) # possible but unlikely
        if iszero(δ_0)
            last_milestone = AtSaturationOverLiquidMilestone # if δ_0 is zero, we are at liquid saturation, so we start there
        elseif iszero(δ_0i)
            last_milestone = AtSaturationOverIceMilestone # if δ_0i is zero, we are at ice saturation, so we start there
        end
        dδdt_0 = get_dδdt_0(δ_0, δ_0i, q.liq, q.ice, τ_liq, τ_ice, dδdt_no_S, below_freezing)
        regime = get_saturation_regime(δ_0, δ_0i, q.liq, q.ice, below_freezing; dδdt = dδdt_0) # use this version to break ties and make sure we start going the right direction.
    else
        regime = get_saturation_regime(δ_0, δ_0i, q.liq, q.ice, below_freezing)
    end



    while (Δt_left > FT(0)) && (depth <= 10 + 5)
        # @debug("------------------- q_liq > 0 = $((q_liq > FT(0))) ($q_liq); q_ice > 0 = $((q_ice > FT(0))) ($q_ice); δ_0 > 0 $(δ_0 > 0) ($δ_0); δ_0i > 0  $(δ_0i > 0) ($δ_0i); -----------------------\n\n")
        milestone_t, milestone, S_ql_addit, S_qi_addit, δ_eq, δi_eq = calculate_next_standard_milestone_time(
            regime,
            q_eq,
            q_liq,
            q_ice,
            δ_0,
            δ_0i,
            T < T_freeze,
            τ_liq,
            τ_ice;
            dδdt_no_S = dδdt_no_S,
            Γ_l = Γ_l,
            Γ_i = Γ_i,
            at_δ_eq_point = (last_milestone == AtSupersaturationStationaryPointMilestone),
        ) # this is the time to the next milestone, i.e. the time to hit the next source limit
        dt_here = min(Δt_left, milestone_t) # this is the time to the next milestone, i.e. the time to hit the next source limit

        # S_ql_addit, S_qi_addit = clamp_S(S_ql_addit, S_qi_addit, regime, δ_0, δ_0i, q_liq, q_ice, dt_here, dδdt_no_S) # for safety... i think this works, clamp to dt_here... [[ i think this is auto clamped by construction though... idk]]



        # sometimes the δ_eq we get out here is not the same as the one input, depending on if we have q_liq, q_ice... so we double check
        if last_milestone == AtSupersaturationStationaryPointMilestone
            if (δ_eq == δ_0) || (δ_eq == δ_0i) # if we hit the eq point, then we are at the stationary point
                at_δ_eq_point = true
            else
                at_δ_eq_point = false # the eq point seems like it can change based on the boundary based on if we calculate
            end
        else
            at_δ_eq_point = false
        end

        q_liq, q_ice, δ_0, δ_0i, new_regime_enum_type = step(
            regime,
            moisture_sources_limiter,
            dt_here,
            q_liq,
            q_ice,
            δ_0,
            δ_0i,
            δ_eq,
            δi_eq,
            q_eq,
            S_ql_addit,
            S_qi_addit,
            ((milestone_t < Δt_left) ? milestone : NotAtSupersaturationMilestone);
            at_δ_eq_point = at_δ_eq_point,
            dδdt_no_S = dδdt_no_S,
            Γ_l = Γ_l,
            Γ_i = Γ_i,
        ) # this is the step to the next milestone, i.e. the time to hit the next source limit
        new_regime_type = get_saturation_regime_type(Val(new_regime_enum_type))
        regime = new_regime_type{q_liq > 0, q_ice > 0, below_freezing}() # reconstruct regime with new q_liq and q_ice

        S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, dt, S_ql_addit, S_qi_addit, dt_here, dt + dt_here)

        # @debug "dt_here = $dt_here; milestone_t = $milestone_t; milestone = $milestone; last_milestone = $last_milestone; Δt_left = $Δt_left; Δt = $Δt; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; δ_eq = $δ_eq; δi_eq = $δi_eq; q_liq = $q_liq; q_ice = $q_ice; δ_0 = $δ_0; δ_0i = $δ_0i; Γ_l = $Γ_l; Γ_i = $Γ_i; dδdt_no_S = $dδdt_no_S; S_ql = $S_ql; S_qi = $S_qi; q_eq = $q_eq; T = $T; T_freeze = $T_freeze; p = $p; ρ = $ρ; area = $area; ts = $ts; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; dTdt = $dTdt; dqvdt = $dqvdt"

        dt += dt_here # add the time to the next milestone to the total time
        Δt_left -= dt_here # reduce the time left by the time we just stepped

        if !(milestone == NotAtSupersaturationMilestone) && (milestone == last_milestone) # if we hit the same milestone twice, we have a problem
            error(
                "Hit the same milestone twice! milestone = $milestone, last_milestone = $last_milestone; dt_here = $dt_here; Δt_left = $Δt_left; S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq, q_ice = $q_ice; δ_eq = $δ_eq; δi_eq = $δi_eq; Γ_l = $Γ_l; Γ_i = $Γ_i; dδdt_no_S = $dδdt_no_S; q_eq = $q_eq; T = $T; p = $p; ρ = $ρ; area = $area; ts = $ts; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; dTdt = $dTdt; dqvdt = $dqvdt",
            )
        end

        last_milestone = milestone # update the last milestone time


        if depth == 10
            @error (
                "Hit the maximum depth of 10 iterations, this is likely a bug in the code, please report it Or your τ values could just be very small, try upgrading to BigFloat or limiting them: milestone = $milestone; last_milestone = $last_milestone; dt = $dt; dt_here = $dt_here; Δt_left = $Δt_left; S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq, q_ice = $q_ice; δ_eq = $δ_eq; δi_eq = $δi_eq; Γ_l = $Γ_l; Γ_i = $Γ_i; dδdt_no_S = $dδdt_no_S; q_eq = $q_eq; T = $T; p = $p; ρ = $ρ; area = $area; ts = $ts; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; dTdt = $dTdt; dqvdt = $dqvdt"
            )
            dTdt = FT(0) # reset the temperature tendency to zero
            dqvdt = FT(0) # reset the vapor tendency to zero
            w = FT(0) # reset the vertical velocity to zero
            dδdt_no_S = FT(0) # reset the supersaturation tendency to zero
        # @debug "setting to 0 and. continuing for 5 iterations"

        elseif depth > (10 + 5)
            error(
                "Hit the maximum depth of 15 iterations, this is likely a bug in the code, please report it! Or your τ values could just be very small, try upgrading to BigFloat or limiting them: milestone = $milestone; last_milestone = $last_milestone; dt = $dt; dt_here = $dt_here; Δt_left = $Δt_left; S_ql = $S_ql; S_qi = $S_qi; δ_0 = $δ_0; δ_0i = $δ_0i; q_liq = $q_liq, q_ice = $q_ice; δ_eq = $δ_eq; δi_eq = $δi_eq; Γ_l = $Γ_l; Γ_i = $Γ_i; dδdt_no_S = $dδdt_no_S; q_eq = $q_eq; T = $T; p = $p; ρ = $ρ; area = $area; ts = $ts; S_ql_addit = $S_ql_addit; S_qi_addit = $S_qi_addit; dTdt = $dTdt; dqvdt = $dqvdt",
            )
        end

        depth += 1 # increase the depth of the loop

    end

    return S_ql, S_qi

end


# [[ Deprecated the old version of StandardSupersaturationLimiter -- we deprecated this version as it can't handle external forcings from w, dqvdt, dTdt, etc. ]]
