"""
Code for handling supersaturation regimes
"""

# === Regime Types  ======================================================================================================================================================================================== #


abstract type AbstractSaturationRegime end

struct Supersaturated{HL, HI, BF} <: AbstractSaturationRegime
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool
end
Supersaturated(has_liquid::Bool, has_ice::Bool, below_freezing::Bool) = Supersaturated{has_liquid, has_ice, below_freezing}(has_liquid, has_ice, below_freezing)

struct WBF{HL, HI, BF} <: AbstractSaturationRegime
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool # allow this type in lieu of created another type for when we're in WBF regime but not below freezing
end
WBF(has_liquid::Bool, has_ice::Bool, below_freezing::Bool) = WBF{has_liquid, has_ice, below_freezing}(has_liquid, has_ice, below_freezing)
WBF(has_liquid::Bool, has_ice::Bool) = WBF{has_liquid, has_ice, true}(has_liquid, has_ice, true)

# [[ Deprecated :: This is problematic, because you really should never be actually on th line. When you are, it's hard to write what the rate of evap from one phase should be since the slower one sets rate controls. We've replaced with just calculating the eq point ]]
# struct LowerSatLine{HL, HI, BF} <: AbstractSaturationRegime # you can end up here if you've reached ice sat but liq effective timescale is slower than ices so you're stuck.
#     has_liquid::Bool
#     has_ice::Bool
#     below_freezing::Bool # below freezing you're stuck on the ice sat line, above freezing you're stuck on the liq sat line but not really bc there's no ice formation
# end

struct Subsaturated{HL, HI, BF}  <: AbstractSaturationRegime
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool # does this really matter for subsaturated? (using it bc whether q_liq_sat or q_ice_sat is lower depends on it)
end
Subsaturated(has_liquid::Bool, has_ice::Bool, below_freezing::Bool) = Subsaturated{has_liquid, has_ice, below_freezing}(has_liquid, has_ice, below_freezing)


# -=-=- Milestones -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
abstract type AbstractSupersaturationLimiterMilestone end
struct NotAtSupersaturationMilestone <: AbstractSupersaturationLimiterMilestone end
struct OutOfLiquid <: AbstractSupersaturationLimiterMilestone end
struct OutOfIce <: AbstractSupersaturationLimiterMilestone end
struct AtSaturationOverLiquid <: AbstractSupersaturationLimiterMilestone end
struct AtSaturationOverIce <: AbstractSupersaturationLimiterMilestone end
struct AtSupersaturationStationaryPoint <: AbstractSupersaturationLimiterMilestone end
struct OutofVapor <: AbstractSupersaturationLimiterMilestone end


# ==== Methods ============================================================================================================================================================================================== #

δi_from_δ(δ::FT, q_sl::FT, q_si::FT) where {FT} = δ + (q_sl - q_si)
δ_from_δi(δi::FT, q_sl::FT, q_si::FT) where {FT} = δi - (q_sl - q_si)

# ----------- #

"""
For getting the initial regime, we err to growth, unless dδ/dt suggests otherwise (i.e. it is explicitly negative)
"""
function get_regime_type(δ::FT, δi::FT, BF::Bool; dδdt::FT = FT(0) ) where {FT} # 2 FT 1 Bool
    if BF
        if δ < FT(0)
            if δi > FT(0) # ≥ bc we are usually heading towards WBF so ties go to WBF
                return WBF
            elseif iszero(δi) # if δi is zero, we are on the WBF/Subsaturated boundary, so we should be WBF
                return (dδdt < FT(0)) ? Subsaturated : WBF # if dδ/dt < 0, we are heading towards subsaturation, so we should be there. otherwise we err to WBF (dδdt == 0 shouldn't cause problems since then nothing is changing for the condensate at saturation)
            else
                return Subsaturated
            end
        else
            if iszero(δ) # if δ is zero, we are on the Supersaturated/WBF boundary, so we should be WBF
                return (dδdt < FT(0)) ? WBF : Supersaturated # if dδ/dt < 0, we are heading towards WBF, so we should be there, otherwise we err to Supersaturated (dδdt == 0 shouldn't cause problems since then nothing is changing for the condensate at saturation)
            end
            return Supersaturated
        end
    else
        if δi < FT(0)
            if δ > FT(0) # ≥ bc we are usually heading towards WBF so ties go to WBF
                return WBF
            elseif iszero(δ) # if δ is zero, we are on the WBF/Subsaturated boundary, so we should be WBF
                return (dδdt < FT(0)) ? Subsaturated : WBF # if dδ/dt < 0, we are heading towards subsaturation, so we should be there ( otherwise we err to WBF, dδdt == 0 shouldn't cause problems since then nothing is changing for the condensate at saturation)
            else
                return Subsaturated
            end
        else
            if iszero(δi) # if δi is zero, we are on the Supersaturated/WBF boundary, so we should be WBF
                return (dδdt < FT(0)) ? WBF : Supersaturated # if dδ/dt < 0, we are heading towards WBF, so we should be there (otherwise we err to Supersaturated, dδdt == 0 shouldn't cause problems since then nothing is changing for the condensate at saturation)
            end
            return Supersaturated # If you don't use T_freeze = T_triple, between 273.15 and 273.16, you can get sent here in this code despite being subsaturated for liquid.... bad things ensue...
        end
    end
end
get_regime_type(δ::FT, δi::FT, T::FT, T_freeze::FT) where {FT} = get_regime_type(δ, δi, T < T_freeze) # 4 FT

# helper to use q_vap
get_regime_type(q_vap::FT, q_sl::FT, q_si::FT, T::FT, T_freeze::FT; dδdt::FT = FT(0)) where {FT} = get_regime_type(q_vap-q_sl, q_vap-q_si, T < T_freeze; dδdt = dδdt) # 5 FT

#=
 Note -- we really should be making a distinction between T_freeze and T_triple/T_0... the latter is actually where the supersaturations cross over...
 Idk what the easiest fix is... we've hardcoded assumptions about the order of q_vap_sat_liq and q_vap_sat_ice in each state, and what you can `hit` next...

 We could just swap T_freeze for T_triple/T_0 and just allow ice to form up to 273.16? might be easier than a full rewrite of all our assumptions...

Because T doesn't really impact the exponential_part_only, we could also round [T_freeze = 273.15] < T < [T_triple = T_0 = 273.16] up to 273.16 and pretend that's all fine and dandy....
That wouldn't solve everything tho bc we'll still have the wrong saturations...
=#




function get_regime(δ::FT, δi::FT, q_liq::FT, q_ice::FT, BF::Bool; dδdt::FT = FT(0)) where {FT} # 4 FT 1 Bool
    regime_type = get_regime_type(δ, δi, BF; dδdt = dδdt)
    return add_regime_parameters(regime_type, q_liq, q_ice, BF)
end

# Helpers to dispatch back to supersaturation defined methods
get_regime(q_vap::FT, q_liq::FT, q_ice::FT, q_sl::FT, q_si::FT, below_freezing::Bool; dδdt::FT = FT(0)) where {FT} = get_regime(q_vap - q_sl, q_vap - q_si, q_liq, q_ice, below_freezing; dδdt = dδdt) # 5 FT 1 Bool
get_regime(q_vap::FT, q_liq::FT, q_ice::FT, q_sl::FT, q_si::FT, T::FT, T_freeze::FT; dδdt::FT = FT(0)) where {FT} = get_regime(q_vap, q_liq, q_ice, q_sl, q_si, T < T_freeze; dδdt = dδdt) # 7 FT
get_regime(q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, below_freezing::Bool; dδdt::FT = FT(0)) where {FT} = get_regime(q_vap, q.liq, q.ice, q_eq.liq, q_eq.ice, below_freezing; dδdt = dδdt) # 1 FT 2 PhasePartition 1 Bool
# get_regime(q::TD.PhasePartition, q_eq::TD.PhasePartition, below_freezing::Bool; dδdt::FT = FT(0)) where {FT} = get_regime(TD.vapor_specific_humidity(q), q.liq, q.ice, q_eq.liq, q_eq.ice, below_freezing; dδdt=dδdt) # 1 FT 2 PhasePartition 1 Bool

get_regime(δ::FT, δi::FT, q_liq::FT, q_ice::FT, T::FT, T_freeze::FT; dδdt::FT = FT(0)) where {FT} = get_regime(δ, δi, q_liq, q_ice, T < T_freeze; dδdt = dδdt) # 6 FT [ can't use BF here bc already have a 5 FT 1 Bool method]


# ----------- #

add_regime_parameters(regime_type::Type{<:AbstractSaturationRegime}, q_liq::FT, q_ice::FT, T::FT, T_freeze::FT) where {FT} = (q_liq, q_ice, T<T_freeze, regime_type)
function add_regime_parameters(regime_type::Type{<:AbstractSaturationRegime}, q_liq::FT, q_ice::FT, BF::Bool) where {FT}
    has_liq = q_liq > FT(0)
    has_ice = q_ice > FT(0)
    regime_type{has_liq, has_ice, BF}(has_liq, has_ice, BF)
end
add_regime_parameters(regime_type::Type{<:AbstractSaturationRegime}, q_liq::FT, has_ice::Bool, BF::Bool) where {FT} = regime_type{q_liq > FT(0), has_ice, BF}(q_liq > FT(0), has_ice, BF) # convenience if you're sure about one
add_regime_parameters(regime_type::Type{<:AbstractSaturationRegime}, has_liq::Bool, q_ice::FT, BF::Bool) where {FT} = regime_type{has_liq, q_ice > FT(0), BF}(has_liq, q_ice > FT(0), BF) # convenience if you're sure about one


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

"""
Get the initial direction of dδdt. Helpful in the very unlikely case you're starting at saturation.
We used to just err to WBF but we can be more precise about it.

You can take the outcome of this and use it in get_regime() to ensure you start in the right regime get_regime_type() in the (very unlikely) case your initial call is exactly at saturation.
"""
function get_dδdt_0(δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, τ_liq::FT, τ_ice::FT, dδdt_no_S::FT, below_freezing::Bool) where {FT}

    if below_freezing
        if δ_0 ≥ FT(0) # Supersaturated or borderline (borderline will return 0 anyway)
            S_ql_δ = (δ_0 / τ_liq)
            S_qi_δ = (δ_0i / τ_ice)
        elseif δ_0i ≥ FT(0) # WBF or borderline (borderline will return 0 anyway)
            S_ql_δ = (q_liq > FT(0)) ? (δ_0 / τ_liq) : FT(0)
            S_qi_δ = (δ_0i / τ_ice)
        else # subsaturated
            S_ql_δ = (q_liq > FT(0)) ? (δ_0 / τ_liq) : FT(0)
            S_qi_δ = (q_ice > FT(0)) ? (δ_0i / τ_ice) : FT(0)
        end
    else # above freezing
        if δ_0i ≥ FT(0) # Supersaturated or borderline ( borderline will return 0 anyway)
            S_ql_δ = δ_0 / τ_liq
            S_qi_δ = FT(0)
        elseif δ_0 ≥ FT(0) # WBF or borderline (borderline will return 0 anyway)
            S_ql_δ = δ_0 / τ_liq
            S_qi_δ = (q_ice > FT(0)) ? (δ_0i / τ_ice) : FT(0)
        else  # subsaturated
            S_ql_δ = (q_liq > FT(0)) ? (δ_0 / τ_liq) : FT(0)
            S_qi_δ = (q_ice > FT(0)) ? (δ_0i / τ_ice) : FT(0)
        end
    end

   return dδdt_no_S - (S_ql_δ + S_qi_δ)

end

# -------------------------------------------------------------------------------------------------------------------------- #


comparable_sign(x::FT, y::FT) where {FT} = (sign(x) == sign(y)) || iszero(x) || iszero(y)
function is_same_supersaturation_regime(δ_candidate::FT, δ_0::FT, δ_0i::FT) where {FT}
    """
    Check if the two supersaturations are in the same regime with respect to the equilibrium point.
    note, zeros are considered to be in the same regime, so we need to be careful with them.
    """
    return comparable_sign(δ_candidate, δ_0) && comparable_sign(δ_candidate + (δ_0i - δ_0), δ_0i) # this is the new solution, it handles zeros correctly and is more robust
end
