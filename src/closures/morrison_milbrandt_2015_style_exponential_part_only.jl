"""
# Milbrandt's equations are in mixing ratio not specific humidity, to convert we need to use

# q_spe = q_mix / (1 + q_mix)   [ in reality this is more like q_spe_c = q_mix_c / (1 + q_tot_mix) , bc 1+q_tot = total ratio all water / dry air...]
# q_mix = q_spe / (1 - q_spe)   [ in reality this is more like q_mix_c = q_spe_c / (1 - q_tot_spe) , bc 1-q_tot = dry mixing ratio...]

# [NOTE: q_tot appears instead of just q because we have more than one water species so we need to be careful...]

# note then d/dy (q_spe) = d/dy(q_mix / (1 + q_mix)) = d/dy(q_mix) / (1 + q_mix)^2  [ really it's more like d/dy(q_spe_c) = d/dy(q_mix_c / (1 + q_tot_mix)) and q_tot_mix should be a constant under phase transformation... ]
# note then d/dy (q_mix) = d/dy(q_spe / (1 - q_spe)) = d/dy(q_spe) / (1 - q_spe)^2  [ really it's more like d/dy(q_mix_c) = d/dy(q_spe_c / (1 - q_tot_spe)) and q_tot_spe should be a constant under phase transformation... ]

NOTE: These mixing ratio, specific humidity, and their derivatives are essentially identical at earthlike conditions...
"""





function get_params_and_go_to_mixing_ratio_exponential_part_only(param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::Real,
    ts::TD.ThermodynamicState;
    use_fix::Bool = false, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it),
    ) where {FT}

    # We have chosen not to bother with mixing ratios in this func... It shouldn't be neceessary.. it's just a scaling and qt is constant... when used as a fallback in morrison_milbrandt_2015_style() we just adjust accordingly....
    
    q_vap = TD.vapor_specific_humidity(q) # ensure it's not in mixing ratio units... [ note this  can cause a problem if q_vap doesn't equal the q_vap passed in due to floating point errors... ]
    q_liq = q.liq
    q_ice = q.ice
    q_sl = q_eq.liq
    q_si = q_eq.ice

    # CAN WE GET AWAY WITH NOT USING MIXING RATIO ANYWHERE? WE HAVE RETURN_MIXING_RATIO = FALSE EVERYWHERE! It's all just a scaling factor of 1/(1+q_tot) or 1/(1-q_tot)... don't think it changes the equations at all...
    # q_vap = TD.shum_to_mixing_ratio(TD.vapor_specific_humidity(q), q.tot) # convert to mixing ratio (use this version instead of q_vap bc q_vap passed from other functions into this fcn can be in mix form already than )
    # q_liq = TD.shum_to_mixing_ratio(q.liq, q.tot)
    # q_ice = TD.shum_to_mixing_ratio(q.ice, q.tot)
    # q_sl = TD.shum_to_mixing_ratio(q_eq.liq, q.tot) # q_eq is the equilibrium mixing ratio, q_sl is the saturation mixing ratio, q_eq we made to contain only saturation vapor pressure values...
    # q_si = TD.shum_to_mixing_ratio(q_eq.ice, q.tot)


    T_freeze = TCP.T_freeze(param_set)

    δ_0 = q_vap - q_sl # supersaturation over liquid
    # δ_0i = δ_0 + (q_sl - q_si) # supersaturation over ice
    δ_0i = q_vap - q_si # =  δ_0 + (q_sl - q_si) but fewer operations, less floating point precise ? # supersaturation over ice

    # @debug "δ_0 = $δ_0; δ_0i = $δ_0i"

    δ_0 = limit_δ(δ_0, q_vap, q_liq, q_ice) # because we're always heading towards WBF, if we're `close enough` we'll just call it WBF to avoid misses.
    δ_0i = limit_δ(δ_0i, q_vap, q_liq, q_ice) # because we're always heading towards WBF, if we're `close enough` we'll just call it WBF to avoid misses.

    # @debug "δ_0 = $δ_0; δ_0i = $δ_0i"

    return (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i)
end

# τ_func_exponential_part_only(τ_liq::FT, τ_ice::FT) where {FT} = 1 / (1 / τ_liq +  1 / τ_ice) # Eqn C2
function τ_func_exponential_part_only(τ_liq::FT, τ_ice::FT) where {FT}
    τ = min( 1 / (1 / τ_liq +  1 / τ_ice), τ_liq, τ_ice) # Eqn C2 [ Better for Floating point safety... ]
    return iszero(τ) ? min(τ_liq, τ_ice) : τ # if τ is 0, we'll just use the minimum of the two... Probably we faced an underflow...
end

const τ_func_EPA = τ_func_exponential_part_only

# THESE ARE WRITTEN ASSUMING WE AREN'T USING MIXING RATIO ANYWHERE, IS THAT BAD? OTHERWISE WE CAN JUST USE THE REGULAR Q_NEW FUNCTIONS FROM morrison_milbrandt_2015_style.jl
q_new_exponential_part_only(q::TD.PhasePartition, S_ql_Δt::FT, S_qi_Δt::FT) where {FT} = TD.PhasePartition(q.tot, q.liq + S_ql_Δt, q.ice + S_qi_Δt) # Eqn C3
q_new_exponential_part_only(q::TD.PhasePartition, S_ql::FT, S_qi::FT, Δt::FT) where {FT} = TD.PhasePartition(q.tot, q.liq + S_ql * Δt, q.ice + S_qi * Δt) # Eqn C3
q_new_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql_Δt::FT, S_qi_Δt::FT) where {FT} = TD.PhasePartition(q.tot, q_liq + S_ql_Δt, q_ice + S_qi_Δt) # Eqn C3
q_new_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql::FT, S_qi::FT, Δt::FT) where {FT} = TD.PhasePartition(q.tot, q_liq + S_ql * Δt, q_ice + S_qi * Δt) # Eqn C3
const q_new_EPA = q_new_exponential_part_only


function morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q::TD.PhasePartition)
    new_q_vap = TD.vapor_specific_humidity(new_q)
    return new_q, new_q_vap
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, S_ql_Δt::FT, S_qi_Δt::FT) where {FT}
    new_q = q_new_exponential_part_only(q, S_ql_Δt, S_qi_Δt)
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, S_ql::FT, S_qi::FT, Δt::FT) where {FT}
    new_q = q_new_exponential_part_only(q, S_ql, S_qi, Δt)
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql_Δt::FT, S_qi_Δt::FT) where {FT}
    # This version should be more floating point stable -- whatever you calculate in mixing ratio you should get exactly out.
    new_q = q_new_exponential_part_only(q, q_liq, q_ice, S_ql_Δt, S_qi_Δt) # use safe version
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end

function morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q::TD.PhasePartition, q_liq::FT, q_ice::FT, S_ql::FT, S_qi::FT, Δt::FT) where {FT}
    # This version should be more floating point stable -- whatever you calculate in mixing ratio you should get exactly out.
    new_q = q_new_exponential_part_only(q, q_liq, q_ice, S_ql, S_qi, Δt) # use safe version
    return morrison_milbrandt_2015_get_new_status_helper_helper_exponential_part_only(new_q)
end
const morrison_milbrandt_2015_get_new_status_helper_EPA = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only


A_c_func_no_WBF_exponential_part_only(FT) = FT(0) # don't use this just use FT(0) where this comes up...
A_c_func_exponential_part_only(τ_ice::FT, q_sl::FT, q_si::FT) where {FT} =  - (q_sl - q_si) / τ_ice # Eq C4
const A_c_func_no_WBF_EPA = A_c_func_no_WBF_exponential_part_only
const A_c_func_EPA = A_c_func_exponential_part_only



function δ_func_exponential_part_only(A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT}
    term_1 = (δ_0 - A_c * τ)
    term_2 = exp(-Δt / τ)
    prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
    return  A_c * τ  + prod # EQN C5
end
const δ_func_EPA = δ_func_exponential_part_only

δi_func_exponential_part_only(δ_0::FT, A_c::FT, τ::FT, Δt::FT, q_sl::FT, q_si::FT) where {FT} = δ_func_exponential_part_only(A_c, τ, δ_0, Δt) + (q_sl - q_si) # EQN C5
const δi_func_EPA = δi_func_exponential_part_only
# indiv for if we need this for just ice with δi and no liquid change
const δi_func_indiv_exponential_part_only = δ_func_exponential_part_only
const δi_func_indiv_EPA = δi_func_indiv_exponential_part_only


# Could be useful for floating point accuracy when dδ < nextfloat(δ) maybe -- seems not, prod_1 and prod_2 differ by neg sign and small amount can lead to prod = nextfloat(prod_1)-prod_1 rather than true value
function dδ_func_exponential_part_only(A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT} 
    term = 1 - exp(-Δt / τ)
    prod_1 = iszero(term) ? FT(0) : ((A_c * τ) * term) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
    prod_2 = iszero(term) ? FT(0) : (-δ_0 * term) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
    prod = prod_1 + prod_2
    return iszero(term) ? FT(0) : prod # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
end
const dδ_func_EPA = dδ_func_exponential_part_only

""" One species τ_species = τ, no WBF """
function S_func_indiv_no_WBF_exponential_part_only( A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT}
    if isfinite(τ) # avoid NaN problems at τ = Inf
        term_1 = (δ_0 - A_c * τ) / (Δt) # can't be nan bc it's finite * finite, neither of which should be 0
        term_2 = (1 - exp(-Δt / τ))
        prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem

        # if !isfinite(term_1) || !isfinite(term_2) || !isfinite(prod) # if any of these are NaN, we'll just return 0
        # end

        return A_c + prod # QCCON EQN C6 # Use this form bc we can keep the first term if prod was gonna be NaN
        # return A_c + (δ_0 - A_c * τ)  / (Δt) * (1 - exp(-Δt / τ)) # QCCON EQN C6 [ did I deprecate this to create the lower structure? need to document better ] This way is ad bc if prod is NaN you also lose the first term.
    else
        return FT(0) # I know it looks like it should be A_c bc τ and τ_c canceled, but infinite τ means no change... Really what this is saying is τ_c == τ == τ_<other> = ∞ , the limit at τ -> ∞ of τ/τ_c goes to 0 not 1. It's a simplification of S_func_no_WBF_exponential_part_only().
    end
end
const S_func_indiv_no_WBF_EPA = S_func_indiv_no_WBF_exponential_part_only
#
const S_ql_func_indiv_exponential_part_only = S_func_indiv_no_WBF_exponential_part_only
const S_qi_func_indiv_exponential_part_only = S_func_indiv_no_WBF_exponential_part_only
const S_ql_func_indiv_EPA = S_ql_func_indiv_exponential_part_only
const S_qi_func_indiv_EPA = S_qi_func_indiv_exponential_part_only



""" One of two species i.e. τ_species ≠ τ, no WBF """
function S_func_no_WBF_exponential_part_only(A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT) where {FT}
    # if τ is Inf, then τ_c and the other τ must have both been Inf
        # if (!isfinite(τ) && isfinite(τ_c)); error("τ is Inf so τ_c and the other τ must have both been Inf but got τ = Inf and τ_c = $τ_c"); end

    # if τ_c is Inf but τ is not, then τ other must not have been Inf. Either way, S is 0.
    if !isfinite(τ_c) # avoid NaN problems at τ = Inf
        return FT(0)
    else
        # @debug "S_func_no_WBF_exponential_part_only: A_c = $A_c; τ = $τ; τ_c = $τ_c; δ_0 = $δ_0; Δt = $Δt"
        # term_1 = (δ_0 - A_c * τ) *  τ / (Δt * τ_c) # can't be nan bc it's finite * finite, neither of which should be 0
        term_1 = (δ_0 - A_c * τ) *  (τ / τ_c)/(Δt) # can't be nan bc it's finite * finite, neither of which should be 0
        if isinf(term_1) #
            # We'll just return the largest number we can... it's still smaller than the correct answer so hopefully we don't end up in any infinite loops... [the tests should hopefully catch problems?]
            # return floatmax(FT) * sign(term_1) # this is a problem, we need to fix it. 
            alt = SΔt_func_no_WBF_exponential_part_only(A_c, τ, τ_c, δ_0, Δt)  / Δt # this is a problem, we need to fix it.
            return isinf(alt) ? floatmax(FT) * sign(alt) : alt
        end

        term_2 = (1 - exp(-Δt / τ)) # this term can be imprecise for small exp(-Δt / τ) leading to floating point problems...
        prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
        S =  A_c * τ / τ_c + prod # QCCON EQN C6 # Use this form bc we can keep the first term if prod was gonna be NaN
        # return (sign(S) == sign(δ_0)) ? S : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
        return S
    end
end

const S_func_no_WBF_EPA = S_func_no_WBF_exponential_part_only
#
const S_ql_func_no_WBF_exponential_part_only = S_func_no_WBF_exponential_part_only
const S_ql_func_EPA = S_ql_func_no_WBF_exponential_part_only
const S_qi_func_no_WBF_exponential_part_only = S_func_no_WBF_exponential_part_only
const S_qi_func_no_WBF_EPA = S_qi_func_no_WBF_exponential_part_only

function SΔt_func_no_WBF_exponential_part_only(A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT) where {FT}
    # if τ_c is Inf but τ is not, then τ other must not have been Inf. Either way, S is 0.
    if !isfinite(τ_c) # avoid NaN problems at τ = Inf
        return FT(0)
    else
        # term_1 = (δ_0 - A_c * τ) *  τ / (Δt * τ_c) # can't be nan bc it's finite * finite, neither of which should be 0
        term_1 = (δ_0 - A_c * τ) *  (τ / τ_c) # can't be nan bc it's finite * finite, neither of which should be 0
        if isinf(term_1) #
            return floatmax(FT) * sign(term_1) # this is a problem, we need to fix it. We'll just return the largest number we can... it's still smaller than the correct answer so hopefully we don't end up in any infinite loops... [the tests should hopefully catch problems?]
        end

        term_2 = (1 - exp(-Δt / τ))
        prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
        SΔt =  (A_c * τ * Δt) / τ_c + prod # QCCON EQN C6 # Use this form bc we can keep the first term if prod was gonna be NaN

        # @error("I think this is bad bc WBF is in A so you can get a different sign than δ_0 from WBF")
        # return (sign(SΔt) == sign(δ_0)) ? SΔt : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
        return SΔt
    end
end

const SΔt_func_no_WBF_EPA = SΔt_func_no_WBF_exponential_part_only
const SΔt_ql_func_no_WBF_exponential_part_only = SΔt_func_no_WBF_exponential_part_only
const SΔt_ql_func_EPA = SΔt_ql_func_no_WBF_exponential_part_only
const SΔt_qi_func_no_WBF_exponential_part_only = SΔt_func_no_WBF_exponential_part_only
const SΔt_qi_func_no_WBF_EPA = SΔt_qi_func_no_WBF_exponential_part_only


""" Both species - default is WBF on ice, but can call the other way with appropriate replacements"""
function S_func_exponential_part_only(A_c::FT, τ::FT, τ_ice::FT, δ_0::FT, Δt::FT, q_sl::FT, q_si::FT) where {FT}
    S = S_func_no_WBF_exponential_part_only( A_c, τ, τ_ice, δ_0, Δt) + (q_sl - q_si) / (τ_ice) #  # QICON Eqn C7 # if τ_ice is Inf this works fine,  no issues...
    # @error("I think this line is bad bc if you're passing in δ_0 and operating on WBF sign can flip")
    # return (sign(S) == sign(δ_0)) ? S : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
    return S
end

const S_func_EPA = S_func_exponential_part_only
const S_qi_func_exponential_part_only = S_func_exponential_part_only
const S_qi_func_EPA = S_qi_func_exponential_part_only

""" Both species - default is WBF on ice, but can call the other way with appropriate replacements"""
function SΔt_func_exponential_part_only(A_c::FT, τ::FT, τ_ice::FT, δ_0::FT, Δt::FT, q_sl::FT, q_si::FT) where {FT}
    SΔt = SΔt_func_no_WBF_exponential_part_only( A_c, τ, τ_ice, δ_0, Δt) + (q_sl - q_si) * Δt / (τ_ice) #  # QICON Eqn C7 # if τ_ice is Inf this works fine,  no issues...
    # @error("I think this line is bad bc if you're passing in δ_0 and operating on WBF sign can flip")
    # return (sign(SΔt) == sign(δ_0)) ? SΔt : FT(0) # probably is a small floating point problem -- just return 0. Honestly we should probably have done this for all large τ, τ_c from the beginning... Now that we have exact/prognostic δ tracking instead of backing out, maybe we should do it for small δ too... edge cases are killing us
    return SΔt
end
const SΔt_func_EPA = SΔt_func_exponential_part_only
const SΔt_qi_func_exponential_part_only = SΔt_func_exponential_part_only
const SΔt_qi_func_EPA = SΔt_qi_func_exponential_part_only



get_t_out_of_q_no_WBF_EPA(δ_0::FT, A_c::FT, τ::FT, τ_c::FT, q_c::FT) where{FT} = get_t_out_of_q_no_WBF(δ_0, A_c, τ, τ_c, q_c, FT(1)) # easier to use this than rewrite, just fill in Γ = 1
const get_t_out_of_q_liq_EPA = get_t_out_of_q_no_WBF_EPA
const get_t_out_of_q_ice_no_WBF_EPA = get_t_out_of_q_no_WBF_EPA

get_t_out_of_q_WBF_EPA(δ_0::FT, A_c::FT, τ::FT, τ_c::FT, q_ice::FT, q_sl::FT, q_si::FT) where {FT} = get_t_out_of_q_WBF(δ_0, A_c, τ, τ_c, q_ice, FT(1), q_sl, q_si) # easier to use this than rewrite, just fill in Γ = 1
const get_t_out_of_q_ice_EPA = get_t_out_of_q_WBF_EPA

# Deprecated for now, I don't think you can actually reach the exact equilibrium point in q_new_exponential_part_only bc we're ignoring external forcings...
# see https://www.wolframalpha.com/input?i2d=true&i=solve+-Divide%5Bx%2Cc%5D%3DDivide%5B%5C%2840%29x%2Ba%5C%2841%29%2Ck%5D+for+x
# function calc_WBF_δ_equilibrium(q_sl::FT, q_si::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT) where {FT}
#     # see writeup
#     return -(q_sl - q_si) * (Γ_l * τ_liq) / ( (Γ_i * τ_ice + Γ_l * τ_liq))
# end

# # see https://www.wolframalpha.com/input?i2d=true&i=solve+-Divide%5Bx%2Cc%5D%3DDivide%5B%5C%2840%29x%2Ba%5C%2841%29%2Ck%5D+for+x
# function calc_WBF_δ_equilibrium_exponential_part_only(q_sl::FT, q_si::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT) where {FT}
#     # see writeup
#     return -(q_sl - q_si) * (τ_liq) /  (τ_ice + τ_liq)
# end
# calc_WBF_δ_equilibrium_EPA = calc_WBF_δ_equilibrium_exponential_part_only


# S_func_exponential_part_only() can have nextfloat problem for large numbers [e.g. nextfloat(S_func_no_WBF_exponential_part_only( A_c, τ, τ_ice, δ_0, Δt)) can be arbitrarily large as S_func_no_WBF_exponential_part_only grows, sometimes leading to S_ql, S_qi that are too big...]
# we could switch to BigFloats which would solve thigns but they're oh so slow and allocate like crazy. Just use this stopgap for now...
""" 
    Unlike S_above_floor(), this does rely on the fact that δ_0, δ_0i are not affected by external forcings and so are only consumed by cond/evap sub/dep.
    In the face of w, mixing, and other forcings, maybe you would jus want to fall back to big float. [or maybe something like DoubleFloats.jl for speed]

    `q_other` is the other species... q_liq if q_ice is being updated and vice versa.

    We also need to decrease gains but not increase losses, so we use switch state
"""
function S_below_ceiling_EPA(S_qc::FT, q_other::FT, δ_0::FT, Δt::FT) where {FT} # Don't exceed consuming all supersat + q_other in a timestep
    (S_qc > FT(0)) ? min(S_qc, (δ_0 + q_other) / Δt, floatmax(FT)) : max(S_qc, -floatmax(FT)) # Assuming gains are limited by first removing any subsaturation then either consuming other species or consuming supersaturation. if not gaining, then don't use limit.
   
    # temp form for debugging problems. The line above should work fine... once all the kinks (elsewhere not in this func) are worked out
    if S_qc > FT(0)
        if (δ_0 + q_other) / Δt < FT(0)
            error("δ_0 + q_other < 0 while S_qc > 0")
        end
        return min(S_qc, (δ_0 + q_other) / Δt, floatmax(FT)) # Assuming gains are limited by first removing any subsaturation then either consuming other species or consuming supersaturation.
    else
        return max(S_qc, -floatmax(FT)) # if not gaining, then don't use limit.
    end
end
@inline S_below_ceiling_EPA(S_ql::FT, S_qi::FT, q_liq::FT, q_ice::FT, δ_0::FT, δ_0i::FT, Δt::FT) where {FT} = S_below_ceiling_EPA(S_ql, q_ice, δ_0, Δt), S_below_ceiling_EPA(S_qi, q_liq, δ_0i, Δt) 


"""
From a list of numbers, return the number with the smallest magnitude
"""
function smallest_magnitude(numbers::Vararg{FT}) where {FT}
    abs_values = abs.(collect(numbers)) # Collect to Vector for argmin
    idx = argmin(abs_values)
    return numbers[idx]
end



function clamp_δ(δ_0::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, q_sl::FT, q_si::FT) where {FT}
    if below_freezing
        if regime_type <: Subsaturated
            return min(δ_0, -(q_sl - q_si))
        elseif regime_type <: WBF
            return clamp(δ_0, -(q_sl - q_si), FT(0))
        elseif regime_type <: Supersaturated
            return max(δ_0, FT(0))
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if regime_type <: Subsaturated
            return min(δ_0, FT(0))
        elseif regime_type <: WBF
            return clamp(δ_0, FT(0), q_si - q_sl)
        elseif regime_type <: Supersaturated
            return max(δ_0, q_si - q_sl)
        else
            error("Unknown regime type: $regime_type")
        end
    end
end         
clamp_δ(δ_0::FT, regime_type::Type{<:AbstractSaturationRegime}, q_sl::FT, q_si::FT) where {FT} = clamp_δ(δ_0, regime_type, regime_type.parameters[3], q_sl, q_si) # don't use this one
clamp_δ(δ_0::FT, regime::AbstractSaturationRegime, q_sl::FT, q_si::FT) where {FT} = clamp_δ(δ_0, typeof(regime), regime.below_freezing, q_sl, q_si)


function clamp_δi(δ_0i::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, q_sl::FT, q_si::FT) where {FT}
    if below_freezing
        if regime_type <: Subsaturated
            return min(δ_0i, FT(0))
        elseif regime_type <: WBF
            return clamp(δ_0i, FT(0), q_sl - q_si)
        elseif regime_type <: Supersaturated
            return max(δ_0i, q_sl - q_si)
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if regime_type <: Subsaturated
            return min(δ_0i, -(q_si - q_sl))
        elseif regime_type <: WBF
            return clamp(δ_0i, -(q_si - q_sl), FT(0))
        elseif regime_type <: Supersaturated
            return max(δ_0i, FT(0))
        else
            error("Unknown regime type: $regime_type")
        end
    end
end
clamp_δi(δ_0i::FT, regime_type::Type{<:AbstractSaturationRegime}, q_sl::FT, q_si::FT) where {FT} = clamp_δi(δ_0i, regime_type, regime_type.parameters[3], q_sl, q_si) # don't use this one
clamp_δi(δ_0i::FT, regime::AbstractSaturationRegime, q_sl::FT, q_si::FT) where {FT} = clamp_δi(δ_0i, typeof(regime), regime.below_freezing, q_sl, q_si)
    

function clamp_S_ql(S_ql::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, δ_0::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT}
    if below_freezing
        if (regime_type <: Subsaturated) || (regime_type <: WBF)
            return clamp(S_ql, max(-q_liq / Δt, -floatmax(FT)), FT(0))
        elseif regime_type <: Supersaturated
            return clamp(S_ql, FT(0), min((δ_0 + q_ice) / Δt, floatmax(FT)))
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if (regime_type <: Subsaturated)
            return clamp(S_ql, FT(0), q_liq / Δt)
        elseif (regime_type <: WBF) || (regime_type <: Supersaturated)
            return clamp(S_ql, FT(0), min((δ_0 + q_ice) / Δt, floatmax(FT)))
        else
            error("Unknown regime type: $regime_type")
        end
    end
end
clamp_S_ql(S_ql::FT, regime_type::Type{<:AbstractSaturationRegime}, δ_0::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S_ql(S_ql, regime_type, regime_type.parameters[3], δ_0, q_liq, q_ice, Δt) # don't use this one
clamp_S_ql(S_ql::FT, regime::AbstractSaturationRegime, δ_0::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S_ql(S_ql, typeof(regime), regime.below_freezing, δ_0, q_liq, q_ice, Δt)

function clamp_S_qi(S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT}
    if below_freezing
        if (regime_type <: Subsaturated)
            return clamp(S_qi, max(-q_ice / Δt, -floatmax(FT)), FT(0))
        elseif (regime_type <: WBF) || (regime_type <: Supersaturated)
            return clamp(S_qi, FT(0), min((δ_0i + q_liq) / Δt, floatmax(FT)))
        else
            error("Unknown regime type: $regime_type")
        end
    else
        if (regime_type <: Subsaturated) || (regime_type <: WBF)
            return clamp(S_qi, max(-q_ice / Δt, -floatmax(FT)), FT(0))
        elseif (regime_type <: Supersaturated)
            return clamp(S_qi, FT(0), min((δ_0i + q_liq) / Δt, floatmax(FT)))
        else
            error("Unknown regime type: $regime_type")
        end
    end
end
clamp_S_qi(S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S_qi(S_qi, regime_type, regime_type.parameters[3], δ_0i, q_liq, q_ice, Δt) # don't use this one
clamp_S_qi(S_qi::FT, regime::AbstractSaturationRegime, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S_qi(S_qi, typeof(regime), regime.below_freezing, δ_0i, q_liq, q_ice, Δt)

clamp_S(S_ql::FT, S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, below_freezing::Bool, δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S_ql(S_ql, regime_type, below_freezing, δ_0, q_liq, q_ice, Δt), clamp_S_qi(S_qi, regime_type, below_freezing, δ_0i, q_liq, q_ice, Δt)
clamp_S(S_ql::FT, S_qi::FT, regime_type::Type{<:AbstractSaturationRegime}, δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S(S_ql, S_qi, regime_type, regime_type.parameters[3], δ_0, δ_0i, q_liq, q_ice, Δt) # don't use this one
clamp_S(S_ql::FT, S_qi::FT, regime::AbstractSaturationRegime, δ_0::FT, δ_0i::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = clamp_S(S_ql, S_qi, typeof(regime), regime.below_freezing, δ_0, δ_0i, q_liq, q_ice, Δt)
            



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

function morrison_milbrandt_2015_style_exponential_part_only(
    param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::Real,
    ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    emit_warnings::Bool = true,
) where {FT}

    # TODO: Track supersaturation directly... should fix floating point problems...

    if area > FT(0)

        if emit_warnings && (Δt < eps(FT))
            @debug "Timestep $(Δt) is very small (smaller than eps(FT) = $(eps(FT))), may cause numerical issues..."
        end

        δ_0 = q_vap - q_eq.liq # supersaturation over liquid
        δ_0i = δ_0 + (q_eq.liq - q_eq.ice) # supersaturation over ice

        # has_liq::Bool = q.liq > FT(0)
        # has_ice::Bool = q.ice > FT(0)
        below_freezing::Bool = T < TCP.T_freeze(param_set)
        regime = get_regime(q_vap, q, q_eq, below_freezing)
        return morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts; use_fix = use_fix, return_mixing_ratio = false, δ_0=δ_0, δ_0i=δ_0i)
        
    else
        return FT(0), FT(0)
    end
end

# ====================================================================================================================================================================================================== #
# ====================================================================================================================================================================================================== #
# ====================================================================================================================================================================================================== #

"""
    Supersaturated (below freezing)
    - Both can grow!
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Supersaturated{true, true, true}, Supersaturated{true, false, true}, Supersaturated{false, true, true}, Supersaturated{false, false, true}}, # these should all work the same, right? you'll end up with some liq/ice at the end no matter what
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end

    @debug "Calling Supersaturated{$(q.liq > FT(0)), $(q.ice > FT(0)), true}..."

    # --- Thermo  constants ------------------------------------------------------------------------------------ #
    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)


    # Both are growing..., no explicit need for WBF right
    A_c = A_c_func_EPA(τ_ice, q_sl, q_si)
    τ = τ_func_EPA(τ_liq, τ_ice)

    t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ) # Eq C5

    min_t, i_min_t = find_min_t([t_hit_liq_sat])  # find_min_t helps resolve if min_t is 0 for example, don't skip this call
    if min_t < Δt
        @debug "will hit liq sat before timestep is over... will transition to wbf at t = $(t_hit_liq_sat)..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t) # This includes the wbf part though...
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, q_sl, q_si) #

        if isinf((S_ql+S_qi)*min_t) # scale down just to hit WBF. This can happen when the timescale is too short to calculate. This will crash the model when you go to calculate new_q for example
            liq_frac = τ/τ_liq
            S_ql = δ_0 * liq_frac
            S_qi = δ_0 * (1-liq_frac)
        end

        # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, min_t) # don't consume more than you have...
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t)

        S_ql_Δt = S_ql * min_t
        S_qi_Δt = S_qi * min_t
        # S_tot = S_ql_Δt + S_qi_Δt
        
        Δt_left = Δt - min_t


        # dδ_from_S = -(S_ql * min_t + S_qi * min_t)
        # dδ_from_theory = dδ_func_EPA(A_c, τ, δ_0, min_t)
        dδ = smallest_magnitude(-(S_ql * min_t + S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.


        new_δ_0 = FT(0) # hit liq sat from above, so δ_0 = 0
        new_δ_0i = δ_0i + dδ
        new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)


        new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q, q_liq, q_ice, S_ql_Δt, S_qi_Δt)
        # not sure if we can have underflow problems here... don't think so
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(WBF{true, true, true}(true, true, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
        S_ql *= (min_t / Δt) # rescale to the timestep
        S_qi *= (min_t / Δt) # rescale to the timestep
        S_ql_addit *= Δt_left / Δt # rescale to the remaining time
        S_qi_addit *= Δt_left / Δt # rescale to the remaining time
        return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
    else
        @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, q_sl, q_si)
        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, Δt) # don't consume more than you have...
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Supersaturated (above freezing)
    - No ice growth, liquid only.
    - Existing ice remains, do not enforce melting here.
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Supersaturated{true, true, false}, Supersaturated{true, false, false}, Supersaturated{false, true, false}, Supersaturated{false, false, false}}, # these should all work the same, right? you'll end up with some liq/ice at the end no matter what
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end
    
    @debug "Calling Supersaturated{$(q.liq > FT(0)), $(q.ice > FT(0)), false}..."


    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    # A_c = A_c_func_EPA(τ_ice, q_sl, q_si)
    A_c = A_c_func_no_WBF_EPA(FT)
    # τ = τ_liq # supersaturated, so no ice decay, but above freezing, so no ice growth


    # t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ_liq) # Eq C5 (can only happen bc of A_c) we're above freezing so we hit ice sat first?
    # min_t, i_min_t = find_min_t([t_hit_liq_sat])  # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    # t_hit_ice_sat = t_δ_hit_value(q_sl - q_si, δ_0, A_c, τ_liq) # Eq C5 (can only happen bc of A_c) [ we're above freezing so we hit ice sat first?]
    t_hit_ice_sat = t_δ_hit_value(q_si - q_si, δ_0, A_c, τ_liq) # Eq C5 (can only happen bc of A_c) [ we're above freezing so we hit ice sat first?]
    min_t, i_min_t = find_min_t([t_hit_ice_sat])  # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    if min_t < Δt # bc of A_c

        @debug "will hit ice sat before timestep is over... will transition to wbf at t = $(min_t)..."
        S_ql = S_ql_func_indiv_EPA( A_c, τ_liq, δ_0, min_t) # This includes the wbf part though...
        S_qi = FT(0)

        # S_ql = S_above_floor(S_ql, q_liq, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, min_t) # don't consume more than you have...
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t)

        new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t)

        dδ = smallest_magnitude(-(S_ql * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
        new_δ_0 = δ_0 + dδ
        new_δ_0i = FT(0) # hit ice sat
        new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
        
        Δt_left = Δt - min_t
        # maybe we should check that new_q actually has values, bc if τ or δ or Δt is really small it might not...
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(WBF{true, false, false}(true, false, false), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
        S_ql *= (min_t / Δt) # rescale to the timestep
        S_qi *= (min_t / Δt) # rescale to the timestep
        S_ql_addit *= Δt_left / Δt # rescale to the remaining time
        S_qi_addit *= Δt_left / Δt # rescale to the remaining time
        return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
    else
        # exponential decay, so will never hit liq sat
        S_ql = S_ql_func_indiv_EPA( A_c, τ_liq, δ_0, Δt)
        S_qi = FT(0)
        # S_ql = S_above_floor(S_ql, q_liq, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, Δt) # don't consume more than you have...
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF has liquid (below freezing)
    - liq can shrink, ice can grow

    - should asymptote to equilibrium, you shouldn't be able to hit sat
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{true, true, true}, WBF{true, false, true}}, # has liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end


    @debug "Calling WBF{true, $(q.ice > FT(0)), true}..."


    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)    
    


    # One is growing, one is shrinking...
    A_c = A_c_func_EPA(τ_ice, q_sl, q_si)
    τ = τ_func_EPA(τ_liq, τ_ice)


    t_out_of_liq = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq) # Eq C6

    # δ_equil = calc_WBF_δ_equilibrium_EPA(q_sl, q_si, τ_liq, τ_ice)
    # t_hit_δ_equil = t_δ_hit_value(δ_equil, δ_0, A_c, τ) # I think this should always be Inf


    # can run out of liq or reach equilibrium (but you'll never reach equilibrium right?)
    # can't reach ice sat bc liq evap would pull you up, can't hit liq sat bc ice would pull you down. equil is guaranteed to be between 0 and (q_sl - q_si) [see writeup]

    min_t, i_min_t = find_min_t([t_out_of_liq]) # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    # @info "min_t = $min_t; t_out_of_liq = $t_out_of_liq; Δt = $Δt"

    if min_t < Δt
        if i_min_t == 1
            @debug "liq will run out first before timestep is over... will transition at t = $(min_t) to just ice growth"
            S_ql = -q_liq / min_t
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, q_sl, q_si)  # however much happens in that time, rescaled to the timestep [ if this underflows, then what? ]
            Δt_left = Δt - min_t

            # S_qi = S_above_floor(S_qi, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_qi = S_below_ceiling_EPA(S_qi, q_liq, δ_0i, min_t) # don't consume more than you have...
            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t)

            dδ = smallest_magnitude(-(-q_liq + S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            # if S_qi is reduced to floatmax due to overflow, then dδ will be larger for the actual δ_0 than for the theoretical one... so that could cause a problem... still that should be a very rare case...?


            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)
            
            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q, q_liq, q_ice, -q_liq, S_qi*min_t) # use multiplied form for floating point accuracy
            # new_regime = get_regime_err_WBF(new_q_vap, new_q.liq, new_q.ice, q_sl, q_si, T, T_freeze) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, T<T_freeze)
            # new_regime = get_regime(δ_0, δ_0i, new_q.liq, new_q.ice, true)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            # S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(WBF{false, true, true}(false, true, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1)  # just ice growth until sat


            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep 
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        # else # i_min_t == 2 # if you hit equil
        # Don't think you should be able to reach here...
        end
        
    else
        @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, q_sl, q_si)
        # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, Δt)
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF has liquid amd ice (above freezing)
    - liq can grow, no ice growth, but ice can shrink
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::WBF{true, true, false}, # has liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end

    @debug "Calling WBF{true, true, false}..."

    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    # No ice growth, only liq shrink
    A_c = A_c_func_EPA(τ_ice, q_sl, q_si)
    τ = τ_func_EPA(τ_liq, τ_ice)

    # t_out_of_liq = t_dδ_hit_value(q_liq, δ_0, A_c, τ) # We only have liq, so this works, don't need the product-log/lambertw stuff (wrong, need to use the product-log stuff bc supersat changing is not 1:1 with q_liq changing)
    t_out_of_liq = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq)
    t_out_of_ice = get_t_out_of_q_ice_EPA(δ_0, A_c, τ, τ_ice, q_ice, q_sl, q_si) # Eq C6
    t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ_liq)
    # hitting ice sat does nothing bc we're above freezing so we ignore

    min_t, i_min_t = find_min_t([t_out_of_liq, t_out_of_ice, t_hit_liq_sat]) # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    # @info "min_t = $min_t; i_min_t = $i_min_t; Δt = $Δt; t_out_of_liq = $t_out_of_liq; t_out_of_ice = $t_out_of_ice; t_hit_liq_sat = $t_hit_liq_sat"


    if min_t < Δt
        if i_min_t == 1 # out of liq, we're done
            S_ql = -q_liq / Δt # sacle to entire timestep (there's no addit)
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, q_sl, q_si) # however much happens in that time, rescaled to the timestep
            # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, min_t)
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t)
            S_qi *= min_t / Δt # rescale to the timestep
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        elseif i_min_t == 2 # out of ice, transition
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t)
            S_qi = -q_ice / min_t
            Δt_left = Δt - min_t
            # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, min_t)
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t)

            dδ = smallest_magnitude(-(S_ql * min_t - q_ice), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)

            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q, q_liq, q_ice, S_ql*min_t, -q_ice) # use multiplied form for floating point accuracy
            # new_regime = get_regime_err_WBF(new_q_vap, new_q.liq, new_q.ice, q_sl, q_si, T, T_freeze) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, false)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        elseif i_min_t == 3 # hit_liq_sat, we're done (stable)
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t)
            S_qi = FT(0)
            # S_ql = S_above_floor(S_ql, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, min_t)
            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t)

            S_ql *= min_t / Δt # rescale to the timestep
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        end

    else # 
        @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt)
        S_qi = FT(0)
        # S_ql = S_above_floor(S_ql, q_liq, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, Δt)
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

"""
    WBF has liq (above freezing)
    - liq can grow
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::WBF{true, false, false}, # has liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end

    @debug "Calling WBF{true, false, false}..."
    
    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    #  liq can grow
    τ = τ_liq # only liq
    A_c = A_c_func_no_WBF_EPA(FT) # just FT(0)

    t_out_of_liq = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq) # Eq C6

    t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ_liq) # Eq C5

    min_t, i_min_t = find_min_t([t_out_of_liq, t_hit_liq_sat]) # find_min_t helps resolve if min_t is 0 for example, don't skip this call

    if min_t < Δt
        if i_min_t == 1 # out of liq, we're done
            S_ql = q_liq / Δt # sacle to entire timestep (there's no addit)
            S_qi = FT(0)
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        else # i_min_t == 2 # hit liq sat, we're done (stable)
            S_ql = S_ql_func_indiv_EPA( A_c, τ_liq, δ_0, min_t) # This includes the wbf part though...
            S_qi = FT(0)
            # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, min_t)
            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t)

            S_ql *= min_t / Δt # rescale to the timestep
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        end

    else
        @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_indiv_EPA( A_c, τ_liq, δ_0, Δt)
        S_qi = FT(0)
        # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, Δt)
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt)

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF no liquid (below freezing)
    - no liq, just ice growth

    -- you should never hit saturation right?
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{false, true, true}, WBF{false, false, true}}, # no liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end

    @debug "Calling WBF{false, $(q.ice > FT(0)), true}..."

    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)


    # no liq so this is just ice growth.
    τ = τ_ice # we have no liq
    A_c = A_c_func_no_WBF_EPA(FT) # just FT(0)
     # We only have ice growth so it's like QCCON but with just ice 

    # can reach ice sat
    t_hit_ice_sat = t_δ_hit_value(FT(0), δ_0i, A_c, τ) # switch to using δ_0i and let that hit 0
    t_hit_liq_sat = t_δ_hit_value(q_sl-q_si, δ_0i, A_c, τ) # (only possible if A_c does it...)

    min_t, i_min_t = find_min_t([t_hit_ice_sat, t_hit_liq_sat])
    # @info "min_t = $min_t, i_min_t = $i_min_t | t_hit_ice_sat = $t_hit_ice_sat, t_hit_liq_sat = $t_hit_liq_sat"
    min_t = FT(Inf)

    if min_t < Δt
        # if i_min_t == 1
        #     @error("should be unreachable bc we're doing this w/ no external forcing")
        #     @debug "will hit ice sat before timestep is over... will transition to wbf at t = $(t_hit_ice_sat)"
        #     S_ql = FT(0)
        #     S_qi = FT(δ_0i / min_t)
        #     Δt_left = Δt - min_t

        #     dδ = smallest_magnitude(-(S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
        #     new_δ_0 = δ_0 + dδ
        #     new_δ_0i = FT(0) # we're at ice sat so δ_0i = 0
        #     new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)

            
        #     new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t)
        #     S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(LowerSatLine{false, true, true}(false, true, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i) # should have  no liq still, no ice above freezing
        #     S_qi *= t_hit_ice_sat / Δt # rescale to the timestep
        #     S_ql_addit *= Δt_left / Δt # rescale to the remaining time
        #     S_qi_addit *= Δt_left / Δt # rescale to the remaining time

        #     return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        # else # i_min_t == 2 # A_c has raised us up to liq_sat, transition to supersaturated? or stay on WBF line?
        #     @error("should be unreachable bc we're doing this w/ no external forcing")
        #     @debug "Reached liq sat, will terminate unless T > T_freeze"
        #     S_ql = FT(0)
        #     # S_qi = FT(δ_0i-δ_0 / min_t) # is only
        #     S_qi = FT(-((q_sl - q_si) - δ_0i) / min_t) #   Δδ_0i = (q_sl - q_si) - δ_0i --> Δq_ice = -Δδ_0i
        #     Δt_left = Δt - min_t

        #     dδ = smallest_magnitude(-(S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
        #     new_δ_0 = FT(0) # we're at liq sat so δ_0 = 0
        #     new_δ_0i = δ_0i + dδ
        #     new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)



        #     new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t)
        #     S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(Supersaturated{false, true, true}(false, true, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i) # should have  no liq still, no ice above freezing
        #     S_qi_addit = FT(0)
        #     S_ql_addit = FT(0)
        #     S_ql *= t_hit_liq_sat / Δt # rescale to the timestep
        #     S_qi *= t_hit_liq_sat / Δt # rescale to the timestep
        #     S_ql_addit *= Δt_left / Δt # rescale to the remaining time
        #     S_qi_addit *= Δt_left / Δt # rescale to the remaining time
        #     return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        # end

    else

        @debug "nothing of note through end of timestep..."
        # just regular return
        S_ql = FT(0)
        S_qi = S_func_indiv_no_WBF_EPA( A_c, τ, δ_0i, Δt) # use the no WBF version w/ just exponential decay of δ_0i
        # S_qi = S_below_ceiling_EPA(S_qi, q_liq, δ_0i, Δt)
        S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, Δt)  

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
    
"""
    WBF no liquid (above freezing)
    - no liq, no ice growth
"""

function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{WBF{false, true, false}, WBF{false, false, false}}, # no liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT} 
    return FT(0), FT(0)
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated with liq and ice
    - both shrink
    - if below freezing, would hit ice_sat first and stop
    - if above freezing, would hit liq_sat first and stop
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{true, true, true}, Subsaturated{true, true, false}}, # can run out of either first
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end


    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    BF::Bool = T < T_freeze
    @debug "Calling Subsaturated{true, true, $BF}..."

    # (!iszero(q_liq)) || error("liq should not be zero")
    # (!iszero(q_ice)) || error("ice should not be zero")

    A_c = A_c_func_EPA(τ_ice, q_sl, q_si)
    τ = τ_func_EPA(τ_liq, τ_ice)

    t_hit_sat = BF ? t_δ_hit_value(q_si-q_sl, δ_0, A_c, τ) : t_δ_hit_value(FT(0), δ_0, A_c, τ) # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower

    t_out_of_liq = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq)
    t_out_of_ice = get_t_out_of_q_ice_EPA(δ_0, A_c, τ, τ_ice, q_ice, q_sl, q_si)

    t_hit_sat = BF ? t_δ_hit_value(q_si-q_sl, δ_0, A_c, τ) : t_δ_hit_value(FT(0), δ_0, A_c, τ) # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower
    
    min_t, i_min_t = find_min_t([t_out_of_liq, t_out_of_ice, t_hit_sat])

    # @info "min_t = $min_t; i_min_t = $i_min_t; Δt = $Δt | t_out_of_liq = $t_out_of_liq; t_out_of_ice = $t_out_of_ice; t_hit_sat = $t_hit_sat"

    if min_t < Δt
        if i_min_t == 1
            @debug "liq will run out first before timestep is over... will transition at t = $(min_t) to just ice decay if ice is present"
            S_ql = -q_liq / min_t
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, q_sl, q_si)
            Δt_left = Δt - min_t

            # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, min_t)
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t)

            dδ = smallest_magnitude(-(-q_liq + S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)

            
            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q, q_liq, q_ice, -q_liq, S_qi * min_t) # use multiplied form for floating point accuracy
            # new_regime = get_regime_err_WBF(new_q_vap, new_q.liq, new_q.ice, q_sl, q_si, T, T_freeze) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        elseif i_min_t == 2
            @debug "ice will run out first before timestep is over... will transition at t = $(min_t) to just liq decay if liq is present"
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t)
            S_qi =  -q_ice / min_t
            Δt_left = Δt - min_t

            # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, min_t)
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t)

            dδ = smallest_magnitude(-(S_ql * min_t - q_ice), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            new_δ_0 = δ_0 + dδ
            new_δ_0i = δ_0i + dδ
            new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)


            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q, q_liq, q_ice, S_ql * min_t, -q_ice) # use multiplied form for floating point accuracy
            # new_regime = get_regime_err_WBF(new_q_vap, new_q.liq, new_q.ice, q_sl, q_si, T, T_freeze) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            new_regime = add_regime_parameters(Subsaturated, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            # S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(Subsaturated{true, false, BF}(true, false, BF), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1)
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        else # i_min_t == 3
            @debug "will hit sat before timestep is over... will transition to wbf at t = $(min_t) if liq is present otherwise we're stuck..."
            S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t)
            S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, q_sl, q_si)
            Δt_left = Δt - min_t

            # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, min_t)
            S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, min_t)

            dδ = smallest_magnitude(-(S_ql * min_t + S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            if BF
                new_δ_0 = δ_0 + dδ
                new_δ_0i = FT(0) # hit ice sat from below so δ_0i = 0
                new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            else
                new_δ_0 = FT(0) # hit liq sat from below so δ_0 = 0
                new_δ_0i = δ_0i + dδ
                new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)
            end
            

            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_exponential_part_only(q, q_liq, q_ice, S_ql, S_qi, min_t)
            # new_regime = get_regime_err_WBF(new_q_vap, new_q.liq, new_q.ice, q_sl, q_si, T, T_freeze) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            # S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(WBF{true, true, BF}(true, true, BF), param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1)
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end
    else
        @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt)
        S_qi = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, Δt, q_sl, q_si)
        # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql, S_qi = S_below_ceiling_EPA(S_ql, S_qi, q_liq, q_ice, δ_0, δ_0i, Δt)
        S_ql, S_qi = clamp_S(S_ql, S_qi, regime, δ_0, δ_0i, q_liq, q_ice, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end


end


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated liq but no ice
    - no ice so liq shrink
    - if below freezing, would hit ice_sat first and stop
    - if above freezing, would hit liq_sat first and stop
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{true, false, true}, Subsaturated{true, false, false}}, # can run out of liq first. if not and we make it to ice sat, transition to WBF
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    if depth ≥ 10
        @info "Failed on inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts"
        error("Failed to converge after 10 iterations")
    end

    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    BF::Bool = T < T_freeze
    @debug "Calling Subsaturated{true, false, $BF}..."


    A_c = A_c_func_no_WBF_EPA(FT) # just FT(0)
    τ = τ_liq # we have no ice so ice can't grow or shrink while subsaturated

    t_out_of_liq = get_t_out_of_q_liq_EPA(δ_0, A_c, τ, τ_liq, q_liq) # Eq C6
    t_hit_sat = BF ? t_δ_hit_value(q_si-q_sl, δ_0, A_c, τ) : t_δ_hit_value(FT(0), δ_0, A_c, τ) # below freezing, stop at ice_sat which is lower, above freezing, stop at liq_sat which is lower


    min_t, i_min_t = find_min_t([t_out_of_liq, t_hit_sat])

    if min_t < Δt
        if i_min_t == 1
            @debug "liq will run out first at t = $(min_t) before timestep is over..."
            S_ql = -q_liq / Δt # can't continue once have nothing
            return return_mixing_ratio ? (S_ql, FT(0)) : (S_mixing_ratio_to_shum(S_ql, q.tot), FT(0))
        else # i_min_t == 2
            @debug "will hit sat before timestep is over... will transition at t = $(min_t) to wbf"
            S_ql = S_ql_func_indiv_EPA( A_c, τ_liq, δ_0, min_t)
            S_qi = FT(0)
            Δt_left = Δt - min_t

            # S_ql = S_above_floor(S_ql, q_liq, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, min_t)
            S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, min_t)

            dδ = smallest_magnitude(-(S_ql * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            if BF
                new_δ_0 = δ_0 + dδ
                new_δ_0i = FT(0) # hit ice sat from below so δ_0i = 0
                new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)
            else
                new_δ_0 = FT(0) # hit liq sat from below so δ_0 = 0
                new_δ_0i = δ_0i + dδ
                new_δ_0 = clamp_δ(new_δ_0, regime, q_sl, q_si)
            end

            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t)
            # new_regime = get_regime_err_WBF(new_q_vap, new_q.liq, new_q.ice, q_sl, q_si, T, T_freeze) # more robust against floating point errors and such, if slightly slower... e.g. underflow in calculating get_q_out_of_q that leads to larger dt than truly needed
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, BF)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            S_ql *= (min_t / Δt) # rescale to the timestep
            S_qi *= (min_t / Δt) # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end
    else
        @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_indiv_EPA( A_c, τ_liq, δ_0, Δt)
        S_qi = FT(0)
        # S_ql = S_above_floor(S_qi, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_ql = S_below_ceiling_EPA(S_ql, q_ice, δ_0, Δt)
        S_ql = clamp_S_ql(S_ql, regime, δ_0, q_liq, q_ice, Δt)

        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated ice but no liq 
    - no ice so liq shrink
    - below freezing, would hit ice_sat first and stop
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{false, true, true}}, # can run out of ice first. then we're stuck.
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    # (q.liq < eps(FT)) || error("shouldn't have any liq here (up to floating point error)..., got q_liq = $(q.liq)")
    # (q.ice > FT(0)) || error("should have some ice here (up to floating point error)..., got q_ice = $(q.ice)")

    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    @debug "Calling Subsaturated{false, true, true}..."

    # (q_vap < (BF ? q_sl : q_si)) || error("should be subsaturated, q_vap = $q_vap, q_vap_sat = $(BF ? q_sl : q_si)")

    τ = τ_ice # we have no liq

    A_c = A_c_func_no_WBF_EPA(FT) # just FT(0)

    # can either lose all the ice or hit ice sat first, either way we're done
    t_out_of_ice = get_t_out_of_q_no_WBF_EPA(δ_0i, A_c, τ, τ_ice, q_ice) # just like QCCON but w/ only ice
    t_hit_sat = t_δ_hit_value(FT(0), δ_0i, A_c, τ) # below freezing, stop at ice_sat which is lower, above freezing, stop at liq_sat which is lower (just like QCCON but now w/ only ice)

    min_t, i_min_t = find_min_t([t_out_of_ice, t_hit_sat])

    # @info "min_t = $min_t, i_min_t = $i_min_t | t_out_of_ice = $t_out_of_ice, t_hit_sat = $t_hit_sat"

    if min_t < Δt # either way it's the same up to 
        if i_min_t == 1
            @debug "ice will run out first before timestep is over at t = $(min_t)..."
            S_ql = FT(0)
            S_qi = -q_ice / Δt # we're done so just scale to the entire timestep
        else
            @debug "will hit sat before timestep is over... Then we're done bc we don't have liq in BF and if above freezing, we can't form ice."
            S_ql = FT(0)
            S_qi = S_qi_func_indiv_EPA( A_c, τ_ice, δ_0i, min_t)
            # S_qi = S_above_floor(S_qi, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_qi = S_below_ceiling_EPA(S_qi, q_liq, δ_0i, min_t)
            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t)
            S_qi *= min_t / Δt # rescale to the timestep
        end
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    else
        @debug "nothing of note through end of timestep..."
        S_ql = FT(0)
        S_qi = S_qi_func_indiv_EPA( A_c, τ_ice, δ_0i, Δt)
        # S_qi = S_above_floor(S_qi, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_qi = S_below_ceiling_EPA(S_qi, q_liq, δ_0i, Δt)
        S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end


"""
    Subsaturated ice but no liq 
    - no ice so liq shrink
    - above freezing, would hit liq_sat first and then go to WBF, losing ice while gaining liq
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Subsaturated{false, true, false}, # can run out of ice first. then we're stuck.
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    # (q.liq < eps(FT)) || error("shouldn't have any liq here (up to floating point error)..., got q_liq = $(q.liq)")
    # (q.ice > FT(0)) || error("should have some ice here (up to floating point error)..., got q_ice = $(q.ice)")

    # (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)
    (; q_sl, q_si, q_vap, q_liq, q_ice, T_freeze) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

    @debug "Calling Subsaturated{false, true, false}..."

    τ = τ_ice # we have no liq

    A_c = A_c_func_no_WBF_EPA(FT) # just FT(0)

    # can either lose all the ice or hit ice sat first, either way we're done
    t_out_of_ice = get_t_out_of_q_no_WBF_EPA(δ_0i, A_c, τ, τ_ice, q_ice) # just like QCCON but w/ only ice
    t_hit_sat =  t_δ_hit_value(q_sl-q_si, δ_0i, A_c, τ) # below freezing, stop at ice_sat which is lower, above freezing, stop at liq_sat which is lower (just like QCCON but now w/ only ice)

    min_t, i_min_t = find_min_t([t_out_of_ice, t_hit_sat])

    # @info "min_t = $min_t, i_min_t = $i_min_t | t_out_of_ice = $t_out_of_ice, t_hit_sat = $t_hit_sat"

    if min_t < Δt # either way it's the same up to 
        if i_min_t == 1
            @debug "ice will run out first before timestep is over at t = $(min_t)..."
            S_ql = FT(0)
            S_qi = -q_ice / Δt # we're done so just scale to the entire timestep
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        else
            @debug "will hit liq sat before timestep is over... Then we're done bc we don't have liq in BF and if above freezing, we can't form ice."
            S_ql = FT(0)
            S_qi = S_qi_func_indiv_EPA( A_c, τ_ice, δ_0i, min_t)
            # S_qi = S_above_floor(S_qi, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            # S_qi = S_below_ceiling_EPA(S_qi, q_liq, δ_0i, min_t)
            S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, min_t)

            Δt_left = Δt - min_t

            dδ = smallest_magnitude(-(S_qi * min_t), dδ_func_EPA(A_c, τ, δ_0, min_t)) # dδ_func_EPA can be bad bc of nextfloat problems w/ A_c, τ creation etc.
            new_δ_0 = FT(0) # hit liq sat from below so δ_0 = 0
            new_δ_0i = δ_0i + dδ
            new_δ_0i = clamp_δi(new_δ_0i, regime, q_sl, q_si)

            new_q, new_q_vap = morrison_milbrandt_2015_get_new_status_helper_EPA(q, q_liq, q_ice, S_ql, S_qi, min_t)
            new_regime = add_regime_parameters(WBF, new_q.liq, new_q.ice, false)
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style_exponential_part_only(new_regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, new_q_vap, new_q, q_eq, Δt_left, ts; use_fix = use_fix, return_mixing_ratio = false, depth = depth+1, δ_0 = new_δ_0, δ_0i = new_δ_0i)
            S_qi *= (min_t / Δt) # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end
    else
        S_ql = FT(0)
        S_qi = S_qi_func_indiv_EPA( A_c, τ_ice, δ_0i, Δt)
        @debug "nothing of note through end of timestep..."
        # S_qi = S_above_floor(S_qi, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        # S_qi = S_below_ceiling_EPA(S_qi, q_liq, δ_0i, Δt)
        S_qi = clamp_S_qi(S_qi, regime, δ_0i, q_liq, q_ice, Δt)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
"""
    Subsaturated with nothing
    - return 0, 0
"""
function morrison_milbrandt_2015_style_exponential_part_only(
    regime::Union{Subsaturated{false, false, true}, Subsaturated{false, false, false}}, # kind of null, you have nothing and can't make anything
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::Real, ts::TD.ThermodynamicState;
    use_fix::Bool = true, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it)
    return_mixing_ratio::Bool = false,
    depth::Int = 0,
    δ_0::FT,
    δ_0i::FT,
    ) where {FT}

    return FT(0), FT(0) 
end


