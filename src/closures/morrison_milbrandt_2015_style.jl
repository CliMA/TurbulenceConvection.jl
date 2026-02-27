"""
# Milbrandt's equations are in mixing ratio not specific humidity, to convert we need to use

If q_tot = all_water / (all_water + dry_air) then 1 - q_tot =  dry_air/(all_water + dry_air) 

So then anything that is  [[ x / (all_water + dry_air) ]]  can be converted to  [[x / dry_air]]  by multiplying by [[ (all_water + dry_air) / dry_air ]] = [[ 1 / (1 - q_tot) ]]

    So anything that is specific humidity [[ x / (all_water + dry_air) ]]  can be converted using [[ x / (1 - q_tot) ]]
    And, then anything that is mixing ratio [[ x / dry_air ]]  can be converted back using [[ x * (1 - q_tot) ]]

this is good for us though bc q_tot is a constant, so derivatives are easier...
So while you might think 
    d/dy(q_spe) = d/dy(q_mix / (1 + q_mix)) = d/dy(q_mix) / (1 + q_mix)^2 unlesl q_spe = q_tot this isn't true...
Instead we can just say d/dy(q_spe_cond) = d/dy(q_spe_cond / (1 - q_tot)) = d/dy(q_mix_cond) / (1 - q_tot) bc q_tot is a constant.
Likewise d/dy(q_mix_cond) = d/dy(q_spe_cond * (1 - q_tot)) = d/dy(q_spe_cond) * (1 - q_tot)
    where _spe = specific humidity, _mix = mixing ratio, _cond = condensate q, _tot = total q

If what you have is w_tot = all_water / dry_air, observe that  (1 + w_tot) = [[ (dry_air + all_water) / dry_air]] = 1 / (1 - q_tot)

    So anything that is specific humidity [[ x / (all_water + dry_air) ]] can be converted using [[ x * (1 + w_tot) ]]
    And, then anything that is mixing ratio [[ x / dry_air ]] can be converted back using [[ x / (1 + w_tot) ]]

# q_spe = q_mix / (1 + q_mix)   [ in reality this is more like q_spe_c = q_mix_c / (1 + q_tot_mix) , bc 1+q_tot = total ratio all water / dry air...]
# q_mix = q_spe / (1 - q_spe)   [ in reality this is more like q_mix_c = q_spe_c / (1 - q_tot_spe) , bc 1-q_tot = dry mixing ratio...]

# [NOTE: q_tot appears instead of just q because we have more than one water species so we need to be careful...]

# note then d/dy (q_spe) = d/dy(q_mix / (1 + q_mix)) = d/dy(q_mix) / (1 + q_mix)^2  [ really it's more like d/dy(q_spe_c) = d/dy(q_mix_c / (1 + q_tot_mix)) and q_tot_mix should be a constant under phase transformation... ]
# note then d/dy (q_mix) = d/dy(q_spe / (1 - q_spe)) = d/dy(q_spe) / (1 - q_spe)^2  [ really it's more like d/dy(q_mix_c) = d/dy(q_spe_c / (1 - q_tot_spe)) and q_tot_spe should be a constant under phase transformation... ]

NOTE: These mixing ratio, specific humidity, and their derivatives are essentially identical at earthlike conditions...
"""



"""
    Because q_vap is calculated from qt, q_liq, and q_ice, it is possibly smaller
    Similarly, q_vap - q_sl, q_vap - q_si can be smaller than floating point precision on q_liq, q_ice, q_vap which we do track...
    Thus, you can construct states that are indistinguishabe from being at liq/ice equilibrium in that you cannot adjust q_liq or q_ice to get a different q_vap... 
    In those cases, we just take δ_0 or δ_0i to be 0 as necessary...
    This would somewhat be alleviated by truncating q_liq, q_ice, Δt, and τs to be >>  eps(FT) but that's not a universal fix...
    It would always be possible that (qt-qc)-q_sc < nextfloat(qc), as the former is at smallest nextfloat(qt)-q_sc which is at smallest O(nextfloat(nextfloat(qc))).
    The solution would be to track supersat explicitly and let ql, qi do what they want but still it would be impossible to construct q states that conserve qt and match this.

    [[ This can probably be deprecated nowadays? idk... in MM2015() it still might be needed because of our `close enough` transitions with tolerance]]
"""
function limit_δ(δ::FT, q_vap::FT, q_liq::FT, q_ice::FT) where {FT}
    if δ > FT(0)
        if (δ < (nextfloat(q_liq) - q_liq)) || (δ < (nextfloat(q_ice) - q_ice)) || (-δ > (prevfloat(q_vap) - q_vap)) # no amount of condensate gain or vapor loss can create δ=0
            return FT(0)
        end
    elseif δ < FT(0)
        if (δ > (prevfloat(q_liq) - q_liq)) || (δ > (prevfloat(q_ice) - q_ice)) || (-δ < (nextfloat(q_vap) - q_vap)) # no amount of condensate loss or vapor gain can create δ=0
            return FT(0)
        end
    end
    return δ
end




function get_params_and_go_to_mixing_ratio(param_set::APS,
    area::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_vap::FT, # should be shum...
    dqvdt::FT,
    dTdt::FT,
    δ_0_shum::FT,
    δ_0i_shum::FT,
    q::TD.PhasePartition,
    q_eq::TD.PhasePartition,
    Δt::FT,
    ts::TD.ThermodynamicState;
    use_fix::Bool = false, # i think something is wrong with the forumula for large timesteps... is less essential now w/ the limiter but still needed... (unless I just don't understand the physics of it),
    ) where {FT}

    g = TCP.grav(param_set) # acceleration of gravity

    thermo_params = TCP.thermodynamics_params(param_set) # currently repeated in both places, pare down later

    L_i = TD.latent_heat_sublim(thermo_params, ts) # Latent heat for ice (L_i)
    L_l = TD.latent_heat_vapor(thermo_params, ts)  # Latent heat for water
    c_p = TD.TP.cp_d(thermo_params) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?
    # c_p = TD.cp_m(thermo_params, q) # specific heat of air, they just say `Specific heat of air at constant pressure` so i assume moist?

    e_sl = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) # saturation vapor pressure, or TD.partial_pressure_vapor(thermo_params, p, q_eq) 
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice()) 

    # analytical forms for change in saturation vapor presure with temperature | Note these are concave down functions when ρ is not constant but p is instead [ See Fig 3: https://doi.org/10.1029/2009RG000302 ]
    dqsl_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(1), T, q_eq.liq) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix... is it a good enough approx?
    dqsi_dT = TD.∂q_vap_sat_∂T(thermo_params, FT(0), T, q_eq.ice) # how do we make this work differently for ice/liquid? or do we need to? currently it's based on the  current condensate mix...

    # dqsl_dT /= (1 - q_eq.liq)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.liq = q_sl at T
    # dqsi_dT /= (1 - q_eq.ice)^2 # convert derivative from specific humidity to mixing_ratio, q_eq.ice = q_si at T

    dqsl_dT = TD.shum_to_mixing_ratio(dqsl_dT, q.tot) # convert to mixing ratio (testinggg)
    dqsi_dT = TD.shum_to_mixing_ratio(dqsi_dT, q.tot) # convert to mixing ratio (testinggg)

    # So, we can do everything in mixing ratio and then convert back at the end instead of trying to convert everythig in their derivation to figure out what q_tot_mix
    # q_vap = TD.shum_to_mixing_ratio(q_vap, q.tot)
    q_vap = TD.shum_to_mixing_ratio(TD.vapor_specific_humidity(q), q.tot) # convert to mixing ratio (use this version instead of q_vap bc q_vap passed from other functions into this fcn can be in mix form already than )
    q_liq = TD.shum_to_mixing_ratio(q.liq, q.tot)
    q_ice = TD.shum_to_mixing_ratio(q.ice, q.tot)
    q_sl = TD.shum_to_mixing_ratio(q_eq.liq, q.tot) # q_eq is the equilibrium mixing ratio, q_sl is the saturation mixing ratio, q_eq we made to contain only saturation vapor pressure values...
    q_si = TD.shum_to_mixing_ratio(q_eq.ice, q.tot)

    # T_freeze = TCP.T_freeze(param_set)
    T_freeze = TCP.T_triple(param_set) # In the code this is actually where the saturation vapor pressures are equal... I think it's wrong though

    # δ_0 = q_vap - q_sl # supersaturation over liquid
    # δ_0i = δ_0 + (q_sl - q_si) # supersaturation over ice

    δ_0 = TD.shum_to_mixing_ratio(δ_0_shum, q.tot) # convert to mixing ratio (testinggg)
    δ_0i = TD.shum_to_mixing_ratio(δ_0i_shum, q.tot) # convert to mixing ratio (testinggg)

    δ_0 = limit_δ(δ_0, q_vap, q_liq, q_ice)
    δ_0i = limit_δ(δ_0i, q_vap, q_liq, q_ice)


    # Balanced at large timestep, but needs limiters (otherwise, may consume more liq/ice than exists... but handles WBF directions fine and most stuff for short timesteps
    if use_fix # setting these makes the asymptote ok, not sure if it breaks other things....
        # i think their equation balances liq evap <-> ice growth, but doesn't limit it properly over the timestep, technically q_sl-q_si addition/subtraction should decay as you approach saturation...
        L_i = L_l
        dqsi_dT = dqsl_dT
        # Does evaporating water releasing more latent heat than forming ice takeup risk evaporating a cloud completely?
    end

    Γ_l = 1 + (L_l / c_p) * dqsl_dT  # Eqn C3
    Γ_i = 1 + (L_i / c_p) * dqsi_dT  # Eqn C3

    # dTdt = FT(0) # ignore for now
    # dqvdt = FT(0) # ignore for now

    return (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i, dqvdt, dTdt)
end


τ_func(τ_liq::FT, τ_ice::FT, L_i::FT, c_p::FT, dqsl_dT::FT, Γ_i::FT) where {FT} = 1 / (1 / τ_liq + (1 + (L_i / c_p) * dqsl_dT) * ((1 / τ_ice) / (Γ_i))) # Eqn C2

A_c_func_no_WBF(q_sl::FT, g::FT, w::FT, c_p::FT, e_sl::FT, dqsl_dT::FT, dqvdt::FT, dTdt::FT, p::FT, ρ::FT) where {FT} =
    dqvdt - (q_sl * ρ * g * w) / (p - e_sl) - dqsl_dT * (dTdt - (w * g) / c_p)

A_c_func(τ_ice::FT, Γ_l::FT, q_sl::FT, q_si::FT, g::FT, w::FT, c_p::FT, e_sl::FT, L_i::FT, dqsl_dT::FT, dqvdt::FT, dTdt::FT, p::FT, ρ::FT) where {FT} =
    A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) - (q_sl - q_si) / (τ_ice * Γ_l) * (1 + (L_i / c_p) * dqsl_dT) # Eq C4

# save some cpu cycles
function A_c_func_with_and_without_WBF(τ_ice::FT, Γ_l::FT, q_sl::FT, q_si::FT, g::FT, w::FT, c_p::FT, e_sl::FT, L_i::FT, dqsl_dT::FT, dqvdt::FT, dTdt::FT, p::FT, ρ::FT, CF_mp::FT = one(FT)) where {FT}
    A_c_no_WBF = A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_c = A_c_no_WBF - (q_sl - q_si) / (τ_ice * Γ_l) * (1 + (L_i / c_p) * dqsl_dT) * CF_mp # save some computation
    return (; A_c, A_c_no_WBF) # return both so we can use the no WBF version if we want to, e.g. for large timesteps where it doesn't matter as much
end

"""
https://www.wolframalpha.com/input?i2d=true&i=solve+A+%CF%84%2B%5C%2840%29Subscript%5Bx%2C0%5D-A%CF%84%5C%2841%29+Power%5Be%2C%5C%28123%29+-Divide%5Bt%2C%CF%84%5D%5C%28125%29%5D+%3D+c+for+t
"""
function t_δ_hit_value( value::FT, δ_0::FT, A_c::FT, τ::FT) where {FT}
    logand = ( δ_0 - A_c * τ ) / ( value - A_c * τ )
    # # @debug "logand = $logand; τ = $τ; log(logand) = $(log(logand))"
    return ((logand ≤ 1) || !isfinite(logand)) ? FT(Inf) : τ * log(logand) # solution. if logand is less than 1, the solution is in the past or otherwise unreachable so we return t = Inf
end
"""
solving our ODE
https://www.wolframalpha.com/input?i2d=true&i=solve+A+%CF%84%2B%5C%2840%29Subscript%5Bx%2C0%5D-A%CF%84%5C%2841%29+Power%5Be%2C%5C%28123%29+-Divide%5Bt%2C%CF%84%5D%5C%28125%29%5D+-+Subscript%5Bx%2C0+%5D%3D+c+for+t
"""
function t_dδ_hit_value(value::FT, δ_0::FT, A_c::FT, τ::FT) where {FT}
    logand =  ( A_c * τ - δ_0 ) / ( A_c * τ - value - δ_0 )
    # # @debug "logand = $logand; τ = $τ; log(logand) = $(log(logand))"
    return ((logand ≤ 1) || !isfinite(logand)) ? FT(Inf) : τ * log(logand) # solution. if logand is less than 1, the solution is in the past or otherwise unreachable so we return t = Inf
    # e.g. say A_c is 0 and value > -δ_0, then after removing all sub/supersat you'll be stuck at sat and will never reach value...
end

function δ_func(δ_0::FT, A_c::FT, τ::FT, t::FT) where {FT} 
    # return A_c * τ + (δ_0 - A_c * τ) * exp(-t/ τ) # Eqn C1
    # term = (1-exp(-t/τ))
    term = -expm1(-t/τ) # this is the same as 1 - exp(-t/τ) but avoids the 0 * Inf problem, which is what we want to avoid
    return  δ_0 + (A_c * τ * term) # more precise, as the change in δ will be accurate even when (-Δt / τ) is small. important if δ_0 is small as well.
end


# Go from mixing ratio world back to specific humidity world
mixing_ratio_to_shum(q_mix::FT, q_tot::FT) where{FT} = q_mix * (1 - q_tot) # undo shum_to_mixing_ratio(), bc it's q_tot, it's not w/(1+w)
# function S_mixing_ratio_to_shum(S_mix::FT, q_tot::FT) where {FT} # q_tot has no difference bc it's tot
#     # q_tot_mix = q_tot / (1 - q_tot) # maybe faster in case it doesnt get inlined... idk.
#     # @inline q_tot_mix = TD.shum_to_mixing_ratio(q_tot, q_tot) # probably gets inlined
#     # return S_mix / (1 + q_tot_mix) # bc we showed  d/dt(q_spe) = d/dt(q_mix) / (1+q_tot_mix) under phase transformation, which makes sense, it's smaller
#     return mixing_ratio_to_shum(S_mix, q_tot) # bc we showed  d/dt(q_spe) = d/dt(q_mix) / (1+q_tot_mix) under phase transformation, which makes sense, it's smaller
# end
const S_mixing_ratio_to_shum = mixing_ratio_to_shum # this is actually equivalent to the above definiton. see /home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES/scaling_funcs.jl for explanation... it's bc q_tot is a constant
const S_shum_to_mixing_ratio = TD.shum_to_mixing_ratio # make it a const to avoid recompilation

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# maybe we should change the default to false... it's an extra slow check and I think we've worked out all the bugs... Furthermore, if you use say get_t_out_of_q[_WBF] but the answer is too small, you may just get t=0 and have to discard that solution. We have no way to fix that, so just disgard the negative... not ideal but... shouldn't really come up unless numbers are very small
function q_new(q_old::TD.PhasePartition, S_ql_mix_Δt::FT, S_qi_mix_Δt::FT, dqvdt_Δt::FT; error_on_q_neg::Bool = true) where {FT}
    new_q_liq = q_old.liq + S_mixing_ratio_to_shum(S_ql_mix_Δt, q_old.tot)
    new_q_ice = q_old.ice + S_mixing_ratio_to_shum(S_qi_mix_Δt, q_old.tot)

    if error_on_q_neg && (new_q_liq < 0 || new_q_ice < 0)
        error("Negative q detected, q_old = $q_old, S_ql_mix_Δt = $S_ql_mix_Δt, S_qi_mix_Δt = $S_qi_mix_Δt, new_q_liq = $new_q_liq, new_q_ice = $new_q_ice")
    end

    return TD.PhasePartition(q_old.tot + dqvdt_Δt, new_q_liq, new_q_ice)
end
    
function q_new(q_old::TD.PhasePartition, S_ql_mix::FT, S_qi_mix::FT, Δt::FT, dqvdt::FT; error_on_q_neg::Bool = true) where {FT}
    new_q_liq = q_old.liq + S_mixing_ratio_to_shum(S_ql_mix*Δt, q_old.tot)
    new_q_ice = q_old.ice + S_mixing_ratio_to_shum(S_qi_mix*Δt, q_old.tot)

    if error_on_q_neg && (new_q_liq < 0 || new_q_ice < 0)
        error("Negative q detected, q_old = $q_old, S_ql_mix = $S_ql_mix, S_qi_mix = $S_qi_mix, Δt = $Δt, new_q_liq = $new_q_liq, new_q_ice = $new_q_ice")
    end

    return TD.PhasePartition(q_old.tot + dqvdt * Δt, new_q_liq, new_q_ice)
end

function q_new(q_old::TD.PhasePartition, q_liq_mix::FT, q_ice_mix::FT, S_ql_mix_Δt::FT, S_qi_mix_Δt::FT, dqvdt_Δt::FT; error_on_q_neg::Bool = true) where {FT}
    new_q_liq = mixing_ratio_to_shum(q_liq_mix + S_ql_mix_Δt, q_old.tot)
    new_q_ice = mixing_ratio_to_shum(q_ice_mix + S_qi_mix_Δt, q_old.tot)

    if error_on_q_neg && (new_q_liq < 0 || new_q_ice < 0)
        error("Negative q detected, q_old = $q_old, q_liq_mix = $q_liq_mix, q_ice_mix = $q_ice_mix, S_ql_mix_Δt = $S_ql_mix_Δt, S_qi_mix_Δt = $S_qi_mix_Δt, new_q_liq = $new_q_liq, new_q_ice = $new_q_ice")
    end

    return TD.PhasePartition(q_old.tot + dqvdt_Δt, new_q_liq, new_q_ice)
end

function q_new(q_old::TD.PhasePartition, q_liq_mix::FT, q_ice_mix::FT, S_ql_mix::FT, S_qi_mix::FT, Δt::FT, dqvdt::FT; error_on_q_neg::Bool = true) where {FT}
    new_q_liq = mixing_ratio_to_shum(q_liq_mix + S_ql_mix*Δt, q_old.tot)
    new_q_ice = mixing_ratio_to_shum(q_ice_mix + S_qi_mix*Δt, q_old.tot)

    if error_on_q_neg && (new_q_liq < 0 || new_q_ice < 0)
        error("Negative q detected, q_old = $q_old, q_liq_mix = $q_liq_mix, q_ice_mix = $q_ice_mix, S_ql_mix = $S_ql_mix, S_qi_mix = $S_qi_mix, Δt = $Δt, new_q_liq = $new_q_liq, new_q_ice = $new_q_ice")
    end

    return TD.PhasePartition(q_old.tot + dqvdt * Δt, new_q_liq, new_q_ice)
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

function morrison_milbrandt_2015_get_new_status_helper_helper(param_set::APS, p::FT, new_q::TD.PhasePartition, q_eq::TD.PhasePartition, ts::TD.ThermodynamicState, Δt::FT, w::FT) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set) # currently repeated in both places, pare down later

    new_q_vap = TD.vapor_specific_humidity(new_q)
    θ = TD.liquid_ice_pottemp(thermo_params, ts) # should be conserved, as is p
    ρ = TD.air_density(thermo_params, ts)
    g = TCP.grav(param_set)

    new_p = p - w * ρ * g * Δt # new pressure [] # Use w = - ω/(ρ*g) --> ω = - w * ρ * g | new_p = p + ω * Δt = p - w * ρ * g * Δt ]
    new_ts = TD.PhaseNonEquil_pθq(thermo_params, p, θ, new_q) # new starting state (sould we conserve ρ or p? probably p bc it's fixed on the grid...)
    new_T = TD.air_temperature(thermo_params, new_ts) # -- note this should really account for the latent heating and dTdt terms...
    new_ρ = TD.air_density(thermo_params, new_ts)
    new_q_sl = TD.q_vap_saturation_generic(thermo_params, new_T, new_ρ, TD.Liquid())
    new_q_si = TD.q_vap_saturation_generic(thermo_params, new_T, new_ρ, TD.Ice())
    new_q_eq = TD.PhasePartition(new_q.tot, new_q_sl, new_q_si)

    # do our conversions after creating everything
    new_q_vap = TD.shum_to_mixing_ratio(new_q_vap, new_q.tot) # convert back to specific humidity
    new_q_sl = TD.shum_to_mixing_ratio(new_q_sl, new_q.tot) # convert back to specific humidity
    new_q_si = TD.shum_to_mixing_ratio(new_q_si, new_q.tot) # convert back to specific humidity

    return (; new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts)
end

function morrison_milbrandt_2015_get_new_status_helper_full_source(param_set::APS, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, S_ql_mix_Δt::FT, S_qi_mix_Δt::FT, Δt::FT, dqvdt_Δt::FT, ts::TD.ThermodynamicState, w::FT; error_on_q_neg::Bool = true) where {FT}
    new_q = q_new(q, S_ql_mix_Δt, S_qi_mix_Δt, dqvdt_Δt; error_on_q_neg = error_on_q_neg)
    return morrison_milbrandt_2015_get_new_status_helper_helper(param_set, p, new_q, q_eq, ts, Δt, w)
end

function morrison_milbrandt_2015_get_new_status_helper(param_set::APS, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, S_ql_mix::FT, S_qi_mix::FT, Δt::FT, dqvdt::FT, ts::TD.ThermodynamicState, w::FT; error_on_q_neg::Bool = true) where {FT}
    new_q = q_new(q, S_ql_mix, S_qi_mix, Δt, dqvdt, dqvdt; error_on_q_neg = error_on_q_neg)
    return morrison_milbrandt_2015_get_new_status_helper_helper(param_set, p, new_q, q_eq, ts, Δt, w)
end

function morrison_milbrandt_2015_get_new_status_helper_full_source(param_set::APS, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq_mix::FT, q_ice_mix::FT, S_ql_mix_Δt::FT, S_qi_mix_Δt::FT, Δt::FT, dqvdt_Δt::FT, ts::TD.ThermodynamicState, w::FT; error_on_q_neg::Bool = true) where {FT} # This version should be more floating point stable -- whatever you calculate in mixing ratio you should get exactly out.
    new_q = q_new(q, q_liq_mix, q_ice_mix, S_ql_mix_Δt, S_qi_mix_Δt, dqvdt_Δt; error_on_q_neg = error_on_q_neg) # use safe version
    return morrison_milbrandt_2015_get_new_status_helper_helper(param_set, p, new_q, q_eq, ts, Δt, w)
end

function morrison_milbrandt_2015_get_new_status_helper(param_set::APS, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq_mix::FT, q_ice_mix::FT, S_ql_mix::FT, S_qi_mix::FT, Δt::FT, dqvdt::FT, ts::TD.ThermodynamicState, w::FT; error_on_q_neg::Bool = true) where {FT} # This version should be more floating point stable -- whatever you calculate in mixing ratio you should get exactly out.
    new_q = q_new(q, q_liq_mix, q_ice_mix, S_ql_mix, S_qi_mix, Δt, dqvdt; error_on_q_neg = error_on_q_neg) # use safe version
    return morrison_milbrandt_2015_get_new_status_helper_helper(param_set, p, new_q, q_eq, ts, Δt, w)
end


""" One species τ_species = τ, no WBF and no competition """
function S_func_indiv_no_WBF( A_c::FT, τ::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    if isfinite(τ)
        term_1 = (δ_0 - A_c * τ) / (Δt * Γ)
        # term_2 =  (1 - exp(-Δt / τ))
        term_2 = -expm1(-Δt / τ) # this is the same as 1 - exp(-Δt / τ) but avoids the 0 * Inf problem, which is what we want to avoid
        if isinf(term_1) || isinf(term_2)
            return floatmax(FT) * sign(term_1) * sign(term_2)
        end
        prod_ = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem
        return A_c / Γ + prod_
        # return A_c / Γ + (δ_0 - A_c * τ)  / (Δt * Γ) * (1 - exp(-Δt / τ))   # QCCON EQN C6 [ did I deprecate this to create the lower structure? need to document better ]
    else
        return FT(0) # I know it looks like it should be A_c bc τ and τ_c canceled, but infinite τ means no change... Really what this is saying is τ_c == τ == τ_<other> = ∞ , the limit at τ -> ∞ of τ/τ_c goes to 0 not 1. It's a simplification of S_func_no_WBF().
    end
end


""" One of two species i.e. τ_species ≠ τ, no WBF """
function S_func_no_WBF(A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    if isfinite(τ_c) # can't have τ_c finite and τ infinite, so we don't need to check τ (at least from properly constructed τ values.)
        term_1 = (δ_0 - A_c * τ) * (τ /(τ_c * Γ))/(Δt)
        if isinf(term_1)
            return floatmax(FT) * sign(term_1)
        end
        # term_2 =  (1 - exp(-Δt / τ))
        term_2 = -expm1(-Δt / τ) # for precision
        if isinf(term_2) # can happen if Δt is neg e.g. from secant method
            return floatmax(FT) * sign(term_2)
        end
        prod = iszero(term_2) ? FT(0) : (term_1 * term_2) # term_1 * term_2 returning Inf is ok bc we avoid the 0 * Inf problem

        preterm = A_c * τ / (τ_c * Γ)

        # if !isfinite(term_1) || !isfinite(term_2) || !isfinite(preterm) || !isfinite(prod)
        #     # @debug "A value is bad, got: term_1 = $term_1; term_2 = $term_2; preterm = $preterm; prod = $prod"
        #     # @debug "A_c = $A_c; τ = $τ; τ_c = $τ_c; δ_0 = $δ_0; Δt = $Δt; Γ = $Γ"
        #     error("Bad value detected: term_1 = $term_1; term_2 = $term_2; preterm = $preterm; prod = $prod")
        # end

        return A_c * τ / (τ_c * Γ) + prod
        # return A_c * τ / (τ_c * Γ) + (δ_0 - A_c * τ) * τ / (Δt * τ_c * Γ) * (1 - exp(-Δt / τ))   # QCCON EQN C6 [ did I deprecate this to create the lower structure? need to document better ]
    else  # if τ_c is infinite, then either τ is infinite or it's some other value from the other condensate, either way this output is 0 all the same.
        return FT(0) # avoid NaN problems at τ_c = ∞
    end
end

const S_ql_func = S_func_no_WBF
const S_ql_func_indiv = S_func_indiv_no_WBF

S_func_WBF(A_c::FT, τ::FT, τ_ice::FT, δ_0::FT, Δt::FT, Γ_i::FT, q_sl::FT, q_si::FT, CF_mp::FT = one(FT)) where {FT} =
    S_func_no_WBF( A_c, τ, τ_ice, δ_0, Δt, Γ_i)  + 
    ((q_sl - q_si) / (τ_ice * Γ_i)) * CF_mp #  # QICON Eqn C7 # if τ_ice is infinite, this works just fine

const S_qi_func = S_func_WBF
const S_qi_func_indiv = S_func_indiv_no_WBF
const S_qi_func_no_WBF = S_func_no_WBF  

# issmallt(Δt::FT) where {FT} = abs(Δt) < eps(FT) # is Δt small? we only use this on Δt, maybe we should just use Δt < eps(FT) instead to prevent putting in any negative numbers
# issmallt(Δt::FT) where {FT} = Δt < eps(FT) # is Δt small? we only use this on Δt, maybe we should just use Δt < eps(FT) instead to prevent putting in any negative numbers
issmallt(Δt::FT) where {FT} = iszero(Δt) # if issmall(Δt), then x is small and in the safe version of the functions, we shouldn't calculate S which might be unreasonably large ... up to Inf if Δt is 0 is 0
limit_large_Δt(Δt::FT) where {FT} = isinf(Δt) ? floatmax(FT) : Δt # if Δt is infinite, then we can't calculate S, so we return the largest possible value
limit_large_Δt(Δt::FT, limit::FT) where {FT} = min(Δt, limit) # if Δt is infinite, then we can't calculate S, so we return the largest possible value

# These methods won't fail and return NaN when Δt is 0 |  I think with const we don't have to make the unions here... 
S_func_zero_Δt_safe(S_func::Union{typeof(S_func_indiv_no_WBF), typeof(S_ql_func_indiv), typeof(S_qi_func_indiv)},
    A_c::FT, τ::FT, δ_0::FT, Δt::FT, Γ::FT; small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt) where {FT} = small_func(Δt) ? FT(0) : S_func(A_c, τ, δ_0, Δt, Γ)
    
S_func_zero_Δt_safe(S_func::Union{typeof(S_func_no_WBF), typeof(S_ql_func), typeof(S_qi_func_no_WBF)},
    A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT, Γ::FT; small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt) where {FT} = small_func(Δt) ? FT(0) : S_func(A_c, τ, τ_c, δ_0, Δt, Γ)

S_func_zero_Δt_safe(S_func::Union{typeof(S_func_WBF), typeof(S_qi_func)},
    A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT, Γ::FT, q_sl::FT, q_si::FT; small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt) where {FT}  = small_func(Δt) ? FT(0) : S_func(A_c, τ, τ_c, δ_0, Δt, Γ, q_sl, q_si)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

"""
see https://www.wolframalpha.com/input?i2d=true&i=solve+A*t+%2B+%5C%2840%29Subscript%5Bx%2C0%5D-A%CF%84%5C%2841%29+%5C%2840%291-Power%5Be%2C%5C%28123%29+-Divide%5Bt%2C%CF%84%5D%5C%28125%29%5D%5C%2841%29+%3D+c+for+t

e.g. 
    Rearrange -q_liq = QCCON * Δt to get: A_c * Δt + (δ_0 - A_c * τ) * (1 - exp(-Δt / τ)) = c where c = -q_liq * (τ_liq * Γ_l / τ)  --> RHS = q_liq

    If you have no liquid processes at all, you can also use this function, replacing δ_0 with δ_0i and Γ and τ with their ice counterparts and using A_c_no_WBF. This is because the QICON formluation in say subsat w/o liquid is just exponential decay w/ no WBF term.

    if the sol'n t is too close to zero, due to floating point summation error we won't be able to find it...
"""
function get_t_out_of_q_no_WBF(δ_0::FT, A_c::FT, τ::FT, τ_c::FT, q_c::FT, Γ::FT, exit_if_fail::Bool=true) where {FT}
    c = -q_c * (τ_c * Γ / τ) # Rearrange QCCON * Δt to get: A_c * Δt + (δ_0 - A_c * τ) * (1 - exp(-Δt / τ)) = c and solve for Δt

    t_out_of_q = FT(Inf)
    # sol'n 1, 2
    if !iszero(A_c) # required for both
        if (δ_0 ≠ (A_c * τ)) # sol'n 1 is valid [almost a surety]
            # if (lambert_argument = (exp(-(c + A_c * τ - δ_0) / (A_c * τ)) * (δ_0 - A_c * τ) / (A_c * τ))) > -exp(-1)
            # if (lambert_argument =  exp(-(c + A_c * τ - δ_0) / (A_c * τ) + log(abs(δ_0 - A_c * τ)) - log(abs(A_c * τ))) * sign(δ_0 - A_c * τ) * sign(A_c * τ)) > -exp(-1) # the log subtract loses accuracy...
            if (lambert_argument =  exp(-(c + A_c * τ - δ_0) / (A_c * τ) +  log(abs((δ_0 - A_c * τ)/(A_c * τ)))) * sign(δ_0 - A_c * τ) * sign(A_c * τ)) > -exp(-1) # this is more accurate...
                # this can underflow... note that LambertW.lambertw.(-(big(10)^-1e10), (0,-1)) is (-1.000...004e-10000000000, -2.3025...e+10)! so the exponentiation will go wrong long before the lambertw function... thus getting a time of zero can happen even if the time is only something like .0003
                if exit_if_fail && (iszero(lambert_argument) || isinf(lambert_argument)) && !iszero(q_c)
                    # note `exit_if_fail` gives us an out, BigFloat can still undreflow, ( see [ https://github.com/JuliaLang/julia/issues/17893 ] and [ https://stackoverflow.com/a/38825652/13331585 ] ), so we want an out then
                    # right now, if we retry with BigFloat, we are setting exit_if_fail to false... You could make an argument we should just pass the nothing up the chain (i.e. return nothing immediately) and then fallback to :standard or something in microphysics_coupling_limiters.jl that will converge.
                    # return nothing # Note, the time is probably not zero. when we go through find_min_t() it'll get raised to eps(FT) but the time is almost certainly longer than that! But fast is fast so hopefully it doesnt cause too many issues. We don't want to have to go to BigFloats (or even ArbFloats https://github.com/JeffreySarnoff/ArbNumerics.jl)
                    return (; sol = t_out_of_q, success = false)
                else
                    for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                        t_candidate = (A_c * τ * LambertW.lambertw(lambert_argument, LW_branch) +  A_c*τ + c - δ_0) / A_c # for some reason this seemed more stable in paractice....
                        # t_candidate = τ * (LambertW.lambertw(lambert_argument, LW_branch) + 1) + (c-δ_0) / A_c # true if δ_0 ≠ A_c [dividing by A_c usually makes things much bigger, so we separate the tau scale stuff from the potentially much smaller ], we done do (τ + 1) though in case τ is large [[ updte wasn't actually more stable in testing... idk why... ]]
                        if t_candidate > FT(0)
                            # # @debug ("q here 1, t_candidate = $t_candidate, lambert_argument = $lambert_argument")
                            t_out_of_q = min(t_out_of_q, t_candidate)
                        elseif iszero(t_candidate) # could be an underflow problem... (i.e. nextfloat is not big enough for the addition and subtraction above to work...)
                            # t_out_of_q = FT(0) # if we're too close to zero, we can't find the solution
                        end
                    end
                end

            end
        else # if δ_0 = A_c * τ, sol'n 2 is valid
            # # @debug "here"
            t_candidate = c/A_c
            if (t_candidate = c/A_c) > FT(0)
                t_out_of_q = min(t_out_of_q, c / A_c)  # invalid unless A * τ = δ_0.. [A_c ≠ 0] ... unlikely, need A_c = δ_0 / τ which is generally not true
            end
        end
    end

    # sol'n 3
    if iszero(δ_0) # only valid if δ_0 = 0, which is possible right at a perfect transition
        if (lambert_argument = -exp(-c/(A_c*τ) - 1)) > -exp(-1)
            if exit_if_fail && iszero(lambert_argument) && !iszero(q_c)
                return (; sol = t_out_of_q, success = false)
            else
                for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                    t_candidate = τ * LambertW.lambertw(lambert_argument, LW_branch) + c/(A_c) + τ # true if δ_0 = 0 which is possible at transition
                    if t_candidate > FT(0)
                        # # @debug("q here 2, t_candidate = $t_candidate, lambert_argument = $lambert_argument")
                        t_out_of_q = min(t_out_of_q, t_candidate)
                    elseif iszero(t_candidate) # could be an underflow problem... (i.e. nextfloat is not big enough for the addition and subtraction above to work...)
                        # t_out_of_q = FT(0) # if we're too close to zero, we can't find the solution
                    end
                end
            end
        end
    end

    # sol'n 4
    # there is another solution that is also valid if A_c * τ = δ_0, but that requires c = 0 implying q_liq == 0, which would also mean we're not in this regime.

    # sol'n 5
    if iszero(A_c) && ((δ_0 * c) > c^2) # only for [δ_0 > c if c > 0] or [δ_0 < c if c < 0] (i.e. δ_0 * c > c^2) do we always get a positive answer. I believe our c is always negative btw but this is most general
        # # @debug("q here 3")
        t_out_of_q = min(t_out_of_q, τ * log(δ_0/(δ_0-c)) )  # invalid unless A_c = 0, unlikely...
    end

    return (; sol = t_out_of_q, success = true)
end
const get_t_out_of_q_liq = get_t_out_of_q_no_WBF
const get_t_out_of_q_ice_no_WBF = get_t_out_of_q_no_WBF

"""
see https://www.wolframalpha.com/input?i2d=true&i=solve+A*t+%2B+%5C%2840%29Subscript%5Bx%2C0%5D-A%CF%84%5C%2841%29+%5C%2840%291-Power%5Be%2C%5C%28123%29+-Divide%5Bt%2C%CF%84%5D%5C%28125%29%5D%5C%2841%29+%2B+B*t+%3D+c+for+t
e.g. 
    Rearrange -q_ice = QICON * Δt to get: A_c * Δt +  B * Δt +  (δ_0 - A_c * τ) * (1 - exp(-Δt / τ)) = c where c = (-q_ice)) * (τ_ice * Γ_i / τ) --> RHS = q_ice
    where B = (q_sl - q_si) / (τ_ice * Γ_i) * (τ_ice * Γ_i / τ) 
"""
# function get_t_out_of_q_WBF(δ_0::FT, A_c::Union{FT, FT2}, τ::FT, τ_c::FT, q_ice::FT, Γ::FT, q_sl::FT, q_si::FT, exit_if_fail::Bool=true) where {FT, FT2}
function get_t_out_of_q_WBF(δ_0::FT, A_c::FT2, τ::FT, τ_c::FT, q_ice::FT, Γ::FT, q_sl::FT, q_si::FT, exit_if_fail::Bool=true, CF_mp::FT = one(FT)) where {FT, FT2}
    B = ((q_sl - q_si) /  τ) * CF_mp # (q_sl - q_si) / (τ_c * Γ) * (τ_c * Γ / τ) #
    c = -q_ice * (τ_c * Γ / τ) # Rearrange QICON * Δt = -q_ice to get: A_c * Δt + (δ_0 - A_c * τ) * (1 - exp(-Δt / τ)) + B * Δt = c and solve for Δt


    t_out_of_q = FT(Inf)

    # sol'n 1,2
    if !iszero(A_c + B)
        if (δ_0 ≠ (A_c * τ)) # sol'n 1 is valid [almost a surety]
            # if (lambert_argument = exp(-(c + A_c * τ - δ_0) / ((A_c + B)*τ)) * (δ_0 - A_c * τ) / ((A_c + B)*τ)) > -exp(-1) # within lambertw domain
            # lambert_argument_big = exp(-(c_big + A_c_big * τ_big - big(δ_0)) / ((A_c_big + B_big)*τ_big)) * (big(δ_0) - A_c_big * τ_big) / ((A_c_big + B_big)*τ_big) # for BigFloat stability
            if (lambert_argument = exp(-(c + A_c*τ - δ_0)/((A_c + B)*τ) + log(abs((δ_0 - A_c*τ)/((A_c + B)*τ))) * sign(δ_0 - A_c*τ) * sign((A_c + B)*τ))) > -exp(-1) # this is more accurate...
                if exit_if_fail && iszero(lambert_argument) && !iszero(q_ice)
                    # return nothing # signifiy a problem
                    return (; sol = t_out_of_q, success = false)
                else    
                    for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                        t_candidate = (τ*(A_c + B) * LambertW.lambertw(lambert_argument, LW_branch) +  A_c*τ + c - δ_0) / (A_c + B) # # for some reason this seemed more stable in paractice....
                        # t_candidate = τ * LambertW.lambertw(lambert_argument, LW_branch) +  (A_c*τ + c - δ_0) / (A_c + B) # trial version
                        # @debug("ice here 1, t_candidate = $t_candidate, lambert_argument = $lambert_argument")
                        if t_candidate > FT(0)
                            # # @debug("ice here 1, t_candidate = $t_candidate, lambert_argument = $lambert_argument")
                            t_out_of_q = min(t_out_of_q, t_candidate)
                        end
                    end
                end
            end
            return (; sol = t_out_of_q, success = true)


        else # if δ_0 = A_c * τ, sol'n 2 is valid
            # @debug("ice here 2")
            if (t_candidate = c / (A_c + B)) > FT(0)
                t_out_of_q = min(t_out_of_q, t_candidate)  # invalid unless A * τ = δ_0.. [A_c ≠ 0] ... unlikely, need A_c = δ_0 / τ which is generally not true
            end
        end
    end


    # sol'n 3, 4
    if iszero(δ_0)
        if iszero(A_c) 
            if !iszero(B*τ) # really is just iszero(B) bc τ shouldnt be 0
                # @debug("ice here 3")
                if (t_candidate = c / B) > FT(0)
                    t_out_of_q = min(t_out_of_q, c / B) 
                end
            end

        end

        if !iszero(A_c) && !iszero(A_c + B) && !iszero(c)
            if (lambert_argument = - A_c * exp(-(c + A_c * τ) / ((A_c + B)*τ)) / (A_c + B)) > -exp(-1) # within lambertw domain
                if exit_if_fail && iszero(lambert_argument) && !iszero(q_ice)
                    # return nothing # signifiy a problem (we can't know sol'n for sure, you can try fallback to BigFloat)
                    return (; sol = t_out_of_q, success = false)
                else
                    for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                        t_candidate = (τ*(A_c + B) * LambertW.lambertw(lambert_argument, LW_branch) + A_c * τ + c) / (A_c + B)
                        if t_candidate > FT(0)
                            # # @debug("ice here 5, t_candidate = $t_candidate")
                            t_out_of_q = min(t_out_of_q, t_candidate)
                        end
                    end
                end
            end
        end
    end

    # sol'n 5
    if (B == -A_c) && !iszero(A_c * τ) && ((logand = (A_c * τ - δ_0) / (A_c * τ + c - δ_0)) > FT(1)) # only for logand > 1 do we get a positive answer
        # # @debug("ice here 6")
        t_out_of_q = min(t_out_of_q, τ * log(logand)) # invalid unless A_c = 0, unlikely...
    end

    # return t_out_of_q
    return (; sol = t_out_of_q, success = true)
end

const get_t_out_of_q_ice = get_t_out_of_q_WBF

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

"""
    dT/dt = L_l/c_p * δ/(τ_L*Γ_l) + L_i/c_p * (δ+(q_sl-q_si))/((τ_I*Γ_I)) + (-g*w/c_p) +(dT/dt)_mix + (dT/dt)_rad)
    We can use QCCON for  δ/(τ_L*Γ_l) and QICON for (δ+(q_sl-q_si))/((τ_I*Γ_I)) , then dT/dt = L_l/c_p * QCCON + L_i/c_p * QICON + (-g*w/c_p) +(dT/dt)_mix + (dT/dt)_rad
    The solution then looks like the solution for q_vap gone, but now:
    K =  τ / (τ_liq * Γ_l) * L_l/c_p + τ / (τ_ice * Γ_i) * L_i/c_p
    B = (q_sl - q_si) /  τ * L_i/c_p + (-g*w/c_p) + (dT/dt)_mix + (dT/dt)_rad
    c = T_freeze - T

    WARNING - This function doesn't working  -- just use iterative solver.... 
    I think it's because unlike with q_liq, q_ice, saturation, temperature has no analytic solution.
    There's second order effects we can't get precisely in dqsl/dT,  dqsi/dT as well as the fact that the equations are writen from liquid perspectiv,e favoring dqsl/dT.
"""

function get_t_T_hit_T_freeze(δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, T::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, g::FT, w::FT, c_p::FT, L_l::FT, L_i::FT, T_freeze::FT, dTdt::FT, exit_if_fail::Bool=true, CF_mp::FT = one(FT)) where {FT}
    B = ((q_sl - q_si) / τ * L_i/c_p) * CF_mp + -g*w/c_p + dTdt
    K = τ / (τ_liq * Γ_l) * L_l/c_p + τ / (τ_ice * Γ_i) * L_i/c_p
    c = T_freeze - T 

    t_out_of_q = FT(Inf)
    # sol'n 1
    if !iszero(A_c * K + B)
        if (δ_0 ≠ (A_c * τ)) # sol'n 1 is valid [almost a surety]
            if (lambert_argument = exp(-(c + A_c*K*τ - K*δ_0) / ((A_c*K + B)*τ)) * K*(δ_0 - A_c * τ) / ((A_c*K + B)*τ)) > -exp(-1) # within lambertw domain
                if exit_if_fail && iszero(lambert_argument) && !iszero(c)
                    # return nothing # signifiy a problem (we can't know sol'n for sure, you can try fallback to BigFloat)
                    return (; sol = t_out_of_q, success = false)
                else
                    for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                        t_candidate = (τ*(A_c*K + B) * LambertW.lambertw(lambert_argument, LW_branch) +  A_c*K*τ + c - K*δ_0) / (A_c*K + B)
                        # @debug "branch is $LW_branch, t_candidate = $t_candidate"
                        if t_candidate > FT(0)
                            # @debug("vap here 1, t_candidate = $t_candidate, lambert_argument = $lambert_argument")
                            t_out_of_q = min(t_out_of_q, t_candidate)
                        end
                    end
                end
            end
        else # if δ_0 = A_c * τ, sol'n 2 is valid
            # @debug("vap here 2")
            t_out_of_q = min(t_out_of_q, c / (A_c * K + B))  # invalid unless A * τ = δ_0.. [A_c ≠ 0] ... unlikely, need A_c = δ_0 / τ which is generally not true
        end    
    end

    # sol'n 3
    if iszero(δ_0) && !iszero(A_c*K)
        if !iszero(A_c*K + B) 
            if (lambert_argument = -A_c*K * exp(-(c + A_c*K*τ)/((B + A_c*K)*τ)) / (B + A_c*K)) > -exp(-1) # within lambertw domain
                if exit_if_fail && iszero(lambert_argument) && !iszero(c)
                    # return nothing # signifiy a problem (we can't know sol'n for sure, you can try fallback to BigFloat)
                    return (; sol = t_out_of_q, success = false)
                else
                    for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                        t_candidate = (τ*(A_c*K + B) * LambertW.lambertw(lambert_argument, LW_branch) + A_c*K*τ + c)/(A_c*K + B)
                        if t_candidate > FT(0)
                            # @debug("vap here 3, t_candidate = $t_candidate")
                            t_out_of_q = min(t_out_of_q, t_candidate)
                        end
                    end
                end
            end
        end

        if iszero(c) 
            if (lambert_argument = -A_c*K * exp(-(A_c*K)/(A_c*K + B)) / (A_c*K + B)) > -exp(-1) # within lambertw domain
                if exit_if_fail && iszero(lambert_argument) && !iszero(c)
                    # return nothing # signifiy a problem (we can't know sol'n for sure, you can try fallback to BigFloat)
                    return (; sol = t_out_of_q, success = false)
                else
                    for LW_branch in ((lambert_argument < FT(0)) ? (0, -1) : (0,))
                        t_candidate = τ*(LambertW.lambertw(lambert_argument, LW_branch) + ((A_c*K) / (A_c*K + B)))
                        if t_candidate > FT(0)
                            # @debug("vap here 4, t_candidate = $t_candidate")
                            t_out_of_q = min(t_out_of_q, t_candidate)
                        end
                    end
                end
            end
        end
    end

    # sol'n 4
    if ((A_c*B*τ) ≠ (B*δ_0)) && ((A_c*c + B*δ_0) ≠ (A_c*B*τ)) && (K == B/A_c) && !iszero(A_c)
        # @debug("vap here 5")
        t_out_of_q = min(t_out_of_q, τ * log(B*(A_c*τ - δ_0)/(A_c*B*τ - A_c*c - B*δ_0))) # invalid unless A_c = 0, unlikely...
    end

    # sol'n 5
    if iszero(A_c) && iszero(B) && !iszero(K*δ_0)
        # @debug("vap here 6")
        t_out_of_q = min(t_out_of_q, τ * log((K*δ_0)/(K*δ_0 - c))) # invalid unless A_c = 0, unlikely...
    end

    # @debug "get_t_T_hit_T_freeze is gonna return $t_out_of_q"
    # return t_out_of_q
    return (; sol = t_out_of_q, success = true)
    
end
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

residual_tol(x::FT; factor::FT = 1.0) where{FT} = max(eps(FT), abs(nextfloat(x) - x)) * factor
target(value::FT, add_or_sub::Union{typeof(+),typeof(-)}; factor::FT = 1.0) where{FT} = add_or_sub(value , residual_tol(value, factor = factor))
# const TF = FT(1e0) # the default target / tolerance factor we're using throughout the code.
const TF = 0 # now that we are explicitly setting δ_0, δ_0i, I don't think we need to make use of the target factor anymore... They're bad anyway if you can't reach that point e.g. for very samll inputs.
# We don't need a buffer really anymore, now we just solve for t from get_t_var_hit_value and then just set δ_0, δ_0i based on that value

"""
    NOTES:   
    
    Consider making value = target_value ± tolerance to ensure you pass your target value including with floating point error.

    Given we have the bounds, and often a good guess, the order from fastest to slowest should be: Newton's Method (fastest, but requires derivative) -> NewtonsMethodAD (haven't tested) -> Secant Method -> Regula Falsi 

    We do not permit Newton's Methods for now, as I don't have derivates and don't want to use autodiff for AD
    Technically, Secant is not and Regula-Falsi is guaranteed to converge, but we're constrained enough I trust Secant (Thermodynamics also defaults to it)

    also returns new_T for quick BF checking... can maybe remove this since we're gonna start checkign T switches...

    Note that while sum of S_ql and S_qi is not necessarily monotonic, We can only cross a saturation line once in a timestep.
    This system as written is stateless -- it only depends on our parameters and q, q_eq, etc... The external forcings are constant and not time dependent.
    Thus the derivative at any state is deterministic, so there is no way to cross a line twice, which would require different derivatives for the same inputs.
    [NOTE: while dδ_0/dt, QCCON = dql_dt, and QICON are all monotonic, they can cross 0 and thus S_ql * Δt, S_qi * Δt may not be monotonic]

    However, q_sl - q_si does change, so that can change if we hit ice saturation (that wouldn't matter if we used ice perspective.... but we'd have the same problem again in WBF)
    The best we can do is just try and hope for convergence, if we don't converge just assume we never hit sat

    Our attempted methodology is:
    - Check the endpoint. if the sign is the same as the start, there may not be a root in the region. If using Regula-Falsi, fallback to secant.
    - Check the predicted point. If we know have a root, we can bisect our search region; if not we've learned nothing.
    - If using Regula-Falsi, bisect the search region N times to see if we find a sign change. This is slow and not ideal but it works surprisingly often... Need a better way around this
    - Search for the root w/ requested method, unless that method was Regula-Falsi and we have no sign change, in which case we fallback to Secant.

    We use residual tolerances so you can specifcy precisely say your e.g. saturation target is eps(FT) and residual tolerance is eps(FT) so that you're guaranteed a positive result, here between [0, 2eps(FT)]

    Note again we do NOT use this method for finding when  we run out of liquid/ice bc S_ql, S_qi are precise, so the t we get out are analytically correct -- it's just predicting saturation/temperature that's hard due to estimating dqsl_dT, dqsi_dT, w effects etc.

        NOTE: For all the work we do, one might consider checking if something like IntervalRootFinding.jl is faster. It certainly should be more guaranteed to find a root if it exists, though it also searchces for all roots which may slow it down.
"""

function get_t_var_hit_value(var::Union{Val{:δ}, Val{:δi}, Val{:δip}, Val{:δipl}, Val{:T}}, param_set::APS, value::FT, min_t::FT, max_t::FT, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, w::FT, p::FT, q_vap::FT, dqvdt::FT, dTdt::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, ts::TD.ThermodynamicState,
    ;
    solver_method::SM = TD.RS.RegulaFalsiMethod, # something like RS.SecantMethod or TD.RS.RegulaFalsiMethod, default to TD.RS.RegulaFalsiMethod bc if it converges it converges inbounds
    residual_tol_factor::FT = FT(TF), # default to using this factor bc it's the same we set the value targets with, guaranteeing a range from [0 to 2*tolerance] with a target at 1xtolerance ensuring correct sign.
    max_iter::Int = 1000,
    num_check_bisect_endpoints::Int = 10, # how many times cut the endpoint in half and retry before quitting
    check_predicted_t::Bool = true, # could cut down the region a lot e.g. if the point is very early..., maybe saves only one iter.
    still_try_secant_method_if_endpoint_signs_match::Bool = true,
    allow_secant_fallback_to_regula_falsi::Bool = true, # it would be nice to allow regula-falsi to fallback to secant as well since it only works when signs are different at the end, but you risk an infinite loop...
    allow_regula_falsi_fallback_to_secant::Bool = true,
    use_WBF::Bool = true, # if false, we use the no WBF version of everything. Very useful for when you're say in the single species regime.
    ) where {FT, SM <: Type{X} where {X <: TD.RS.RootSolvingMethod{RSFT} where RSFT} }

    # It would be nice to store all the return values from our function, but we can only return f(t) for rootfinding, so we will just have to live with those values getting recalculated (i know, wasteful.)

    # dont use this, instead we use let block to avoid penalty from capturing variables  (though it had seemed ok without let block using anonymous (NOT named) func, naming the func really killed our inference times, up from 10s to 5mins)
    # function tester_func(t::FT) where {FT}
    #     out = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF) - value
    #     # isnan(out) && error("out is nan at t = $t for $var")
    #     return out
    # end

    known_root::Bool = false

    # check endpoint to know if to continue (only continue if signs match -- root may still exist even if they do but the only way we can search is secant which can diverge outside the range and waste our time...
    if min_t > FT(0)
        var_at_min_t = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, min_t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF)
    else # shorcut to save some evaluation at t = 0, these better be correc though or you'll have bugs (as in return the same thing get_t_var_hit_value_helper_Δt_safe would return)
        var_at_min_t = begin
            if (var isa Val{:δ}) || (var isa Val{:δip}) # supersat over liquid or supersat over ice bc we're doing thigns from ice perspective.
                δ_0
            elseif var isa Val{:δi}
                δ_0 + (q_sl - q_si) # supersat over ice...
            elseif var isa Val{:δipl}
                δ_0 - (q_sl - q_si) # supersat over <f ill in here >
            elseif var isa Val{:T}
                T
            else
                error("var must be Val{:δ}(), Val{:δi}(), Val{:δip}(), Val{:δipl}(), or Val{:T}() for now")
            end
        end
        # (var_at_min_t == get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, min_t, ts, w; error_on_q_neg = false, use_WBF = use_WBF)) || error("var_at_min_t != get_t_var_hit_value_helper_Δt_safe(var, min_t)") # bug checking
    end

    var_at_max_t = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, max_t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF)
    # if the endpoint is at the value, we can't find the root or the endpoint is the root (if we're already at the root at min_t, this is also a fine thing to just continue to the end...)
    if iszero(var_at_max_t - value)
        return max_t, true
    end
    
    # if the starting point is 0, that's fine, we could be asking for a return to that point... Note 0 has signbit false (same as positive numbers), so going from 0 to negative would claim to be a root... however 0*x would still fail RegulaFalsiMethod()'s negative answer check, Secant would just return 0 again.
    if iszero(var_at_min_t - value)
        min_t += max(eps(FT), nextfloat(min_t) - min_t) # Add a little and re-solve, this can help us ensure we find the next 0 and not the one we're already at
        var_at_min_t = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, min_t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF)
    end
    
    # sometimes, there is a root even if the endpoints match. This just means there are multiple (an even number) roots. Bisection search can help us find some of these cases.
    if signbit(var_at_max_t - value) == signbit(var_at_min_t - value) # there's no root in the interval, could also use and 
        # # @debug(" possibly no root in interval, got var_at_max_t = $var_at_max_t, var_at_min_t = $var_at_min_t, value = $value | (var_at_max_t - value) = $(var_at_max_t - value), (var_at_min_t - value) = $(var_at_min_t - value), min_t = $min_t, max_t = $max_t")
        # # @debug "attempting to bisect range to ensure root doesn't exist"
        # bisect search bc sometimes the root is just hard to find...
        new_max_t = max_t
        bisection_failed::Bool = true
        while num_check_bisect_endpoints > 0
            # # @debug "$num_check_bisect_endpoints bisect iters remaining | new_max_t = $new_max_t ==> $(new_max_t/2)"
            new_max_t /= FT(2)
            var_at_new_max_t = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, new_max_t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF)
            # # @debug "(var_at_new_max_t - value) = $(var_at_new_max_t - value) | new_max_t = $new_max_t, var_at_new_max_t = $var_at_new_max_t, value = $value"
            if signbit(var_at_new_max_t - value) != signbit(var_at_min_t - value) # there's a root in the interval, could also use ⊻ (\xor)
                # # @debug " bisecting endpoint worked, found root in interval, var_at_max_t $new_max_t = $var_at_new_max_t, var_at_min_t $min_t = $var_at_min_t, value = $value"
                max_t = new_max_t
                var_at_max_t = var_at_new_max_t
                bisection_failed  = false
                known_root = true
                break
            end
            num_check_bisect_endpoints -= 1
        end

        if bisection_failed
            if (still_try_secant_method_if_endpoint_signs_match && (solver_method == TD.RS.SecantMethod)) # if we're not using secant with adjustment attempt even after endpoint signs don't match, then we're done.
                # do nothing
            elseif (allow_regula_falsi_fallback_to_secant && (solver_method == TD.RS.RegulaFalsiMethod)) # if we're not using secant with adjustment attempt even after endpoint signs don't match, then we're done.
                solver_method = TD.RS.SecantMethod # switch to secant bc endpoint signs match so we can't use regula falsi
                allow_secant_fallback_to_regula_falsi = false # don't allow regula falsi to fallback to secant again
            else # we are using secant but not continuing if endpoint signs match, or we're using regula falsi and we can't fallback to secant
                return max_t, false
            end
        end # otherwise the bisection succeeded and we will continue with whatever method we've chosen.

    else
        known_root = true
    end
    # @debug "get_t_var_hit_value() debug point: $var | var_at_max_t - value = $(var_at_max_t - value), var_at_min_t - value = $(var_at_min_t - value), var_at_max_t = $var_at_max_t, var_at_min_t = $var_at_min_t, value = $value, min_t = $min_t, max_t = $max_t"


    # check predicted t to maybe narrow our range faster (our predicted t should be quite close to the answer normally)
    if check_predicted_t # we calculat these outside anyway, so we should probably turn this off normally
        if (var isa Val{:δ}) || (var isa Val{:δip})
            t_guess = t_δ_hit_value(value, δ_0, A_c, τ)  # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower
        elseif var isa Val{:δi} # ice saturation hits value, so liq hits value + (q_sl - q_si)
            t_guess = t_δ_hit_value(value + (q_sl - q_si), δ_0, A_c, τ)  # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower
        elseif var isa Val{:δipl} # liquid saturation hits value (from ice perspective) so liq its value - (q_sl - q_si) bc δ_0 is liq persepctive's δ_0i
            t_guess = t_δ_hit_value(value - (q_sl - q_si), δ_0, A_c, τ)  # below freezing, stop at ice sat which is lower, above freezing, stop at liq sat which is lower
        elseif (var isa Val{:T})
            t_guess = FT(Inf) # deprecated using get_t_T_hit_T_freeze() bc it was highly inaccurate
        end

        # evaluate what the outcome is when we actually use the prediction
        if t_guess < max_t
            var_at_t_guess = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, t_guess, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF) # this should be an exact function
            if iszero(var_at_t_guess - value) # unlikely but have to heck
                return t_guess, true
            end
            if signbit(var_at_t_guess - value) != signbit(var_at_min_t - value) # there's a root in the interval, could also use ⊻ (\xor)
                # # @debug "var_at_t_guess = $var_at_t_guess, var_at_min_t = $var_at_min_t, value = $value"
                max_t = t_guess # the sign change we found between beginning and end must be between min_t and t_guess
            else # probably not a root could still search w/ secant...
                min_t = t_guess # the sign change we found between beginning and end must be between t_guess and the end
            end
        end # if t_guess > max_t, we could check it but we don't know how far away it is and can't guarantee it's near enough to the root for the methods to be valid...
    end

    # @debug "get_t_var_hit_value() debug point: min_t = $min_t, max_t = $max_t, solver_method = $solver_method"


    # use let block to avoid penalty from capturing variables  (though it had seemed ok without let block using anonymous (NOT named) func, naming the func really killed our inference times, up from 10s to 5mins)
    val_min_t = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, min_t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF) - value
    val_max_t = get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, max_t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF) - value
    if !(val_min_t * val_max_t < FT(0)) && (solver_method == TD.RS.RegulaFalsiMethod)
        # @debug "signs are the same at endpoints, we can't use regula falsi, got val_min_t = $val_min_t, val_max_t = $val_max_t, min_t = $min_t, max_t = $max_t, solver_method = $solver_method"
        # # @debug "var = $var; δ_0 = $δ_0; A_c = $A_c; τ = $τ; τ_liq = $τ_liq; τ_ice = $τ_ice; Γ_l = $Γ_l; Γ_i = $Γ_i; q_sl = $q_sl; q_si = $q_si; T = $T; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; q_liq = $q_liq; q_ice = $q_ice; ts = $ts; w = $w"
        return max_t, false
    end
    sol = let var = var, param_set = param_set, δ_0 = δ_0, A_c = A_c, τ = τ, τ_liq = τ_liq, τ_ice = τ_ice, Γ_l = Γ_l, Γ_i = Γ_i, q_sl = q_sl, q_si = q_si, T = T, p = p, q_vap = q_vap, q = q, q_eq = q_eq, q_liq = q_liq, q_ice = q_ice, ts = ts, w = w
        TD.RS.find_zero(
            t::FT -> get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, t, ts, w; error_on_q_neg = false, use_WBF = use_WBF) - value,
            solver_method(min_t, max_t), # valid form for both secant and regula falsi
            TD.RS.CompactSolution(),
            TD.RS.ResidualTolerance(residual_tol(value; factor = residual_tol_factor)), # Using residual tolerance ensures our eps(FT) roots are safe -- also it ensures we shouldn't need issmallt() I think..., no guarantee i suppose.
            max_iter
        )
    end

    # handle convergence out of bounds caused by secant method
    if !(min_t ≤ sol.root ≤ max_t)
        # @debug "sol.root = $(sol.root) is out of bounds min_t = $min_t, max_t = $max_t, solver_method = $solver_method"
        if (solver_method == TD.RS.SecantMethod) && allow_secant_fallback_to_regula_falsi && known_root
            # @debug "secant method converged out of bounds but we know a root exists (endpoint signs are different...). Falling back to RegulaFalsiMethod() which is guaranteed to converge inbounds..."

            # use let block to avoid penalty from capturing variables  (though it had seemed ok without let block using anonymous (NOT named) func, naming the func really killed our inference times, up from 10s to 5mins)
            sol = let var = var, param_set = param_set, δ_0 = δ_0, A_c = A_c, τ = τ, τ_liq = τ_liq, τ_ice = τ_ice, Γ_l = Γ_l, Γ_i = Γ_i, q_sl = q_sl, q_si = q_si, T = T, p = p, q_vap = q_vap, q = q, q_eq = q_eq, q_liq = q_liq, q_ice = q_ice, ts = ts, w = w, dqvdt = dqvdt, dTdt = dTdt
                TD.RS.find_zero(
                    t::FT -> get_t_var_hit_value_helper_Δt_safe(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q_vap, q, q_eq, q_liq, q_ice, t, ts, w; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = false, use_WBF = use_WBF) - value,
                    TD.RS.RegulaFalsiMethod(min_t, max_t), # valid form for both secant and regula falsi
                    TD.RS.CompactSolution(),
                    TD.RS.ResidualTolerance(residual_tol(value; factor = residual_tol_factor)), # we're not going to get a perfect answer, but we can get close enough
                    max_iter
                )
            end
        else
            # @debug("secant method converged out of bounds -- the root may still exist. Consider setting allow_secant_fallback_to_regula_falsi = true or switching to regula_falsi alltogether. returning Non-Convergence and sol'n at max_t")
            # could use S_func_zero_Δt_safe here but also want new_T for BF checking...
            return max_t, false
        end
    end

    if sol.converged

        # @debug "returning 2 $var| sol.root = $(sol.root)"
        return sol.root, true
    else
        # @debug "returning 3 $var | max_t = $max_t"
        return max_t, false
    end

end



# -------------- #
# If not use_WBF, we call the no_WBf fcn directly on δ_0i if :δ/δ_0i, and on δ_0ipl if :δip/δipl
# -------------- #


function get_t_var_hit_value_helper(::Val{:δ}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    
    if !isfinite(Δt)
        # error("Val{:δ} | Δt is not finite, got Δt = $Δt")
        Δt = floatmax(FT)
    end
    S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
    S_qi = use_WBF ? S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) : S_qi_func_no_WBF( A_c, τ, τ_ice, δ_0 + (q_sl - q_si), Δt, Γ_i)

    (; new_q_vap, new_q_sl) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)

    return new_q_vap - new_q_sl
end

function get_t_var_hit_value_helper_return_outputs(::Val{:δ}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
    S_qi = use_WBF ? S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) : S_qi_func_no_WBF( A_c, τ, τ_ice, δ_0 + (q_sl - q_si), Δt, Γ_i)
    (; new_q_vap, new_q_sl, new_T) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)

    return new_q_vap - new_q_sl, S_ql, S_qi, new_T
end

# ----- #

# from ice perspective
function get_t_var_hit_value_helper(::Val{:δip}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_qi = S_func_no_WBF( A_c, τ, τ_ice, δ_0, Δt, Γ_i)
    S_ql = use_WBF ? S_func_WBF( A_c, τ, τ_liq, δ_0, Δt, Γ_l, q_si, q_sl, CF_mp) : S_func_no_WBF( A_c, τ, τ_liq, δ_0 - (q_si - q_sl), Δt, Γ_l)
    (; new_q_vap, new_q_si) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_q_vap - new_q_si
end

function get_t_var_hit_value_helper_return_outputs(::Val{:δip}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_qi = S_func_no_WBF( A_c, τ, τ_ice, δ_0, Δt, Γ_i)
    S_ql = use_WBF ? S_func_WBF( A_c, τ, τ_liq, δ_0, Δt, Γ_l, q_si, q_sl, CF_mp) : S_func_no_WBF( A_c, τ, τ_liq, δ_0 - (q_si - q_sl), Δt, Γ_l)
    (; new_q_vap, new_q_si, new_T) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_q_vap - new_q_si, S_ql, S_qi, new_T
end

# ----- #

function get_t_var_hit_value_helper(::Val{:δi}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
    S_qi = use_WBF ? S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) : S_qi_func_no_WBF( A_c, τ, τ_ice, δ_0 + (q_sl - q_si), Δt, Γ_i)
    (; new_q_vap, new_q_si) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_q_vap - new_q_si
end

function get_t_var_hit_value_helper_return_outputs(::Val{:δi}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
    S_qi = use_WBF ? S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) : S_qi_func_no_WBF( A_c, τ, τ_ice, δ_0 + (q_sl - q_si), Δt, Γ_i)
    (; new_q_vap, new_q_si, new_T) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_q_vap - new_q_si, S_ql, S_qi, new_T
end

# ----- #

# from ice perspective, so (δ_0l, δ_0) are equivalent to liquid perspective (δ_0, δ_0i)
function get_t_var_hit_value_helper(::Val{:δipl}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = use_WBF ? S_func_WBF( A_c, τ, τ_liq, δ_0, Δt, Γ_l, q_si, q_sl, CF_mp) : S_func_no_WBF( A_c, τ, τ_liq, δ_0 - (q_si - q_sl), Δt, Γ_l) # swap l for i
    S_qi = S_func_no_WBF( A_c, τ, τ_ice, δ_0, Δt, Γ_i)
    (; new_q_vap, new_q_sl) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_q_vap - new_q_sl
end

# from ice perspective, so (δ_0l, δ_0) are equivalent to liquid perspective (δ_0, δ_0i)
function get_t_var_hit_value_helper_return_outputs(::Val{:δipl}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = use_WBF ? S_func_WBF( A_c, τ, τ_liq, δ_0, Δt, Γ_l, q_si, q_sl, CF_mp) : S_func_no_WBF( A_c, τ, τ_liq, δ_0 - (q_si - q_sl), Δt, Γ_l) # swap l for i
    S_qi = S_func_no_WBF( A_c, τ, τ_ice, δ_0, Δt, Γ_i)
    (; new_q_vap, new_q_sl, new_T) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_q_vap - new_q_sl, S_ql, S_qi, new_T
end

# ----- #

function get_t_var_hit_value_helper(::Val{:T}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
    S_qi = use_WBF ? S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) : S_qi_func_no_WBF( A_c, τ, τ_ice, δ_0 + (q_sl - q_si), Δt, Γ_i)
    (; new_T) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_T
end

function get_t_var_hit_value_helper_return_outputs(::Val{:T}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT; dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, use_WBF::Bool = true) where {FT}
    S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
    S_qi = use_WBF ? S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si) : S_qi_func_no_WBF( A_c, τ, τ_ice, δ_0 + (q_sl - q_si), Δt, Γ_i)
    (; new_T) = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, Δt, dqvdt, ts, w; error_on_q_neg = error_on_q_neg)
    return new_T, S_ql, S_qi, new_T # second new_T is for consistency, it's used in other fcns for checking BF automtically
end

# ------- #

# using small_func  = issmallt here is bad bc it doesn't actually fix S_ql, S_qi, you need to also use small_func in the S_ql, S_qi functions
# so if you take the t here and call S_ql, S_qi, you'll get the wrong answer... so we need to fix that in the S_ql, S_qi functions as well

# we could write versions to replace S_ql_func w/ zero safe version or adjust S_ql_func to be zero safe but for now we'll just do this here to skip call to morrison_milbrandt_2015_get_new_status_helper()
# --- #
function get_t_var_hit_value_helper_Δt_safe(var::Val{:δ}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return (q_vap - q_sl)
    else
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

function get_t_var_hit_value_helper_return_outputs_Δt_safe(var::Val{:δ}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return ((q_vap - q_sl), FT(0), FT(0), T)
    else
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper_return_outputs(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

# --- ice persepctive "
function get_t_var_hit_value_helper_Δt_safe(var::Val{:δip}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return (q_vap - q_si)
    else
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

function get_t_var_hit_value_helper_return_outputs_Δt_safe(var::Val{:δip}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
   if small_func(Δt)
        return ((q_vap - q_si), FT(0), FT(0), T)
   else
        Δt = large_func(Δt)
        get_t_var_hit_value_helper_return_outputs(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
   end
end

# --- #
function get_t_var_hit_value_helper_Δt_safe(var::Val{:δi}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return (q_vap - q_si)
    else 
        Δt = large_func(Δt)
        get_t_var_hit_value_helper(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

function get_t_var_hit_value_helper_return_outputs_Δt_safe(var::Val{:δi}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return ((q_vap - q_si), FT(0), FT(0), T)
    else 
        Δt = large_func(Δt)
        get_t_var_hit_value_helper_return_outputs(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

# ---  ice perspective liq
function get_t_var_hit_value_helper_Δt_safe(var::Val{:δipl}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return (q_vap - q_sl)
    else 
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

function get_t_var_hit_value_helper_return_outputs_Δt_safe(var::Val{:δipl}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return ((q_vap - q_sl), FT(0), FT(0), T)
    else
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper_return_outputs(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

# --- #
function get_t_var_hit_value_helper_Δt_safe(var::Val{:T}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return T
    else
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

function get_t_var_hit_value_helper_return_outputs_Δt_safe(var::Val{:T}, param_set::APS, δ_0::FT, A_c::FT, τ::FT, τ_liq::FT, τ_ice::FT, Γ_l::FT, Γ_i::FT, q_sl::FT, q_si::FT, T::FT, p::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, q_liq::FT, q_ice::FT, Δt::FT, ts::TD.ThermodynamicState, w::FT;
    dqvdt::FT = FT(0), dTdt::FT = FT(0), error_on_q_neg::Bool = true, small_func::Union{typeof(issmallt), typeof(iszero)} = issmallt, large_func::typeof(limit_large_Δt) = limit_large_Δt, use_WBF::Bool = true
    ) where {FT}
    if small_func(Δt)
        return (T, FT(0), FT(0), T)
    else
        Δt = large_func(Δt)
        return get_t_var_hit_value_helper_return_outputs(var, param_set, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, p, q, q_eq, q_liq, q_ice, Δt, ts, w::FT; dqvdt=dqvdt, dTdt=dTdt, error_on_q_neg = error_on_q_neg, use_WBF = use_WBF)
    end
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

function morrison_milbrandt_2015_style(
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
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
)::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, fallback_to_standard_supersaturation_limiter, emit_warnings, time_tolerance, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    if area > 0

        if emit_warnings && (Δt < eps(FT))
            # @debug "Timestep $(Δt) is very small (smaller than eps(FT) = $(eps(FT))), may cause numerical issues..."
        end

        δ_0 = q_vap - q_eq.liq # supersaturation over liquid
        δ_0i = δ_0 + (q_eq.liq - q_eq.ice) # supersaturation over ice

        # has_liq::Bool = q.liq > FT(0)
        # has_ice::Bool = q.ice > FT(0)


        # below_freezing::Bool = T < TCP.T_freeze(param_set)
        below_freezing::Bool = T < TCP.T_triple(param_set) # In the code this is actually where the saturation vapor pressures are equal... I think it's wrong though
        regime = get_saturation_regime(q_vap, q, q_eq, below_freezing)

        return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, δ_0_shum = δ_0, δ_0i_shum = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))


        # if below_freezing # [ s_qi < s_ql --> δ_0i > δ_0 ]
        #     if δ_0i < 0 # always evap for both... (only lasts until run out of one...)
        #         return morrison_milbrandt_2015_style(Subsaturated{has_liq, has_ice, true}(has_liq, has_ice, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, use_fix = use_fix, return_mixing_ratio = false, δ_0 = δ_0, δ_0i = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)
        #     elseif (δ_0 < 0 < δ_0i) # WBF (evap and deposition)
        #         return morrison_milbrandt_2015_style(WBF{has_liq, has_ice, true}(has_liq, has_ice, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, use_fix = use_fix, return_mixing_ratio = false, δ_0 = δ_0, δ_0i = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)
        #     else
        #         return morrison_milbrandt_2015_style(Supersaturated{has_liq, has_ice, true}(has_liq, has_ice, true), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, use_fix = use_fix, return_mixing_ratio = false, δ_0 = δ_0, δ_0i = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)
        #     end

        # else # [ s_qi > s_ql --> δ_0i < δ_0 ]
        #     if δ_0 < 0 # always evap for both... (only lasts until run out of one...)
        #         return morrison_milbrandt_2015_style(Subsaturated{has_liq, has_ice, false}(has_liq, has_ice, false), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, use_fix = use_fix, return_mixing_ratio = false, δ_0 = δ_0, δ_0i = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)
        #     elseif (δ_0i < 0 < δ_0) # WBF [ if were symmetric, liq would grow and ice would shrink, so that's fine... don't enforce melting here... ]
        #         return morrison_milbrandt_2015_style(WBF{has_liq, has_ice, false}(has_liq, has_ice, false), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, use_fix = use_fix, return_mixing_ratio = false, δ_0 = δ_0, δ_0i = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)
        #     else
        #         return morrison_milbrandt_2015_style(Supersaturated{has_liq, has_ice, false}(has_liq, has_ice, false), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts, use_fix = use_fix, return_mixing_ratio = false, δ_0 = δ_0, δ_0i = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)
        #     end

        # end
        
    else
        return FT(0), FT(0)
    end
end

# ====================================================================================================================================================================================================== #
# ====================================================================================================================================================================================================== #
# ====================================================================================================================================================================================================== #
@inline min_if_converged(x::FT, x_candidate::FT, x_candidate_converged::Bool) where{FT} = x_candidate_converged ? min(x, x_candidate) : x # helper function to update the shortest time only if converged
@inline S_above_floor(S_qc::FT, q_c::FT, Δt::FT) where {FT} = max(S_qc, -q_c/Δt, -floatmax(FT)) # ensure S_qc does not deplete q_c in floating point
@inline S_above_floor(S_ql::FT, S_qi::FT, q_liq::FT, q_ice::FT, Δt::FT) where {FT} = S_above_floor(S_ql, q_liq, Δt), S_above_floor(S_qi, q_ice, Δt) # ensure S_ql, S_qi do not deplete q_liq, q_ice in floating point
"""
If get_t_out_of liq/ice ... functions return 0, it could be because of underflow, i.e the number was too small to calculate with floating point accuracy.
However, a 0 is bad because it means q/Δt is infinite, which makes it impossible to compose with anything else.
Here we remedy that - given a vector of Δts we return the minimum unless it's zero then we return the next smallest but keep the index the same. 
    - WARNING: This depends on, in the code, you having a separate path for the out of liq or out of ice category! In that path you should set S_ql = q_liq/Δt and S_qi = q_ice/Δt and not use the fucntions that calculate them.
If the next smallest is < eps(FT), we just return eps(FT) in the spirit of the time being very fast
"""
# function find_min_t(Δts::AbstractVector{FT}) where {FT}
#     min_t, i_min_t = findmin(Δts)
#     if iszero(min_t)
#         min_t = minimum(Δts -> (iszero(Δts) ? FT(Inf) : Δts), Δts)  # min of all values that are not zero. Don't change index though. it should be fine as the new_regime will send us to the right place
#         min_t = min(min_t, eps(FT)) # but not further than eps(FT)
#     end
#     return min_t, i_min_t
# end


# function find_min_t(Δts::AbstractVector{FT}) where {FT}
#     min_t, i_min_t = findmin(Δts)

#     if min_t == zero(FT)
#         # Find the smallest non-zero Δt
#         best = FT(typemax(FT))
#         @inbounds for x in Δts
#             if x != zero(FT) && x < best
#                 best = x
#             end
#         end

#         # If everything was zero, best is typemax → clamp to eps(FT)
#         if best == FT(typemax(FT))
#             min_t = eps(FT)
#         else
#             min_t = min(best, eps(FT))
#         end
#     end

#     return min_t, i_min_t
# end



# function find_min_t(Δts::NTuple{_, FT}) where {_, FT} # This should be much more type stable
#     # 1. Initial findmin (highly optimized for NTuple)
#     min_t, i_min_t = findmin(Δts)

#     if iszero(min_t)
#         # 2. Find the minimum non-zero value using the direct, fast loop
        
#         # Initialize the variable to search for the smallest non-zero time
#         min_t_non_zero = FT(Inf) 
        
#         for val in Δts
#             # Compiler unrolls this loop completely for speed
#             if !iszero(val) && val < min_t_non_zero
#                 min_t_non_zero = val
#             end
#         end
        
#         # Use the non-zero minimum, but clamp it at eps(FT)
#         min_t = min(min_t_non_zero, eps(FT))
#     end
    
#     # Returns the potentially adjusted min_t and the original index.
#     return min_t, i_min_t
# end


# function find_min_t(x::FT...) where {FT}
function find_min_t(x::Vararg{FT, _}) where {FT, _} # force specialization [ faster than tuple for small N, and we never go past 4-5]
    # find minimum and its index
    min_t = x[1]
    i_min_t = 1
    for (i, val) in enumerate(x)
        if val < min_t
            min_t = val
            i_min_t = i
        end
    end

    if min_t == zero(FT)
        # smallest non-zero
        best = FT(typemax(FT))
        for val in x
            if val != zero(FT) && val < best
                best = val
            end
        end

        if best == FT(typemax(FT))
            min_t = eps(FT)
        else
            min_t = min(best, eps(FT))
        end
    end

    return min_t, i_min_t
end


# ====================================================================================================================================================================================================== #


"""
    Supersaturated (below freezing)
    - Both can grow!
"""
function morrison_milbrandt_2015_style(
    regime::Union{Supersaturated{true, true, true}, Supersaturated{true, false, true}, Supersaturated{false, true, true}, Supersaturated{false, false, true}}, # these should all work the same, right? you'll end up with some liq/ice at the end no matter what
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT

    # @debug "Calling Supersaturated{$(q.liq > FT(0)), $(q.ice > FT(0)), true}"

    # --- Thermo  constants ------------------------------------------------------------------------------------ #
    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts) 

    # (δ_0 > FT(0)) && (δ_0i > FT(0)) || error(" Should be supersaturated over both phases, got δ_0 = $δ_0, δ_0i = $δ_0i")


    # Both are growing..., no explicit need for WBF right

    A_c = A_c_func(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τ = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i) # Eq C2


    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # don't calculate if we're starting with no liquid, calculation would fail
    (t_out_of_ice, t_out_of_ice_valid) = get_t_out_of_q_ice(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si) # don't calculate if we're starting with no ice, calculation would fail  # call no matter what..., apparently this can be 0 even if you start w/ no ice and gain some eventually you can lose it again? had no idea... evem w no updraft or other external stuff i mean...
    
    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end

        # upgrade to BigFloat Call
        A_c_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
        t_out_of_liq = FT(get_t_out_of_q_liq(big(δ_0), A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_liq = $t_out_of_liq"
    end
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        # upgrade to BigFloat Call
        A_c_big = (@isdefined A_c_big) ? A_c_big : A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        τ_big = (@isdefined τ_big) ? τ_big : τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
        t_out_of_ice = FT(get_t_out_of_q_ice(big(δ_0), A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_ice = $t_out_of_ice"
    end


    
    max_t = min(max_t, t_out_of_liq, t_out_of_ice) #
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(Val{:δ}(), param_set, target(FT(0), -; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(Val{:δ}(), param_set, target(FT(0), -; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)    
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)

    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 
    # @debug "t_out_of_liq = $t_out_of_liq, t_out_of_ice = $t_out_of_ice, t_hit_δ_sat = $t_hit_δ_sat, t_hit_T_freeze = $t_hit_T_freeze | Δt = $Δt"

    min_t, i_min_t = find_min_t(t_out_of_liq, t_out_of_ice, t_hit_δ_sat, t_hit_T_freeze)
    # @debug "min_t = $min_t, i_min_t = $i_min_t"

    if (min_t < Δt)
        if i_min_t == 1 # out of liquid
            S_ql = -q_liq / min_t
            S_qi = S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)
            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, S_qi*min_t, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # floating point safe

            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liquid before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= min_t / Δt
            S_qi *= min_t / Δt
            S_ql_addit *= Δt_left / Δt
            S_qi_addit *= Δt_left / Δt
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        elseif i_min_t == 2 # out of ice
            S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = - q_ice / min_t
            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, S_ql*min_t, -q_ice, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # floating point safe
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of ice before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        
            S_ql *= min_t / Δt
            S_qi *= min_t / Δt
            S_ql_addit *= Δt_left / Δt
            S_qi_addit *= Δt_left / Δt
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        else # i_min_t == 3,4 #
            S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)
            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            (i_min_t == 3) && # @debug "Hit liq saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            (i_min_t == 4) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            if i_min_t == 3 # hit liq saturation
                new_δ_0_shum = FT(0)
                new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice # not sure what temp is so need to be sure, could use new_q_sl - new_q_si but that's not always correct depending on new T
            else # hit freezing
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are equal at freezing
            end
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end
    else
        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        S_qi = S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si)
        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))

    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Supersaturated (above freezing)
    - No ice growth, liquid only.
    - Existing ice remains, do not enforce melting here.
"""
function morrison_milbrandt_2015_style(
    regime::Union{Supersaturated{true, true, false} , Supersaturated{true, false, false} , Supersaturated{false, true, false} , Supersaturated{false, false, false} },# these should all work the same, right? you'll end up with some liq/ice at the end no matter what
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    # @debug "Calling Supersaturated{$(q.liq > FT(0)), $(q.ice > FT(0)), false}"


    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT


    τ_ice_here = FT(Inf) # no ice/decay above freezing
    τ = τ_liq # no ice/decay above freezing
    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice_here, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)

    A_c = A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) # supersaturated so no loss and now growth because above freezing


    # no ice changing above freezing
    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # don't calculate if we're starting with no liquid, calculation would fail (actually call it anyway, i think we ignore the t = 0 sol'n)
    
    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            q_vap = TD.vapor_specific_humidity(q)
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        # upgrade to BigFloat Call
        A_c_big = A_c_func_no_WBF(big(q_sl), big(g), big(w), big(c_p), big(e_sl), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # no WBF above freezing
        t_out_of_liq = FT(get_t_out_of_q_liq(big(δ_0), A_c_big, big(τ), big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_liq = $t_out_of_liq"
    end

    
    max_t = min(max_t, t_out_of_liq) #
    # t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(Val{:δ}(), param_set, target(FT(0), +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    # above freezing so ice sat is first right?
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(Val{:δi}(), param_set, target(FT(0), +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, -; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    max_t = min_if_converged(max_t, t_hit_T_freeze, t_hit_T_freeze_converged)


    # @error(" If we're above freezing, shouldn't we be hitting ice sat first if we're coming from supersaturated")
    
    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 
    
    min_t, i_min_t = find_min_t(t_out_of_liq, t_hit_δ_sat, t_hit_T_freeze)

    if (min_t < Δt)
        if i_min_t == 1 # out of liquid
            S_ql = - q_liq / min_t
            S_qi = FT(0)
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, FT(0), min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # floating point safe
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liquid before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt
            S_qi *= min_t / Δt
            S_ql_addit *= Δt_left / Δt
            S_qi_addit *= Δt_left / Δt
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        else # i_min_t == 2,3 hit sat or freezing
            S_ql = S_ql_func_indiv( A_c, τ_liq, δ_0, min_t, Γ_l)
            S_qi = FT(0) # no ice growth above freezing

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            (i_min_t==2) && # @debug "Hit saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            (i_min_t==3) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            if i_min_t == 2 # hit ice saturation
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = FT(0)
            else # hit freezing
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are equal at freezing
            end
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= (min_t / Δt) # rescale to the timestep
            S_qi *= (min_t / Δt) # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end        
    else
        # exponential decay, so will never hit liq sat ( should we switch to just normal func with τ_liq = Inf?)
        S_ql = S_ql_func_indiv( A_c, τ_liq, δ_0, Δt, Γ_l)
        S_qi = FT(0)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF has liquid (below freezing)
    - liq can shrink, ice can grow
"""
function morrison_milbrandt_2015_style(
    regime::Union{WBF{true, true, true} ,  WBF{true, false, true}},# has liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    # @debug "Calling WBF{true, $(q.ice > FT(0)), true}"


    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT
    
    local use_ice_perspective::Bool    

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)    
    
    # (δ_0i > FT(0) > δ_0) || error(" WBF but got δ_0 = $δ_0, δ_0i = $δ_0i, q_vap = $q_vap, q_sl = $q_sl, q_si = $q_si")

    A_cL = A_c_func(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τL = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i) # Eq C2

    A_cI = A_c_func(τ_liq, Γ_i, q_si, q_sl, g, w, c_p, e_si, L_l, dqsi_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τI = τ_func(τ_ice, τ_liq, L_l, c_p, dqsi_dT, Γ_l)

    if (A_cL ≤ FT(0)) && (A_cI ≥ FT(0))
        # @debug "You're probably staying in the WBF regime... unclear which perspective is better. Maybe if A_cL * τL > δ_0  (i.e. moving towards liq sat) use liq perspective? and otherwise use ice?"
        if A_cL > δ_0
            use_ice_perspective = false # we're moving towards liq sat, use liq perspective
        else
            use_ice_perspective = true # we're moving towards ice sat, use ice perspective
        end
    elseif (A_cL ≥ FT(0)) && (A_cI ≤ FT(0))
        error("Not sure how compatible asymtotes of supersat over liquid and subsat over ice are...")
    elseif (A_cL ≤ FT(0)) && (A_cI ≤ FT(0))
        # @debug "You 100% are ending up subsaturated... since below freezing maybe use ice perspective?"
        use_ice_perspective = true
    elseif (A_cL ≥ FT(0)) && (A_cI ≥ FT(0))
        # @debug "You 100 % are ending up supersaturated... since above freezing maybe use liq perspective"
        use_ice_perspective = false
    end

    #= 
        Here we shouldn't set use_WBF to false bc we're not moving down to a no WBF (single changing variable) setup, we're just swapping ice and liquid perspective completely 
        In ice perspective we should, however, use δi and δipl inplace of δ and δi. we set δ and δo accordingly
    =#
    if use_ice_perspective
        # @debug "Using ice perspective"
        A_c = A_cI
        τ = τ_func(τ_ice, τ_liq, L_l, c_p, dqsi_dT, Γ_l)  # flip into ice perspective
        δ_0o = δ_0  # other 
        δ_0 = δ_0i
        self = :ice
        other = :liq # just for printing, can deprecate
    else
        # @debug "Using liq perspective"
        A_c = A_cL
        τ = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i) 
        δ_0o = δ_0i # other
        self = :liq
        other = :ice # just for printing, can deprecate
    end

    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = use_ice_perspective ? get_t_out_of_q_WBF(δ_0, A_c, τ, τ_liq, q_liq, Γ_l, q_si, q_sl, true, CF_mp) : get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # don't calculate if we're starting with no liquid, calculation would fail
    (t_out_of_ice, t_out_of_ice_valid) = use_ice_perspective ? get_t_out_of_q_ice_no_WBF(δ_0, A_c, τ, τ_ice, q_ice, Γ_i) : get_t_out_of_q_ice(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si) # don't calculate if we're starting with no ice, calculation would fail
    
    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        # upgrade to BigFloat Call
        A_cL_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        τL_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) # Eq C2
        A_cI_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τI_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l)) 

        if use_ice_perspective
            A_c_big = A_cI_big
            τ_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))  # flip into ice perspective
        else
            A_c_big = A_cL_big
            τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) 
        end

        (t_out_of_liq, t_out_of_liq_vaid) = use_ice_perspective ? FT(get_t_out_of_q_WBF(big(δ_0), A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), big(q_si), big(q_sl), false).sol) : FT(get_t_out_of_q_liq(big(δ_0), A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't exit again if fail...
        # @debug "After upgrading to BigFloat, t_out_of_liq = $t_out_of_liq, t_out_of_ice = $t_out_of_ice"
    end
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        # upgrade to BigFloat Call

        if !(@isdefined A_c_big)
            A_cL_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
            τL_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) # Eq C2
            A_cI_big = (@isdefined A_cI_big) ? A_cI_big : A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
            τI_big = (@isdefined τI_big) ? τI_big : τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))

            if use_ice_perspective
                A_c_big = A_cI_big
                τ_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))  # flip into ice perspective
            else
                A_c_big = A_cL_big
                τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) 
            end
        end
        t_out_of_ice = use_ice_perspective ? FT(get_t_out_of_q_ice_no_WBF(big(δ_0), A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), false).sol) : FT(get_t_out_of_q_ice(big(δ_0), A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol)
        # @debug "After upgrading to BigFloat, t_out_of_liq = $t_out_of_liq, t_out_of_ice = $t_out_of_ice"

    end 
    
    
    
    
    max_t = min(max_t, t_out_of_liq, t_out_of_ice)
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δip}() : Val{:δ}(), param_set, target(FT(0), use_ice_perspective ? (-) : (+) ; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_δo_sat, t_hit_δo_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δipl}() : Val{:δ}() , param_set, target(FT(0), use_ice_perspective ? (+) : (-); factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δo_sat, t_hit_δo_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    
    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_δo_sat = t_hit_δo_sat_converged ? t_hit_δo_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 

    min_t, i_min_t = find_min_t(t_out_of_liq, t_out_of_ice, t_hit_δ_sat, t_hit_δo_sat, t_hit_T_freeze)
    # @debug "t_out_of_liq = $t_out_of_liq, t_out_of_ice = $t_out_of_ice, t_hit_δ_sat = $t_hit_δ_sat, t_hit_δo_sat = $t_hit_δo_sat, t_hit_T_freeze = $t_hit_T_freeze"

    if (min_t < Δt)
        if i_min_t == 1 # out of liquid
            S_ql = -q_liq / min_t
            S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, S_qi*min_t, min_t, dqvdt*min_t,ts, w; error_on_q_neg = false) # use multiplied form for floating point accuracy
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liquid before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep 
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        elseif i_min_t == 2  # out of ice
            # S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
            S_qi = -q_ice / min_t

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, S_ql*min_t, -q_ice, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # use multiplied form for floating point accuracy
            Δt_left = Δt - min_t

            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of ice before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        else # i_min_t == 3,4,5 
            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
            S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            (i_min_t==3) && # @debug "Hit $self saturation before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==4) && # @debug "Hit $other saturation before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==5) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            if i_min_t == 3 # hit saturation
                if use_ice_perspective  # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                else # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                end
            elseif i_min_t == 4 # hit other saturation
                if use_ice_perspective  # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                else # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                end
            elseif i_min_t == 5 # hit freezing
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are equal at freezing
            end

            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end

    else

        # @debug "nothing of note through end of timestep..."
        S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, Δt, Γ_l, q_si, q_sl) : S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l) # however much happens in that time, rescaled to the timestep
        S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, Δt, Γ_i) : S_qi_func( A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si)

        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
"""
    WBF has liquid (above freezing)
    - liq can shrink/grow, no ice growth, but ice can shrink
"""
function morrison_milbrandt_2015_style(
    regime::Union{WBF{true, false, false}}, # above freezing, so liq sat is lower. so we can lose liq. can't gain ice above freezing and we're subsat anyway
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    # @debug "Calling WBF{true, false, false}"
    # although we'd be moving towards liquid sat, with outside forcing could still be either, so still need to choose perspectives
    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT

    τ_ice_here = FT(Inf) # no ice growth above freezing, no shrink bc WBF is supersaturated over ice. Do this instead of just using A_c_func_no_WBF bc it allows for perspective changing.
    τ = τ_liq # no ice growth, only liq shrink

    # No ice growth, only liq shrink
    # A_c = A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) # supersaturated so no loss and now growth because above freezing

    # search locator 2
    # ======================================================================================== #
    local use_ice_perspective::Bool    

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice_here, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)

    
    A_cL = A_c_func(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τL = τ_func(τ_liq, τ_ice_here, L_i, c_p, dqsl_dT, Γ_i) # τ_ice_here essentially makes this A_c_func_no_WBF

    A_cI = A_c_func(τ_liq, Γ_i, q_si, q_sl, g, w, c_p, e_si, L_l, dqsi_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τI = τ_func(τ_ice_here, τ_liq, L_l, c_p, dqsi_dT, Γ_l)

    if (A_cI ≤ FT(0)) && (A_cL ≥ FT(0)) # swap criteeria from below freezing
        # @debug "You're probably staying in the WBF regime... unclear which perspective is better. Maybe if A_cL * τL > δ_0  (i.e. moving towards liq sat) use liq perspective? and otherwise use ice?"
        if A_cI > δ_0
            use_ice_perspective = true # we're moving towards ice sat, use ice perspective
        else
            use_ice_perspective = false # we're moving towards liq sat, use liq perspective
        end
    elseif (A_cI ≥ FT(0)) && (A_cL ≤ FT(0))
        error("Not sure how compatible asymtotes of supersat over liquid and subsat over ice are...")
    elseif (A_cI ≤ FT(0)) && (A_cL ≤ FT(0))
        # @debug "You 100% are ending up subsaturated... since above freezing maybe use liq perspective?"
        use_ice_perspective = false
    elseif (A_cI ≥ FT(0)) && (A_cL ≥ FT(0))
        # @debug "You 100 % are ending up supersaturated... since below freezing maybe use ice perspective"
        use_ice_perspective = true
    end

    #= 
        Here we shouldn't set use_WBF to false bc we're not moving down to a no WBF (single changing variable) setup, we're just swapping ice and liquid perspective completely 
        In ice perspective we should, however, use δi and δipl inplace of δ and δi. we set δ and δo accordingly
    =#
    if use_ice_perspective
        # @debug "Using ice perspective"
        A_c = A_cI
        τ = τ_func(τ_ice, τ_liq, L_l, c_p, dqsi_dT, Γ_l)  # flip into ice perspective
        δ_0o = δ_0  # other 
        δ_0 = δ_0i
        self = :ice
        other = :liq # just for printing, can deprecate
    else
        # @debug "Using liq perspective"
        A_c = A_cL
        τ = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i) 
        δ_0o = δ_0i # other
        self = :liq
        other = :ice # just for printing, can deprecate
    end 


    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = use_ice_perspective ? get_t_out_of_q_WBF(δ_0, A_c, τ, τ_liq, q_liq, Γ_l, q_si, q_sl, true, CF_mp) : get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # don't calculate if we're starting with no liquid, calculation would fail

    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        
        # upgrade to BigFloat Call
        A_cL_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        τL_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) # Eq C2
        A_cI_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τI_big = τ_func(big(τ_ice_here), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l)) 

        if use_ice_perspective
            A_c_big = A_cI_big
            τ_big = τ_func(big(τ_ice_here), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))  # flip into ice perspective
            δ_0o_big = big(δ_0)
            δ_0_big = big(δ_0i)
        else
            A_c_big = A_cL_big
            τ_big = τ_func(big(τ_liq), big(τ_ice_here), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
            δ_0_big = big(δ_0)
            δ_0o_big = big(δ_0i)
        end

        t_out_of_liq = use_ice_perspective ? FT(get_t_out_of_q_WBF(δ_0_big, A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), big(q_si), big(q_sl), false, big(CF_mp)).sol) : FT(get_t_out_of_q_liq(δ_0_big, A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't calculate if we're starting with no liquid, calculation would fail
        # @debug "After upgrading to BigFloat, t_out_of_liq = $t_out_of_liq"
    end


    max_t = min(max_t, t_out_of_liq)
    # flip target signs from above freezing, vars remain the same
    # below freezing, we're moving towards either liq or ice sat
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δip}() : Val{:δ}(), param_set, target(FT(0), use_ice_perspective ? (+) : (-) ; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_δo_sat, t_hit_δo_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δipl}() : Val{:δi}(), param_set, target(FT(0), use_ice_perspective ? (-) : (+); factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δo_sat, t_hit_δo_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, -; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)

    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_δo_sat = t_hit_δo_sat_converged ? t_hit_δo_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 

    min_t, i_min_t = find_min_t(t_out_of_liq, t_hit_δ_sat, t_hit_δo_sat, t_hit_T_freeze)
    # @debug "t_out_of_liq = $t_out_of_liq, t_hit_δ_sat = $t_hit_δ_sat, t_hit_δo_sat = $t_hit_δo_sat, t_hit_T_freeze = $t_hit_T_freeze"

    # ======================================================================================== #
    # ======================================================================================== #

    if (min_t < Δt)
        if i_min_t == 1 # out of liq, we're done
            S_ql = -q_liq / min_t
            S_qi = FT(0)

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, S_qi*min_t, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # use multiplied form for floating point accuracy
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liquid before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        else # i_min_t == 2,3,4

            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
            S_qi = FT(0)

            if (S_qi > FT(0))
                # @debug "S_qi > 0 in WBF above freezing... can't happen. Remedying by recalling with τ_ice = Inf"
                return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
            end


            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)

            (i_min_t==2) && # @debug "Hit $self saturation before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==3) && # @debug "Hit $other saturation before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==4) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."

            if i_min_t == 2 # hit saturation
                if use_ice_perspective  # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                else # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                end
            elseif i_min_t == 3 # hit other saturation
                if use_ice_perspective  # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                else # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                end
            else # i_min_t == 4
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
            end

            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end

    else 
        # @debug "nothing of note through end of timestep..."
        S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
        S_qi = FT(0)

        if (S_qi > FT(0))
            # @debug "S_qi > 0 in WBF above freezing... can't happen. Remedying by recalling with τ_ice = Inf"
            return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
        end

        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


"""
    WBF has liquid (above freezing)
    - liq can shrink/grow, no ice growth, but ice can shrink
"""
function morrison_milbrandt_2015_style(
    regime::WBF{true, true, false}, # above freezing so we're subsat over ice, supersat over liq. ice can shrink, never grow bc above freeezing.
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    # @debug "Calling WBF{true, true, false}"

    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts) 

    # No ice growth, only liq shrink
    # A_c = A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) # supersaturated so no loss and now growth because above freezing

    # WBF above freezing is subsat for ice, supersat for liq. So here ice can shrink. But hitting ice sat does nothing because ice can't grow... should we still use ice perspective? may stop oscillations... hmmm
    # Need to reject any positive ice tendencies though... can do the method of recalling w/ timescale inf (preferable to just zeroing out ice tendencies bc then we can recalculate tendencies on other things)

    # search_locator
    # ======================================================================================== #
    local use_ice_perspective::Bool    

    
    A_cL = A_c_func(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τL = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i) # Eq C2

    A_cI = A_c_func(τ_liq, Γ_i, q_si, q_sl, g, w, c_p, e_si, L_l, dqsi_dT, dqvdt, dTdt, p, ρ) # Eq C4
    τI = τ_func(τ_ice, τ_liq, L_l, c_p, dqsi_dT, Γ_l)

    if (A_cI ≤   FT(0)) && (A_cL ≥ FT(0)) # swap criteeria from below freezing
        # @debug "You're probably staying in the WBF regime... unclear which perspective is better. Maybe if A_cL * τL > δ_0  (i.e. moving towards liq sat) use liq perspective? and otherwise use ice?"
        if A_cI > δ_0
            use_ice_perspective = true # we're moving towards ice sat, use ice perspective
        else
            use_ice_perspective = false # we're moving towards liq sat, use liq perspective
        end
    elseif (A_cI ≥ FT(0)) && (A_cL ≤ FT(0))
        error("Not sure how compatible asymtotes of supersat over liquid and subsat over ice are...")
    elseif (A_cI ≤ FT(0)) && (A_cL ≤ FT(0))
        # @debug "You 100% are ending up subsaturated... since above freezing maybe use liq perspective?"
        use_ice_perspective = false
    elseif (A_cI ≥ FT(0)) && (A_cL ≥ FT(0))
        # @debug "You 100 % are ending up supersaturated... since below freezing maybe use ice perspective"
        use_ice_perspective = true
    end

    #= 
        Here we shouldn't set use_WBF to false bc we're not moving down to a no WBF (single changing variable) setup, we're just swapping ice and liquid perspective completely 
        In ice perspective we should, however, use δi and δipl inplace of δ and δi. we set δ and δo accordingly
    =#
    if use_ice_perspective
        # @debug "Using ice perspective"
        A_c = A_cI
        τ = τ_func(τ_ice, τ_liq, L_l, c_p, dqsi_dT, Γ_l)  # flip into ice perspective
        δ_0o = δ_0  # other 
        δ_0 = δ_0i
        self = :ice
        other = :liq # just for printing, can deprecate
    else
        # @debug "Using liq perspective"
        A_c = A_cL
        τ = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i) 
        δ_0o = δ_0i # other
        self = :liq
        other = :ice # just for printing, can deprecate
    end 


    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = use_ice_perspective ? get_t_out_of_q_WBF(δ_0, A_c, τ, τ_liq, q_liq, Γ_l, q_si, q_sl, true, CF_mp) : get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # don't calculate if we're starting with no liquid, calculation would fail
    (t_out_of_ice, t_out_of_ice_valid) = use_ice_perspective ? get_t_out_of_q_ice_no_WBF(δ_0, A_c, τ, τ_ice, q_ice, Γ_i) : get_t_out_of_q_ice(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si) # don't calculate if we're starting with no ice, calculation would fail

    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        # upgrade to BigFloat Call
        A_cL_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        τL_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) # Eq C2
        A_cI_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
        τI_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l)) # Eq C4

        if use_ice_perspective
            A_c_big = A_cI_big
            τ_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))  # flip into ice perspective
            δ_0o_big = big(δ_0)
            δ_0_big = big(δ_0i)
        else
            A_c_big = A_cL_big
            τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
            δ_0o_big = big(δ_0i)
        end

        t_out_of_liq = use_ice_perspective ? FT(get_t_out_of_q_WBF(δ_0_big , A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), big(q_si), big(q_sl), false, big(CF_mp)).sol) : FT(get_t_out_of_q_liq(δ_0_big , A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol)
        # @debug "After upgrading to BigFloat, t_out_of_liq = $t_out_of_liq"
    end
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        if !(@isdefined A_c_big)
            A_cL_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
            τL_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) # Eq C2
            A_cI_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
            τI_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l)) # Eq C4

            if use_ice_perspective
                A_c_big = A_cI_big
                τ_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))  # flip into ice perspective
                δ_0o_big = big(δ_0)
                δ_0_big = big(δ_0i)
            else
                A_c_big = A_cL_big
                τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
                δ_0o_big = big(δ_0i)
            end
        end
        t_out_of_ice = use_ice_perspective ? FT(get_t_out_of_q_ice_no_WBF(δ_0_big, A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), false).sol) : FT(get_t_out_of_q_ice(δ_0_big, A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol)
        # @debug "After upgrading to BigFloat, t_out_of_ice = $t_out_of_ice"
    end






    max_t = min(max_t, t_out_of_liq, t_out_of_ice)
    # flip target signs from above freezing, vars remain the same
    # below freezing, we're moving towards either liq or ice sat
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δip}() : Val{:δ}(), param_set, target(FT(0), use_ice_perspective ? (+) : (-) ; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_δo_sat, t_hit_δo_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δipl}() : Val{:δi}(), param_set, target(FT(0), use_ice_perspective ? (-) : (+); factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δo_sat, t_hit_δo_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, -; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)

    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_δo_sat = t_hit_δo_sat_converged ? t_hit_δo_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 

    min_t, i_min_t = find_min_t(t_out_of_liq, t_out_of_ice, t_hit_δ_sat, t_hit_δo_sat, t_hit_T_freeze)
    # @debug "t_out_of_liq = $t_out_of_liq, t_out_of_ice = $t_out_of_ice, t_hit_δ_sat = $t_hit_δ_sat, t_hit_δo_sat = $t_hit_δo_sat, t_hit_T_freeze = $t_hit_T_freeze"

    # ======================================================================================== #

    if (min_t < Δt)
        if i_min_t == 1 # out of liq, we're done
            S_ql = -q_liq / min_t
            S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            if (S_qi > FT(0))
                # @debug "S_qi > 0 in WBF above freezing... can't happen. Remedying by recalling with τ_ice = Inf"
                return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
            end

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, S_qi*min_t, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # use multiplied form for floating point accuracy
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liquid before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
        elseif i_min_t == 2 # out of ice
            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
            S_qi = -q_ice / min_t

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, S_ql*min_t, -q_ice, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # use multiplied form for floating point accuracy
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liquid before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))

        else # i_min_t == 2,3,4
            # S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            # S_qi = FT(0)

            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
            S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            if (S_qi > FT(0))
                # @debug "S_qi > 0 in WBF above freezing... can't happen. Remedying by recalling with τ_ice = Inf"
                return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
            end


            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)

            (i_min_t==3) && # @debug "Hit $self before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==4) && # @debug "Hit $other saturation before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==5) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            if i_min_t == 3 # hit self saturation
                if use_ice_perspective  # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                else # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                end
            elseif i_min_t == 4 # hit other saturation
                if use_ice_perspective  # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                else # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                end
            else # i_min_t == 5
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
            end

            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end

    else 
        # @debug "nothing of note through end of timestep..."
        # S_ql = S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        # S_qi = FT(0)
        S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l) # however much happens in that time, rescaled to the timestep
        S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

        if (S_qi > FT(0))
            # @debug "S_qi > 0 in WBF above freezing... can't happen. Remedying by recalling with τ_ice = Inf"
            return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, fallback_to_standard_supersaturation_limiter=fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
        end

        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    WBF no liquid (below freezing)
    - no liq, just ice growth
"""
function morrison_milbrandt_2015_style(
    regime::Union{WBF{false, true, true}, WBF{false, false, true}}, # no liq, no matter what we'll end up with some ice
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    # @debug "Calling WBF{false, $(q.ice > FT(0)), true}"

    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT
    local use_ice_perspective::Bool

    τ_liq_here = FT(Inf) # we have no liq and can't immediately grow it while in WBF mode (don't actually set in case need to pass to another func)
    τ = τ_ice # no liq growth, only ice growth

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq_here, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)

    # iszero(q.liq) ||  error("we have liq $(q.liq), this is wrong")


    
    A_cL = A_c_func_no_WBF( q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_cI = A_c_func_no_WBF( q_si, g, w, c_p, e_si, dqsi_dT, dqvdt, dTdt, p, ρ)

    if (A_cL ≤ FT(0)) && (A_cI ≥ FT(0))
        # @debug "You're probably staying in the WBF regime... unclear which perspective is better. Maybe if A_cL * τL > δ_0  (i.e. moving towards liq sat) use liq perspective? and otherwise use ice?"
        if (A_cL) > δ_0
            use_ice_perspective = false # we're moving towards liq sat, use liq perspective
        else
            use_ice_perspective = true # we're moving towards ice sat, use ice perspective
        end
    elseif (A_cL ≥ FT(0)) && (A_cI ≤ FT(0))
        # # @debug "Not sure how compatible asymtotes of supersat over liquid and subsat over ice are..."
        error("Not sure how compatible asymtotes of supersat over liquid and subsat over ice are...")
    elseif (A_cL ≤ FT(0)) && (A_cI ≤ FT(0))
        # @debug "You 100% are ending up subsaturated... since below freezing maybe use ice perspective?"
        use_ice_perspective = true
    elseif (A_cL ≥ FT(0)) && (A_cI ≥ FT(0))
        # @debug "You 100 % are ending up supersaturated... since above freezing maybe use liq perspective"
        use_ice_perspective = false
    end
    use_ice_perspective = true
    #= 
        Here we shouldn't set use_WBF to false bc we're not moving down to a no WBF (single changing variable) setup, we're just swapping ice and liquid perspective completely 
        We should, however, use δ instead of δi and δipl instead of δ
    =#
    if use_ice_perspective
        # @debug "Using ice perspective"
        A_c = A_cI
        τ = τ_func(τ_ice, τ_liq_here, L_l, c_p, dqsi_dT, Γ_l) 
        δ_0o = δ_0  # other 
        δ_0 = δ_0i
        self = :ice
        other = :liq # just for printing, can deprecate
    else
        # @debug "Using liq perspective"
        A_c = A_cL
        τ = τ_func(τ_liq_here, τ_ice, L_l, c_p, dqsl_dT, Γ_i) # treat as if is only ice cause no liq <-- can't use pure liq perspective bc no wbf, but Can have A_c defined from liq perspective...
        δ_0o = δ_0i # other
        # δ_0 = δ_0 # already set
        self = :liq
        other = :ice # just for printing, can deprecate
    end

    max_t = Δt
    (t_out_of_ice, t_out_of_ice_valid) = use_ice_perspective ? get_t_out_of_q_ice_no_WBF(δ_0, A_c, τ, τ_ice, q_ice, Γ_i) : get_t_out_of_q_ice(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si) # no liq growth while subsat, only ice changes
    
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        # upgrade to BigFloat Call
        if use_ice_perspective
            A_cI_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
            A_c_big = A_cI_big
            τ_big = τ_func(big(τ_ice), big(τ_liq_here), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l)) # Eq C4
            δ_0o_big = big(δ_0)  # other
            δ_0_big = big(δ_0i)
        else
            A_cL_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C2
            A_c_big = A_cL_big
            τ_big = τ_func(big(τ_liq_here), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i)) # Eq C2
            δ_0o_big = big(δ_0i) # other
            δ_0_big = big(δ_0) # already set
        end

        t_out_of_ice = use_ice_perspective ? FT(get_t_out_of_q_ice_no_WBF(δ_0_big, A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), false).sol) : FT(get_t_out_of_q_ice(δ_0_big, A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol)
    end
    
    
    
    
    max_t = min(max_t, t_out_of_ice)
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δip}() : Val{:δ}(), param_set, target(FT(0), use_ice_perspective ? (-) : (+) ; factor = FT(TF)), FT(0), max_t, δ_0, A_c, τ, τ_liq_here, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_δo_sat, t_hit_δo_sat_converged  = get_t_var_hit_value(use_ice_perspective ? Val{:δipl}() : Val{:δi}() , param_set, target(FT(0), use_ice_perspective ? (+) : (-); factor = FT(TF)), FT(0), max_t, δ_0, A_c, τ, τ_liq_here, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_δo_sat, t_hit_δo_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq_here, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    
    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_δo_sat = t_hit_δo_sat_converged ? t_hit_δo_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 

    min_t, i_min_t = find_min_t(t_out_of_ice, t_hit_δ_sat, t_hit_δo_sat, t_hit_T_freeze)

    if (min_t < Δt)
        if i_min_t == 1 # out of ice
            S_ql = FT(0) # no liq gain during WBF and started with none
            S_qi = -q_ice / min_t

            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, FT(0), -q_ice, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # use multiplied form for floating point accuracy
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of ice before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))

        else # i_min_t == 2,3,4 
            S_ql = FT(0) # no liq gain during WBF and started with none
            S_qi = use_ice_perspective ? S_func_indiv_no_WBF( A_c, τ, δ_0, min_t, Γ_i) : S_func_indiv_no_WBF( A_c, τ, δ_0i, min_t, Γ_i) 

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            (i_min_t==2) && # @debug "Hit $self before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==3) && # @debug "Hit $other saturation before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            (i_min_t==4) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(regime)) at t = $(min_t)..."
            if i_min_t == 2 # hit liq saturation
                if use_ice_perspective  # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                else # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                end
            elseif i_min_t == 3 # hit ice saturation
                if use_ice_perspective  # hit liq saturation
                    new_δ_0_shum = FT(0)
                    new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
                else # hit ice saturation
                    new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                    new_δ_0i_shum = FT(0)
                end
            else # i_min_t == 4
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
            end

            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end

    else

        # @debug "nothing of note through end of timestep..."
        # just regular return
        S_ql = FT(0)
        S_qi = use_ice_perspective ? S_func_indiv_no_WBF( A_c, τ, δ_0, Δt, Γ_i) : S_func_indiv_no_WBF( A_c, τ, δ_0i, Δt, Γ_i) 
        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
    
"""
    WBF no liquid (above freezing)
    - no liq, no ice growth
"""

function morrison_milbrandt_2015_style(
    regime::Union{WBF{false, true, false}, WBF{false, false, false}}, # above freezing, no liq, ice can't shrink.
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)

    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT 

    HI::Bool = (q.ice > FT(0))
    # @debug "Calling WBF{false, $HI, false}"


    τ_liq_here = FT(Inf) # we have no liq and can't immediately grow it while in WBF mode (don't actually set in case need to pass to another func)
    τ_ice_here = FT(Inf) # no ice growth above freezing
    τ_here = FT(Inf)

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq_here, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)


    A_c = A_c_func_no_WBF( q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)

    # no ice changing above freezing so no growth. No shrinking because is WBF so supersaturated
    # no liq growth because subsat. No liq shrink because we have no liq.

    min_t = Δt
    t_hit_δ_sat, t_hit_δ_sat_converged = get_t_var_hit_value(Val{:δ}(), param_set, target(FT(0), -; factor = FT(TF)), FT(0), min_t,  δ_0, A_c, τ_here, τ_liq_here, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts; use_WBF = false) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    min_t = min_if_converged(min_t, t_hit_δ_sat, t_hit_δ_sat_converged)
    t_hit_δi_sat, t_hit_δi_sat_converged = get_t_var_hit_value(Val{:δi}(), param_set, target(FT(0), +; factor = FT(TF)), FT(0), min_t,  δ_0, A_c, τ_here, τ_liq_here, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts; use_WBF = false) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    min_t = min_if_converged(min_t, t_hit_δi_sat, t_hit_δi_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, -; factor = FT(TF)), FT(0), min_t,  δ_0, A_c, τ_here, τ_liq_here, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts; use_WBF = false)
    # @debug "t_hit_δ_sat = $t_hit_δ_sat, t_hit_δi_sat = $t_hit_δi_sat, t_hit_T_freeze = $t_hit_T_freeze | BF = false"
    
    t_hit_δ_sat = t_hit_δ_sat_converged ? t_hit_δ_sat : FT(Inf) 
    t_hit_δi_sat = t_hit_δi_sat_converged ? t_hit_δi_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 
    min_t, i_min_t = find_min_t(t_hit_δ_sat, t_hit_δi_sat, t_hit_T_freeze)
    
    
    if (min_t < Δt)
        # all are the same
        S_ql = FT(0) # no liq gain during WBF and started with none
        S_qi = FT(0) # no ice growth above freezing

        new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
        Δt_left = Δt - min_t
        new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
        (i_min_t == 1) && # @debug "Hit liquid saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
        (i_min_t == 2) && # @debug "Hit ice saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
        (i_min_t == 3) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."

        if i_min_t == 1 # hit liq saturation
            new_δ_0_shum = FT(0)
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
        elseif i_min_t == 2 # hit ice saturation
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = FT(0)
        else # i_min_t == 3
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
        end

        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

        S_ql *= min_t / Δt # rescale to the timestep
        S_qi *= min_t / Δt # rescale to the timestep
        S_ql_addit *= Δt_left / Δt # rescale to the remaining time
        S_qi_addit *= Δt_left / Δt # rescale to the remaining time
        return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))              
            
    else
        # @debug "nothing of note through end of timestep..."
        # just regular return
        S_ql = FT(0)
        S_qi = FT(0)
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))

    end
end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated with liq and ice
    - both shrink
    - if below freezing, would hit ice_sat first and stop
    - if above freezing, would hit liq_sat first and stop
"""
function morrison_milbrandt_2015_style(
    regime::Union{Subsaturated{true, true, true}, Subsaturated{true, true, false}}, # can run out of either first
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter, liq_fraction, ice_fraction, cld_fraction) = opts
    CF_mp = min(liq_fraction, ice_fraction)


    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts) 

    local BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{true, true, $BF}"

    # (δ_0 < FT(0)) && (δ_0i < FT(0)) || error("δ_0 and δ_0i should be negative, got δ_0 = $δ_0, δ_0i = $δ_0i")
    # !iszero(q_liq) || error("liq should not be zero")
    # !iszero(q_ice) || error("ice should not be zero")

    local use_ice_perspective::Bool = BF # if we're below freezing, we're heading to ice sat. use ice perspective. if above freezing, we're heading to liq sat. use liq perspective

    #= 
        Here we shouldn't set use_WBF to false bc we're not moving down to a no WBF (single changing variable) setup, we're just swapping ice and liquid perspective completely 
        We should, however, use δ instead of δi and δipl instead of δ
    =#
    if use_ice_perspective
        # @debug "Using ice perspective"
        A_c = A_c_func(τ_liq, Γ_i, q_si, q_sl, g, w, c_p, e_si, L_l, dqsi_dT, dqvdt, dTdt, p, ρ) # Eq C4
        τ = τ_func(τ_ice, τ_liq, L_l, c_p, dqsi_dT, Γ_l)
        δ_0o = δ_0  # other 
        δ_0 = δ_0i
        other = :liq # just for printing, can deprecate
    else
        # @debug "Using liq perspective"
        A_c = A_c_func(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
        τ = τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)
        δ_0o = δ_0i # other
        other = :ice # just for printing, can deprecate
    end

    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = use_ice_perspective ? get_t_out_of_q_WBF(δ_0, A_c, τ, τ_liq, q_liq, Γ_l, q_si, q_sl, true, CF_mp) : get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) # don't calculate if we're starting with no liquid, calculation would fail
    (t_out_of_ice, t_out_of_ice_valid) = use_ice_perspective ? get_t_out_of_q_ice_no_WBF(δ_0, A_c, τ, τ_ice, q_ice, Γ_i) : get_t_out_of_q_ice(δ_0, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si) # don't calculate if we're starting with no ice, calculation would fail   
    

    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end


        if use_ice_perspective
            A_c_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
            τ_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))
            δ_0o_big = big(δ_0)  # other
            δ_0_big = big(δ_0i)
        else
            A_c_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
            τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
            δ_0o_big = big(δ_0i) # other
            δ_0_big = big(δ_0) #
        end

        t_out_of_liq = use_ice_perspective ? FT(get_t_out_of_q_WBF(δ_0_big, A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), big(q_si), big(q_sl), false, big(CF_mp)).sol) : FT(get_t_out_of_q_liq(δ_0_big, A_c_big, τ_big, big(τ_liq), big(q_liq), big(Γ_l), false).sol) # don't calculate if we're starting with no liquid, calculation would fail
        # @debug "After upgrading to big, t_out_of_liq = $t_out_of_liq"
    end
    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        if !(@isdefined A_c_big)
            if use_ice_perspective
                A_c_big = A_c_func(big(τ_liq), big(Γ_i), big(q_si), big(q_sl), big(g), big(w), big(c_p), big(e_si), big(L_l), big(dqsi_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
                τ_big = τ_func(big(τ_ice), big(τ_liq), big(L_l), big(c_p), big(dqsi_dT), big(Γ_l))
                δ_0o_big = big(δ_0)  # other
                δ_0_big = big(δ_0i)
            else
                A_c_big = A_c_func(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
                τ_big = τ_func(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
                δ_0o_big = big(δ_0i) # other
                δ_0_big = big(δ_0) #
            end
        end
        t_out_of_ice = use_ice_perspective ? FT(get_t_out_of_q_ice_no_WBF(δ_0_big, A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), false).sol) : FT(get_t_out_of_q_ice(δ_0_big, A_c_big, τ_big, big(τ_ice), big(q_ice), big(Γ_i), big(q_sl), big(q_si), false).sol) # don't calculate if we're starting with no ice, calculation would fail
        # @debug "After upgrading to big, t_out_of_ice = $t_out_of_ice"
    end
    



    max_t = min(t_out_of_liq, t_out_of_ice, max_t)
    if BF # below freezing, heading for ice sat first, from ice perspective (δip) or from liq perspective (δi)
        t_hit_sat, t_hit_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δip}() : Val{:δi}(), param_set, target(FT(0), + ; factor = FT(TF)), FT(0), max_t, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    else # above freezing we're heading for liq sat first, using ice perspective (δipl) or liq perspective (δ)
        t_hit_sat, t_hit_sat_converged = get_t_var_hit_value(use_ice_perspective ? Val{:δipl}() : Val{:δ}(), param_set, target(FT(0), + ; factor = FT(TF)), FT(0), max_t, δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the
    end
    max_t = min_if_converged(max_t, t_hit_sat, t_hit_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged  = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, BF ? (+) : (-); factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    
    t_hit_sat = t_hit_sat_converged ? t_hit_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 
    min_t, i_min_t = find_min_t(t_out_of_liq, t_out_of_ice, t_hit_sat, t_hit_T_freeze)
    
    if (min_t < Δt)
        if i_min_t == 1 # out of liq
            S_ql = -q_liq / min_t
            S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func(A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            if !BF && (S_qi > FT(0))
                # @debug "S_qi > 0 in Subsaturated{false, true, false}... setting to zero by making ice time scale infinite"
                return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth)) # we had ice growth instead of loss, set τ_ice to infinity and recall
            end


            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, S_qi*min_t, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liq before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        elseif i_min_t == 2 # out of ice
            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = -q_ice / min_t

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, S_ql*min_t, -q_ice, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of ice before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        else # i_min_t == 3,4 
            S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, min_t, Γ_l, q_si, q_sl, CF_mp) : S_ql_func( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, min_t, Γ_i) : S_qi_func(A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            if !BF && (S_qi > FT(0))
                # @debug "S_qi > 0 in Subsaturated{false, true, false}... setting to zero by making ice time scale infinite"
                return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, δ_0_shum = δ_0, δ_0i_shum = δ_0i, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
            end

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, new_T, T_freeze)
            (i_min_t==3) && # @debug "Hit ice saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            (i_min_t==4) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."

            if i_min_t == 3 # hit sat
                new_δ_0_shum = BF ? TD.vapor_specific_humidity(new_q) - new_q_eq.liq : FT(0)
                new_δ_0i_shum = BF ? FT(0) : TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            else # i_min_t == 4 # hit freezing
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
            end

            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end
    else

        # @debug "nothing of note through end of timestep..."
        S_ql = use_ice_perspective ? S_func_WBF(A_c, τ, τ_liq, δ_0, Δt, Γ_l, q_si, q_sl) : S_ql_func( A_c, τ, τ_liq, δ_0, Δt, Γ_l)
        S_qi = use_ice_perspective ? S_func_no_WBF(A_c, τ, τ_ice, δ_0, Δt, Γ_i) : S_qi_func(A_c, τ, τ_ice, δ_0, Δt, Γ_i, q_sl, q_si)

        if !BF && (S_qi > FT(0))
            # @debug "S_qi > 0 in Subsaturated{false, true, false}... setting to zero by making ice time scale infinite"
            return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, fallback_to_standard_supersaturation_limiter=fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
        end

        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
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
function morrison_milbrandt_2015_style(
    regime::Union{Subsaturated{true, false, true}, Subsaturated{true, false, false}}, # can run out of liq first. if not and we make it to ice sat, transition to WBF
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter) = opts


    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end



    # no ice and can't gain any bc subsat
    τ_ice_here = FT(Inf) # we have no ice, and ice can't grow bc subsat
    τ = τ_liq # liq can decay
    
    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq, τ_ice_here, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)

    
    local BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{true, false, $BF}"

    # iszero(q_ice) || error("ice should be zero, got q_ice = $(q_ice)")

    A_c = A_c_func_no_WBF( q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)

    max_t = Δt
    (t_out_of_liq, t_out_of_liq_valid) = get_t_out_of_q_liq(δ_0, A_c, τ, τ_liq, q_liq, Γ_l) 

    if !t_out_of_liq_valid
        if fallback_to_standard_supersaturation_limiter
            q_vap = TD.vapor_specific_humidity(q)
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        A_c_big = A_c_func_no_WBF(big(q_sl), big(g), big(w), big(c_p), big(e_sl), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        t_out_of_liq = FT(get_t_out_of_q_liq(big(δ_0), A_c_big, big(τ), big(τ_liq), big(q_liq), big(Γ_l), false).sol)
        # @debug "After upgrading to big, t_out_of_liq = $t_out_of_liq"
    end


    max_t = min(t_out_of_liq, max_t)
    t_hit_sat, t_hit_sat_converged = get_t_var_hit_value(BF ? Val{:δi}() : Val{:δ}(), param_set, target(FT(0), +; factor = FT(TF)), FT(0), max_t, δ_0, A_c, τ, τ_liq, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_sat, t_hit_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, (BF ? (+) : (-)); factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    # @debug "t_out_of_liq = $t_out_of_liq, t_hit_sat = $t_hit_sat, t_hit_T_freeze = $t_hit_T_freeze"
    
    t_hit_sat = t_hit_sat_converged ? t_hit_sat : FT(Inf)
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 

    min_t, i_min_t = find_min_t(t_out_of_liq, t_hit_sat, t_hit_T_freeze)

    if (min_t < Δt)
        if i_min_t == 1 # out of liq
            S_ql = -q_liq / min_t
            S_qi = FT(0) # no ice growth above freezing

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, -q_liq, S_qi*min_t, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # floating point safe
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, T, T_freeze)
            # @debug "Ran out of liq before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt - min_t, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        else # i_min_t == 2,3 
            S_ql = S_ql_func_indiv( A_c, τ_liq, δ_0, min_t, Γ_l)
            S_qi = FT(0)

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, new_T, T_freeze)
            sat_hit = (T < T_freeze) ? "ice" : "liq"
            (i_min_t==2) && # @debug "Hit $sat_hit saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            (i_min_t==3) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            if i_min_t == 2 # hit sat
                new_δ_0_shum = BF ? TD.vapor_specific_humidity(new_q) - new_q_eq.liq : FT(0)
                new_δ_0i_shum = BF ? FT(0) : TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            else # i_min_t == 3 # hit freezing
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
            end
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))

            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        end

    else

        # @debug "nothing of note through end of timestep..."
        S_ql = S_ql_func_indiv( A_c, τ_liq, δ_0, Δt, Γ_l)
        S_qi = FT(0)
        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

"""
    Subsaturated ice but no liq 
    - no ice so liq shrink
    - if below freezing, would hit ice_sat first and stop
    - if above freezing, would hit liq_sat first and stop
"""
function morrison_milbrandt_2015_style(
    regime::Union{Subsaturated{false, true, true}, Subsaturated{false, true, false}}, # can run out of ice first. then we're stuck.
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter) = opts

    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    local S_ql_addit::FT
    local S_qi_addit::FT

    τ_liq_here = FT(Inf) # no liq growth and no liq to evap
    τ = τ_ice # we have no liq


    # (q.liq < eps(FT)) || error("shouldn't have any liq here (up to floating point error)..., got q_liq = $(q.liq)")
    # (q.ice > eps(FT)) || error("should have some ice here (up to floating point error)..., got q_ice = $(q.ice)")


    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq_here, τ_ice, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)


    (q_vap < q_si)  || error("should be subsaturated, got q_vap = $q_vap, q_si = $q_si, (q_vap - q_si) = $(q_vap - q_si)")

    A_c = A_c_func_no_WBF( q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)

    # need to treet like is liquid...
    local BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{false, true, $BF}"

    max_t = Δt
    (t_out_of_ice, t_out_of_ice_valid) = get_t_out_of_q_ice_no_WBF(δ_0i, A_c, τ, τ_ice, q_ice, Γ_i) # calcualte as if is liquid since  there's no WBF

    if !t_out_of_ice_valid
        if fallback_to_standard_supersaturation_limiter
            # @debug "Falling back to StandardSupersaturationMoistureSourcesLimiter due to t_out_of_ice being nothing"
            return standard_supersaturation_sources(StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts) # use this version with δ_0, δ_0i to keep our up to date supersaturation
        end
        A_c_big = A_c_func_no_WBF(big(q_sl), big(g), big(w), big(c_p), big(e_sl), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ))
        t_out_of_ice = FT(get_t_out_of_q_ice_no_WBF(big(δ_0i), A_c_big, big(τ), big(τ_ice), big(q_ice), big(Γ_i), false).sol) # calcualte as if is liquid since  there's no WBF
        # @debug "After upgrading to big, t_out_of_ice = $t_out_of_ice"
    end

    max_t = min(max_t, t_out_of_ice)
    t_hit_sat, t_hit_sat_converged = get_t_var_hit_value(BF ? Val{:δi}() : Val{:δ}(), param_set, target(FT(0), +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq_here, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts; use_WBF = false) # if didn't hit sat and t_hit_sat is the least, it'll just equal Δt, then we can use these S_ql and S_qi and shouldn't have to recalculate them at the end...
    max_t = min_if_converged(max_t, t_hit_sat, t_hit_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, target(T_freeze, (BF ? (+) : (-)); factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ, τ_liq_here, τ_ice, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts; use_WBF = false)
    
    t_hit_sat = t_hit_sat_converged ? t_hit_sat : FT(Inf) 
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 
    min_t, i_min_t = find_min_t(t_out_of_ice, t_hit_sat, t_hit_T_freeze)

    if (min_t < Δt) # either way it's the same up to 
        if i_min_t == 1 # out of ice
            S_ql = FT(0) # no liq to lose, cant grow in subsat
            S_qi = -q_ice / min_t

            # S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi

            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper_full_source(param_set, p, q, q_eq, q_liq, q_ice, S_ql*min_t, -q_ice, min_t, dqvdt*min_t, ts, w; error_on_q_neg = false) # floating point safe
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, new_T, T_freeze)
            # @debug "Ran out of ice before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

        else i_min_t == 2,3 #
            S_ql = FT(0) # no liq to lose, cant grow in subsat
            S_qi = S_qi_func_indiv( A_c, τ_ice, δ_0i, Δt, Γ_i)

            if !BF && (S_qi > FT(0))
                # @debug "S_qi > 0 in Subsaturated{false, true, false}... setting to zero by making ice time scale infinite"
                return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter ) # we had ice growth instead of loss, set τ_ice to infinity and recall
            end

            S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, min_t) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
            new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
            Δt_left = Δt - min_t
            new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, new_T, T_freeze)
            (i_min_t==2) && # @debug "Hit ice saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(t_hit_sat)..."
            (i_min_t==3) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(t_hit_T_freeze)..."
            if i_min_t == 2 # hit sat
                new_δ_0_shum = BF ? TD.vapor_specific_humidity(new_q) - new_q_eq.liq : FT(0)
                new_δ_0i_shum = BF ? FT(0) : TD.vapor_specific_humidity(new_q) - new_q_eq.ice
            else # i_min_t == 3 # hit freezing
                new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
                new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
            end
            S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
            S_ql *= min_t / Δt # rescale to the timestep
            S_qi *= min_t / Δt # rescale to the timestep
            S_ql_addit *= Δt_left / Δt # rescale to the remaining time
            S_qi_addit *= Δt_left / Δt # rescale to the remaining time
            return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))
        
        end
        
    else
        # @debug "nothing of note through end of timestep..."
        S_ql = FT(0)
        S_qi = S_qi_func_indiv( A_c, τ_ice, δ_0i, Δt, Γ_i)
        if !BF && (S_qi > FT(0))
            # @debug "S_qi > 0 in Subsaturated{false, true, false}... setting to zero by making ice time scale infinite"
            return morrison_milbrandt_2015_style(regime, param_set, area, ρ, p, T, w, τ_liq, FT(Inf), q_vap, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = return_mixing_ratio, max_depth = max_depth, depth = depth, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter)) # we had ice growth instead of loss, set τ_ice to infinity and recall
        end
        S_ql, S_qi = S_above_floor(S_ql, S_qi, q_liq, q_ice, Δt) # if t was too small for t_out_of_q for example, this makes sure we put a floor on S_ql and S_qi
        return return_mixing_ratio ? (S_ql, S_qi) : (S_mixing_ratio_to_shum(S_ql, q.tot), S_mixing_ratio_to_shum(S_qi, q.tot))
    end

end


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
"""
    Subsaturated with nothing
    - return 0, 0
"""
function morrison_milbrandt_2015_style(
    regime::Union{Subsaturated{false, false, false}, Subsaturated{false, false, true}}, # kind of null, you have nothing and can't make anything
    param_set::APS, area::FT, ρ::FT, p::FT, T::FT, w::FT, τ_liq::FT, τ_ice::FT, q_vap::FT, q::TD.PhasePartition, q_eq::TD.PhasePartition, Δt::FT, ts::TD.ThermodynamicState;
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
    )::Tuple{FT, FT} where {FT}
    (; use_fix, return_mixing_ratio, max_depth, depth, δ_0_shum, δ_0i_shum, dqvdt, dTdt, fallback_to_standard_supersaturation_limiter) = opts

    if depth == max_depth
        # @debug "Max depth reached, returning exponential part only fallback (no thermo adjustments)"
        q_vap = TD.vapor_specific_humidity(q)
        S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = false, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        return return_mixing_ratio ? (S_shum_to_mixing_ratio(S_ql, q.tot), S_shum_to_mixing_ratio(S_qi, q.tot)) : (S_ql, S_qi) # exponential only doesn't bother with mixing ratios at the moment
    end

    τ_liq_here = FT(Inf) # no liq growth and no liq to evap
    τ_ice_here = FT(Inf) # no ice growth and no ice to sublimate
    τ_here = FT(Inf) # no liq or ice

    (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_vap, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i) = get_params_and_go_to_mixing_ratio(param_set, area, ρ, p, T, w, τ_liq_here, τ_ice_here, q_vap, dqvdt, dTdt, δ_0_shum, δ_0i_shum, q, q_eq, Δt, ts)


    local BF::Bool = T < T_freeze
    # @debug "Calling Subsaturated{false, false, $BF}"

    # (q_vap < (true ? q_si : q_sl)) || error("should be subsaturated, q_vap = $q_vap, q_vap_sat = $(true ? q_si : q_sl)")

    A_c = A_c_func_no_WBF( q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)

    max_t = Δt
    # no ice or liq, ice/liq growth/decay
    t_hit_sat, t_hit_sat_converged = get_t_var_hit_value(BF ? Val{:δi}() : Val{:δ}(), param_set, target(FT(0), +; factor = FT(TF)), FT(0), max_t,  δ_0, A_c, τ_here, τ_liq_here, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    max_t = min_if_converged(max_t, t_hit_sat, t_hit_sat_converged)
    t_hit_T_freeze, t_hit_T_freeze_converged = get_t_var_hit_value(Val{:T}(), param_set, T_freeze, target(T_freeze, (BF ? (+) : (-)); factor = FT(TF)), max_t,  δ_0, A_c, τ_here, τ_liq_here, τ_ice_here, Γ_l, Γ_i, q_sl, q_si, T, w, p, q_vap, dqvdt, dTdt, q, q_eq, q_liq, q_ice, ts)
    # @debug " = $t_hit_sat = $t_hit_sat, t_hit_T_freeze = $t_hit_T_freeze"
    t_hit_sat = t_hit_sat_converged ? t_hit_sat : FT(Inf)
    t_hit_T_freeze = t_hit_T_freeze_converged ? t_hit_T_freeze : FT(Inf) 

    min_t, i_min_t = find_min_t(t_hit_sat, t_hit_T_freeze)

    if (min_t < Δt)
        # both are same
        S_ql = FT(0) # no liq to lose, cant grow in subsat
        S_qi = FT(0) # no ice to sublimate, can't grow in subsat

        new_ρ, new_p, new_T, new_q_vap, new_q, new_q_eq, new_q_sl, new_q_si, new_ts = morrison_milbrandt_2015_get_new_status_helper(param_set, p, q, q_eq, q_liq, q_ice, S_ql, S_qi, min_t, dqvdt, ts, w; error_on_q_neg = false)
        Δt_left = Δt - min_t
        new_regime = get_saturation_regime(new_q_vap, new_q.liq, new_q.ice, new_q_sl, new_q_si, new_T, T_freeze)
        (i_min_t == 1) && # @debug "Hit ice saturation before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
        (i_min_t == 2) && # @debug "Hit freezing before timestep is over... transitioning to $(typeof(new_regime)) at t = $(min_t)..."
        if i_min_t == 1 # hit sat
            new_δ_0_shum = BF ? TD.vapor_specific_humidity(new_q) - new_q_eq.liq : FT(0)
            new_δ_0i_shum = BF ? FT(0) : TD.vapor_specific_humidity(new_q) - new_q_eq.ice
        else # i_min_t == 2 # hit freezing
            new_δ_0_shum = TD.vapor_specific_humidity(new_q) - new_q_eq.liq
            new_δ_0i_shum = new_δ_0_shum # these are the same at freezing
        end
        S_ql_addit, S_qi_addit = morrison_milbrandt_2015_style(new_regime, param_set, area, new_ρ, new_p, new_T, w, τ_liq, τ_ice, new_q_vap, new_q, new_q_eq, Δt_left, new_ts; opts = MM2015Opts{FT}(use_fix = use_fix, return_mixing_ratio = true, depth = depth + 1, δ_0_shum = new_δ_0_shum, δ_0i_shum = new_δ_0i_shum, dqvdt=dqvdt, dTdt=dTdt, fallback_to_standard_supersaturation_limiter = fallback_to_standard_supersaturation_limiter))
        S_ql *= min_t / Δt # rescale to the timestep
        S_qi *= min_t / Δt # rescale to the timestep
        S_ql_addit *= Δt_left / Δt # rescale to the remaining time
        S_qi_addit *= Δt_left / Δt # rescale to the remaining time
        return return_mixing_ratio ? (S_ql + S_ql_addit, S_qi + S_qi_addit) : (S_mixing_ratio_to_shum(S_ql + S_ql_addit, q.tot), S_mixing_ratio_to_shum(S_qi + S_qi_addit, q.tot))

    else
        return FT(0), FT(0) 
    end
end

