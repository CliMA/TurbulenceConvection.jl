


const FT = Float64

using Pkg
using Statistics: mean

thisdir = dirname(@__FILE__)

reload_TC = false
if reload_TC || !isdefined(Main, :TurbulenceConvection) || !isdefined(Main, :param_set) || !isdefined(Main, :thermo_params)
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/"))
    # Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
    tc_dir = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl")
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/"))
    # include(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/src/TurbulenceConvection.jl")) 
    using ForwardDiff # for netcdf.io
    include(joinpath(tc_dir, "driver", "NetCDFIO.jl"))
    include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
    include(joinpath(tc_dir, "driver", "Cases.jl"))
    include(joinpath(tc_dir, "driver", "parameter_set.jl"))
    include(joinpath(tc_dir, "driver", "dycore.jl"))
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    namelist = NameList.default_namelist("SOCRATES_RF09_Obs_data"; write = false) # we don't need a directory
    namelist["root"] = thisdir
    namelist["output"]["output_root"] = thisdir
    namelist["meta"] = Dict("simname" => "test_mm_2015", "uuid" => "")

    const param_set = create_parameter_set(namelist, toml_dict, FT) # this creates an override file in  a directory we don' need...
    const thermo_params = TCP.thermodynamics_params(param_set)
end

call_from_TC = false
reload_MM2015 = true
if reload_MM2015 || !isdefined(Main, :morrison_milbrandt_2015_style) || !isdefined(Main, :morrison_milbrandt_2015_style_exponential_part_only) || !isdefined(Main, :limiters)
    if call_from_TC
        morrison_milbrandt_2015_style = TC.morrison_milbrandt_2015_style
        morrison_milbrandt_2015_style_exponential_part_only = TC.morrison_milbrandt_2015_style_exponential_part_only
    else
        # Things to load if we're calling directly form the script rather than TC methods
        const AbstractSaturationRegime = TC.AbstractSaturationRegime
        const Subsaturated = TC.Subsaturated
        const Supersaturated = TC.Supersaturated
        const WBF = TC.WBF
        # const LowerSatLine = TC.LowerSatLine # deprecated -- was supposed to be used for the WBF but in WBF you never actually get to the line and it's really hard to estimate the WBF rate on the line.
        const LambertW = TC.LambertW
        const CMT = TC.CMT
        const CMP = TC.CMP

        const AbstractSupersaturationLimiterMilestone = TC.AbstractSupersaturationLimiterMilestone
        const NotAtSupersaturationMilestone = TC.NotAtSupersaturationMilestone
        const OutOfLiquid = TC.OutOfLiquid
        const OutOfIce = TC.OutOfIce
        const AtSaturationOverLiquid = TC.AtSaturationOverLiquid
        const AtSaturationOverIce = TC.AtSaturationOverIce
        const AtSupersaturationStationaryPoint = TC.AtSupersaturationStationaryPoint

        const get_saturation_regime_type = TC.get_saturation_regime_type
        const get_saturation_regime = TC.get_saturation_regime
        const add_regime_parameters = TC.add_regime_parameters
        const get_new_saturation_regime_type_from_milestone = TC.get_new_regime_type_from_milestone
        const δi_from_δ = TC.δi_from_δ
        const δ_from_δi = TC.δ_from_δi
        const get_dδdt_0 = TC.get_dδdt_0
        const is_same_supersaturation_regime = TC.is_same_supersaturation_regime

        const safe_clamp = TC.safe_clamp

        include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style.jl"))
        include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style_exponential_part_only.jl"))

        const NoMoistureSourcesLimiter = TC.NoMoistureSourcesLimiter
        const BasicMoistureSourcesLimiter = TC.BasicMoistureSourcesLimiter
        const StandardSupersaturationMoistureSourcesLimiter = TC.StandardSupersaturationMoistureSourcesLimiter
        const MorrisonMilbrandt2015MoistureSourcesLimiter = TC.MorrisonMilbrandt2015MoistureSourcesLimiter
        const MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter = TC.MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter
        # const limiters = (NoMoistureSourcesLimiter, BasicMoistureSourcesLimiter, StandardSupersaturationMoistureSourcesLimiter, MorrisonMilbrandt2015MoistureSourcesLimiter, MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter)
        # limiters = (MorrisonMilbrandt2015MoistureSourcesLimiter,) # the most annoying one
        limiters = (MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter,) # the most annoying one
        # limiters = (NoMoistureSourcesLimiter, BasicMoistureSourcesLimiter, StandardSupersaturationMoistureSourcesLimiter, MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter)
        
        include(joinpath(tc_dir, "src", "microphysics_coupling_limiters.jl"))

        const Thermodynamics = TD
        Thermodynamics.PhasePartition{Float64}(x::Float64, y::Float64, z::Float64) = TD.PhasePartition(x,y,z)

        const area = FT(1.0)
    end
end


Ts = FT.((260, 273.1, 273.2, 275))
q_small = 10 .^ range(-200, stop = -50, length = 3)
q_large = 10 .^ range(-10, stop = -3, length = 5)
qs = (FT(0), q_small..., q_large...)
# qs = (FT(0),)
q_liqs = qs
q_ices = qs

ρs = FT(.7):FT(.3):FT(1.3)

τs_small = 10.0 .^ range(log10(floatmin(FT))-1, stop = -50, length = 3)
τs_medium = 10.0 .^ range(-10, stop = 10, length = 6)
τs_large = 10.0 .^ range(50, stop = log10(floatmax(FT))-1, length = 3)
τs = (τs_small..., τs_medium..., τs_large..., FT(Inf))
τ_liqs = τs
τ_ices = τs

Δts_small = 10.0 .^ range(-50, stop = -20, length = 3)
Δts_medium = 10.0 .^ range(-10, stop = -5, length = 3)
Δts_large = 10.0 .^ range(-2, stop = 0, length = 3)
Δts = (Δts_small..., Δts_medium..., Δts_large...)

signs = (-1, 1)
sat_factors = FT(10.0) .^ [-100, -10, -5, -1]
sat_factors = collect(sat_factor * sign for sign in signs for sat_factor in sat_factors) # collect so length is known
sat_factors = sort(collect((FT(0), sat_factors...)))
# sat_factors = (FT(0),)

ws = FT.((.1, 1.0, 10.0))
ws = collect(w * sign for sign in signs for w in ws) # collect so length is known
ws = sort(collect((FT(0), ws...)))
# ws = (FT(0),)

dqvdts = FT[0, 1e-7, -1e-7]

dTdts = FT[0, .002, -0.002]

inputs = Iterators.product(Ts, ρs, q_liqs, q_ices, limiters, τ_liqs, τ_ices, Δts, sat_factors, (:liqsat, :icesat), ws, dqvdts, dTdts)

using ProgressMeter

# inputs = collect(Iterators.take(inputs, 10000))

local q_si::FT
local q_sl::FT
local dδ::FT
local q_vap::FT
local qt::FT
local q::TD.PhasePartition{FT}
local q_eq::TD.PhasePartition{FT}
local ts::TD.ThermodynamicState{FT}
local p::FT
local S_ql::FT
local S_qi::FT
local new_q_liq::FT
local new_q_ice::FT
local δ_ql::FT
local δ_qi::FT
local δ_0::FT
local δ_0i::FT
local sat::Symbol

check_all_w = true

@showprogress showspeed=true for (T, ρ, q_liq, q_ice, limiter, τ_liq, τ_ice, Δt, sat_factor, sat, w, dqvdt, dTdt) in inputs
    if limiter ∈ (MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter, StandardSupersaturationMoistureSourcesLimiter, BasicMoistureSourcesLimiter)
        if !iszero(w) && !check_all_w # these limiters don't use w...
            continue
        end
    end
    local q_sl = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
    local q_si = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
    local dδ = q_sl - q_si
    if sat === :liqsat
        local q_vap = q_sl + (dδ * sat_factor)
    elseif sat === :icesat
        local q_vap = q_si + (dδ * sat_factor)
    else
        error("Invalid saturation regime: $sat")
    end
    local q_tot = q_vap + q_liq + q_ice
    local q = TD.PhasePartition(q_tot, q_liq, q_ice)
    local q_eq = TD.PhasePartition(q_tot, q_sl, q_si)
    local ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
    local p = TD.air_pressure(thermo_params, ts)
    try
        
        q_vap = TD.vapor_specific_humidity(q) # for exactness
        (q_vap == TD.vapor_specific_humidity(q)) || error("q_vap != vapor_specific_humidity(q), got $q_vap = $q_vap; vapor_specific_humidity(q) = $(TD.vapor_specific_humidity(q))")

        limiter_here =
        if limiter === MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter
            MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter(false) # don't fallback_to_standard_supersaturation_limiter bc BigFloat works most of the time and for super fast things we mostly fallback by default...
        else
            limiter()
        end

        
        S_ql, S_qi = calculate_timestep_limited_sources(limiter_here, param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, q, q_eq, Δt, ts)

        δ_ql = S_ql * Δt
        δ_qi = S_qi * Δt
        new_q_liq = q_liq + δ_ql
        new_q_ice = q_ice + δ_qi
        isfinite(new_q_liq) || begin (@error "new_q_liq $new_q_liq is not finite"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area"); continue; end
        isfinite(new_q_ice) || begin (@error "new_q_ice $new_q_ice is not finite"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area"); continue; end
        (new_q_liq >= FT(0)) || ( ((abs(new_q_liq/q_liq) < 1e-8) || (q_liq < eps(FT))) || begin (@error "new_q_liq $new_q_liq is negative"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area"); continue; end)
        (new_q_ice >= FT(0)) || ( ((abs(new_q_ice/q_ice) < 1e-8) || (q_ice < eps(FT))) || begin (@error "new_q_ice $new_q_ice is negative"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area"); continue; end)

        (δ_ql ≥ (-q_liq - 1e-8 * q_liq)) || ( ((abs(new_q_liq/q_liq) < 1e-8) || (q_liq < eps(FT))) || begin (@error "δ_ql $δ_ql is too negative for q_liq = $q_liq"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area"); continue; end)
        (δ_qi ≥ (-q_ice - 1e-8 * q_ice)) || ( ((abs(new_q_ice/q_ice) < 1e-8) || (q_ice < eps(FT))) || begin (@error "δ_qi $δ_qi is too negative for q_ice = $q_ice"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area"); continue; end)

        if limiter ∈ (MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter, StandardSupersaturationMoistureSourcesLimiter, BasicMoistureSourcesLimiter)
            δ_0 = q_vap - q_sl
            δ_0i = q_vap - q_si
            # You can test these if you want but not you're probably always gonna get rounding errors on extreme cases
            # ((δ_ql > FT(0)) ? (δ_ql ≤ (δ_0 + q_ice)) : true) || begin (@error "δ_ql is too large"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat_factor = $sat_factor; w = $w"); continue; end
            # ((δ_qi > FT(0)) ? (δ_qi ≤ (δ_0i + q_liq)) : true) || begin (@error "δ_qi is too large"); (@error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat_factor = $sat_factor; w = $w"); continue; end

            # dont let sources exceed available vapor
            if ((S_ql + S_qi)*Δt ≥ (q_vap+eps(FT))) && (!iszero(q_vap)) # consume all vapor (should never happen, really we should allow for liq and ice to contribute but all vapor gone is ridiculous, qv is like 2 orders of magnitude larger than q_liq, q_ice)
                error("Source calculations returned a total source greater than the available vapor. Got S_ql = $S_ql; S_qi = $S_qi; S_ql*Δt=$δ_ql; S_qi*Δt=$δ_qi while q_liq = $q_liq; q_ice = $q_ice; from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dqvdt = $dqvdt; dTdt = $dTdt; δ_ql = $δ_ql; δ_qi = $δ_qi; δ_0 = $δ_0; δ_0i = $δ_0i")
            elseif (S_ql + S_qi)*Δt < -(q.liq + q.ice + eps(FT)) # consume more than all liquid and ice
                error("Source calculations returned a total source less than the available liquid and ice. Got S_ql = $S_ql; S_qi = $S_qi from inputs area = $area; ρ = $ρ; p = $p; T = $T; w = $w; τ_liq = $τ_liq; τ_ice = $τ_ice; q_vap = $q_vap; q = $q; q_eq = $q_eq; Δt = $Δt; ts = $ts; dqvdt = $dqvdt; dTdt = $dTdt; δ_ql = $δ_ql; δ_qi = $δ_qi; δ_0 = $δ_0; δ_0i = $δ_0i")
            end
        end

    catch e
        if isa(e, InterruptException)
            throw(e)
        end

        if occursin("Timescales", sprint(showerror, e)) && occursin("are too fast", sprint(showerror, e)) && ((τ_liq < eps(FT)) || (τ_ice < eps(FT))) # ignore timescales too fast error in Standard() 
            # @warn "Ignoring known timescale error: $(sprint(showerror, e))"
            continue # keep going
        end

        δ_0 = q_vap - q_sl
        δ_0i = q_vap - q_si

        @error "Failed on T = $T; ρ = $ρ; q_liq = $q_liq; q_ice = $q_ice; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w; dqvdt = $dqvdt; dTdt = $dTdt; ts = $ts; p = $p; q_vap = $q_vap; q = $q; q_eq = $q_eq; area = $area; δ_0 = $δ_0; δ_0i = $δ_0i; q_sl = $q_sl; q_si = $q_si"
        throw(e)
    end


end


# q_sl = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_si = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
# dδ = q_sl - q_si
# if sat === :liqsat
#     q_vap = q_sl + (dδ * sat_factor)
# elseif sat === :icesat
#     q_vap = q_si + (dδ * sat_factor)
# else
#     error("Invalid saturation regime: $sat")
# end
# # q_vap = q_sl + (dδ * sat_factor)
# qt = q_vap + q_liq + q_ice
# q = TD.PhasePartition(qt, q_liq, q_ice)
# q_vap = TD.vapor_specific_humidity(q) # for exactness
# q_eq = TD.PhasePartition(qt, q_sl, q_si)
# ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
# p = TD.air_pressure(thermo_params, ts)
# δ_0 = q_vap - q_sl
# δ_0i = q_vap - q_si

