


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
        const LowerSatLine = TC.LowerSatLine
        const LambertW = TC.LambertW
        const CMT = TC.CMT
        const CMP = TC.CMP
        include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style_exponential_part_only.jl"))
        include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style.jl"))

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
qls = qs
qis = qs

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

inputs = Iterators.product(Ts, ρs, qls, qis, limiters, τ_liqs, τ_ices, Δts, sat_factors, (:liqsat, :icesat), ws)

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
local new_ql::FT
local new_qi::FT
local δ_ql::FT
local δ_qi::FT
local δ_0::FT
local δ_0i::FT
local sat::Symbol

check_all_w = false

@showprogress showspeed=true for (T, ρ, ql, qi, limiter, τ_liq, τ_ice, Δt, sat_factor, sat, w) in inputs
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
    local qt = q_vap + ql + qi
    local q = TD.PhasePartition(qt, ql, qi)
    local q_eq = TD.PhasePartition(qt, q_sl, q_si)
    local ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
    local p = TD.air_pressure(thermo_params, ts)
    try
        
        q_vap = TD.vapor_specific_humidity(q) # for exactness
        (q_vap == TD.vapor_specific_humidity(q)) || error("q_vap != vapor_specific_humidity(q), got $q_vap = $q_vap; vapor_specific_humidity(q) = $(TD.vapor_specific_humidity(q))")
        S_ql, S_qi = calculate_timestep_limited_sources(limiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, q, q_eq, Δt, ts)

        δ_ql = S_ql * Δt
        δ_qi = S_qi * Δt
        new_ql = ql + δ_ql
        new_qi = qi + δ_qi
        isfinite(new_ql) || begin (@error "new_ql is not finite"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"); continue; end
        isfinite(new_qi) || begin (@error "new_qi is not finite"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"); continue; end
        (new_ql >= FT(0)) || begin (@error "new_ql is negative"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"); continue; end
        (new_qi >= FT(0)) || begin (@error "new_qi is negative"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"); continue; end

        (δ_ql ≥ -ql) || begin (@error "δ_ql is too negative"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"); continue; end
        (δ_qi ≥ -qi) || begin (@error "δ_qi is too negative"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"); continue; end

        if limiter ∈ (MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter, StandardSupersaturationMoistureSourcesLimiter, BasicMoistureSourcesLimiter)
            δ_0 = q_vap - q_sl
            δ_0i = q_vap - q_si
            # You can test these if you want but not you're probably always gonna get rounding errors on extreme cases 
            # ((δ_ql > FT(0)) ? (δ_ql ≤ (δ_0 + qi)) : true) || begin (@error "δ_ql is too large"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat_factor = $sat_factor; w = $w"); continue; end
            # ((δ_qi > FT(0)) ? (δ_qi ≤ (δ_0i + ql)) : true) || begin (@error "δ_qi is too large"); (@error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat_factor = $sat_factor; w = $w"); continue; end
        end

    catch e
        if isa(e, InterruptException)
            throw(e)
        end
        @error "Failed on T = $T; ρ = $ρ; ql = $ql; qi = $qi; limiter = $limiter; τ_liq = $τ_liq; τ_ice = $τ_ice; Δt = $Δt; sat = :$sat; sat_factor = $sat_factor; w = $w"
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
# qt = q_vap + ql + qi
# q = TD.PhasePartition(qt, ql, qi)
# q_vap = TD.vapor_specific_humidity(q) # for exactness
# q_eq = TD.PhasePartition(qt, q_sl, q_si)
# ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
# p = TD.air_pressure(thermo_params, ts)
# δ_0 = q_vap - q_sl
# δ_0i = q_vap - q_si

