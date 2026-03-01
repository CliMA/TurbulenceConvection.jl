

# reload_SSCF = false
# if reload_SSCF | !isdefined(Main, :SOCRATESSingleColumnForcings)
#     using Pkg
#     using Revise
#     Pkg.activate(expanduser("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test"))
#     import CLIMAParameters as CP # use CLIMAParameters = "0.7, 0.8, 0.9, 0.10"
#     # import ClimaParams as CPP # would using this trouble w/ TC.jl? it's a different uuid technically..., use ClimaParams = "0.10"
#     import Thermodynamics as TD
#     import Thermodynamics.Parameters as TDP
#     # FT = Float64

#     toml_dict = CP.create_toml_dict(FT; dict_type = "alias") # CP 0.7 and below, Thermodynamics 0.11 and above
#     aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
#     param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
#     thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
#     Pkg.activate(expanduser("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/"))
#     include(
#         "/home/jbenjami/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/src/SOCRATESSingleColumnForcings.jl",
#     )
# end

# const FT = Float64
# FT = BigFloat
FT = Float64


using Pkg
using Statistics: mean
using Printf: @sprintf

thisdir = dirname(@__FILE__)

reload_TC = false
if reload_TC ||
   !isdefined(Main, :TurbulenceConvection) ||
   !isdefined(Main, :param_set) ||
   !isdefined(Main, :thermo_params) ||
   (typeof(thermo_params).parameters[1] != FT)
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

    param_set = create_parameter_set(namelist, toml_dict, FT) # this creates an override file in  a directory we don' need...
    thermo_params = TCP.thermodynamics_params(param_set)
end

call_from_TC = false
if call_from_TC
    morrison_milbrandt_2015_style = TC.morrison_milbrandt_2015_style
    morrison_milbrandt_2015_style_exponential_part_only = TC.morrison_milbrandt_2015_style_exponential_part_only
else
    # Things to load if we're calling directly form the script rather than TC methods
    const AbstractSaturationRegime = TC.AbstractSaturationRegime
    const CMP = TC.CMP
    const Subsaturated = TC.Subsaturated
    const Supersaturated = TC.Supersaturated
    const WBF = TC.WBF

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

    const LambertW = TC.LambertW
    const CMT = TC.CMT
    const safe_clamp = TC.safe_clamp
    include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style.jl"))
    include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style_exponential_part_only.jl"))

    const NoMoistureSourcesLimiter = TC.NoMoistureSourcesLimiter
    const BasicMoistureSourcesLimiter = TC.BasicMoistureSourcesLimiter
    const StandardSupersaturationMoistureSourcesLimiter = TC.StandardSupersaturationMoistureSourcesLimiter
    const MorrisonMilbrandt2015MoistureSourcesLimiter = TC.MorrisonMilbrandt2015MoistureSourcesLimiter
    const MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter =
        TC.MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter
    include(joinpath(tc_dir, "src", "microphysics_coupling_limiters.jl"))
end




T = FT(260)
# T = FT(272)
ρ = FT(0.5)



# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) + .5 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) - 0.5 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) - 5.0 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))


# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) + 0.05 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
q_vap_0 =
    TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) .* 1.0 +
    (
        TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) -
        TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
    ) .* 0.5

q_vap_0 *= 1


q_liq = FT(1e-8) * 1 * 25000
q_ice = FT(1e-8) * 1 * 10 * 100
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
q = TD.PhasePartition(q_vap_0 + q_liq + q_ice, q_liq, q_ice)


q_eq = TD.PhasePartition(
    q.tot,
    TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
    TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
)


dqvdt = FT(1e-5)
dTdt = FT(1e-5)

use_fix = true
# w = FT(+1000000.5)
w = FT(+0.0)

τ_liq = FT(1e-2)
τ_ice = FT(.935e-4)
# τ_ice = FT(.95e-4)

# τ_ice = FT(1e-2)
# τ_liq = FT(.935e-4)

τ_ice = FT(100.0)
τ_liq = FT(1.0)

τ_ice = FT(1.0)
τ_liq = FT(10.0)


# q_vap_0 = 0.004211891174058872
# q_tot = 0.004194225552472401
# q_liq = 
# q_ice = 0.
# T = 276.01926130123644
# τ_liq = 5.684678316404151


# τ_liq = FT(1e-20)
# τ_ice = FT(1e-13)

# τ_liq = FT(Inf)
# τ_ice = FT(Inf)

# q_vap_0, q_liq, q_ice = (0.0027351234431657536, 4.4465037826352686e-61, 2.094458880577658e-19)
# q_vap_0, q_liq, q_ice = (0.0027351234431657536, 4.4465037826352686e-4, 2.094458880577658e-4)
# q = TD.PhasePartition(q_vap_0, q_liq, q_ice)


# τ_liq = FT(Inf)
# τ_ice = FT(Inf)


# Δt = FT(0.0005)
# Δt = FT(1e-17)

area = FT(1.0)
# p = FT(750 .* 100)

# ##### An old failure case I had
Thermodynamics = TD
TD.PhasePartition{FT}(x::FT, y::FT, z::FT) where {FT} = TD.PhasePartition(x, y, z)



# convert everything to FT
q = TD.PhasePartition(FT(q.tot), FT(q.liq), FT(q.ice))
# q_eq = TD.PhasePartition(FT(q_eq.tot), FT(q_eq.liq), FT(q_eq.ice))

q_liq = q.liq
q_ice = q.ice
ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
q_vap = TD.vapor_specific_humidity(thermo_params, ts)
q_vap_0 = q_vap


#one line version
# (area, ρ, p, T, w, τ_liq, τ_ice, q_vap_0) = (FT(area), FT(ρ), FT(p), FT(T), FT(w), FT(τ_liq), FT(τ_ice), FT(q_vap_0)); q = TD.PhasePartition(FT(q.tot), FT(q.liq), FT(q.ice)); q_eq = TD.PhasePartition(FT(q_eq.tot), FT(q_eq.liq), FT(q_eq.ice)); ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q); q_vap = TD.vapor_specific_humidity(thermo_params, ts); q_vap_0 = q_vap


# #####

ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
p = TD.air_pressure(thermo_params, ts)
ρ = TD.air_density(thermo_params, ts)

(area, ρ, p, T, w, τ_liq, τ_ice, q_vap_0) = (FT(area), FT(ρ), FT(p), FT(T), FT(w), FT(τ_liq), FT(τ_ice), FT(q_vap_0))


q_vap = TD.vapor_specific_humidity(thermo_params, ts)
# @info(q_vap_0 - q_vap)


# ===================================================================== #
# S_ql = -1.722723211535336e-15; S_qi = -2.4881238809763347e-40
# area = 0.023915617886319285; ρ = 1.2324726009797524; p = 97835.32551757492; T = 275.7521203563053; w = 0.23720546040942375; τ_liq = 3.052296704e9; τ_ice = 4.503599627370496e15; q_vap = 0.004757075963271857; q = Thermodynamics.PhasePartition{Float64}(0.004757075963271857, 0.0, 1.244061940488168e-40); q_eq = Thermodynamics.PhasePartition{Float64}(0.004757075963271857, 0.0046979006207082, 0.004817189091832905); Δt = 0.5000000000000002; ts = Thermodynamics.PhaseNonEquil{Float64}(13165.01988635991, 1.232641626301441, Thermodynamics.PhasePartition{Float64}(0.004757075963271857, 0.0, 1.244061940488168e-40)); dqvdt = 0.0; dTdt = 0.0

S_ql = 1.5740282460754777e-6;
S_qi = 8.837601096511839e-8;
area = 0.049966362669001334;
ρ = 1.4607008779607433;
p = 74651.29711651805;
T = 178.1058319793084;
w = 0.018706561325269898;
τ_liq = 6.656910134239524e-3;
τ_ice = 0.14858242869377136;
q_vap = 5.095529671220612e-10;
q = Thermodynamics.PhasePartition{Float64}(0.00020172170899797012, 5.7675685444915575e-9, 0.0002017154318764585);
q_eq = Thermodynamics.PhasePartition{Float64}(0.00020172170899797012, 7.545234665120191e-8, 2.9644619882828738e-8);
Δt = 0.007696511708622092;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    -68295.73889316698,
    1.4607008779607433,
    Thermodynamics.PhasePartition{Float64}(0.00020172170899797012, 5.7675685444915575e-9, 0.0002017154318764585),
);
dqvdt = 1.2898707234822461e-6;
dTdt = 0.01384903619503847;

# S_ql = 0.24484948201431436; S_qi = 0.0002609807690189025
# area = 0.04998396610455697; ρ = 1.1004823628116984; p = 81881.18225822448; T = 247.41072508045403; w = 0.01157821665104497; τ_liq = 0.0007316029514186084; τ_ice = 1.1408493518829346; q_vap = 0.07880207008012662; q = Thermodynamics.PhasePartition{Float64}(0.07887520480324778, 0.0, 7.313472312116094e-5); q_eq = Thermodynamics.PhasePartition{Float64}(0.07887520480324778, 0.0006020599833018769, 0.0004676899808379902); Δt = 20.0; ts = Thermodynamics.PhaseNonEquil{Float64}(167251.6888354985, 1.1004823628116984, Thermodynamics.PhasePartition{Float64}(0.07887520480324778, 0.0, 7.313472312116094e-5)); dqvdt = -2.4576585297496547e-7; dTdt = -0.0005850420073265253

T = 273.1;
ρ = 0.7;
q_liq = 1.778279410038923e-5;
q_ice = 0.0;
τ_liq = 2.22507385850719e-309;
τ_ice = 2.22507385850719e-309;
Δt = 1.0e-50;
sat = :liqsat;
sat_factor = -0.1;
w = -10.0;
dqvdt = 0.0;
dTdt = 0.0;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    16347.138620801497,
    0.7,
    Thermodynamics.PhasePartition{Float64}(0.006919818397398777, 1.778279410038923e-5, 0.0),
);
p = 55095.54778139178;
q_vap = 0.006902035603298388;
q = Thermodynamics.PhasePartition{Float64}(0.006919818397398777, 1.778279410038923e-5, 0.0);
q_eq = Thermodynamics.PhasePartition{Float64}(0.006919818397398777, 0.006902436691973973, 0.006898425805218117);
area = 1.0;
# T = 273.1; ρ = 1.3; q_liq = 1.778279410038923e-5; q_ice = 0.0; τ_liq = 4.717068855239647e-180; τ_ice = 2.22507385850719e-309; Δt = 1.0e-50; sat = :liqsat; sat_factor = -0.1; w = -10.0; dqvdt = 0.0; dTdt = 0.0; ts = Thermodynamics.PhaseNonEquil{Float64}(8782.43430804145, 1.3, Thermodynamics.PhasePartition{Float64}(0.003734263503568751, 1.778279410038923e-5, 0.0)); p = 102122.9262908164; q_vap = 0.003716480709468362; q = Thermodynamics.PhasePartition{Float64}(0.003734263503568751, 1.778279410038923e-5, 0.0); q_eq = Thermodynamics.PhasePartition{Float64}(0.003734263503568751, 0.0037166966802936775, 0.003714536972040524); area = 1.0
# T = 273.2; ρ = 1.3; q_liq = 0.0; q_ice = 1.778279410038923e-5;  τ_liq = 4.717068855239647e-180; τ_ice = 2.22507385850719e-309; Δt = 1.0e-50; sat = :liqsat; sat_factor = -0.1; w = -10.0; dqvdt = 0.0; dTdt = 0.0; ts = Thermodynamics.PhaseNonEquil{Float64}(8910.447241624219, 1.3, Thermodynamics.PhasePartition{Float64}(0.0037603450205353195, 0.0, 1.778279410038923e-5)); p = 102161.93685917695; q_vap = 0.0037425622264349303; q = Thermodynamics.PhasePartition{Float64}(0.0037603450205353195, 0.0, 1.778279410038923e-5); q_eq = Thermodynamics.PhasePartition{Float64}(0.0037603450205353195, 0.003742417187122792, 0.0037438675802441753); area = 1.0

T = 273.2;
ρ = 1.0;
q_liq = 0.0;
q_ice = 1.0e-10;
τ_liq = 1.0e-50;
τ_ice = 1.0e-10;
Δt = 0.010000000000000002;
sat = :liqsat;
sat_factor = -0.1;
w = -10.0;
dqvdt = 0.0;
dTdt = 0.0;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    11582.682072221132,
    1.0,
    Thermodynamics.PhasePartition{Float64}(0.00486533099436541, 0.0, 1.0e-10),
);
p = 78641.0320095916;
q_vap = 0.00486533089436541;
q = Thermodynamics.PhasePartition{Float64}(0.00486533099436541, 0.0, 1.0e-10);
q_eq = Thermodynamics.PhasePartition{Float64}(0.00486533099436541, 0.00486514234325963, 0.004867027854317428);
area = 1.0;

area = 1.0;
ρ = 0.7;
p = 55512.5292119373;
T = 275.0;
w = -10.0;
τ_liq = 0.010000000000000002;
τ_ice = 100.0;
q_vap = 0.007875108050448738;
q = Thermodynamics.PhasePartition{Float64}(0.007875108150448737, 0.0, 1.0e-10);
q_eq = Thermodynamics.PhasePartition{Float64}(0.007875108150448737, 0.00786098494726133, 0.00800221597913541);
Δt = 1.0;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    20031.323030729105,
    0.7,
    Thermodynamics.PhasePartition{Float64}(0.007875108150448737, 0.0, 1.0e-10),
);
dqvdt = 0.0;
dTdt = 0.0;
dqsl_dT = 0.0005336975485260993;
dqsi_dT = 0.0005432859978098328;

area = 1.0;
ρ = 1.3;
p = 102161.92786926191;
T = 273.2;
w = -10.0;
τ_liq = 2.22507385850719e-309;
τ_ice = 2.22507385850719e-309;
q_vap = 0.003742417187122792;
q = Thermodynamics.PhasePartition{Float64}(0.0037601999812231814, 0.0, 1.778279410038923e-5);
q_eq = Thermodynamics.PhasePartition{Float64}(0.0037601999812231814, 0.003742417187122792, 0.0037438675802441753);
Δt = 1.0e-50;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    8910.102808405085,
    1.3,
    Thermodynamics.PhasePartition{Float64}(0.0037601999812231814, 0.0, 1.778279410038923e-5),
);
dqvdt = 0.0;
dTdt = 0.0;
δ_ql = 3.969202437239124e252;
δ_qi = -8.831768582223793e-7;
δ_0 = 0.0;
δ_0i = -1.4503931213831651e-6;

area = 1.0;
ρ = 0.7;
p = 55118.30555399315;
T = 273.2;
w = -10.0;
τ_liq = 1.0e-6;
τ_ice = 100.0;
q_vap = 0.006950203347513758;
q = Thermodynamics.PhasePartition{Float64}(0.006950203347513758, 1.0e-200, 1.0e-200);
q_eq = Thermodynamics.PhasePartition{Float64}(0.006950203347513758, 0.006950203347513758, 0.006952896934739183);
Δt = 3.162277660168379e-8;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    16533.748957298405,
    0.7,
    Thermodynamics.PhasePartition{Float64}(0.006950203347513758, 1.0e-200, 1.0e-200),
);
dqvdt = 0.0;
dTdt = 0.0;
dqsl_dT = 0.0004791120152533068;
dqsi_dT = 0.00047929769759083086;

area = 1.0;
ρ = 0.7;
p = 55512.054696767365;
T = 275.0;
w = -1.0;
τ_liq = 2.22507385850719e-309;
τ_ice = 0.010000000000000002;
q_vap = 0.00786098353495101;
q = Thermodynamics.PhasePartition{Float64}(0.00786098363495101, 0.0, 1.0e-10);
q_eq = Thermodynamics.PhasePartition{Float64}(0.00786098363495101, 0.00786098494726133, 0.00800221597913541);
Δt = 1.0e-5;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    19997.763441520052,
    0.7,
    Thermodynamics.PhasePartition{Float64}(0.00786098363495101, 0.0, 1.0e-10),
);
dqvdt = 0.0;
dTdt = 0.0;
dqsl_dT = 0.0005336975485260993;
dqsi_dT = 0.0005432859978098328;

area = 1.0;
ρ = 0.7;
p = 55118.314538419516;
T = 273.2;
w = -10.0;
τ_liq = 2.22507385850719e-309;
τ_ice = 100.0;
q_vap = 0.0069504727062363;
q = Thermodynamics.PhasePartition{Float64}(0.0069504728062363, 0.0, 1.0e-10);
q_eq = Thermodynamics.PhasePartition{Float64}(0.0069504728062363, 0.006950203347513758, 0.006952896934739183);
Δt = 0.010000000000000002;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    16534.388585636614,
    0.7,
    Thermodynamics.PhasePartition{Float64}(0.0069504728062363, 0.0, 1.0e-10),
);
dqvdt = 0.0;
dTdt = 0.0;
dqsl_dT = 0.0004791120152533068;
dqsi_dT = 0.00047929769759083086;

T = 273.2;
ρ = 0.7;
q_liq = 0.0;
q_ice = 1.778279410038923e-5;
τ_liq = 2.22507385850719e-309;
τ_ice = 2.22507385850719e-309;
Δt = 1.0e-50;
sat = :liqsat;
sat_factor = -0.1;
w = 0.0;
dqvdt = 0.0;
dTdt = 0.0;
ts = Thermodynamics.PhaseNonEquil{Float64}(
    16528.457262263426,
    0.7,
    Thermodynamics.PhasePartition{Float64}(0.00696825550033669, 0.0, 1.778279410038923e-5),
);
p = 55117.33851141129;
q_vap = 0.0069504727062363;
q = Thermodynamics.PhasePartition{Float64}(0.00696825550033669, 0.0, 1.778279410038923e-5);
q_eq = Thermodynamics.PhasePartition{Float64}(0.00696825550033669, 0.006950203347513758, 0.006952896934739183);
area = 1.0;

# -------------------------------------- #
q_liq = q.liq;
q_ice = q.ice;
q_vap_0 = q_vap;
# ===================================================================== #



# Dict(
#     "Normal" => morrison_milbrandt_2015_style(param_set, area, ρ,p, T, w, τ_liq, τ_ice, q_vap_0, q, q_eq, Δt, ts),
#     "Exp Only" => morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, q_eq, Δt,)
#     )



function symlog(x, n = -3)
    result = zeros(size(x, 1))
    for i in eachindex(x)
        result[i] = sign(x[i]) * (log10(1 + abs(x[i]) / (10^n)))
    end
    result
end

function symlogformatter(x, n; ndigits = 20)
    # if sign(x) == 0
    #     "10^$(n)"
    # else
    s = sign(x) == 1 ? "+" : "-"
    nexp = sign(x) * (abs(x) + n)
    if sign(x) == -1
        nexp = -nexp
    end
    n_exp_print = round(nexp, digits = ndigits)
    # get minimal
    # s * "10^$(n_exp_print)"
    num_str = @sprintf("%.*f", ndigits, nexp)  # ✅ dynamic precision
    return s * "10^" * num_str

    # end
end

# Δts = 10 .^ range(-7, stop = 5, length = 100)
# Δts = FT(10) .^ range(FT(-4), stop = FT(2), length = 100); @warn("temp test")
# Δts = 10 .^ range(log10(5), stop = log10(10), length = 100); @warn("temp test")

# Δts = FT(10) .^ range(-17, stop = 7, length = 1000); @warn("Really temp test")

Δts = 10 .^ range(-7, stop = 5, length = 1000)


qls = zeros(length(Δts))
qis = zeros(length(Δts))
qls_exp = zeros(length(Δts))
qis_exp = zeros(length(Δts))
qls_standard = zeros(length(Δts))
qis_standard = zeros(length(Δts))
δ_ΔTs = zeros(length(Δts))
using Plots
ENV["GKSwstype"] = "nul"
using ProgressMeter
@showprogress dt = 1 desc = "Computing..." for (i, Δt) in enumerate(Δts)
    # qls[i], qis[i], δ_ΔTs[i] = morrison_milbrandt_2015_style(param_set, area, ρ,p, T, w, τ_liq, τ_ice, q_vap_0, q, q_eq, Δt, ts; use_fix=use_fix)
    # qls_exp[i], qis_exp[i] = morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q, q_eq, Δt,)

    println(
        "============================================================================================================================================================== Δt = $Δt",
    )

    qls[i], qis[i] = (q_vap_0 .* FT(NaN), q_vap_0 * FT(NaN))
    # morrison_milbrandt_2015_style( param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q, q_eq, Δt, ts; use_fix = use_fix)

    println(
        "----------------------------------------------------------------------------------------------------------------------------------------------------- Δt = $Δt (exp only)",
    )

    qls_exp[i], qis_exp[i] =
    # (q_vap_0 .* FT(NaN), q_vap_0 * FT(NaN))
        morrison_milbrandt_2015_style_exponential_part_only(
            param_set,
            area,
            ρ,
            p,
            T,
            w,
            τ_liq,
            τ_ice,
            q_vap_0,
            dqvdt,
            dTdt,
            q,
            q_eq,
            Δt,
            ts;
            opts = MM2015Opts{Float64}(use_fix = use_fix, fallback_to_standard_supersaturation_limiter = false),
        )

    println(
        "----------------------------------------------------------------------------------------------------------------------------------------------------- Δt = $Δt (exp only)",
    )

    qls_standard[i], qis_standard[i] =
    # (q_vap_0 .* FT(NaN), q_vap_0 * FT(NaN))
        calculate_timestep_limited_sources(
            TC.StandardSupersaturationMoistureSourcesLimiter(),
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

    qls[i], qis[i] = qls_exp[i], qis_exp[i]

end
# plot the results, x scale log
linthresh = FT(-5)
valid = abs.(qls .+ qis) .> 0
linthresh = mean(log10.((abs.((qls .+ qis) .* Δts))[valid])) - 1
yformatter = x -> symlogformatter.(x, linthresh, ndigits = 4)
yscale = :symlog
# yscale = :linear
if yscale === :symlog
    plot(
        Δts,
        symlog(qls .* Δts, linthresh),
        label = "S_ql",
        dpi = 600,
        size = (1000, 400),
        xaxis = :log,
        legend = :outertopright,
        color = :green,
        yformatter = yformatter,
        minorgrid = true,
        linewidth = 2,
        # left_margin=40Plots.PlotMeasures.mm, # sometimes we need it sometimes not
    )
    plot!(Δts, symlog(qis .* Δts, linthresh), label = "S_qi", color = :blue, yformatter = yformatter, linewidth = 2)
    plot!(
        Δts,
        symlog((qls .+ qis) .* Δts, linthresh),
        label = "sum",
        color = :red,
        yformatter = yformatter,
        linewidth = 1,
    )
    #
    plot!(
        Δts,
        symlog(qls_exp .* Δts, linthresh),
        label = "S_ql_exp",
        color = :green,
        linestyle = :dash,
        yformatter = yformatter,
        linewidth = 2,
    )
    plot!(
        Δts,
        symlog(qis_exp .* Δts, linthresh),
        label = "S_qi_exp",
        color = :blue,
        linestyle = :dash,
        yformatter = yformatter,
        linewidth = 2,
    )
    plot!(
        Δts,
        symlog((qls_exp .+ qis_exp) .* Δts, linthresh),
        label = "sum_exp",
        color = :red,
        linestyle = :dash,
        yformatter = yformatter,
        linewidth = 1,
    )
    #
    plot!(
        Δts,
        symlog(qls_standard .* Δts, linthresh),
        label = "S_ql_standard",
        color = :green,
        linestyle = :dot,
        yformatter = yformatter,
        linewidth = 2,
    )
    plot!(
        Δts,
        symlog(qis_standard .* Δts, linthresh),
        label = "S_qi_standard",
        color = :blue,
        linestyle = :dot,
        yformatter = yformatter,
        linewidth = 2,
    )
    plot!(
        Δts,
        symlog((qls_standard .+ qis_standard) .* Δts, linthresh),
        label = "sum_standard",
        color = :red,
        linestyle = :dot,
        yformatter = yformatter,
        linewidth = 1,
    )


    hline!([0], color = :black, linestyle = :dash, yformatter = yformatter, label = "")
    hline!(
        symlog([(q_vap_0 - q_eq.liq)], linthresh),
        color = :darkgreen,
        linestyle = :dashdot,
        yformatter = yformatter,
        label = "δ_l",
        alpha = 0.5,
    )
    hline!(
        symlog([q_vap_0 - q_eq.ice], linthresh),
        color = :darkblue,
        linestyle = :dashdot,
        yformatter = yformatter,
        label = "δ_i",
        alpha = 0.5,
    )
    # plot!(Δts, symlog(-δ_ΔTs .* Δts, linthresh), label = "δ_ΔT", color=:orange, linestyle=:dashdot, yformatter=yformatter)
    plot!(
        Δts,
        symlog((q_vap_0 .- q_eq.ice) ./ τ_ice .* Δts, linthresh),
        label = "S_i base",
        color = :darkblue,
        linestyle = :solid,
        yformatter = yformatter,
    )
    plot!(
        Δts,
        symlog((q_vap_0 .- q_eq.liq) ./ τ_liq .* Δts, linthresh),
        label = "S_l base",
        color = :darkgreen,
        linestyle = :solid,
        yformatter = yformatter,
    )

    plot!(
        Δts,
        symlog(fill(-q.liq, size(Δts)), linthresh),
        label = "-ql/Δt base",
        color = :darkgreen,
        linestyle = :solid,
        yformatter = yformatter,
    )
    plot!(
        Δts,
        symlog(fill(-q.ice, size(Δts)), linthresh),
        label = "-qi/Δt base",
        color = :darkblue,
        linestyle = :solid,
        yformatter = yformatter,
    )
elseif yscale === :linear
    plot(
        Δts,
        qls .* Δts,
        label = "S_ql",
        dpi = 600,
        size = (1000, 400),
        xaxis = :log,
        yaxis = :linear,
        legend = :outertopright,
        color = :green,
    )
    plot!(Δts, qis .* Δts, label = "S_qi", color = :blue)
    plot!(Δts, (qls .+ qis) .* Δts, label = "sum", color = :red)
    #
    plot!(Δts, qls_exp .* Δts, label = "ql_exp", color = :green, linestyle = :dash)
    plot!(Δts, qis_exp .* Δts, label = "qi_exp", color = :blue, linestyle = :dash)
    plot!(Δts, (qls_exp .+ qis_exp) .* Δts, label = "sum_exp", color = :red, linestyle = :dash)
    #
    plot!(Δts, qls_standard .* Δts, label = "ql_standard", color = :green, linestyle = :dot)
    plot!(Δts, qis_standard .* Δts, label = "qi_standard", color = :blue, linestyle = :dot)
    plot!(Δts, (qls_standard .+ qis_standard) .* Δts, label = "sum_standard", color = :red, linestyle = :dot)
    #
    hline!([0], color = :black, linestyle = :dash, label = "")
    hline!([q_vap_0 - q_eq.liq], color = :gray, linestyle = :dashdot, yformatter = yformatter, label = "δ_l")
    plot!(Δts, (q_vap_0 .- q_eq.ice) ./ τ_ice .* Δts, label = "S_i base", color = :pink, linestyle = :solid)


end

nanminimum = x -> minimum(x[isfinite.(x)])
nanmaximum = x -> maximum(x[isfinite.(x)])
nanextrema = x -> (nanminimum(x), nanmaximum(x))

ymin_l =
    yscale === :symlog ?
    nanminimum(symlog(((qls .* Δts)..., (qls_exp .* Δts)..., (qls_standard .* Δts)...), linthresh)) :
    nanminimum(((qls .* Δts)..., (qls_exp .* Δts)..., (qls_standard .* Δts)...))
ymin_i =
    yscale === :symlog ?
    nanminimum(symlog(((qis .* Δts)..., (qis_exp .* Δts)..., (qis_standard .* Δts)...), linthresh)) :
    nanminimum(((qis .* Δts)..., (qis_exp .* Δts)..., (qis_standard .* Δts)...))
ymax_l =
    yscale === :symlog ?
    nanmaximum(symlog(((qls .* Δts)..., (qls_exp .* Δts)..., (qls_standard .* Δts)...), linthresh)) :
    nanmaximum(((qls .* Δts)..., (qls_exp .* Δts)..., (qls_standard .* Δts)...))
ymax_i =
    yscale === :symlog ?
    nanmaximum(symlog(((qis .* Δts)..., (qis_exp .* Δts)..., (qis_standard .* Δts)...), linthresh)) :
    nanmaximum(((qis .* Δts)..., (qis_exp .* Δts)..., (qis_standard .* Δts)...))
ymin = min(ymin_l, ymin_i)
ymax = max(ymax_l, ymax_i)

if yscale === :symlog
    dy = ymax - ymin
    factor = 0.2
    ymin -= factor * dy
    ymax += factor * dy
    # ymin = ymin > 0 ? ymin * 0.9 : ymin * 1.1
    # ymax = ymax > 0 ? ymax * 1.1 : ymax * 0.9
elseif yscale == :linear
    dy = ymax - ymin
    ymin -= 0.1 * dy
    ymax += 0.1 * dy
end

ylims!(ymin, ymax)


vline!([Δts[1]], color = :purple, linestyle = :dash, label = "")
vline!([Δts[end]], color = :purple, linestyle = :dash, label = "")

# hline!([-q_liq], color=:green, linestyle=:dashdot, label="-q_liq(t=0)")
# hline!([-q_ice], color=:blue, linestyle=:dashdot, label="-q_ice(t=0)")

# save pdf and png
for fmt in [:pdf, :png]
    savefig(joinpath(thisdir, "test_moisture_sources_limiters.$fmt"))
end
# savefig(joinpath(thisdir, "test_moisture_sources_limiters.png"))
