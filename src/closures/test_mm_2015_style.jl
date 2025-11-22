

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

if !isdefined(Main, :FT)
    const FT = Float64
end

using Pkg
using Statistics: mean
using Printf: @sprintf
using ProgressMeter
using Plots
ENV["GKSwstype"]="nul"

thisdir = dirname(@__FILE__)

reload_TC = false
if reload_TC || !isdefined(Main, :TurbulenceConvection) || !isdefined(Main, :param_set) || !isdefined(Main, :thermo_params)
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/"))
    using OrderedCollections: OrderedCollections # is only available from /integration_tests/ so we have to load it here
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
reload_MM2015 = false

if reload_MM2015 || !isdefined(Main, :morrison_milbrandt_2015_style) || !isdefined(Main, :morrison_milbrandt_2015_style_exponential_part_only)
    if call_from_TC
        morrison_milbrandt_2015_style = TC.morrison_milbrandt_2015_style
        morrison_milbrandt_2015_style_exponential_part_only = TC.morrison_milbrandt_2015_style_exponential_part_only
    else
        # Things to load if we're calling directly form the script rather than TC methods
        const AbstractSaturationRegime = TC.AbstractSaturationRegime
        const Subsaturated = TC.Subsaturated
        const Supersaturated = TC.Supersaturated
        const WBF = TC.WBF
        # const LowerSatLine = TC.LowerSatLine
        const LambertW = TC.LambertW
        const CMT = TC.CMT
        const CMP = TC.CMP
        const safe_clamp = TC.safe_clamp
        const get_qv_eq_point = TC.get_qv_eq_point
        const calculate_timestep_limited_sources = TC.calculate_timestep_limited_sources
        const calculate_next_standard_milestone_time_given_milestones = TC.calculate_next_standard_milestone_time_given_milestones
        const step = TC.step
        const calculate_next_standard_milestone_time = TC.calculate_next_standard_milestone_time
        const StandardSupersaturationMoistureSourcesLimiter = TC.StandardSupersaturationMoistureSourcesLimiter
        const t_out_of_x = TC.t_out_of_x

        const AbstractSupersaturationLimiterMilestone = TC.AbstractSupersaturationLimiterMilestone
        const NotAtSupersaturationMilestone = TC.NotAtSupersaturationMilestone
        const OutOfLiquid = TC.OutOfLiquid
        const OutOfIce = TC.OutOfIce
        const AtSaturationOverLiquid = TC.AtSaturationOverLiquid
        const AtSaturationOverIce = TC.AtSaturationOverIce
        const AtSupersaturationStationaryPoint = TC.AtSupersaturationStationaryPoint

        
        const resolve_S_S_addit = TC.resolve_S_S_addit
        const get_saturation_regime_type = TC.get_saturation_regime_type
        const get_saturation_regime = TC.get_saturation_regime
        const add_regime_parameters = TC.add_regime_parameters
        const get_new_saturation_regime_type_from_milestone = TC.get_new_regime_type_from_milestone
        const δi_from_δ = TC.δi_from_δ
        const δ_from_δi = TC.δ_from_δi
        const get_dδdt_0 = TC.get_dδdt_0
        const is_same_supersaturation_regime = TC.is_same_supersaturation_regime

        include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style.jl"))
        include(joinpath(tc_dir, "src", "closures", "morrison_milbrandt_2015_style_exponential_part_only.jl"))

        const Thermodynamics = TD
        Thermodynamics.PhasePartition{Float64}(x::Float64, y::Float64, z::Float64) = TD.PhasePartition(x,y,z)
    end
end




# T = FT(260.8)
# T = FT(275)
# ρ = FT(0.5)
# ρ = FT(0.99)
# RH = FT(1.01)
# RH = FT(1.12)

# T = FT(254.58)
# ρ = FT(1.0195)
# p = FT(74651.)
# RH = FT(1.7031)

# # T = FT(258.62)
# # ρ = FT(1.0597)
# # p = FT(78720.)
# # RH = FT(1.3249)


# dq = 0.000099182 * 5
# # dq = 0.000085828 * 5
# # dq = 0.0006537095128834309

# dT = dq * L_l/c_p
# dqls = dqsl_dT * dq * L_l/c_p

# new_RH = (q_vap_0 - dq) / TD.q_vap_saturation_generic(thermo_params, T + dT, ρ, TD.Liquid())



# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) + .5 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) - 0.5 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) - 5.0 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))


# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) + 0.05 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 =
#     TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) .* (RH) +
#     (
#         TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) -
#         TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
#     ) .* 0.0005 * 0

# q_vap_0 *= 1


# q_liq = FT(5e-5) * 1 * 1
# q_ice = FT(7e-8) * 1 * 1

# q_ice = 1.748459986808224e-110
# q_ice = 1e-15
# q_liq = 7.868468755570234e-150

# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
# q = TD.PhasePartition(q_vap_0 + q_liq + q_ice, q_liq, q_ice)

use_fix = false
# w = FT(+1000000.5)
# w = FT(-0.0)
# dqvdt = FT(-2e-9)
# dTdt = FT(-3e-5)
# dTdt = FT(-3e-2)


# τ_ice = FT(1e8)
# τ_liq = FT(1e2)


# area = FT(1.0)
# p = FT(750 .* 100)

# ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
# p = TD.air_pressure(thermo_params, ts)
# ρ = TD.air_density(thermo_params, ts)

# q_vap = TD.vapor_specific_humidity(thermo_params, ts)
# # @info(q_vap_0 - q_vap)
# q_eq = TD.PhasePartition(
#     q.tot,
#     TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
#     TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
# )



test_loop = true
# S_ql = 2.6648728745568e-7; S_qi = -5.4216454135895356e-5; area = 0.9345002260562765; ρ = 1.2178189520267437; p = 96153.27798188626; T = 274.3873889742498; w = -0.06924515117191868; τ_liq = 0.017190976068377495; τ_ice = 0.5335972309112549; q_vap = 0.004334269514558791; q = Thermodynamics.PhasePartition{Float64}(0.004360670191519954, 2.6400676961162774e-5, 6.462348535570529e-27); q_eq = Thermodynamics.PhasePartition{Float64}(0.004360670191519954, 0.004333936405449472, 0.004385739554245602); Δt = 1.2500000000000002; ts = Thermodynamics.PhaseNonEquil{Float64}(11177.110439755896, 1.2178189520267437, Thermodynamics.PhasePartition{Float64}(0.004360670191519954, 2.6400676961162774e-5, 6.462348535570529e-27)); dqvdt = -4.8235384433029e-8; dTdt = 0.00012646893809203843
# S_ql = -7.2292665105809514e-9; S_qi = 6.615638565144957e-6; area = 0.5851509522828972; ρ = 1.6513423208935223; p = 99997.59915358007; T = 211.09592611590728; w = 0.0; τ_liq = 587.9485473632812; τ_ice = 4.0064773202175274e-5; q_vap = 4.038688636945147e-5; q = Thermodynamics.PhasePartition{Float64}(0.0005551210324224197, 8.202420379429755e-7, 0.0005139139040150252); q_eq = Thermodynamics.PhasePartition{Float64}(0.0005551210324224197, 9.287251731221621e-6, 5.023674409734277e-6); Δt = 10.0; ts = Thermodynamics.PhaseNonEquil{Float64}(-44652.856741360054, 1.6513423208935223, Thermodynamics.PhasePartition{Float64}(0.0005551210324224197, 8.202420379429755e-7, 0.0005139139040150252)); dqvdt = 3.157688201160764e-6; dTdt = 0.05699837957444119
# area = 0.04999715734476376; ρ = 1.2020954015651961; p = 94829.59591704236; T = 274.18050776832496; w = 0.010000033432521083; τ_liq = 1.4250951529959366e-9; τ_ice = 2015.6162109375; δ_0 = 0.0; δ_0i = 0.043247192979746156; q = Thermodynamics.PhasePartition{Float64}(0.004455928315161662, 0.0001292893116379244, 0.0); q_eq = Thermodynamics.PhasePartition{Float64}(0.004455928315161662, 0.004329017736538598, 0.004372002122178419); Δt = 1.2499999461209619; ts = Thermodynamics.PhaseNonEquil{Float64}(11002.616314194878, 1.2020954015651961, Thermodynamics.PhasePartition{Float64}(0.004455928315161662, 0.00013252088726586428, 0.0))
# area = 0.04999715538960087; ρ = 1.2047613974582667; p = 94965.32251252233; T = 273.99443472463236; w = 0.00999980070533047; τ_liq = 1.0265231997763635e-9; τ_ice = 2295.97998046875; δ_0 = 0.0; δ_0i = 0.04261337201761381; q = Thermodynamics.PhasePartition{Float64}(0.0044539389202048005, 0.00019285925715926635, 0.0); q_eq = Thermodynamics.PhasePartition{Float64}(0.0044539389202048005, 0.004264796979900834, 0.004299394761295365); Δt = 1.2499999525146464; ts = Thermodynamics.PhaseNonEquil{Float64}(10708.462730601326, 1.2047613974582667, Thermodynamics.PhasePartition{Float64}(0.0044539389202048005, 0.00019797642446978285, 0.0))
# area = 0.043961462006736535; ρ = 0.9876874310242214; p = 94433.03794505185; T = 337.2494197169845; w = 1.1397681035854128; τ_liq = 9.578226232986875e-12; τ_ice = 7.7581320192e11; δ_0 = 0.0; δ_0i = 1.4277671685089102; q = Thermodynamics.PhasePartition{Float64}(0.016795466635614025, 0.004537123505354369, 0.0); q_eq = Thermodynamics.PhasePartition{Float64}(0.016795466635614025, 0.15545313779022987, 0.28221734718361857); Δt = 9.999999998499758; ts = Thermodynamics.PhaseNonEquil{Float64}(55985.072457444214, 0.9876874310242214, Thermodynamics.PhasePartition{Float64}(0.016795466635614025, 0.013939627662109452, 0.0)); dTdt = FT(0); dqvdt = FT(0.0)
S_ql = 1.5740282460754777e-6; S_qi = 8.837601096511839e-8; area = 0.049966362669001334; ρ = 1.4607008779607433; p = 74651.29711651805; T = 178.1058319793084; w = 0.018706561325269898; τ_liq = 6.656910134239524e-9; τ_ice = 0.14858242869377136; q_vap = 5.095529671220612e-10; q = Thermodynamics.PhasePartition{Float64}(0.00020172170899797012, 5.7675685444915575e-9, 0.0002017154318764585); q_eq = Thermodynamics.PhasePartition{Float64}(0.00020172170899797012, 7.545234665120191e-8, 2.9644619882828738e-8); Δt = 0.007696511708622092; ts = Thermodynamics.PhaseNonEquil{Float64}(-68295.73889316698, 1.4607008779607433, Thermodynamics.PhasePartition{Float64}(0.00020172170899797012, 5.7675685444915575e-9, 0.0002017154318764585)); dqvdt = 1.2898707234822461e-5; dTdt = 0.01384903619503847

area = 0.09999999999999998; ρ = 1.1297764844315694; p = 72678.5848761201; T = 224.23736891366676; w = -0.1271288806541458; τ_liq = 8.70519108187029e-12; τ_ice = 2.0788375361402434e-12; δ_0 = -2.3831539958081173e-5; δ_0i = 0.0; q = Thermodynamics.PhasePartition{Float64}(0.00046794210973515136, 0.0, 0.0004296824485967858); q_eq = Thermodynamics.PhasePartition{Float64}(0.00046794210973515136, 6.216785765933966e-5, 3.833631770125848e-5); Δt = 9.999999999480808; ts = Thermodynamics.PhaseNonEquil{Float64}(-35197.11200601141, 1.1297764844315694, Thermodynamics.PhasePartition{Float64}(0.00046794210973515136, 0.0, 0.0004340611316873781)); dTdt = 0.00045330389447046477; dqvdt = 8.245975226821409e-10;
q = Thermodynamics.PhasePartition{Float64}(0.001976663686494517, 0.0, 2.5519121725500658e-20); T = 262.691521171828; p = 84043.30717958807; ρ = 1.1133964727328423; w = -0.0005480114811913959; relaxation_timescale = TurbulenceConvection.GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale{Float64}(5048.808821939158, 0.1821615879079399, 7.334264855778188, 4.0874493700288115, 0.4645823732473287, 0.10566549047690477, -7.215527259640249, 0.0030042933101672126, -0.8520013654401164, true, true, TurbulenceConvection.RelaxationTimescaleArgs{Float64}(1.0e-6, 1.0e-6, 4.503599627370496e15, 4.503599627370496e15, 2.220446049250313e-16, 2.220446049250313e-16, 4.503599627370496e15, 4.503599627370496e15))

# q = TD.PhasePartition(q_tot, q_liq, q_ice)
q_tot, q_liq, q_ice = q.tot, q.liq, q.ice
q_vap = TD.vapor_specific_humidity(q)
q_vap_0 = q_vap

# area = FT(1.0)
# ts = TD.PhaseNonEquil_pTq(thermo_params, p_c, T, q)
# ρ = TD.air_density(thermo_params, ts)
# p = p_c
# ρ = FT(ρ_c) # 
# q_eq = TD.PhasePartition(
#     q.tot,
#     TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
#     TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
# )
# q_eq = TD.PhasePartition(
#     q.tot,
#     TD.q_vap_saturation_generic(thermo_params, T, TD.air_density(thermo_params, ts), TD.Liquid()),
#     TD.q_vap_saturation_generic(thermo_params, T, TD.air_density(thermo_params, ts), TD.Ice()),
# )
RH_liq = q_vap / q_eq.liq
RH_ice = q_vap / q_eq.ice

# δ_0 = q_vap - q_eq.liq
# δ_0i = q_vap - q_eq.ice
      
if test_loop


    # ρ_c = FT(1.0195120197862853); p_c = FT(74651.29711651804)
    # nc_file = "/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/test/debug/Output.SOCRATES_RF09_obs_data.bHES/stats/Stats.SOCRATES_RF09_obs_data.nc"
    # nc_file = "/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/test/debug/Output.SOCRATES_RF09_obs_data.rFUT/stats/Stats.SOCRATES_RF09_obs_data.nc"
    # nc_file = "/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/test/debug/Output.SOCRATES_RF09_obs_data.o4u0/stats/Stats.SOCRATES_RF09_obs_data.nc"
    nc_file = "/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/test/debug/Output.SOCRATES_RF09_obs_data.He0F/stats/Stats.SOCRATES_RF09_obs_data.nc"
    z_ind = 11

    p_c = NCDatasets.NCDataset(nc_file, "r") do data; data.group["reference"]["p_c"][z_ind+1]; end
    ρ_c = NCDatasets.NCDataset(nc_file, "r") do data; data.group["reference"]["ρ_c"][z_ind+1]; end

    get_mm_2015_style_verification_data_from_nc_file(t_ind::Int, env_or_updraft::Symbol; nc_file::String = nc_file, thermo_params = thermo_params, param_set = param_set, p_c = p_c, ρ_c = ρ_c, z_ind=z_ind) = begin
        T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt = NCDatasets.NCDataset(nc_file, "r") do data
            T = data.group["profiles"][string(env_or_updraft) * "_temperature"][z_ind+1, t_ind+1]
            w = data.group["profiles"][string(env_or_updraft) * "_w"][z_ind+1, t_ind+1]
            area = data.group["profiles"][string(env_or_updraft) * "_area"][z_ind+1, t_ind+1]
            q_tot = data.group["profiles"][string(env_or_updraft) * "_qt"][z_ind+1, t_ind+1]
            q_liq = data.group["profiles"][string(env_or_updraft) * "_ql"][z_ind+1, t_ind+1]
            q_ice = data.group["profiles"][string(env_or_updraft) * "_qi"][z_ind+1, t_ind+1]
            dTdt = data.group["profiles"][string(env_or_updraft) * "_dTdt"][z_ind+1, t_ind+1]
            dqvdt = data.group["profiles"][string(env_or_updraft) * "_dqvdt"][z_ind+1, t_ind+1]
            τ_ice = data.group["profiles"][string(env_or_updraft) * "_τ_ice"][z_ind+1, t_ind+1]
            τ_liq = data.group["profiles"][string(env_or_updraft) * "_τ_liq"][z_ind+1, t_ind+1]
            Δt = data.group["timeseries"]["t"][t_ind+2] - data.group["timeseries"]["t"][t_ind+1]
            return (T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt)
        end


        q = TD.PhasePartition(q_tot, q_liq, q_ice)
        ts = TD.PhaseNonEquil_pTq(thermo_params, p_c, T, q)
        ρ = TD.air_density(thermo_params, ts)
        # ρ = FT(ρ_c) # 
        q_vap = TD.vapor_specific_humidity(q)
        q_vap_0 = q_vap
        q_eq = TD.PhasePartition(
            q.tot,
            TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
            TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
        )

        return (; T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt, q, ts, ρ, q_vap, q_vap_0, q_eq)

    end
    #
    S = OrderedCollections.OrderedDict{Symbol, OrderedCollections.OrderedDict}(
        # :env => OrderedCollections.OrderedDict(),
        :updraft => OrderedCollections.OrderedDict()
    )
    Δt_factor = 2
    # Δt_factor = (1/4) * 2

    # # ------ #
    add_to_existing_plot = false
    linthresh = FT(NaN)
    # plot_only = :all
    plot_only = :liq
    # plot_only = :ice


    t_inds = 3950:3951

    for (i, t_ind) in enumerate(t_inds)
        @info "t_ind = $t_ind, i = $i/$(length(t_inds))"
        for updraft_or_env in filter(x -> x in [:updraft, :env], keys(S))
            global S_ql, S_qi, T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt, q, ts, ρ, q_vap, q_vap_0, q_eq, p, δ_0, δ_0i

            (; T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt, q, ts, ρ, q_vap, q_vap_0, q_eq) = get_mm_2015_style_verification_data_from_nc_file(t_ind, updraft_or_env)
            # S_ql, S_qi = morrison_milbrandt_2015_style(param_set, area, ρ, p_c, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q, q_eq, Δt*Δt_factor, ts; use_fix = use_fix)
            S_ql, S_qi = morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, p_c, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q, q_eq, Δt*Δt_factor, ts; use_fix = true)
            S[updraft_or_env][t_ind] = OrderedCollections.OrderedDict(:area => area, :S_ql => S_ql, :S_qi => S_qi)
            println("------------------------------------------------")
            p = p_c
            δ_0 = q_vap - q_eq.liq
            δ_0i = q_vap - q_eq.ice
            (; g, L_i, L_l, c_p, e_sl, e_si, dqsl_dT, dqsi_dT, q_sl, q_si, q_liq, q_ice, T_freeze, δ_0, δ_0i, Γ_l, Γ_i, dqvdt, dTdt) = get_params_and_go_to_mixing_ratio_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; use_fix = true)
            A_c = A_c_func_EPA(τ_ice, Γ_l, q_sl, q_si, g, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ) # Eq C4
            τ = τ_func_EPA(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)
            t_hit_liq_sat = t_δ_hit_value(FT(0), δ_0, A_c, τ)

            A_c_big = A_c_func_EPA(big(τ_ice), big(Γ_l), big(q_sl), big(q_si), big(g), big(w), big(c_p), big(e_sl), big(L_i), big(dqsl_dT), big(dqvdt), big(dTdt), big(p), big(ρ)) # Eq C4
            τ_big = τ_func_EPA(big(τ_liq), big(τ_ice), big(L_i), big(c_p), big(dqsl_dT), big(Γ_i))
            t_hit_liq_sat_big = t_δ_hit_value(big(0.), big(δ_0), A_c_big, τ_big)
            dδ = dδ_func_EPA(A_c, τ, δ_0, t_hit_liq_sat)

            S_ql_big = S_ql_func_EPA( A_c_big, τ_big, big(τ_liq), big(δ_0), big(Δt*Δt_factor), big(Γ_l))
            S_ql2 = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, Δt*Δt_factor, Γ_l)
            min_t = isinf(t_hit_liq_sat) ? 7.8657389978492365 : t_hit_liq_sat
            S_ql3 = S_ql_func_EPA( A_c, τ, τ_liq, δ_0, min_t, Γ_l)
            S_qi3 = S_qi_func_EPA( A_c, τ, τ_ice, δ_0, min_t, Γ_i, q_sl, q_si)

            # dδ3 = S_ql3*min_t + S_qi3*min_t
            dδ_tot = dδ_func_EPA(A_c, τ, δ_0,  Δt*Δt_factor)
            dδ_tot_big = dδ_func_EPA(A_c_big, τ_big, big(δ_0), big(Δt*Δt_factor))
            @info "dδ = $dδ; δ_0 + dδ = $(δ_0 + dδ); big(δ_0) + big(dδ) = $(big(δ_0) + dδ_tot_big)"

            @info "S_ql = $S_ql; S_ql2 = $S_ql2; S_ql3 = $S_ql3; Δt = $(Δt*Δt_factor); min_t = $min_t; t_hit_liq_sat = $t_hit_liq_sat"
            @info "t_hit_liq_sat = $t_hit_liq_sat; t_hit_liq_sat_big = $t_hit_liq_sat_big; A_c = $A_c; A_c_big = $A_c_big; τ = $τ; τ_big = $τ_big; dδ=$dδ; δ_0 = $δ_0; δ_0i = $δ_0i; Γ_l = $Γ_l; Γ_i = $Γ_i; q_sl = $q_sl; q_si = $q_si; q_liq = $q_liq; q_ice = $q_ice; T_freeze = $T_freeze; dqvdt = $dqvdt; dTdt = $dTdt"
            @info "S_ql = $S_ql; S_ql_big = $S_ql_big"
            @info "dTdT = $dTdt; dqvdt = $dqvdt; dTdt * dqsl_dT = $(dTdt * dqsl_dT); dTdt * dqsi_dT = $(dTdt * dqsi_dT) || dqvdt + (dTdt * dqsl_dT) = $(dqvdt - (dTdt * dqsl_dT)); dqvdt + (dTdt * dqsi_dT) = $(dqvdt - (dTdt * dqsi_dT))"
            println("==================================================================================================================================================================================================================================== #")
            println("")
        end
    end


    S[:mean] = OrderedCollections.OrderedDict(
        t_ind => OrderedCollections.OrderedDict(
            :area => sum(S[r][t_ind][:area] for r in keys(S) if r != :mean),
            :S_ql => sum(S[r][t_ind][:S_ql] * S[r][t_ind][:area] for r in keys(S) if r != :mean) / sum(S[r][t_ind][:area] for r in keys(S) if r != :mean),
            :S_qi => sum(S[r][t_ind][:S_qi] * S[r][t_ind][:area] for r in keys(S) if r != :mean) / sum(S[r][t_ind][:area] for r in keys(S) if r != :mean)
            ) for t_ind in sort(collect(keys(S[first(keys(S))])))
    )
    updraft_or_envs = filter(x -> x in [:updraft, :env], keys(S))
else
    t_inds = 1:1
    updraft_or_envs = [Symbol()]
    q_tot = q.tot
    q_liq = q.liq
    q_ice = q.ice
    global add_to_existing_plot = false
    plot_only = :all

end


#### PLOT ####
plot_rates = true

for t_ind in t_inds
    for updraft_or_env in updraft_or_envs
        if test_loop
            (; T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt, q, ts, ρ, q_vap, q_vap_0, q_eq) = get_mm_2015_style_verification_data_from_nc_file(t_ind, updraft_or_env)
        else
            global T, w, area, q_tot, q_liq, q_ice, dTdt, dqvdt, τ_ice, τ_liq, Δt, q, ts, ρ, q_vap, q_vap_0, q_eq, q
        end


        # ==================================================================================================================================================================================================================================== #
        # ==================================================================================================================================================================================================================================== #
        # ==================================================================================================================================================================================================================================== #
        global add_to_existing_plot
        global linthresh


        function symlog(x, n = -3)
            result = zeros(size(x, 1))
            for i in eachindex(x)
                result[i] = sign(x[i]) * (log10(1 + abs(x[i]) / (10^n)))
            end
            result
        end

        function symlogformatter(x, n; ndigits=20)
            # if sign(x) == 0
            #     "10^$(n)"
            # else
                s = sign(x) == 1 ? "+" : "-"
                nexp = sign(x) * (abs(x) + n)
                if sign(x) == -1
                    nexp = -nexp
                end
                n_exp_print = round(nexp, digits=ndigits)
                # get minimal
                # s * "10^$(n_exp_print)"
                num_str = @sprintf("%.*f", ndigits, nexp)  # ✅ dynamic precision
                return s * "10^" * num_str

            # end
        end

        Δts = 10 .^ range(-1, stop = 2, length = 1000)
        # Δts = 10 .^ range(-7, stop = 0, length = 100); @warn("temp test")
        # Δts = 10 .^ range(-20, stop = -10, length = 100); @warn("temp test")

        # Δts = 10 .^ range(log10(5), stop = log10(10), length = 100); @warn("temp test")

        # Δts = 10 .^ range(-18, stop = -10, length = 1000)

        qls = zeros(length(Δts))
        qis = zeros(length(Δts))
        qls_exp = zeros(length(Δts))
        qis_exp = zeros(length(Δts))
        δ_ΔTs = zeros(length(Δts))

        @showprogress dt=1 desc="Computing..." for (i, Δt) in enumerate(Δts)
            # qls[i], qis[i], δ_ΔTs[i] = morrison_milbrandt_2015_style(param_set, area, ρ,p, T, w, τ_liq, τ_ice, q_vap_0, q, q_eq, Δt, ts; use_fix=use_fix)
            # qls_exp[i], qis_exp[i] = morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, q_eq, Δt,)

            println("============================================================================================================================================================== Δt = $Δt")

            qls[i], qis[i] = 
                # (q_vap_0 .* FT(NaN), q_vap_0 * FT(NaN))
                # morrison_milbrandt_2015_style( param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q, q_eq, Δt, ts; use_fix = use_fix) 
                calculate_timestep_limited_sources(TC.StandardSupersaturationMoistureSourcesLimiter(), param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap, dqvdt, dTdt, q, q_eq, Δt, ts)


            println("----------------------------------------------------------------------------------------------------------------------------------------------------- Δt = $Δt (exp only)")

            qls_exp[i], qis_exp[i] = 
                # (q_vap_0 .* FT(NaN), q_vap_0 * FT(NaN))
                morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, p, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q, q_eq, Δt, ts; use_fix = use_fix)

            # qls[i] = qls_exp[i]
            # qis[i] = qis_exp[i]

            # qls_exp[i] = qls[i]
            # qis_exp[i] = qis[i]




            # if plot_rates
            #     qls[i] /= Δts[i]; qis[i] /= Δts[i]; qls_exp[i] /= Δts[i]; qis_exp[i] /= Δts[i];
            # end

            if plot_only == :ice
                qls[i] = qis[i]; qls_exp[i] = qis_exp[i]
            elseif plot_only == :liq
                qis[i] = qls[i]; qis_exp[i] = qls_exp[i]
            end




        end
        # plot the results, x scale log
        # linthresh = FT(-5)

        @info "τ_liq = $τ_liq; τ_ice = $τ_ice; T = $T; w = $w; ρ = $ρ; area = $area; p = $p; q_vap_0 = $q_vap_0; dqvdt = $dqvdt; dTdt = $dTdt; q_tot = $q_tot; q_liq = $q_liq; q_ice = $q_ice; q_eq.liq = $(q_eq.liq); q_eq.ice = $(q_eq.ice)"


        if add_to_existing_plot && Base.@isdefined(linthresh) # isdefined(:linthresh)
            initial_plot_func = Plots.plot!
        else
            valid  = (abs.(qls) .> 0) .|| (abs.(qis) .> 0) .|| (abs.(qls_exp) .> 0) .|| (abs.(qis_exp) .> 0)
            linthresh = mean(log10.((abs.((qls .+ qis .+ qls_exp .+ qis_exp) .* Δts))[valid])) - 1
            initial_plot_func = Plots.plot
        end

        yformatter = x -> symlogformatter.(x, linthresh)
        yscale = :symlog
        # yscale = :linear
        #


        factors = (plot_rates ? 1 : Δts)

        if yscale === :symlog
            initial_plot_func(
                Δts,
                symlog(qls .* factors, linthresh),
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
            plot!(Δts, symlog(qis .* factors, linthresh), label = "S_qi", color = :blue, yformatter = yformatter, linewidth = 2)
            plot!(Δts, symlog((qls .+ qis) .* factors, linthresh), label = "sum", color = :red, yformatter = yformatter, linewidth = 1)
            plot!(
                Δts,
                symlog(qls_exp .* factors, linthresh),
                label = "S_ql_exp",
                color = :green,
                linestyle = :dash,
                yformatter = yformatter,
                linewidth = 2,
            )
            plot!(
                Δts,
                symlog(qis_exp .* factors, linthresh),
                label = "S_qi_exp",
                color = :blue,
                linestyle = :dash,
                yformatter = yformatter,
                linewidth = 2,
            )
            plot!(
                Δts,
                symlog((qls_exp .+ qis_exp) .* factors, linthresh),
                label = "sum_exp",
                color = :red,
                linestyle = :dash,
                yformatter = yformatter,
                linewidth = 1,
            )
            hline!([0], color = :black, linestyle = :dash, yformatter = yformatter, label = "")

            if !plot_rates
                hline!(
                    symlog([(q_vap_0 - q_eq.liq) ], linthresh),
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
            end
            # plot!(Δts, symlog(-δ_ΔTs .* Δts, linthresh), label = "δ_ΔT", color=:orange, linestyle=:dashdot, yformatter=yformatter)
            plot!(
                Δts,
                symlog((q_vap_0 .- q_eq.ice) ./ τ_ice .* (plot_rates ? fill(1, size(Δts)) : Δts), linthresh),
                label = "S_i base",
                color = :darkblue,
                linestyle = :solid,
                yformatter = yformatter,
            )
            plot!(
                Δts,
                symlog((q_vap_0 .- q_eq.liq) ./ τ_liq .* (plot_rates ? fill(1, size(Δts)) : Δts), linthresh),
                label = "S_l base",
                color = :darkgreen,
                linestyle = :solid,
                yformatter = yformatter,
            )

            plot!(
                Δts,
                symlog( fill(-q.liq, size(Δts)) ./ (plot_rates ? Δts : 1) , linthresh),
                label = "-ql/Δt base",
                color = :darkgreen,
                linestyle = :solid,
                yformatter = yformatter,
            )
            plot!(
                Δts,
                symlog( fill(-q.ice, size(Δts)) ./ (plot_rates ? Δts : 1) , linthresh),
                label = "-qi/Δt base",
                color = :darkblue,
                linestyle = :solid,
                yformatter = yformatter,
            )
        elseif yscale === :linear
            initial_plot_func(
                Δts,
                qls .* Δts,
                label = "S_ql",
                dpi = 600,
                size = (1000, 400),
                xaxis = :log,
                yaxis = :linear,
                legend = :outertopright,
                color = :green,
                minorgrid = true,
                linewidth = 2,
            )
            plot!(Δts, qis .* Δts, label = "S_qi", color = :blue, linewidth = 2)
            plot!(Δts, (qls .+ qis) .* Δts, label = "sum", color = :red, linewidth = 1)
            #
            plot!(Δts, qls_exp .* Δts, label = "ql_exp", color = :green, linestyle = :dash, linewidth = 2)
            plot!(Δts, qis_exp .* Δts, label = "qi_exp", color = :blue, linestyle = :dash, linewidth = 2)
            plot!(Δts, (qls_exp .+ qis_exp) .* Δts, label = "sum_exp", color = :red, linestyle = :dash, linewidth = 1)
            #
            hline!([0], color = :black, linestyle = :dash, label = "")
            hline!([q_vap_0 - q_eq.liq], color = :gray, linestyle = :dashdot, yformatter = yformatter, label = "δ_l", alpha = 0.5)
            hline!([q_vap_0 - q_eq.ice], color = :blue, linestyle = :dashdot, yformatter = yformatter, label = "δ_i", alpha = 0.5)
            # plot!(Δts, -δ_ΔTs .* Δts, label = "δ_ΔT", color=:orange, linestyle=:dashdot, yformatter=yformatter)
            plot!(Δts, (q_vap_0 .- q_eq.ice) ./ τ_ice .* Δts, label = "S_i base", color = :darkblue, linestyle = :solid)
            plot!(Δts, (q_vap_0 .- q_eq.liq) ./ τ_liq .* Δts, label = "S_l base", color = :darkgreen, linestyle = :solid)

            plot!(Δts, fill(-q.liq, size(Δts)), label = "-ql/Δt base", color = :darkgreen, linestyle = :solid)
            plot!(Δts, fill(-q.ice, size(Δts)), label = "-qi/Δt base", color = :darkblue, linestyle = :solid)
        else
            error("Unsupported yscale: $yscale")

        end

        factors = (plot_rates ? 1 : Δts)
        ymin_l = yscale === :symlog ? minimum(symlog(((qls .* factors)..., (qls_exp .* factors)...), linthresh)) : minimum(((qls .* factors)..., (qls_exp .* factors)...))
        ymin_i = yscale === :symlog ? minimum(symlog(((qis .* factors)..., (qis_exp .* factors)...), linthresh)) : minimum(((qis .* factors)..., (qis_exp .* factors)...))
        ymax_l = yscale === :symlog ? maximum(symlog(((qls .* factors)..., (qls_exp .* factors)...), linthresh)) : maximum(((qls .* factors)..., (qls_exp .* factors)...))
        ymax_i = yscale === :symlog ? maximum(symlog(((qis .* factors)..., (qis_exp .* factors)...), linthresh)) : maximum(((qis .* factors)..., (qis_exp .* factors)...))
        ymin = min(ymin_l, ymin_i)
        ymax = max(ymax_l, ymax_i)

        if yscale === :symlog
            dy = ymax - ymin
            factor = 0.2
            ymin -=  factor * dy
            ymax +=  factor * dy
            # ymin = ymin > 0 ? ymin * 0.9 : ymin * 1.1
            # ymax = ymax > 0 ? ymax * 1.1 : ymax * 0.9
        elseif yscale ==:linear
            dy = ymax - ymin
            ymin -=  0.1 * dy
            ymax +=  0.1 * dy
        end

        if add_to_existing_plot
            current_ylims = Plots.ylims()
            ylims!(min(current_ylims[1], ymin), max(current_ylims[2], ymax))
        else
            ylims!(ymin, ymax)
        end


        vline!([Δts[1]], color = :purple, linestyle = :dash, label = "")
        vline!([Δts[end]], color = :purple, linestyle = :dash, label = "")

        # hline!([-q_liq], color=:green, linestyle=:dashdot, label="-q_liq(t=0)")
        # hline!([-q_ice], color=:blue, linestyle=:dashdot, label="-q_ice(t=0)")


        global add_to_existing_plot = true

    end
end
# save
savefig(joinpath(thisdir, "test_morrison_milbrandt_2015_style.png"))
