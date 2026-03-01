

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

const FT = Float64

using Pkg
using Statistics: mean
using Printf: @sprintf

thisdir = dirname(@__FILE__)
@info thisdir

reload_TC = false
if reload_TC ||
   !isdefined(Main, :TurbulenceConvection) ||
   !isdefined(Main, :param_set) ||
   !isdefined(Main, :thermo_params)
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
    namelist["meta"] = Dict("simname" => "test_standard_supersat", "uuid" => "")

    param_set = create_parameter_set(namelist, toml_dict, FT) # this creates an override file in  a directory we don' need...
    thermo_params = TCP.thermodynamics_params(param_set)
end

call_from_TC = false
if call_from_TC
    # morrison_milbrandt_2015_style = TC.morrison_milbrandt_2015_style
    # morrison_milbrandt_2015_style_exponential_part_only = TC.morrison_milbrandt_2015_style_exponential_part_only
    calculate_timestep_limited_sources = TC.calculate_timestep_limited_sources
else
    # Things to load if we're calling directly form the script rather than TC methods
    # need these types defined before including
    const NoMoistureSourcesLimiter = TC.NoMoistureSourcesLimiter
    const BasicMoistureSourcesLimiter = TC.BasicMoistureSourcesLimiter
    const StandardSupersaturationMoistureSourcesLimiter = TC.StandardSupersaturationMoistureSourcesLimiter
    const MorrisonMilbrandt2015MoistureSourcesLimiter = TC.MorrisonMilbrandt2015MoistureSourcesLimiter
    const MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter =
        TC.MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter
    include("../microphysics_coupling_limiters.jl")
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
    ) .* 1.005

q_vap_0 *= 1


q_liq = FT(1e-8) * 1 * 25
q_ice = FT(1e-8) * 1 * 10 * 1000
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_liq = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
# q_ice = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
q = TD.PhasePartition(q_vap_0 + q_liq + q_ice, q_liq, q_ice)

use_fix = false
# w = FT(+1000000.5)
w = FT(+0.0)

τ_liq = FT(1e-2)
τ_ice = FT(.935e-4)
# τ_ice = FT(.95e-4)

# τ_ice = FT(1e-2)
# τ_liq = FT(.935e-4)

τ_ice = FT(10.0)
τ_liq = FT(10.0)

# τ_ice = FT(1.0)
# τ_liq = FT(100.0)


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

ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
p = TD.air_pressure(thermo_params, ts)
ρ = TD.air_density(thermo_params, ts)

q_vap = TD.vapor_specific_humidity(thermo_params, ts)
q_eq = TD.PhasePartition(
    q.tot,
    TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()),
    TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()),
)



# Dict(
#     "Normal" => morrison_milbrandt_2015_style(param_set, area, ρ,p, T, w, τ_liq, τ_ice, q_vap_0, q, q_eq, Δt, ts),
#     "Exp Only" => morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, q_eq, Δt,)
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
    # s * "10^$(n_exp_print)" # this doenst work bc stuff like BigFloats dont use minimal printing and wont be truncated
    num_str = @sprintf("%.*f", ndigits, nexp)  # ✅ dynamic precision
    return s * "10^" * num_str
    # end
end

# Δts = 10 .^ range(-7, stop = 5, length = 100)
Δts = 10 .^ range(-4, stop = 4, length = 100);
@warn("temp test");
# Δts = 10 .^ range(log10(5), stop = log10(10), length = 100); @warn("temp test")

# Δts = 10 .^ range(-18, stop = -10, length = 1000)

qls = zeros(length(Δts))
qis = zeros(length(Δts))
qls_exp = zeros(length(Δts))
qis_exp = zeros(length(Δts))
δ_ΔTs = zeros(length(Δts))
using Plots
ENV["GKSwstype"] = "nul"
using ProgressMeter
@showprogress dt = 1 desc = "Computing..." for (i, Δt) in enumerate(Δts)
    # qls[i], qis[i], δ_ΔTs[i] = morrison_milbrandt_2015_style(param_set, area, ρ,p, T, w, τ_liq, τ_ice, q_vap_0, q, q_eq, Δt, ts; opts = MM2015Opts{Float64}(use_fix=use_fix))
    # qls_exp[i], qis_exp[i] = morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, dqvdt, dTdt, q_eq, Δt,)

    # println("============================================================================================================================================================== Δt = $Δt")

    qls[i], qis[i] =
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




end
# plot the results, x scale log
linthresh = FT(-5)
valid = abs.(qls .+ qis) .> 0
linthresh = mean(log10.((abs.((qls .+ qis) .* Δts))[valid])) - 1
yformatter = x -> symlogformatter.(x, linthresh)
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
        # left_margin=40Plots.PlotMeasures.mm, # can't tell when it seems needed lol
    )
    plot!(Δts, symlog(qis .* Δts, linthresh), label = "S_qi", color = :blue, yformatter = yformatter, linewidth = 2)
    plot!(Δts, symlog((qls .+ qis) .* Δts, linthresh), label = "sum", color = :red, yformatter = yformatter)
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
    hline!([0], color = :black, linestyle = :dash, label = "")
    hline!([q_vap_0 - q_eq.liq], color = :gray, linestyle = :dashdot, yformatter = yformatter, label = "δ_l")
    plot!(Δts, (q_vap_0 .- q_eq.ice) ./ τ_ice .* Δts, label = "S_i base", color = :pink, linestyle = :solid)


end

ymin_l = yscale === :symlog ? minimum(symlog(qls .* Δts, linthresh)) : minimum(qls .* Δts)
ymin_i = yscale === :symlog ? minimum(symlog(qis .* Δts, linthresh)) : minimum(qis .* Δts)
ymax_l = yscale === :symlog ? maximum(symlog(qls .* Δts, linthresh)) : maximum(qls .* Δts)
ymax_i = yscale === :symlog ? maximum(symlog(qis .* Δts, linthresh)) : maximum(qis .* Δts)
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

# save
savefig(joinpath(thisdir, "test_standard_supersat.png"))
