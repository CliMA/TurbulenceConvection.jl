

# reload_SSCF = false
# if reload_SSCF | !isdefined(Main, :SOCRATESSingleColumnForcings)
#     using Pkg
#     using Revise
#     Pkg.activate(expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl/test"))
#     import CLIMAParameters as CP # use CLIMAParameters = "0.7, 0.8, 0.9, 0.10"
#     # import ClimaParams as CPP # would using this trouble w/ TC.jl? it's a different uuid technically..., use ClimaParams = "0.10"
#     import Thermodynamics as TD
#     import Thermodynamics.Parameters as TDP
#     # FT = Float64

#     toml_dict = CP.create_toml_dict(FT; dict_type = "alias") # CP 0.7 and below, Thermodynamics 0.11 and above
#     aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
#     param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
#     thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
#     Pkg.activate(expanduser("~/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl/"))
#     include(
#         "/home/jbenjami/Research_Schneider/CliMa/SOCRATESSingleColumnForcings.jl/src/SOCRATESSingleColumnForcings.jl",
#     )
# end

FT = Float64

using Pkg

reload_TC = false
if reload_TC | !isdefined(Main, :TurbulenceConvection)
    Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/integration_tests/"))
    # Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/"))
    tc_dir = expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl")
    Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/integration_tests/"))
    # include(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/src/TurbulenceConvection.jl")) 
    using ForwardDiff # for netcdf.io
    include(joinpath(tc_dir, "driver", "NetCDFIO.jl"))
    include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
    include(joinpath(tc_dir, "driver", "Cases.jl"))
    include(joinpath(tc_dir, "driver", "parameter_set.jl"))
    include(joinpath(tc_dir, "driver", "dycore.jl"))
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    namelist = NameList.default_namelist("SOCRATES_RF09_Obs_data")
    param_set = create_parameter_set(namelist, toml_dict, FT)
    thermo_params = TCP.thermodynamics_params(param_set)
end

include("./morrison_milbrandt_2015_style.jl")




T = FT(260)
ρ = FT(0.5)

# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) + .5 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) - 0.5 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()) - 5.0 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))


# q_vap_0 = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) + 0.05 .* ( TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) - TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice()))
q_vap_0 =
    TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) +
    2.0 .* (
        TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid()) -
        TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())
    )


q_liq = FT(1e-5) * 10
q_ice = FT(1e-5) * 1
q = TD.PhasePartition(q_vap_0, q_liq, q_ice)

use_fix = true
w = FT(0)
τ_liq = FT(1.1)
τ_ice = FT(0.1)

τ_liq = FT(1e-20)
τ_ice = FT(1e-13)
# q_vap_0, q_liq, q_ice = (0.0027351234431657536, 4.4465037826352686e-61, 2.094458880577658e-19)
# q_vap_0, q_liq, q_ice = (0.0027351234431657536, 4.4465037826352686e-4, 2.094458880577658e-4)
# q = TD.PhasePartition(q_vap_0, q_liq, q_ice)


# τ_liq = FT(Inf)
# τ_ice = FT(Inf)

# τ_liq = FT(1000)
# τ_ice = FT(1000)

Δt = FT(0.0005)
area = FT(1.0)
# p = FT(750 .* 100)

ts = TD.PhaseNonEquil_ρTq(thermo_params, ρ, T, q)
p = TD.air_pressure(thermo_params, ts)
ρ = TD.air_density(thermo_params, ts)

q_vap = TD.vapor_specific_humidity(thermo_params, ts)
# @info(q_vap_0 - q_vap)
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

function symlogformatter(x, n)
    if sign(x) == 0
        "10^$(n)"
    else
        s = sign(x) == 1 ? "+" : "-"
        nexp = sign(x) * (abs(x) + n)
        if sign(x) == -1
            nexp = -nexp
        end
        s * "10^$(nexp)"
    end
end

Δts = 10 .^ range(-5, stop = 3, length = 100)
qls = zeros(length(Δts))
qis = zeros(length(Δts))
qls_exp = zeros(length(Δts))
qis_exp = zeros(length(Δts))
δ_ΔTs = zeros(length(Δts))
using Plots
for (i, Δt) in enumerate(Δts)
    # qls[i], qis[i], δ_ΔTs[i] = morrison_milbrandt_2015_style(param_set, area, ρ,p, T, w, τ_liq, τ_ice, q_vap_0, q, q_eq, Δt, ts; use_fix=use_fix)
    # qls_exp[i], qis_exp[i] = morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, q_eq, Δt,)

    qls[i], qis[i] = morrison_milbrandt_2015_style(
        param_set,
        area,
        ρ,
        p,
        T,
        w,
        τ_liq,
        τ_ice,
        q_vap_0,
        q,
        q_eq,
        Δt,
        ts;
        use_fix = use_fix,
    )
    qls_exp[i], qis_exp[i] =
        morrison_milbrandt_2015_style_exponential_part_only(param_set, area, ρ, T, w, τ_liq, τ_ice, q_vap_0, q_eq, Δt)
end
# plot the results, x scale log
linthresh = FT(-5)
linthresh = mean(log10.(abs.((qls .+ qis) .* Δts))) - 1
yformatter = x -> symlogformatter.(x, linthresh)
yscale = :symlog
# yscale = :linear
if yscale == :symlog
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
    )
    plot!(Δts, symlog(qis .* Δts, linthresh), label = "S_qi", color = :blue, yformatter = yformatter, linewidth = 2)
    plot!(Δts, symlog((qls .+ qis) .* Δts, linthresh), label = "sum", color = :red, yformatter = yformatter)
    plot!(
        Δts,
        symlog(qls_exp .* Δts, linthresh),
        label = "S_ql_exp",
        color = :green,
        linestyle = :dash,
        yformatter = yformatter,
    )
    plot!(
        Δts,
        symlog(qis_exp .* Δts, linthresh),
        label = "S_qi_exp",
        color = :blue,
        linestyle = :dash,
        yformatter = yformatter,
    )
    plot!(
        Δts,
        symlog((qls_exp .+ qis_exp) .* Δts, linthresh),
        label = "sum_exp",
        color = :red,
        linestyle = :dash,
        yformatter = yformatter,
    )
    hline!([0], color = :black, linestyle = :dash, yformatter = yformatter, label = "")
    hline!(
        symlog([q_vap_0 - q_eq.liq], linthresh),
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
elseif yscale == :linear
    plot(
        Δts,
        ql .* Δts,
        label = "ql",
        dpi = 600,
        size = (1000, 400),
        xaxis = :log,
        yaxis = :linear,
        legend = :outertopright,
        color = :green,
    )
    plot!(Δts, qis .* Δts, label = "qi", color = :blue)
    plot!(Δts, (ql .+ qi) .* Δts, label = "sum", color = :red)
    plot!(Δts, qls_exp .* Δts, label = "ql_exp", color = :green, linestyle = :dash)
    plot!(Δts, qis_exp .* Δts, label = "qi_exp", color = :blue, linestyle = :dash)
    plot!(Δts, (qls_exp .+ qis_exp) .* Δts, label = "sum_exp", color = :red, linestyle = :dash)
    hline!([0], color = :black, linestyle = :dash, label = "")
    hline!([q_vap_0 - q_eq.liq], color = :gray, linestyle = :dashdot, yformatter = yformatter, label = "δ_l")

end

ymin_l = yscale == :symlog ? minimum(symlog(qls .* Δts, linthresh)) : minimum(qls .* Δts)
ymin_i = yscale == :symlog ? minimum(symlog(qis .* Δts, linthresh)) : minimum(qis .* Δts)
ymax_l = yscale == :symlog ? maximum(symlog(qls .* Δts, linthresh)) : maximum(qls .* Δts)
ymax_i = yscale == :symlog ? maximum(symlog(qis .* Δts, linthresh)) : maximum(qis .* Δts)
ymin = min(ymin_l, ymin_i)
ymax = max(ymax_l, ymax_i)
ymin = ymin > 0 ? ymin * 0.9 : ymin * 1.1
ymax = ymax > 0 ? ymax * 1.1 : ymax * 0.9
ylims!(ymin, ymax)

# hline!([-q_liq], color=:green, linestyle=:dashdot, label="-q_liq(t=0)")
# hline!([-q_ice], color=:blue, linestyle=:dashdot, label="-q_ice(t=0)")

# save
thisfile = "/home/jbenjami/Research_Schneider/CliMa/TurbulenceConvection.jl/src/closures/test_mm_2015_style.jl"
thisdir = dirname(thisfile)
savefig(joinpath(thisdir, "morrison_milbrandt_2015_style_test.png"))
