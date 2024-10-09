# test vertical velocity


using Revise
using Pkg
Pkg.activate(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/"))

includet(expanduser("~/Research_Schneider/CliMa/TurbulenceConvection.jl/src/TurbulenceConvection.jl"))


import CLIMAParameters as CP
# # import Thermodynamics as TD
TC = TurbulenceConvection # set if loaded driver/main files
CM = TC.CM
TD = TC.TD
TDP = TC.TD.Parameters

FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)

aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
# append :thermo_params => thermo_params to the param_pairs
param_pairs = vcat(param_pairs, :thermo_params => thermo_params)
TP = typeof(thermo_params)
microphys_params = CM.Parameters.CloudMicrophysicsParameters{FT, TP}(; param_pairs..., thermo_params)

ρ_i = FT(916.7) # kg/m^3
r_is = FT(62.5e-6)
r_0 = FT(1e-6)

TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), 1.5, 1e-7)
TC.my_terminal_velocity(microphys_params, CM.CommonTypes.SnowType(), CM.CommonTypes.Chen2022Type(), 1.5, 1e-7)

# TC.my_terminal_velocity(microphys_params, CM.CommonTypes.SnowType(), CM.CommonTypes.Blk1MVelType(), 1.5, 1e-7)
# TC.my_terminal_velocity(microphys_params, CM.CommonTypes.RainType(), CM.CommonTypes.Chen2022Type(), 1.5, 1e-7)


# q = collect(0:.01:5) * 1e-3
q = range(0, stop = 5e-4, length = 5000)

D_trans = 0.625e-3 # the default
# ρ_air = 1.2 # the documentation plots
ρ_air = 0.8 # more realistic
v_rain = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.RainType(), CM.CommonTypes.Chen2022Type(), ρ_air, q)
v_ice =
    TC.my_terminal_velocity.(
        microphys_params,
        CM.CommonTypes.IceType(),
        CM.CommonTypes.Chen2022Type(),
        ρ_air,
        q;
        D_transition = FT(Inf),
    )
v_ice_large =
    TC.my_terminal_velocity.(
        microphys_params,
        CM.CommonTypes.IceType(),
        CM.CommonTypes.Chen2022Type(),
        ρ_air,
        q;
        D_transition = FT(0),
    )
v_snow =
    TC.my_terminal_velocity.(
        microphys_params,
        CM.CommonTypes.SnowType(),
        CM.CommonTypes.Chen2022Type(),
        ρ_air,
        q;
        D_transition = FT(Inf),
    )


v_ice_dmax =
    TC.my_terminal_velocity.(
        microphys_params,
        CM.CommonTypes.IceType(),
        CM.CommonTypes.Chen2022Type(),
        ρ_air,
        q;
        D_transition = D_trans,
        Dmax = FT(62.5e-6),
    )
v_ice_dmax2 =
    TC.my_terminal_velocity.(
        microphys_params,
        CM.CommonTypes.IceType(),
        CM.CommonTypes.Chen2022Type(),
        ρ_air,
        q;
        D_transition = D_trans,
        Dmax = FT(0.001),
    )

v_ice_dmax3 =
    (
        (q, N) -> TC.my_terminal_velocity(
            microphys_params,
            CM.CommonTypes.IceType(),
            CM.CommonTypes.Chen2022Type(),
            ρ_air,
            q;
            D_transition = D_trans,
            Dmax = FT(62.5e-6 * 2),
            Nt = FT(N),
        )
    ).(q, q ./ (4 / 3 * ρ_i * π * r_is^3))
v_ice_dmax4 =
    (
        (q, N) -> TC.my_terminal_velocity(
            microphys_params,
            CM.CommonTypes.IceType(),
            CM.CommonTypes.Chen2022Type(),
            ρ_air,
            q;
            D_transition = D_trans,
            Dmax = FT(62.5e-6 * 2),
            Nt = FT(N),
        )
    ).(q, q ./ (4 / 3 * ρ_i * π * r_0^3))

# make plot
plot(
    q,
    v_rain,
    label = "Rain",
    xlabel = "q [kg/kg]",
    ylabel = "v [m/s]",
    title = "Terminal velocity vs q",
    color = :blue,
    linewidth = 2,
    dpi = 600,
    xaxis = :log10,
    yaxis = :symlog,
)
plot!(q, v_snow, label = "Snow", color = :red, linewidth = 2)
plot!(q, v_ice, label = "Ice", color = :pink, linewidth = 2)
plot!(q, v_ice_large, label = "large ice", color = :pink, linestyle = :dot, linewidth = 2)


plot!(q, v_ice_dmax, label = "Ice Dmax", color = :purple, linewidth = 2)
plot!(q, v_ice_dmax2, label = "Ice Dmax2", color = :purple, linestyle = :dot, linewidth = 2)
plot!(q, v_ice_dmax3, label = "Ice Dmax3 - low N", color = :purple, linestyle = :dash, linewidth = 2)
plot!(q, v_ice_dmax4, label = "Ice Dmax4 - high N", color = :purple, linestyle = :dashdot, linewidth = 2)



# set ylim
ylims!(-0.5, 5)
# ylims!(.00005, 5)

# set xlim
xlims!(q[2], q[end])

# save
thisdir = dirname(@__FILE__)
@info(thisdir)
savefig(joinpath(thisdir, "terminal_velocity_vs_q.png"))
