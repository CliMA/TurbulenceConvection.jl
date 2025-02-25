# test vertical velocity


using Revise
using Pkg
using Plots
using Statistics: mean
Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))

includet(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/src/TurbulenceConvection.jl"))


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
# ρ_air = 1.2 # the documentation plots
ρ_air = 0.8 # more realistic
# ρ_air = ρ_i


r_is = 0.00011799138385768925

q_is = 4/3 * ρ_i *  π * r_is^3
q_i0 = 4/3 * ρ_i *  π * r_0^3




# TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), 1.5, 1e-7)
# TC.my_terminal_velocity(microphys_params, CM.CommonTypes.SnowType(), CM.CommonTypes.Chen2022Type(), 1.5, 1e-7)

# TC.my_terminal_velocity(microphys_params, CM.CommonTypes.SnowType(), CM.CommonTypes.Blk1MVelType(), 1.5, 1e-7)
# TC.my_terminal_velocity(microphys_params, CM.CommonTypes.RainType(), CM.CommonTypes.Chen2022Type(), 1.5, 1e-7)


CloudMicrophysics = TC.CM
Thermodynamics = TC.TD
microphys_params = CloudMicrophysics.Parameters.CloudMicrophysicsParameters{Float64, Thermodynamics.Parameters.ThermodynamicsParameters{Float64}}(0.024, 2.26e-5, 1000.0, 0.072, 10.0, 10.0, 1000.0, 5.0e-6, 0.02, 0.55, 1.6e-5, 0.00011799138385768925, 2.0e7, 1.0e-5, 3.0, 1.0, 0.0, 1.5, 0.53, 1.6e7, 0.001, 3.0, 2.0, 0.5, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 43.58563716234225, 6.149368373672493e-7, 8.210052874283406e-8, 151437.07675216204, 0.65, 0.44, 0.63, 4.36e9, 0.001, 2.0, 2.0, 0.25, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.8, 0.1, 1.0, 0.1, 1.0, 916.7, Thermodynamics.Parameters.ThermodynamicsParameters{Float64}(273.16, 100000.0, 100000.0, 1859.0, 4181.0, 2100.0, 2.5008e6, 2.8344e6, 611.657, 273.16, 273.15, 150.0, 1000.0, 298.15, 6864.8, 10513.6, 0.28571428571, 8.3144598, 0.02897, 0.01801528, 290.0, 220.0, 9.81, 233.0, 1.0e7), 3268.0, 2.3333333333333335, -0.3333333333333333, 7.0e-6, 4.7, 3.0, 3.0e34, -1.7, 4.7, -3.3, 3.9, 9.9, 2.0e8, 6.0, 7.42e13, 2.47, -1.79, -1.47, 67.0, 1.15, -1.3, 7.5, 1.08e10, 1.0, 2.0, 4.44e9, 5.25, 7.12, 60.7, 2.6e-10, 5.0e-6, 2.0, 1.225, 400.0, 0.7, 3.0, 5.0e-5, 4.0, -5.0, 0.0009, 0.00035, 1000.0, 2300.0, 9.65, 10.3, 600.0, 0.78, 0.308, 159.0, 0.266, 250000.0, 2.0e7, 1000.0, 10000.0, 0.115231, 0.044612, -0.263166, 4.7178, -0.47335, 2.2955, 2.2955, 1.1451, 0.038465, 0.0, 0.184325, 0.184325, 0.00174079, 0.0378769, 0.263503, 0.575231, 0.0909307, 0.515579, -0.345387, 0.177362, -0.000427794, 0.00419647, -0.156593, 0.0189334, 0.1377817, -3.35641, 0.0156199, 0.765337, -0.0309715, 1.55054, 0.518349, 1.35, 220.0, 1.03, 1.07, 4.7, 9.2, 1.17, 1.03, 0.43, 2.35, 22.62, -1.35, 54.58834, -10.54758, 54.48075, -10.66873, 0.26, 0.34, -906.7, 8502.0, 26924.0, 29180.0, 235.0, 185.0, 1.4408, 23.306, 5.3465, 12.0, 8.19, -5814.0, 928.9, 1876.7, 0.058443, 2170.0, 0.9, 2.0, 1.0, 1.12, 0.132, 1770.0, 1.0, 3.0, 1.0, 0.53, 7.38e-11, 1.9, 1.88, 0.2285, 3.95451, 3.373738, 9.702973, -11.48166, 12.62259, 25.49469, -0.007066146, 0.1810722, 2.891024, 3.138719, 182.4495, -23.8002, 1.203451, 37.03029, -4.188065, 0.227413, 8.003471, 3.071246, 1.5703478e-6, 0.0048314, 0.0400097, 1.84826, 0.00136641, 1.56588, 0.186303, 0.029, 0.012, 8.05e-16, 1.2e-11, -640.0, 440.0, 3.27e-21)




# q = collect(0:.01:5) * 1e-3
# q = range(0, stop = 5e-4, length = 5000)
# q = 10 .^ range(-300, stop = -3, length = 5000)
q = 10 .^ range(-8, stop = -3, length = 5000)



# N0 = 1.1784787739629964e-150; q0 = 7.43342e-159; w0 = 0.868177
# N0 = 1.3538690695057672e-134; q0 = 8.53972e-143; w0 = 1283.26
N0 = 7.378643617682342e-131; q0 = 5.00278e-127; w0 = 6005.38
# N0 = Nt; q = qt;  q_is = q_ist; w0 = wt
# valid = q .> 0
# q = q[valid]; N0 = N0[valid]; w0 = w0[valid]

D_trans = 0.625e-3 # the default
# D_trans = Inf # test

v_rain      = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.RainType(), CM.CommonTypes.Chen2022Type(), ρ_air, q)
v_ice       = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.IceType() , CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=FT(Inf))
v_ice_large = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.IceType() , CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=FT(0))
v_snow      = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.SnowType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=FT(Inf))


v_ice_dmax       = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.IceType() , CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(62.5e-6))
v_ice_dmax2      = TC.my_terminal_velocity.(microphys_params, CM.CommonTypes.IceType() , CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(.001))

v_ice_dmax3 = ((q, N) -> TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(62.5e-6 * 2), Nt=FT(N))).(q, q./(q_is) )
# v_ice_dmax4 = ((q, N) -> TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(62.5e-6 * 2), Nt=FT(N))).(q, q./q_i0 )
v_ice_dmax4 = ((q, N) -> TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(Inf), Nt=FT(N))).(q, q./q_i0 )
# v_ice_dmax4 = ((q, N) -> TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(Inf), Nt=FT(N))).(q, N0 )

v_ice_dmax5 = ((q, N) -> TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(Inf), Nt=FT(N))).(q, max.(N0, q/q_is) )

N = NaN; v_ice_dmax6 = ((q, N) -> TC.my_terminal_velocity(microphys_params, CM.CommonTypes.IceType(), CM.CommonTypes.Chen2022Type(), ρ_air, q; D_transition=D_trans, Dmax=FT(Inf), Nt=FT(N))).(q, N )


# -------------------------------------------------------------------- #

function symlog(x; n = -3)
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


# linthresh = FT(-100)
# linthreshy = mean(log10.(abs.(v_ice_dmax[2500:end]))) - 1
linthreshy = -3.
yformatter = x -> symlogformatter.(x, linthreshy)
yscale = :symlog

use_symlog = true
if use_symlog
    scale_func = x -> symlog(x; n=linthreshy)
else
    scale_func = x -> x
end

# -------------------------------------------------------------------- #

# make plot
plot()
# plot!(q, scale_func(v_rain), label = "Rain", xlabel = "q [kg/kg]", ylabel = "v [m/s]", title = "Terminal velocity vs q", color = :blue, linewidth=2, dpi=1200, xaxis=:log10, yaxis=:symlog)
# plot!(q, scale_func(v_snow), label = "Snow", color = :red, linewidth=2)
# plot!(q, scale_func(v_ice), label = "Ice", color = :pink, linewidth=2)
# plot!(q, scale_func(v_ice_large), label = "large ice", color=:pink, linestyle=:dot, linewidth=2)


plot!(q, scale_func(v_ice_dmax), label = "Ice Dmax", color = :purple, linewidth=2)
plot!(q, scale_func(v_ice_dmax2), label = "Ice Dmax2", color = :purple, linestyle=:dot, linewidth = 2)
plot!(q, scale_func(v_ice_dmax3), label = "Ice Dmax3 - low N", color = :purple, linestyle=:dash, linewidth = 2)
plot!(q, scale_func(v_ice_dmax4), label = "Ice Dmax4 -  tiny N", color = :green, linestyle=:dashdot, linewidth = 2)
plot!(q, scale_func(v_ice_dmax5), label = "Ice Dmax5 - truncate tiny N", color = :lime, linestyle=:dashdot, linewidth = 2)
plot!(q, scale_func(v_ice_dmax6), label = "Ice Dmax6 - assume N", color = :orange, linestyle=:dashdot, linewidth = 2)

plot!([q[1], q[end]], scale_func([0, 0]), color = :black, linestyle = :dash, linewidth = 1, dpi=1200,  label=nothing)
plot!([q[1], q[end]], scale_func(+10 .^ [linthreshy, linthreshy]), color = :gray, linestyle = :dash, linewidth = 1, dpi=1200,  label=nothing)
plot!([q[1], q[end]], scale_func(-10 .^ [linthreshy, linthreshy]), color = :gray, linestyle = :dash, linewidth = 1, dpi=1200,  label=nothing)
plot!([q[1], q[end]], scale_func([1, 1]), color = :black, linestyle = :dot, linewidth = 1, dpi=1200,  label=nothing)


plot!([q0, q0], scale_func([- 10 .^ 1, 10 .^ 10]), color = :red, linestyle = :dot, linewidth = 1, dpi=1200,  label=nothing)

plot!([q[1], q[end]], scale_func([w0, w0]), color = :red, linestyle = :dot, linewidth = 1, dpi=1200,  label=nothing)
# plot!(q, scale_func(w0), color = :pink, linestyle = :dot, linewidth = 1, dpi=1200,  label="truth")


plot!(yformatter = yformatter, xaxis=:log,  legend=:outertopright, minorgrid=false, linewidth=2)

# set ylim
# ylims!(-.5, 5)
# ylims!(.00005, 5)

# set xlim
xlims!(q[2], q[end])



# save
thisdir = dirname(@__FILE__)
@info(thisdir)
savefig(joinpath(thisdir, "terminal_velocity_vs_q.png"))


