using NCDatasets
using CairoMakie
using Glob

case_names = ["DYCOMS_RF02","Rico","TRMM_LBA"]
rates = ["autocon", "accretion", "depsub", "evap", "melt"]
w_paths = ["lwp_mean", "iwp_mean", "rwp_mean", "swp_mean"]
clouds = ["cloud_cover", "cloud_base", "cloud_top"]
updrafts_qw = ["updraft_qt", "updraft_ql", "updraft_qi", "updraft_w"]
env_qw = ["env_qt", "env_ql", "env_qi"]
q_water = ["qt_mean","ql_mean","qi_mean","qr_mean","qs_mean"]
w_flux = ["rain_flux", "snow_flux"]


dir = "C:\\Users\\khanh\\SURF\\TurbulenceConvection.jl\\"

# plotting function for wp
function plot_wp(path::String, plot::Axis, ds::Dataset)
    t = ds.group["timeseries"]["t"]/3600
    p = ds.group["timeseries"][path]
    l =lines!(plot, t, p)
    return l
end

# plotting function for clouds
function plot_clouds(type::String, plot::Axis, ds::Dataset)
    t = ds.group["timeseries"]["t"]/3600
    p = ds.group["timeseries"][type]*1e-3
    l =lines!(plot, t, p)
    return l
end

# plotting heatmap for qw
function plot_qw_timeseries(name::String, fig::GridPosition, ds::Dataset, type::String)
    height = ds.group["profiles"]["zc"]*1e-3
    time = ds.group["timeseries"]["t"]/3.6e3
    qw = transpose(ds.group["profiles"][name]*1e3)
    ax, hm =heatmap(fig[1,1], time, height, qw)
    ax.xlabel = "hours(h)"
    ax.ylabel = "height(km)"
    ax.title = name*" time evolution "*type 
    Colorbar(fig[1,2], hm, label = "g/kg")
end


function plot_all_wp(case::String, type::String)
    fig = Figure(resolution =(1200,1200))
    plots = [Axis(fig[j,i], 
        xlabel = "hours(h)", 
        ylabel = "density(kg/m^2)",
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for i in 1:2, j in 1:2] 
    
    wp_lines = [[] for i in 1:4]
    rates = 0:0.5:2.0
    # deposition/sublimation rates of 2.0 occur on shorter time spans and thus are not considered
    if type == "depsub"
        rates = 0:0.5:1.5
    end

    for v in rates
        ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")
        for (i, path) in zip(1:4, w_paths)
            plots[i].title = case*" "*path*" variable "*type
            push!(wp_lines[i], plot_wp(path, plots[i], ds))
        end
    end

    for i in 1:4
        axislegend(plots[i], wp_lines[i], [string(j) for j in rates], type*" rates", position = :lt)
    end
    save(dir*case*"\\wp_plots\\"*case*"_wp_"*type*".png", fig)

    return fig
end

function plot_all_clouds(case::String, type::String)
    fig = Figure(resolution =(1800,600))
    plots = [Axis(fig[j,i], 
        xlabel = "hours(h)", 
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for i in 1:3, j in 1:1] 
    
    c_lines = [[] for i in 1:3]
    rates = 0:0.5:2.0
    # deposition/sublimation rates of 2.0 occur on shorter time spans and thus are not considered
    if type == "depsub"
        rates = 0:0.5:1.5
    end

    for v in rates
        ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")
        for (i, cloud) in zip(1:4, clouds)
            plots[i].title = case*" "*cloud*" variable "*type
            if cloud == "cloud_cover"
                plots[i].ylabel = "coverage"
                push!(c_lines[i], plot_wp(cloud*"_mean", plots[i], ds))
            else
                plots[i].ylabel = "height(km)"
                push!(c_lines[i], plot_clouds(cloud*"_mean", plots[i], ds))
            end
        end
    end

    for i in 1:3
        axislegend(plots[i], c_lines[i], [string(j) for j in rates], type*" rates", position = :lt)
    end
    save(dir*case*"\\cloud_plots\\"*case*"_clouds_"*type*".png", fig)

    return fig
end

function plot_qw(case::String, type::String)
    # initialize
    print(type*"\n")
    figs = [Figure(resolution=(1400, 1600)) for i in 1:4]
    base = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_1.0.nc")
    mm = size(base.group["profiles"]["t"])[1]

    # plotting normal rates
    plots = [[figs[i][j,k] for  k in 1:2, j in 2:3] for i in 1:4]
    for (qw,i) in zip(["ql","qi", "qr", "qs"],1:4) 
        plot_qw_timeseries(qw*"_mean", figs[i][1,1], base, "1.0")
    end

    values = ["0.0", "0.5", "1.5", "2.0"]
    count = 4
    if type == "depsub"
        values = ["0.0", "0.5", "1.5"]
        count = 3
    end
    #getting differences
    diffs = [[] for i in 1:4]
    for (v, j) in zip(values, 1:count)
        print(v*"\n")
        ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*v*".nc")
        for (qw,i) in zip(["ql","qi", "qr", "qs"], 1:4)
            height = ds.group["profiles"]["zc"]*1e-3
            time = ds.group["profiles"]["t"]*3600
            q = transpose(ds.group["profiles"][qw*"_mean"][:,1:mm]-base.group["profiles"][qw*"_mean"])*1e3
            push!(diffs[i], (time, height, q))
        end
    end

    #plotting differences
    for (qw,i) in zip(["ql","qi", "qr", "qs"], 1:4)
        print(qw*"\n")
        extremas = map(extrema, [diffs[3][i][3] for i in 1:count])
        global_min = minimum(t->first(t), extremas)
        global_max = maximum(t->last(t), extremas)
        m = max(abs(global_min), abs(global_max))
        range = (-m, m)
        for (v, j) in zip(values, 1:count)
            ax, hm = heatmap(plots[i][j], diffs[i][j][1][1:mm], diffs[i][j][2], diffs[i][j][3], colormap = :bwr, colorrange = range)
            ax.xlabel = "hours(h)"
            ax.ylabel = "height(km)"
            ax.title = qw*"_mean time evolution difference "*type*"="*v        
        end
        cb = Colorbar(figs[i][2:3, 3], label="g/kg", limits=range, colormap = :bwr)
        save(dir*case*"\\qw_plots\\"*case*"_"*type*"_"*qw*"_mean.png", figs[i]) 
    end
    return(figs)
end

function plot_wp_condensed(case::String)
    fig = Figure(resolution =(1200,1200))
    plots = [Axis(fig[i,j], 
        xlabel = "rate", 
        ylabel = "density(kg/m^2)",
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for j in 1:2, i in 1:2]
    c_lines = [[] for i in 1:4]
    for type in rates
        l = [[] for i in 1:4]
        values = 0.0:0.5:2.0
        if type == "depsub"
            values = 0.0:0.5:1.5
        end
        for v in values
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")
            t_max = size(ds.group["timeseries"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["timeseries"]["t"][t_max]/3600/6
            for (wp,i) in zip(w_paths, 1:4)
                push!(l[i], sum(ds.group["timeseries"][wp][t_min:t_max])/(t_max-t_min+1))
            end
        end
        for (wp,i) in zip(w_paths, 1:4)
            plots[i].title = wp*" last x hour(s)"
            r_l = lines!(plots[i], values, convert(Array{Float32}, l[i]))
            push!(c_lines[i], r_l)            
        end
    end 
    for i in 1:4
        axislegend(plots[i], c_lines[i], rates, "rates", position = :lt)
    end
    save(dir*case*"\\wp_condensed\\"*case*"_wp_condensed.png", fig)
    return fig
end

function plot_cloud_condensed(case::String)
    fig = Figure(resolution =(1800,600))
    plots = [Axis(fig[1,j], 
        xlabel = "rate", 
        ylabel = "",
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for j in 1:3]
    c_lines = [[] for i in 1:3]
    for type in rates
        l = [[] for i in 1:3]
        values = 0.0:0.5:2.0
        if type == "depsub"
            values = 0.0:0.5:1.5
        end
        for v in values
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")
            t_max = size(ds.group["timeseries"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["timeseries"]["t"][t_max]/3600/6
            for (cloud,i) in zip(clouds, 1:3)
                push!(l[i], sum(ds.group["timeseries"][cloud*"_mean"][t_min:t_max])/(t_max-t_min+1))
            end
        end
        for (cloud,i) in zip(clouds, 1:3)
            plots[i].title = cloud*" last x hour(s)"
            if cloud == "cloud_cover"
                plots[i].ylabel = "coverage"
                r_l = lines!(plots[i], values, convert(Array{Float32}, l[i]))
                push!(c_lines[i], r_l)
            else
                plots[i].ylabel = "height(km)"
                r_l = lines!(plots[i], values, convert(Array{Float32}, l[i]/1e3))
                push!(c_lines[i], r_l)
            end            
        end
    end 
    for i in 1:3
        axislegend(plots[i], c_lines[i], rates, "rates", position = :lb)
    end
    save(dir*case*"\\cloud_condensed\\"*case*"_cloud_condensed.png", fig)
    return fig
end

function plot_qw_condensed(case::String, type::String)
    fig = Figure(resolution =(1800,1200))
    values = 0.0:0.5:2.0
    if type == "depsub"
        values = 0.0:0.5:1.5
    end
    subplots = [Axis(fig[i,j], 
        xlabel = "g/kg", 
        ylabel = "height(km)",
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for j in 1:3, i in 1:2]
    qw_lines = [[] for i in 1:5]
    for r in values
        ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(r)*".nc")
        height = ds.group["profiles"]["zc"]*1e-3
        t_max = size(ds.group["profiles"]["t"])[1]
        t_min = floor(Int, 5/6 * t_max)
        period = ds.group["profiles"]["t"][t_max]/3600/6
        for (qw,i) in zip(q_water, 1:5)
            q=0
            for x in t_min:t_max
                q = ds.group["profiles"][qw][:, x]
            end
            q = q/(t_max-t_min+1)*1e3
            subplots[i].title = qw*" last "*string(period)*" hour(s)"
            l = lines!(subplots[i], q, height)
            push!(qw_lines[i], l)
        end
    end
    for i in 1:5
        axislegend(subplots[i], qw_lines[i], [string(j) for j in values], type, position = :rt)
    end
    
    save(dir*case*"\\qw_condensed\\"*case*"_qw_mean_condensed_"*type*".png", fig)    
    return fig
end


function plot_flux_condensed(case::String)
    fig = Figure(resolution =(1200,600))
    plots = [Axis(fig[1,j], 
        xlabel = "rate", 
        ylabel = "",
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for j in 1:2]
    c_lines = [[] for i in 1:2]
    for type in rates
        l = [[] for i in 1:2]
        values = 0.0:0.5:2.0
        if type == "depsub"
            values = 0.0:0.5:1.5
        end
        for v in values
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")
            t_max = size(ds.group["timeseries"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["timeseries"]["t"][t_max]/3600/6
            for (flux,i) in zip(w_flux, 1:2)
                push!(l[i], sum(ds.group["profiles"][flux][1,t_min:t_max])/(t_max-t_min+1))
            end
        end
        for (flux,i) in zip(w_flux, 1:2)
            plots[i].title = flux*" last x hour(s)"
            r_l = lines!(plots[i], values, convert(Array{Float32}, l[i]))
            push!(c_lines[i], r_l)
        end
    end 
    for i in 1:2
        axislegend(plots[i], c_lines[i], rates, "rates", position = :lt)
    end
    save(dir*case*"\\flux_condensed\\"*case*"_flux_condensed.png", fig)
    return fig
end

function plot_env_updraft_condensed(case::String, type::String)
    fig = Figure(resolution =(1800,1200))
    values = 0.0:0.5:2.0
    if type == "depsub"
        values = 0.0:0.5:1.5
    end
    subplots = [Axis(fig[i,j], 
        xlabel = "", 
        ylabel = "height(km)",
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for j in 1:3, i in 1:2]
    qw_lines = [[] for i in 1:6]
    for r in values
        ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(r)*".nc")
        height = ds.group["profiles"]["zc"]*1e-3
        t_max = size(ds.group["profiles"]["t"])[1]
        t_min = floor(Int, 5/6 * t_max)
        period = ds.group["profiles"]["t"][t_max]/3600/6        
        for (qw,i) in zip(cat(env_qw, updrafts_qw, dims=1), 1:6)
            q=0
            for x in t_min:t_max
                q = ds.group["profiles"][qw][:, x]
            end
            q = q/(t_max-t_min+1)
            subplots[i].title = qw*" last "*string(period)*" hour(s)"
            l = lines!(subplots[i], q, height)
            push!(qw_lines[i], l)
        end
    end
    for i in 1:5
        axislegend(subplots[i], qw_lines[i], [string(j) for j in values], type, position = :rt)
    end
    
    save(dir*case*"\\qw_condensed\\"*case*"_env_updraft_condensed_"*type*".png", fig)    
    return fig
end

function plot_w_updraft_condensed(case::String)
    fig = Figure(resolution =(1800,1200))
    subplots = [Axis(fig[i,j], 
        palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
        for j in 1:3, i in 1:2]
    qw_lines = [[] for i in 1:6]
    
    for (type, i) in zip(rates, 1:5)
        values = 0.0:0.5:2.0
        if type == "depsub"
            values = 0.0:0.5:1.5
        end
        for r in values
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(r)*".nc")
            height = ds.group["profiles"]["zf"]*1e-3
            t_max = size(ds.group["profiles"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["profiles"]["t"][t_max]/3600/6        
            q=0
            for x in t_min:t_max
                q = ds.group["profiles"]["updraft_w"][:, x]
            end
            q = q/(t_max-t_min+1)
            subplots[i].title = type*" updraft_w last "*string(period)*" hour(s)"
            l = lines!(subplots[i], q, height)
            push!(qw_lines[i], l)
        end

        axislegend(subplots[i], qw_lines[i], [string(j) for j in values], type, position = :rt)
    end        
       
    save(dir*case*"\\qw_condensed\\"*case*"_updraft_w_condensed.png", fig)    
    return fig
end

for case in case_names
    plot_flux_condensed(case)
    plot_wp_condensed(case)
    plot_cloud_condensed(case)
    plot_w_updraft_condensed(case)
    for type in rates
        if case!="Rico"
            plot_qw(case, type)
        end
        plot_all_wp(case, type)
        plot_all_clouds(case, type)
        plot_env_updraft_condensed(case, type)
        plot_qw_condensed(case, type)
    end
end

