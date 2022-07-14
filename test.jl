using NCDatasets
using CairoMakie
using Glob
using ColorSchemes

# cases  
case_names = ["DYCOMS_RF02","Rico","TRMM_LBA"]
# rates (no autocon currently)
rates = ["autocon", "accretion", "depsub", "evap", "melt"]
# rates = ["autocon"]


dir = "C:\\Users\\khanh\\SURF\\TurbulenceConvection.jl\\"


#wp plotting function
function plot_wp(path::String, plot::Axis, ds::Dataset)
    l =lines!(plot,
        ds.group["timeseries"]["t"]/3600, 
        ds.group["timeseries"][path])
    return l
end

#plotting all wp
case_figs = []
for case in case_names
    type_figs = []
    for type in rates
        fig = Figure(resolution =(1600,1200))
        subplots = [Axis(fig[i,j], 
            xlabel = "hours(h)", 
            ylabel = "density(kg/m^3)",
            palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
            for i in 1:2, j in 1:2]
        wp_lines = [[],[],[],[],[]]
        subplots[1].title = case*" lwp_mean variable "*type
        subplots[2].title = case*" iwp_mean variable "*type
        subplots[3].title = case*" rwp_mean variable "*type
        subplots[4].title = case*" swp_mean variable "*type
        for v in 0.0:0.5:2.0
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")      
            push!(wp_lines[1], plot_wp("lwp_mean",subplots[1],ds))
            push!(wp_lines[2], plot_wp("iwp_mean",subplots[2],ds))
            push!(wp_lines[3], plot_wp("rwp_mean",subplots[3],ds))
            push!(wp_lines[4], plot_wp("swp_mean",subplots[4],ds))
        end
        for i in 1:4
            axislegend(subplots[i], wp_lines[i], [string(j) for j in 0.0:0.5:2.0], type, position = :lt)
        end
        push!(type_figs, fig)
        save(dir*case*"\\wp_plots\\"*case*"_wp_"*type*".png", fig)
    end
    push!(case_figs, type_figs)
end


function plot_clouds(path::String, plot::Axis, ds::Dataset)
    l =lines!(plot,
        ds.group["timeseries"]["t"]/3600, 
        ds.group["timeseries"][path])
    return l

end

#plotting cloud functions (fix cloud cover)
for case in case_names
    type_figs = []
    for type in rates
        fig = Figure(resolution =(1600,1200))
        subplots = [Axis(fig[i,j], 
            xlabel = "hours(h)", 
            ylabel = "height(km)",
            palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
            for i in 1:2, j in 1:2]
        wp_lines = [[],[],[],[],[]]
        subplots[1].title = case*" cloud_cover_mean variable "*type
        subplots[2].title = case*" cloud_top_mean variable "*type
        subplots[3].title = case*" cloud_base_mean variable "*type
        for v in 0.0:0.5:2.0
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(v)*".nc")      
            push!(wp_lines[1], plot_wp("cloud_cover_mean",subplots[1],ds))
            push!(wp_lines[2], plot_wp("cloud_top_mean",subplots[2],ds))
            push!(wp_lines[3], plot_wp("cloud_base_mean",subplots[3],ds))
        end
        for i in 1:3
            axislegend(subplots[i], wp_lines[i], [string(j) for j in 0.0:0.5:2.0], type, position = :lt)
        end
        push!(type_figs, fig)
        save(dir*case*"\\cloud_plots\\"*case*"_clouds_"*type*".png", fig)
    end
    push!(case_figs, type_figs)
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

# setting case (Rico cannot be graphed)
case = "TRMM_LBA"
# setting rate type (For TRMM_LBA, depsub should not be graphed for 2.0)
# note that the size may differ between different rate values and thus might need to be specified
type = "accretion"

for type in rates
    # plottiing
    print(type*"\n")
    figs = [Figure(resolution=(1400, 1600)) for i in 1:4]
    base = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_1.0.nc")
    
    # plotting normal rates
    plots = [[figs[i][j,k] for j in 2:3, k in 1:2] for i in 1:4]
    for (qw,i) in zip(["ql","qi", "qr", "qs"],1:4) 
        plot_qw_timeseries(qw*"_mean", figs[i][1,1], base, "1.0")
    end


    #getting differences
    diffs = [[] for i in 1:4]
    for (v, j) in zip(["0.0", "0.5", "1.5", "2.0"], 1:4)
        print(v*"\n")
        ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*v*".nc")
        for (qw,i) in zip(["ql","qi", "qr", "qs"], 1:4)
            height = ds.group["profiles"]["zc"]*1e-3
            time = ds.group["profiles"]["t"]*3600
            q = transpose(ds.group["profiles"][qw*"_mean"][:,1:334]-base.group["profiles"][qw*"_mean"])*1e3
            push!(diffs[i], (time, height, q))
        end
    end

    #plotting differences
    for (qw,i) in zip(["ql","qi", "qr", "qs"], 1:4)
        print(qw*"\n")
        extremas = map(extrema, [diffs[3][i][3] for i in 1:4])
        global_min = minimum(t->first(t), extremas)
        global_max = maximum(t->last(t), extremas)
        m = max(abs(global_min), abs(global_max))
        range = (-m, m)
        for (v, j) in zip(["0.0", "0.5", "1.5", "2.0"], 1:4)
            ax, hm = heatmap(plots[i][j], diffs[i][j][1][1:334], diffs[i][j][2], diffs[i][j][3], colormap = :bwr, colorrange = range)
            ax.xlabel = "hours(h)"
            ax.ylabel = "height(km)"
            ax.title = qw*"_mean time evolution difference "*type*"="*v        
        end
        cb = Colorbar(figs[i][2:3, 3], limits=range, colormap = :bwr)
        save(dir*case*"\\qw_plots\\"*case*"_"*type*"_"*qw*"_mean.png", figs[i]) 
    end
end

base = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_1.0.nc")

base.group["profiles"]["updraft_w"]


for case in case_names
    for type in rates
        print(case*"\n")
        fig = Figure(resolution =(600,600))
        subplots = [Axis(fig[i,j], 
            xlabel = "", 
            ylabel = "height(km)",
            palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
            for j in 1:1, i in 1:1]
        qw_lines = [[] for i in 1:5]
        for r in 0.0:0.5:2.0
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(r)*".nc")
            height = ds.group["profiles"]["zc"]*1e-3
            t_max = size(ds.group["profiles"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["profiles"]["t"][t_max]/3600/6
            for (qw,i) in zip(["v","qt","ql","qi","qr","qs"], 1:1)
                q=0
                for x in t_min:t_max
                    q = ds.group["profiles"][qw*"_mean"][:, x]
                end
                q = q/(t_max-t_min+1)
                subplots[i].title = qw*" last "*string(period)*" hour(s)"
                l = lines!(subplots[i], q, height)
                push!(qw_lines[i], l)
            end
        end
        for i in 1:1
            axislegend(subplots[i], qw_lines[i], [string(j) for j in 0.0:0.5:2.0], type, position = :lt)
        end

        save(dir*case*"\\qw_condensed\\"*case*"_v_condensed_"*type*".png", fig)
    end
end

base.group["timeseries"]
base.group["profiles"]


w_paths = ["lwp", "iwp", "rwp", "swp"]
clouds = ["cloud_cover", "cloud_base", "cloud_top"]


for case in case_names
    for type in rates
        paths = [[] for i in 1:4]
        for r in 0.0:0.5:2.0
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(r)*".nc")
            t_max = size(ds.group["timeseries"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["timeseries"]["t"][t_max]/3600/6
            for (wp,i) in zip(w_paths, 1:4)
                p = ds.group["timeseries"][wp*"_mean"][t_min:t_max]
                push!(paths[i],sum(p)/(t_max-t_min+1))
            end
        end

        fig = Figure(resolution =(600,600))
        subplots = [Axis(fig[i,j], 
                    xlabel = "rate", 
                    ylabel = "",
                    palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
                    for j in 1:2, i in 1:2]

        for (wp,i) in zip(clouds, 1:4)
            l = lines!(subplots[i],0.0:0.5:2.0, convert(Array{Float32}, paths[i]))
            subplots[i].title = type*" "*wp*"_mean "
        end
        save(dir*case*"\\cloud_condensed\\"*case*"_cloud_condensed_"*type*".png", fig)
    end
end


for case in case_names
    for type in rates
        paths = [[] for i in 1:4]
        for r in 0.0:0.5:2.0
            ds = Dataset(dir*"\\"*case*"\\"*case*"_"*type*"_rate_"*string(r)*".nc")
            t_max = size(ds.group["timeseries"]["t"])[1]
            t_min = floor(Int, 5/6 * t_max)
            period = ds.group["timeseries"]["t"][t_max]/3600/6
            for (wp,i) in zip(["rain_flux","snow_flux"], 1:2)
                p = ds.group["profiles"][wp][1,t_min:t_max]
                push!(paths[i],sum(p)/(t_max-t_min+1))
            end
        end

        fig = Figure(resolution =(600,600))
        subplots = [Axis(fig[i,j], 
                    xlabel = "rate", 
                    ylabel = "",
                    palette = (color = [(:red, 0.5), (:orange, 0.5), (:green, 0.5), (:blue, 0.5), (:cyan, 0.5)],))
                    for j in 1:2, i in 1:2]

        for (wp,i) in zip(["rain_flux","snow_flux"], 1:2)
            l = lines!(subplots[i],0.0:0.5:2.0, convert(Array{Float32}, paths[i]))
            subplots[i].title = type*" "*wp
        end
        save(dir*case*"\\flux_condensed\\"*case*"_flux_condensed_"*type*".png", fig)
    end
end





