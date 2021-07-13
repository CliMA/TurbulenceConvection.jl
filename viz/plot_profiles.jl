using NCDatasets
using Plots
using Statistics
ENV["GKSwstype"] = "nul"

function plot_profiles(ds_filename::String, group; kwargs...)
    Dataset(ds_filename) do ds
        plot_profiles(ds, dirname(ds_filename), group; kwargs...)
    end
end

function plot_profiles(ds,
    output_dir,
    group;
    skip_fields = [],
    fields_only = nothing,
    profiles_only = false,
    skip_contours = true,
    skip_all_zeros_data = true,
)
    mkpath(output_dir)
    @assert any((
        group == "profiles",
        group == "reference",
        group == "timeseries",
    ))
    z = ds.group["profiles"]["z"][:]
    z_half = ds.group["profiles"]["z"][:]
    time = ds.group[group]["t"]
    t_end = length(time)
    t_start = ((time .- 3600) .> 0)[1][1]

    @show keys(ds.group[group])
    for var in keys(ds.group[group])
        var=="z" && continue
        var=="z_half" && continue
        var=="t" && continue
        var in skip_fields && continue

        values = ds.group[group][var][:]
        figname = "profile_"*var
        yvals = mean(values[t_start+1:t_end,:], dims = 1)
        yvals = reshape(yvals, length(yvals))
        if skip_all_zeros_data && all(values .≈ 0)
            @warn "Skipping $var, since all data is 0"
            continue
        end
        plot(
            yvals,
            z;
            xlabel = var,
            ylabel = "height (m)"
        )
        xlabel!(var)
        ylabel!("height (m)")
        fname = joinpath(output_dir,figname*".png")
        isfile(fname) && rm(fname; force=true)
        savefig(fname)

        skip_contours && continue
        figname = "contours_"*var
        contourf(
            time[:], z, values';
            xlabel = var,
            ylabel = "height (m)",
            c = :viridis
        )
        xlabel!(var)
        ylabel!("height (m)")
        savefig(joinpath(output_dir,figname*".png"))
    end
end

