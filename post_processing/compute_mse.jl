import Plots
import OrderedCollections
using Test

import Dates
import JSON

import NCDatasets
const NC = NCDatasets

import StatsBase
import Dierckx
import PrettyTables
import Printf

ENV["GKSwstype"] = "nul"

include(joinpath(@__DIR__, "..", "integration_tests", "artifact_funcs.jl"))

function get_time(ds, var)
    if haskey(ds, var)
        return ds[var][:]
    elseif haskey(ds.group["timeseries"], var)
        return ds.group["profiles"][var][:]
    end
    error("No key for $var found in the nc file.")
end

append_dict(::Nothing, sym::Symbol, dict = Dict()) = dict
function append_dict(filename::AbstractString, sym::Symbol, dict = Dict())
    if isfile(filename)
        dict[sym] = NC.Dataset(filename, "r")
    end
    return dict
end

function find_latest_dataset_folder(; dir = pwd())
    matching_paths = String[]
    for file in readdir(dir)
        !ispath(joinpath(dir, file)) && continue
        push!(matching_paths, joinpath(dir, file))
    end
    isempty(matching_paths) && return ""
    # sort by timestamp
    sorted_paths = sort(matching_paths; by = f -> Dates.unix2datetime(stat(f).mtime))
    return pop!(sorted_paths)
end

function compute_mse_wrapper(
    case_name,
    best_mse,
    ds_tc_filename,
    ds_tc_main_filename = nothing;
    ds_les_filename = joinpath(pycles_output_dataset_folder(), "$case_name.nc"),
    ds_scm_filename = joinpath(scampy_output_dataset_folder(), "$case_name.nc"),
    plot_dir = joinpath(dirname(ds_tc_filename), "comparison"),
    kwargs...,
)

    # Note: cluster_data_prefix is also defined in utils/move_output.jl
    if haskey(ENV, "BUILDKITE_COMMIT")
        cluster_data_prefix = "/resnick/scratch/esm/slurm-buildkite/turbulenceconvection-main"
        path = find_latest_dataset_folder(; dir = cluster_data_prefix)

        # TODO: make this more robust in case folder/file changes
        folder_name = joinpath("Output.$case_name.01", "stats")
        ds_tc_main_filename = joinpath(path, folder_name, "Stats.$case_name.nc")
        @info "TC.jl main dataset:`$ds_tc_main_filename`"
    end
    all_args = (case_name, best_mse, plot_dir)

    ds_dict = append_dict(ds_les_filename, :ds_pycles)
    ds_dict = append_dict(ds_scm_filename, :ds_scampy, ds_dict)
    ds_dict = append_dict(ds_tc_filename, :ds_tc, ds_dict)
    ds_dict = append_dict(ds_tc_main_filename, :ds_tc_main, ds_dict)

    computed_mse_success = false

    local computed_mse
    try
        computed_mse = compute_mse(all_args...; ds_dict, kwargs...)
        computed_mse_success = true
    finally
        for ds in values(ds_dict)
            close(ds)
        end
    end
    if computed_mse_success
        return computed_mse
    else
        error("MSE computation unsuccessful.")
    end
end


function compute_mse(
    case_name,
    best_mse,
    plot_dir;
    ds_dict,
    plot_comparison = true,
    group_figs = true,
    t_start,
    t_stop,
    skip_comparison = false,
)

    TC = TurbulenceConvection
    ds_scampy = haskey(ds_dict, :ds_scampy) ? ds_dict[:ds_scampy] : nothing
    ds_tc_main = haskey(ds_dict, :ds_tc_main) ? ds_dict[:ds_tc_main] : nothing
    ds_pycles = haskey(ds_dict, :ds_pycles) ? ds_dict[:ds_pycles] : nothing
    ds_tc = haskey(ds_dict, :ds_tc) ? ds_dict[:ds_tc] : nothing
    # Make pycles optional
    have_pycles_ds = ds_pycles ≠ nothing
    have_scampy_ds = ds_scampy ≠ nothing
    have_tc_main = ds_tc_main ≠ nothing
    if !have_pycles_ds
        ds_pycles = ds_tc
        @warn "No PyCLES data, using TC.jl instead"
    end
    if !have_scampy_ds
        ds_scampy = ds_tc
        @warn "No SCAMPy data, using TC.jl instead"
    end
    if !have_tc_main
        ds_tc_main = ds_tc
        @warn "No TC.jl main data, using TC.jl instead"
    end
    mse = OrderedCollections.OrderedDict()
    time_tcc = get_time(ds_tc, "t")
    time_les = get_time(ds_pycles, "t")
    time_tcm = get_time(ds_tc_main, "t")
    time_scm = get_time(ds_scampy, "t")
    n_time_points = length(time_tcc)

    mkpath(plot_dir)
    # Ensure domain matches:
    z_les = Array(TC.get_nc_data(ds_pycles, "zc"))
    z_tcc_c = Array(TC.get_nc_data(ds_tc, "zc"))
    z_tcc_f = Array(TC.get_nc_data(ds_tc, "zf"))
    z_tcm_c = Array(TC.get_nc_data(ds_tc_main, "zc"))
    z_tcm_f = Array(TC.get_nc_data(ds_tc_main, "zf"))
    z_scm = Array(TC.get_nc_data(ds_scampy, "zc"))
    n_grid_points = length(z_tcc_c)
    @info "z extrema (les,scm,tcm,tcc): $(extrema(z_les)), $(extrema(z_scm)), $(extrema(z_tcm_c)), $(extrema(z_tcc_c))"
    @info "n-grid points (les,scm,tcm,tcc): $(length(z_les)), $(length(z_scm)), $(length(z_tcm_c)), $(length(z_tcc_c))"

    @info "time extrema (les,scm,tcm,tcc): $(extrema(time_les)), $(extrema(time_scm)), $(extrema(time_tcm)), $(extrema(time_tcc))"
    @info "n-time points (les,scm,tcm,tcc): $(length(time_les)), $(length(time_scm)), $(length(time_tcm)), $(length(time_tcc))"

    # Find the nearest matching final time:
    t_cmp = min(time_tcc[end], time_tcm[end], time_les[end], time_scm[end], t_stop)
    @info "time compared: $t_cmp"

    if length(time_tcc) == 1
        msg = "IO was collected for a single point in time."
        msg *= "This script (compute_mse) requires comparing "
        msg *= "space-time fields with multiple timesteps."
        error(msg)
    end

    # Accidentally running a short simulation
    # could improve MSE. So, let's test that
    # we run for at least t_compare. We should
    # increase this as we can reach higher CFL.
    # @test t_cmp >= t_compare

    # Ensure z_tcc and fields are consistent lengths:
    @test length(z_tcc_c) == n_grid_points
    tcc_variables = []
    computed_mse = []
    table_best_mse = []
    mse_reductions = []
    data_scales_scm = []
    data_scales_les = []
    data_scales_tcc = []
    data_scales_tcm = []
    pycles_weight = []
    plots_dict = Dict()
    plots_dict["profiles"] = Dict()
    plots_dict["contours"] = Dict()
    plot_attr = Dict()
    plot_attr[true] = Dict()
    plot_attr[false] = Dict()

    local fig_height

    for tc_var in keys(best_mse)
        data_les_arr = TC.get_nc_data(ds_pycles, tc_var)
        data_tcm_arr = TC.get_nc_data(ds_tc_main, tc_var)
        data_tcc_arr = TC.get_nc_data(ds_tc, tc_var)
        data_scm_arr = TC.get_nc_data(ds_scampy, tc_var)
        # Only compare fields that exist in the nc files
        missing_les_var = isnothing(data_les_arr)
        missing_tcm_var = isnothing(data_tcm_arr)
        missing_tcc_var = isnothing(data_tcc_arr)
        missing_scm_var = isnothing(data_scm_arr)

        coord_name = first(NC.dimnames(data_tcc_arr))
        coord_name == "zf" || coord_name == "zc" || error("Bad coord_name")

        if coord_name == "zf"
            z_tcc = z_tcc_f
        elseif coord_name == "zc"
            z_tcc = z_tcc_c
        end
        if have_tc_main
            if coord_name == "zf"
                z_tcm = z_tcm_f
            elseif coord_name == "zc"
                z_tcm = z_tcm_c
            end
        else
            z_tcm = z_tcc
        end
        if !have_scampy_ds
            z_scm = z_tcc
        end
        if !have_pycles_ds
            z_les = z_tcc
        end

        if missing_les_var
            @warn "les missing data for variable $tc_var, filling with zeros"
            data_les_arr = zeros(length(z_les), length(time_les))
            warn_msg_les = " Warning: missing data"
        else
            warn_msg_les = ""
        end
        if missing_tcm_var
            @warn "tcm missing data for variable $tc_var, filling with zeros"
            data_tcm_arr = zeros(length(z_tcm), length(time_tcm))
            warn_msg_tcm = " Warning: missing data"
        else
            warn_msg_tcm = ""
        end
        if missing_tcc_var
            @warn "tcc missing data for variable $tc_var, filling with zeros"
            data_tcc_arr = zeros(length(z_tcc), length(time_tcc))
            warn_msg_tcc = " Warning: missing data"
        else
            warn_msg_tcc = ""
        end
        if missing_scm_var
            @warn "scm missing data for variable $tc_var, filling with zeros"
            data_scm_arr = zeros(length(z_scm), length(time_scm))
            warn_msg_scm = " Warning: missing data"
        else
            warn_msg_scm = ""
        end

        # TC and TC main are the same locally
        if haskey(ENV, "BUILDKITE_COMMIT")
            new_variable = all((missing_les_var, missing_tcm_var, missing_scm_var))
        else
            new_variable = all((missing_les_var, missing_scm_var))
        end

        data_les_arr = Array(data_les_arr)'
        data_tcm_arr = Array(data_tcm_arr)'
        data_tcc_arr = Array(data_tcc_arr)'
        data_scm_arr = Array(data_scm_arr)'
        push!(tcc_variables, tc_var)

        @info "Assembling plots for variable:`$tc_var`"

        @debug "     Data sizes (les,scm,tcc,tcm): $(size(data_les_arr)), $(size(data_scm_arr)), $(size(data_tcc_arr)), $(size(data_tcm_arr))"
        # Scale the data for comparison
        push!(pycles_weight, "1")

        # Interpolate data
        steady_data = size(data_les_arr, 1) == 1
        if steady_data
            data_les_cont = Dierckx.Spline1D(z_les, vec(data_les_arr); k = 1)
            data_tcc_cont = Dierckx.Spline1D(z_tcc, vec(data_tcc_arr); k = 1)
            data_tcm_cont = Dierckx.Spline1D(z_tcm, vec(data_tcm_arr); k = 1)
            data_scm_cont = Dierckx.Spline1D(z_scm, vec(data_scm_arr); k = 1)
            data_les_cont_mapped = map(z -> data_les_cont(z), z_tcc)
            data_tcm_cont_mapped = map(z -> data_tcm_cont(z), z_tcc)
            data_tcc_cont_mapped = map(z -> data_tcc_cont(z), z_tcc)
            data_scm_cont_mapped = map(z -> data_scm_cont(z), z_tcc)
        else # unsteady data
            data_les_cont = Dierckx.Spline2D(time_les, z_les, data_les_arr; kx = 1, ky = 1)
            data_tcm_cont = Dierckx.Spline2D(time_tcm, z_tcm, data_tcm_arr; kx = 1, ky = 1)
            data_tcc_cont = Dierckx.Spline2D(time_tcc, z_tcc, data_tcc_arr; kx = 1, ky = 1)
            data_scm_cont = Dierckx.Spline2D(time_scm, z_scm, data_scm_arr; kx = 1, ky = 1)
            R = range(t_start, t_cmp; length = 50)
            data_les_cont_mapped = map(z_tcc) do z
                StatsBase.mean(map(t -> data_les_cont(t, z), R))
            end
            data_tcm_cont_mapped = map(z_tcc) do z
                StatsBase.mean(map(t -> data_tcm_cont(t, z), R))
            end
            data_tcc_cont_mapped = map(z_tcc) do z
                StatsBase.mean(map(t -> data_tcc_cont(t, z), R))
            end
            data_scm_cont_mapped = map(z_tcc) do z
                StatsBase.mean(map(t -> data_scm_cont(t, z), R))
            end

        end

        # Plot comparison
        if plot_comparison
            p = Plots.plot()
            if have_pycles_ds
                Plots.plot!(
                    data_les_cont_mapped,
                    z_tcc ./ 10^3,
                    title = tc_var,
                    ylabel = "z [km]",
                    label = "PyCLES$warn_msg_les",
                )
            end
            Plots.plot!(
                data_scm_cont_mapped,
                z_tcc ./ 10^3,
                title = tc_var,
                ylabel = "z [km]",
                label = "SCAMPy$warn_msg_scm",
            )
            if have_tc_main
                Plots.plot!(
                    data_tcm_cont_mapped,
                    z_tcc ./ 10^3,
                    title = tc_var,
                    ylabel = "z [km]",
                    label = "TC.jl (main)$warn_msg_tcm",
                )
            end
            plots_dict["profiles"][tc_var] = Plots.plot!(
                data_tcc_cont_mapped,
                z_tcc ./ 10^3,
                title = tc_var,
                ylabel = "z [km]",
                label = "TC.jl$warn_msg_tcc",
            )

            clims_min = min(minimum(data_scm_arr), minimum(data_tcc_arr))
            clims_max = max(maximum(data_scm_arr), maximum(data_tcc_arr))
            clims = (clims_min, clims_max)

            width_to_height_ratio = have_tc_main ? 15 / 10 : 15 / 10
            fig_height = 1500

            plot_attr[true]["contour_1_kwargs"] = (
                bottom_margin = 0 * Plots.PlotMeasures.px,
                left_margin = 40 * Plots.PlotMeasures.px,
                right_margin = 0 * Plots.PlotMeasures.px,
                top_margin = 0 * Plots.PlotMeasures.px,
                xticks = false,
                colorbar = false,
                size = (width_to_height_ratio * fig_height, fig_height),
                ylabel = "height (km)",
            )
            plot_attr[true]["contour_2_kwargs"] = (
                bottom_margin = 0 * Plots.PlotMeasures.px,
                left_margin = 0 * Plots.PlotMeasures.px,
                right_margin = 0 * Plots.PlotMeasures.px,
                top_margin = 0 * Plots.PlotMeasures.px,
                xticks = false,
                yticks = false,
                colorbar = false,
                size = (width_to_height_ratio * fig_height, fig_height),
            )
            plot_attr[true]["contour_3_kwargs"] = (
                bottom_margin = 0 * Plots.PlotMeasures.px,
                left_margin = 0 * Plots.PlotMeasures.px,
                right_margin = -20 * Plots.PlotMeasures.px,
                top_margin = 0 * Plots.PlotMeasures.px,
                xticks = false,
                yticks = false,
                colorbar = true,
                size = (width_to_height_ratio * fig_height, fig_height), # extra space for colorbar
            )
            plot_attr[false]["contour_1_kwargs"] = (ylabel = "height (km)",)
            plot_attr[false]["contour_2_kwargs"] = (ylabel = "height (km)",)
            plot_attr[false]["contour_3_kwargs"] = (ylabel = "height (km)",)

            clim_kwarg = clims_min == clims_max ? () : (; clims = clims)

            p1 = Plots.contourf(
                time_scm ./ 3600,
                z_scm ./ 10^3,
                data_scm_arr';
                c = :viridis,
                clim_kwarg...,
                title = "$tc_var (SCAMPy)$warn_msg_scm",
                plot_attr[group_figs]["contour_1_kwargs"]...,
            )
            if have_tc_main
                p2 = Plots.contourf(
                    time_tcm ./ 3600,
                    z_tcm ./ 10^3,
                    data_tcm_arr';
                    c = :viridis,
                    clim_kwarg...,
                    title = "$tc_var (TC.jl main)$warn_msg_tcm",
                    plot_attr[group_figs]["contour_2_kwargs"]...,
                )
            end
            p3 = Plots.contourf(
                time_tcc ./ 3600,
                z_tcc ./ 10^3,
                data_tcc_arr';
                c = :viridis,
                clim_kwarg...,
                title = have_tc_main ? "$tc_var (TC.jl PR)$warn_msg_tcc" : "$tc_var (TC.jl)$warn_msg_tcc",
                plot_attr[group_figs]["contour_3_kwargs"]...,
            )
            if group_figs
                if have_tc_main
                    widths = [0.29, 0.31]
                    push!(widths, 1 - sum(widths))
                    plots_dict["contours"][tc_var] = Plots.plot(p1, p2, p3; layout = Plots.grid(1, 3, widths = widths))
                else
                    widths = [0.45]
                    push!(widths, 1 - sum(widths))
                    plots_dict["contours"][tc_var] = Plots.plot(p1, p3; layout = Plots.grid(1, 2, widths = widths))
                end
            else
                if have_tc_main
                    plots_dict["contours"][tc_var] = Plots.plot(p1, p2, p3; layout = (3, 1))
                else
                    plots_dict["contours"][tc_var] = Plots.plot(p1, p3; layout = (2, 1))
                end
            end
        end

        # Compute data scale
        data_scale_tcm = sum(abs.(data_tcm_arr)) / length(data_tcm_arr)
        data_scale_tcc = sum(abs.(data_tcc_arr)) / length(data_tcc_arr)
        data_scale_scm = sum(abs.(data_scm_arr)) / length(data_scm_arr)
        data_scale_les = sum(abs.(data_les_arr)) / length(data_les_arr)
        push!(data_scales_tcm, data_scale_tcm) # TODO: should we add this to the mse table?
        push!(data_scales_tcc, data_scale_tcc)
        push!(data_scales_scm, data_scale_scm)
        push!(data_scales_les, data_scale_les) # TODO: should we add this to the mse table?

        # Compute mean squared error (mse)
        if have_pycles_ds # LES takes first precedence
            mse_single_var = sum((data_les_cont_mapped .- data_tcc_cont_mapped) .^ 2)
            data_scale_used = data_scale_les
            scaled_mse = mse_single_var / data_scale_used^2 # Normalize by data scale
        elseif have_tc_main # TC.jl main takes second precedence
            mse_single_var = sum((data_tcm_cont_mapped .- data_tcc_cont_mapped) .^ 2)
            data_scale_used = data_scale_tcm
            scaled_mse = mse_single_var / data_scale_used^2 # Normalize by data scale
        elseif have_scampy_ds # SCAMPy takes third precedence
            mse_single_var = sum((data_scm_cont_mapped .- data_tcc_cont_mapped) .^ 2)
            data_scale_used = data_scale_scm
            scaled_mse = mse_single_var / data_scale_used^2 # Normalize by data scale
        elseif skip_comparison
            mse[tc_var] = "NA"
            continue
        else
            error("No dataset to compute MSE")
        end

        haskey(best_mse, tc_var) || error("No key found in best_mse for variable $tc_var")

        if new_variable || best_mse[tc_var] isa String
            mse[tc_var] = "NA"
            push!(mse_reductions, "NA")
            push!(table_best_mse, "NA")
            push!(computed_mse, "NA")
        else
            new_variable && @warn "$tc_var must be added to the MSE tables."
            mse[tc_var] = scaled_mse
            push!(mse_reductions, (best_mse[tc_var] - mse[tc_var]) / best_mse[tc_var] * 100)
            push!(table_best_mse, best_mse[tc_var])
            push!(computed_mse, mse[tc_var])
        end
    end

    save_plots(plot_dir, plots_dict; group_figs, have_tc_main, fig_height, case_name)

    if skip_comparison
        return mse
    end

    # Tabulate output
    header = (
        ["Variable", "Weight", "Data scale", "Data scale", "MSE", "MSE", "MSE"],
        ["TC.jl (EDMF)", "PyCLES", "tcc", "scm", "Computed", "Best", "Reduction (%)"],
    )
    table_data = hcat(
        tcc_variables,
        pycles_weight,
        data_scales_tcc,
        data_scales_scm,
        computed_mse,
        table_best_mse,
        mse_reductions,
    )

    @info Printf.@sprintf("case_name comparison: %s at time t=%s\n", case_name, t_cmp)
    hl_worsened_mse = PrettyTables.Highlighter(
        (data, i, j) -> !sufficient_mse(data[i, 5], data[i, 6]) && j == 5,
        PrettyTables.crayon"red bold",
    )
    hl_worsened_mse_reduction = PrettyTables.Highlighter(
        (data, i, j) -> !sufficient_mse(data[i, 5], data[i, 6]) && j == 7,
        PrettyTables.crayon"red bold",
    )
    hl_improved_mse = PrettyTables.Highlighter(
        (data, i, j) -> sufficient_mse(data[i, 5], data[i, 6]) && j == 7,
        PrettyTables.crayon"green bold",
    )
    PrettyTables.pretty_table(
        table_data;
        header,
        formatters = PrettyTables.ft_printf("%.16e", 5:6),
        header_crayon = PrettyTables.crayon"yellow bold",
        subheader_crayon = PrettyTables.crayon"green bold",
        highlighters = (hl_worsened_mse, hl_improved_mse, hl_worsened_mse_reduction),
        crop = :none,
    )

    return mse
end

function save_plots(plot_dir, plots_dict; group_figs = true, have_tc_main, fig_height, case_name)
    vars_to_skip = [
        "thetal_mean",
        "updraft_thetal",
        "u_mean",
        "v_mean",
        "qt_mean",
        "updraft_qt",
        "temperature_mean",
        "total_flux_h",
        "total_flux_qt",
        "massflux",
    ]
    all_contours = plots_dict["contours"]
    filter!(p -> !occursin("var", p.first), all_contours) # skip higher order moment contours
    filter!(p -> !any(map(x -> occursin(x, p.first), vars_to_skip)), all_contours)
    if group_figs
        all_contours = collect(values(all_contours))
        Plots.plot!(all_contours[end], xticks = true, xlabel = "Time (hr)", bottom_margin = 20 * Plots.PlotMeasures.px)

        n_plots = length(all_contours)
        n_cols = 1
        n_rows = ceil(Int, n_plots / n_cols)
        @info "     Saving file:`$(joinpath(plot_dir, "contours.png"))`"
        title = Plots.plot(
            title = case_name,
            grid = false,
            showaxis = false,
            xticks = false,
            yticks = false,
            bottom_margin = -20 * Plots.PlotMeasures.px,
        )
        layout = Plots.@layout [a{0.01h}; Plots.grid(n_rows, n_cols)]

        Plots.plot(
            title,
            all_contours...;
            layout = layout,
            framestyle = :box,
            margin = 20 * Plots.PlotMeasures.px,
            dpi = 200,
        )
        Plots.savefig(joinpath(plot_dir, "contours.png"))

        width_to_height_ratio = 15 / 10
        all_profiles = collect(values(plots_dict["profiles"]))
        n_plots = length(all_profiles)
        n_cols = 5
        n_rows = ceil(Int, n_plots / n_cols)
        @info "     Saving file:`$(joinpath(plot_dir, "profiles.png"))`"
        n_empty_plots = n_rows * n_cols - n_plots
        all_profiles =
            [all_profiles..., ntuple(i -> Plots.plot(; framestyle = :none, xticks = false), n_empty_plots)...]
        for (k, I) in enumerate(CartesianIndices((1:n_cols, 1:n_rows)))
            k > length(all_profiles) && continue
            j, i = (Tuple(I)...,)
            if j == 1 # left plots
                Plots.plot!(
                    all_profiles[k],
                    left_margin = 50 * Plots.PlotMeasures.px,
                    bottom_margin = -2 * Plots.PlotMeasures.px,
                    right_margin = -2 * Plots.PlotMeasures.px,
                    top_margin = -2 * Plots.PlotMeasures.px,
                )
            elseif j == n_cols # right plots
                Plots.plot!(
                    all_profiles[k],
                    yticks = false,
                    ylabel = "",
                    bottom_margin = -2 * Plots.PlotMeasures.px,
                    left_margin = -2 * Plots.PlotMeasures.px,
                    right_margin = -2 * Plots.PlotMeasures.px,
                    top_margin = -2 * Plots.PlotMeasures.px,
                )
            else # middle plots
                Plots.plot!(
                    all_profiles[k],
                    yticks = false,
                    ylabel = "",
                    bottom_margin = -2 * Plots.PlotMeasures.px,
                    left_margin = -2 * Plots.PlotMeasures.px,
                    right_margin = -2 * Plots.PlotMeasures.px,
                    top_margin = -2 * Plots.PlotMeasures.px,
                )
            end
        end
        title = Plots.plot(
            title = case_name,
            grid = false,
            showaxis = false,
            xticks = false,
            yticks = false,
            bottom_margin = -20 * Plots.PlotMeasures.px,
        )
        layout = Plots.@layout [a{0.01h}; Plots.grid(n_rows, n_cols)]
        Plots.plot(
            title,
            all_profiles...;
            layout = layout,
            size = (width_to_height_ratio * fig_height, fig_height),
            framestyle = :box,
            margin = 20 * Plots.PlotMeasures.px,
            dpi = 200,
        )
        Plots.savefig(joinpath(plot_dir, "profiles.png"))
    else
        for tc_var in keys(all_contours)
            @info "     Saving file:`$(joinpath(plot_dir, "contours_$tc_var.png"))`"
            Plots.plot(all_contours[tc_var])
            Plots.savefig(joinpath(plot_dir, "contours_$tc_var.png"))
        end

        all_profiles = plots_dict["profiles"]
        for tc_var in keys(all_profiles)
            @info "     Saving file:`$(joinpath(plot_dir, "profiles_$tc_var.png"))`"
            Plots.plot(all_profiles[tc_var])
            Plots.savefig(joinpath(plot_dir, "profiles_$tc_var.png"))
        end
    end
end

sufficient_mse(computed_mse::Real, best_mse::Real) = computed_mse <= best_mse + sqrt(eps())
# New variables report string MSEs (NA), so those should always be sufficient.
sufficient_mse(computed_mse::String, best_mse::String) = true
sufficient_mse(computed_mse::Real, best_mse::String) = true
sufficient_mse(computed_mse::String, best_mse::Real) = true

function test_mse(computed_mse, best_mse, key)
    mse_not_regressed = sufficient_mse(computed_mse[key], best_mse[key])
    mse_not_regressed || @show key
    @test mse_not_regressed
end
