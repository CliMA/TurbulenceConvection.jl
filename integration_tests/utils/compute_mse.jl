import Plots
using OrderedCollections
using Test
using NCDatasets
import StatsBase
using Dierckx
using PrettyTables
using Printf
using ArtifactWrappers

ENV["GKSwstype"] = "nul"

# Get PyCLES_output dataset folder:
#! format: off
PyCLES_output_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "PyCLES_output",
    ArtifactFile[
    ArtifactFile(url = "https://caltech.box.com/shared/static/johlutwhohvr66wn38cdo7a6rluvz708.nc", filename = "Rico.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/zraeiftuzlgmykzhppqwrym2upqsiwyb.nc", filename = "Gabls.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/toyvhbwmow3nz5bfa145m5fmcb2qbfuz.nc", filename = "DYCOMS_RF01.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/ivo4751camlph6u3k68ftmb1dl4z7uox.nc", filename = "TRMM_LBA.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/4osqp0jpt4cny8fq2ukimgfnyi787vsy.nc", filename = "ARM_SGP.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/jci8l11qetlioab4cxf5myr1r492prk6.nc", filename = "Bomex.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/pzuu6ii99by2s356ij69v5cb615200jq.nc", filename = "Soares.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/7upt639siyc2umon8gs6qsjiqavof5cq.nc", filename = "Nieuwstadt.nc",),
    ],
)
PyCLES_output_dataset_path = get_data_folder(PyCLES_output_dataset)

SCAMPy_output_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "SCAMPy_output",
    ArtifactFile[
    ArtifactFile(url = "https://caltech.box.com/shared/static/1dzpydqiagjvzfpyv9lbic3atvca93hl.nc", filename = "Rico.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/agkycoum93bd6xjduyaeo3oqg6asoru5.nc", filename = "GABLS.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/fqpq1q74uxfh1e8018hwhqogw1htn2lq.nc", filename = "DYCOMS_RF01.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/wevi0rqiwo6sgkqdhcddr72u5ylt0tqp.nc", filename = "TRMM_LBA.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/fis6n0g9x9lts70m0zmve5ullqnw0pzq.nc", filename = "ARM_SGP.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/t6qq6plt2oxcmmy40r1szgokahykqggp.nc", filename = "Bomex.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/wp8k4m7ta1hs0c6e4j2fpsp3kj05wdip.nc", filename = "Soares.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/vbuxzg85scwy9mg0ziiuzy6fimo0e389.nc", filename = "Nieuwstadt.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/72t6fr1gq10tg3jjputtp35nfzex0o4k.nc", filename = "DryBubble.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/7axeussneeg8g3k0ndvagsn0pkmbij3e.nc", filename = "life_cycle_Tan2018.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/r6t7dk6g35bmbvc86h006yusb6rd8vkw.nc", filename = "SP.nc",),
    ],
)
SCAMPy_output_dataset_path = get_data_folder(SCAMPy_output_dataset)
#! format: on

include("variable_map.jl")

function get_data(ds, var)
    if haskey(ds, var)
        return ds[var][:]
    elseif haskey(ds.group["profiles"], var)
        return ds.group["profiles"][var][:]
    elseif haskey(ds.group["reference"], var)
        return ds.group["reference"][var][:]
    end
    error("No key for $var found in the nc file.")
end

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
        dict[sym] = Dataset(filename, "r")
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
    return pop!(sort(matching_paths))
end

function compute_mse_wrapper(
    case_name,
    best_mse,
    ds_tc_filename,
    ds_tc_main_filename = nothing;
    ds_les_filename = joinpath(PyCLES_output_dataset_path, "$case_name.nc"),
    ds_scm_filename = joinpath(SCAMPy_output_dataset_path, "$case_name.nc"),
    plot_dir = joinpath(dirname(ds_tc_filename), "comparison"),
    kwargs...,
)
    # Note: cluster_data_prefix is also defined in pipeline.yml
    # if haskey(ENV, "BUILDKITE_BUILD_PATH") && haskey(ENV, "BUILDKITE_PIPELINE_SLUG")
    #     bkp = ENV["BUILDKITE_BUILD_PATH"]
    #     bks = ENV["BUILDKITE_PIPELINE_SLUG"]
    #     cluster_data_prefix = joinpath(bkp, bks, "main_results")
    #     path = find_latest_dataset_folder(; dir = cluster_data_prefix)

    #     # TODO: make this more robust in case folder/file changes
    #     folder_name = joinpath("Output.$case_name.01", "stats")
    #     ds_tc_main_filename = joinpath(path, folder_name, "Stats.$case_name.nc")
    # end
    # Note that we're using a closure over
    #  - PyCLES_output_dataset_path
    #  - SCAMPy_output_dataset_path
    all_args = (case_name, best_mse, plot_dir)

    ds_dict = append_dict(ds_les_filename, :ds_pycles)
    ds_dict = append_dict(ds_scm_filename, :ds_scampy, ds_dict)
    ds_dict = append_dict(ds_tc_filename, :ds_tc, ds_dict)
    ds_dict = append_dict(ds_tc_main_filename, :ds_tc_main, ds_dict)

    try
        computed_mse = compute_mse(all_args...; ds_dict, kwargs...)
    finally
        for ds in values(ds_dict)
            close(ds)
        end
    end
end


function compute_mse(experiment, best_mse, plot_dir; ds_dict, plot_comparison = true, t_start, t_stop)

    ds_scampy = haskey(ds_dict, :ds_scampy) ? ds_dict[:ds_scampy] : nothing
    ds_tc_main = haskey(ds_dict, :ds_tc_main) ? ds_dict[:ds_tc_main] : nothing
    ds_pycles = haskey(ds_dict, :ds_pycles) ? ds_dict[:ds_pycles] : nothing
    ds_tc = haskey(ds_dict, :ds_tc) ? ds_dict[:ds_tc] : nothing
    # Make pycles optional
    have_pycles_ds = ds_pycles ≠ nothing
    have_tc_main = ds_tc_main ≠ nothing
    if !have_pycles_ds
        ds_pycles = ds_scampy
    end
    if !have_tc_main
        ds_tc_main = ds_tc
    end
    mse = Dict()
    time_tcc = get_time(ds_tc, "t")
    time_les = get_time(ds_pycles, "t")
    time_tcm = get_time(ds_tc_main, "t")
    time_scm = get_time(ds_scampy, "t")
    n_time_points = length(time_tcc)

    mkpath(plot_dir)
    # Ensure domain matches:
    z_les = get_data(ds_pycles, "z_half")
    z_tcc = get_data(ds_tc, "z_half")
    z_tcm = get_data(ds_tc_main, "z_half")
    z_scm = get_data(ds_scampy, "z_half")
    n_grid_points = length(z_tcc)
    @info "z extrema (les,scm,tcm,tcc): $(extrema(z_les)), $(extrema(z_scm)), $(extrema(z_tcm)), $(extrema(z_tcc))"
    @info "n-grid points (les,scm,tcm,tcc): $(length(z_les)), $(length(z_scm)), $(length(z_tcm)), $(length(z_tcc))"

    @info "time extrema (les,scm,tcm,tcc): $(extrema(time_les)), $(extrema(time_scm)), $(extrema(time_tcm)), $(extrema(time_tcc))"
    @info "n-time points (les,scm,tcm,tcc): $(length(time_les)), $(length(time_scm)), $(length(time_tcm)), $(length(time_tcc))"

    # Find the nearest matching final time:
    t_cmp = min(time_tcc[end], time_tcm[end], time_les[end], time_scm[end], t_stop)
    @info "time compared: $t_cmp"

    # Accidentally running a short simulation
    # could improve MSE. So, let's test that
    # we run for at least t_compare. We should
    # increase this as we can reach higher CFL.
    # @test t_cmp >= t_compare

    # Ensure z_tcc and fields are consistent lengths:
    @test length(z_tcc) == n_grid_points
    tcc_variables = []
    computed_mse = []
    table_best_mse = []
    mse_reductions = []
    pycles_variables = []
    data_scales_scm = []
    data_scales_les = []
    data_scales_tcc = []
    data_scales_tcm = []
    pycles_weight = []
    for tc_var in keys(best_mse)

        # Only compare fields defined for var_map_les
        if have_pycles_ds
            les_var = var_map_les(tc_var)
        else
            les_var = var_map_scampy(tc_var)
        end
        scm_var = var_map_scampy(tc_var)
        les_var == nothing && continue
        les_var isa String || continue
        @info "Assembling plots for $tc_var"
        push!(tcc_variables, tc_var)
        push!(pycles_variables, les_var)

        data_les_arr = get_data(ds_pycles, les_var)'
        data_tcm_arr = get_data(ds_tc_main, tc_var)'
        data_tcc_arr = get_data(ds_tc, tc_var)'
        data_scm_arr = get_data(ds_scampy, scm_var)'
        @info "     Data sizes (les,scm,tcc,tcm): $(size(data_les_arr)), $(size(data_scm_arr)), $(size(data_tcc_arr)), $(size(data_tcm_arr))"
        # Scale the data for comparison
        push!(pycles_weight, "1")

        # Interpolate data
        steady_data = length(size(data_les_arr)) == 1
        if steady_data
            data_les_cont = Spline1D(z_les, data_les_arr)
            data_tcc_cont = Spline1D(z_tcc, data_tcc_arr)
            data_tcm_cont = Spline1D(z_tcm, data_tcm_arr)
            data_scm_cont = Spline1D(z_scm, data_scm_arr)
            data_les_cont_mapped = map(z -> data_les_cont(z), z_tcc)
            data_tcm_cont_mapped = map(z -> data_tcm_cont(z), z_tcc)
            data_tcc_cont_mapped = map(z -> data_tcc_cont(z), z_tcc)
            data_scm_cont_mapped = map(z -> data_scm_cont(z), z_tcc)
        else # unsteady data
            data_les_cont = Spline2D(time_les, z_les, data_les_arr)
            data_tcm_cont = Spline2D(time_tcm, z_tcm, data_tcm_arr)
            data_tcc_cont = Spline2D(time_tcc, z_tcc, data_tcc_arr)
            data_scm_cont = Spline2D(time_scm, z_scm, data_scm_arr)
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

        # Compute data scale
        data_scale_tcm = sum(abs.(data_tcm_arr)) / length(data_tcm_arr)
        data_scale_tcc = sum(abs.(data_tcc_arr)) / length(data_tcc_arr)
        data_scale_scm = sum(abs.(data_scm_arr)) / length(data_scm_arr)
        data_scale_les = sum(abs.(data_les_arr)) / length(data_les_arr)
        push!(data_scales_tcm, data_scale_tcm)
        push!(data_scales_tcc, data_scale_tcc)
        push!(data_scales_scm, data_scale_scm)
        push!(data_scales_les, data_scale_les)

        # Plot comparison
        if plot_comparison
            p = Plots.plot()
            if have_pycles_ds
                Plots.plot!(data_les_cont_mapped, z_tcc ./ 10^3, xlabel = tc_var, ylabel = "z [km]", label = "PyCLES")
            end
            Plots.plot!(data_scm_cont_mapped, z_tcc ./ 10^3, xlabel = tc_var, ylabel = "z [km]", label = "SCAMPy")
            if have_tc_main
                Plots.plot!(
                    data_tcm_cont_mapped,
                    z_tcc ./ 10^3,
                    xlabel = tc_var,
                    ylabel = "z [km]",
                    label = "TC.jl (main)",
                )
            end
            Plots.plot!(data_tcc_cont_mapped, z_tcc ./ 10^3, xlabel = tc_var, ylabel = "z [km]", label = "TC.jl")
            @info "     Saving $(joinpath(plot_dir, "$tc_var.png"))"
            Plots.savefig(joinpath(plot_dir, "profile_$tc_var.png"))

            if tc_var == "updraft_thetal"
                # TODO: remove this if-else when artifacts are updated
                #       So that updraft_thetal is initialized before IO
                clims_min = minimum(data_tcc_arr)
                clims_max = maximum(data_tcc_arr)
                clims = (clims_min, clims_max)
            else
                clims_min = min(minimum(data_scm_arr), minimum(data_tcc_arr))
                clims_max = max(maximum(data_scm_arr), maximum(data_tcc_arr))
                clims = (clims_min, clims_max)
            end

            p1 = Plots.contourf(
                time_scm,
                z_scm,
                data_scm_arr';
                xlabel = tc_var,
                ylabel = "height (m)",
                c = :viridis,
                clims = clims,
                title = "SCAMPy",
            )
            if have_tc_main
                p2 = Plots.contourf(
                    time_tcc,
                    z_tcc,
                    data_tcm_arr';
                    xlabel = tc_var,
                    ylabel = "height (m)",
                    c = :viridis,
                    clims = clims,
                    title = "TC.jl (main)",
                )
            end
            p3 = Plots.contourf(
                time_tcc,
                z_tcc,
                data_tcc_arr';
                xlabel = tc_var,
                ylabel = "height (m)",
                c = :viridis,
                clims = clims,
                title = have_tc_main ? "TC.jl (PR)" : "TC.jl",
            )
            if have_tc_main
                Plots.plot(p1, p2, p3; layout = (3, 1), size = (400, 600))
            else
                Plots.plot(p1, p3; layout = (2, 1))
            end
            Plots.savefig(joinpath(plot_dir, "contours_" * tc_var * ".png"))
        end

        # Compute mean squared error (mse)
        mse_single_var = sum((data_les_cont_mapped .- data_tcc_cont_mapped) .^ 2)
        # Normalize by data scale
        mse[tc_var] = mse_single_var / data_scale_tcc^2

        push!(mse_reductions, (best_mse[tc_var] - mse[tc_var]) / best_mse[tc_var] * 100)
        push!(computed_mse, mse[tc_var])
        push!(table_best_mse, best_mse[tc_var])
    end

    # Tabulate output
    header = [
        "Variable" "Variable" "Weight" "Data scale" "Data scale" "MSE" "MSE" "MSE"
        "TC.jl (EDMF)" "PyCLES" "PyCLES" "tcc" "scm" "Computed" "Best" "Reduction (%)"
    ]
    table_data = hcat(
        tcc_variables,
        pycles_variables,
        pycles_weight,
        data_scales_tcc,
        data_scales_scm,
        computed_mse,
        table_best_mse,
        mse_reductions,
    )

    @info @sprintf("Experiment comparison: %s at time t=%s\n", experiment, t_cmp)
    hl_worsened_mse = Highlighter((data, i, j) -> !sufficient_mse(data[i, 6], data[i, 7]) && j == 6, crayon"red bold")
    hl_worsened_mse_reduction =
        Highlighter((data, i, j) -> !sufficient_mse(data[i, 6], data[i, 7]) && j == 8, crayon"red bold")
    hl_improved_mse = Highlighter((data, i, j) -> sufficient_mse(data[i, 6], data[i, 7]) && j == 8, crayon"green bold")
    pretty_table(
        table_data,
        header,
        formatters = ft_printf("%.16e", 6:7),
        header_crayon = crayon"yellow bold",
        subheader_crayon = crayon"green bold",
        highlighters = (hl_worsened_mse, hl_improved_mse, hl_worsened_mse_reduction),
        crop = :none,
    )

    return mse
end

sufficient_mse(computed_mse, best_mse) = computed_mse <= best_mse + sqrt(eps())

function test_mse(computed_mse, best_mse, key)
    mse_not_regressed = sufficient_mse(computed_mse[key], best_mse[key])
    mse_not_regressed || @show key
    @test mse_not_regressed
end
