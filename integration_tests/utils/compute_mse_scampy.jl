using Plots
using OrderedCollections
using Test
using NCDatasets
using Dierckx
using PrettyTables
using Printf
using ArtifactWrappers

ENV["GKSwstype"] = "nul"

# Get scampy_output dataset folder:
#! format: off
scampy_output_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "scampy_output",
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
scampy_output_dataset_path = get_data_folder(scampy_output_dataset)
#! format: on

include("variable_map.jl")

function compute_mse_scampy(
    ds,
    ds_scampy,
    experiment,
    best_mse,
    foldername;
    plot_comparison=true,
    interpolate_in_time = false
)
    mse = Dict()
    time_tcc = ds.group["timeseries"]["t"][:]
    time_scm = ds_scampy.group["timeseries"]["t"][:]
    n_time_points = length(time_tcc)
    n_time_points = length(time_scm)

    # Ensure domain matches:
    z_scm = ds_scampy.group["profiles"]["z_half"][:]
    z_tcc = ds.group["profiles"]["z_half"][:]
    n_grid_points = length(z_tcc)
    @info "Z extent for LES vs CLIMA:"
    @show extrema(z_tcc)
    @show extrema(z_scm)

    # Find the nearest matching final time:
    t_cmp = min(time_tcc[end], time_scm[end])

    # Accidentally running a short simulation
    # could improve MSE. So, let's test that
    # we run for at least t_compare. We should
    # increase this as we can reach higher CFL.
    # @test t_cmp >= t_compare

    # Ensure z_tcc and fields are consistent lengths:
    @test length(z_tcc) == n_grid_points

    data_tcc = Dict()
    dons_cont = Dict()
    tcc_variables = []
    computed_mse = []
    table_best_mse = []
    mse_reductions = []
    scampy_variables = []
    data_scales_tcc = []
    data_scales_scm = []
    scampy_weight = []
    for tc_var in keys(best_mse)

        # Only compare fields defined for var_map
        # les_var = var_map(tc_var)
        les_var = tc_var
        les_var == nothing && continue
        les_var isa String || continue

        push!(tcc_variables, tc_var)
        push!(scampy_variables, les_var)
        if haskey(ds_scampy.group["profiles"], les_var)
            data_scm = ds_scampy.group["profiles"][les_var][:]
        elseif haskey(ds_scampy.group["reference"], les_var)
            data_scm = ds_scampy.group["reference"][les_var][:]
        else
            error("No key for $tc_var found (les) in the nc file.")
        end
        # Scale the data for comparison
        push!(scampy_weight, "1")

        # Interpolate data
        steady_data = length(size(data_scm)) == 1
        if interpolate_in_time
            if steady_data
                data_scm_cont = Spline1D(z_scm, data_scm)
            else # unsteady data
                data_scm_cont = Spline2D(time_scm, z_scm, data_scm')
            end
        else
            data_scm_cont = data_scm
        end

        if haskey(ds.group["profiles"], tc_var)
            data_tcc_arr_ = ds.group["profiles"][tc_var][:]
        elseif haskey(ds.group["reference"], tc_var)
            data_tcc_arr_ = ds.group["reference"][tc_var][:]
        else
            error("No key for $tc_var found in the nc file.")
        end
        if interpolate_in_time
            data_tcc_arr = reshape(data_tcc_arr_, (n_time_points, length(z_tcc)))
            data_tcc[tc_var] = Spline2D(time_tcc, z_tcc, data_tcc_arr)
        else
            data_tcc_arr = reshape(data_tcc_arr_, (n_time_points, length(z_tcc)))
            # @show size(data_tcc_arr)
            data_tcc[tc_var] = data_tcc_arr
        end

        # Compute data scale
        data_scale_scm = sum(abs.(data_scm)) / length(data_scm)
        data_scale_tcc = sum(abs.(data_tcc[tc_var])) / length(data_tcc[tc_var])
        push!(data_scales_scm, data_scale_scm)
        push!(data_scales_tcc, data_scale_tcc)

        # Plot comparison
        if plot_comparison
            p = plot()
            if interpolate_in_time
                if steady_data
                    data_scm_cont_mapped = map(z -> data_scm_cont(z), z_tcc)
                else
                    data_scm_cont_mapped = map(z -> data_scm_cont(t_cmp, z), z_tcc)
                end
            else
                data_scm_cont_mapped = data_scm_cont
            end
            if interpolate_in_time
                data_tcc_cont_mapped = map(z -> data_tcc[tc_var](t_cmp, z), z_tcc)
            else
                data_tcc_cont_mapped = reshape(data_tcc[tc_var], length(data_tcc[tc_var]), 1)
            end
            # @show size(data_scm_cont_mapped)
            # @show size(data_tcc_cont_mapped)
            # @show size(z_tcc)
            plot!(
                data_scm_cont_mapped,
                z_tcc ./ 10^3,
                xlabel = tc_var,
                ylabel = "z [km]",
                label = "scampy",
            )
            plot!(
                data_tcc_cont_mapped,
                z_tcc ./ 10^3,
                xlabel = tc_var,
                ylabel = "z [km]",
                label = "TC.jl",
            )
            mkpath(foldername)
            savefig(joinpath(foldername, "$tc_var.png"))
        end

        # Compute mean squared error (mse)
        if interpolate_in_time
            if steady_data
                mse_single_var = sum(map(z_tcc) do z
                    (data_scm_cont(z) - data_tcc[tc_var](t_cmp, z))^2
                end)
            else
                mse_single_var = sum(map(z_tcc) do z
                    (data_scm_cont(t_cmp, z) - data_tcc[tc_var](t_cmp, z))^2
                end)
            end
        else
            mse_single_var = sum((data_scm_cont_mapped .- data_tcc_cont_mapped).^2)
        end
        # Normalize by data scale
        mse[tc_var] = mse_single_var / data_scale_scm^2

        push!(mse_reductions, (best_mse[tc_var] - mse[tc_var]) / best_mse[tc_var] * 100)
        push!(computed_mse, mse[tc_var])
        push!(table_best_mse, best_mse[tc_var])
    end

    # Tabulate output
    header = [
        "Variable" "Variable" "Weight" "Data scale" "Data scale" "MSE" "MSE" "MSE"
        "TC.jl (EDMF)" "scampy" "scampy" "tcc" "scm" "Computed" "Best" "Reduction (%)"
    ]
    table_data = hcat(
        tcc_variables,
        scampy_variables,
        scampy_weight,
        data_scales_tcc,
        data_scales_scm,
        computed_mse,
        table_best_mse,
        mse_reductions,
    )

    @info @sprintf(
        "Experiment comparison: %s at time t=%s\n",
        experiment,
        t_cmp
    )
    hl_worsened_mse = Highlighter(
        (data, i, j) -> !sufficient_mse(data[i, 6], data[i, 7]) && j == 6,
        crayon"red bold",
    )
    hl_worsened_mse_reduction = Highlighter(
        (data, i, j) -> !sufficient_mse(data[i, 6], data[i, 7]) && j == 8,
        crayon"red bold",
    )
    hl_improved_mse = Highlighter(
        (data, i, j) -> sufficient_mse(data[i, 6], data[i, 7]) && j == 8,
        crayon"green bold",
    )
    pretty_table(
        table_data,
        header,
        formatters = ft_printf("%.16e", 6:7),
        header_crayon = crayon"yellow bold",
        subheader_crayon = crayon"green bold",
        highlighters = (
            hl_worsened_mse,
            hl_improved_mse,
            hl_worsened_mse_reduction,
        ),
        crop = :none,
    )

    return mse
end

sufficient_mse(computed_mse, best_mse) = computed_mse <= best_mse + sqrt(eps())

function test_mse(computed_mse, best_mse, key)
    mse_not_regressed = sufficient_mse(computed_mse[key], best_mse[key])
    @test mse_not_regressed
    mse_not_regressed || @show key
end
