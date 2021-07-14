using Plots
using OrderedCollections
using Test
using NCDatasets
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
#! format: on

include("variable_map.jl")

function compute_mse(
    ds,
    ds_pycles,
    experiment,
    best_mse,
    foldername;
    plot_comparison=true,
)
    mse = Dict()
    time_tcc = ds.group["timeseries"]["t"][:]
    time_les = ds_pycles["t"][:]
    n_time_points = length(time_tcc)

    # Ensure domain matches:
    z_les = ds_pycles["z_half"][:]
    z_tcc = ds.group["profiles"]["z_half"][:]
    n_grid_points = length(z_tcc)
    @info "Z extent for LES vs CLIMA:"
    @show extrema(z_tcc)
    @show extrema(z_les)

    # Find the nearest matching final time:
    t_cmp = min(time_tcc[end], time_les[end])

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
    pycles_variables = []
    data_scales = []
    pycles_weight = []
    for tc_var in keys(best_mse)

        # Only compare fields defined for var_map
        les_var = var_map(tc_var)
        les_var == nothing && continue
        les_var isa String || continue

        push!(tcc_variables, tc_var)
        push!(pycles_variables, les_var)
        ds_les = ds_pycles.group["profiles"]
        data_les = ds_les[les_var][:]

        # Scale the data for comparison
        push!(pycles_weight, "1")

        # Interpolate data
        steady_data = length(size(data_les)) == 1
        if steady_data
            data_les_cont = Spline1D(z_les, data_les)
        else # unsteady data
            data_les_cont = Spline2D(time_les, z_les, data_les')
        end

        if haskey(ds.group["profiles"], tc_var)
            data_tcc_arr_ = ds.group["profiles"][tc_var][:]
        elseif haskey(ds.group["reference"], tc_var)
            data_tcc_arr_ = ds.group["reference"][tc_var][:]
        else
            error("No key for $tc_var found in the nc file.")
        end
        data_tcc_arr = reshape(data_tcc_arr_, (n_time_points, length(z_tcc)))
        data_tcc[tc_var] = Spline2D(time_tcc, z_tcc, data_tcc_arr)

        # Compute data scale
        data_scale = sum(abs.(data_les)) / length(data_les)
        push!(data_scales, data_scale)

        # Plot comparison
        if plot_comparison
            p = plot()
            if steady_data
                data_les_cont_mapped = map(z -> data_les_cont(z), z_tcc)
            else
                data_les_cont_mapped = map(z -> data_les_cont(t_cmp, z), z_tcc)
            end
            plot!(
                data_les_cont_mapped,
                z_tcc ./ 10^3,
                xlabel = tc_var,
                ylabel = "z [km]",
                label = "PyCLES",
            )
            plot!(
                map(z -> data_tcc[tc_var](t_cmp, z), z_tcc),
                z_tcc ./ 10^3,
                xlabel = tc_var,
                ylabel = "z [km]",
                label = "TC.jl",
            )
            mkpath(foldername)
            savefig(joinpath(foldername, "$tc_var.png"))
        end

        # Compute mean squared error (mse)
        if steady_data
            mse_single_var = sum(map(z_tcc) do z
                (data_les_cont(z) - data_tcc[tc_var](t_cmp, z))^2
            end)
        else
            mse_single_var = sum(map(z_tcc) do z
                (data_les_cont(t_cmp, z) - data_tcc[tc_var](t_cmp, z))^2
            end)
        end
        # Normalize by data scale
        mse[tc_var] = mse_single_var / data_scale^2

        push!(mse_reductions, (best_mse[tc_var] - mse[tc_var]) / best_mse[tc_var] * 100)
        push!(computed_mse, mse[tc_var])
        push!(table_best_mse, best_mse[tc_var])
    end

    # Tabulate output
    header = [
        "Variable" "Variable" "Weight" "Data scale" "MSE" "MSE" "MSE"
        "TC.jl (EDMF)" "PyCLES" "PyCLES" "" "Computed" "Best" "Reduction (%)"
    ]
    table_data = hcat(
        tcc_variables,
        pycles_variables,
        pycles_weight,
        data_scales,
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
        (data, i, j) -> !sufficient_mse(data[i, 5], data[i, 6]) && j == 5,
        crayon"red bold",
    )
    hl_worsened_mse_reduction = Highlighter(
        (data, i, j) -> !sufficient_mse(data[i, 5], data[i, 6]) && j == 7,
        crayon"red bold",
    )
    hl_improved_mse = Highlighter(
        (data, i, j) -> sufficient_mse(data[i, 5], data[i, 6]) && j == 7,
        crayon"green bold",
    )
    pretty_table(
        table_data,
        header,
        formatters = ft_printf("%.16e", 5:6),
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
