import Plots
using OrderedCollections
using Test
using NCDatasets
using Dierckx
using PrettyTables
using Printf
using Statistics
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
    ArtifactFile(url = "https://caltech.box.com/shared/static/7upt639siyc2umon8gs6qsjiqavof5cq.nc", filename = "Nieuwstadt.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/72t6fr1gq10tg3jjputtp35nfzex0o4k.nc", filename = "DryBubble.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/7axeussneeg8g3k0ndvagsn0pkmbij3e.nc", filename = "life_cycle_Tan2018.nc",),
    ArtifactFile(url = "https://caltech.box.com/shared/static/r6t7dk6g35bmbvc86h006yusb6rd8vkw.nc", filename = "SP.nc",),
    ],
)
SCAMPy_output_dataset_path = get_data_folder(SCAMPy_output_dataset)
#! format: on

include("variable_map.jl")

function has_data(ds, var)
    if haskey(ds.group["profiles"], var)
        return true
    elseif haskey(ds.group["reference"], var)
        return true
    else
        return false
    end
end

function get_data(ds, var)
    if haskey(ds.group["profiles"], var)
        return ds.group["profiles"][var][:]
    elseif haskey(ds.group["reference"], var)
        return ds.group["reference"][var][:]
    else
        error("No key for $var found in the nc file.")
    end
end

function compute_mse(
    experiment,
    best_mse,
    foldername;
    ds_scampy=nothing,
    ds_pycles=nothing,
    ds_turb_conv=nothing,
    plot_comparison=true,
    skip_contours = true,
)
    mse = Dict()
    time_tcc = ds_turb_conv.group["timeseries"]["t"][:]
    # time_les = ds_pycles["t"][:]
    time_scm = ds_scampy.group["timeseries"]["t"][:]
    n_time_points = length(time_tcc)

    mkpath(foldername)
    # Ensure domain matches:
    # z_les = ds_pycles["z_half"][:]
    z_tcc = ds_turb_conv.group["profiles"]["z_half"][:]
    z_scm = ds_scampy.group["profiles"]["z_half"][:]
    n_grid_points = length(z_tcc)
    @info "Z extent for LES vs CLIMA:"
    @show extrema(z_tcc)
    # @show extrema(z_les)
    @show extrema(z_scm)
    @info "n-grid points"
    @show length(z_tcc)
    # @show length(z_les)
    @show length(z_scm)

    # Find the nearest matching final time:
    t_cmp = min(time_tcc[end], time_scm[end])

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
    # data_scales_les = []
    data_scales_tcc = []
    pycles_weight = []
    max_mse = 0
    mismatch_found = false

    all_keys = unique([keys(ds_scampy.group["profiles"])..., keys(best_mse)...])
    # for tc_var in keys(best_mse)
    # @show keys(ds_scampy.group["profiles"])
    for tc_var in all_keys

        # Only compare fields defined for var_map_les
        # les_var = var_map_les(tc_var)
        scm_var = var_map_scampy(tc_var)
        tc_var == "z" && continue
        tc_var == "z_half" && continue
        tc_var == "t" && continue
        tc_var == "rho0" && continue
        tc_var == "rho0_half" && continue
        tc_var == "p0" && continue
        tc_var == "p0_half" && continue
        tc_var == "alpha0" && continue
        tc_var == "alpha0_half" && continue
        # les_var == nothing && continue
        # les_var isa String || continue

        # push!(pycles_variables, les_var)

        # data_les_arr = get_data(ds_pycles, les_var)'
        data_tcc_arr = has_data(ds_turb_conv, tc_var)
        data_scm_arr = has_data(ds_scampy, scm_var)
        if !(data_tcc_arr && data_scm_arr)
            continue
        end
        data_tcc_arr = get_data(ds_turb_conv, tc_var)'
        data_scm_arr = get_data(ds_scampy, scm_var)'

        push!(tcc_variables, tc_var)

        # Scale the data for comparison
        # push!(pycles_weight, "1")

        # Interpolate data
        steady_data = size(data_tcc_arr,1) == 1
        # @show tc_var, steady_data, size(data_tcc_arr), size(data_scm_arr)
        # if steady_data
        #     # data_les_cont = Spline1D(z_les, data_les_arr)
        #     data_tcc_cont = Spline1D(z_tcc, data_tcc_arr)
        #     data_scm_cont = Spline1D(z_scm, data_scm_arr)
        #     data_les_cont_mapped = map(z -> data_les_cont(z), z_tcc)
        #     data_tcc_cont_mapped = map(z -> data_tcc_cont(z), z_tcc)
        #     data_scm_cont_mapped = map(z -> data_scm_cont(z), z_tcc)
        # else # unsteady data
        #     # data_les_cont = Spline2D(time_les, z_les, data_les_arr)
        #     data_tcc_cont = Spline2D(time_tcc, z_tcc, data_tcc_arr)
        #     data_scm_cont = Spline2D(time_scm, z_scm, data_scm_arr)
        #     # data_les_cont_mapped = map(z -> data_les_cont(t_cmp, z), z_tcc)
        #     data_tcc_cont_mapped = map(z -> data_tcc_cont(t_cmp, z), z_tcc)
        #     data_scm_cont_mapped = map(z -> data_scm_cont(t_cmp, z), z_tcc)
        # end

        # if steady_data
        #     # data_les_cont = Spline1D(z_les, data_les_arr)
        #     data_tcc_cont_mapped = data_tcc_arr[:]
        #     data_scm_cont_mapped = data_scm_arr[:]
        # else # unsteady data
        #     data_tcc_cont_mapped = data_tcc_arr[:,end]
        #     data_scm_cont_mapped = data_scm_arr[:,end]
        # end

        data_tcc_cont_mapped = data_tcc_arr[:,:]
        data_scm_cont_mapped = data_scm_arr[:,:]

        profile_data_tcc_cont_mapped = Statistics.mean(data_tcc_arr[:,:]; dims=1)
        profile_data_scm_cont_mapped = Statistics.mean(data_scm_arr[:,:]; dims=1)
        profile_data_tcc_cont_mapped = reshape(profile_data_tcc_cont_mapped, length(z_tcc))
        profile_data_scm_cont_mapped = reshape(profile_data_scm_cont_mapped, length(z_scm))

        # Compute data scale
        data_scale_tcc = sum(abs.(data_tcc_arr)) / length(data_tcc_arr)
        data_scale_scm = sum(abs.(data_scm_arr)) / length(data_scm_arr)
        # data_scale_les = sum(abs.(data_les_arr)) / length(data_les_arr)
        push!(data_scales_tcc, data_scale_tcc)
        push!(data_scales_scm, data_scale_scm)
        # push!(data_scales_les, data_scale_les)
        # Plot comparison
        if plot_comparison
            p = Plots.plot()
            # plot!(
            #     data_les_cont_mapped,
            #     z_tcc ./ 10^3,
            #     xlabel = tc_var,
            #     ylabel = "z [km]",
            #     label = "PyCLES",
            # )
            Plots.plot!(
                profile_data_scm_cont_mapped,
                z_tcc ./ 10^3,
                xlabel = tc_var,
                ylabel = "z [km]",
                label = "SCAMPy",
                markershape = :circle,
                markersize = 6,
            )
            Plots.plot!(
                profile_data_tcc_cont_mapped,
                z_tcc ./ 10^3,
                xlabel = tc_var,
                ylabel = "z [km]",
                label = "TC.jl",
                markershape = :diamond,
            )
            @info "Saving $(joinpath(foldername, "$tc_var.png"))"
            Plots.savefig(joinpath(foldername, "profiles_$tc_var.png"))

            if !skip_contours
                clims_min = min(minimum(data_tcc_arr), minimum(data_tcc_arr))
                clims_max = max(maximum(data_tcc_arr), maximum(data_tcc_arr))
                clims = (clims_min, clims_max)
                p1 = Plots.contourf(
                    time_scm, z_scm, data_scm_arr';
                    xlabel = tc_var,
                    ylabel = "height (m)",
                    c = :viridis,
                    clims = clims,
                    title="SCAMPy",
                )

                p2 = Plots.contourf(
                    time_tcc, z_tcc, data_tcc_arr';
                    xlabel = tc_var,
                    ylabel = "height (m)",
                    c = :viridis,
                    clims = clims,
                    title="TC.jl",
                )
                Plots.plot(p1, p2; layout= (2, 1))

                Plots.savefig(joinpath(foldername,"contours_"*tc_var*".png"))
            end
        end

        # Compute mean squared error (mse)
        # if steady_data
        #     mse_single_var = sum(map(z_tcc) do z
        #         (data_scm_cont(z) - data_tcc_cont(t_cmp, z))^2
        #     end)
        # else
        mse_single_var = sum((data_scm_cont_mapped .- data_tcc_cont_mapped).^2)
        # end
        # Normalize by data scale
        mse[tc_var] = mse_single_var / data_scale_tcc^2
        if !(data_scale_tcc ≈ 0 && data_scale_scm ≈ 0)
            if !isnan(mse[tc_var])
                max_mse = max(max_mse, mse[tc_var])
            elseif (isnan(data_scale_tcc) && !isnan(data_scale_scm)) ||
                   (isnan(data_scale_scm) && !isnan(data_scale_tcc))
                mismatch_found = true
            elseif ((data_scale_tcc ≈ 0) && !(data_scale_scm ≈ 0)) ||
                   ((data_scale_scm ≈ 0) && !(data_scale_tcc ≈ 0))
                mismatch_found = true
            end
        end

        push!(mse_reductions, (best_mse[tc_var] - mse[tc_var]) / best_mse[tc_var] * 100)
        push!(computed_mse, mse[tc_var])
        push!(table_best_mse, best_mse[tc_var])
    end

    # Tabulate output
    header = [
        "Variable" "Data scale" "Data scale" "MSE" "MSE" "MSE"
        "" "tcc" "scm" "Computed" "Best" "Reduction (%)"
    ]

    @show size(tcc_variables)
    @show size(data_scales_tcc)
    @show size(data_scales_scm)
    @show size(computed_mse)
    @show size(table_best_mse)
    @show size(mse_reductions)

    table_data = hcat(
        tcc_variables,
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
        (data, i, j) -> !sufficient_mse(data[i, 4], data[i, 5]) && j == 4,
        crayon"red bold",
    )
    hl_worsened_mse_reduction = Highlighter(
        (data, i, j) -> !sufficient_mse(data[i, 4], data[i, 5]) && j == 6,
        crayon"red bold",
    )
    hl_improved_mse = Highlighter(
        (data, i, j) -> sufficient_mse(data[i, 4], data[i, 5]) && j == 6,
        crayon"green bold",
    )
    pretty_table(
        table_data,
        header,
        formatters = ft_printf("%.16e", 4:5),
        header_crayon = crayon"yellow bold",
        subheader_crayon = crayon"green bold",
        highlighters = (
            hl_worsened_mse,
            hl_improved_mse,
            hl_worsened_mse_reduction,
        ),
        crop = :none,
    )
    @info "max_mse = $max_mse"
    if mismatch_found
        @warn "Mismatch found (one sim has NaNs, while the other does not)"
    end

    return mse
end

sufficient_mse(computed_mse, best_mse) = computed_mse <= best_mse + sqrt(eps())

function test_mse(computed_mse, best_mse, key)
    mse_not_regressed = sufficient_mse(computed_mse[key], best_mse[key])
    @test mse_not_regressed
    mse_not_regressed || @show key
end
