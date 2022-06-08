import Plots

space_string(::CC.Spaces.FaceExtrudedFiniteDifferenceSpace) = "(Face field)"
space_string(::CC.Spaces.CenterExtrudedFiniteDifferenceSpace) = "(Center field)"

function process_name(s::AbstractString)
    # "c_ρ", "c_ρe", "c_uₕ_1", "c_uₕ_2", "f_w_1"
    s = replace(s, "components_data_" => "")
    s = replace(s, "ₕ" => "_h")
    s = replace(s, "ρ" => "rho")
    return s
end
processed_varname(pc::Tuple) = process_name(join(pc, "_"))

function plot_profiles(Y, output_dir; filter_prop_chain = pc -> true)
    # Column animations
    mkpath(output_dir)
    n_columns = TC.number_of_columns(Y)
    @info "Creating profiles with `n_columns` = $n_columns"
    prop_chains = filter(filter_prop_chain, CC.Fields.property_chains(Y))
    for prop_chain in prop_chains
        var_name = processed_varname(prop_chain)
        var_space = axes(CC.Fields.single_field(Y, prop_chain))
        @info "    Creating profile for `$var_name`"
        var = CC.Fields.single_field(Y, prop_chain)
        temporary = ClimaCore.column(var, 1, 1, 1)
        ϕ_col_ave = deepcopy(vec(temporary))
        ϕ_col_std = deepcopy(vec(temporary))
        ϕ_col_ave .= 0
        ϕ_col_std .= 0
        local_geom = CC.Fields.local_geometry_field(axes(var))
        z_f = ClimaCore.column(local_geom, 1, 1, 1)
        z_f = z_f.coordinates.z
        z = vec(z_f)
        for inds in TC.iterate_columns(Y)
            ϕ_col = ClimaCore.column(var, inds...)
            ϕ_col_ave .+= vec(ϕ_col) ./ n_columns
        end
        for inds in TC.iterate_columns(Y)
            ϕ_col = ClimaCore.column(var, inds...)
            ϕ_col_std .+= sqrt.((vec(ϕ_col) .- ϕ_col_ave) .^ 2 ./ n_columns)
        end

        Plots.plot(ϕ_col_ave, z ./ 1000; label = "Mean & Std", grid = false, xerror = ϕ_col_std, fillalpha = 0.5)
        Plots.plot!(; ylabel = "z [km]", xlabel = "$var_name", markershape = :circle)
        Plots.title!("$(space_string(var_space))")
        Plots.savefig(joinpath(output_dir, "$var_name.png"))
    end
end

function test_zero_horizontal_variance(Y; filter_prop_chain = pc -> true)
    # Column animations
    n_columns = TC.number_of_columns(Y)
    prop_chains = filter(filter_prop_chain, CC.Fields.property_chains(Y))
    std_per_var = map(CC.Fields.property_chains(Y)) do prop_chain
        var_name = processed_varname(prop_chain)
        var_space = axes(CC.Fields.single_field(Y, prop_chain))
        var = CC.Fields.single_field(Y, prop_chain)
        temporary = ClimaCore.column(var, 1, 1, 1)
        ϕ_col_ave = deepcopy(vec(temporary))
        ϕ_col_std = deepcopy(vec(temporary))
        ϕ_col_ave .= 0
        ϕ_col_std .= 0
        local_geom = CC.Fields.local_geometry_field(axes(var))
        z_f = ClimaCore.column(local_geom, 1, 1, 1)
        z_f = z_f.coordinates.z
        z = vec(z_f)
        for inds in TC.iterate_columns(Y)
            ϕ_col = ClimaCore.column(var, inds...)
            ϕ_col_ave .+= vec(ϕ_col) ./ n_columns
        end
        for inds in TC.iterate_columns(Y)
            ϕ_col = ClimaCore.column(var, inds...)
            ϕ_col_std .+= sqrt.((vec(ϕ_col) .- ϕ_col_ave) .^ 2 ./ n_columns)
        end
        (prop_chain, ϕ_col_std)
    end
    success_dict = Dict(
        (:cent, :ρ) => false,
        (:cent, :uₕ, :components, :data, 1) => false,
        (:cent, :uₕ, :components, :data, 2) => false,
        (:cent, :ρe_tot) => false,
        (:cent, :ρq_tot) => false,
        (:cent, :turbconv, :en, :ρatke) => false,
        (:cent, :turbconv, :up, 1, :ρarea) => false,
        (:cent, :turbconv, :up, 1, :ρaq_tot) => false,
        (:cent, :turbconv, :up, 1, :ρae_tot) => false,
        (:cent, :turbconv, :pr, :q_rai) => true,
        (:cent, :turbconv, :pr, :q_sno) => true,
        (:face, :w, :components, :data, 1) => true,
        (:face, :turbconv, :up, 1, :ρaw) => false,
    )
    @testset "Identical columns" begin
        for (prop_chain, ϕ_col_std) in std_per_var
            bool = all(ϕ_col_std .== 0)
            if success_dict[prop_chain]
                bool || @show prop_chain
                @test bool
            else
                bool && @show prop_chain
                @test_broken bool
            end
        end
    end

end
