import ClimaCore
const CC = ClimaCore

function construct_grid(namelist; FT = Float64)

    truncated_gcm_mesh = TC.parse_namelist(namelist, "grid", "stretch", "flag"; default = false)

    if Cases.get_case(namelist) == Cases.LES_driven_SCM()
        Δz = get(namelist["grid"], "dz", nothing)
        nz = get(namelist["grid"], "nz", nothing)
        @assert isnothing(Δz) ⊻ isnothing(nz) string(
            "LES_driven_SCM supports nz or Δz, not both.",
            "The domain height is enforced to be the same as in LES.",
        )

        les_filename = namelist["meta"]["lesfile"]
        TC.valid_lespath(les_filename)
        zmax = NC.Dataset(les_filename, "r") do data
            Array(TC.get_nc_data(data, "zf"))[end]
        end
        nz = isnothing(nz) ? Int(zmax ÷ Δz) : Int(nz)
        Δz = isnothing(Δz) ? FT(zmax ÷ nz) : FT(Δz)
    else
        Δz = FT(namelist["grid"]["dz"])
        nz = namelist["grid"]["nz"]
    end

    z₀, z₁ = FT(0), FT(nz * Δz)
    if truncated_gcm_mesh
        nzₛ = namelist["grid"]["stretch"]["nz"]
        Δzₛ_surf = FT(namelist["grid"]["stretch"]["dz_surf"])
        Δzₛ_top = FT(namelist["grid"]["stretch"]["dz_toa"])
        zₛ_toa = FT(namelist["grid"]["stretch"]["z_toa"])
        stretch = CC.Meshes.GeneralizedExponentialStretching(Δzₛ_surf, Δzₛ_top)
        domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{FT}(z₀),
            CC.Geometry.ZPoint{FT}(zₛ_toa),
            boundary_tags = (:bottom, :top),
        )
        gcm_mesh = CC.Meshes.IntervalMesh(domain, stretch; nelems = nzₛ)
        mesh = TC.TCMeshFromGCMMesh(gcm_mesh; z_max = z₁)
    else
        CC.Meshes.Uniform()
        domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{FT}(z₀),
            CC.Geometry.ZPoint{FT}(z₁),
            boundary_tags = (:bottom, :top),
        )
        mesh = CC.Meshes.IntervalMesh(domain, nelems = nz)
    end
    @info "z heights" mesh.faces
    return TC.Grid(mesh)
end
