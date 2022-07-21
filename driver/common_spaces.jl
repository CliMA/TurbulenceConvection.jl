import ClimaCore
const CC = ClimaCore

function periodic_line_mesh(; x_max, x_elem)
    domain = CC.Domains.IntervalDomain(CCG.XPoint(zero(x_max)), CCG.XPoint(x_max); periodic = true)
    return CC.Meshes.IntervalMesh(domain; nelems = x_elem)
end

function periodic_rectangle_mesh(; x_max, y_max, x_elem, y_elem)
    x_domain = CC.Domains.IntervalDomain(CCG.XPoint(zero(x_max)), CCG.XPoint(x_max); periodic = true)
    y_domain = CC.Domains.IntervalDomain(CCG.YPoint(zero(y_max)), CCG.YPoint(y_max); periodic = true)
    domain = CC.Domains.RectangleDomain(x_domain, y_domain)
    return CC.Meshes.RectilinearMesh(domain, x_elem, y_elem)
end

# h_elem is the number of elements per side of every panel (6 panels in total)
function cubed_sphere_mesh(; radius, h_elem)
    domain = CC.Domains.SphereDomain(radius)
    return CC.Meshes.EquiangularCubedSphere(domain, h_elem)
end

function make_horizontal_space(mesh, quad)
    if mesh isa CC.Meshes.AbstractMesh1D
        topology = CC.Topologies.IntervalTopology(mesh)
        space = CC.Spaces.SpectralElementSpace1D(topology, quad)
    elseif mesh isa CC.Meshes.AbstractMesh2D
        topology = CC.Topologies.Topology2D(mesh)
        space = CC.Spaces.SpectralElementSpace2D(topology, quad)
    end
    return space
end

function make_hybrid_spaces(h_space, z_max, z_mesh)
    FT = eltype(z_max)

    @info "z heights" z_mesh.faces
    z_topology = CC.Topologies.IntervalTopology(z_mesh)
    z_space = CC.Spaces.CenterFiniteDifferenceSpace(z_topology)
    center_space = CC.Spaces.ExtrudedFiniteDifferenceSpace(h_space, z_space)
    face_space = CC.Spaces.FaceExtrudedFiniteDifferenceSpace(center_space)

    svpc_domain =
        CC.Domains.IntervalDomain(CC.Geometry.ZPoint{FT}(0), CC.Geometry.ZPoint{FT}(1), boundary_tags = (:bottom, :top))
    svpc_mesh = CC.Meshes.IntervalMesh(svpc_domain, nelems = 1)
    svpc_space = CC.Spaces.CenterFiniteDifferenceSpace(svpc_mesh)

    return center_space, face_space, svpc_space
end

function construct_mesh(namelist; FT = Float64)

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

    z_max = FT(nz * Δz)

    if truncated_gcm_mesh
        nz_s= namelist["grid"]["stretch"]["nz"]
        Δz_s_surf = FT(namelist["grid"]["stretch"]["dz_surf"])
        Δz_s_top = FT(namelist["grid"]["stretch"]["dz_toa"])
        z_s_toa = FT(namelist["grid"]["stretch"]["z_toa"])

        z_stretch = CC.Meshes.GeneralizedExponentialStretching(Δz_s_surf, Δz_s_top)
        z_domain = CC.Domains.IntervalDomain(
            CCG.ZPoint(zero(z_s_toa)),
            CCG.ZPoint(z_s_toa);
            boundary_tags = (:bottom, :top),
        )
        gcm_mesh = CC.Meshes.IntervalMesh(z_domain, z_stretch; nelems = nz_s)
        z_mesh = TC.TCMeshFromGCMMesh(gcm_mesh; z_max = z_max)
    else
        z_stretch = CC.Meshes.Uniform()
        z_domain = CC.Domains.IntervalDomain(
            CCG.ZPoint(zero(z_max)),
            CCG.ZPoint(z_max);
            boundary_tags = (:bottom, :top)
        )
        z_mesh = CC.Meshes.IntervalMesh(z_domain, z_stretch; nelems = nz)
    end

    return (; z_mesh = z_mesh, z_max = z_max)
end

function get_spaces(namelist, param_set, FT)
    center_space, face_space, svpc_space = if namelist["config"] == "sphere"

        h_elem = 1
        quad = CC.Spaces.Quadratures.GLL{2}()
        horizontal_mesh = cubed_sphere_mesh(; radius = FT(TCP.planet_radius(param_set)), h_elem)
        h_space = make_horizontal_space(horizontal_mesh, quad)

        (; z_mesh, z_max) = construct_mesh(namelist; FT = FT)

        center_space, face_space, svpc_space = make_hybrid_spaces(h_space, z_max, z_mesh)
        return center_space, face_space, svpc_space

    elseif namelist["config"] == "column" # single column (default)

        Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
        quad = CC.Spaces.Quadratures.GL{1}()
        horizontal_mesh = periodic_rectangle_mesh(; x_max = Δx, y_max = Δx, x_elem = 1, y_elem = 1)
        h_space = make_horizontal_space(horizontal_mesh, quad)

        (; z_mesh, z_max) = construct_mesh(namelist; FT = FT)

        center_space, face_space, svpc_space = make_hybrid_spaces(h_space, z_max, z_mesh)
        return center_space, face_space, svpc_space
    end
end
