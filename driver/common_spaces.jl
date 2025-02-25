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

function make_hybrid_spaces(h_space, z_mesh)
    FT = CC.Geometry.float_type(z_mesh.domain)
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

    elseif typeof(Cases.get_case(namelist)) <: Cases.SOCRATES # Maybe we dont wanna keep doing this so move the data part to SOCRATES data and just pass in the arg to all the surface_ref_state construcors instead...
        flight_number = namelist["meta"]["flight_number"]
        forcing_type = namelist["meta"]["forcing_type"]
        new_zc = vec(Cases.SSCF.get_default_new_z(flight_number))[:] # redundant constructor here, can we evade this somehow somewhere by passing in the object? this gets called before Cases stuff tho in main.

        # convert zc to zf (this is not guaranteed to always work but it does for the grids they gave us)
        new_z = FT[0] # zf
        for zc in new_zc
            append!(new_z, 2 * zc - new_z[end]) # new_z[end] + 2*(zc - zf_data[end]) = 2*zc - new_z[end] # This happens to work for the grids they gave us and not yield any negative numbers...
        end
        @assert all(diff(new_z) .> 0) "calulated face points zf from new_z is not monotonically increasing, try passing in a new_z that is a a valid center zc for some zf grid or adding functionality to interpret input grid as zf instead of as zc"

        # # probably better to change this to have some min dz since this reduction could be really overkill on a stretched/squeezed grid
        # new_z = new_z[begin:namelist["grid"]["z_reduction_factor"]:end] #  a way to reduce the resolution of the LES z so we can for example avoid CFL errors when running with a larger timestep

        old_z = new_z # zf from old zc
        # new_z = FT[]
        new_z = FT[0] # we can have a face at 0 right? then zc will be too close to ground but i think that's fine? for CFL it's on faces... but maybe it messes w/ BCs? But the meshes for non-Socrates use z₀ = 0 so it should be fine... I think the problem was CFL on ice sed since that's not on faces... but that's just going into the ground? i dont see/recall why it would be a problem... (actually it shouldn't be because we go C2F w/ wsed anyway...)
        current_z = FT(0) # first level can't be too close to sfc or else CFL still fails for faces
        dz_min = namelist["grid"]["dz_min"]


        dz_min_old = minimum(diff(old_z))
        if dz_min_old < dz_min
            @info(
                "minimum Δz $dz_min_old is smaller than the minimum allowed Δz of  $dz_min between z levels in the LES file, reducing..."
            )
            @info("old_z for $flight_number, $forcing_type: ", old_z)
            for (i, val) in enumerate(old_z)
                if (val - current_z) < dz_min # hopefully this is goin in the right direction...
                    continue
                else
                    append!(new_z, val)
                    current_z = val
                end
            end
            @info("new_z for $flight_number, $forcing_type: ", new_z)
        else
            new_z = old_z # we could interpolate to higher res but we won't for now...
        end

        z_mesh = CC.Geometry.ZPoint{FT}.(new_z) # added 0 to beginning? copy from the file #Array(TC.get_nc_data(data, "zc")) also idk what to do about paths like this
        # z_mesh = CC.Geometry.ZPoint{FT}.([FT(0), new_z...]) # added 0 to beginning? copy from the file #Array(TC.get_nc_data(data, "zc")) also idk what to do about paths like this (not needed if we do zc to zf conversion)
        nz = length(z_mesh)
        z₀, z₁ = z_mesh[1], z_mesh[end]
        zmax = z_mesh[end]
        domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{FT}(z₀),
            CC.Geometry.ZPoint{FT}(z₁),
            boundary_tags = (:bottom, :top),
        )
        z_mesh = CC.Meshes.IntervalMesh(domain, z_mesh)
        return (; z_mesh)
    else
        Δz = FT(namelist["grid"]["dz"])
        nz = namelist["grid"]["nz"]
    end

    z₀, z₁ = FT(0), FT(nz * Δz)

    z_mesh = if truncated_gcm_mesh
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
        TC.TCMeshFromGCMMesh(gcm_mesh; z_max = z₁)
    else
        stretch = CC.Meshes.Uniform()
        domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{FT}(z₀),
            CC.Geometry.ZPoint{FT}(z₁),
            boundary_tags = (:bottom, :top),
        )
        CC.Meshes.IntervalMesh(domain, stretch; nelems = nz)
    end
    return (; z_mesh)
end


function get_spaces(namelist, param_set, FT)
    center_space, face_space, svpc_space = if namelist["config"] == "sphere"
        h_elem = 1
        quad = CC.Spaces.Quadratures.GLL{2}()
        horizontal_mesh = cubed_sphere_mesh(; radius = FT(TCP.planet_radius(param_set)), h_elem)
        h_space = make_horizontal_space(horizontal_mesh, quad)
        (; z_mesh) = construct_mesh(namelist; FT = FT)
        center_space, face_space, svpc_space = make_hybrid_spaces(h_space, z_mesh)
        center_space, face_space, svpc_space
    elseif namelist["config"] == "column" # single column (default)
        Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
        quad = CC.Spaces.Quadratures.GL{1}()
        horizontal_mesh = periodic_rectangle_mesh(; x_max = Δx, y_max = Δx, x_elem = 1, y_elem = 1)
        h_space = make_horizontal_space(horizontal_mesh, quad)
        (; z_mesh) = construct_mesh(namelist; FT = FT)
        center_space, face_space, svpc_space = make_hybrid_spaces(h_space, z_mesh)
        center_space, face_space, svpc_space
    end
end
