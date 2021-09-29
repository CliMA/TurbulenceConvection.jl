
# ClimaCore hacks
Base.isless(zp1::CC.Geometry.ZPoint, zp2::CC.Geometry.ZPoint) = Base.isless(zp1.z, zp2.z)
Base.isless(zp::CC.Geometry.ZPoint, val::Number) = Base.isless(zp.z, val)
Base.isless(val::Number, zp::CC.Geometry.ZPoint) = Base.isless(val, zp.z)

Base.:+(zp::CC.Geometry.ZPoint, val) = Base.:+(zp.z, val)
Base.:+(val, zp::CC.Geometry.ZPoint) = Base.:+(val, zp.z)
Base.:/(val, zp::CC.Geometry.ZPoint) = Base.:/(val, zp.z)
Base.:+(zp1::CC.Geometry.ZPoint, zp2::CC.Geometry.ZPoint) = Base.:+(zp1.z, zp2.z)

Base.:*(zp1::CC.Geometry.ZPoint, zp2::CC.Geometry.ZPoint) = Base.:*(zp1.z, zp2.z)

Base.log(zp::CC.Geometry.ZPoint) = Base.log(zp.z)
Base.:-(zp::CC.Geometry.ZPoint) = Base.:-(zp.z)
Base.:-(zp::CC.Geometry.ZPoint, val) = Base.:-(zp.z, val)
Base.:-(val, zp::CC.Geometry.ZPoint) = Base.:-(val, zp.z)
Base.:-(zp1::CC.Geometry.ZPoint, zp2::CC.Geometry.ZPoint) = Base.:-(zp1.z, zp2.z)

Base.convert(::Type{Float64}, zp::CC.Geometry.ZPoint) = zp.z


struct Grid{FT, CS, FS, SC, SF}
    zmin::FT
    zmax::FT
    Δz::FT
    Δzi::FT
    nz::Int
    cs::CS
    fs::FS
    zc::SC
    zf::SF
    function Grid(namelist)

        #Get the grid spacing
        Δz = namelist["grid"]["dz"]
        nz = namelist["grid"]["nz"]
        FT = Float64
        z₀, z₁ = FT(0), FT(nz * Δz)

        domain = CC.Domains.IntervalDomain(
            CC.Geometry.ZPoint{FT}(z₀),
            CC.Geometry.ZPoint{FT}(z₁),
            boundary_tags = (:bottom, :top),
        )

        mesh = CC.Meshes.IntervalMesh(domain, nelems = nz)

        cs = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
        fs = CC.Spaces.FaceFiniteDifferenceSpace(cs)
        zc = CC.Fields.coordinate_field(cs)
        zf = CC.Fields.coordinate_field(fs)

        #Set the inverse grid spacing
        Δzi = 1.0 / Δz

        zmin = minimum(parent(zf))
        zmax = maximum(parent(zf))
        CS = typeof(cs)
        FS = typeof(fs)
        SC = typeof(zc)
        SF = typeof(zf)
        FT = typeof(zmax)
        return new{FT, CS, FS, SC, SF}(zmin, zmax, Δz, Δzi, nz, cs, fs, zc, zf)
    end
end

# Index of the first interior cell above the surface
kc_surface(grid::Grid) = 1
kf_surface(grid::Grid) = 1
kc_top_of_atmos(grid::Grid) = grid.nz
kf_top_of_atmos(grid::Grid) = grid.nz + 1

is_surface_center(grid::Grid, k::Int) = k == kc_surface(grid)
is_toa_center(grid::Grid, k::Int) = k == kc_top_of_atmos(grid)
is_surface_face(grid::Grid, k::Int) = k == kf_surface(grid)
is_toa_face(grid::Grid, k::Int) = k == kf_top_of_atmos(grid)

zc_surface(grid::Grid) = grid.zc[kc_surface(grid)]
zf_surface(grid::Grid) = grid.zf[kf_surface(grid)]
zc_toa(grid::Grid) = grid.zc[kc_top_of_atmos(grid)]
zf_toa(grid::Grid) = grid.zf[kf_top_of_atmos(grid)]

real_center_indices(grid::Grid) = kc_surface(grid):kc_top_of_atmos(grid)
real_face_indices(grid::Grid) = kf_surface(grid):kf_top_of_atmos(grid)

face_space(grid::Grid) = grid.fs
center_space(grid::Grid) = grid.cs


Base.eltype(::Grid{FT}) where {FT} = FT
