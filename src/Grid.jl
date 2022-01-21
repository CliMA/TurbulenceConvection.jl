
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


struct Grid{FT, NZ, CS, FS, SC, SF}
    zmin::FT
    zmax::FT
    Δz::FT
    cs::CS
    fs::FS
    zc::SC
    zf::SF
    function Grid(mesh)

        nz = length(mesh.faces) - 1
        cs = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
        fs = CC.Spaces.FaceFiniteDifferenceSpace(cs)
        zc = CC.Fields.coordinate_field(cs)
        zf = CC.Fields.coordinate_field(fs)
        Δz = zf[CCO.PlusHalf(2)] - zf[CCO.PlusHalf(1)]

        FT = eltype(parent(zf))
        zmin = minimum(parent(zf))
        zmax = maximum(parent(zf))
        CS = typeof(cs)
        FS = typeof(fs)
        SC = typeof(zc)
        SF = typeof(zf)
        return new{FT, nz, CS, FS, SC, SF}(zmin, zmax, Δz, cs, fs, zc, zf)
    end
end

function Grid(Δz::FT, nz::Int) where {FT <: AbstractFloat}
    z₀, z₁ = FT(0), FT(nz * Δz)

    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(z₀),
        CC.Geometry.ZPoint{FT}(z₁),
        boundary_tags = (:bottom, :top),
    )

    mesh = CC.Meshes.IntervalMesh(domain, nelems = nz)
    return Grid(mesh)
end

n_cells(::Grid{FT, NZ}) where {FT, NZ} = NZ

# Index of the first interior cell above the surface
kc_surface(grid::Grid) = Cent(1)
kf_surface(grid::Grid) = CCO.PlusHalf(1)
kc_top_of_atmos(grid::Grid) = Cent(n_cells(grid))
kf_top_of_atmos(grid::Grid) = CCO.PlusHalf(n_cells(grid) + 1)

is_surface_center(grid::Grid, k) = k == kc_surface(grid)
is_toa_center(grid::Grid, k) = k == kc_top_of_atmos(grid)
is_surface_face(grid::Grid, k) = k == kf_surface(grid)
is_toa_face(grid::Grid, k) = k == kf_top_of_atmos(grid)

zc_surface(grid::Grid) = grid.zc[kc_surface(grid)]
zf_surface(grid::Grid) = grid.zf[kf_surface(grid)]
zc_toa(grid::Grid) = grid.zc[kc_top_of_atmos(grid)]
zf_toa(grid::Grid) = grid.zf[kf_top_of_atmos(grid)]

real_center_indices(grid::Grid) = CenterIndices(grid)
real_face_indices(grid::Grid) = FaceIndices(grid)

struct FaceIndices{Nstart, Nstop, G}
    grid::G
    function FaceIndices(grid::G) where {G <: Grid}
        Nstart, Nstop = kf_surface(grid).i, kf_top_of_atmos(grid).i
        new{Nstart, Nstop, G}(grid)
    end
end

struct CenterIndices{Nstart, Nstop, G}
    grid::G
    function CenterIndices(grid::G) where {G <: Grid}
        Nstart, Nstop = kc_surface(grid).i, kc_top_of_atmos(grid).i
        new{Nstart, Nstop, G}(grid)
    end
end

Base.keys(ci::CenterIndices) = 1:length(ci)
Base.keys(fi::FaceIndices) = 1:length(fi)

n_start(::CenterIndices{Nstart}) where {Nstart} = Nstart
n_start(::FaceIndices{Nstart}) where {Nstart} = Nstart
n_stop(::CenterIndices{Nstart, Nstop}) where {Nstart, Nstop} = Nstop
n_stop(::FaceIndices{Nstart, Nstop}) where {Nstart, Nstop} = Nstop

Base.getindex(ci::CenterIndices, i::Int) = Cent(Base.getindex(n_start(ci):n_stop(ci), i))
Base.getindex(fi::FaceIndices, i::Int) = CCO.PlusHalf(Base.getindex(n_start(fi):n_stop(fi), i))

Base.length(::FaceIndices{Nstart, Nstop}) where {Nstart, Nstop} = Nstop - Nstart + 1
Base.length(::CenterIndices{Nstart, Nstop}) where {Nstart, Nstop} = Nstop - Nstart + 1

Base.iterate(fi::CenterIndices{Nstart, Nstop}, state = Nstart) where {Nstart, Nstop} =
    state > Nstop ? nothing : (Cent(state), state + 1)

Base.iterate(fi::FaceIndices{Nstart, Nstop}, state = Nstart) where {Nstart, Nstop} =
    state > Nstop ? nothing : (CCO.PlusHalf(state), state + 1)

Base.iterate(fi::Base.Iterators.Reverse{T}, state = Nstop) where {Nstart, Nstop, T <: CenterIndices{Nstart, Nstop}} =
    state < Nstart ? nothing : (Cent(state), state - 1)

Base.iterate(fi::Base.Iterators.Reverse{T}, state = Nstop) where {Nstart, Nstop, T <: FaceIndices{Nstart, Nstop}} =
    state < Nstart ? nothing : (CCO.PlusHalf(state), state - 1)

face_space(grid::Grid) = grid.fs
center_space(grid::Grid) = grid.cs

#=
    findfirst_center
    findlast_center
    findfirst_face
    findlast_face

Grid-aware find-first / find-last indices with
surface/toa (respectively) as the default index
=#

function findfirst_center(f::Function, grid::Grid)
    RI = real_center_indices(grid)
    k = findfirst(f, RI)
    return RI[isnothing(k) ? kc_surface(grid).i : k]
end
function findlast_center(f::Function, grid::Grid)
    RI = real_center_indices(grid)
    k = findlast(f, RI)
    return RI[isnothing(k) ? kc_top_of_atmos(grid).i : k]
end
z_findfirst_center(f::F, grid::Grid) where {F} = grid.zc[findfirst_center(f, grid)].z
z_findlast_center(f::F, grid::Grid) where {F} = grid.zc[findlast_center(f, grid)].z

function findfirst_face(f::F, grid::Grid) where {F}
    RI = real_face_indices(grid)
    k = findfirst(f, RI)
    return RI[isnothing(k) ? kf_surface(grid).i : k]
end
function findlast_face(f::F, grid::Grid) where {F}
    RI = real_face_indices(grid)
    k = findlast(f, RI)
    return RI[isnothing(k) ? kf_top_of_atmos(grid).i : k]
end
z_findfirst_face(f::F, grid::Grid) where {F} = grid.zf[findfirst_face(f, grid)].z
z_findlast_face(f::F, grid::Grid) where {F} = grid.zf[findlast_face(f, grid)].z


Base.eltype(::Grid{FT}) where {FT} = FT
