struct Grid{A1, FT}
    zmin::FT
    zmax::FT
    cinterior::UnitRange
    finterior::UnitRange
    dz::FT
    dzi::FT
    gw::Int
    nz::Int
    nzg::Int
    z_half::A1
    z::A1
    function Grid(namelist)

        #Get the grid spacing
        dz = namelist["grid"]["dz"]

        #Set the inverse grid spacing

        dzi = 1.0 / dz

        #Get the grid dimensions and ghost points
        gw = namelist["grid"]["gw"]
        nz = namelist["grid"]["nz"]
        nzg = nz + 2 * gw

        cinterior = (gw + 1):(nzg - gw)
        finterior = (gw + 1):(nzg - gw)

        # TODO: make cell centers and cell faces different sizes
        z_half = zeros(nz + 2 * gw)
        z = zeros(nz + 2 * gw)
        count = 1
        @inbounds for i in range(1 - gw, stop = nz + gw)
            z[count] = i * dz
            z_half[count] = (i - 0.5) * dz
            count += 1
        end
        zmin = z[gw]
        zmax = z[nzg - gw]
        A1 = typeof(z)
        FT = typeof(zmax)
        return new{A1, FT}(zmin, zmax, cinterior, finterior, dz, dzi, gw, nz, nzg, z_half, z)
    end
end

# Index of the first interior cell above the surface
kc_surface(grid::Grid) = grid.gw + 1
kf_surface(grid::Grid) = grid.gw
kc_top_of_atmos(grid::Grid) = grid.nzg - grid.gw
kf_top_of_atmos(grid::Grid) = grid.nzg - grid.gw

is_surface_center(grid::Grid, k::Int) = k == kc_surface(grid)
is_toa_center(grid::Grid, k::Int) = k == kc_top_of_atmos(grid)
is_surface_face(grid::Grid, k::Int) = k == kf_surface(grid)
is_toa_face(grid::Grid, k::Int) = k == kf_top_of_atmos(grid)

zc_surface(grid::Grid) = grid.z_half[kc_surface(grid)]
zf_surface(grid::Grid) = grid.z[kf_surface(grid)]
zc_toa(grid::Grid) = grid.z_half[kc_top_of_atmos(grid)]
zf_toa(grid::Grid) = grid.z[kf_top_of_atmos(grid)]

real_center_indices(grid::Grid) = kc_surface(grid):kc_top_of_atmos(grid)
real_face_indices(grid::Grid) = kf_surface(grid):kf_top_of_atmos(grid)


Base.eltype(::Grid{A1, FT}) where {A1, FT} = FT
