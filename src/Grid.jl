struct Grid{A1, FT}
    zmin::FT
    zmax::FT
    dz::FT
    dzi::FT
    nz::Int
    nzg::Int
    z_half::A1
    z::A1
    function Grid(namelist)

        #Get the grid spacing
        dz = namelist["grid"]["dz"]

        #Set the inverse grid spacing

        dzi = 1.0 / dz

        # Get the grid dimensions
        nz = namelist["grid"]["nz"]
        nzg = nz

        # TODO: make cell centers and cell faces different sizes
        z_half = zeros(nz)
        z = zeros(nz + 1)

        z .= map(1:(nz + 1)) do k
            (k - 1.0) * dz
        end
        z_half .= map(1:nz) do k
            (z[k] + z[k + 1]) / 2.0
        end

        zmin = z[1]
        zmax = z[nzg + 1]
        A1 = typeof(z)
        FT = typeof(zmax)
        return new{A1, FT}(zmin, zmax, dz, dzi, nz, nzg, z_half, z)
    end
end

# Index of the first interior cell above the surface
kc_surface(grid::Grid) = 1
kf_surface(grid::Grid) = 1
kc_top_of_atmos(grid::Grid) = grid.nzg
kf_top_of_atmos(grid::Grid) = grid.nzg + 1

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
