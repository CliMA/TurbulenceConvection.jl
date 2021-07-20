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


        cinterior = gw:((nzg - gw) - 1)
        finterior = gw:((nzg - gw) - 1)

        # TODO: make cell centers and cell faces different sizes
        z_half = pyzeros(nz + 2 * gw)
        z = pyzeros(nz + 2 * gw)
        count = 0
        @inbounds for i in xrange(-gw, nz + gw)
            z[count] = (i + 1) * dz
            z_half[count] = (i + 0.5) * dz
            count += 1
        end
        zmin = z[gw - 1]
        zmax = z[nzg - gw - 1]
        A1 = typeof(z)
        FT = typeof(zmax)
        return new{A1, FT}(zmin, zmax, cinterior, finterior, dz, dzi, gw, nz, nzg, z_half, z)
    end
end

center_indicies(grid::Grid) = xrange(grid.nzg)
face_indicies(grid::Grid) = xrange(grid.nzg)
real_center_indicies(grid::Grid) = xrange(grid.gw, grid.nzg - grid.gw)
real_face_indicies(grid::Grid) = xrange(grid.gw, grid.nzg - grid.gw)

Base.eltype(::Grid{A1, FT}) where {A1, FT} = FT
