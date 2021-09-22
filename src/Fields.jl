
function center_field(grid::Grid)
    return zeros(grid.nz)
end

function center_field(grid::Grid, nu::Int)
    return zeros(nu, grid.nz)
end

# TODO: change to use unequal length to center_field
function face_field(grid::Grid)
    return zeros(grid.nz + 1)
end

function face_field(grid::Grid, nu::Int)
    return zeros(nu, grid.nz + 1)
end

function field(grid::Grid, loc::String)
    @assert loc == "full" || loc == "half"
    if loc == "full"
        return face_field(grid)
    else
        return center_field(grid)
    end
end

function field(grid::Grid, loc::String, nu::Int)
    @assert loc == "full" || loc == "half"
    if loc == "full"
        return face_field(grid, nu)
    else
        return center_field(grid, nu)
    end
end
