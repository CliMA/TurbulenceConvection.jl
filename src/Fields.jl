
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

surface_value(f::AbstractVector, grid::Grid) = f[kf_surface(grid)]

Base.length(space::CC.Spaces.FiniteDifferenceSpace) = length(CC.Spaces.coordinates_data(space))

Base.@propagate_inbounds function Base.getindex(field::CC.Fields.FiniteDifferenceField, i::Integer)
    Base.getindex(CC.Fields.field_values(field), i)
end

Base.@propagate_inbounds function Base.setindex!(field::CC.Fields.FiniteDifferenceField, v, i::Integer)
    Base.setindex!(CC.Fields.field_values(field), v, i)
end

function FieldFromNamedTuple(space, nt::NamedTuple)
    cmv(z) = nt
    return cmv.(CC.Fields.coordinate_field(space))
end
