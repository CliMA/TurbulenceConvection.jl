
function center_field(grid::Grid)
    return zeros(grid.nz)
end

function center_field(grid::Grid, nu::Int)
    return zeros(nu, grid.nz)
end

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

# A complementary struct to ClimaCore's `PlusHalf` type.
struct Cent{I <: Integer}
    i::I
end

Base.:+(h::Cent) = h
Base.:-(h::Cent) = Cent(-h.i - one(h.i))
Base.:+(i::Integer, h::Cent) = Cent(i + h.i)
Base.:+(h::Cent, i::Integer) = Cent(h.i + i)
Base.:+(h1::Cent, h2::Cent) = h1.i + h2.i + one(h1.i)
Base.:-(i::Integer, h::Cent) = Cent(i - h.i - one(h.i))
Base.:-(h::Cent, i::Integer) = Cent(h.i - i)
Base.:-(h1::Cent, h2::Cent) = h1.i - h2.i

Base.:<=(h1::Cent, h2::Cent) = h1.i <= h2.i
Base.:<(h1::Cent, h2::Cent) = h1.i < h2.i
Base.max(h1::Cent, h2::Cent) = Cent(max(h1.i, h2.i))
Base.min(h1::Cent, h2::Cent) = Cent(min(h1.i, h2.i))


Base.@propagate_inbounds Base.getindex(field::CC.Fields.FiniteDifferenceField, i::Integer) = Base.getproperty(field, i)

Base.@propagate_inbounds Base.getindex(field::CC.Fields.CenterFiniteDifferenceField, i::Cent) =
    Base.getindex(CC.Fields.field_values(field), i.i)
Base.@propagate_inbounds Base.setindex!(field::CC.Fields.CenterFiniteDifferenceField, v, i::Cent) =
    Base.setindex!(CC.Fields.field_values(field), v, i.i)

Base.@propagate_inbounds Base.getindex(field::CC.Fields.FaceFiniteDifferenceField, i::CCO.PlusHalf) =
    Base.getindex(CC.Fields.field_values(field), i.i)
Base.@propagate_inbounds Base.setindex!(field::CC.Fields.FaceFiniteDifferenceField, v, i::CCO.PlusHalf) =
    Base.setindex!(CC.Fields.field_values(field), v, i.i)

Base.@propagate_inbounds Base.getindex(field::CC.Fields.FaceFiniteDifferenceField, ::Cent) =
    error("Attempting to getindex with a center index (Cent) into a Face field")
Base.@propagate_inbounds Base.getindex(field::CC.Fields.CenterFiniteDifferenceField, ::CCO.PlusHalf) =
    error("Attempting to getindex with a face index (PlusHalf) into a Center field")

Base.@propagate_inbounds Base.setindex!(field::CC.Fields.FaceFiniteDifferenceField, v, ::Cent) =
    error("Attempting to setindex with a center index (Cent) into a Face field")
Base.@propagate_inbounds Base.setindex!(field::CC.Fields.CenterFiniteDifferenceField, v, ::CCO.PlusHalf) =
    error("Attempting to setindex with a face index (PlusHalf) into a Center field")

# TODO: deprecate, we should not overload getindex/setindex for ordinary arrays.
Base.@propagate_inbounds Base.getindex(arr::AbstractArray, i::Cent) = Base.getindex(arr, i.i)
Base.@propagate_inbounds Base.setindex!(arr::AbstractArray, v, i::Cent) = Base.setindex!(arr, v, i.i)
Base.@propagate_inbounds Base.getindex(arr::AbstractArray, i::CCO.PlusHalf) = Base.getindex(arr, i.i)
Base.@propagate_inbounds Base.setindex!(arr::AbstractArray, v, i::CCO.PlusHalf) = Base.setindex!(arr, v, i.i)
Base.@propagate_inbounds Base.getindex(arr::AbstractArray, i::Int, j::Cent) = Base.getindex(arr, i, j.i)
Base.@propagate_inbounds Base.setindex!(arr::AbstractArray, v, i::Int, j::Cent) = Base.setindex!(arr, v, i, j.i)
Base.@propagate_inbounds Base.getindex(arr::AbstractArray, i::Int, j::CCO.PlusHalf) = Base.getindex(arr, i, j.i)
Base.@propagate_inbounds Base.setindex!(arr::AbstractArray, v, i::Int, j::CCO.PlusHalf) = Base.setindex!(arr, v, i, j.i)

function FieldFromNamedTuple(space, nt::NamedTuple)
    cmv(z) = nt
    return cmv.(CC.Fields.coordinate_field(space))
end
