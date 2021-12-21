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

# https://github.com/CliMA/ClimaCore.jl/issues/275
transform_broadcasted(bc::Base.Broadcast.Broadcasted{CC.Fields.FieldVectorStyle}, symb, axes) =
    Base.Broadcast.Broadcasted(bc.f, map(arg -> transform_broadcasted(arg, symb, axes), bc.args), axes)
transform_broadcasted(fv::CC.Fields.FieldVector, symb, axes) = parent(getproperty(fv, symb))
transform_broadcasted(x, symb, axes) = x
@inline function Base.copyto!(dest::CC.Fields.FieldVector, bc::Base.Broadcast.Broadcasted{CC.Fields.FieldVectorStyle})
    for symb in propertynames(dest)
        p = parent(getproperty(dest, symb))
        Base.copyto!(p, transform_broadcasted(bc, symb, axes(p)))
    end
    return dest
end

function Base.cumsum(field::CC.Fields.FiniteDifferenceField)
    Base.cumsum(parent(CC.Fields.weighted_jacobian(field)) .* vec(field), dims = 1)
end
function Base.cumsum!(fieldout::CC.Fields.FiniteDifferenceField, fieldin::CC.Fields.FiniteDifferenceField)
    Base.cumsum!(fieldout, parent(CC.Fields.weighted_jacobian(fieldin)) .* vec(fieldin), dims = 1)
end
