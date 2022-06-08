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

toscalar(x::CCG.Covariant3Vector) = x.u₃

const FDFields = Union{CC.Fields.ExtrudedFiniteDifferenceField, CC.Fields.FiniteDifferenceField}

const FaceFields = Union{CC.Fields.FaceExtrudedFiniteDifferenceField, CC.Fields.FaceFiniteDifferenceField}

const CenterFields = Union{CC.Fields.CenterExtrudedFiniteDifferenceField, CC.Fields.CenterFiniteDifferenceField}

Base.@propagate_inbounds Base.getindex(field::FDFields, i::Integer) = Base.getproperty(field, i)

Base.@propagate_inbounds Base.getindex(field::CenterFields, i::Cent) = Base.getindex(CC.Fields.field_values(field), i.i)
Base.@propagate_inbounds Base.setindex!(field::CenterFields, v, i::Cent) =
    Base.setindex!(CC.Fields.field_values(field), v, i.i)

Base.@propagate_inbounds Base.getindex(field::FaceFields, i::CCO.PlusHalf) =
    Base.getindex(CC.Fields.field_values(field), i.i)
Base.@propagate_inbounds Base.setindex!(field::FaceFields, v, i::CCO.PlusHalf) =
    Base.setindex!(CC.Fields.field_values(field), v, i.i)

Base.@propagate_inbounds Base.getindex(field::FaceFields, ::Cent) =
    error("Attempting to getindex with a center index (Cent) into a Face field")
Base.@propagate_inbounds Base.getindex(field::CenterFields, ::CCO.PlusHalf) =
    error("Attempting to getindex with a face index (PlusHalf) into a Center field")

Base.@propagate_inbounds Base.setindex!(field::FaceFields, v, ::Cent) =
    error("Attempting to setindex with a center index (Cent) into a Face field")
Base.@propagate_inbounds Base.setindex!(field::CenterFields, v, ::CCO.PlusHalf) =
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

# Constant field
function FieldFromNamedTuple(space, nt::NamedTuple)
    cmv(z) = nt
    return cmv.(CC.Fields.coordinate_field(space))
end
# Non-constant field
function FieldFromNamedTuple(space, initial_conditions::Function, ::Type{FT}, params...) where {FT}
    local_geometry = CC.Fields.local_geometry_field(space)
    return initial_conditions.(FT, local_geometry, params...)
end

function Base.cumsum(field::CC.Fields.FiniteDifferenceField)
    Base.cumsum(parent(CC.Fields.weighted_jacobian(field)) .* vec(field), dims = 1)
end
function Base.cumsum!(fieldout::CC.Fields.FiniteDifferenceField, fieldin::CC.Fields.FiniteDifferenceField)
    Base.cumsum!(fieldout, parent(CC.Fields.weighted_jacobian(fieldin)) .* vec(fieldin), dims = 1)
end

get_Δz(field::CC.Fields.FiniteDifferenceField) = parent(CC.Fields.weighted_jacobian(field))

# TODO: move these things into ClimaCore

isa_center_space(space) = false
isa_center_space(::CC.Spaces.CenterFiniteDifferenceSpace) = true
isa_center_space(::CC.Spaces.CenterExtrudedFiniteDifferenceSpace) = true

isa_face_space(space) = false
isa_face_space(::CC.Spaces.FaceFiniteDifferenceSpace) = true
isa_face_space(::CC.Spaces.FaceExtrudedFiniteDifferenceSpace) = true

function first_center_space(fv::CC.Fields.FieldVector)
    for prop_chain in CC.Fields.property_chains(fv)
        f = CC.Fields.single_field(fv, prop_chain)
        space = axes(f)
        isa_center_space(space) && return space
    end
    error("Unfound space")
end

function first_face_space(fv::CC.Fields.FieldVector)
    for prop_chain in CC.Fields.property_chains(fv)
        f = CC.Fields.single_field(fv, prop_chain)
        space = axes(f)
        isa_face_space(space) && return space
    end
    error("Unfound space")
end

"""
    set_z!(field, f::Function)

Populate field `field` with values computed
from function `f(z)`
"""
function set_z!(field::CC.Fields.Field, f::F) where {F::Function}
    space = axes(field)
    z = CC.Fields.coordinate_field(space).z
    lg = CC.Fields.local_geometry_field(space)
    @. field = f(z)
    return nothing
end

"""
    set_z!(field, u::Function, v::Function)

Populate field `field` with values computed
from function `Covariant12Vector(UVVector(u(z), v(z)))`
"""
function set_z!(field::CC.Fields.Field, u::U, v::V) where {U::Function, V::Function}
    space = axes(field)
    z = CC.Fields.coordinate_field(space).z
    lg = CC.Fields.local_geometry_field(space)
    @. field = CCO.Covariant12Vector(CCG.UVVector(u(z), v(z)), lg)
    return nothing
end

"""
    z_function(field)

Create a callable function of `z` from field `field`
"""
function z_function(field::CC.Fields.Field)
    space = axes(field)
    z = CC.Fields.coordinate_field(space).z
    return Dierckx.Spline1D(z, vec(field); k = 1)
end
