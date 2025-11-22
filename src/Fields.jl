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
toscalar(x::CCG.Contravariant3Vector) = CC.Geometry.WVector(x).components.data.:1 # from dennis., not sure precisely why it works... or why it needs to wvector... [ think method might not exist properly like this]

@inline function zero_field!(field::CC.Fields.FiniteDifferenceField)
    fill!(parent(field), zero(eltype(parent(field))))
end

const FDFields = Union{CC.Fields.ExtrudedFiniteDifferenceField, CC.Fields.FiniteDifferenceField}

const FaceFields = Union{CC.Fields.FaceExtrudedFiniteDifferenceField, CC.Fields.FaceFiniteDifferenceField}

const CenterFields = Union{CC.Fields.CenterExtrudedFiniteDifferenceField, CC.Fields.CenterFiniteDifferenceField}

Base.@propagate_inbounds Base.getindex(field::FDFields, i::Integer) = Base.getproperty(field, i)

# Not sure if this is a ClimaCore or a julia 1.11 problem, seemed 1.10.8 with CC 0.14.11 still didn't work tho
if VERSION ≥ v"1.11.0-beta"
# if VERSION ≥ v"1.10.8"

    """
    My hotfix attempt to define this method for julia 1.11
        Leads to segfault though lol
    """
    function Base.getindex(field_values::CC.DataLayouts.VIJFH{S}, v::Integer) where {S}
        i,j,_,_,h = size(field_values) # f will be 1, and we replace i... (our) (actually we wanna replace v)
        # i, j, _, v, h = I.I
        @inbounds CC.DataLayouts.get_struct(parent(field_values), S, Val(4), CartesianIndex(v, i, j, 1, h))
    end

    function Base.setindex!(field_values::CC.DataLayouts.VIJFH{S}, val, v::Integer) where {S}
        i,j,_,_,h = size(field_values) # f will be 1, and we replace i...
        # i, j, _, v, h = I.I
        @inbounds CC.DataLayouts.set_struct!(parent(field_values), convert(S, val), Val(4), CartesianIndex(v, i, j, 1, h))
    end
end

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
function FieldFromNamedTuple(space, initial_conditions::Function, ::Type{FT}, params...; calibrate_io::Bool) where {FT}
    local_geometry = CC.Fields.local_geometry_field(space)
    return initial_conditions.(FT, local_geometry, params..., Val{calibrate_io}()) # Val so it can be dispatched on at compile time
end

function Base.cumsum(field::CC.Fields.FiniteDifferenceField)
    Base.cumsum(parent(CC.Spaces.weighted_jacobian(axes(field))) .* vec(field), dims = 1)
end
function Base.cumsum!(fieldout::CC.Fields.FiniteDifferenceField, fieldin::CC.Fields.FiniteDifferenceField)
    Base.cumsum!(fieldout, parent(CC.Spaces.weighted_jacobian(axes(fieldin))) .* vec(fieldin), dims = 1)
end

get_Δz(field::CC.Fields.FiniteDifferenceField) = parent(CC.Spaces.weighted_jacobian(axes(field)))

# TODO: move these things into ClimaCore

isa_center_space(space) = false
isa_center_space(::CC.Spaces.CenterFiniteDifferenceSpace) = true
isa_center_space(::CC.Spaces.CenterExtrudedFiniteDifferenceSpace) = true

isa_face_space(space) = false
isa_face_space(::CC.Spaces.FaceFiniteDifferenceSpace) = true
isa_face_space(::CC.Spaces.FaceExtrudedFiniteDifferenceSpace) = true

const CallableZType = Union{Function, Dierckx.Spline1D} # SSCF.Fast1DLinearInterpolant is only defined in driver rn... I dont think we need it though

function set_z!(field::CC.Fields.Field, u::CallableZType = identity, v::CallableZType = identity)
    z = CC.Fields.coordinate_field(axes(field)).z
    @. field = CCG.Covariant12Vector(CCG.UVVector(u(z), v(z)))
end

function set_z!(field::CC.Fields.Field, u::FT, v::FT) where {FT <: Real}
    lg = CC.Fields.local_geometry_field(axes(field))
    uconst(_) = u # uconst(coord) = u
    vconst(_) = v # vconst(coord) = v
    @. field = CCG.Covariant12Vector(CCG.UVVector(uconst(lg), vconst(lg)))
end

function set_z!(k::Cent, field::CC.Fields.Field, u::FT, v::FT) where {FT <: Real}
    lg = CC.Fields.local_geometry_field(axes(field))
    field[k] = CCG.Covariant12Vector(CCG.UVVector(u, v), lg[k])
end
