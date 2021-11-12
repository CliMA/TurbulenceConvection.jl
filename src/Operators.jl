
abstract type AbstractBCTag end
struct BottomBCTag <: AbstractBCTag end
struct TopBCTag <: AbstractBCTag end
struct InteriorTag <: AbstractBCTag end

abstract type AbstractBC end
struct UseBoundaryValue <: AbstractBC end
struct NoBCGivenError <: AbstractBC end
struct SetZeroGradient <: AbstractBC end
struct SetValue{FT} <: AbstractBC
    value::FT
end
struct SetGradient{FT} <: AbstractBC
    value::FT
end
struct Extrapolate <: AbstractBC end
struct FreeBoundary <: AbstractBC end # when no BC is used (one-sided derivative at surface that takes first and second interior points)

∇f2c(f_dual::SA.SVector, grid::Grid, k; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    ∇f2c(f_dual, grid, k, bottom, top)

function ∇f2c(f_dual::SA.SVector, grid::Grid, k, bottom::AbstractBC, top::AbstractBC)
    if is_surface_center(grid, k)
        return ∇f2c(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return ∇f2c(f_dual, grid, TopBCTag(), top)
    else
        return ∇f2c(f_dual, grid, InteriorTag())
    end
end
∇f2c(f::SA.SVector, grid::Grid, k, ::UseBoundaryValue, top::UseBoundaryValue) = ∇f2c(f, grid, InteriorTag())
∇f2c(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.Δzi
∇f2c(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * grid.Δzi
∇f2c(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
∇f2c(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (f[2] - bc.value) * grid.Δzi
∇f2c(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

∇c2f(f_dual::SA.SVector, grid::Grid, k; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    ∇c2f(f_dual, grid, k, bottom, top)

function ∇c2f(f_dual::SA.SVector, grid::Grid, k, bottom::AbstractBC, top::AbstractBC)
    if is_surface_face(grid, k)
        return ∇c2f(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return ∇c2f(f_dual, grid, TopBCTag(), top)
    else
        return ∇c2f(f_dual, grid, InteriorTag())
    end
end
∇c2f(f::SA.SVector, grid::Grid, ::Int, ::UseBoundaryValue, top::UseBoundaryValue) = ∇c2f(f, grid, InteriorTag())
∇c2f(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.Δzi
∇c2f(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * grid.Δzi * 2.0
∇c2f(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
∇c2f(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (f[1] - bc.value) * grid.Δzi * 2.0
∇c2f(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

# Actually, this method is needed for rain (upwind is in reverse direction due to rain)
function c∇_downwind(f_dual::SA.SVector, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_center(grid, k)
        return c∇_downwind(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return c∇_downwind(f_dual, grid, TopBCTag(), top)
    else
        return c∇_downwind(f_dual, grid, InteriorTag())
    end
end
c∇_downwind(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.Δzi
c∇_downwind(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * (grid.Δzi * 2)
# TODO: this is a crud approximation, as we're specifying what should be the derivative
# at the boundary, and we're taking this as the derivative at the first interior at the
# top of the domain.
c∇_downwind(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
c∇_downwind(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::FreeBoundary) = (f[2] - f[1]) * grid.Δzi # don't use BC info
# TODO: this is a crud approximation, as we're specifying what should be the derivative
# at the boundary, and we're taking this as the derivative at the first interior at the
# top of the domain.
c∇_downwind(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

c∇_upwind(f, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    c∇_upwind(ccut_upwind(f, grid, k), grid, k; bottom, top)

function c∇_upwind(f_dual::SA.SVector, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_center(grid, k)
        return c∇_upwind(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return c∇_upwind(f_dual, grid, TopBCTag(), top)
    else
        return c∇_upwind(f_dual, grid, InteriorTag())
    end
end
c∇_upwind(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.Δzi
c∇_upwind(f::SA.SVector, grid::Grid, ::TopBCTag, bc::FreeBoundary) = (f[2] - f[1]) * grid.Δzi
c∇_upwind(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
c∇_upwind(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (f[1] - bc.value) * (grid.Δzi * 2)
c∇_upwind(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

function f∇_onesided(f_dual::SA.SVector, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_face(grid, k)
        return f∇_onesided(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return f∇_onesided(f_dual, grid, TopBCTag(), top)
    else
        return f∇_onesided(f_dual, grid, InteriorTag())
    end
end
f∇_onesided(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.Δzi
f∇_onesided(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * (grid.Δzi * 2)
f∇_onesided(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
f∇_onesided(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::FreeBoundary) = (f[2] - f[1]) * grid.Δzi # don't use BC info
f∇_onesided(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

# Used when traversing cell faces

interpc2f(f, grid::Grid, k::CCO.PlusHalf; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    interpc2f(dual_centers(f, grid, k), grid, k; bottom, top)

interpc2f(f, grid::Grid, k::CCO.PlusHalf, i_up::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    interpc2f(dual_centers(f, grid, k, i_up), grid, k; bottom, top)

function interpc2f(f_dual::SA.SVector, grid::Grid, k::CCO.PlusHalf; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_face(grid, k)
        return interpc2f(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return interpc2f(f_dual, grid, TopBCTag(), top)
    else
        return interpc2f(f_dual, grid, InteriorTag())
    end
end
interpc2f(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[1] + f[2]) / 2
interpc2f(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = bc.value
interpc2f(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetZeroGradient) = f[1]
interpc2f(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = bc.value

# Used when traversing cell centers
interpf2c(f, grid::Grid, k; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(dual_faces(f, grid, k), grid, k, bottom, top)
interpf2c(f, grid::Grid, k, i_up::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(dual_faces(f, grid, k, i_up), grid, k, bottom, top)

interpf2c(f::SA.SVector, grid::Grid, k; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(f, grid, k, bottom, top)
interpf2c(f::SA.SVector, grid::Grid, k, i_up::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(f, grid, k, bottom, top)

function interpf2c(f_dual::SA.SVector, grid::Grid, k, bottom::AbstractBC, top::AbstractBC)
    if is_surface_center(grid, k)
        return interpf2c(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return interpf2c(f_dual, grid, TopBCTag(), top)
    else
        return interpf2c(f_dual, grid, InteriorTag())
    end
end
interpf2c(f::SA.SVector, grid::Grid, k, ::UseBoundaryValue, top::UseBoundaryValue) = interpf2c(f, grid, InteriorTag())
interpf2c(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[1] + f[2]) / 2
interpf2c(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (f[1] + bc.value) / 2
interpf2c(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (bc.value + f[2]) / 2

#####
##### ∇(center data)
#####

c∇(f, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError()) = c∇(ccut(f, grid, k), grid, k; bottom, top)

function c∇(f_cut::SA.SVector, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_center(grid, k)
        return c∇(f_cut, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return c∇(f_cut, grid, TopBCTag(), top)
    else
        return c∇(f_cut, grid, InteriorTag())
    end
end
c∇(f::SA.SVector, grid::Grid, ::AbstractBCTag, ::NoBCGivenError) = error("No BC given")
function c∇(f::SA.SVector, grid::Grid, ::InteriorTag)
    @assert length(f) == 3
    f_dual⁺ = SA.SVector(f[2], f[3])
    f_dual⁻ = SA.SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = cg+ci => cg = 2fb-ci
    f_dual⁺ = SA.SVector(f[2], 2 * bc.value - f[2])
    f_dual⁻ = SA.SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = cg+ci => cg = 2fb-ci
    f_dual⁺ = SA.SVector(f[1], f[2])
    f_dual⁻ = SA.SVector(2 * bc.value - f[1], f[1])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient)
    @assert length(f) == 2
    f_dual⁻ = SA.SVector(f[1], f[2])
    return (bc.value + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient)
    @assert length(f) == 2
    f_dual⁺ = SA.SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + bc.value) / 2
end
function c∇(f::SA.SVector, grid::Grid, ::TopBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2ci = cg+cii => cg = 2ci-cii. Note: f[3] not used
    f_dual⁺ = SA.SVector(f[2], 2 * f[2] - f[1])
    f_dual⁻ = SA.SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SA.SVector, grid::Grid, ::BottomBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2ci = cg+cii => cg = 2ci-cii. Note: f[1] not used
    f_dual⁺ = SA.SVector(f[1], f[2])
    f_dual⁻ = SA.SVector(2 * f[1] - f[2], f[1])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end

#####
##### ∇(face data)
#####

f∇(f, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError()) = f∇(fcut(f, grid, k), grid, k; bottom, top)

function f∇(f_cut::SA.SVector, grid::Grid, k; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_face(grid, k)
        return f∇(f_cut, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return f∇(f_cut, grid, TopBCTag(), top)
    else
        return f∇(f_cut, grid, InteriorTag())
    end
end
f∇(f::SA.SVector, grid::Grid, ::AbstractBCTag, ::NoBCGivenError) = error("No BC given")
function f∇(f::SA.SVector, grid::Grid, ::InteriorTag)
    @assert length(f) == 3
    f_dual⁺ = SA.SVector(f[2], f[3])
    f_dual⁻ = SA.SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SA.SVector(bc.value, 2 * bc.value - f[1])
    f_dual⁻ = SA.SVector(f[1], bc.value)
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SA.SVector(bc.value, f[2])
    f_dual⁻ = SA.SVector(2 * bc.value - f[2], bc.value)
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient)
    @assert length(f) == 2
    return bc.value
end
function f∇(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient)
    @assert length(f) == 2
    return bc.value
end
function f∇(f::SA.SVector, grid::Grid, ::TopBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SA.SVector(f[2], 2 * f[2] - f[1])
    f_dual⁻ = SA.SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SA.SVector, grid::Grid, ::BottomBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SA.SVector(f[1], f[2])
    f_dual⁻ = SA.SVector(2 * f[1] - f[2], f[1])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end

#####
##### Generic functions
#####

function ∇_staggered(f::SA.SVector, grid::Grid)
    @assert length(f) == 2
    return (f[2] - f[1]) * grid.Δzi
end

"""
    ccut

Used when
     - traversing cell centers
     - grabbing centered center stencils
"""
function ccut(f, grid, k::Cent)
    if is_surface_center(grid, k)
        return SA.SVector(f[k], f[k + 1])
    elseif is_toa_center(grid, k)
        return SA.SVector(f[k - 1], f[k])
    else
        return SA.SVector(f[k - 1], f[k], f[k + 1])
    end
end

"""
    fcut

Used when
     - traversing cell faces
     - grabbing centered face stencils
"""
function fcut(f, grid, k::CCO.PlusHalf)
    if is_surface_face(grid, k)
        return SA.SVector(f[k], f[k + 1])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[k - 1], f[k])
    else
        return SA.SVector(f[k - 1], f[k], f[k + 1])
    end
end

"""
    ccut_downwind

Used when
     - traversing cell centers
     - grabbing one-sided (downwind) stencil of cell center `k` and cell center `k+1`

This is needed for "upwinding" rain, which travels _down_ (hence the direction change).
"""
function ccut_downwind(f, grid, k::Cent)
    if is_toa_center(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k], f[k + 1])
    end
end

"""
    ccut_upwind

Used when
     - traversing cell centers
     - grabbing one-sided (upwind) stencil of cell center `k` and cell center `k-1`
"""
function ccut_upwind(f, grid, k::Cent)
    if is_surface_center(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k - 1], f[k])
    end
end

"""
    fcut_upwind

Used when
     - traversing cell faces
     - grabbing one-sided (upwind) stencil of cell face `k` and cell face `k-1`
"""
function fcut_upwind(f, grid, k)
    if is_surface_face(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k - 1], f[k])
    end
end

"""
    daul_c2f_upwind

Used when
     - traversing cell faces
     - grabbing _interpolated_ one-sided (upwind) stencil of cell face `k` and cell face `k-1`
"""
function daul_c2f_upwind(f, grid, k::CCO.PlusHalf; bottom::SetValue, top::SetZeroGradient)
    kc = Cent(k.i)
    if is_toa_face(grid, k)
        return SA.SVector((f[kc - 2] + f[kc - 1]) / 2, (f[kc - 2] + f[kc - 1]) / 2)
    elseif is_surface_face(grid, k) # never actually called
        error("Uncaught case")
    elseif is_surface_face(grid, k - 1)
        return SA.SVector(bottom.value, (f[kc - 1] + f[kc]) / 2)
    else
        return SA.SVector((f[kc - 2] + f[kc - 1]) / 2, (f[kc - 1] + f[kc]) / 2)
    end
end

"""
    daul_f2c_upwind

Used when
     - traversing cell centers
     - grabbing _interpolated_ one-sided (upwind) stencil of cell center `k` and cell center `k-1`
"""
function daul_f2c_upwind(f, grid, k::Cent)
    kf = CCO.PlusHalf(k.i)
    if is_surface_center(grid, k)
        return SA.SVector((f[kf] + f[kf + 1]) / 2)
    else
        return SA.SVector((f[kf - 1] + f[kf]) / 2, (f[kf] + f[kf + 1]) / 2)
    end
end

"""
    dual_faces

Used when
     - traversing cell centers
     - grabbing stencil of 2 neighboring cell faces
"""
function dual_faces(f, grid, k::Cent)
    kf = CCO.PlusHalf(k.i)
    SA.SVector(f[kf], f[kf + 1])
end
function dual_faces(f, grid, k::Cent, i_up::Int)
    kf = CCO.PlusHalf(k.i)
    SA.SVector(f[i_up, kf], f[i_up, kf + 1])
end

"""
    dual_centers

Used when
     - traversing cell faces
     - grabbing stencil of 2 neighboring cell centers
"""
function dual_centers(f, grid, k::CCO.PlusHalf)
    if is_surface_face(grid, k)
        return SA.SVector(f[Cent(k.i)])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[Cent(k.i - 1)])
    else
        return SA.SVector(f[Cent(k.i - 1)], f[Cent(k.i)])
    end
end
function dual_centers(f, grid, k::CCO.PlusHalf, i_up::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[i_up, Cent(k.i)])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[i_up, Cent(k.i - 1)])
    else
        return SA.SVector(f[i_up, Cent(k.i - 1)], f[i_up, Cent(k.i)])
    end
end


#####
##### Masked operators
#####

"""
    SubMasks

A container for collecting a tuple of
index ranges corresponding to a mask
where `mask == true`.
"""
struct SubMasks{N, T, M, R <: NTuple{N, T}}
    mask::M
    ranges::R
    function SubMasks(mask, min_range_len = 2)
        ranges = UnitRange{Int64}[]
        iter = 0
        mask_offset = 0
        mask_len = length(mask)
        # TODO: this algorithm can be improved, at least, by
        # searching through mask_remaining, instead of the
        # entire mask for every submask.
        while true
            iter > mask_len && break # safety net
            mask_start = findfirst(i -> mask[i] == 1 && i > mask_offset, 1:mask_len)
            isnothing(mask_start) && break
            mask_stop = findfirst(i -> mask[i] == 0 && i > mask_start, 1:mask_len)
            if isnothing(mask_stop)
                R = mask_start:mask_len
                if length(R) >= min_range_len
                    push!(ranges, R)
                end
                break
            else
                R = mask_start:(mask_stop - 1)
                if length(R) >= min_range_len
                    push!(ranges, R)
                end
            end
            mask_offset = mask_stop
            iter += 1
        end
        ranges = Tuple(ranges)
        M = typeof(mask)
        N = length(ranges)
        T = UnitRange{Int64}
        R = typeof(ranges)
        return new{N, T, M, R}(mask, ranges)
    end
end

function shrink_mask(mask)
    return map(enumerate(mask)) do (i, m)
        if i == 1 || i == length(mask)
            m
        elseif m == 1 && mask[i - 1] == 0
            m = 0
        elseif m == 1 && mask[i + 1] == 0
            m = 0
        else
            m
        end
    end
end

Base.eltype(::Type{SDM}) where {N, T, SDM <: SubMasks{N, T}} = T
Base.length(sm::SubMasks{N}) where {N} = N
Base.iterate(sm::SubMasks{N}, state = 1) where {N} = state > N ? nothing : (sm.ranges[state], state + 1)

"""
    subdomain_field(
        field_in,
        space_in::CC.Spaces.CenterFiniteDifferenceSpace,
        ind
    )

Allocates and returns a cell-centered field in the
index range `ind` which is equal to the input field
`field_in`.
"""
function subdomain_field(field_in, space_in::CC.Spaces.CenterFiniteDifferenceSpace, ind)
    FT = eltype(field_in)
    fs_in = CC.Spaces.FaceFiniteDifferenceSpace(space_in)
    zf_in = CC.Fields.coordinate_field(fs_in)
    zc_in = CC.Fields.coordinate_field(fs_in)
    nelems = length(vec(zc_in)[ind])
    zmin = minimum(vec(zf_in)[ind])
    zmax = maximum(vec(zf_in)[ind])
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(zmin),
        CC.Geometry.ZPoint{FT}(zmax);
        boundary_tags = (:bottom, :top),
    )
    mesh = CC.Meshes.IntervalMesh(domain; nelems = nelems)
    cs = CC.Spaces.CenterFiniteDifferenceSpace(mesh)
    cf = CC.Fields.coordinate_field(cs)
    sub_field = sin.(cf.z)
    parent(sub_field) .= parent(field_in)[ind]
    return sub_field
end

"""
    masked_interpolate!(
        face::CC.Fields.FaceFiniteDifferenceField,
        center::CC.Fields.CenterFiniteDifferenceField,
        grid::Grid,
        sub_masks::SubMasks,
    )

Interpolate a cell-centered to cell faces, while avoiding
the use of the cell-centered field where `mask == 0`.
"""
function masked_interpolate!(
    face::CC.Fields.FaceFiniteDifferenceField,
    center::CC.Fields.CenterFiniteDifferenceField,
    grid::Grid,
    sub_masks::SubMasks,
)
    If = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    for sm in sub_masks
        sub_field = subdomain_field(center, grid.cs, sm)
        parent(face)[(sm.start):(sm.stop + 1)] .= vec(If.(sub_field))

        # Linearly extrapolate, instead of zero-th order extrapolate:
        # 2ci = fb+fi, => fb = 2ci - fi
        # TODO: add linear extrapolation support in ClimaCore
        parent(face)[sm.start] = 2 * parent(sub_field)[1] - parent(face)[sm.start + 1]
        parent(face)[sm.stop + 1] = 2 * parent(sub_field)[end] - parent(face)[sm.stop]
    end
end
