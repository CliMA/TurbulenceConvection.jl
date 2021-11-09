
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

# ∇(center data) for possibly vanishing subdomains.

"""
    c∇_vanishing_subdomain

Used when the vertical gradient of a field must be computed over a conditional subdomain which may vanish
below or above the current cell. If the subdomains vanishes, a default user-defined gradient is returned.

Inputs:
 - f_cut: Slice of field.
 - sd_cut: Slice of subdomain volume fraction.
 - default∇: Gradient used in vanishing subdomains.
"""
function c∇_vanishing_subdomain(
    f_cut::SA.SVector,
    sd_cut::SA.SVector,
    default∇::FT,
    grid::Grid,
    k;
    bottom = NoBCGivenError(),
    top = NoBCGivenError(),
) where {FT <: Real}
    if is_surface_center(grid, k)
        return c∇(f_cut, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return c∇(f_cut, grid, TopBCTag(), top)
    else
        return c∇_vanishing_subdomain(f_cut, sd_cut, default∇, grid, InteriorTag())
    end
end
function c∇_vanishing_subdomain(
    f::SA.SVector,
    sd::SA.SVector,
    default∇::FT,
    grid::Grid,
    ::InteriorTag,
) where {FT <: Real}
    @assert length(f) == 3
    @assert length(sd) == 3
    if sd[1] * sd[3] ≈ FT(0)
        return default∇
    else
        f_dual⁺ = SA.SVector(f[2], f[3])
        f_dual⁻ = SA.SVector(f[1], f[2])
        return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
    end
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
