
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
##### Generic functions
#####

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
##### Helpers for masked operators
#####

"""
    shrink_mask(mask)

Shrinks the subdomains where `mask == 1`.

Example:

```julia
using Test
mask = Bool[0, 0, 0, 1, 1, 1, 0, 0, 1, 1]
shrunken_mask = Bool[0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
@test shrink_mask(mask) == shrunken_mask
```
"""
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
