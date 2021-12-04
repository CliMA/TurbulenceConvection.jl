
abstract type AbstractBCTag end
struct BottomBCTag <: AbstractBCTag end
struct TopBCTag <: AbstractBCTag end
struct InteriorTag <: AbstractBCTag end

abstract type AbstractBC end
struct UseBoundaryValue <: AbstractBC end
struct SetValue{FT} <: AbstractBC
    value::FT
end

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
