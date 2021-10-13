
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
##### advection operators
#####

function upwind_advection_area(ρ0_half, a_up, w_up, grid, k)
    ρ_0_cut = ccut_upwind(ρ0_half, grid, k)
    a_up_cut = ccut_upwind(a_up, grid, k)
    w_up_cut = daul_f2c_upwind(w_up, grid, k)
    m_cut = ρ_0_cut .* a_up_cut .* w_up_cut
    ∇m = c∇_upwind(m_cut, grid, k; bottom = SetValue(0), top = SetGradient(0))
    return -∇m / ρ0_half[k]
end

function upwind_advection_velocity(ρ0, a_up, w_up, grid, k; a_up_bcs)
    a_dual = daul_c2f_upwind(a_up, grid, k; a_up_bcs...)
    ρ_0_dual = fcut_upwind(ρ0, grid, k)
    w_up_dual = fcut_upwind(w_up, grid, k)
    adv_dual = a_dual .* ρ_0_dual .* w_up_dual .* w_up_dual
    ∇ρaw = f∇_onesided(adv_dual, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
    return ∇ρaw
end

function upwind_advection_scalar(ρ0_half, a_up, w_up, var, grid, k)
    ρ_0_cut = ccut_upwind(ρ0_half, grid, k)
    a_up_cut = ccut_upwind(a_up, grid, k)
    w_up_cut = daul_f2c_upwind(w_up, grid, k)
    var_cut = ccut_upwind(var, grid, k)
    m_cut = ρ_0_cut .* a_up_cut .* w_up_cut .* var_cut
    ∇m = c∇_upwind(m_cut, grid, k; bottom = SetValue(0), top = SetGradient(0))
    return ∇m
end

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
function ccut(f, grid, k)
    if is_surface_center(grid, k)
        return SA.SVector(f[k], f[k + 1])
    elseif is_toa_center(grid, k)
        return SA.SVector(f[k - 1], f[k])
    else
        return SA.SVector(f[k - 1], f[k], f[k + 1])
    end
end
function ccut(f, grid, k, i_up::Int)
    if is_surface_center(grid, k)
        return SA.SVector(f[i_up, k], f[i_up, k + 1])
    elseif is_toa_center(grid, k)
        return SA.SVector(f[i_up, k - 1], f[i_up, k])
    else
        return SA.SVector(f[i_up, k - 1], f[i_up, k], f[i_up, k + 1])
    end
end

"""
    fcut

Used when
     - traversing cell faces
     - grabbing centered face stencils
"""
function fcut(f, grid, k)
    if is_surface_face(grid, k)
        return SA.SVector(f[k], f[k + 1])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[k - 1], f[k])
    else
        return SA.SVector(f[k - 1], f[k], f[k + 1])
    end
end
function fcut(f, grid, k, i_up::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[i_up, k], f[i_up, k + 1])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[i_up, k - 1], f[i_up, k])
    else
        return SA.SVector(f[i_up, k - 1], f[i_up, k], f[i_up, k + 1])
    end
end

"""
    ccut_downwind

Used when
     - traversing cell centers
     - grabbing one-sided (downwind) stencil of cell center `k` and cell center `k+1`

This is needed for "upwinding" rain, which travels _down_ (hence the direction change).
"""
function ccut_downwind(f, grid, k)
    if is_toa_center(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k], f[k + 1])
    end
end
function ccut_downwind(f, grid, k, i_up::Int)
    if is_toa_center(grid, k)
        return SA.SVector(f[i_up, k])
    else
        return SA.SVector(f[i_up, k], f[i_up, k + 1])
    end
end

"""
    ccut_upwind

Used when
     - traversing cell centers
     - grabbing one-sided (upwind) stencil of cell center `k` and cell center `k-1`
"""
function ccut_upwind(f, grid, k)
    if is_surface_center(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k - 1], f[k])
    end
end
function ccut_upwind(f, grid, k, i_up::Int)
    if is_surface_center(grid, k)
        return SA.SVector(f[i_up, k])
    else
        return SA.SVector(f[i_up, k - 1], f[i_up, k])
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
function fcut_upwind(f, grid, k, i_up::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[i_up, k])
    else
        return SA.SVector(f[i_up, k - 1], f[i_up, k])
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
function daul_c2f_upwind(f, grid, k::CCO.PlusHalf, i_up::Int; bottom::SetValue, top::SetZeroGradient)
    kc = Cent(k.i)
    if is_toa_face(grid, k)
        return SA.SVector((f[i_up, kc - 2] + f[i_up, kc - 1]) / 2, (f[i_up, kc - 2] + f[i_up, kc - 1]) / 2)
    elseif is_surface_face(grid, k) # never actually called
        error("Uncaught case")
    elseif is_surface_face(grid, kc - 1)
        return SA.SVector(bottom.value, (f[i_up, kc - 1] + f[i_up, kc]) / 2)
    else
        return SA.SVector((f[i_up, kc - 2] + f[i_up, kc - 1]) / 2, (f[i_up, kc - 1] + f[i_up, kc]) / 2)
    end
end

"""
    daul_f2c_upwind

Used when
     - traversing cell centers
     - grabbing _interpolated_ one-sided (upwind) stencil of cell center `k` and cell center `k-1`
"""
function daul_f2c_upwind(f, grid, k)
    if is_surface_center(grid, k)
        return SA.SVector((f[k] + f[k + 1]) / 2)
    else
        return SA.SVector((f[k - 1] + f[k]) / 2, (f[k] + f[k + 1]) / 2)
    end
end
function daul_f2c_upwind(f, grid, k, i_up::Int)
    if is_surface_center(grid, k)
        return SA.SVector((f[i_up, k] + f[i_up, k + 1]) / 2)
    else
        return SA.SVector((f[i_up, k - 1] + f[i_up, k]) / 2, (f[i_up, k] + f[i_up, k + 1]) / 2)
    end
end

"""
    dual_faces

Used when
     - traversing cell centers
     - grabbing stencil of 2 neighboring cell faces
"""
dual_faces(f, grid, k) = SA.SVector(f[k], f[k + 1])
dual_faces(f, grid, k, i_up::Int) = SA.SVector(f[i_up, k], f[i_up, k + 1])

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
##### Implicit operators
#####

#=
We consider the time-discretized equation:

``
(ρaeϕⁿ⁺¹ - ρaeϕⁿ)/Δt = D ϕⁿ⁺¹ + Sⁿ
D = ∂_z (ρaeK ∂_z)
``

Which leads to:

``
ρaeϕⁿ⁺¹ - ρaeϕⁿ = Δt D ϕⁿ⁺¹ + Δt Sⁿ
ρaeϕⁿ⁺¹ - Δt D ϕⁿ⁺¹ = (ϕⁿ + Δt Sⁿ)
``

We relax this to

``
(I - Δt D/ρaeⁿ) ϕⁿ⁺¹ = (ϕⁿ + Δt Sⁿ)/ρaeⁿ
``

Let `A = (I - Δt D)` and `b = ϕⁿ + Δt Sⁿ`, and we have

A ϕⁿ⁺¹ = b

This function constructs and returns the matrix `A`,
where `D = ∂_z (ρaeK ∂_z)`

=#
function construct_tridiag_diffusion_gm(grid::Grid, dt, ρ_ae_K, ρ_0, ae)
    a = center_field(grid) # for tridiag solver
    b = center_field(grid) # for tridiag solver
    c = center_field(grid) # for tridiag solver
    Δzi = grid.Δzi
    @inbounds for k in real_center_indices(grid)
        X = ρ_0[k] * ae[k] / dt
        Y = ρ_ae_K[k + 1] * Δzi * Δzi
        Z = ρ_ae_K[k] * Δzi * Δzi
        if is_surface_center(grid, k)
            Z = 0.0
        elseif is_toa_center(grid, k)
            Y = 0.0
        end
        a[k] = -Z / X
        b[k] = 1.0 + Y / X + Z / X
        c[k] = -Y / X
    end
    A = LinearAlgebra.Tridiagonal(a[2:end], vec(b), c[1:(end - 1)])
    return A
end

tridiag_solve(b_rhs, A) = A \ b_rhs

#= TODO: clean this up somehow!

construct_tridiag_diffusion_en configures the tridiagonal
matrix for solving 2nd order environment variables where
some terms are treated implicitly. Here is a list of all
the terms in the matrix (and the overall equation solved):

    N_a terms = 1 (diffusion)
    N_diagonal terms = 6 (unsteady+upwind_advection+2diffusion+entr_detr+dissipation)
    N_c terms = 2 (diffusion+upwind_advection)

    ∂_t ρa*tke                              ✓ unsteady
      + ∂_z(ρaw*tke) =                      ✓ advection
      + ρaK*(∂²_z(u)+∂²_z(v)+∂²_z(w̄))
      + ρawΣᵢ(εᵢ(wⱼ-w₀)²-δ₀*tke)            ✓ entr_detr
      + ∂_z(ρa₀K ∂_z(tke))                  ✓ diffusion
      + ρa₀*w̅₀b̅₀
      - a₀(u⋅∇p)
      + ρa₀D                                ✓ dissipation
=#
function construct_tridiag_diffusion_en(
    grid::Grid,
    param_set::APS,
    state,
    TS,
    KM,
    KH,
    a_up_bulk,
    up_prog,
    w_up,
    w_en,
    tke_en,
    n_updrafts::Int,
    minimum_area::Float64,
    pressure_plume_spacing::Vector,
    frac_turb_entr,
    entr_sc,
    mixing_length,
    is_tke,
)

    c_d = CPEDMF.c_d(param_set)
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    Δzi = grid.Δzi
    dti = TS.dti
    a = center_field(grid)
    b = center_field(grid)
    c = center_field(grid)
    ρ0_c = center_ref_state(state).ρ0
    ρ0_f = face_ref_state(state).ρ0

    ae = 1 .- a_up_bulk
    rho_ae_K_m = face_field(grid)
    w_en_c = center_field(grid)
    D_env = 0.0

    aeK = is_tke ? ae .* KM : ae .* KH
    aeK_bcs = (; bottom = SetValue(aeK[kc_surf]), top = SetValue(aeK[kc_toa]))

    @inbounds for k in real_face_indices(grid)
        rho_ae_K_m[k] = interpc2f(aeK, grid, k; aeK_bcs...) * ρ0_f[k]
    end

    @inbounds for k in real_center_indices(grid)
        w_en_c[k] = interpf2c(w_en, grid, k)
    end

    @inbounds for k in real_center_indices(grid)
        D_env = 0.0

        @inbounds for i in 1:n_updrafts
            if up_prog[i].area[k] > minimum_area
                turb_entr = frac_turb_entr[i, k]
                R_up = pressure_plume_spacing[i]
                w_up_c = interpf2c(w_up, grid, k, i)
                D_env += ρ0_c[k] * up_prog[i].area[k] * w_up_c * (entr_sc[i, k] + turb_entr)
            else
                D_env = 0.0
            end
        end

        # TODO: this tridiagonal matrix needs to be re-verified, as it's been pragmatically
        #       modified to not depend on ghost points, and these changes have not been
        #       carefully verified.
        if is_surface_center(grid, k)
            a[k] = 0.0
            b[k] = 1.0
            c[k] = 0.0
        else
            a[k] = (-rho_ae_K_m[k] * Δzi * Δzi)
            b[k] = (
                ρ0_c[k] * ae[k] * dti - ρ0_c[k] * ae[k] * w_en_c[k] * Δzi +
                rho_ae_K_m[k + 1] * Δzi * Δzi +
                rho_ae_K_m[k] * Δzi * Δzi +
                D_env +
                ρ0_c[k] * ae[k] * c_d * sqrt(max(tke_en[k], 0)) / max(mixing_length[k], 1)
            )
            if is_toa_center(grid, k)
                c[k] = 0.0
            else
                c[k] = (ρ0_c[k + 1] * ae[k + 1] * w_en_c[k + 1] * Δzi - rho_ae_K_m[k + 1] * Δzi * Δzi)
            end
        end
    end

    A = LinearAlgebra.Tridiagonal(a[2:end], vec(b), c[1:(end - 1)])
    return A
end
