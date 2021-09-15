
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

∇f2c(f_dual::SA.SVector, grid::Grid, k::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    ∇f2c(f_dual, grid, k, bottom, top)

function ∇f2c(f_dual::SA.SVector, grid::Grid, k::Int, bottom::AbstractBC, top::AbstractBC)
    if is_surface_center(grid, k)
        return ∇f2c(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return ∇f2c(f_dual, grid, TopBCTag(), top)
    else
        return ∇f2c(f_dual, grid, InteriorTag())
    end
end
∇f2c(f::SA.SVector, grid::Grid, ::Int, ::UseBoundaryValue, top::UseBoundaryValue) = ∇f2c(f, grid, InteriorTag())
∇f2c(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.Δzi
∇f2c(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * grid.Δzi
∇f2c(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
∇f2c(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (f[2] - bc.value) * grid.Δzi
∇f2c(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

# TODO: this should be changed to `c∇` or `c∇_upwind`
# Actually, this method is needed for rain (upwind is in reverse direction due to rain)
function c∇_downwind(f_dual::SA.SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
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

c∇_upwind(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    c∇_upwind(ccut_upwind(f, grid, k), grid, k; bottom, top)

function c∇_upwind(f_dual::SA.SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
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

function f∇_onesided(f_dual::SA.SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
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

interpc2f(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    interpc2f(dual_centers(f, grid, k), grid, k; bottom, top)

interpc2f(f, grid::Grid, k::Int, i_up::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    interpc2f(dual_centers(f, grid, k, i_up), grid, k; bottom, top)

function interpc2f(f_dual::SA.SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
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
interpf2c(f, grid::Grid, k::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(dual_faces(f, grid, k), grid, k, bottom, top)
interpf2c(f, grid::Grid, k::Int, i_up::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(dual_faces(f, grid, k, i_up), grid, k, bottom, top)

interpf2c(f::SA.SVector, grid::Grid, k::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(f, grid, k, bottom, top)
interpf2c(f::SA.SVector, grid::Grid, k::Int, i_up::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(f, grid, k, bottom, top)

function interpf2c(f_dual::SA.SVector, grid::Grid, k::Int, bottom::AbstractBC, top::AbstractBC)
    if is_surface_center(grid, k)
        return interpf2c(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return interpf2c(f_dual, grid, TopBCTag(), top)
    else
        return interpf2c(f_dual, grid, InteriorTag())
    end
end
interpf2c(f::SA.SVector, grid::Grid, ::Int, ::UseBoundaryValue, top::UseBoundaryValue) =
    interpf2c(f, grid, InteriorTag())
interpf2c(f::SA.SVector, grid::Grid, ::InteriorTag) = (f[1] + f[2]) / 2
interpf2c(f::SA.SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (f[1] + bc.value) / 2
interpf2c(f::SA.SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (bc.value + f[2]) / 2

#####
##### advection operators
#####

function upwind_advection_area(ρ0_half::Vector{Float64}, a_up::Vector{Float64}, w_up::Vector{Float64}, grid, k)
    ρ_0_cut = ccut_upwind(ρ0_half, grid, k)
    a_up_cut = ccut_upwind(a_up, grid, k)
    w_up_cut = daul_f2c_upwind(w_up, grid, k)
    m_cut = ρ_0_cut .* a_up_cut .* w_up_cut
    ∇m = c∇_upwind(m_cut, grid, k; bottom = SetValue(0), top = SetGradient(0))
    return -∇m / ρ0_half[k]
end

function upwind_advection_velocity(ρ0::Vector{Float64}, a_up::Vector{Float64}, w_up::Vector{Float64}, grid, k; a_up_bcs)
    a_dual = daul_c2f_upwind(a_up, grid, k; a_up_bcs...)
    ρ_0_dual = fcut_upwind(ρ0, grid, k)
    w_up_dual = fcut_upwind(w_up, grid, k)
    adv_dual = a_dual .* ρ_0_dual .* w_up_dual .* w_up_dual
    ∇ρaw = f∇_onesided(adv_dual, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
    return ∇ρaw
end

function upwind_advection_scalar(
    ρ0_half::Vector{Float64},
    a_up::Vector{Float64},
    w_up::Vector{Float64},
    var::Vector{Float64},
    grid,
    k,
)
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

c∇(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    c∇(ccut(f, grid, k), grid, k; bottom, top)

function c∇(f_cut::SA.SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
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

f∇(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    f∇(fcut(f, grid, k), grid, k; bottom, top)

function f∇(f_cut::SA.SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
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
function ccut(f::AbstractVector, grid, k::Int)
    if is_surface_center(grid, k)
        return SA.SVector(f[k], f[k + 1])
    elseif is_toa_center(grid, k)
        return SA.SVector(f[k - 1], f[k])
    else
        return SA.SVector(f[k - 1], f[k], f[k + 1])
    end
end
function ccut(f::AbstractMatrix, grid, k::Int, i_up::Int)
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
function fcut(f::AbstractVector, grid, k::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[k], f[k + 1])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[k - 1], f[k])
    else
        return SA.SVector(f[k - 1], f[k], f[k + 1])
    end
end
function fcut(f::AbstractMatrix, grid, k::Int, i_up::Int)
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
function ccut_downwind(f::AbstractVector, grid, k::Int)
    if is_toa_center(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k], f[k + 1])
    end
end
function ccut_downwind(f::AbstractMatrix, grid, k::Int, i_up::Int)
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
function ccut_upwind(f::AbstractVector, grid, k::Int)
    if is_surface_center(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k - 1], f[k])
    end
end
function ccut_upwind(f::AbstractMatrix, grid, k::Int, i_up::Int)
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
function fcut_upwind(f::AbstractVector, grid, k::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[k])
    else
        return SA.SVector(f[k - 1], f[k])
    end
end
function fcut_upwind(f::AbstractMatrix, grid, k::Int, i_up::Int)
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
function daul_c2f_upwind(f::AbstractVector, grid, k::Int; bottom::SetValue, top::SetZeroGradient)
    if is_toa_face(grid, k)
        return SA.SVector((f[k - 2] + f[k - 1]) / 2, (f[k - 2] + f[k - 1]) / 2)
    elseif is_surface_face(grid, k) # never actually called
        error("Uncaught case")
    elseif is_surface_face(grid, k - 1)
        return SA.SVector(bottom.value, (f[k - 1] + f[k]) / 2)
    else
        return SA.SVector((f[k - 2] + f[k - 1]) / 2, (f[k - 1] + f[k]) / 2)
    end
end
function daul_c2f_upwind(f::AbstractVector, grid, k::Int, i_up::Int; bottom::SetValue, top::SetZeroGradient)
    if is_toa_face(grid, k)
        return SA.SVector((f[i_up, k - 2] + f[i_up, k - 1]) / 2, (f[i_up, k - 2] + f[i_up, k - 1]) / 2)
    elseif is_surface_face(grid, k) # never actually called
        error("Uncaught case")
    elseif is_surface_face(grid, k - 1)
        return SA.SVector(bottom.value, (f[i_up, k - 1] + f[i_up, k]) / 2)
    else
        return SA.SVector((f[i_up, k - 2] + f[i_up, k - 1]) / 2, (f[i_up, k - 1] + f[i_up, k]) / 2)
    end
end

"""
    daul_f2c_upwind

Used when
     - traversing cell centers
     - grabbing _interpolated_ one-sided (upwind) stencil of cell center `k` and cell center `k-1`
"""
function daul_f2c_upwind(f::AbstractVector, grid, k::Int)
    if is_surface_center(grid, k)
        return SA.SVector((f[k] + f[k + 1]) / 2)
    else
        return SA.SVector((f[k - 1] + f[k]) / 2, (f[k] + f[k + 1]) / 2)
    end
end
function daul_f2c_upwind(f::AbstractMatrix, grid, k::Int, i_up::Int)
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
dual_faces(f::AbstractVector, grid, k::Int) = SA.SVector(f[k], f[k + 1])
dual_faces(f::AbstractMatrix, grid, k::Int, i_up::Int) = SA.SVector(f[i_up, k], f[i_up, k + 1])

"""
    dual_centers

Used when
     - traversing cell faces
     - grabbing stencil of 2 neighboring cell centers
"""
function dual_centers(f::AbstractVector, grid, k::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[k])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[k - 1])
    else
        return SA.SVector(f[k - 1], f[k])
    end
end
function dual_centers(f::AbstractMatrix, grid, k::Int, i_up::Int)
    if is_surface_face(grid, k)
        return SA.SVector(f[i_up, k])
    elseif is_toa_face(grid, k)
        return SA.SVector(f[i_up, k - 1])
    else
        return SA.SVector(f[i_up, k - 1], f[i_up, k])
    end
end
