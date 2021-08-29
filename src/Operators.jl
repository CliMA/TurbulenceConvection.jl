
abstract type AbstractBCTag end
struct BottomBCTag <: AbstractBCTag end
struct TopBCTag <: AbstractBCTag end
struct InteriorTag <: AbstractBCTag end

abstract type AbstractBC end
struct UseBoundaryValue <: AbstractBC end
struct NoBCGivenError <: AbstractBC end
struct SetValue{FT} <: AbstractBC
    value::FT
end
struct SetGradient{FT} <: AbstractBC
    value::FT
end
struct Extrapolate <: AbstractBC end
struct FreeBoundary <: AbstractBC end # when no BC is used (one-sided derivative at surface that takes first and second interior points)

∇f2c(f_dual::SVector, grid::Grid, k::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    ∇f2c(f_dual, grid, k, bottom, top)

function ∇f2c(f_dual::SVector, grid::Grid, k::Int, bottom::AbstractBC, top::AbstractBC)
    if is_surface_face(grid, k - 1)
        return ∇f2c(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return ∇f2c(f_dual, grid, TopBCTag(), top)
    else
        return ∇f2c(f_dual, grid, InteriorTag())
    end
end
∇f2c(f::SVector, grid::Grid, ::Int, ::UseBoundaryValue, top::UseBoundaryValue) = ∇f2c(f, grid, InteriorTag())
∇f2c(f::SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.dzi
∇f2c(f::SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * grid.dzi
∇f2c(f::SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
∇f2c(f::SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (f[2] - bc.value) * grid.dzi
∇f2c(f::SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

function ∇_onesided(f_dual::SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_center(grid, k)
        return ∇_onesided(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return ∇_onesided(f_dual, grid, TopBCTag(), top)
    else
        return ∇_onesided(f_dual, grid, InteriorTag())
    end
end
∇_onesided(f::SVector, grid::Grid, ::InteriorTag) = (f[2] - f[1]) * grid.dzi
∇_onesided(f::SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (bc.value - f[1]) * (grid.dzi / 2)
# TODO: this is a crud approximation, as we're specifying what should be the derivative
# at the boundary, and we're taking this as the derivative at the first interior at the
# top of the domain.
∇_onesided(f::SVector, grid::Grid, ::TopBCTag, bc::SetGradient) = bc.value
∇_onesided(f::SVector, grid::Grid, ::BottomBCTag, bc::FreeBoundary) = (f[2] - f[1]) * grid.dzi # don't use BC info
# TODO: this is a crud approximation, as we're specifying what should be the derivative
# at the boundary, and we're taking this as the derivative at the first interior at the
# top of the domain.
∇_onesided(f::SVector, grid::Grid, ::BottomBCTag, bc::SetGradient) = bc.value

# Used when traversing cell faces

interpc2f(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    interpc2f(dual_centers(f, grid, k), grid, k; bottom, top)

interpc2f(f, grid::Grid, k::Int, i_up::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    interpc2f(dual_centers(f, grid, k, i_up), grid, k; bottom, top)

function interpc2f(f_dual::SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_face(grid, k)
        return interpc2f(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return interpc2f(f_dual, grid, TopBCTag(), top)
    else
        return interpc2f(f_dual, grid, InteriorTag())
    end
end
interpc2f(f::SVector, grid::Grid, ::InteriorTag) = (f[1] + f[2]) / 2
interpc2f(f::SVector, grid::Grid, ::TopBCTag, bc::SetValue) = bc.value
interpc2f(f::SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = bc.value

# Used when traversing cell centers
interpf2c(f, grid::Grid, k::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(dual_faces(f, grid, k), grid, k, bottom, top)

interpf2c(f, grid::Grid, k::Int, i_up::Int; bottom = UseBoundaryValue(), top = UseBoundaryValue()) =
    interpf2c(dual_faces(f, grid, k, i_up), grid, k, bottom, top)

function interpf2c(f_dual::SVector, grid::Grid, k::Int, bottom::AbstractBC, top::AbstractBC)
    if is_surface_center(grid, k)
        return interpf2c(f_dual, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return interpf2c(f_dual, grid, TopBCTag(), top)
    else
        return interpf2c(f_dual, grid, InteriorTag())
    end
end
interpf2c(f::SVector, grid::Grid, ::Int, ::UseBoundaryValue, top::UseBoundaryValue) = interpf2c(f, grid, InteriorTag())
interpf2c(f::SVector, grid::Grid, ::InteriorTag) = (f[1] + f[2]) / 2
interpf2c(f::SVector, grid::Grid, ::TopBCTag, bc::SetValue) = (f[1] + bc.value) / 2
interpf2c(f::SVector, grid::Grid, ::BottomBCTag, bc::SetValue) = (bc.value + f[2]) / 2


# To be deprecated
interp2pt(val1, val2) = 0.5 * (val1 + val2)

#####
##### advection operators
#####

function upwind_advection_area(ρ0_half::Vector{Float64}, a_up::Vector{Float64}, w_up::Vector{Float64}, grid, k)
    whalf_kp = interp2pt(w_up[k - 1], w_up[k])
    whalf_k = interp2pt(w_up[k - 2], w_up[k - 1])

    m_kp = (ρ0_half[k] * a_up[k] * whalf_kp)
    m_k = (ρ0_half[k - 1] * a_up[k - 1] * whalf_k)
    return -grid.dzi * (m_kp - m_k) / ρ0_half[k]
end

function upwind_advection_velocity(ρ0::Vector{Float64}, a_up::Vector{Float64}, w_up::Vector{Float64}, grid, k)
    a_k = interp2pt(a_up[k], a_up[k + 1])
    a_km = interp2pt(a_up[k - 1], a_up[k])
    adv = (ρ0[k] * a_k * w_up[k] * w_up[k] * grid.dzi - ρ0[k - 1] * a_km * w_up[k - 1] * w_up[k - 1] * grid.dzi)
    return adv
end

function upwind_advection_scalar(
    ρ0_half::Vector{Float64},
    a_up::Vector{Float64},
    w_up::Vector{Float64},
    var::Vector{Float64},
    grid,
    k,
)
    m_k = (ρ0_half[k] * a_up[k] * interp2pt(w_up[k - 1], w_up[k]))
    m_km = (ρ0_half[k - 1] * a_up[k - 1] * interp2pt(w_up[k - 2], w_up[k - 1]))
    return (m_k * var[k] - m_km * var[k - 1]) * grid.dzi
end

#####
##### ∇(center data)
#####

c∇(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    c∇(ccut(f, grid, k), grid, k; bottom, top)

function c∇(f_cut::SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_center(grid, k)
        return c∇(f_cut, grid, BottomBCTag(), bottom)
    elseif is_toa_center(grid, k)
        return c∇(f_cut, grid, TopBCTag(), top)
    else
        return c∇(f_cut, grid, InteriorTag())
    end
end
c∇(f::SVector, grid::Grid, ::AbstractBCTag, ::NoBCGivenError) = error("No BC given")
function c∇(f::SVector, grid::Grid, ::InteriorTag)
    @assert length(f) == 3
    f_dual⁺ = SVector(f[2], f[3])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::TopBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = cg+ci => cg = 2fb-ci
    f_dual⁺ = SVector(f[2], 2 * bc.value - f[2])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::BottomBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = cg+ci => cg = 2fb-ci
    f_dual⁺ = SVector(f[1], f[2])
    f_dual⁻ = SVector(2 * bc.value - f[1], f[1])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::TopBCTag, bc::SetGradient)
    @assert length(f) == 2
    f_dual⁻ = SVector(f[1], f[2])
    return (bc.value + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::BottomBCTag, bc::SetGradient)
    @assert length(f) == 2
    f_dual⁺ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + bc.value) / 2
end
function c∇(f::SVector, grid::Grid, ::TopBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2ci = cg+cii => cg = 2ci-cii. Note: f[3] not used
    f_dual⁺ = SVector(f[2], 2 * f[2] - f[1])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::BottomBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2ci = cg+cii => cg = 2ci-cii. Note: f[1] not used
    f_dual⁺ = SVector(f[1], f[2])
    f_dual⁻ = SVector(2 * f[1] - f[2], f[1])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end

#####
##### ∇(face data)
#####

f∇(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) =
    f∇(fcut(f, grid, k), grid, k; bottom, top)

function f∇(f_cut::SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
    if is_surface_face(grid, k)
        return f∇(f_cut, grid, BottomBCTag(), bottom)
    elseif is_toa_face(grid, k)
        return f∇(f_cut, grid, TopBCTag(), top)
    else
        return f∇(f_cut, grid, InteriorTag())
    end
end
f∇(f::SVector, grid::Grid, ::AbstractBCTag, ::NoBCGivenError) = error("No BC given")
function f∇(f::SVector, grid::Grid, ::InteriorTag)
    @assert length(f) == 3
    f_dual⁺ = SVector(f[2], f[3])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SVector, grid::Grid, ::TopBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SVector(bc.value, 2 * bc.value - f[1])
    f_dual⁻ = SVector(f[1], bc.value)
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SVector, grid::Grid, ::BottomBCTag, bc::SetValue)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SVector(bc.value, f[2])
    f_dual⁻ = SVector(2 * bc.value - f[2], bc.value)
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SVector, grid::Grid, ::TopBCTag, bc::SetGradient)
    @assert length(f) == 2
    return bc.value
end
function f∇(f::SVector, grid::Grid, ::BottomBCTag, bc::SetGradient)
    @assert length(f) == 2
    return bc.value
end
function f∇(f::SVector, grid::Grid, ::TopBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SVector(f[2], 2 * f[2] - f[1])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function f∇(f::SVector, grid::Grid, ::BottomBCTag, ::Extrapolate)
    @assert length(f) == 2
    # 2fb = fg+fi => fg = 2fb-fi
    f_dual⁺ = SVector(f[1], f[2])
    f_dual⁻ = SVector(2 * f[1] - f[2], f[1])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end

#####
##### Generic functions
#####

function ∇_staggered(f::SVector, grid::Grid)
    @assert length(f) == 2
    return (f[2] - f[1]) * grid.dzi
end

# A 3-point field stencil for ordinary and updraft variables
function ccut(f::AbstractVector, grid, k::Int)
    if is_surface_center(grid, k)
        return SVector(f[k], f[k + 1])
    elseif is_toa_center(grid, k)
        return SVector(f[k - 1], f[k])
    else
        return SVector(f[k - 1], f[k], f[k + 1])
    end
end
function ccut(f::AbstractMatrix, grid, k::Int, i_up::Int)
    if is_surface_center(grid, k)
        return SVector(f[i_up, k], f[i_up, k + 1])
    elseif is_toa_center(grid, k)
        return SVector(f[i_up, k - 1], f[i_up, k])
    else
        return SVector(f[i_up, k - 1], f[i_up, k], f[i_up, k + 1])
    end
end

function fcut(f::AbstractVector, grid, k::Int)
    if is_surface_face(grid, k)
        return SVector(f[k], f[k + 1])
    elseif is_toa_face(grid, k)
        return SVector(f[k - 1], f[k])
    else
        return SVector(f[k - 1], f[k], f[k + 1])
    end
end
function fcut(f::AbstractMatrix, grid, k::Int, i_up::Int)
    if is_surface_face(grid, k)
        return SVector(f[i_up, k], f[i_up, k + 1])
    elseif is_toa_face(grid, k)
        return SVector(f[i_up, k - 1], f[i_up, k])
    else
        return SVector(f[i_up, k - 1], f[i_up, k], f[i_up, k + 1])
    end
end

# A 2-point field stencil for ordinary and updraft variables
function ccut_onesided(f::AbstractVector, grid, k::Int)
    if is_toa_center(grid, k)
        return SVector(f[k])
    else
        return SVector(f[k], f[k + 1])
    end
end
function ccut_onesided(f::AbstractMatrix, grid, k::Int, i_up::Int)
    if is_toa_center(grid, k)
        return SVector(f[i_up, k])
    else
        return SVector(f[i_up, k], f[i_up, k + 1])
    end
end

# A 2-point field stencil for ordinary and updraft variables
dual_faces(f::AbstractVector, grid, k::Int) = SVector(f[k - 1], f[k])
dual_faces(f::AbstractMatrix, grid, k::Int, i_up::Int) = SVector(f[i_up, k - 1], f[i_up, k])

function dual_centers(f::AbstractVector, grid, k::Int)
    if is_surface_face(grid, k)
        return SVector(f[k + 1])
    elseif is_toa_face(grid, k)
        return SVector(f[k])
    else
        return SVector(f[k], f[k + 1])
    end
end
function dual_centers(f::AbstractMatrix, grid, k::Int, i_up::Int)
    if is_surface_face(grid, k)
        return SVector(f[i_up, k + 1])
    elseif is_toa_face(grid, k)
        return SVector(f[i_up, k])
    else
        return SVector(f[i_up, k], f[i_up, k + 1])
    end
end
