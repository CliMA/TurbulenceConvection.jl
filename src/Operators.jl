
abstract type AbstractBCTag end
struct BottomBCTag <: AbstractBCTag end
struct TopBCTag <: AbstractBCTag end
struct InteriorTag <: AbstractBCTag end

struct NoBCGivenError end
struct SetValue{FT}
    value::FT
end
struct SetGradient{FT}
    value::FT
end
struct Extrapolate end

function ∇f2c(f, grid::Grid, k::Int)
    return (f[k] - f[k - 1]) * grid.dzi
end

function ∇c2f(f, grid::Grid, k::Int)
    return (f[k + 1] - f[k]) * grid.dzi
end

function ∇_upwind(f, grid::Grid, k::Int)
    return (f[k + 1] - f[k]) * grid.dzi
end

function interpc2f(f, grid::Grid, k::Int)
    return 0.5 * (f[k + 1] + f[k])
end

function interpc2f(f, grid::Grid, k::Int, i_up::Int)
    return 0.5 * (f[i_up, k + 1] + f[i_up, k])
end

function interpf2c(f, grid::Grid, k::Int)
    return 0.5 * (f[k] + f[k - 1])
end

function interpf2c(f, grid::Grid, k::Int, i_up::Int)
    return 0.5 * (f[i_up, k] + f[i_up, k - 1])
end

#####
##### ∇(center data)
#####

c∇(f, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError()) = c∇(cut(f, k), grid, k; bottom, top)

function c∇(f_cut::SVector, grid::Grid, k::Int; bottom = NoBCGivenError(), top = NoBCGivenError())
    gw = grid.gw
    if k == gw # NOTE: 0-based indexing (k = 3 for 1-based indexing) bottom boundary
        return c∇(f_cut, grid, BottomBCTag(), bottom)
    elseif k == grid.nzg - gw - 1 # NOTE: 0-based indexing
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
    @assert length(f) == 3
    # 2fb = cg+ci => cg = 2fb-ci
    f_dual⁺ = SVector(f[2], 2 * bc.value - f[2])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::BottomBCTag, bc::SetValue)
    @assert length(f) == 3
    # 2fb = cg+ci => cg = 2fb-ci
    f_dual⁺ = SVector(f[2], f[3])
    f_dual⁻ = SVector(2 * bc.value - f[2], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::TopBCTag, bc::SetGradient)
    @assert length(f) == 3
    f_dual⁺ = SVector(f[2], f[3])
    f_dual⁻ = SVector(f[1], f[2])
    return (bc.value + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::BottomBCTag, bc::SetGradient)
    @assert length(f) == 3
    f_dual⁺ = SVector(f[2], f[3])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + bc.value) / 2
end
function c∇(f::SVector, grid::Grid, ::TopBCTag, ::Extrapolate)
    @assert length(f) == 3
    # 2ci = cg+cii => cg = 2ci-cii. Note: f[3] not used
    f_dual⁺ = SVector(f[2], 2 * f[2] - f[1])
    f_dual⁻ = SVector(f[1], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end
function c∇(f::SVector, grid::Grid, ::BottomBCTag, ::Extrapolate)
    @assert length(f) == 3
    # 2ci = cg+cii => cg = 2ci-cii. Note: f[1] not used
    f_dual⁺ = SVector(f[2], f[3])
    f_dual⁻ = SVector(2 * f[2] - f[3], f[2])
    return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
end

#####
##### Generic functions
#####

function ∇_staggered(f::SVector, grid::Grid)
    @assert length(f) == 2
    return (f[2] - f[1]) * grid.dzi
end

# A 3-point field stencil
function cut(f::AbstractVector, k::Int)
    return SVector(f[k - 1], f[k], f[k + 1])
end

# A 3-point field stencil for updraft variables
function cut(f::AbstractMatrix, k::Int, i_up::Int)
    return SVector(f[i_up, k - 1], f[i_up, k], f[i_up, k + 1])
end

function ∇_collocated(f::SVector, grid::Grid)
    @assert length(f) == 3
    return (f[3] - f[1]) * 0.5 * grid.dzi
end

# TODO: use this implementation
# function ∇_collocated(f::SVector, grid::Grid)
#     @assert length(f) == 3
#     f_dual⁺ = SVector(f[2], f[3])
#     f_dual⁻ = SVector(f[1], f[2])
#     return (∇_staggered(f_dual⁺, grid) + ∇_staggered(f_dual⁻, grid)) / 2
# end
