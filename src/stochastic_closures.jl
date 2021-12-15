
abstract type AbstractEntDet end
struct Entrainment <: AbstractEntDet end
struct Detrainment <: AbstractEntDet end

stochastic_closure(::APS, ::sde_struct{NoneClosureType}, ::AbstractEntDet) = 1
function stochastic_closure(param_set::APS, ::sde_struct{LogNormalClosureType}, term::AbstractEntDet)
    u = lognormal_closure(param_set, term)
    if u < 0
        @warn("Negative stochastic parameter ⟹ negative $(typeof(term)). Setting to 0")
        u = 0
    end
    return u
end
function stochastic_closure(param_set::APS, sde_model::sde_struct{SDEClosureType}, term::AbstractEntDet)
    u = sde_closure(param_set, sde_model, term)
    if u < 0
        @warn("Negative stochastic parameter ⟹ negative $(typeof(term)). Setting to 0")
        u = 0
    end
    return u
end

#
## Lognormal closure functions
#
stoch_lognormal_var(param_set, ::Entrainment) = ICP.stoch_ε_lognormal_var(param_set)
stoch_lognormal_var(param_set, ::Detrainment) = ICP.stoch_δ_lognormal_var(param_set)

function lognormal_closure(param_set::APS, term::AbstractEntDet)
    lognormal_var = stoch_lognormal_var(param_set, term)
    return lognormal_sampler(1.0, lognormal_var)
end

function lognormal_sampler(m, var)
    μ = log(m^2 / √(m^2 + var))
    σ = √(log(1 + var / m^2))
    return rand(Distributions.LogNormal(μ, σ))
end


# Stochastic differential equation

sde_θ(param_set, ::Entrainment) = ICP.sde_ϵ_θ(param_set)
sde_σ(param_set, ::Entrainment) = ICP.sde_ϵ_σ(param_set)
sde_θ(param_set, ::Detrainment) = ICP.sde_δ_θ(param_set)
sde_σ(param_set, ::Detrainment) = ICP.sde_δ_σ(param_set)

function sde_closure(param_set::APS, sde_model::sde_struct{SDEClosureType}, term::AbstractEntDet)
    θ = sde_θ(param_set, term)
    σ = sde_σ(param_set, term)
    u0 = sde_model.u0
    dt = sde_model.dt
    u = sde(θ, σ, u0, dt)
    sde_model.u0 = u
    return u
end

"""
    Solve the Ornstein-Uhlenbeck process `du = f(u,p,t)⋅dt + g(u,p,t)⋅dW` numerically

This formulation solves the Cox-Ingersoll-Ross model, a variant of the Vasicek model
which ensures that `u` remains positive so long as `θ, μ > 0`,

    `du = θ⋅(μ - u)⋅dt + σ⋅√u⋅dW`.

ref: https://en.wikipedia.org/wiki/Cox–Ingersoll–Ross_model
"""
function sde(θ::FT, σ::FT, u0::FT, dt::FT) where {FT <: Real}
    μ = FT(1)                       # μ :: long-term mean (fixed to 1)
    u_pos(u) = max(u, FT(0))        # ensure u remains non-negative
    f(u, p, t) = θ * (μ - u)        # θ :: speed of reversion
    g(u, p, t) = σ * √(u_pos(u))    # σ :: standard deviation
    tspan = (0.0, dt)
    prob = SDE.SDEProblem(f, g, u0, tspan)
    sol = SDE.solve(prob, SDE.SOSRI())
    return u_pos(sol[end])
end
