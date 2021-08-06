import Distributions

struct Entrainment end
struct Detrainment end

function stochastic_closure(param_set::APS, term::EntDet) where {EntDet}
    closure = ICP.stoch_closure(param_set)
    if closure == "none"
        # By default, unity scaling (i.e. entr/detr are unaffected)
        return 1.0
    elseif closure == "lognormal"
        return lognormal_closure(param_set, term)
    else
        error("Unknown stochastic closure $closure.")
    end
end

#
## Lognormal closure functions
#
stoch_lognormal_var(param_set, ::Entrainment) = ICP.stoch_ε_lognormal_var(param_set)
stoch_lognormal_var(param_set, ::Detrainment) = ICP.stoch_δ_lognormal_var(param_set)

function lognormal_closure(param_set::APS, term::EntDet) where {EntDet}
    lognormal_var = stoch_lognormal_var(param_set, term)
    return lognormal_sampler(1.0, lognormal_var)
end

function lognormal_sampler(m, var)
    μ = log(m^2 / √(m^2 + var))
    σ = √(log(1 + var / m^2))
    return rand(Distributions.LogNormal(μ, σ))
end

#
## Stochastic differential equations
#
# function sde(a, b, u0, dt)
#     f(u,p,t) = a*u
#     g(u,p,t) = b*u
#     sub_dt = dt // 8
#     tspan = (0.0,dt)
#     prob = SDEProblem(f,g,u0,tspan)
#     sol = solve(prob,EM(),dt=sub_dt)
# end
