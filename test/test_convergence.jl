include(joinpath(dirname(@__DIR__), "perf", "common.jl"))

import ClimaCore
import LinearAlgebra
const CC = ClimaCore

non_dimensionalize!(u, v) = nothing
function non_dimensionalize!(u::FV, v::FV) where {FV <: CC.Fields.Field}
    if isempty(propertynames(u))
        ∑u = sum(abs.(parent(u))) / length(parent(u))
        ∑v = sum(abs.(parent(v))) / length(parent(v))
        ∑uv = (∑u + ∑v) / 2
        if !(∑uv ≈ 0)
            u .= u ./ ∑uv
            v .= v ./ ∑uv
        end
    else
        for sym in propertynames(u)
            uvar = getproperty(u, sym)
            vvar = getproperty(v, sym)
            non_dimensionalize!(uvar, vvar)
        end
    end
    return nothing
end

function non_dimensionalize!(u::FV, v::FV) where {FV <: CC.Fields.FieldVector}
    if all(sym -> sym in propertynames(u), (:cent, :face))
        non_dimensionalize!(u.cent, v.cent)
        non_dimensionalize!(u.face, v.face)
        return nothing
    end
    for sym in propertynames(u)
        uvar = getproperty(u, sym)
        vvar = getproperty(v, sym)
        non_dimensionalize!(uvar, vvar)
    end
    return nothing
end

function compute_err_norm(u, v)
    # We need to normalize individual fields so that,
    # for example, temperature doesn't dominate the error
    # over velocity/momentum.
    non_dimensionalize!(u, v)
    LinearAlgebra.norm(parent(u) .- parent(v))
end

sim = init_sim("Bomex"; single_timestep = false)
sim.TS.t_max = sim.TS.t_max / 5
(prob, alg, kwargs) = solve_args(sim)
integrator = ODE.init(prob, alg; kwargs...);

import ODEConvergenceTester
ODEConvergenceTester.refinement_study(integrator; refinement_range = 11:16, compute_err_norm = compute_err_norm)
nothing
