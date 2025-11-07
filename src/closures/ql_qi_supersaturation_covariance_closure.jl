

"""
    get_qs_from_saturation_excesses(
        S::AbstractMatrix{FT},
        weights::AbstractVector{FT},
        q_target::FT;
        q′S′::FT = one(FT),
        q′q′::Union{Nothing, FT} = nothing,
        ε::Union{Nothing, AbstractMatrix{FT}} = nothing
    ) where {FT}

Generates a log-normally distributed field `q` that is correlated with a
normally distributed input field `S`.

This function constructs `q` to meet a target mean (`q_target`) and a target
covariance with `S` (`q′S′`).

# Arguments
- `S::AbstractMatrix{FT}`: The input field (supersaturation), assumed to be normally distributed.
- `weights::AbstractVector{FT}`: Quadrature weights for calculating moments.
- `q_target::FT`: The desired weighted mean of the output log-normal field `q`.
- `q′S′::FT`: The desired weighted covariance between `q` and `S`, i.e., `⟨(q-μ_q)(S-μ_S)⟩`.
- `q′q′::Union{Nothing, FT}`: (Optional) The variance of `q`. If `nothing`, it's estimated
  using a physically consistent proxy derived from the Cauchy-Schwarz inequality.
- `ε::Union{Nothing, AbstractMatrix{FT}}`: (Optional) An uncorrelated random field to source
  the residual (uncorrelated) variance in `q`.

# Returns
- `q_ens::Matrix{FT}`: The generated log-normal field `q`.
"""
function get_qs_from_saturation_excesses(
    S::AbstractMatrix{FT},
    weights::AbstractVector{FT},
    q_target::FT;
    q′S′::FT = one(FT),
    q′q′::Union{Nothing, FT} = nothing,
    ε::Union{Nothing, AbstractMatrix{FT}} = nothing
) where {FT}

    if iszero(q_target) 
        return zeros(FT, size(S))
    end

    if all(iszero, S)
        # If S has no variance, q must be constant to have the correct mean
        return fill(q_target, size(S))
    end

    # --- 1. Calculate Statistics for Gaussian S ---
    W = weights .* weights'
    Wtot = sum(W)
    eps_FT = eps(FT)
    eps_q = iszero(q_target) ? eps_FT : abs(q_target) * eps_FT

    S̄ = sum(W .* S) / Wtot      # Mean of S (μ_S)
    δS = S .- S̄
    var_S = sum(W .* δS .^ 2) / Wtot # Variance of S (σ_S²)

    if var_S < eps_FT
        # If S has no variance, q must be constant to have the correct mean
        return fill(q_target, size(S))
    end

    # --- 2. Calculate Parameters for Log-Normal q ---
    # Determine the variance of q (q′q′)
    q_var = if isnothing(q′q′)
        # Stable proxy from Cauchy-Schwarz: q′q′ ≥ (q′S′)² / var_S
        # (q′S′)^2 / var_S
        (q′S′^2) / max(var_S, eps_q)
    else
        q′q′
    end

    # Convert linear-space stats (q_target, q_var) to log-space stats (ν_q, s_q²)
    # for the underlying Gaussian variable Y = log(q).
    s_q_sq = log(q_var / max(q_target^2, eps_FT) + 1.0)
    ν_q = log(q_target) - s_q_sq / 2.0 # This is `c` in the derivation


    # --- THIS IS THE FIX ---
    s_q = sqrt(s_q_sq) # <--- ADD THIS LINE
    # Enforce the physical limit on the covariance to prevent numerical overflow.
    # The requested covariance q′S′ cannot require more log-space variance than is available (s_q_sq).
    max_q_S_cov = s_q * q_target * sqrt(var_S)
    # q_S_cov_limited = max(min(q′S′, max_q_S_cov), -max_q_S_cov)
    q_S_cov_limited = safe_clamp(q′S′, -max_q_S_cov, max_q_S_cov)
    q′S′ = q_S_cov_limited
    # --- END FIX ---

    # --- 3. Link the Two Distributions via Covariance ---
    # This is the core of the hybrid model, calculating the regression coefficient `b`
    # that maps the variation in S to variation in log(q).
    # b = Cov(log(q), S) / Var(S) ≈ (Cov(q, S)/E[q]) / Var(S)
    b = (q′S′ / q_target) / max(var_S, eps_FT)

    # --- 4. Partition Variance and Handle Residuals ---
    log_q_var_explained_by_S = b^2 * var_S
    # The residual variance must be non-negative.
    log_q_var_residual = max(0.0, s_q_sq - log_q_var_explained_by_S)

    residual_term = zeros(FT, size(S))
    if !isnothing(ε) && log_q_var_residual > eps_FT
        # Standardize the optional input ε to create the uncorrelated residual field
        με = sum(W .* ε) / Wtot
        δ_ε = ε .- με
        var_ε = sum(W .* δ_ε .^ 2) / Wtot
        if var_ε > eps_FT
            # Scale the residual to have the correct variance
            residual_term = δ_ε ./ sqrt(var_ε) .* sqrt(log_q_var_residual)
        end
    end

    # --- 5. Generate the Correlated Log-Normal Field ---
    # Construct log(q) as a linear function of S, making it Gaussian.
    log_q_ens = @. ν_q + b * δS + residual_term
    # Exponentiate to get the final log-normal field q.
    q_ens = exp.(log_q_ens)

    return q_ens
end



# using Statistics   # for mean
# using LinearAlgebra  # for dot product if needed
# using Random
# Random.seed!(1234)  # reproducible results

# # Helper to compute weighted covariance
# function weighted_cov(W, S, Q)
#     S̄ = sum(W .* S) / sum(W)
#     Q̄ = sum(W .* Q) / sum(W)
#     return sum(W .* (S .- S̄) .* (Q .- Q̄)) / sum(W)
# end

# # Run a single test
# function run_test_old(S::AbstractMatrix{FT}, weights::AbstractVector{FT}; B::FT=FT(1.0), tol::FT=FT(1e-12)) where FT
#     W = weights * weights'
#     q_target = rand() * 10
#     q_ens = get_qs_from_saturation_excesses_old(S, weights, q_target; B=B)

#     weighted_mean = sum(W .* q_ens) / sum(W)
#     @assert abs(weighted_mean - q_target) < tol "Weighted mean does not match q_target!"

#     @assert all(q_ens .>= 0) "q_ens contains negative values!"

#     if B > 0
#         cov_s_q = weighted_cov(W, S, q_ens)
#         @assert cov_s_q >= 0 "Weighted covariance with s should be positive for B>0!"
#     end
#     @info "Test passed" q_target=q_target weighted_mean=weighted_mean all_non_negative=all(q_ens .>= 0) weighted_covariance=cov_s_q
#     @info "S" S
#     @info "q_ens" q_ens
#     println("")
# end

# function run_test(S::AbstractMatrix{FT}, weights::AbstractVector{FT}; q′S′::FT=FT(1.0), tol::FT=FT(1e-12)) where FT
#     W = weights * weights'
#     q_target = rand() * 10
#     # q_ens_old = get_qs_from_saturation_excesses_old(s, weights, q_target; B=B)
#     q_ens = get_qs_from_saturation_excesses(S, weights, q_target, q′S′=q′S′)

#     weighted_mean = sum(W .* q_ens) / sum(W)
#     # @assert abs(weighted_mean - q_target) < tol "Weighted mean does not match q_target! Got $weighted_mean, expected $q_target"
#     @assert abs(1-abs(weighted_mean/ q_target)) < 0.2 "Weighted mean does not match q_target! Got $weighted_mean, expected $q_target"

#     @assert all(q_ens .>= 0) "q_ens contains negative values!"

#     cov_s_q = weighted_cov(W, S, q_ens)
#     # @assert abs(cov_s_q - q′S′) < tol "Weighted covariance with s does not match target q′S′! Got $cov_s_q, expected $q′S′"
#     @assert abs(1-abs(cov_s_q/q′S′)) < 0.2 "Weighted covariance with s does not match target q′S′! Got $cov_s_q, expected $q′S′"
#     @info "Test passed" q_target=q_target weighted_mean=weighted_mean all_non_negative=all(q_ens .>= 0) q′S′=q′S′ weighted_covariance=cov_s_q
#     @info "S" S
#     @info "q_ens" q_ens
#     println("")
# end

# runtests = true
# if runtests; println(" ")
#     FT = Float64
#     # Edge cases
#     # weights = FT[1.0, 1.5, 0.5, 2.0]
#     weights = FT[1.0, 1.0, 1.0, 1.0, 3.4, 1.1, 0.9, 2.2, 1.3, 0.7, 1.6, 1.4, 0.8, 2.5, 1.2, 0.6]
#     N = length(weights)

#     @info "--- All positive s ---"
#     s = rand(N,N) .* 2.0
#     # run_test_old(s, weights; B=0.8)
#     run_test(s, weights; q′S′=0.5)

#     @info "--- All negative s ---"
#     s = -rand(N,N) .* 2.0
#     # run_test_old(s, weights; B=0.8)
#     run_test(s, weights; q′S′=0.5)

#     @info "--- Mixed positive/negative s ---"
#     s = rand(N,N) .* 2.0 .- 1.0
#     # run_test_old(s, weights; B=0.8)
#     run_test(s, weights; q′S′=0.5)

#     # @info "--- All zero s ---"
#     # s = zeros(N,N)
#     # # run_test_old(s, weights; B=0.8)
#     # run_test(s, weights; q′S′=0.5)
# end
