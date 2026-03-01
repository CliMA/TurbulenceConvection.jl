import Thermodynamics as TD
import StaticArrays as SA



"""
    get_quadrature_saturation_excess_matrices(
        param_set,
        thermo_params,
        quadrature_type::AbstractSGSQuadratureType,
        χ::SA.SVector{N, FT},
        qt_mean::FT,
        qt′qt′::FT,
        θl_mean::FT,
        θl′θl′::FT,
        θl′qt′::FT,
        p::FT,
        q_liq_mean::FT,
        q_ice_mean::FT,
    ) where {N, FT}

Computes the matrices of saturation excess (S_liq, S_ice) for all quadrature points 
defined by the given mean state and variances.
Returns just the required SMatrix objects to avoid heap allocations and massive memory bloat.
"""
function get_quadrature_saturation_excess_matrices(
    param_set,
    thermo_params,
    quadrature_type,
    χ::SA.SVector{N, FT},
    qt_mean::FT,
    qt′qt′::FT,
    θl_mean::FT,
    θl′θl′::FT,
    θl′qt′::FT,
    p::FT,
    q_liq_mean::FT,
    q_ice_mean::FT,
) where {N, FT}

    sqrt2 = FT(sqrt(2))
    eps_q = (qt_mean ≈ FT(0)) ? eps(FT) : (eps(FT) * qt_mean)
    eps_θ = eps(FT)

    # Initialize variables for the distribution parameters
    ν_q, ν_θ, s_q, s_θ, s2_θq, s_c, corr, σ_q, σ_θ, σ_c = (FT(0) for _ in 1:10)

    if quadrature_type isa LogNormalQuad
        # Lognormal parameters (ν, s) from mean and variance
        ν_q = log(qt_mean^2 / max(sqrt(qt_mean^2 + qt′qt′), eps_q))
        ν_θ = log(θl_mean^2 / sqrt(θl_mean^2 + θl′θl′))
        s_q = sqrt(log(qt′qt′ / max(qt_mean, eps_q)^2 + 1))
        s_θ = sqrt(log(θl′θl′ / θl_mean^2 + 1))

        # Enforce Cauchy-Schwarz inequality, numerically stable compute
        corr = θl′qt′ / max(sqrt(qt′qt′), eps_q)
        corr = max(min(corr / max(sqrt(θl′θl′), eps_θ), 1), -1)

        # Conditionals
        s2_θq = log(corr * sqrt(θl′θl′ * qt′qt′) / θl_mean / max(qt_mean, eps_q) + 1)
        s_c = sqrt(max(s_θ^2 - s2_θq^2 / max(s_q, eps_q)^2, 0))

    elseif quadrature_type isa GaussianQuad
        # limit σ_q to prevent negative qt_hat
        σ_q_lim = -qt_mean / (sqrt2 * χ[1])
        σ_q = min(sqrt(qt′qt′), σ_q_lim)
        σ_θ = sqrt(θl′θl′)

        # Enforce Cauchy-Schwarz inequality, numerically stable compute
        corr = θl′qt′ / max(σ_q, eps_q)
        corr = max(min(corr / max(σ_θ, eps_θ), 1), -1)

        # Conditionals
        σ_c = sqrt(max(1 - corr^2, 0)) * σ_θ
    end

    # Use purely allocation-free StaticArrays broadcasting
    # Broadcast a column vector (m_q) against a row vector (m_h)
    S_tup_M = broadcast(χ, χ') do χ_q, χ_h
        qt_hat = FT(0)
        h_hat = FT(0)
        
        if quadrature_type isa LogNormalQuad
            qt_hat = exp(ν_q + sqrt2 * s_q * χ_q)
            ν_c = ν_θ + s2_θq / max(s_q, eps_q)^2 * (log(qt_hat) - ν_q)
            h_hat = exp(ν_c + sqrt2 * s_c * χ_h)
        elseif quadrature_type isa GaussianQuad
            qt_hat = qt_mean + sqrt2 * σ_q * χ_q
            μ_c = θl_mean + sqrt2 * corr * σ_θ * χ_q
            h_hat = μ_c + sqrt2 * σ_c * χ_h
        end

        ts_hat = thermo_state_pθq(param_set, p, h_hat, qt_hat, q_liq_mean, q_ice_mean)
        T_hat = TD.air_temperature(thermo_params, ts_hat)
        ρ_hat = TD.air_density(thermo_params, ts_hat)

        q_vap_sat_liq_hat = TD.q_vap_saturation_generic(thermo_params, T_hat, ρ_hat, TD.Liquid())
        q_vap_sat_ice_hat = TD.q_vap_saturation_generic(thermo_params, T_hat, ρ_hat, TD.Ice())

        # Return a Tuple for the SMatrix element
        (qt_hat - q_vap_sat_liq_hat, qt_hat - q_vap_sat_ice_hat)
    end

    # Return exactly the required SMatrix objects
    return (map(t -> t[1], S_tup_M), map(t -> t[2], S_tup_M))
end

"""
    partition_condensate_into_quadrature_fractions(
        S_liq_quad,
        S_ice_quad,
        W,
        q_liq_mean,
        q_ice_mean,
        ql′ql′,
        qi′qi′,
        closure_type
    )

Returns the matrices `q_liq_quad` and `q_ice_quad` across quadrature points, partitioned 
to satisfy the mean targets and variance targets using `get_quadrature_condensate_state`.
"""
function partition_condensate_into_quadrature_fractions(
    S_liq_quad,
    S_ice_quad,
    W,
    q_liq_mean::FT,
    q_ice_mean::FT,
    ql′ql′::FT,
    qi′qi′::FT,
    closure_type::Val;
    max_boost_factor_liq::FT = FT(2.0),
    max_boost_factor_ice::FT = FT(1.1),
) where {FT}
    # Create the condensate matrix across quadrature points using our closure
    q_liq_quad = get_quadrature_condensate_state(
        S_liq_quad,
        W,
        q_liq_mean,
        ql′ql′,
        closure_type;
        max_boost_factor = max_boost_factor_liq,
    )
    q_ice_quad = get_quadrature_condensate_state(
        S_ice_quad,
        W,
        q_ice_mean,
        qi′qi′,
        closure_type;
        max_boost_factor = max_boost_factor_ice,
    )

    return q_liq_quad, q_ice_quad
end

# Define the user-facing wrapper
function get_quadrature_condensate_state(
    S::AbstractMatrix{FT},
    weights::AbstractVector{FT},
    q_target::FT,
    q′q′::FT,
    closure_type::Val; # <-- Forces compile-time dispatch
    q′S′_target::FT = zero(FT),
    max_boost_factor::FT = FT(2.0),
) where {FT}

    # Delegate to the type-stable internal methods
    return _compute_sgs_condensate(closure_type, S, weights, q_target, q′q′, q′S′_target, max_boost_factor)
end

# Compile-time route for Linear
function _compute_sgs_condensate(::Val{:linear}, S, weights, q_target, q′q′, q′S′_target, max_boost_factor)
    return get_qs_linear(S, weights, q_target, q′q′; max_boost_factor = max_boost_factor)
end

# Compile-time route for Lognormal
function _compute_sgs_condensate(::Val{:lognormal}, S, weights, q_target, q′q′, q′S′_target, max_boost_factor)
    return get_qs_lognormal(S, weights, q_target, q′q′, q′S′_target; max_boost_factor = max_boost_factor)
end

# Compile-time route for Proportional (Fallback)
function _compute_sgs_condensate(::Val{:proportional}, S, weights, q_target, q′q′, q′S′_target, max_boost_factor)
    return get_qs_proportional(S, weights, q_target; max_boost_factor = max_boost_factor)
end

function get_qs_lognormal(
    S::AbstractMatrix{FT},
    weights::AbstractVector{FT},
    q_target::FT,
    q′q′::FT,
    q′S′::FT;
    max_boost_factor::FT = FT(2.0),
) where {FT}

    if iszero(q_target) || all(iszero, S)
        return fill(q_target, size(S))
    end

    W = weights * weights'
    Wtot = sum(W)
    eps_FT = eps(FT)

    S̄ = sum(W .* S) / Wtot
    δS = S .- S̄
    var_S = sum(W .* δS .^ 2) / Wtot

    if var_S < eps_FT
        return fill(q_target, size(S))
    end

    # 1. Baseline Phase Locking
    f_sat = sum(W[S .> 0]) / Wtot
    f_sat_safe = clamp(f_sat, eps_FT, one(FT) - eps_FT)

    q_supersat = min(q_target / f_sat_safe, max_boost_factor * q_target)
    q_subsat = max(zero(FT), (q_target - (q_supersat * f_sat_safe)) / (1.0 - f_sat_safe))

    μ_q = similar(S)
    @. μ_q = ifelse(S > 0, q_supersat, q_subsat)

    # 2. Covariance Match
    cov_base = sum(W .* μ_q .* δS) / Wtot

    # If the target covariance is lower than what the strict phase-locked baseline provides,
    # and the sub-saturated regions have been tightly constrained to zero mass, it is mathematically 
    # impossible to lower the covariance. In this strict physical edge-case, we gracefully soften 
    # the phase-locking requirement to uniformly distribute mass, granting the exponential stretch 
    # full access to the sub-saturated space to accurately reach the necessary covariance target.
    if q′S′ < cov_base && q_subsat < eps_FT
        fill!(μ_q, q_target)
    end

    q_S_cov_residual = q′S′ - cov_base
    var_S_weighted = sum(W .* μ_q .* δS .^ 2) / max(Wtot, eps_FT)

    b = q_S_cov_residual / max(var_S_weighted, eps_FT)

    max_abs_δS = maximum(abs, δS)
    max_b = 4.0 * sqrt(max(q′q′, eps_FT)) / sqrt(var_S)
    max_b = max(max_b, FT(10.0)) # Fallback relax if strictly constrained

    # Prevent Float64 overflow in exp(b * δS) for highly skewed grids with tiny var_S
    if max_abs_δS > eps_FT
        max_b = min(max_b, FT(80.0) / max_abs_δS)
    end

    b = clamp(b, -max_b, max_b)

    # EXACT Newton solve for b to overcome the linear approximation breakdown
    for _ in 1:7
        # Moments of the exponentially stretched distribution
        D = sum(W .* μ_q .* exp.(b .* δS))
        N_val = sum(W .* μ_q .* exp.(b .* δS) .* δS)
        N_prime = sum(W .* μ_q .* exp.(b .* δS) .* δS .^ 2)

        D_safe = max(D, eps_FT)

        # Mean-invariant covariance: (N / D) * q_target
        f_b = q_target * (N_val / D_safe) - q′S′
        f_prime_b = q_target * (D_safe * N_prime - N_val^2) / (D_safe^2)

        if abs(f_b) < eps_FT * 100
            break
        end

        step = f_b / max(f_prime_b, eps_FT)
        # damp the Newton step to avoid wild oscillations
        b = clamp(b - 0.5 * step, -max_b, max_b)
    end

    q_cov = similar(S)
    @. q_cov = μ_q * exp(b * δS)

    q_mean_cov = sum(W .* q_cov) / Wtot
    if q_mean_cov > eps_FT
        @. q_cov = q_cov * (q_target / q_mean_cov)
    end

    # 3. EXACT Discrete Variance Injection
    q_ens = copy(q_cov)
    var_base = sum(W .* (q_cov .- q_target) .^ 2) / Wtot
    var_needed = max(zero(FT), q′q′ - var_base)

    if var_needed > eps_FT
        # We find a deterministic basis vector P_i = q_cov_i * E_i
        # where E_i is a quadratic polynomial of δS_i: E_i = δS_i^2 - a * δS_i - b
        # that is analytically orthogonalized to both 1 and δS under the q_cov weighting
        W_tilde = W .* q_cov
        W_tilde_sum = max(sum(W_tilde), eps_FT)

        μ1 = sum(W_tilde .* δS) / W_tilde_sum
        μ2 = sum(W_tilde .* δS .^ 2) / W_tilde_sum
        μ3 = sum(W_tilde .* δS .^ 3) / W_tilde_sum

        var_tilde = max(μ2 - μ1^2, eps_FT)
        a = (μ3 - μ1 * μ2) / var_tilde
        b = μ2 - a * μ1

        P = similar(S)
        @. P = q_cov * (δS^2 - a * δS - b)

        var_P = sum(W .* P .^ 2) / Wtot
        cov_q_P = sum(W .* (q_cov .- q_target) .* P) / Wtot

        A = var_P
        B = 2 * cov_q_P
        C_quad = -var_needed

        if A > eps_FT
            disc = B^2 - 4 * A * C_quad
            if disc >= 0
                c = (-B + sqrt(disc)) / (2 * A)
                @. q_ens = q_cov + c * P
            end
        end
    end

    # Enforce strict positivity
    @. q_ens = max(zero(FT), q_ens)

    # EXACT Final Mass Rescaling
    q_mean_final = sum(W .* q_ens) / Wtot

    if q_mean_final > eps_FT
        # Multiplicative scaling preserves positivity but alters variance.
        # We shift it additively where possible to preserve variance better.
        shift = q_target - q_mean_final
        # Distribute shift only where q_ens > 0 to maintain phase locking
        f_active = sum(W[q_ens .> 0]) / Wtot
        if f_active > eps_FT
            @. q_ens = q_ens + ifelse(q_ens > 0, shift / f_active, zero(FT))
        else
            @. q_ens = q_ens * (q_target / q_mean_final)
        end
        # Final safety clamp
        @. q_ens = max(zero(FT), q_ens)

        # If additive shift broke mass again because of clamping, fallback to multiplicative
        q_mean_final_2 = sum(W .* q_ens) / Wtot
        if abs(q_mean_final_2 - q_target) > eps_FT * 10 && q_mean_final_2 > eps_FT
            @. q_ens = q_ens * (q_target / q_mean_final_2)
        end
    end

    return q_ens
end


function get_qs_linear(
    S::AbstractMatrix{FT},
    weights::AbstractVector{FT},
    q_target::FT,
    q′q′::FT;
    max_boost_factor::FT = FT(2.0),
) where {FT}

    if iszero(q_target) || all(iszero, S)
        return fill(q_target, size(S))
    end

    # Guarantee 2D outer product for the grid
    W = weights * weights'
    Wtot = sum(W)
    eps_FT = eps(FT)

    f_sat = sum(W[S .> 0]) / Wtot
    f_sat_safe = clamp(f_sat, eps_FT, one(FT) - eps_FT)

    q_supersat_max = min(q_target / f_sat_safe, max_boost_factor * q_target)
    q_subsat_min = max(zero(FT), (q_target - (q_supersat_max * f_sat_safe)) / (1.0 - f_sat_safe))

    q_base = similar(S)
    if f_sat > eps_FT
        S_super_mean = sum(W[S .> 0] .* S[S .> 0]) / f_sat
        @. q_base = ifelse(S > 0, q_supersat_max * (S / max(S_super_mean, eps_FT)), q_subsat_min)
    else
        @. q_base = q_subsat_min
    end

    # Enforce exact mean on the baseline
    q_mean_base = sum(W .* q_base) / Wtot
    if q_mean_base > eps_FT
        @. q_base = q_base * (q_target / q_mean_base)
    end

    # --- EXACT Discrete Variance Injection ---
    q_ens = copy(q_base)
    var_base = sum(W .* (q_base .- q_target) .^ 2) / Wtot
    var_needed = max(zero(FT), q′q′ - var_base)

    if var_needed > eps_FT
        # Center the deterministic feature δS under the q_base weighting to perfectly preserve the mean
        S̄ = sum(W .* S) / Wtot
        δS = S .- S̄

        q_mean_sum = sum(W .* q_base)
        k = sum(W .* q_base .* δS) / max(q_mean_sum, eps_FT)

        P = similar(S)
        @. P = q_base * (δS - k)

        # Calculate exactly how P impacts the discrete variance
        var_P = sum(W .* P .^ 2) / Wtot
        cov_q_P = sum(W .* (q_base .- q_target) .* P) / Wtot

        # Solve the exact quadratic: var(q + cP) = q'q'
        A = var_P
        B = 2 * cov_q_P
        C_quad = -var_needed

        if A > eps_FT
            disc = B^2 - 4 * A * C_quad
            if disc >= 0
                c = (-B + sqrt(disc)) / (2 * A)
                @. q_ens = q_base + c * P
            end
        end
    end

    # Enforce strict positivity
    @. q_ens = max(zero(FT), q_ens)

    # EXACT Final Mass Rescaling
    q_mean_final = sum(W .* q_ens) / Wtot
    if q_mean_final > eps_FT
        shift = q_target - q_mean_final
        f_active = sum(W[q_ens .> 0]) / Wtot
        if f_active > eps_FT
            @. q_ens = q_ens + ifelse(q_ens > 0, shift / f_active, zero(FT))
        else
            @. q_ens = q_ens * (q_target / q_mean_final)
        end
        @. q_ens = max(zero(FT), q_ens)
        q_mean_final_2 = sum(W .* q_ens) / Wtot
        if abs(q_mean_final_2 - q_target) > eps_FT * 10 && q_mean_final_2 > eps_FT
            @. q_ens = q_ens * (q_target / q_mean_final_2)
        end
    end

    return q_ens
end







"""
    get_qs_proportional(
        S::AbstractMatrix{FT},
        weights::AbstractVector{FT},
        q_target::FT;
        max_boost_factor::FT = FT(2.0),
        kwargs... 
    ) where {FT}

Generates an SGS condensate field `q` that is perfectly continuous, strictly 
monotonic, and rigorously phase-locked.

### The Continuous Geometric Closure
To prevent unphysical "cliffs" or jumps at the `S = 0` phase boundary, the 
condensate is constructed as a linear combination of two continuous shape functions:
1. `H_base = S - S_min` (A continuous background tracer representing evaporation).
2. `H_wet = max(S, 0)` (Strictly phase-locked active condensation).

The coefficients `a` and `b` for `q = a*H_base + b*H_wet` are solved analytically. 
The system mathematically maximizes `b` (packing mass into the wet region) while 
strictly respecting the `max_boost_factor` limits and ensuring `a, b >= 0` to 
guarantee a monotonically increasing distribution.
"""
function get_qs_proportional(
    S::AbstractMatrix{FT},
    weights::AbstractVector{FT},
    q_target::FT;
    max_boost_factor::FT = FT(2.0),
    kwargs..., # Safely absorbs legacy variance/covariance args
) where {FT}

    if iszero(q_target) || all(iszero, S)
        return fill(q_target, size(S))
    end

    W = weights .* weights'
    Wtot = sum(W)
    eps_FT = eps(FT)

    S_min = minimum(S)
    S_max = maximum(S)

    # Fallback for perfectly uniform moisture field
    if (S_max - S_min) < eps_FT
        return fill(q_target, size(S))
    end

    # --- 1. Define the Continuous Shape Functions ---
    H_base = similar(S)
    H_wet = similar(S)
    @. H_base = max(S - S_min, zero(FT))
    @. H_wet = max(S, zero(FT))

    # --- 2. Calculate the Shape Integrals ---
    I_base = sum(W .* H_base)
    I_wet = sum(W .* H_wet)

    # How much of the base shape naturally falls into the subsaturated region?
    I_base_sub = sum(W .* H_base .* (S .<= 0))

    # --- 3. Solve for Coefficients (a, b) ---
    f_sat = sum(W[S .> 0]) / Wtot
    M_limit = f_sat * max_boost_factor * q_target * Wtot

    a = zero(FT)
    b = zero(FT)

    if I_base_sub > eps_FT
        # We want to minimize `a`, but if the target mass exceeds what we are legally 
        # allowed to pack into the wet bucket (M_limit), `a` must increase to spill 
        # the remaining mass into the subsaturated evaporating region.
        a_ideal = (q_target * Wtot - M_limit) / I_base_sub

        # Clamp ensures we never create negative slopes (a >= 0) and never exceed 
        # the total mass available (a <= q_target / I_base)
        a = clamp(a_ideal, zero(FT), (q_target * Wtot) / max(I_base, eps_FT))
    end

    if I_wet > eps_FT
        # Whatever mass is not consumed by the base shape is packed into the wet shape
        b_ideal = (q_target * Wtot - a * I_base) / I_wet
        b = max(zero(FT), b_ideal)
    end

    # --- 4. Assemble the Final Continuous Field ---
    q_ens = similar(S)
    @. q_ens = a * H_base + b * H_wet

    # Final scalar sweep to guarantee exact machine-precision mass conservation
    q_mean_raw = sum(W .* q_ens) / Wtot
    if q_mean_raw > eps_FT
        @. q_ens = q_ens * (q_target / q_mean_raw)
    end

    return q_ens
end
