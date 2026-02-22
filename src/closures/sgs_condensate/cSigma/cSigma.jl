import StaticArrays
import LinearAlgebra

#=
--------------------------------------------------------------------------------
METHODOLOGY: The "cSigma" Parameterization
Reference: Larson et al. (2011), "Parameterizing correlations between hydrometeor 
species in mixed-phase Arctic clouds", J. Geophys. Res.

THE PROBLEM:
Microphysical processes (like accretion) depend heavily on the subgrid correlation 
between species. Arbitrarily defining these correlations often leads to "impossible" 
matrices that are not positive semidefinite.

THE SOLUTION:
1. Spherical Parameterization: We predict the Cholesky factor (L) of the correlation 
   matrix (Sigma = L * L'), guaranteeing positive semidefiniteness.
   
2. "cSigma" Closure: We estimate the angles of L using:
   a) Known correlations with vertical velocity (w).
   b) A closure based on subgrid variability (std/mean).

   c_ij = c_1i * c_1j + f_ij * s_1i * s_1j   [Eq. 15]
   f_ij = alpha * S_i * S_j * sgn(c_1i * c_1j)  [Eq. 16]

EXTENSION (Prognostic Covariances):
If specific covariances (e.g., H-QT) are known/prognosed, they can be passed in 
to override the cSigma prediction for those specific elements.
--------------------------------------------------------------------------------
=#

"""
    hydrometeor_covariances(mean_ql, ..., mean_h, tke; kwargs...)

Computes the full 7x7 **Covariance Matrix** for vertical velocity and hydrometeors.
This is the primary entry point for physical models.

# Arguments
- `mean_...`: Prognostic mean values (e.g., kg/kg). Used to normalize fluctuations.
- `tke`: Turbulent Kinetic Energy (m²/s²). Used to estimate `var_w` if not provided.

# Keyword Arguments
- `var_...`: Variance of each species (e.g., kg²/kg²). 
  **Fallback:** If `NaN`, assumes `sigma = |mean|` (i.e., coefficient of variation = 1).
- `flux_...`: Turbulent flux <w'x'>. Constrains the first row of the correlation matrix.
- `cov_qt_h`: **Prognostic Covariance** between total water (qt) and enthalpy (h).
   If provided, this overrides the cSigma prediction for `Cov[6,7]`.
- `alpha`: **Decorrelation Parameter** (default 0.19).

# Returns
- `SMatrix{7,7,FT}`: The covariance matrix. 
  Indices: 1=w, 2=ql, 3=qr, 4=qi, 5=qs, 6=qt, 7=h
"""
function hydrometeor_covariances(
    mean_ql::FT, mean_qr::FT, mean_qi::FT, mean_qs::FT, mean_qt::FT, mean_h::FT, tke::FT;
    
    var_w::FT = FT(NaN),

    # Variances
    var_ql::FT=FT(NaN), var_qr::FT=FT(NaN), var_qi::FT=FT(NaN), 
    var_qs::FT=FT(NaN), var_qt::FT=FT(NaN), var_h::FT=FT(NaN),

    # Fluxes <w'x'>
    flux_ql::FT=FT(NaN), flux_qr::FT=FT(NaN), flux_qi::FT=FT(NaN), 
    flux_qs::FT=FT(NaN), flux_qt::FT=FT(NaN), flux_h::FT=FT(NaN),

    # Prognostic Covariances (Optional overrides)
    cov_qt_h::FT=FT(NaN),

    alpha::FT = FT(0.19)
) where {FT <: AbstractFloat}

    # 1. Estimate vertical velocity variance
    w_variance_val = if isnan(var_w) 
        if !isnan(tke)
            (FT(2)/FT(3) * tke) 
        else
            error("Cannot compute w variance: both var_w and tke are NaN.")
        end
    else 
        var_w
    end
    sigma_w = sqrt(w_variance_val)

    # 2. Compute Sigmas (with fallback)
    _sigma(v, m) = isnan(v) ? abs(m) : sqrt(v)

    sigmas = StaticArrays.SVector{6, FT}(
        _sigma(var_ql, mean_ql), _sigma(var_qr, mean_qr),
        _sigma(var_qi, mean_qi), _sigma(var_qs, mean_qs),
        _sigma(var_qt, mean_qt), _sigma(var_h, mean_h)
    )

    # 3. Compute Normalized Inputs (S_x = sigma / mean)
    _norm(s, m) = abs(m) <= sqrt(eps(FT)) ? one(FT) : s / m

    std_over_mean = StaticArrays.SVector{6, FT}(
        _norm(sigmas[1], mean_ql), _norm(sigmas[2], mean_qr),
        _norm(sigmas[3], mean_qi), _norm(sigmas[4], mean_qs),
        _norm(sigmas[5], mean_qt), _norm(sigmas[6], mean_h)
    )

    # 4. Compute Correlations with w
    # rho = Flux / (sigma_w * sigma_x)
    _rho(flux, s, sw) = (isnan(flux) || s <= sqrt(eps(FT))) ? zero(FT) : flux / (sw * s)

    rhos_w = StaticArrays.SVector{6, FT}(
        _rho(flux_ql, sigmas[1], sigma_w), _rho(flux_qr, sigmas[2], sigma_w),
        _rho(flux_qi, sigmas[3], sigma_w), _rho(flux_qs, sigmas[4], sigma_w),
        _rho(flux_qt, sigmas[5], sigma_w), _rho(flux_h, sigmas[6], sigma_w)
    )

    # 5. Core cSigma: Predict Correlation Matrix
    CorrMatrix = _csigma_core(rhos_w, std_over_mean, alpha)

    # 6. Scale to Covariance
    # Broadcasting (Outer Product) is zero-allocation for SMatrix
    full_sigmas = StaticArrays.SVector{7, FT}(sigma_w, sigmas[1], sigmas[2], sigmas[3], sigmas[4], sigmas[5], sigmas[6])
    Cov = CorrMatrix .* (full_sigmas * full_sigmas')

    # 7. Apply Prognostic Covariance Override
    # If a specific covariance (qt-h) is prognosed, we overwrite the cSigma prediction.
    if !isnan(cov_qt_h)
        # We must convert to MMatrix briefly to mutate, then freeze back to SMatrix.
        # This branch is only taken if override is present.
        Cov_M = StaticArrays.MMatrix(Cov)
        Cov_M[6, 7] = cov_qt_h
        Cov_M[7, 6] = cov_qt_h
        return StaticArrays.SMatrix(Cov_M)
    end

    return Cov
end

"""
    hydrometeor_correlations_from_covariance(Cov)

Extract unique correlations from a 7x7 covariance matrix.
"""
function hydrometeor_correlations_from_covariance(Cov::StaticArrays.SMatrix{7,7,FT}) where {FT}
    @inline function _corr(Cov, i, j)
        denom = sqrt(Cov[i,i] * Cov[j,j])
        # Use relative threshold based on the covariance scale, not absolute
        rel_tol = max(abs(Cov[i,j]), abs(Cov[i,i]), abs(Cov[j,j])) * sqrt(eps(FT))
        return denom > rel_tol ? Cov[i,j] / denom : zero(FT)
    end
    
    return (
        corr_w_ql = _corr(Cov, 1, 2), corr_w_qr = _corr(Cov, 1, 3),
        corr_w_qi = _corr(Cov, 1, 4), corr_w_qs = _corr(Cov, 1, 5),
        corr_w_qt = _corr(Cov, 1, 6), corr_w_h = _corr(Cov, 1, 7),
        corr_ql_qr = _corr(Cov, 2, 3), corr_ql_qi = _corr(Cov, 2, 4),
        corr_ql_qs = _corr(Cov, 2, 5), corr_ql_qt = _corr(Cov, 2, 6),
        corr_ql_h = _corr(Cov, 2, 7), corr_qr_qi = _corr(Cov, 3, 4),
        corr_qr_qs = _corr(Cov, 3, 5), corr_qr_qt = _corr(Cov, 3, 6),
        corr_qr_h = _corr(Cov, 3, 7), corr_qi_qs = _corr(Cov, 4, 5),
        corr_qi_qt = _corr(Cov, 4, 6), corr_qi_h = _corr(Cov, 4, 7),
        corr_qs_qt = _corr(Cov, 5, 6), corr_qs_h = _corr(Cov, 5, 7),
        corr_qt_h = _corr(Cov, 6, 7),
    )
end

"""
    hydrometeor_correlations(::Type{FT}; kwargs...)

Low-level wrapper for direct access to the correlation matrix.
"""
function hydrometeor_correlations(::Type{FT};
    rho_w_ql = FT(0), rho_w_qr = FT(0), rho_w_qi = FT(0),
    rho_w_qs = FT(0), rho_w_qt = FT(0), rho_w_h = FT(0),
    std_over_mean_ql = FT(1), std_over_mean_qr = FT(1), std_over_mean_qi = FT(1),
    std_over_mean_qs = FT(1), std_over_mean_qt = FT(1), std_over_mean_h = FT(1),
    alpha = FT(0.19)
) where {FT <: AbstractFloat}

    rhos = StaticArrays.SVector{6, FT}(
        rho_w_ql, rho_w_qr, rho_w_qi, rho_w_qs, rho_w_qt, rho_w_h
    )
    stds = StaticArrays.SVector{6, FT}(
        std_over_mean_ql, std_over_mean_qr, std_over_mean_qi, 
        std_over_mean_qs, std_over_mean_qt, std_over_mean_h
    )
    
    return _csigma_core(rhos, stds, alpha)
end

"""
    _csigma_core(rhos_w, std_over_mean, alpha)

Internal kernel implementing the geometric closure.
"""
@inline function _csigma_core(
    rhos_w::StaticArrays.SVector{N, T}, 
    std_over_mean::StaticArrays.SVector{N, T}, 
    alpha::T
) where {N, T}
    
    Dim = N + 1
    # Use MMatrix for manual iterative construction (fastest compilation)
    c_param = StaticArrays.MMatrix{Dim, Dim, T}(undef)
    s_param = StaticArrays.MMatrix{Dim, Dim, T}(undef)

    # --- Step 1: Initialize Row 1 ---
    c_param[1, 1] = one(T); s_param[1, 1] = zero(T) 
    @inbounds for j in 1:N
        val = clamp(rhos_w[j], -one(T) + eps(T), one(T) - eps(T))
        c_param[1, j+1] = val
        s_param[1, j+1] = sqrt(one(T) - val^2)
    end

    # --- Step 2: Inner Parameters ---
    rows = axes(c_param, 1)
    @inbounds for i in 2:last(rows)
        for j in (i+1):last(rows)
            f_ij = clamp(alpha * std_over_mean[i-1] * std_over_mean[j-1] * sign(c_param[1, i] * c_param[1, j]), T(-0.50), T(0.50))
            val = c_param[1, i] * c_param[1, j] + f_ij * s_param[1, i] * s_param[1, j]
            c_param[i, j] = clamp(val, T(-0.95), T(0.95))
            s_param[i, j] = sqrt(one(T) - c_param[i, j]^2)
        end
    end

    # --- Step 3: Construct L^T ---
    L_T = StaticArrays.MMatrix{Dim, Dim, T}(undef)
    fill!(L_T, zero(T)) 

    @inbounds for i in axes(L_T, 1)
        for j in i:last(axes(L_T, 2))
            if i == 1
                L_T[1, j] = (j == 1) ? one(T) : c_param[1, j]
            else
                term = one(T)
                for k in 1:(i-1); term *= s_param[k, j]; end
                L_T[i, j] = term * (i == j ? one(T) : c_param[i, j])
            end
        end
    end

    # --- Step 4: Compute Sigma = L * L^T ---
    L_S = StaticArrays.SMatrix(L_T)
    Sigma = L_S' * L_S

    # Enforce unit diagonal using a Generator Expression (Zero Allocation)
    result = StaticArrays.SMatrix{Dim, Dim, T}(
        (i == j ? one(T) : Sigma[i, j]) for i in 1:Dim, j in 1:Dim
    )
    
    return result
end

# # ====================================================================================
# # TESTS
# # ====================================================================================
# using Test

# @testset "CSigma Logic Tests" begin
#     # 1. Test Float64
#     Sigma = hydrometeor_correlations(Float64;
#         rho_w_ql = 0.5, rho_w_qr = -0.5,
#         std_over_mean_ql = 1.0, std_over_mean_qr = 1.0, alpha = 0.19
#     )
#     @test size(Sigma) == (7, 7)
#     @test eltype(Sigma) == Float64
#     @test Sigma[1, 1] ≈ 1.0
    
#     # 2. Test Float32 Fallback
#     mean_val = 2.0f0
#     tke = 1.5f0
#     Cov = hydrometeor_covariances(mean_val, mean_val, mean_val, mean_val, mean_val, mean_val, tke)
#     @test eltype(Cov) == Float32
#     @test isapprox(Cov[2, 2], mean_val^2, atol=1e-5)

#     # 3. Test Prognostic Covariance Override
#     # Pass a known covariance for QT-H (index 6, 7)
#     known_cov = 0.123f0
#     Cov_O = hydrometeor_covariances(
#         mean_val, mean_val, mean_val, mean_val, mean_val, mean_val, tke;
#         cov_qt_h = known_cov
#     )
#     @test Cov_O[6, 7] == known_cov
#     @test Cov_O[7, 6] == known_cov
# end

# # Global function to verify zero allocations
# function perform_allocation_check()
#     println("Running Allocation Check...")
#     # Warmup
#     hydrometeor_correlations(Float32; rho_w_ql=0.5f0)
#     # Measure
#     allocs = @allocated hydrometeor_correlations(Float32; rho_w_ql=0.5f0)
    
#     if allocs == 0
#         println("✅ Zero allocations confirmed.")
#     else
#         error("❌ Allocations detected: $allocs bytes")
#     end
# end

# perform_allocation_check()