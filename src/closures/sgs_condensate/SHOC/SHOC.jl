"""
SHOC (Simplified Higher-Order Closure) Turbulence Parameterization

Adapted from Bogenschutz & Krueger (2013) JAMES for TurbulenceConvection.jl
This implementation uses ClimaCore fields and TC infrastructure instead of matrices.

References:
- Bogenschutz & Krueger (2013) https://doi.org/10.1002/jame.20018  
- SHOC-MF: Chinita et al (2022) https://doi.org/10.5194/gmd-16-1909-2023
"""

import SpecialFunctions: SpecialFunctions

# ============================================================
# Physical constants matching SHOC Fortran
# ============================================================
const SHOC_VK = 0.4         # von Karman constant
const SHOC_MAXLEN = 20000.0 # Maximum mixing length [m]
const SHOC_MINLEN = 20.0    # Minimum mixing length [m]  
const SHOC_MAXTKE = 50.0    # Maximum TKE [m²/s²]
const SHOC_MINTKE = 0.0004  # Minimum TKE [m²/s²]
const SHOC_LENGTH_FAC = 0.5 # Length scale factor
const SHOC_TROPPRES = 80000.0 # Tropopause pressure [Pa]
const SHOC_THL2TUNE = 1.0
const SHOC_QW2TUNE = 1.0
const SHOC_QWTHL2TUNE = 1.0
const SHOC_W3CLIP = 1.2
const SHOC_C_DIAG_3RD_MOM = 7.0
const SHOC_LARGENEG = -9.9999999e7

# ============================================================
# Brunt-Väisälä frequency (using TC thermodynamics)
# ============================================================
"""
    compute_brunt_vaisala!(N²::CC.Fields.Field, state::State, param_set::APS)

Compute the squared Brunt-Väisälä frequency (buoyancy frequency).

Physics: The Brunt-Väisälä frequency measures atmospheric stability. It describes oscillations of a
parcel displaced vertically in a stably stratified fluid. This is essential for the SHOC mixing length
scheme, which uses stability to modulate turbulent length scales.

Formula: N² = (g/θ_v) · ∂θ_v/∂z

where:
  - g is gravitational acceleration
  - θ_v is virtual potential temperature (accounts for moisture effects)
  - ∂θ_v/∂z is computed using ClimaCore's vertical gradient operator

Physical meaning:
  - N² > 0: Stable layer (displacing parcel oscillates)
  - N² ≈ 0: Neutral layer (displacing parcel is indifferent)
  - N² < 0: Unstable layer (parcel accelerates further from equilibrium)

Arguments:
  - N²: Output field to store squared Brunt-Väisälä frequency (s⁻²)
  - state: Model state containing environment thermodynamic variables
  - param_set: Physical parameters including gravity and thermodynamics info
"""
function compute_brunt_vaisala!(N²::CC.Fields.Field, state::State, param_set::APS) where {APS}
    g = TCP.grav(param_set)
    aux_en = center_aux_environment(state)
    θ_virt = aux_en.θ_virt

    # Compute N² = (g/θ_v) * ∂θ_v/∂z using ClimaCore operators
    @. N² = (g / θ_virt) * ∇c(wvec(Ifx(θ_virt)))

    return nothing
end

# ============================================================
# L_inf: Integral length scale
# ============================================================
"""
    compute_l_inf(tke::CC.Fields.Field, grid::Grid{FT}, Δz)

Compute the integral length scale characterizing the largest turbulent eddies.

Physics: The integral length scale represents a bulk estimate of eddy size in the domain. It is used to
scale the mixing length formulation in stable to neutral conditions. The definition integrates the
vertical structure of turbulent kinetic energy, weighting it by height.

Formula: L_inf = 0.1 * ∫(√TKE · z · dz) / ∫(√TKE · dz)

where:
  - √TKE is proportional to turbulent velocity scale
  - z is height above surface
  - Integral is computed over the full domain
  - Factor 0.1 is an empirical normalization from SHOC tuning

Physical meaning: Typical eddy size in meters. Larger L_inf → larger turbulent eddies.
Typical range: 100-2000 m in boundary layer layers.

Arguments:
  - tke: Turbulent kinetic energy field (m²/s²)
  - grid: Vertical grid information for vertical indices and heights
  - Δz: Vertical grid spacing array (m)

Returns: Scalar length scale in meters
"""
function compute_l_inf(tke::CC.Fields.Field, grid::Grid{FT}, Δz) where {FT}
    numer = FT(0)
    denom = FT(0)

    @inbounds for k in real_center_indices(grid)
        sqrt_tke = sqrt(max(tke[k], SHOC_MINTKE))
        z_k = grid.zc[k].z
        dz_k = Δz[k.i]
        numer += sqrt_tke * z_k * dz_k
        denom += sqrt_tke * dz_k
    end

    return FT(0.1) * (numer / max(denom, eps(FT)))
end

# ============================================================
# Convective velocity scale
# ============================================================
"""
    compute_conv_velocity(wthv_flux::CC.Fields.Field, θ_virt::CC.Fields.Field, 
                         pblh::FT, grid::Grid{FT}, g::FT, Δz)

Compute the convective velocity scale (w*) characterizing strength of convection.

Physics: In convective boundary layers, updrafts and downdrafts transport energy and moisture. The
convective velocity scale represents a characteristic overturning speed and is deeply connected to
buoyancy-driven motion. It scales nonlinearly with the integrated buoyancy flux over the PBL depth.

Formula: w* = (∫₀^pblh 2.5 · (g/θ_v) · w'θ_v' dz)^(1/3)

where:
  - g/θ_v accounts for density variations
  - w'θ_v' is the virtual potential temperature flux (proportional to buoyancy flux)
  - Integration is over planetary boundary layer depth (pblh)
  - 2.5 factor is an empirical constant from SHOC
  - Cubic root converts integrated buoyancy to velocity scale

Physical meaning: Typical vertical velocity in convective updrafts (m/s). Larger w* → stronger
convection with deeper mixing.

Arguments:
  - wthv_flux: Vertical flux of virtual potential temperature (K·m/s)
  - θ_virt: Virtual potential temperature (K)
  - pblh: Planetary boundary layer height (m)
  - grid: Vertical grid information
  - g: Gravitational acceleration (m/s²)
  - Δz: Vertical grid spacing array (m)

Returns: (wstar, tscale) where:
  - wstar: Convective velocity scale (m/s)
  - tscale: Convective time scale = pblh/wstar (s), or default 100 s if wstar is very small
"""
function compute_conv_velocity(
    wthv_flux::CC.Fields.Field,
    θ_virt::CC.Fields.Field,
    pblh::FT,
    grid::Grid{FT},
    g::FT,
    Δz,
) where {FT}
    conv_vel3 = FT(0)

    @inbounds for k in real_center_indices(grid)
        if grid.zc[k].z < pblh
            conv_vel3 += FT(2.5) * Δz[k.i] * (g / θ_virt[k]) * wthv_flux[k]
        end
    end

    wstar = cbrt(max(conv_vel3, FT(0)))
    tscale = wstar > 0 ? pblh / wstar : FT(100)

    return wstar, tscale
end

# ============================================================
# Mixing length computation
# ============================================================
"""
    compute_shoc_mixing_length!(ℓ_mix::CC.Fields.Field, tke::CC.Fields.Field,
                                N²::CC.Fields.Field, pblh::FT, l_inf::FT,
                                tscale::FT, grid::Grid{FT}, vk::FT)

Compute SHOC mixing length: the characteristic size of turbulent eddies that transport momentum,
heat, and moisture.

Physics: The mixing length is the fundamental parameter in K-theory (eddy diffusivity). It bridges
macroscopic turbulence properties (TKE, N²) to actual transport coefficients. SHOC uses a
stability-dependent formulation: in stable layers, stratification suppresses eddies (mixing length
reduces); in unstable layers, buoyancy amplifies eddies (mixing length grows with height).

Formulation:
  Stable (N² > 0):    ℓ = L_inf · min(1, √[TKE / (N² · z²)])
  Unstable (N² ≤ 0): ℓ = κ · z  for z < pblh
                      ℓ = L_inf  for z ≥ pblh

where:
  - L_inf is the integral length scale (maximum eddies can be)
  - κ is von Kármán constant (≈ 0.4)
  - z is height above surface
  - Min function prevents ℓ from exceeding domain-scale L_inf
  - √(TKE / (N²·z²)) is the Froude number (~velocity fluctuation / stratification)

Physical meaning: In stable layers, strong N² suppresses turbulence, reducing ℓ. Near surface in unstable
layers, ℓ grows linearly with height (wall-bounded turbulence). Above PBL, ℓ is capped by L_inf.

Arguments:
  - ℓ_mix: Output mixing length field (m)
  - tke: Turbulent kinetic energy (m²/s²)
  - N²: Squared Brunt-Väisälä frequency (s⁻²)
  - pblh: Planetary boundary layer height (m)
  - l_inf: Integral length scale (m)
  - tscale: Convective time scale (s) [unused in current implementation]
  - grid: Vertical grid information
  - vk: von Kármán constant (≈ 0.4)

"""
function compute_shoc_mixing_length!(
    ℓ_mix::CC.Fields.Field,
    tke::CC.Fields.Field,
    N²::CC.Fields.Field,
    pblh::FT,
    l_inf::FT,
    tscale::FT,
    grid::Grid{FT},
    vk::FT,
) where {FT}
    @inbounds for k in real_center_indices(grid)
        z_k = grid.zc[k].z
        tke_k = max(tke[k], SHOC_MINTKE)

        # Stability-dependent length scale
        if N²[k] > eps(FT)  # Stable
            L = l_inf * min(FT(1), sqrt(tke_k / (N²[k] * z_k + eps(FT))))
        else  # Unstable
            if z_k < pblh
                L = vk * z_k
            else
                L = l_inf
            end
        end

        # Apply bounds and length_fac tuning
        ℓ_mix[k] = SHOC_LENGTH_FAC * clamp(L, SHOC_MINLEN, SHOC_MAXLEN)
    end

    return nothing
end

# ============================================================
# Shear production
# ============================================================
"""
    compute_shear_production!(shear_prod::CC.Fields.Field, u::CC.Fields.Field,
                             v::CC.Fields.Field, grid::Grid{FT}, Ck_sh::FT, Δz)

Compute shear production: the rate at which wind shear converts kinetic energy into turbulence.

Physics: When horizontal winds differ with height (wind shear), turbulent eddies extract energy from
the mean flow and convert it to turbulent kinetic energy. This is a key source term in the TKE budget,
especially in stable layers and near the surface where winds vary rapidly.

Formula: S = Ck_sh · [(∂u/∂z)² + (∂v/∂z)²]

where:
  - ∂u/∂z, ∂v/∂z are vertical gradients of horizontal wind components
  - Ck_sh is a proportionality constant (typically 0.1 for unstable, 0.16 for stable)
  - Both horizontal components contribute to shear production
  - Units: m²/s³ (power per unit mass)

Physical meaning: Stronger wind shear (rapid wind change with height) produces more turbulence.
Typically important above 100 m height and in stably-stratified regions.

Arguments:
  - shear_prod: Output shear production rate field (m²/s³)
  - u: Zonal (E-W) wind component (m/s)
  - v: Meridional (N-S) wind component (m/s)
  - grid: Vertical grid information
  - Ck_sh: Shear production coefficient (dimensionless, ~0.1)
  - Δz: Vertical grid spacing array (m) [currently unused, kept for future use]

"""
function compute_shear_production!(
    shear_prod::CC.Fields.Field,
    u::CC.Fields.Field,
    v::CC.Fields.Field,
    grid::Grid{FT},
    Ck_sh::FT,
    Δz,
) where {FT}
    # Compute velocity gradients using ClimaCore operators
    # Use inline computation to avoid allocating dudz and dvdz separately
    @. shear_prod = Ck_sh * ((∇c(wvec(Ifx(u))))^2 + (∇c(wvec(Ifx(v))))^2)

    return nothing
end

# ============================================================
# Eddy diffusivities
# ============================================================
"""
    compute_eddy_diffusivities!(Km::CC.Fields.Field, Kh::CC.Fields.Field,
                                tke::CC.Fields.Field, ℓ_mix::CC.Fields.Field,
                                Pr::FT = 1.0)

Compute eddy diffusivity coefficients for momentum (Km) and heat/moisture (Kh).

Physics: K-theory assumes turbulent fluxes are proportional to mean gradients with diffusivity
coefficients Km and Kh. These depend on turbulent kinetic energy (energy available for mixing) and
mixing length (eddy size). The ratio Kh/Km is the Prandtl number, relating momentum diffusion to
scalar diffusion.

Formula:
  K_m = ℓ_mix · √(TKE_eff)
  K_h = K_m / Pr

where:
  - ℓ_mix is the mixing length (m)
  - √(TKE_eff) = √(max(TKE, TKE_min)) is the effective turbulent velocity scale (m/s)
  - TKE_min prevents division by zero when TKE becomes very small
  - Pr is the turbulent Prandtl number (typically 1.0, meaning Kh ≈ Km)

Physical meaning: Larger mixing length or stronger turbulence → larger diffusivity → stronger mixing.
K_m ≈ 1-100 m²/s typical in boundary layers.

Arguments:
  - Km: Output momentum diffusivity field (m²/s)
  - Kh: Output heat/moisture diffusivity field (m²/s)
  - tke: Turbulent kinetic energy (m²/s²)
  - ℓ_mix: Mixing length (m)
  - Pr: Turbulent Prandtl number (dimensionless, default 1.0)

Usage: Output Km and Kh are used to compute turbulent fluxes:
  - Momentum flux: τ = -ρ·Km·(∂u/∂z)
  - Heat flux: Q = -ρ·cp·Kh·(∂T/∂z)

"""
function compute_eddy_diffusivities!(
    Km::CC.Fields.Field,
    Kh::CC.Fields.Field,
    tke::CC.Fields.Field,
    ℓ_mix::CC.Fields.Field,
    Pr::FT = 1.0,
) where {FT}
    @. Km = ℓ_mix * sqrt(max(tke, SHOC_MINTKE))
    @. Kh = Km / Pr

    return nothing
end

# ============================================================
# Third moment of vertical velocity (Canuto et al. 2001)
# ============================================================
@inline function shoc_w3_diag_third_moment(aa0::FT, aa1::FT, x0::FT, x1::FT, f5::FT) where {FT}
    denom = FT(SHOC_C_DIAG_3RD_MOM) - FT(1.2) * x0 + aa0
    return (aa1 - FT(1.2) * x1 - FT(1.5) * f5) / denom
end
"""Helper function to compute return-to-isotropy timescale from TKE and mixing length.
Δz should be from get_Δz(field) - contains vertical spacing weights (weighted jacobian parent)."""
function compute_shoc_isotropy!(
    isotropy::CC.Fields.Field,
    brunt::CC.Fields.Field,
    tke::CC.Fields.Field,
    ℓ_mix::CC.Fields.Field,
    p_c::CC.Fields.Field,
    Δz::AbstractArray,
    grid::Grid{FT},
    g::FT,
) where {FT}
    lambda_low = FT(0.001)
    lambda_high = FT(0.04)
    lambda_slope = FT(0.65)
    brunt_low = FT(0.02)
    maxiso = FT(20000)

    Ck = FT(0.1)
    Cs = FT(0.15)
    Ce = (Ck^3) / (Cs^4)
    Ce1 = Ce / FT(0.7) * FT(0.19)
    Ce2 = Ce / FT(0.7) * FT(0.51)
    Cee = Ce1 + Ce2

    # Compute domain-averaged Brunt frequency to determine lambda parameter
    brunt_int = zero(FT)
    @inbounds for k in real_center_indices(grid)
        if p_c[k] > FT(SHOC_TROPPRES)
            brunt_int += Δz[k.i] * brunt[k]
        end
    end

    lambda = lambda_low + ((brunt_int / g) - brunt_low) * lambda_slope
    lambda = clamp(lambda, lambda_low, lambda_high)

    # Compute return-to-isotropy timescale for each level
    @inbounds for k in real_center_indices(grid)
        a_diss = Cee / max(ℓ_mix[k], eps(FT)) * (tke[k]^(FT(1.5)))
        tscale = (FT(2) * tke[k]) / max(a_diss, eps(FT))
        eff_lambda = brunt[k] <= FT(0) ? FT(0) : lambda
        isotropy[k] = min(maxiso, tscale / (FT(1) + eff_lambda * brunt[k] * tscale^2))
    end

    return nothing
end

"""Helper function to compute w³ (third moment) via Canuto et al. closure.
Δz should be from get_Δz(field) - contains vertical spacing weights (weighted jacobian parent).

Use pre-allocated temporaries when possible to reduce memory allocations."""
function compute_shoc_w3!(
    w3_f::CC.Fields.Field,
    w2_c::CC.Fields.Field,
    thl_var_c::CC.Fields.Field,
    wthl_f::CC.Fields.Field,
    tke::CC.Fields.Field,
    ℓ_mix::CC.Fields.Field,
    brunt_c::CC.Fields.Field,
    θl_c::CC.Fields.Field,
    p_c::CC.Fields.Field,
    Δz::AbstractArray,
    grid::Grid{FT},
    g::FT;
    isotropy_c::CC.Fields.Field,
) where {FT}

    compute_shoc_isotropy!(isotropy_c, brunt_c, tke, ℓ_mix, p_c, Δz, grid, g)

    # Face interpolation operators
    If_wsec = CCO.InterpolateC2F(;
        bottom = CCO.SetValue(FT(2 / 3) * FT(SHOC_MINTKE)),
        top = CCO.SetValue(FT(2 / 3) * FT(SHOC_MINTKE)),
    )
    If_iso = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    If_brunt = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(SHOC_LARGENEG)), top = CCO.SetValue(FT(SHOC_LARGENEG)))
    If_thetal = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    If_thl = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))

    # Interpolate to face level
    w_sec_zi = If_wsec.(w2_c)
    isotropy_zi = If_iso.(isotropy_c)
    brunt_zi = If_brunt.(brunt_c)
    thetal_zi = If_thetal.(θl_c)
    thl_sec = If_thl.(thl_var_c)

    # Initialize w3_f to zero
    zero_field!(w3_f)

    # Main Canuto closure computation loop over center levels
    @inbounds for k in real_center_indices(grid)
        # Compute grid spacing inverse for this level
        thedz = FT(1) / Δz[k]
        thedz2 = FT(1) / (Δz[k] + Δz[k - 1])

        # Isotropic return-to-isotropy parameter and buoyancy interaction
        iso = isotropy_c[k]
        isosqrd = iso^2
        buoy_sgs2 = isosqrd * brunt_c[k]
        bet2 = g / max(θl_c[k], eps(FT))

        # Vertical differences in key variables (face values around center k)
        kf = CCO.PlusHalf(k.i)
        kpf = CCO.PlusHalf((k + 1).i)
        thl_sec_diff = thl_sec[kf] - thl_sec[kpf]
        wthl_sec_diff = wthl_f[kf] - wthl_f[kpf]
        wsec_diff = w2_c[k - 1] - w2_c[k]
        tke_diff = tke[k - 1] - tke[k]

        # Intermediate closure terms (Canuto et al. 2001)
        f0 = thedz2 * (bet2^3) * (iso^4) * wthl_f[kf] * thl_sec_diff
        f1 = thedz2 * (bet2^2) * (iso^3) * (wthl_f[kf] * wthl_sec_diff + FT(0.5) * w_sec_zi[kf] * thl_sec_diff)
        f2 =
            thedz * bet2 * isosqrd * wthl_f[kf] * wsec_diff +
            FT(2) * thedz2 * bet2 * isosqrd * w_sec_zi[kf] * wthl_sec_diff
        f3 = thedz2 * bet2 * isosqrd * w_sec_zi[kf] * wthl_sec_diff + thedz * bet2 * isosqrd * (wthl_f[kf] * tke_diff)
        f4 = thedz * iso * w_sec_zi[kf] * (wsec_diff + tke_diff)
        f5 = thedz * iso * w_sec_zi[kf] * wsec_diff

        a4 = FT(2.4) / (FT(3) * FT(SHOC_C_DIAG_3RD_MOM) + FT(5))
        a5 = FT(0.6) / (FT(SHOC_C_DIAG_3RD_MOM) * (FT(3) + FT(5) * FT(SHOC_C_DIAG_3RD_MOM)))
        omega0 = a4 / (FT(1) - a5 * buoy_sgs2)
        omega1 = omega0 / (FT(2) * FT(SHOC_C_DIAG_3RD_MOM))
        omega2 = omega1 * f3 + (FT(5) / FT(4)) * omega0 * f4

        a0 = (FT(0.52) * (FT(1) / (FT(SHOC_C_DIAG_3RD_MOM)^2))) / (FT(SHOC_C_DIAG_3RD_MOM) - FT(2))
        a1 = FT(0.87) / (FT(SHOC_C_DIAG_3RD_MOM)^2)
        a2 = FT(0.5) / FT(SHOC_C_DIAG_3RD_MOM)
        a3 = FT(0.6) / (FT(SHOC_C_DIAG_3RD_MOM) * (FT(SHOC_C_DIAG_3RD_MOM) - FT(2)))

        x0 = (a2 * buoy_sgs2 * (FT(1) - a3 * buoy_sgs2)) / (FT(1) - (a1 + a3) * buoy_sgs2)
        y0 = (FT(2) * a2 * buoy_sgs2 * x0) / (FT(1) - a3 * buoy_sgs2)
        x1 = (a0 * f0 + a1 * f1 + a2 * (FT(1) - a3 * buoy_sgs2) * f2) / (FT(1) - (a1 + a3) * buoy_sgs2)
        y1 = (FT(2) * a2 * (buoy_sgs2 * x1 + (a0 / a1) * f0 + f1)) / (FT(1) - a3 * buoy_sgs2)

        aa0 = omega0 * x0 + omega1 * y0
        aa1 = omega0 * x1 + omega1 * y1 + omega2

        w3_f[kf] = shoc_w3_diag_third_moment(aa0, aa1, x0, x1, f5)
    end

    # Apply clipping to w3 amplitudes at face levels
    @inbounds for k in real_face_indices(grid)
        theterm = max(w_sec_zi[k], FT(0))
        cond = FT(SHOC_W3CLIP) * sqrt(FT(2) * theterm^3)
        tsign = w3_f[k] < FT(0) ? FT(-1) : FT(1)
        if tsign * w3_f[k] > cond
            w3_f[k] = tsign * cond
        end
    end

    return nothing
end

"""
    compute_shoc_pdf_and_ql_qi_fluxes!(edmf, state, surf, param_set, ...)

Compute SHOC PDF-based cloud statistics and partition ql/qi fluxes.
  
  This function extracts the shared logic between compute_SHOC_fluxes! and
  apply_shoc_ql_qi_fluxes! for:
  1. Computing SHOC Analytic Double-Gaussian (ADG1) cloud statistics at each level
  2. Partitioning total condensate flux (wqls) into liquid and ice fluxes
  
  Algorithm:
  1. Compute w² = (2/3) * TKE for each level
  2. Get/diagnose second moments (thl_var, qt_var, thl_qt_cov):
     - From TC fields if use_tc_second_moments=true
     - From SHOC diagnostic formulas otherwise
  3. Get/compute third moment w³:
     - From TC W_third_m if use_tc_wthird=true
     - From Canuto et al. closure otherwise
  4. For each vertical level:
     a. Call shoc_assumed_pdf_point to get:
        - cloud_frac: Cloud fraction (0..1)
        - qc_mean: TOTAL condensate (liquid + ice) [kg/kg]
        - wqls: Vertical condensate flux [m/s * kg/kg]
        - wthv_sec: Buoyancy flux
     b. Update aux_en.cloud_fraction[k] = cloud_frac
     c. **Optionally** (if update_ql_qi_from_pdf=true AND EquilibriumMoisture): 
        Partition total condensate into liquid/ice state variables:
        - aux_en.q_liq[k] = qc_mean * liquid_fraction
        - aux_en.q_ice[k] = qc_mean * (1 - liquid_fraction)
        Default: OFF (only compute fluxes, don't update state)
     d. Store wqls for flux partitioning
  5. Partition wqls flux to ql/qi fluxes using thermodynamic liquid_fraction:
     - diffusive_flux_ql = wqls * liquid_fraction
     - diffusive_flux_qi = wqls * (1 - liquid_fraction)
  
  Phase Partitioning:
  - The SHOC PDF uses mixed-phase saturation: qs = λ*qs_liq + (1-λ)*qs_ice
  - Returns TOTAL condensate (not separated into liquid/ice within PDF)
  - Partitioning happens based on thermodynamic liquid_fraction computed from
    temperature using pow_icenuc parameter in the thermodynamics module
  - **State variables (q_liq, q_ice)**: Only updated if update_ql_qi_from_pdf=true (default: false)
    * For NonEquilibriumMoisture: Never updated (q_liq/q_ice are prognostic)
    * For EquilibriumMoisture: Only if flag is true (usually diagnosed elsewhere)
  - **Fluxes (diffusive_flux_ql, diffusive_flux_qi)**: Always computed and partitioned
  
  Inputs:
  - edmf: EDMF model structure
  - state: TC state with all auxiliary fields
  - wthl_f, wqw_f: Face-level h/qt fluxes
  - tke, ℓ_mix, N², Δz: Turbulence parameters
  - use_tc_second_moments, use_tc_wthird: Flags for using TC diagnostics
  - update_ql_qi_from_pdf: If true, update q_liq/q_ice state (default: false, only compute fluxes)
  - match_qc_to_state: If true, scale PDF qc and fluxes to match existing env ql+qi
  - qc_match_min, qc_match_max: Bounds on the scaling factor (default: 0..2)
  - wthl_sec_c, wqw_sec_c: Scratch fields for second moments
  - w2_c, thl_var_scratch, wqls_c, liq_frac_c: Scratch fields
  
  Outputs (modifies in place):
  - aux_en.cloud_fraction: Updated with SHOC PDF cloud fraction
  - aux_en.q_liq, aux_en.q_ice: Updated with partitioned condensate (only if update_ql_qi_from_pdf=true AND equilibrium mode)
  - aux_en_f.diffusive_flux_ql, diffusive_flux_qi: Updated with partitioned fluxes
  
  Note: By default (update_ql_qi_from_pdf=false), only SGS fluxes are computed, not state variables.
  """
function compute_shoc_pdf_and_ql_qi_fluxes!(
    edmf::EDMFModel,
    state::State,
    surf::SurfaceBase,
    param_set::APS,
    wthl_f::CC.Fields.Field,
    wqw_f::CC.Fields.Field,
    tke::CC.Fields.Field,
    ℓ_mix::CC.Fields.Field,
    N²::CC.Fields.Field,
    Δz::AbstractArray,
    grid::Grid{FT},
    g::FT;
    use_tc_second_moments::Bool,
    use_tc_wthird::Bool,
    update_ql_qi_from_pdf::Bool = false,
    match_qc_to_state::Bool = false,
    qc_match_min::FT = FT(0),
    qc_match_max::FT = FT(2),
    wthl_sec_c::CC.Fields.Field,
    wqw_sec_c::CC.Fields.Field,
    w2_c::CC.Fields.Field,
    thl_var_scratch::CC.Fields.Field,
    wqls_c::CC.Fields.Field,
    liq_frac_c::CC.Fields.Field,
) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)
    aux_en = center_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    # aux_shoc = aux_tc.shoc  # SHOC diagnostic storage
    aux_tc_f = face_aux_turbconv(state)
    Ic = CCO.InterpolateF2C()

    @. w2_c = max((2 / 3) * tke, SHOC_MINTKE)
    # @. aux_shoc.w2 = w2_c  # Store for debugging
    # @. aux_shoc.tke_eff = w2_c * (3 / 2)  # Store effective TKE for debugging

    if use_tc_second_moments
        thl_var_c = aux_en.Hvar
        qt_var_c = aux_en.QTvar
        thl_qt_cov_c = aux_en.HQTcov
    else
        ∂θl∂z_c = aux_tc.∂θl∂z
        ∂qt∂z_c = aux_tc.∂qt∂z
        thl_var_c = thl_var_scratch
        @. thl_var_c = SHOC_THL2TUNE * (ℓ_mix^2) * (∂θl∂z_c^2)
        qt_var_c = @. SHOC_QW2TUNE * (ℓ_mix^2) * (∂qt∂z_c^2)
        thl_qt_cov_c = @. SHOC_QWTHL2TUNE * (ℓ_mix^2) * ∂θl∂z_c * ∂qt∂z_c
    end
    # @. aux_shoc.thl_var = thl_var_c  # Store for debugging
    # @. aux_shoc.qt_var = qt_var_c  # Store for debugging
    # @. aux_shoc.thl_qt_cov = thl_qt_cov_c  # Store for debugging

    if use_tc_wthird
        w3_c = aux_gm.W_third_m
    else
        w3_f = aux_tc_f.temporary_f1
        isotropy_c = aux_tc.temporary_2
        compute_shoc_w3!(
            w3_f,
            w2_c,
            thl_var_c,
            wthl_f,
            tke,
            ℓ_mix,
            N²,
            aux_en.θ_liq_ice,
            aux_en.p,
            Δz,
            grid,
            g;
            isotropy_c = isotropy_c,
        )
        w3_c = Ic.(w3_f)
    end
    # @. aux_shoc.w3 = w3_c  # Store for debugging

    @. wthl_sec_c = Ic.(wthl_f)
    @. wqw_sec_c = Ic.(wqw_f)
    # @. aux_shoc.wthl_sec = wthl_sec_c  # Store for debugging
    # @. aux_shoc.wqw_sec = wqw_sec_c  # Store for debugging
    w_en_c = aux_tc.w_en_c

    zero_field!(wqls_c)

    @inbounds for k in real_center_indices(grid)
        thl_mean_k = aux_en.θ_liq_ice[k]
        qt_mean_k = aux_en.q_tot[k]
        w_mean_k = w_en_c[k]
        w2_k = w2_c[k]
        w3_k = w3_c[k]
        thl_var_k = thl_var_c[k]
        qt_var_k = qt_var_c[k]
        thl_qt_cov_k = thl_qt_cov_c[k]
        wthl_sec_k = wthl_sec_c[k]
        wqw_sec_k = wqw_sec_c[k]
        p_k = aux_en.p[k]

        liq_frac_mean_k = TD.liquid_fraction(thermo_params, aux_en.ts[k])
        # Note: SHOC PDF returns TOTAL CONDENSATE (liquid + ice), not just liquid
        # The variable name ql_mean_k is misleading - it's actually total condensate
        (
            cloud_frac_k,
            qc_mean_k,
            wqls_k,
            _wthv_sec_k,
            _ql_var_k,
            s1_k,
            s2_k,
            thl1_1_k,
            thl1_2_k,
            qw1_1_k,
            qw1_2_k,
            w1_1_k,
            w1_2_k,
            a_k,
            qs1_k,
            qs2_k,
        ) = shoc_assumed_pdf_point(
            thl_mean_k,
            qt_mean_k,
            w_mean_k,
            w2_k,
            w3_k,
            thl_var_k,
            qt_var_k,
            thl_qt_cov_k,
            wthl_sec_k,
            wqw_sec_k,
            p_k,
            thermo_params;
            liq_frac_mean = liq_frac_mean_k,
        )

        # # Store PDF diagnostics for debugging
        # aux_shoc.pdf_qt_mean[k] = qt_mean_k
        # aux_shoc.pdf_thl_mean[k] = thl_mean_k
        # aux_shoc.pdf_qt_var_input[k] = qt_var_k
        # aux_shoc.pdf_thl_var_input[k] = thl_var_k
        # aux_shoc.pdf_thl_qt_cov_input[k] = thl_qt_cov_k

        # # Store plume states
        # aux_shoc.pdf_thl1_1[k] = thl1_1_k
        # aux_shoc.pdf_thl1_2[k] = thl1_2_k
        # aux_shoc.pdf_qw1_1[k] = qw1_1_k
        # aux_shoc.pdf_qw1_2[k] = qw1_2_k
        # aux_shoc.pdf_w1_1[k] = w1_1_k
        # aux_shoc.pdf_w1_2[k] = w1_2_k
        # aux_shoc.pdf_area_frac[k] = a_k

        # # Store saturation values (SIGNED - can be negative for undersaturation)
        # aux_shoc.pdf_qs1[k] = qs1_k
        # aux_shoc.pdf_qs2[k] = qs2_k
        # aux_shoc.pdf_s1_signed[k] = s1_k  # NOT clipped, shows true saturation/undersaturation
        # aux_shoc.pdf_s2_signed[k] = s2_k  # NOT clipped
        # aux_shoc.pdf_qsat_mean[k] = (qs1_k + qs2_k) / FT(2)

        # Optional scaling to match existing env condensate (ql+qi) when state is not updated from PDF
        qc_mean_scaled = qc_mean_k
        wqls_scaled = wqls_k
        ql_var_scaled = _ql_var_k
        if match_qc_to_state
            qc_state_k = aux_en.q_liq[k] + aux_en.q_ice[k]
            if qc_mean_k > eps(FT)
                r_qc = clamp(qc_state_k / qc_mean_k, qc_match_min, qc_match_max)
                qc_mean_scaled = qc_mean_k * r_qc
                wqls_scaled = wqls_k * r_qc
                ql_var_scaled = _ql_var_k * r_qc^2
            end
        end

        aux_en.cloud_fraction[k] = cloud_frac_k

        # aux_shoc.cloud_fraction[k] = cloud_frac_k
        # aux_shoc.ql_mean[k] = qc_mean_scaled * liq_frac_mean_k
        # aux_shoc.qi_mean[k] = qc_mean_scaled * (FT(1) - liq_frac_mean_k)
        # aux_shoc.qc_total[k] = qc_mean_scaled
        # aux_shoc.wqls[k] = wqls_scaled
        # aux_shoc.wthv_sec[k] = _wthv_sec_k
        # aux_shoc.ql_var[k] = ql_var_scaled

        # Optionally partition total condensate into liquid and ice state variables
        # (only in equilibrium mode where q_liq/q_ice are diagnostic)
        # By default, only SGS fluxes are computed, not state variables
        if update_ql_qi_from_pdf && (edmf.moisture_model isa EquilibriumMoisture)
            aux_en.q_liq[k] = qc_mean_scaled * liq_frac_mean_k
            aux_en.q_ice[k] = qc_mean_scaled * (FT(1) - liq_frac_mean_k)
        end
        wqls_c[k] = wqls_scaled
    end

    if (edmf.moisture_model isa EquilibriumMoisture) || (edmf.moisture_model isa NonEquilibriumMoisture)
        wqls_f = Ifx.(wqls_c)
        @. liq_frac_c = TD.liquid_fraction(thermo_params, aux_en.ts)
        liq_frac_f = Ifx.(liq_frac_c)
        ice_frac_f = @. one(FT) - liq_frac_f

        @. aux_tc_f.diffusive_flux_ql = wqls_f * liq_frac_f
        @. aux_tc_f.diffusive_flux_qi = wqls_f * ice_frac_f

        kf_surf = kf_surface(grid)
        kf_toa = kf_top_of_atmos(grid)
        aux_tc_f.diffusive_flux_ql[kf_surf] = surf.ρq_liq_flux
        aux_tc_f.diffusive_flux_qi[kf_surf] = surf.ρq_ice_flux
        aux_tc_f.diffusive_flux_ql[kf_toa] = FT(0)
        aux_tc_f.diffusive_flux_qi[kf_toa] = FT(0)
    end

    return nothing
end

# ============================================================
# SHOC Assumed PDF (Analytic Double Gaussian)
# ============================================================
"""
    shoc_assumed_pdf_point(
        thl_mean, qt_mean, w_mean, w2, w3,
        thl_var, qt_var, thl_qt_cov, wthl, wqw,
        p, thermo_params; dothetal_skew = false, liq_frac_mean = 1,
    )

Compute SHOC Analytic Double-Gaussian (ADG1) cloud statistics and liquid-water flux
for a single level using the full SHOC closure.

Algorithm:
1. Computes bimodal vertical velocity distribution (w1, w2) based on w2 and w3 moments
2. Computes bimodal θl and qt distributions for each velocity mode
3. Computes saturation using mixed-phase saturation: qs = λ*qs_liq + (1-λ)*qs_ice
4. Computes saturation deficit and cloud fraction per plume using ADG1 formulation
5. Computes area-weighted cloud fraction and condensate amount

Important notes:
- **ql_mean is TOTAL CONDENSATE (liquid + ice), not just liquid!** 
  The function computes condensate using mixed-phase saturation based on liq_frac_mean.
  Caller must partition into liquid/ice using liquid_fraction if needed.
- This uses saturation mixing ratio and the ADG1 saturation-deficit formulation; it is
  not a full saturation-adjustment step.
- The closure is identical for Equilibrium and NonEquilibrium moisture models; any
  phase partitioning is handled outside this function.
- `liq_frac_mean` is a linearized liquid fraction (0..1) used to blend liquid/ice
  saturation in the plume thermodynamics.

Returns:
- cloud_frac: Area-weighted cloud fraction (0..1)
- ql_mean: Mean TOTAL condensate (liquid + ice) [kg/kg]
- wqls: Vertical flux of liquid water (should be partitioned to ql/qi by caller)
- wthv_sec: Virtual potential temperature flux (buoyancy flux)
- ql_var: Variance of liquid water
- s1: Signed saturation deficit in plume 1 (can be negative for undersaturation)
- s2: Signed saturation deficit in plume 2 (can be negative for undersaturation)
- thl1_1, thl1_2: Liquid potential temperature in plumes 1a and 1b
- qw1_1, qw1_2: Total water content in plumes 1a and 1b
- w1_1, w1_2: Vertical velocity in plumes 1a and 1b
- a: Bimodal area fraction
- qs1, qs2: Saturation mixing ratio in plumes 1a and 1b

Reference: Bogenschutz & Krueger (2013), ADG1 closure
"""
function shoc_assumed_pdf_point(
    thl_mean::FT,
    qt_mean::FT,
    w_mean::FT,
    w2::FT,
    w3::FT,
    thl_var::FT,
    qt_var::FT,
    thl_qt_cov::FT,
    wthl::FT,
    wqw::FT,
    p::FT,
    thermo_params::TD.Parameters.ThermodynamicsParameters{FT};
    dothetal_skew::Bool = false,
    liq_frac_mean::FT = one(FT),
) where {FT}

    Rd = TD.TP.R_d(thermo_params)
    Rv = TD.TP.R_v(thermo_params)
    cp = TD.TP.cp_d(thermo_params)
    Lv = TD.latent_heat_vapor(thermo_params, thl_mean)
    basepres = FT(100000)
    basetemp = FT(300)
    epsterm = Rd / Rv

    thl_tol = FT(1e-2)
    rt_tol = FT(1e-4)
    w_tol_sqd = FT(2e-2)^2

    # NOTE: No early return here! The Fortran code has no such early return.
    # The helper functions (shoc_assumed_pdf_thl_parameters, etc.) already handle
    # small variances correctly by returning mean values, not zeros.
    # An early return here was causing thl1_1, thl1_2, qw1_1, qw1_2 to be zero
    # instead of ~thl_mean (~290K) and ~qt_mean, breaking cloud formation.

    sqrtw2 = sqrt(max(w2, FT(0)))
    sqrtthl = max(thl_tol, sqrt(max(thl_var, FT(0))))
    sqrtqt = max(rt_tol, sqrt(max(qt_var, FT(0))))

    Skew_w, w1_1, w1_2, w2_1, w2_2, a = shoc_assumed_pdf_vv_parameters(w_mean, w2, w3, w_tol_sqd)

    thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2 =
        shoc_assumed_pdf_thl_parameters(wthl, sqrtw2, sqrtthl, thl_var, thl_mean, w1_1, w1_2, Skew_w, a, dothetal_skew)

    qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2 =
        shoc_assumed_pdf_qw_parameters(wqw, sqrtw2, Skew_w, sqrtqt, qt_var, w1_2, w1_1, qt_mean, a)

    w1_1 = shoc_assumed_pdf_tilde_to_real(w_mean, sqrtw2, w1_1)
    w1_2 = shoc_assumed_pdf_tilde_to_real(w_mean, sqrtw2, w1_2)

    r_qwthl_1 = shoc_assumed_pdf_inplume_correlations(
        sqrtqw2_1,
        sqrtthl2_1,
        a,
        sqrtqw2_2,
        sqrtthl2_2,
        thl_qt_cov,
        qw1_1,
        qt_mean,
        thl1_1,
        thl_mean,
        qw1_2,
        thl1_2,
    )

    Tl1_1 = shoc_assumed_pdf_compute_temperature(thl1_1, basepres, p, Rd, cp)
    Tl1_2 = shoc_assumed_pdf_compute_temperature(thl1_2, basepres, p, Rd, cp)

    qs1, beta1, qs2, beta2 = shoc_assumed_pdf_compute_qs(Tl1_1, Tl1_2, p, Rd, Rv, cp, Lv, thermo_params, liq_frac_mean)

    s1, std_s1, qn1, C1 = shoc_assumed_pdf_compute_s(
        qw1_1,
        qs1,
        beta1,
        p,
        thl2_1,
        qw2_1,
        sqrtthl2_1,
        sqrtqw2_1,
        r_qwthl_1,
        basepres,
        Rd,
        cp,
        Lv,
    )

    if (qw1_1 == qw1_2) && (thl2_1 == thl2_2) && (qs1 == qs2)
        s2, std_s2, qn2, C2 = s1, std_s1, qn1, C1
    else
        s2, std_s2, qn2, C2 = shoc_assumed_pdf_compute_s(
            qw1_2,
            qs2,
            beta2,
            p,
            thl2_2,
            qw2_2,
            sqrtthl2_2,
            sqrtqw2_2,
            r_qwthl_1,
            basepres,
            Rd,
            cp,
            Lv,
        )
    end

    ql1 = min(qn1, qw1_1)
    ql2 = min(qn2, qw1_2)

    cloud_frac = min(FT(1), a * C1 + (FT(1) - a) * C2)
    ql_mean = shoc_assumed_pdf_compute_sgs_liquid(a, ql1, ql2)
    ql_var = shoc_assumed_pdf_compute_cloud_liquid_variance(a, s1, ql1, C1, std_s1, s2, ql2, C2, std_s2, ql_mean)
    wqls = shoc_assumed_pdf_compute_liquid_water_flux(a, w1_1, w_mean, ql1, w1_2, ql2)
    wthv_sec = shoc_assumed_pdf_compute_buoyancy_flux(wthl, epsterm, wqw, p, wqls, basetemp, basepres, Rd, cp, Lv)

    return cloud_frac, ql_mean, wqls, wthv_sec, ql_var, s1, s2, thl1_1, thl1_2, qw1_1, qw1_2, w1_1, w1_2, a, qs1, qs2
end

"""Compute vertical velocity distribution parameters for ADG1 PDF."""
@inline function shoc_assumed_pdf_vv_parameters(w_first::FT, w2::FT, w3::FT, w_tol_sqd::FT) where {FT}
    Skew_w = w3 / sqrt(max(w2, eps(w2))^3)
    if w2 <= w_tol_sqd
        return zero(w2), w_first, w_first, zero(w2), zero(w2), FT(0.5)
    end

    w2_1 = FT(0.4)
    w2_2 = FT(0.4)
    a = FT(0.5) * (one(w2) - Skew_w * sqrt(one(w2) / (4 * (one(w2) - w2_1)^3 + Skew_w^2)))
    a = clamp(a, FT(0.01), FT(0.99))

    sqrtw2t = sqrt(one(w2) - w2_1)
    w1_1 = sqrt((one(w2) - a) / a) * sqrtw2t
    w1_2 = -sqrt(a / (one(w2) - a)) * sqrtw2t

    w2_1 *= w2
    w2_2 *= w2

    return Skew_w, w1_1, w1_2, w2_1, w2_2, a
end

"""Compute liquid potential temperature distribution parameters for ADG1 PDF."""
@inline function shoc_assumed_pdf_thl_parameters(
    wthl::FT,
    sqrtw2::FT,
    sqrtthl::FT,
    thl_var::FT,
    thl_first::FT,
    w1_1::FT,
    w1_2::FT,
    Skew_w::FT,
    a::FT,
    dothetal_skew::Bool,
) where {FT}
    thl_tol = FT(1e-2)
    w_thresh = zero(thl_var)
    corr = clamp(wthl / (sqrtw2 * sqrtthl), FT(-1), FT(1))

    if thl_var <= thl_tol^2 || abs(w1_2 - w1_1) <= w_thresh
        return thl_first, thl_first, zero(thl_var), zero(thl_var), zero(thl_var), zero(thl_var)
    end

    thl1_1 = (-corr) / w1_2
    thl1_2 = (-corr) / w1_1

    Skew_thl = zero(thl_var)
    if dothetal_skew
        tsign = abs(thl1_2 - thl1_1)
        if tsign > FT(0.4)
            Skew_thl = FT(1.2) * Skew_w
        elseif tsign <= FT(0.2)
            Skew_thl = zero(thl_var)
        else
            Skew_thl = (FT(1.2) * Skew_w / FT(0.2)) * (tsign - FT(0.2))
        end
    end

    thl2_1 =
        (3 * thl1_2 * (1 - a * thl1_1^2 - (1 - a) * thl1_2^2) - (Skew_thl - a * thl1_1^3 - (1 - a) * thl1_2^3)) /
        (3 * a * (thl1_2 - thl1_1))
    thl2_2 =
        (-3 * thl1_1 * (1 - a * thl1_1^2 - (1 - a) * thl1_2^2) + (Skew_thl - a * thl1_1^3 - (1 - a) * thl1_2^3)) /
        (3 * (1 - a) * (thl1_2 - thl1_1))

    thl2_1 = min(FT(100), max(zero(thl_var), thl2_1)) * thl_var
    thl2_2 = min(FT(100), max(zero(thl_var), thl2_2)) * thl_var

    thl1_1 = thl1_1 * sqrtthl + thl_first
    thl1_2 = thl1_2 * sqrtthl + thl_first

    return thl1_1, thl1_2, thl2_1, thl2_2, sqrt(thl2_1), sqrt(thl2_2)
end

"""Compute total water distribution parameters for ADG1 PDF."""
@inline function shoc_assumed_pdf_qw_parameters(
    wqw::FT,
    sqrtw2::FT,
    Skew_w::FT,
    sqrtqt::FT,
    qt_var::FT,
    w1_2::FT,
    w1_1::FT,
    qw_first::FT,
    a::FT,
) where {FT}
    rt_tol = FT(1e-4)
    w_thresh = zero(qt_var)
    corr = clamp(wqw / (sqrtw2 * sqrtqt), FT(-1), FT(1))

    if qt_var <= rt_tol^2 || abs(w1_2 - w1_1) <= w_thresh
        return qw_first, qw_first, zero(qt_var), zero(qt_var), zero(qt_var), zero(qt_var)
    end

    qw1_1 = (-corr) / w1_2
    qw1_2 = (-corr) / w1_1

    tsign = abs(qw1_2 - qw1_1)
    if tsign > FT(0.4)
        Skew_qw = FT(1.2) * Skew_w
    elseif tsign <= FT(0.2)
        Skew_qw = zero(qt_var)
    else
        Skew_qw = (FT(1.2) * Skew_w / FT(0.2)) * (tsign - FT(0.2))
    end

    qw2_1 =
        (3 * qw1_2 * (1 - a * qw1_1^2 - (1 - a) * qw1_2^2) - (Skew_qw - a * qw1_1^3 - (1 - a) * qw1_2^3)) /
        (3 * a * (qw1_2 - qw1_1))
    qw2_2 =
        (-3 * qw1_1 * (1 - a * qw1_1^2 - (1 - a) * qw1_2^2) + (Skew_qw - a * qw1_1^3 - (1 - a) * qw1_2^3)) /
        (3 * (1 - a) * (qw1_2 - qw1_1))

    qw2_1 = min(FT(100), max(zero(qt_var), qw2_1)) * qt_var
    qw2_2 = min(FT(100), max(zero(qt_var), qw2_2)) * qt_var

    qw1_1 = qw1_1 * sqrtqt + qw_first
    qw1_2 = qw1_2 * sqrtqt + qw_first

    return qw1_1, qw1_2, qw2_1, qw2_2, sqrt(qw2_1), sqrt(qw2_2)
end

"""Transform normalized coordinates to physical values."""
@inline function shoc_assumed_pdf_tilde_to_real(w_first::FT, sqrtw2::FT, w1::FT) where {FT}
    return w1 * sqrtw2 + w_first
end

"""
Compute in-plume correlation coefficient between theta_l and q_t.

Inputs:
- `sqrtqw2_*`, `sqrtthl2_*`: standard deviations of total water and liquid potential temperature
  within each plume component.
- `a`: area fraction for plume 1.
- `qwthlsec`: total covariance between q_t and theta_l.
- `qw1_*`, `qw_first`, `thl1_*`, `thl_first`: component means and overall means.
"""
@inline function shoc_assumed_pdf_inplume_correlations(
    sqrtqw2_1::FT,
    sqrtthl2_1::FT,
    a::FT,
    sqrtqw2_2::FT,
    sqrtthl2_2::FT,
    qwthlsec::FT,
    qw1_1::FT,
    qw_first::FT,
    thl1_1::FT,
    thl_first::FT,
    qw1_2::FT,
    thl1_2::FT,
) where {FT}
    testvar = a * sqrtqw2_1 * sqrtthl2_1 + (one(a) - a) * sqrtqw2_2 * sqrtthl2_2
    if iszero(testvar)
        return zero(testvar)
    end
    r =
        (
            qwthlsec - a * (qw1_1 - qw_first) * (thl1_1 - thl_first) -
            (one(a) - a) * (qw1_2 - qw_first) * (thl1_2 - thl_first)
        ) / testvar
    return clamp(r, FT(-1), FT(1))
end

"""
Compute temperature from liquid potential temperature and pressure.

Inputs:
- `thl1`: liquid potential temperature for the plume component.
- `basepres`: reference pressure for theta_l.
- `pval`: local pressure.
"""
@inline function shoc_assumed_pdf_compute_temperature(thl1, basepres, pval, Rd, cp)
    return thl1 / (basepres / pval)^(Rd / cp)
end

"""
Compute saturation mixing ratio and Clausius-Clapeyron factor for both plumes.

Inputs:
- `Tl1_1`, `Tl1_2`: plume temperatures derived from theta_l.
- `pval`: local pressure.
- `thermo_params`: thermodynamic parameters for saturation functions.
- `liq_frac_mean`: linearized liquid fraction used to blend liquid/ice saturation.
"""
@inline function shoc_assumed_pdf_compute_qs(
    Tl1_1::FT,
    Tl1_2::FT,
    pval::FT,
    Rd::FT,
    Rv::FT,
    cp::FT,
    Lv::FT,
    thermo_params,
    liq_frac_mean::FT,
) where {FT}
    λ = clamp(liq_frac_mean, zero(FT), one(FT))

    ρ1 = TD.air_density(thermo_params, Tl1_1, pval)
    qs1_liq = TD.q_vap_saturation_generic(thermo_params, Tl1_1, ρ1, TD.Liquid())
    qs1_ice = TD.q_vap_saturation_generic(thermo_params, Tl1_1, ρ1, TD.Ice())
    qs1 = λ * qs1_liq + (one(FT) - λ) * qs1_ice
    beta1_liq = (Rd / Rv) * (Lv / (Rd * Tl1_1)) * (Lv / (cp * Tl1_1))
    beta1_ice = beta1_liq
    beta1 = λ * beta1_liq + (one(FT) - λ) * beta1_ice

    if Tl1_1 == Tl1_2
        return qs1, beta1, qs1, beta1
    end

    ρ2 = TD.air_density(thermo_params, Tl1_2, pval)
    qs2_liq = TD.q_vap_saturation_generic(thermo_params, Tl1_2, ρ2, TD.Liquid())
    qs2_ice = TD.q_vap_saturation_generic(thermo_params, Tl1_2, ρ2, TD.Ice())
    qs2 = λ * qs2_liq + (one(FT) - λ) * qs2_ice
    beta2_liq = (Rd / Rv) * (Lv / (Rd * Tl1_2)) * (Lv / (cp * Tl1_2))
    beta2_ice = beta2_liq
    beta2 = λ * beta2_liq + (one(FT) - λ) * beta2_ice
    return qs1, beta1, qs2, beta2
end

"""
Compute saturation deficit, its standard deviation, condensate qn, and cloud fraction.

Inputs:
- `qw1`, `qw2`: plume-mean total water and its variance.
- `qs1`: saturation mixing ratio for the plume temperature.
- `beta`: Clausius-Clapeyron factor.
- `thl2`, `sqrtthl2`: liquid potential temperature variance and standard deviation.
- `sqrtqw2`: total water standard deviation.
- `r_qwthl`: in-plume correlation between q_t and theta_l.
- `require_cloud_for_condensate`: If true, set qn=0 when cloud fraction C=0. 
  Set to false when computing fluxes with non-equilibrium microphysics.
"""
@inline function shoc_assumed_pdf_compute_s(
    qw1::FT,
    qs1::FT,
    beta::FT,
    pval::FT,
    thl2::FT,
    qw2::FT,
    sqrtthl2::FT,
    sqrtqw2::FT,
    r_qwthl::FT,
    basepres::FT,
    Rd::FT,
    cp::FT,
    Lv::FT;
    require_cloud_for_condensate::Bool = false,
) where {FT}
    sqrt2 = sqrt(FT(2))
    sqrt2pi = sqrt(FT(2) * π)

    s = qw1 - qs1 * ((one(qw1) + beta * qw1) / (one(qw1) + beta * qs1))
    cthl = ((one(qw1) + beta * qw1) / (one(qw1) + beta * qs1)^2) * (cp / Lv) * beta * qs1 * (pval / basepres)^(Rd / cp)
    cqt = one(qw1) / (one(qw1) + beta * qs1)
    tmp_val = max(zero(qw1), cthl^2 * thl2 + cqt^2 * qw2 - 2 * cthl * sqrtthl2 * cqt * sqrtqw2 * r_qwthl)
    std_s = sqrt(tmp_val)

    if std_s > sqrt(eps(std_s)) * 100
        C = FT(0.5) * (one(qw1) + SpecialFunctions.erf(s / (sqrt2 * std_s)))
        if require_cloud_for_condensate && iszero(C)
            qn = zero(qw1)
        else
            qn = s * C + (std_s / sqrt2pi) * exp(-FT(0.5) * (s / std_s)^2)
        end
        return s, std_s, qn, C
    end

    if s > 0
        return s, std_s, s, one(qw1)
    end

    return s, std_s, zero(qw1), zero(qw1)
end

"""Compute subgrid liquid water mean from ADG1 PDF parameters."""
@inline function shoc_assumed_pdf_compute_sgs_liquid(a::FT, ql1::FT, ql2::FT) where {FT}
    return max(zero(ql1), a * ql1 + (one(ql1) - a) * ql2)
end

"""Compute variance of cloud liquid water from double-Gaussian statistics."""
@inline function shoc_assumed_pdf_compute_cloud_liquid_variance(
    a::FT,
    s1::FT,
    ql1::FT,
    C1::FT,
    std_s1::FT,
    s2::FT,
    ql2::FT,
    C2::FT,
    std_s2::FT,
    ql_mean::FT,
) where {FT}
    ql2_val = a * (s1 * ql1 + C1 * std_s1^2) + (one(ql1) - a) * (s2 * ql2 + C2 * std_s2^2) - ql_mean^2
    return max(zero(ql_mean), ql2_val)
end

@inline function shoc_assumed_pdf_compute_liquid_water_flux(
    a::FT,
    w1_1::FT,
    w_first::FT,
    ql1::FT,
    w1_2::FT,
    ql2::FT,
) where {FT}
    return a * (w1_1 - w_first) * ql1 + (one(w_first) - a) * (w1_2 - w_first) * ql2
end

"""Compute buoyancy flux w'θ_v' from sensible and latent heat fluxes and liquid water flux."""
@inline function shoc_assumed_pdf_compute_buoyancy_flux(
    wthlsec::FT,
    epsterm::FT,
    wqwsec::FT,
    pval::FT,
    wqls::FT,
    basetemp::FT,
    basepres::FT,
    Rd::FT,
    cp::FT,
    Lv::FT,
) where {FT}
    return wthlsec +
           ((one(wthlsec) - epsterm) / epsterm) * basetemp * wqwsec +
           ((Lv / cp) * (basepres / pval)^(Rd / cp) - (one(wthlsec) / epsterm) * basetemp) * wqls
end

# ============================================================
# Main SHOC flux computation
# ============================================================
"""
    compute_SHOC_fluxes!(edmf::EDMFModel, state::State, surf::SurfaceBase, param_set::APS)


Physics Overview:
  SHOC is a PDF-based subgrid turbulence closure designed for single-column or coarse 3D models.
  It diagnoses turbulent mixing from mean-state variables (no TKE equations) using:
    - N²: buoyancy frequency → stratification effects
    - TKE: available turbulent energy → mixing strength
    - Wind shear: horizontal wind gradients → mechanical production
    - Double-Gaussian PDF: finite-area assumptions → cloud properties

EDMF Structure:
  This is a Simplified Higher-Order Closure implementation WITHIN the EDMF (Eddy-Diffusivity/
  Mass-Flux) framework. TC has updrafts + environment. SHOC currently provides diffusivity for
  the ENVIRONMENT component (not updrafts). All operations use environment variables from aux_en.

Workflow (this function performs steps 1-12):
  1. Compute Brunt-Väisälä frequency N² from θ_v profile → measures stratification
  2. Get TKE from environment state (aux_en.tke)
  3. Compute integral length scale L_inf from ∫√TKE·z structure
  4. Estimate PBL height (currently hardcoded to 1000 m) [TODO: use diagnostic]
  5. Compute convective velocity w* from integrated buoyancy flux using ENVIRONMENT θ_virt
  6. Compute mixing length ℓ from stability (N²), TKE, and L_inf
  7. Compute eddy diffusivities Km, Kh from ℓ and TKE
  8. Optionally compute shear production (for future TKE budget)
  9. Interpolate diffusivities to vertical faces
  10. Use pre-computed environment gradients (∂θl∂z, ∂qt∂z from update_aux.jl)
  11. Compute turbulent fluxes as -K·(∂var/∂z) and store in auxiliary fields
  12. Run SHOC ADG1 assumed-PDF to diagnose cloud fraction, ql, and wqls

Key Features:
  - Uses ClimaCore fields (not matrices) for all state variables
  - Operates on ENVIRONMENT fields (aux_en), not grid-mean (aux_gm)
  - Uses pre-computed gradients from aux_tc (computed in update_aux.jl) to avoid recomputation
  - Fully diagnostic: all computations from instantaneous fields (no prognostic TKE equation)
  - Full ADG1 double-Gaussian PDF for cloud fraction and liquid water
  - Computes momentum (Km) and heat/moisture (Kh) diffusivities for environment

Arguments:
  - edmf: EDMF model (contains updraft/environment structure info)
  - state: Complete model state containing:
    * center_aux_environment: Environment variables (θ_virt, TKE, velocities, θ_liq_ice, q_tot)
    * center_aux_turbconv: Working arrays and pre-computed gradients (∂θl∂z, ∂qt∂z)
    * center_prog_grid_mean: Prognostic variables (ρ for grid spacing)
    * face_aux_turbconv: Face-level fluxes (output location)
  - surf: Surface model (unused currently)
  - param_set: Physical parameters (gravity, thermodynamics info)

Keyword arguments (defaults defer to TC EDMF values):
  - use_tc_n2: If true, use aux_tc.∂θv∂z to build N²; else compute N² with SHOC
  - use_tc_mixing_length: If true, keep aux_tc.mixing_length; else compute SHOC mixing length
  - use_tc_kh_km: If true, use aux_tc.KH/KM; else compute SHOC diffusivities
  - use_tc_shear: If true, use TC shear production (aux_en_2m.tke.shear); else compute SHOC shear from grid-mean u/v
  - use_tc_second_moments: If true, use aux_en.Hvar/QTvar/HQTcov; else diagnose from SHOC closure
  - use_tc_wthird: If true, use W_third_m from aux_gm; else diagnose w3 via Canuto et al. (2001)

Output:
  Fluxes stored in state.face_aux_turbconv fields:
    * diffusive_flux_h: Environment vertical heat flux (K·m/s) [downgradient]
    * diffusive_flux_qt: Environment vertical moisture flux (kg/kg·m/s) [downgradient]
    * diffusive_flux_ql: Liquid-water flux from SHOC PDF (if NonEq)
    * diffusive_flux_qi: Ice-water flux from SHOC PDF (if NonEq)
    * diffusive_flux_qr: Rain flux using standard TC diffusive formulation
    * diffusive_flux_qs: Snow flux using standard TC diffusive formulation
  These fluxes represent turbulent transport in the ENVIRONMENT component of EDMF.
  They are then integrated vertically by TC's time-stepping to update prognostic variables.

Limitations (Known Issues for Future Work):
  - PBL height hardcoded to 1000 m (should use diagnosis)
  - TKE source must exist in state (aux_en.tke) - currently assumed but not validated
  - PDF skewness set to zero (should depend on buoyancy flux)
  - Buoyancy fluxes computed but set to zero (marked TODO)
  - Surface boundary conditions enforced for h/qt; ql/qi derived from SHOC wqls
  - No TKE equation source terms (shear production, dissipation not integrated)

References:
  - Bogenschutz & Krueger (2013): https://doi.org/10.1002/jame.20018
  - SHOC-MF (2022): https://doi.org/10.5194/gmd-16-1909-2023
"""
function compute_SHOC_fluxes!(
    edmf::EDMFModel,
    state::State,
    surf::SurfaceBase,
    param_set::APS;
    use_tc_n2::Bool = true,
    use_tc_mixing_length::Bool = true,
    use_tc_kh_km::Bool = true,
    use_tc_shear::Bool = true,
    use_tc_second_moments::Bool = true,
    use_tc_wthird::Bool = false,
)
    FT = float_type(state)
    grid = Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    g = TCP.grav(param_set)

    # Access state variables
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_tc = center_aux_turbconv(state)
    # aux_shoc = aux_tc.shoc  # SHOC diagnostic storage
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    aux_tc_f = face_aux_turbconv(state)

    # Get grid spacing
    Δz = get_Δz(prog_gm.ρ) # centers

    # Get or allocate working arrays
    N² = aux_tc.temporary_1  # Reuse temporary field
    ℓ_mix = aux_tc.mixing_length

    # 1. Compute Brunt-Väisälä frequency using TC or SHOC formulation
    if use_tc_n2
        θ_virt = aux_en.θ_virt
        ∂θv∂z = aux_tc.∂θv∂z
        @. N² = (g / θ_virt) * ∂θv∂z # is this the same definition as the N2_dry = ∂b∂θv * aux_tc.∂θv∂z[k] definition in EDMF_functions.jl ?
    else
        compute_brunt_vaisala!(N², state, param_set)
    end
    # @. aux_shoc.N2 = N²  # Store N² for debugging

    # 2. Get TKE from TC environment (SHOC is a diagnostic closure, doesn't compute its own TKE)
    tke = aux_en.tke

    # 3. Compute integral length scale
    l_inf = compute_l_inf(tke, grid, Δz)
    # @. aux_shoc.l_inf = l_inf  # Store for debugging

    # 4. Estimate PBL height (simplified - use PBLH from state if available)
    pblh = FT(1000)  # Default, should use aux_gm.pblh or similar
    # @. aux_shoc.pblh = pblh  # Store for debugging

    # 5. Compute convective velocity scale (using ENVIRONMENT fields, not grid-mean)
    wthv_flux = aux_en.buoy  # Buoyancy flux in environment
    θ_virt = aux_en.θ_virt  # Environment virtual potential temperature
    wstar, tscale = compute_conv_velocity(wthv_flux, θ_virt, pblh, grid, g, Δz)
    # @. aux_shoc.wstar = wstar  # Store for debugging
    # @. aux_shoc.tscale = tscale  # Store for debugging

    # 6. Compute mixing length (or use TC precomputed value)
    if !use_tc_mixing_length
        compute_shoc_mixing_length!(ℓ_mix, tke, N², pblh, l_inf, tscale, grid, FT(SHOC_VK))
    end
    # @. aux_shoc.mixing_length = ℓ_mix  # Store for debugging

    # 7. Compute eddy diffusivities (or use TC precomputed values)
    if use_tc_kh_km
        Km_shoc = aux_tc.KM
        Kh_shoc = aux_tc.KH
    else
        Km_shoc = aux_tc.temporary_2
        Kh_shoc = aux_tc.temporary_3
        compute_eddy_diffusivities!(Km_shoc, Kh_shoc, tke, ℓ_mix, FT(1.0))
    end
    # @. aux_shoc.Km = Km_shoc  # Store for debugging
    # @. aux_shoc.Kh = Kh_shoc  # Store for debugging

    # 8. Compute shear production (optional; TC already computes aux_en_2m.tke.shear)
    if !use_tc_shear
        shear_prod = aux_tc.temporary_4
        u_gm = physical_grid_mean_u(state)
        v_gm = physical_grid_mean_v(state)
        compute_shear_production!(shear_prod, u_gm, v_gm, grid, FT(0.1), Δz)
    end

    # 9. Interpolate diffusivities to faces for flux computation
    Km_f = Ifx.(Km_shoc)
    Kh_f = Ifx.(Kh_shoc)

    # 10. Compute diffusive fluxes using ENVIRONMENT gradients
    # NOTE: Compute gradients directly on faces using GradientC2F to avoid checkerboarding
    # from double interpolation (C→F→C→F pattern). The standard approach in update_aux.jl
    # computes gradients at centers (∇c(Ifx(...))), which when interpolated to faces again
    # causes checkerboarding artifacts.

    # Compute gradients directly on faces from center values
    ∇f = CCO.GradientC2F(; bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))
    dθl_dz_f = @. toscalar(wvec(∇f(aux_en.θ_liq_ice)))
    dqt_dz_f = @. toscalar(wvec(∇f(aux_en.q_tot)))

    # Diffusive fluxes (downgradient) on faces using environment diffusivity
    # NOTE: These represent the ENVIRONMENT turbulent diffusive contribution in the EDMF framework
    @. aux_tc_f.diffusive_flux_h = -Kh_f * dθl_dz_f
    @. aux_tc_f.diffusive_flux_qt = -Kh_f * dqt_dz_f

    # Enforce surface/top boundary conditions for h and qt fluxes (match compute_diffusive_fluxes)
    kf_surf = kf_surface(grid)
    kf_toa = kf_top_of_atmos(grid)
    aux_tc_f.diffusive_flux_h[kf_surf] = surf.ρθ_liq_ice_flux
    aux_tc_f.diffusive_flux_qt[kf_surf] = surf.ρq_tot_flux
    aux_tc_f.diffusive_flux_h[kf_toa] = FT(0)
    aux_tc_f.diffusive_flux_qt[kf_toa] = FT(0)

    compute_shoc_pdf_and_ql_qi_fluxes!(
        edmf,
        state,
        surf,
        param_set,
        aux_tc_f.diffusive_flux_h,
        aux_tc_f.diffusive_flux_qt,
        tke,
        ℓ_mix,
        N²,
        Δz,
        grid,
        g;
        use_tc_second_moments = use_tc_second_moments,
        use_tc_wthird = use_tc_wthird,
        wthl_sec_c = aux_tc.temporary_1,
        wqw_sec_c = aux_tc.temporary_5,
        w2_c = aux_tc.temporary_3,
        thl_var_scratch = aux_tc.temporary_4,
        wqls_c = aux_tc.temporary_6,
        liq_frac_c = aux_tc.temporary_4,
    )

    # 13. Precipitation diffusive fluxes (fallback to standard TC formulation)
    ρ_f = aux_gm_f.ρ
    a_en = aux_en.area
    aeKQ = aux_tc.temporary_1
    @. aeKQ = a_en * aux_tc.KQ
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_surf = kf_surface(grid)
    IfKQ = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKQ[kc_surf]), top = CCO.SetValue(aeKQ[kc_toa]))
    @. aux_tc_f.ρ_ae_KQ = IfKQ(aeKQ) * ρ_f

    aeKQq_rai_bc = FT(0)
    aeKQq_sno_bc = FT(0)
    ∇q_rai_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_rai_bc), top = CCO.SetDivergence(FT(0)))
    ∇q_sno_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_sno_bc), top = CCO.SetDivergence(FT(0)))

    c_KQr = mixing_length_params(edmf).c_KQr
    c_KQs = mixing_length_params(edmf).c_KQs
    @. aux_tc_f.diffusive_flux_qr = -(c_KQr * aux_tc_f.ρ_ae_KQ) * ∇q_rai_en(wvec(prog_pr.q_rai))
    @. aux_tc_f.diffusive_flux_qs = -(c_KQs * aux_tc_f.ρ_ae_KQ) * ∇q_sno_en(wvec(prog_pr.q_sno))

    rain_advection_factor = mixing_length_params(edmf).c_KTKEqr
    snow_advection_factor = mixing_length_params(edmf).c_KTKEqs
    @. aux_tc_f.diffusive_flux_qr +=
        ρ_f * Ifx(
            a_en *
            min(sqrt(aux_en.tke) * rain_advection_factor, sqrt(2 * aux_en.CAPE)) *
            (aux_en.frac_supersat - FT(0.25)) *
            prog_pr.q_rai,
        )
    @. aux_tc_f.diffusive_flux_qs +=
        ρ_f * Ifx(
            a_en *
            min(sqrt(aux_en.tke) * snow_advection_factor, sqrt(2 * aux_en.CAPE)) *
            (aux_en.frac_supersat - FT(0.25)) *
            prog_pr.q_sno,
        )

    return nothing
end

"""
    apply_shoc_ql_qi_fluxes!(edmf::EDMFModel, state::State, surf::SurfaceBase, param_set::APS;
                 use_tc_second_moments::Bool = true,
                 use_tc_wthird::Bool = false)

  Update only ql/qi diffusive fluxes using the SHOC ADG1 PDF while keeping the
  existing diffusive fluxes (h, qt, qr, qs) unchanged.

  Expected workflow:
    1. Call `compute_diffusive_fluxes` to populate standard TC diffusive fluxes.
    2. Call `apply_shoc_ql_qi_fluxes!` to overwrite only ql/qi fluxes.
  """
function apply_shoc_ql_qi_fluxes!(
    edmf::EDMFModel,
    state::State,
    surf::SurfaceBase,
    param_set::APS;
    use_tc_second_moments::Bool = true,
    use_tc_wthird::Bool = false,
)
    FT = float_type(state)
    grid = Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    g = TCP.grav(param_set)

    aux_en = center_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    aux_tc_f = face_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    Δz = get_Δz(prog_gm.ρ)

    tke = aux_en.tke
    ℓ_mix = aux_tc.mixing_length
    N² = aux_tc.temporary_1
    if !use_tc_wthird
        θ_virt = aux_en.θ_virt
        ∂θv∂z = aux_tc.∂θv∂z
        @. N² = (g / θ_virt) * ∂θv∂z
    end

    compute_shoc_pdf_and_ql_qi_fluxes!(
        edmf,
        state,
        surf,
        param_set,
        aux_tc_f.diffusive_flux_h,
        aux_tc_f.diffusive_flux_qt,
        tke,
        ℓ_mix,
        N²,
        Δz,
        grid,
        g;
        use_tc_second_moments = use_tc_second_moments,
        use_tc_wthird = use_tc_wthird,
        wthl_sec_c = aux_tc.temporary_3,
        wqw_sec_c = aux_tc.temporary_5,
        w2_c = aux_tc.temporary_4,
        thl_var_scratch = aux_tc.temporary_6,
        wqls_c = aux_tc.temporary_2,
        liq_frac_c = aux_tc.temporary_4,
    )

    return nothing
end

# ============================================================
# Integration with TC update_aux
# ============================================================
"""
    apply_shoc_diffusion!(edmf::EDMFModel, state::State, surf::SurfaceBase, 
                          param_set::APS, Δt::FT)

Integration wrapper: compute SHOC fluxes for the ENVIRONMENT and store in auxiliary state.

Physics Context (EDMF):
  SHOC diagnoses turbulent diffusivity for the ENVIRONMENT in the EDMF framework. The computed
  fluxes are diagnostic (no explicit update to prognostics here). Instead, fluxes are stored in
  auxiliary fields, which TC's time-stepping routine then uses in vertical divergence calculations.
  The ENVIRONMENT receives this turbulent diffusive contribution while UPDRAFTS have their own
  mass-flux transport.

Workflow:
  1. Call compute_SHOC_fluxes! to diagnose all turbulent diffusivities and fluxes
  2. SHOC uses pre-computed environment gradients from update_aux.jl (already in aux_tc)
  3. Fluxes are automatically written to state.face_aux_turbconv.diffusive_flux_h etc.
  4. TC's time integrator reads these fluxes and computes divergence: ∂ρθ/∂t = -∂(ρ·Flux_θ)/∂z
  5. This applies specifically to the ENVIRONMENT component budget

Arguments:
  - edmf: EDMF model object (contains updraft/environment structure)
  - state: Model state (modified in-place: auxiliary fluxes are filled with environment contribution)
  - surf: Surface model
  - param_set: Physical parameters
  - Δt: Time step (s) [currently unused, but available for future implicit methods]

Usage:
Called from update_aux! in src/update_aux.jl when use_shoc_diffusive_fluxes = true.
  Must be called AFTER environment thermodynamics updates to have valid θ_virt, etc.
  No direct output; all results stored in state.face_aux_turbconv (environment fluxes).

Note on Coupling (EDMF):
  This wrapper cleanly integrates SHOC for the ENVIRONMENT without modifying TC's core time-stepping.
  Gradients used here (∂θl∂z, ∂qt∂z) are pre-computed in update_aux.jl from aux_en fields,
  so SHOC does not perform redundant gradient calculations.

"""
function apply_shoc_diffusion!(
    edmf::EDMFModel,
    state::State,
    surf::SurfaceBase,
    param_set::APS,
    Δt::FT;
    use_tc_n2::Bool = true,
    use_tc_mixing_length::Bool = true,
    use_tc_kh_km::Bool = true,
    use_tc_shear::Bool = true,
    use_tc_second_moments::Bool = true,
    use_tc_wthird::Bool = false,
) where {FT}
    compute_SHOC_fluxes!(
        edmf,
        state,
        surf,
        param_set;
        use_tc_n2 = use_tc_n2,
        use_tc_mixing_length = use_tc_mixing_length,
        use_tc_kh_km = use_tc_kh_km,
        use_tc_shear = use_tc_shear,
        use_tc_second_moments = use_tc_second_moments,
        use_tc_wthird = use_tc_wthird,
    )
    return nothing
end
