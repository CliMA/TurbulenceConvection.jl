"""
    InternalClimaParams

This is a module is an interface to CLIMAParameters.jl.
"""
module InternalClimaParams

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractEarthParameterSet

α_a(ps::APS) = CP.Atmos.EDMF.α_a(ps)
α_b(ps::APS) = CP.Atmos.EDMF.α_b(ps)
α_d(ps::APS) = CP.Atmos.EDMF.α_d(ps)

β(ps::APS) = CP.Atmos.EDMF.β(ps)

c_d(ps::APS) = CP.Atmos.EDMF.c_d(ps)
c_m(ps::APS) = CP.Atmos.EDMF.c_m(ps)
cv_l(ps::APS) = CP.Planet.cv_l(ps)

c_δ(ps::APS) = CP.Atmos.EDMF.c_δ(ps)
c_ε(ps::APS) = CP.Atmos.EDMF.c_ε(ps)
c_λ(ps::APS) = CP.Atmos.EDMF.c_λ(ps)
c_γ(ps::APS) = CP.Atmos.EDMF.c_γ(ps)

grav(ps::APS) = CP.Planet.grav(ps)

H_up_min(ps::APS) = CP.Atmos.EDMF.H_up_min(ps)

χ(ps::APS) = CP.Atmos.EDMF.χ(ps)

κ_star²(ps::APS) = CP.Atmos.EDMF.κ_star²(ps)

LH_v0(ps::APS) = CP.Planet.LH_v0(ps)
LH_s0(ps::APS) = CP.Planet.LH_s0(ps)
LH_f0(ps::APS) = CP.Planet.LH_f0(ps)

microph_scaling(ps::APS) = CP.Atmos.Microphysics.microph_scaling(ps)
microph_scaling_dep_sub(ps::APS) = CP.Atmos.Microphysics.microph_scaling_dep_sub(ps)
microph_scaling_melt(ps::APS) = CP.Atmos.Microphysics.microph_scaling_melt(ps)
molmass_ratio(ps::APS) = CP.Planet.molmass_ratio(ps)

μ_0(ps::APS) = CP.Atmos.EDMF.μ_0(ps)

Pr_n(ps::APS) = CP.Atmos.EDMF.Pr_n(ps)

Ri_c(ps::APS) = CP.Atmos.EDMF.Ri_c(ps)
R_d(ps::APS) = CP.Planet.R_d(ps)
R_v(ps::APS) = CP.Planet.R_v(ps)

smin_ub(ps::APS) = CP.Atmos.EDMF.smin_ub(ps)
smin_rm(ps::APS) = CP.Atmos.EDMF.smin_rm(ps)

T_freeze(ps::APS) = CP.Planet.T_freeze(ps)

von_karman_const(ps::APS) = CP.SubgridScale.von_karman_const(ps)

w_min(ps::APS) = CP.Atmos.EDMF.w_min(ps)
ω_pr(ps::APS) = CP.Atmos.EDMF.ω_pr(ps)

end
