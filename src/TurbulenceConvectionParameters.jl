"""
    TurbulenceConvectionParameters

This is a module is an interface to CLIMAParameters.jl.
"""
module TurbulenceConvectionParameters

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractEarthParameterSet

#####
##### TurbulenceConvection parameters
#####

cv_l(ps::APS) = CP.Planet.cv_l(ps)
grav(ps::APS) = CP.Planet.grav(ps)

LH_v0(ps::APS) = CP.Planet.LH_v0(ps)
LH_s0(ps::APS) = CP.Planet.LH_s0(ps)

microph_scaling(ps::APS) = CP.Atmos.Microphysics.microph_scaling(ps)
microph_scaling_dep_sub(ps::APS) = CP.Atmos.Microphysics.microph_scaling_dep_sub(ps)
microph_scaling_melt(ps::APS) = CP.Atmos.Microphysics.microph_scaling_melt(ps)
molmass_ratio(ps::APS) = CP.Planet.molmass_ratio(ps)

R_d(ps::APS) = CP.Planet.R_d(ps)
R_v(ps::APS) = CP.Planet.R_v(ps)

T_freeze(ps::APS) = CP.Planet.T_freeze(ps)

von_karman_const(ps::APS) = CP.SubgridScale.von_karman_const(ps)

Omega(ps::APS) = CP.Planet.Omega(ps)
T_0(ps::APS) = CP.Planet.T_0(ps)
cp_v(ps::APS) = CP.Planet.cp_v(ps)
cv_v(ps::APS) = CP.Planet.cv_v(ps)
cp_d(ps::APS) = CP.Planet.cp_d(ps)

α_a(ps::APS) = CP.Atmos.EDMF.α_a(ps)
α_b(ps::APS) = CP.Atmos.EDMF.α_b(ps)
α_d(ps::APS) = CP.Atmos.EDMF.α_d(ps)
β(ps::APS) = CP.Atmos.EDMF.β(ps)
c_d(ps::APS) = CP.Atmos.EDMF.c_d(ps)
c_m(ps::APS) = CP.Atmos.EDMF.c_m(ps)
c_λ(ps::APS) = CP.Atmos.EDMF.c_λ(ps)
χ(ps::APS) = CP.Atmos.EDMF.χ(ps)
κ_star²(ps::APS) = CP.Atmos.EDMF.κ_star²(ps)
Pr_n(ps::APS) = CP.Atmos.EDMF.Pr_n(ps)
Ri_c(ps::APS) = CP.Atmos.EDMF.Ri_c(ps)
smin_ub(ps::APS) = CP.Atmos.EDMF.smin_ub(ps)
smin_rm(ps::APS) = CP.Atmos.EDMF.smin_rm(ps)
ω_pr(ps::APS) = CP.Atmos.EDMF.ω_pr(ps)

#####
##### Experimental parameters
#####

#=
    ExperimentalClimaParams

This is a module for staging experimental
clima parameters-- that we are unsure will
be used in the final model. Once they are
no longer experimental, they should be moved
to CLIMAParameters.jl.

TODO: move non-experimental parameters to CLIMAParameters
=#

#= divergence factor for bubble case (zero otherwise) =#
static_stab_coeff(ps::APS) = ps.nt.c_b
l_max(ps::APS) = ps.nt.l_max

end
