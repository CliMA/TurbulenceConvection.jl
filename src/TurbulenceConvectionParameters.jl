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

end
