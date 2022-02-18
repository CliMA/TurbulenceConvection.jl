"""
    ClimaParams

This is a module for staging experimental
clima parameters-- that we are unsure will
be used in the final model. Once they are
no longer experimental, they should be moved
to CLIMAParameters.jl.

TODO: move non-experimental parameters to CLIMAParameters
"""
module ClimaParams

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractEarthParameterSet

""" divergence factor for bubble case (zero otherwise) """
entrainment_massflux_div_factor(ps::APS) = ps.nt.c_div
area_limiter_scale(ps::APS) = ps.nt.γ_lim
area_limiter_power(ps::APS) = ps.nt.β_lim
static_stab_coeff(ps::APS) = ps.nt.c_b
l_max(ps::APS) = ps.nt.l_max

"""Additional entrainment/detrainment parameters."""
c_gen(ps::APS) = ps.nt.c_gen
c_fno(ps::APS) = ps.nt.c_fno
c_rf_fix(ps::APS) = ps.nt.c_rf_fix
c_rf_opt(ps::APS) = ps.nt.c_rf_opt

""" stochastic parameters """
c_gen_stoch(ps::APS) = ps.nt.c_gen_stoch

end
