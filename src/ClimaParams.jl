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
entrainment_sigma(ps::APS) = ps.nt.μ
static_stab_coeff(ps::APS) = ps.nt.c_b
l_max(ps::APS) = ps.nt.l_max

end
