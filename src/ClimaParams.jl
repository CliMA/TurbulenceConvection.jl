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

using CLIMAParameters: AbstractEarthParameterSet
const APS = AbstractEarthParameterSet

""" divergence factor for bubble case (zero otherwise) """
entrainment_massflux_div_factor(ps::APS) = ps.nt.c_div
entrainment_sigma(ps::APS) = ps.nt.μ
static_stab_coeff(ps::APS) = ps.nt.c_b
l_max(ps::APS) = ps.nt.l_max


""" stochastic parameters """
stoch_closure(ps::APS) = ps.nt.stoch_closure
stoch_ε_lognormal_var(ps::APS) = ps.nt.stoch_ε_lognormal_var
stoch_δ_lognormal_var(ps::APS) = ps.nt.stoch_δ_lognormal_var

end
