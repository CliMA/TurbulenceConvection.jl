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
area_limiter_scale(ps::APS) = ps.nt.γ_lim
area_limiter_power(ps::APS) = ps.nt.β_lim
static_stab_coeff(ps::APS) = ps.nt.c_b
l_max(ps::APS) = ps.nt.l_max

"""Additional entrainment/detrainment parameters."""
c_gen(ps::APS) = ps.nt.c_gen

""" stochastic parameters """
stoch_ε_lognormal_var(ps::APS) = ps.nt.stoch_ε_lognormal_var
stoch_δ_lognormal_var(ps::APS) = ps.nt.stoch_δ_lognormal_var

sde_ϵ_θ(ps::APS) = ps.nt.sde_ϵ_θ
sde_ϵ_σ(ps::APS) = ps.nt.sde_ϵ_σ
sde_δ_θ(ps::APS) = ps.nt.sde_δ_θ
sde_δ_σ(ps::APS) = ps.nt.sde_δ_σ

end
