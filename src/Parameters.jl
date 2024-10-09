"""
    Parameters

"""
module Parameters

import Thermodynamics as TD
import SurfaceFluxes as SF
import CloudMicrophysics as CM

abstract type AbstractTurbulenceConvectionParameters end
const ATCP = AbstractTurbulenceConvectionParameters

#####
##### TurbulenceConvection parameters
#####

Base.@kwdef struct TurbulenceConvectionParameters{FT, MP, SFP, NT, NT2} <: ATCP
    Omega::FT
    planet_radius::FT
    microph_scaling::FT
    microph_scaling_dep_sub::FT
    microph_scaling_melt::FT
    microph_scaling_acnv::FT
    microph_scaling_accr::FT
    microphys_params::MP
    surf_flux_params::SFP
    user_args::NT # Let this be a Namedtuple of user arguments passed taken straight from namelist["user_args"]::NamedTuple in the main program
    user_aux::NT2 # Let this be a NamedTuple created from namelist["user_aux"]::Dict in the main program.. to ensure this is `isbits`, the values in the dict should have been be turned to tuples if they were arrays... string keys should have been  converted to symbols (get_parameter_valus() in parameter_set.jl did a similar thing )
end

thermodynamics_params(ps::ATCP) = CM.Parameters.thermodynamics_params(ps.microphys_params)
surface_fluxes_params(ps::ATCP) = ps.surf_flux_params
microphysics_params(ps::ATCP) = ps.microphys_params

Base.eltype(::TurbulenceConvectionParameters{FT}) where {FT} = FT
Omega(ps::ATCP) = ps.Omega
planet_radius(ps::ATCP) = ps.planet_radius
# TODO - microph_scaling is the factor for adjusting evaporation.
# The name will be fixed in CLIMAParameters first.
microph_scaling(ps::ATCP) = ps.microph_scaling
microph_scaling_dep_sub(ps::ATCP) = ps.microph_scaling_dep_sub
microph_scaling_melt(ps::ATCP) = ps.microph_scaling_melt
microph_scaling_acnv(ps::ATCP) = ps.microph_scaling_acnv
microph_scaling_accr(ps::ATCP) = ps.microph_scaling_accr

#####
##### Forwarding parameters
#####

##### Forwarding Thermodynamics.jl

const TDPS = TD.Parameters.ThermodynamicsParameters
for var in fieldnames(TDPS)
    @eval $var(ps::ATCP) = TD.Parameters.$var(thermodynamics_params(ps))
end

# derived parameters
molmass_ratio(ps::ATCP) = TD.Parameters.molmass_ratio(thermodynamics_params(ps))
R_d(ps::ATCP) = TD.Parameters.R_d(thermodynamics_params(ps))
R_v(ps::ATCP) = TD.Parameters.R_v(thermodynamics_params(ps))
cp_d(ps::ATCP) = TD.Parameters.cp_d(thermodynamics_params(ps))
cv_v(ps::ATCP) = TD.Parameters.cv_v(thermodynamics_params(ps))
cv_l(ps::ATCP) = TD.Parameters.cv_l(thermodynamics_params(ps))

##### Forwarding SurfaceFluxes.jl

von_karman_const(ps::ATCP) = SF.Parameters.von_karman_const(surface_fluxes_params(ps))

##### Forwarding CloudMicrophysics.jl

ρ_cloud_liq(ps::ATCP) = CM.Parameters.ρ_cloud_liq(microphysics_params(ps))

end
