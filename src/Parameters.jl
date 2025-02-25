"""
    Parameters

"""
module Parameters

import Thermodynamics as TD
import SurfaceFluxes as SF
import CloudMicrophysics as CM

abstract type AbstractTurbulenceConvectionParameters end
const ATCP = AbstractTurbulenceConvectionParameters


"""
Because we pack our parameters in a named tuple and put symbols in Vals so they're isbits, we use a convenience get fcn.
Put in here so can use here and in driver
"""

 # like SciMLBase._unwrap_val with _unwrap_val(::Val{B}) where {B} = B and _unwrap_val(B) = B
 unwrap_val(::Val{B}) where {B} = B 
 unwrap_val(B) = B # idk if we should keep this definition, seems error prone

 
function get_isbits_nt(named_tuple::NamedTuple, symbol::Symbol, default)

    if symbol in keys(named_tuple)
        val = Base.getproperty(named_tuple, symbol)
    else
        val = default
    end

    # val = get(named_tuple, symbol, default)

    # if isa(val, Val) # extract from it's Val wrapped prison
    #     val = typeof(val).parameters[1] # this is a hack to get the case name from the Val object
    # end
    # return val
    return unwrap_val(val)
end





function get_isbits_nt(named_tuple::NamedTuple, symbol::Symbol)

    val = Base.getproperty(named_tuple, symbol)

    # if isa(val, Val) # extract from it's Val wrapped prison
    #     val = typeof(val).parameters[1] # this is a hack to get the case name from the Val object 
    #     # could also do something like SciMLBase._unwrap_val with _unwrap_val(::Val{B}) where {B} = B and _unwrap_val(B) = B
    # end
    # return val
    return unwrap_val(val)
end



"""
Just lik get_isbits_nt, but warns if using the default (just like parse_namelist() does)
Probably slower, so only use e.g. in types.jl when building our objects.
"""
function parse_isbits_nt(named_tuple::NamedTuple, symbol::Symbol, default)

    if symbol in keys(named_tuple)
        val = Base.getproperty(named_tuple, symbol)
    else
        val = default
        @info "Using default value, $default, for parameter $symbol."
    end

    # if isa(val, Val) # extract from it's Val wrapped prison
    #     val = typeof(val).parameters[1] # this is a hack to get the case name from the Val object
    # end
    # return val
    return unwrap_val(val)
end

# ---------------------------------------------------------------------- #

# # Base.@kwdef struct UserParameters{FT} <: ATCP # Things i wanna have just for me... (so i dont have to keep reaching into user_params for things that may or may not be there)
# #     particle_min_radius::FT
# #     # add more here as needed
# # end

# # If this new UserParameters works, maybe we can get rid of / combine user_params and user_args
# # still have to check on isbits ness... rn user_params / user_args are namedtuples, now this would be a struct containing a named tuple
# # iirc i used the dict/nt mix because of strings or something i found hard to generate? but now it just seems unnecessary... it also allowed us to go through and make sure all the values are isbits

# # Based on https://discourse.julialang.org/t/creating-a-struct-from-a-named-tuple/94586/6?u=jbphyswx 
# Base.@kwdef struct UserParameters{NT <: NamedTuple} <: ATCP # Things i wanna have just for me... (so i dont have to keep reaching into user_params for things that may or may not be there)
#     nt::NT
# end
# Base.getproperty(x::UserParameters, y::Symbol) = get_isbits_nt(getfield(x,:nt), y) # overload getproperty to access the value in the namedtuple. Idk how fragile this is...
# # Base.getproperty(x::UserParameters, y::Symbol) = Base.getproperty(getfield(x,:nt),y) # overload getproperty to access the value in the namedtuple. Idk how fragile this is...
# get_isbits_nt(UP::UserParameters, symbol::Symbol) = get_isbits_nt(UP.nt, symbol)
# get_isbits_nt(UP::UserParameters, symbol::Symbol, default) = get_isbits_nt(UP.nt, symbol, default)

# # ---------------------------------------------------------------------- #


#####
##### TurbulenceConvection parameters
#####

Base.@kwdef struct TurbulenceConvectionParameters{FT, MP, SFP, NT <: NamedTuple, UP <: NamedTuple} <: ATCP
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
    user_params::UP # Let this be a NamedTuple created from namelist["user_params"]::Dict in the main program.. to ensure this is `isbits`, the values in the dict should have been be turned to tuples if they were arrays... string keys should have been  converted to symbols (get_parameter_valus() in parameter_set.jl did a similar thing )
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
