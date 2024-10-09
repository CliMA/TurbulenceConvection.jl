import TurbulenceConvection as TC
import CLIMAParameters as CP
import CloudMicrophysics as CM
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import TurbulenceConvection.Parameters as TCP

#! format: off
function create_parameter_set(
    namelist,
    toml_dict_default::CP.AbstractTOMLDict,
    FTD = CP.float_type(toml_dict_default)
)
    FT = CP.float_type(toml_dict_default)
    _, out_dir = nc_fileinfo(namelist)
    override_file = joinpath(out_dir, "override_dict.toml")

    # Read in data from namelist to overwrite parameters
    τ_precip = TC.parse_namelist(namelist, "microphysics", "τ_precip"; default = 1000.0)
    τ_cond_evap = TC.parse_namelist(namelist, "microphysics", "τ_cond_evap"; default = 10.0)
    τ_sub_dep = TC.parse_namelist(namelist, "microphysics", "τ_sub_dep"; default = 10.0)
    τ_acnv_rai = TC.parse_namelist(namelist, "microphysics", "τ_acnv_rai"; default = 2500.0)
    τ_acnv_sno = TC.parse_namelist(namelist, "microphysics", "τ_acnv_sno"; default = 100.0)
    q_liq_threshold = TC.parse_namelist(namelist, "microphysics", "q_liq_threshold"; default = 0.5e-3)
    q_ice_threshold = TC.parse_namelist(namelist, "microphysics", "q_ice_threshold"; default = 1e-6)
    microph_scaling_acnv = TC.parse_namelist(namelist, "microphysics", "microph_scaling_acnv"; default = 1.0)
    microph_scaling_accr = TC.parse_namelist(namelist, "microphysics", "microph_scaling_accr"; default = 1.0)
    microph_scaling = TC.parse_namelist(namelist, "microphysics", "microph_scaling"; default = 1.0)
    microph_scaling_dep_sub = TC.parse_namelist(namelist, "microphysics", "microph_scaling_dep_sub"; default = 1.0)
    microph_scaling_melt = TC.parse_namelist(namelist, "microphysics", "microph_scaling_melt"; default = 1.0)
    E_liq_rai = TC.parse_namelist(namelist, "microphysics", "E_liq_rai"; default = 0.8)
    E_liq_sno = TC.parse_namelist(namelist, "microphysics", "E_liq_sno"; default = 0.1)
    E_ice_rai = TC.parse_namelist(namelist, "microphysics", "E_ice_rai"; default = 1.0)
    E_ice_sno = TC.parse_namelist(namelist, "microphysics", "E_ice_sno"; default = 0.1)
    E_rai_sno = TC.parse_namelist(namelist, "microphysics", "E_rai_sno"; default = 1.0)
    A_acnv_KK2000 = TC.parse_namelist(namelist, "microphysics", "A_acnv_KK2000"; default = 7.42e13)
    a_acnv_KK2000 = TC.parse_namelist(namelist, "microphysics", "a_acnv_KK2000"; default = 2.47)
    b_acnv_KK2000 = TC.parse_namelist(namelist, "microphysics", "b_acnv_KK2000"; default = -1.79)
    c_acnv_KK2000 = TC.parse_namelist(namelist, "microphysics", "c_acnv_KK2000"; default = -1.47)
    pow_icenuc=TC.parse_namelist(namelist, "microphysics", "pow_icenuc"; default = 1e7) # I think this has to be here to overwrite toml..., and to place it so we can overwrite_namelist it .. picked high value initially to keep ramp but maybe should keep original default... https://github.com/CliMA/CLIMAParameters.jl/blob/2f298661d527c1b9a11188ee2abc92ba4ba0c5ec/src/parameters.toml#L558

    # Override the default files in the toml file
    open(override_file, "w") do io
        println(io, "[mean_sea_level_pressure]")
        println(io, "alias = \"MSLP\"")
        println(io, "value = 100000.0")
        println(io, "type = \"float\"")
        println(io, "[precipitation_timescale]")
        println(io, "alias = \"τ_precip\"")
        println(io, "value = " * string(τ_precip))
        println(io, "type = \"float\"")
        println(io, "[condensation_evaporation_timescale]")
        println(io, "alias = \"τ_cond_evap\"")
        println(io, "value = " * string(τ_cond_evap))
        println(io, "type = \"float\"")
        println(io, "[sublimation_deposition_timescale]")
        println(io, "alias = \"τ_sub_dep\"")
        println(io, "value = " * string(τ_sub_dep))
        println(io, "type = \"float\"")
        println(io, "[rain_autoconversion_timescale]")
        println(io, "alias = \"τ_acnv_rai\"")
        println(io, "value = " * string(τ_acnv_rai))
        println(io, "type = \"float\"")
        println(io, "[snow_autoconversion_timescale]")
        println(io, "alias = \"τ_acnv_sno\"")
        println(io, "value = " * string(τ_acnv_sno))
        println(io, "type = \"float\"")
        println(io, "[cloud_liquid_water_specific_humidity_autoconversion_threshold]")
        println(io, "alias = \"q_liq_threshold\"")
        println(io, "value = " * string(q_liq_threshold))
        println(io, "type = \"float\"")
        println(io, "[cloud_ice_specific_humidity_autoconversion_threshold]")
        println(io, "alias = \"q_ice_threshold\"")
        println(io, "value = " * string(q_ice_threshold))
        println(io, "type = \"float\"")
        println(io, "[microph_scaling_acnv]")
        println(io, "alias = \"microph_scaling_acnv\"")
        println(io, "value = " * string(microph_scaling_acnv))
        println(io, "type = \"float\"")
        println(io, "[microph_scaling_accr]")
        println(io, "alias = \"microph_scaling_accr\"")
        println(io, "value = " * string(microph_scaling_accr))
        println(io, "type = \"float\"")
        println(io, "[microph_scaling]")
        println(io, "alias = \"microph_scaling\"")
        println(io, "value = " * string(microph_scaling))
        println(io, "type = \"float\"")
        println(io, "[microph_scaling_dep_sub]")
        println(io, "alias = \"microph_scaling_dep_sub\"")
        println(io, "value = " * string(microph_scaling_dep_sub))
        println(io, "type = \"float\"")
        println(io, "[microph_scaling_melt]")
        println(io, "alias = \"microph_scaling_melt\"")
        println(io, "value = " * string(microph_scaling_melt))
        println(io, "type = \"float\"")
        println(io, "[cloud_liquid_rain_collision_efficiency]")
        println(io, "alias = \"E_liq_rai\"")
        println(io, "value = " * string(E_liq_rai))
        println(io, "type = \"float\"")
        println(io, "[cloud_liquid_snow_collision_efficiency]")
        println(io, "alias = \"E_liq_sno\"")
        println(io, "value = " * string(E_liq_sno))
        println(io, "type = \"float\"")
        println(io, "[cloud_ice_rain_collision_efficiency]")
        println(io, "alias = \"E_ice_rai\"")
        println(io, "value = " * string(E_ice_rai))
        println(io, "type = \"float\"")
        println(io, "[cloud_ice_snow_collision_efficiency]")
        println(io, "alias = \"E_ice_sno\"")
        println(io, "value = " * string(E_ice_sno))
        println(io, "type = \"float\"")
        println(io, "[rain_snow_collision_efficiency]")
        println(io, "alias = \"E_rai_sno\"")
        println(io, "value = " * string(E_rai_sno))
        println(io, "type = \"float\"")
        println(io, "[KK2000_auctoconversion_coeff_A]")
        println(io, "alias = \"A_acnv_KK2000\"")
        println(io, "value = " * string(A_acnv_KK2000))
        println(io, "type = \"float\"")
        println(io, "[KK2000_auctoconversion_coeff_a]")
        println(io, "alias = \"a_acnv_KK2000\"")
        println(io, "value = " * string(a_acnv_KK2000))
        println(io, "type = \"float\"")
        println(io, "[KK2000_auctoconversion_coeff_b]")
        println(io, "alias = \"b_acnv_KK2000\"")
        println(io, "value = " * string(b_acnv_KK2000))
        println(io, "type = \"float\"")
        println(io, "[KK2000_auctoconversion_coeff_c]")
        println(io, "alias = \"c_acnv_KK2000\"")
        println(io, "value = " * string(c_acnv_KK2000))
        println(io, "type = \"float\"")
        println(io, "[pow_icenuc]") # add in our overwrite above of the original toml from ClimaParameters
        println(io, "alias = \"pow_icenuc\"")
        println(io, "value = " * string(pow_icenuc))
        println(io, "type = \"float\"")
    end

    toml_dict = CP.create_toml_dict(FT; override_file, dict_type="alias")
    isfile(override_file) && rm(override_file; force=true)

    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FTD}(; param_pairs...)
    # logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    # CP.log_parameter_information(toml_dict, logfilepath)
    TP = typeof(thermo_params)

    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    aliases = setdiff(aliases, ["thermo_params"])
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    microphys_params = CM.Parameters.CloudMicrophysicsParameters{FTD, TP}(;
        pairs...,
        thermo_params,
    )
    MP = typeof(microphys_params)

    aliases = ["Pr_0_Businger", "a_m_Businger", "a_h_Businger", "ζ_a_Businger", "γ_Businger"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (; Pr_0 = pairs.Pr_0_Businger, a_m = pairs.a_m_Businger, a_h = pairs.a_h_Businger, ζ_a = pairs.ζ_a_Businger, γ = pairs.γ_Businger)
    ufp = UF.BusingerParams{FTD}(; pairs...)
    UFP = typeof(ufp)

    pairs = CP.get_parameter_values!(toml_dict, ["von_karman_const"], "SurfaceFluxesParameters")
    surf_flux_params = SF.Parameters.SurfaceFluxesParameters{FTD, UFP, TP}(; pairs..., ufp, thermo_params)

    aliases = [
    "microph_scaling_dep_sub",
    "microph_scaling_melt",
    "microph_scaling",
    "microph_scaling_acnv",
    "microph_scaling_accr",
    "Omega",
    "planet_radius"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "TurbulenceConvection")

    SFP = typeof(surf_flux_params)
    namelist["user_args"] = get(namelist,"user_args", (;)) # handle if empty so we don't have to add to namelist defaults
    namelist["user_aux"] = get(namelist,"user_aux", (;)) # handle if empty so we don't have to add to namelist defaults

    user_args = namelist["user_args"] # Let this be a NamedTuple
    user_aux  = namelist["user_aux"] # Let this be a Dict()

    leaves_to_tuples(x) = x
    leaves_to_tuples(x::Array) = Tuple(x)
    leaves_to_tuples(x::Dict) = Dict{String, Any}((k => leaves_to_tuples(v) for (k,v) in x))
    user_aux = leaves_to_tuples(user_aux) # convert any leaves that are Arrays to tuples to preserve isbits
    user_args = leaves_to_tuples(user_args) # convert any leaves that are Arrays to tuples to preserve isbits

    leaves_to_Vals(x) = x # this turns leaves that are string into symbols and then all symbols into Vals so that they can be isbits, thus we are allowed to pass in symbols or strings, to extract, use typeof(valitem).parameters[1]
    leaves_to_Vals(x::Symbol) = Val(x)
    leaves_to_Vals(x::String) = Val(Symbol(x))
    leaves_to_Vals(x::Dict) = Dict{String, Any}((k => leaves_to_Vals(v) for (k,v) in x))
    user_aux = leaves_to_Vals(user_aux)  # convert leaves that are strings or symbols into Vals to preserve isbits
    user_args = leaves_to_Vals(user_args) # convert leaves that are strings or symbols into Vals to preserve isbits

    namedtuple_fromdict(x) = x 
    namedtuple_fromdict(d::Dict) = (; (Symbol(k) => namedtuple_fromdict(v) for (k,v) in d)...) # from https://discourse.julialang.org/t/how-to-make-a-named-tuple-from-a-dictionary/10899/46?u=jbphyswx
    user_aux = namedtuple_fromdict(user_aux) # convert dict to NamedTuple to preserve isbits
    user_args = namedtuple_fromdict(user_args) # convert dict to NamedTuple to preserve isbits

    param_set = TCP.TurbulenceConvectionParameters{FTD, MP, SFP, typeof(user_args), typeof(user_aux)}(; pairs..., microphys_params, surf_flux_params, user_args, user_aux) # `typeof` so we have concrete types such that if user_args and user_aux are isbits, then param_set is isbits
    if !isbits(param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on