import TurbulenceConvection as TC
import CLIMAParameters as CP
import CloudMicrophysics as CM
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import TurbulenceConvection.Parameters as TCP


leaves_to_tuples(x) = x
leaves_to_tuples(x::Array) = Tuple(x)
leaves_to_tuples(x::Dict) = Dict{String, Any}((k => leaves_to_tuples(v) for (k,v) in x))

leaves_to_Vals(x) = x # this turns leaves that are string into symbols and then all symbols into Vals so that they can be isbits, thus we are allowed to pass in symbols or strings, to extract, use typeof(valitem).parameters[1]
leaves_to_Vals(x::Symbol) = Val(x)
leaves_to_Vals(x::String) = Val(Symbol(x))
leaves_to_Vals(x::Dict) = Dict{String, Any}((k => leaves_to_Vals(v) for (k,v) in x))
leaves_to_Vals(x::NamedTuple) = map(leaves_to_Vals, x)

namedtuple_fromdict(x) = x 
namedtuple_fromdict(d::Dict) = (; (Symbol(k) => namedtuple_fromdict(v) for (k,v) in d)...) # from https://discourse.julialang.org/t/how-to-make-a-named-tuple-from-a-dictionary/10899/46?u=jbphyswx

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
    pow_icenuc = TC.parse_namelist(namelist, "microphysics", "pow_icenuc"; default = 1e7) # I think this has to be here to overwrite toml..., and to place it so we can overwrite_namelist it .. picked high value initially to keep ramp but maybe should keep original default... https://github.com/CliMA/CLIMAParameters.jl/blob/2f298661d527c1b9a11188ee2abc92ba4ba0c5ec/src/parameters.toml#L558
    r_ice_snow = TC.parse_namelist(namelist, "microphysics", "r_ice_snow"; default = 62.5e-6) # allow changing this in the namelist/param_set/microphys_params instead of user_params (probably not ideal if you used a 2-moment or something where this mattered more but ...)
    #
    χm_ice = TC.parse_namelist(namelist, "microphysics", "χm_ice"; default = 1.0)
    #
    ν_sno = TC.parse_namelist(namelist, "microphysics", "ν_sno"; default = 0.63)
    μ_sno = TC.parse_namelist(namelist, "microphysics", "μ_sno"; default = 4.36e9)
    me_sno = TC.parse_namelist(namelist, "microphysics", "me_sno"; default = 2.0)
    χm_sno = TC.parse_namelist(namelist, "microphysics", "χm_sno"; default = 1.0)
    #
    χv_rai = TC.parse_namelist(namelist, "microphysics", "χv_rai"; default = 1.0)
    χv_sno = TC.parse_namelist(namelist, "microphysics", "χv_sno"; default = 1.0)

    # Override the default files in the toml file with the values from the namelist -- [[ i think it's plain language name name, alias = code_name ]]
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
        println(io, "[ice_snow_threshold_radius]") # allow chaninging r_ice_snow in the namelist (probably not how you wanna do it in Blk1mVel, etc, but we're not using)
        println(io, "alias = \"r_ice_snow\"")
        println(io, "value = " * string(r_ice_snow))
        println(io, "type = \"float\"")
        println(io, "[cloud_ice_mass_size_relation_coefficient_chim]") # see https://github.com/CliMA/ClimaParams.jl/blob/70c847bd01a2963a1bbddbfc283bbcdc306888d3/src/parameters.toml#L754C2-L754C47
        println(io, "alias = \"χm_ice\"") 
        println(io, "value = " * string(χm_ice))
        println(io, "type = \"float\"")
        println(io, "[snow_flake_size_distribution_coefficient_nu]") # ν_sno
        println(io, "alias = \"ν_sno\"")
        println(io, "value = " * string(ν_sno))
        println(io, "type = \"float\"")
        println(io, "[snow_flake_size_distribution_coefficient_mu]") # μ_sno prefactor
        println(io, "alias = \"μ_sno\"")
        println(io, "value = " * string(μ_sno))
        println(io, "type = \"float\"")
        println(io, "[snow_mass_size_relation_coefficient_me]")
        println(io, "alias = \"me_sno\"")
        println(io, "value = " * string(me_sno))
        println(io, "type = \"float\"")
        println(io, "[snow_mass_size_relation_coefficient_chim]") # χm_sno
        println(io, "alias = \"χm_sno\"")
        println(io, "value = " * string(χm_sno))
        println(io, "type = \"float\"")
        println(io, "[rain_terminal_velocity_size_relation_coefficient_chiv]") # χv_rai
        println(io, "alias = \"χv_rai\"")
        println(io, "value = " * string(χv_rai))
        println(io, "type = \"float\"")
        println(io, "[snow_terminal_velocity_size_relation_coefficient_chiv]") # χv_sno
        println(io, "alias = \"χv_sno\"")
        println(io, "value = " * string(χv_sno))
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

    # CM 0.13
    # see  https://github.com/CliMA/CloudMicrophysics.jl/pull/193/files#diff-2624a4f6514b5d8d4cdf38f495c23a1167682c3accbdda813b80b425a02ab2f1
    # not sure if I should really pivot and change `import CLIMAParameters.Parameters` to `CMP` by using `import CloudMicrophysics.Parameters as CMP``
    # aliases = string.(fieldnames(CM.Parameters.ModalNucleationParameters))
    # pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    # modal_nucleation_params = CM.Parameters.ModalNucleationParameters{FT}(; pairs...)
    # MNP = typeof(modal_nucleation_params)

    # aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    # pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    # microphys_params = CM.Parameters.CloudMicrophysicsParameters{FT, TP, MNP}(;
    # pairs...,
    # thermo_params,
    # modal_nucleation_params,
    # )

    # CM 0.14
    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    microphys_params = CM.Parameters.CloudMicrophysicsParameters{FT, TP}(;
    pairs...,
    thermo_params,
    )

    #
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
    "planet_radius",]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "TurbulenceConvection")

    SFP = typeof(surf_flux_params)



    namelist_user_args = get!(namelist, "user_args", Dict()) # handle if empty so we don't have to add to namelist defaults
    namelist_user_params = get!(namelist, "user_params", Dict())  # handle if empty so we don't have to add to namelist defaults, put in nameslist so it'll be there for storage
    
    # relaxation_timescale_params = get(namelist, "relaxation_timescale_params", Dict()) # handle if empty so we don't have to add to namelist defaults
    # merge relaxation_timescale_params into user_params (in the future we could handle separately , but this is backwards compatible)
    # the relaxation_timescale_params go into relaxation_timescale_types so they shouldn't actually need to be stored but oh well
    # user_params = merge(user_params, relaxation_timescale_params) # it's no longer just a reference into namelist, so be careful
    # namelist["user_params"] = user_params # put it back into namelist so we can use it later


    # honestly now that we are creating object types w/ the relaxation params, we really don't need them to be stored in user_params in the param_set, maybe we should move them all to like namelist["relaxation_timescale_parms"] or something ...
   

    # user_params = get(namelist, "user_params", Dict())
    # overwrites, but this is all local so no need for TOML
    namelist_user_params["particle_min_radius"] =  TC.parse_namelist(namelist, "user_params", "particle_min_radius"; default = 0.2e-6) # this is set no matter what
    # namelist_user_params["r_liq_rain"] =  TC.parse_namelist(namelist, "user_params", "r_liq_rain"; default = 100e-6) # this is set no matter what [ We take the default form the cutoff for raindrop in Chen]
    namelist_user_params["r_liq_rain"] =  TC.parse_namelist(namelist, "user_params", "r_liq_rain"; default = 30e-6) # I think literature likes to stay with between 20-50 microns... https://ntrs.nasa.gov/api/citations/20220009419/downloads/SSantosJAMESLimitationsReprint.pdf used 40
    namelist_user_params["χm_liq"] =  TC.parse_namelist(namelist, "user_params", "χm_liq"; default = 1.0) # this is set no matter what

    namelist_user_params["r_ice_acnv_scaling_factor"] =  TC.parse_namelist(namelist, "user_params", "r_ice_acnv_scaling_factor"; default = 1.0) # this is set no matter what
    namelist_user_params["r_ice_snow_threshold_scaling_factor"] =  TC.parse_namelist(namelist, "user_params", "r_ice_snow_threshold_scaling_factor"; default = 200/125) # this is set no matter what

    namelist_user_params["τ_acnv_sno_threshold"] =  TC.parse_namelist(namelist, "user_params", "τ_acnv_sno_threshold"; default = 100.0) # this is set no matter waht
    namelist_user_params["τ_acnv_liq_thresh"] =  TC.parse_namelist(namelist, "user_params", "τ_acnv_liq_thresh"; default = 100.0) # this is set no matter what

    namelist_user_params["massflux_N_i_boost_factor"] =  TC.parse_namelist(namelist, "user_params", "massflux_N_i_boost_factor"; default = 2.0) # this is set no matter what
    namelist_user_params["sedimentation_N_i_boost_factor"] =  TC.parse_namelist(namelist, "user_params", "sedimentation_N_i_boost_factor"; default = 0.2) # this is set no matter what
    namelist_user_params["apply_massflux_N_i_boost"] =  TC.parse_namelist(namelist, "user_params", "apply_massflux_N_i_boost"; default = false) # this is set no matter what
    namelist_user_params["apply_sedimentation_N_i_boost"] =  TC.parse_namelist(namelist, "user_params", "apply_sedimentation_N_i_boost"; default = false) # this is set no matter what

    namelist_user_params["massflux_N_i_boost_max_ratio"] =  TC.parse_namelist(namelist, "user_params", "massflux_N_i_boost_max_ratio"; default = 0.7) # this is set no matter what
    namelist_user_params["massflux_N_i_boost_progress_fraction"] =  TC.parse_namelist(namelist, "user_params", "massflux_N_i_boost_progress_fraction"; default = 0.5) # this is set no matter what
    
    namelist_user_params["use_ice_mult"] =  TC.parse_namelist(namelist, "user_params", "use_ice_mult"; default = true) # this is set no matter what

    # If we decide to add this to the calibration... 
    # default_r_factor_liq = FT(1) # unknown
    # default_r_factor_ice = FT(((4/3) / 8)^(1/3))
    # namelist_user_params["mean_r_factor_liq"] =  TC.parse_namelist(namelist, "user_params", "mean_r_factor_liq"; default = FT(1)) # this is set no matter what
    # namelist_user_params["mean_r_factor_ice"] =  TC.parse_namelist(namelist, "user_params", "mean_r_factor_ice"; default = FT(1)) # this is set no matter what

    namelist_user_params["q_min"] = TC.parse_namelist(namelist, "user_params", "q_min"; default = zero(FT)) # maybe one day swap to var limiter step
    # user_params = namelist["user_params"] # Let this be a Dict()


    user_params = deepcopy(namelist_user_params) # make a copy of user_params so we can modify it without changing the original
    user_args = deepcopy(namelist_user_args) # make a copy of user_args so we can modify it without changing the original
    # delete things stored in relaxation timescales
    delete!.(Ref(user_params), ["neural_microphysics_relaxation_network", "model_re_location", "model_x_0_characteristic"]) # drop neural_microphysics_relaxation_network, model_re_location, model_x_0_characteristic, etc. from user_param.  deletes from user_params, but not from namelist, so we can use it later, so it doesnt end up in the param_set
    delete!.(Ref(user_args), ["nonequilibrium_moisture_scheme", "adjust_liq_N", "adjust_ice_N", "use_heterogeneous_ice_nucleation",]) # these are now in relaxtion_timescale object
    delete!.(Ref(user_params), ["heterogeneous_ice_nucleation_coefficient", "heterogeneous_ice_nucleation_exponent", "min_τ_liq", "min_τ_ice", "max_τ_liq", "max_τ_ice", "min_N_liq", "min_N_ice", "max_N_liq", "max_N_ice"]) # these are now in relaxtion_timescale object

    # delete things stored in cloud_sedimentation_model (or in its setup)
    delete!.(Ref(user_params), ["liq_sedimentation_Dmax", "ice_sedimentation_Dmax", "liq_sedimentation_scaling_factor", "ice_sedimentation_scaling_factor", "E_liq_ice"])
    delete!.(Ref(user_args), ["use_sedimentation", "liq_terminal_velocity_scheme", "ice_terminal_velocity_scheme", "sedimentation_differencing_scheme", "grid_mean", ]) 

    # delete things stored in precip_model
    delete!.(Ref(user_params), ["rain_sedimentation_scaling_factor", "snow_sedimentation_scaling_factor"])
    delete!.(Ref(user_args), ["rain_terminal_velocity_scheme", "snow_terminal_velocity_scheme",]) 

    # delete things stored in precip_formation model
    delete!.(Ref(user_params), ["ice_dep_acnv_scaling_factor", "ice_dep_acnv_scaling_factor_above", "ice_acnv_power",]) # "r_ice_acnv_scaling_factor", "r_ice_snow_threshold_scaling_factor"

    # delete things stored in tendency limiter set
    delete!.(Ref(user_args), ["truncated_basic_limiter_factor", "default_tendency_limiter_type", "fallback_default_tendency_limiter_type", "nonequilibrium_moisture_sources_limiter_type", "fallback_nonequilibrium_moisture_sources_limiter_type", "fallback_to_standard_supersaturation_limiter", "entr_detr_limiter_type", "fallback_entr_detr_limiter_type", "precipitation_tendency_limiter_type", "fallback_precipitation_tendency_limiter_type", "tendency_resolver_setup"])


    # delete stored area partition things

    # delete things not stored anywhere but that we don't need
    # delete!.(Ref(user_args, [
    #     "", # we make a cloud sedimentation model...
    #     ])
    # )

    user_args = leaves_to_tuples(user_args) # convert any leaves that are Arrays to tuples to preserve isbits
    user_args = leaves_to_Vals(user_args) # convert leaves that are strings or symbols into Vals to preserve isbits
    user_args = namedtuple_fromdict(user_args) # convert dict to NamedTuple to preserve isbits

    user_params = leaves_to_tuples(user_params) # convert any leaves that are Arrays to tuples to preserve isbits
    user_params = leaves_to_Vals(user_params)  # convert leaves that are strings or symbols into Vals to preserve isbits
    user_params = namedtuple_fromdict(user_params) # convert dict to NamedTuple to preserve isbits

    # ------------------------------------------------ #

    # # Do we actually need to wrap the type if it's a nt anyway...
    # # Create user_params, see Parameters.jl for what we store in here...
    # user_params = get(namelist, "user_params", Dict()) # should be a dict...
    # # overwrites, not using TOML
    # # ensure isbitsness
    # user_params = leaves_to_tuples(user_params) # convert leaves that are strings or symbols into Vals to preserve isbits
    # user_params = leaves_to_Vals(user_params) # convert leaves that are strings or symbols into Vals to preserve isbits
    # # set the user_params
    # user_params = namedtuple_fromdict(user_args)
    # user_params = TCP.UserParameters{typeof(user_args)}(user_params) # version that store NamedTuple

    # # I also used user_args and user_params to separate parameters like process_x_parameter_1::FT=0.1 and keyword arguments like use_x_process::Bool = true, etc.


    # Some things in param_set are used in the construction of other objects, is kind of redundant lol. No way to fix at the moment I guess given the namelist isn't passed everywhere...

    # ------------------------------------------------ #

    # ------------------------------------------------ #

    #= [[ deprecate `user_args` ]]
        user_args was originally created to store arguments that might be useful in the models evolution. it still serves that purpose in the namelist.
        However, arguments need not be passed as parameters. The are helpful in construction, in types.jl. 
        If they are important for runtime, either place them into user_params or into an specific object that is part of edmf.

    =#

    param_set = TCP.TurbulenceConvectionParameters{FTD, MP, SFP, typeof(user_params)}(; pairs..., microphys_params, surf_flux_params, user_params) # `typeof` so we have concrete types such that if user_args and user_params are isbits, then param_set is isbits
    if !isbits(param_set)
        @info("param_set", param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on
