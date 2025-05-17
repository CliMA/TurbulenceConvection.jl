FT = Float64

reload_environment = false
if reload_environment || !isdefined(Main, :TurbulenceConvection) || !isdefined(Main, :main1d)
    using Revise
    using TurbulenceConvection
    tc = pkgdir(TurbulenceConvection)
    # include("/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl/src/TurbulenceConvection.jl")
    # tc = "/home/jbenjami/Research_Schneider/CliMA/TurbulenceConvection.jl"
    # includet(joinpath(tc, "driver", "main_revise.jl"))
    # includet(joinpath(tc, "driver", "generate_namelist.jl"))
    includet(joinpath(tc, "driver", "main.jl"))
    includet(joinpath(tc, "driver", "generate_namelist.jl"))
end

reload_param_set = false
if reload_param_set | !isdefined(Main, :param_set)
    using Pkg
    using Revise
    import CLIMAParameters as CP # use CLIMAParameters = "0.7, 0.8, 0.9, 0.10"
    # import ClimaParams as CPP # would using this trouble w/ TC.jl? it's a different uuid technically..., use ClimaParams = "0.10"
    import Thermodynamics as TD
    import Thermodynamics.Parameters as TDP
    # FT = Float64

    toml_dict = CP.create_toml_dict(FT; dict_type = "alias") # CP 0.7 and below, Thermodynamics 0.11 and above
    aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)


    namelist = NameList.default_namelist("SOCRATES_RF09_Obs_data"; write = false) # we don't need a directory
    param_set = create_parameter_set(namelist, toml_dict, FT) # this creates an override file in  a directory we don' need...

    microphys_params = TCP.microphysics_params(param_set)
end

S_qs = 6.5237565853107145e-19;
q_ice = 6.452298441776985e-9;
Dmax = Inf;
Dmin = 0.0;
r_ice_snow = 6.25e-5
N = 10.4128991889143
ρ = 1.0228946958166607

μ_ice = FT(0)


v_term = TC.my_terminal_velocity(microphys_params, TC.ice_type, TC.CMT.Chen2022Type(), ρ, q_ice; Dmin = Dmin, Dmax=Dmax, Nt=N)

v_diff_number_mean = TC.mean_velocity_difference(microphys_params, TC.ice_type, q_ice, ρ, μ_ice; Nt=N) / N^2

v_diff_mass_mean = TC.mass_weighted_velocity_difference(microphys_params, TC.ice_type, q_ice, ρ, μ_ice; Nt=N)


q_int = TC.q_int(microphys_params, TC.ice_type, q_ice, ρ, μ_ice, Nt=N, Dmin=r_ice_snow, Dmax=Inf)

@info "Calcs" v_term v_diff_number_mean v_diff_mass_mean q_int (q_int / q_ice)
