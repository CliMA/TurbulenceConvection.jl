# THIS IS A WAY TO TEST THAT CREATED THERMODYNAMIC STATES ARE VALID AND STABLE..

FT = Float64

reload_TC = false
if reload_TC ||
   !isdefined(Main, :TurbulenceConvection) ||
   !isdefined(Main, :param_set) ||
   !isdefined(Main, :thermo_params)
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/"))
    # Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/"))
    tc_dir = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl")
    Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests/"))
    # include(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/src/TurbulenceConvection.jl")) 
    using ForwardDiff # for netcdf.io
    include(joinpath(tc_dir, "driver", "NetCDFIO.jl"))
    include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
    include(joinpath(tc_dir, "driver", "Cases.jl"))
    include(joinpath(tc_dir, "driver", "parameter_set.jl"))
    include(joinpath(tc_dir, "driver", "dycore.jl"))
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    namelist = NameList.default_namelist("SOCRATES_RF09_Obs_data"; write = false) # we don't need a directory
    namelist["root"] = thisdir
    namelist["output"]["output_root"] = thisdir
    namelist["meta"] = Dict("simname" => "test_mm_2015", "uuid" => "")

    param_set = create_parameter_set(namelist, toml_dict, FT) # this creates an override file in  a directory we don' need...
    thermo_params = TCP.thermodynamics_params(param_set)
end


ρ_c = 1.240569830292043;
area = 0.1;
ρarea = 0.1240569830292043;
w = 0.01;
ρaw = 0.0012380112658666187;
ρaq_tot = 0.0005582662190320437;
ρaθ_liq_ice = 34.495528744250095;
q_liq = 0.0;
q_ice = 0.0;
q_tot = 0.004500078958881517;
θ_liq_ice = 278.061967185914;
T = 277.23428535758245;
p_c = 98960.9332762119;

# θ_liq_ice = ρaθ_liq_ice / ρarea
# q_tot = ρaq_tot / ρarea
# q_liq = FT(0)
# q_ice = FT(0)

q = TD.PhasePartition(q_tot, q_liq, q_ice)


ts1 = TD.PhaseNonEquil_ρθq(thermo_params, ρ_c, θ_liq_ice, q)
println("ts1 = $ts1")
p_c2 = TD.air_pressure(thermo_params, ts)
ts2 = TC.thermo_state_pθq(param_set, p_c2, θ_liq_ice, q_tot, q_liq, q_ice)
println("ts2 = $ts2")

ts = TC.thermo_state_pθq(param_set, p_c, θ_liq_ice, q_tot, q_liq, q_ice)
println("ts = $ts")


ts = ts1
T = TD.air_temperature(thermo_params, ts)
println("T = $T;")

T1 = TD.air_temperature(thermo_params, ts1)
println("T1 = $T1;")

T2 = TD.air_temperature(thermo_params, ts2)
println("T2 = $T2;")




ts4 =
    TD.PhaseNonEquil{Float64}(35818.442110822296, 1.1163156887380072, TD.PhasePartition(0.004500078958881517, 0.0, 0.0));
q = TD.PhasePartition(0.004500078958881517, 0.0, 0.0);
T = 308.0380948417583;
p = 98960.93327621189;
ρ = 1.240569830292043;
w = 0.0;
z = 23.96;



ts5 = TD.PhaseNonEquil(12979.93376715246, 1.2355507018081924, TD.PhasePartition(0.00436904618444261, 0.0, 0.0));
q = TD.PhasePartition(0.00436904618444261, 0.0, 0.0);
T = 276.7751427474668;
p = 98406.88158043097;
ρ = 1.2357025148253538;
w = 0.0;
z = 69.56;
θ_liq_ice5 = TD.liquid_ice_pottemp(thermo_params, ts5)
println("θ_liq_ice5 = $θ_liq_ice5")
