
"""
This is to allow precompilation on HPC so that jobs don't require huge amounts of memory during startup.
"""

# using Pkg: Pkg
# # Import things we need
thisdir = @__DIR__
tc_dir = joinpath(thisdir, "..")


# Check if we are on an HPC system and set namelist path accordingly

is_HPC = occursin(r"(login|node|compute|batch|cluster)", lowercase(gethostname())) # We could also condition this on an ENV var for more control...

@warn("Precompilation workload running. is_HPC = $is_HPC")

if is_HPC # (if we can)
    # Maybe there's some simple methods we can infer here, but we don't have an easy way of say, generating a param_set.

    # include(joinpath(tc_dir, "driver", "parameter_set.jl")) # This imports TC so can't use as is, need to extract out the relevant bits.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
    using PrecompileTools
    import CLIMAParameters as CP
    import Thermodynamics as TD
    import CloudMicrophysics as CM
    import SurfaceFluxes as SF
    import SurfaceFluxes.UniversalFunctions as UF

    # RELATIVE IMPORT: Access the submodule defined in src/Parameters.jl
    import .Parameters as TCP 

    # Helper for the NamedTuple logic
    namedtuple_fromdict(d::Dict) = (; (Symbol(k) => namedtuple_fromdict(v) for (k,v) in d)...)
    namedtuple_fromdict(x) = x

    PrecompileTools.@setup_workload begin
        # include("../integration_tests/overwrite_namelist.jl")
        # overwrite_namelist!(namelist, parsed_args)

        # --- STEP 1: SETUP MOCK DATA ---
        FT = Float64
        
        # Create a base TOML dictionary (Standard Physics Constants)
        # Note: If this fails without a file, use: CP.create_toml_dict(FT; dict_type="alias") 
        # but usually the default works in newer versions.

        # -=- Handle overrides (from parameters_set.jl) -=- #
        override_file = joinpath(thisdir, "override_dict.toml")
        # Read in data from namelist to overwrite parameters
        τ_precip = FT(1000)
        τ_cond_evap = FT(10)
        τ_sub_dep = FT(10)
        τ_acnv_rai = FT(2500)
        τ_acnv_sno = FT(100)
        q_liq_threshold = FT(0.5e-3)
        q_ice_threshold = FT(1e-6)
        microph_scaling_acnv = FT(1.0)
        microph_scaling_accr = FT(1.0)
        microph_scaling_evap = FT(1.0)
        microph_scaling_dep_sub = FT(1.0)
        microph_scaling_melt = FT(1.0)
        E_liq_rai = FT(0.8)
        E_liq_sno = FT(0.1)
        E_ice_rai = FT(1.0)
        E_ice_sno = FT(0.1)
        E_rai_sno = FT(1.0)
        A_acnv_KK2000 = FT(7.42e13)
        a_acnv_KK2000 = FT(2.47)
        b_acnv_KK2000 = FT(-1.79)
        c_acnv_KK2000 = FT(-1.47)
        pow_icenuc = FT(1e7) # I think this has to be here to overwrite toml..., and to place it so we can overwrite_namelist it .. picked high value initially to keep ramp but maybe should keep original default... https://github.com/CliMA/CLIMAParameters.jl/blob/2f298661d527c1b9a11188ee2abc92ba4ba0c5ec/src/parameters.toml#L558
        r_ice_snow = FT(62.5e-6) # allow changing this in the namelist/param_set/microphys_params instead of user_params (probably not ideal if you used a 2-moment or something where this mattered more but ...)
        #
        χm_ice = FT(1.0)
        #
        ν_sno = FT(0.63)
        μ_sno = FT(4.36e9)
        me_sno = FT(2.0)
        χm_sno = FT(1.0)
        #
        χv_rai = FT(1.0)
        χv_sno = FT(1.0)
        open(override_file, "w") do io # create the override file
            # Override the default files in the toml file with the values from the namelist -- [[ i think it's plain language name name, alias = code_name ]]
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
            println(io, "[microph_scaling_evap]")
            println(io, "alias = \"microph_scaling_evap\"")
            println(io, "value = " * string(microph_scaling_evap))
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
        toml_dict = CP.create_toml_dict(FT; override_file = override_file, dict_type="alias")   
        isfile(override_file) && rm(override_file; force=true)

        # --- STEP 2: CONSTRUCT SUB-PARAMETERS ---
        # A. Thermodynamics
        aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
        # FIX: Renamed 'pairs' to 'param_pairs' to avoid conflict with Base.pairs
        param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
        thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; param_pairs...)
        TP = typeof(thermo_params)

        # B. CloudMicrophysics
        aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
        param_pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
        microphys_params = CM.Parameters.CloudMicrophysicsParameters{FT, TP}(; param_pairs..., thermo_params)
        MP = typeof(microphys_params)

        # C. SurfaceFluxes
        # Universal Functions (Businger)
        aliases = ["Pr_0_Businger", "a_m_Businger", "a_h_Businger", "ζ_a_Businger", "γ_Businger"]
        uf_pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
        uf_pairs = (; uf_pairs...) 
        # Map the TOML names to the struct field names manually
        pairs_uf_mapped = (; Pr_0 = uf_pairs.Pr_0_Businger, a_m = uf_pairs.a_m_Businger, a_h = uf_pairs.a_h_Businger, ζ_a = uf_pairs.ζ_a_Businger, γ = uf_pairs.γ_Businger)
        ufp = UF.BusingerParams{FT}(; pairs_uf_mapped...)
        UFP = typeof(ufp)
        
        # Surface Flux Parameters
        sf_pairs = CP.get_parameter_values!(toml_dict, ["von_karman_const"], "SurfaceFluxesParameters")
        surf_flux_params = SF.Parameters.SurfaceFluxesParameters{FT, UFP, TP}(; sf_pairs..., ufp, thermo_params)
        SFP = typeof(surf_flux_params)

        # --- STEP 3: CONSTRUCT USER PARAMETERS ---
        user_params_dict = Dict(
            "particle_min_radius" => 0.2e-6,
            "r_liq_rain" => 30e-6,
            "χm_liq" => 1.0,
            "r_ice_acnv_scaling_factor" => 1.0,
            "r_ice_snow_threshold_scaling_factor" => 1.6, 
            "τ_acnv_sno_threshold" => 100.0,
            "massflux_N_i_boost_factor" => 2.0,
            "sedimentation_N_i_boost_factor" => 0.2,
            "apply_massflux_N_i_boost" => false,
            "apply_sedimentation_N_i_boost" => false,
            "massflux_N_i_boost_max_ratio" => 0.7,
            "massflux_N_i_boost_progress_fraction" => 0.5,
            "q_min" => 0.0
        )
        user_params = namedtuple_fromdict(user_params_dict)

        # --- STEP 4: ASSEMBLE THE FINAL STRUCT ---
        aliases = [
            "microph_scaling_dep_sub", "microph_scaling_melt", "microph_scaling_evap",
            "microph_scaling_acnv", "microph_scaling_accr", "Omega", "planet_radius"
        ]
        tc_pairs = CP.get_parameter_values!(toml_dict, aliases, "TurbulenceConvection")
        
        # Construct the REAL struct from src/Parameters.jl
        param_set = TCP.TurbulenceConvectionParameters{FT, MP, SFP, typeof(user_params)}(; 
            tc_pairs..., 
            microphys_params, 
            surf_flux_params, 
            user_params
        )

        microphysics_params = TCP.microphysics_params(param_set)
        prs = microphysics_params
        thermo_params = TCP.thermodynamics_params(param_set)















        PrecompileTools.@compile_workload begin
            # ======================================================================
            # 1. SETUP LOCAL SCALARS & STATE (Copied from your snippets)
            # ======================================================================

            local FT, area, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q_liq, q_ice, δ_eq, δi_eq, dδdt_no_S, dδdt, Δt, dqvdt, dTdt, Γ_l, Γ_i, qi_tendency_sub_dep, N, ρ, ql, N_INP, q_s, μ, Δt_factor, p_c, ρ_c, _n0, λ, k, regime, milestone_t, milestone, ts, q, q_eq, S_ql, S_qi, new_regime_type, limiter
            FT = Float64

            # Physics variables
            qi_tendency_sub_dep = FT(1e-9); qi = FT(1e-6); N = FT(200.0); ρ = FT(1.0); ql = FT(0); 
            N_INP = FT(300.0); q_s = FT(1e-5); μ = FT(0); Δt_factor = FT(1.0); p_c = FT(1e5); ρ_c = FT(1.0)
            
            # State variables
            S_ql = FT(1.57e-6); S_qi = FT(8.83e-8); area = FT(0.0499); p = FT(74651.0); T = FT(178.1)
            w = FT(0.0187); τ_liq = FT(6.65e-9); τ_ice = FT(0.148); q_vap = FT(5.09e-10)
            Δt = FT(0.00769)
            dqvdt = FT(1.28e-5); dTdt = FT(0.0138)
            Γ_l = FT(1.7); Γ_i = FT(1.6)
            
            # Construct State Objects
            q = TD.PhasePartition(FT(0.00495), ql, qi)
            q_eq = TD.PhasePartition(FT(0.00020), FT(5.76e-9), FT(2.96e-8))
            ts = TD.PhaseNonEquil(FT(-68295.0), ρ, q)
            
            # Derived variables
            q_liq = q.liq; q_ice = q.ice
            δ_0 = q_vap - q_eq.liq
            δ_0i = q_vap - q_eq.ice
            δ_eq = 0.0 # Placeholder
            δi_eq = 0.0 # Placeholder
            dδdt_no_S = dqvdt
            dδdt = dδdt_no_S - (S_ql*Γ_l + S_qi*Γ_i) 
            
            # ======================================================================
            # 2. ICE NUCLEATION / ADJUSTMENT
            # ======================================================================
            adjust_ice_N(param_set, N, N_INP, q.ice; 
                ρ=ρ, S_i=FT(0.05), q_l=ql, q_s=q_s, monodisperse=false, 
                ice_type=ice_type, decrease_N_if_subsaturated=true, N_INP_top=FT(150.)
            )

            # ======================================================================
            # 3. SATURATION REGIMES & MILESTONES
            # ======================================================================
            # A. Get Regime
            below_freezing = (T < TCP.T_freeze(param_set))
            regime = get_saturation_regime(δ_0, δ_0i, q.liq, q.ice, below_freezing)

            # B. Calculate Eq Point
            get_δ_eq_point(q_eq, τ_liq, τ_ice; dδdt_no_S = dδdt_no_S, Γ_l = Γ_l, Γ_i = Γ_i)

            # C. Calculate Milestones
            milestone_t, milestone, _, _, _, _ = calculate_next_standard_milestone_time(
                q_eq, q_liq, q_ice, δ_0, δ_0i, δ_eq, δi_eq, S_ql, S_qi; 
                dδdt=dδdt, dδdt_is_full_tendency=true, Γ_l=Γ_l, Γ_i=Γ_i, 
                at_δ_eq_point=false, allow_δ_eq_point=true
            )

            # D. Time Step Logic
            limiter = StandardSupersaturationMoistureSourcesLimiter()
            step(
                regime, limiter, Δt, q_liq, q_ice, δ_0, δ_0i, δ_eq, δi_eq, q_eq, S_ql, S_qi, 
                ((milestone_t < Δt) ? milestone : NotAtSupersaturationMilestone); 
                dδdt_no_S=dδdt_no_S, Γ_l=Γ_l, Γ_i=Γ_i
            )

            # E. Regime Transitions
            new_regime_type = get_new_saturation_regime_type_from_milestone(milestone, regime, δ_0, δ_0i)
            add_regime_parameters(new_regime_type, q_liq, q_ice, below_freezing)
            
            resolve_S_S_addit(S_ql, S_qi, Δt, FT(0.0), FT(0.0), milestone_t, Δt + milestone_t)

            # ======================================================================
            # 4. MORRISON-MILBRANDT PHYSICS KERNELS
            # ======================================================================
            # A. Exponential Part
            morrison_milbrandt_2015_style_exponential_part_only(
                regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, q, q_eq, Δt, ts; 
                use_fix=false, return_mixing_ratio=false, dqvdt=dqvdt, dTdt=dTdt, 
                fallback_to_standard_supersaturation_limiter=false, time_tolerance=FT(1e-6)
            )

            # B. Standard Fallback
            do_standard_fallback(
                milestone_t, milestone, FT(1e-6), S_ql, S_qi, q_liq, q_ice, δ_eq, δi_eq, 
                dδdt_no_S, Γ_l, Γ_i, regime, param_set, area, ρ, p, T, w, τ_liq, τ_ice, 
                δ_0, δ_0i, q, q_eq, Δt, ts; 
                use_fix=false, return_mixing_ratio=true, depth=0, dqvdt=dqvdt, dTdt=dTdt, 
                fallback_to_standard_supersaturation_limiter=false
            )

            # ======================================================================
            # 5. AUXILIARY / EPA FUNCTIONS
            # ======================================================================
            # Get params helper
            get_params_and_go_to_mixing_ratio_exponential_part_only(
                param_set, area, ρ, p, T, w, τ_liq, τ_ice, δ_0, δ_0i, dqvdt, dTdt, q, q_eq, Δt, ts; 
                use_fix=false
            )
            
            # Helper variables for EPA calls
            # (Assuming standard constants for A_c calculation are mocked or calculated)
            q_sl = FT(0.001); q_si = FT(0.0008); L_i = FT(2.8e6); c_p = FT(1004.0); e_sl = FT(600.0); dqsl_dT = FT(1e-4)
            
            # A_c Calc
            A_c_res = A_c_func_with_and_without_WBF(τ_ice, Γ_l, q_sl, q_si, 9.81, w, c_p, e_sl, L_i, dqsl_dT, dqvdt, dTdt, p, ρ)
            A_c = A_c_res.A_c
            
            # Tau func
            τ_EPA = τ_func_EPA(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)
            
            # Hit values
            t_δ_hit_value(FT(0), δ_0, A_c, τ_EPA)
            get_t_out_of_q_liq_EPA(δ_0, A_c, τ_EPA, τ_liq, q_liq, Γ_l)
            get_t_out_of_q_ice_EPA(δ_0, A_c, τ_EPA, τ_ice, q_ice, Γ_i, q_sl, q_si)

            # ======================================================================
            # 6. INTEGRALS
            # ======================================================================
            aiu, bi, ciu = FT(0.1), FT(0.01), FT(0.001)
            _n0, λ, k = FT(1e6), FT(1e3), FT(2.0)
            int__v_Dk_n__dD(aiu, bi, ciu, _n0, λ, k; Dmin=FT(0.0), Dmax=FT(Inf), μ=FT(0.0))

            int_nav_dr(
                param_set, ice_type, Chen2022Vel, FT(1e-6), FT(1.0), FT(1e-5); 
                Nt=FT(200.0), Dmin=FT(0.0), Dmax=FT(Inf), D_transition=FT((0.625e-3)/2)
            )

            # ======================================================================
            # 7. SEDIMENTATION & TERMINAL VELOCITY
            # ======================================================================
            velo_scheme = Chen2022Vel
            
            # Effective q
            q_eff = q_effective_nan_N_safe(param_set, liq_type, FT(1e-6), FT(200.0); monodisperse=true)
            
            # Terminal Velocity
            calculate_sedimentation_velocity(
                param_set, q_eff, ρ, liq_type, FT(200.0); velo_scheme=velo_scheme, Dmax=FT(Inf)
            )
            
            my_terminal_velocity(
                param_set, rain_type, velo_scheme, ρ, q_eff; Dmax=FT(Inf), Nt=FT(2.5e8)
            )
                
        end
    end    
end
