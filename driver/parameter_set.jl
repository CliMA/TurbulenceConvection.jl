import TurbulenceConvection
import CLIMAParameters
const CP = CLIMAParameters

struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CLIMAParameters.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
CLIMAParameters.Planet.cp_v(ps::EarthParameterSet) = ps.nt.cp_v
CLIMAParameters.Planet.R_d(ps::EarthParameterSet) = ps.nt.R_d
CLIMAParameters.Planet.R_v(ps::EarthParameterSet) = ps.nt.R_v
CLIMAParameters.Planet.molmass_ratio(ps::EarthParameterSet) = ps.nt.molmass_ratio
# microphysics parameters
CLIMAParameters.Atmos.Microphysics_0M.τ_precip(ps::EarthParameterSet) = ps.nt.τ_precip
CLIMAParameters.Atmos.Microphysics.τ_cond_evap(ps::EarthParameterSet) = ps.nt.τ_cond_evap
CLIMAParameters.Atmos.Microphysics.τ_sub_dep(ps::EarthParameterSet) = ps.nt.τ_sub_dep
CLIMAParameters.Atmos.Microphysics.τ_acnv_rai(ps::EarthParameterSet) = ps.nt.τ_acnv_rai
CLIMAParameters.Atmos.Microphysics.τ_acnv_sno(ps::EarthParameterSet) = ps.nt.τ_acnv_sno
CLIMAParameters.Atmos.Microphysics.q_liq_threshold(ps::EarthParameterSet) = ps.nt.q_liq_threshold
CLIMAParameters.Atmos.Microphysics.q_ice_threshold(ps::EarthParameterSet) = ps.nt.q_ice_threshold
CLIMAParameters.Atmos.Microphysics.microph_scaling(ps::EarthParameterSet) = ps.nt.microph_scaling
CLIMAParameters.Atmos.Microphysics.microph_scaling_dep_sub(ps::EarthParameterSet) = ps.nt.microph_scaling_dep_sub
CLIMAParameters.Atmos.Microphysics.microph_scaling_melt(ps::EarthParameterSet) = ps.nt.microph_scaling_melt
CLIMAParameters.Atmos.Microphysics.E_liq_rai(ps::EarthParameterSet) = ps.nt.E_liq_rai
CLIMAParameters.Atmos.Microphysics.E_liq_sno(ps::EarthParameterSet) = ps.nt.E_liq_sno
CLIMAParameters.Atmos.Microphysics.E_ice_rai(ps::EarthParameterSet) = ps.nt.E_ice_rai
CLIMAParameters.Atmos.Microphysics.E_ice_sno(ps::EarthParameterSet) = ps.nt.E_ice_sno
CLIMAParameters.Atmos.Microphysics.E_rai_sno(ps::EarthParameterSet) = ps.nt.E_rai_sno

# entrainment/detrainment parameters
CLIMAParameters.Atmos.EDMF.c_ε(ps::EarthParameterSet) = ps.nt.c_ε # factor multiplier for dry term in entrainment/detrainment
CLIMAParameters.Atmos.EDMF.α_b(ps::EarthParameterSet) = ps.nt.α_b # factor multiplier for pressure buoyancy terms (effective buoyancy is (1-α_b))
CLIMAParameters.Atmos.EDMF.α_a(ps::EarthParameterSet) = ps.nt.α_a # factor multiplier for pressure advection
CLIMAParameters.Atmos.EDMF.α_d(ps::EarthParameterSet) = ps.nt.α_d # factor multiplier for pressure drag
CLIMAParameters.Atmos.EDMF.H_up_min(ps::EarthParameterSet) = ps.nt.H_up_min # minimum updraft top to avoid zero division in pressure drag and turb-entr
CLIMAParameters.Atmos.EDMF.c_δ(ps::EarthParameterSet) = ps.nt.c_δ # factor multiplier for moist term in entrainment/detrainment
CLIMAParameters.Atmos.EDMF.β(ps::EarthParameterSet) = ps.nt.β # sorting power for ad-hoc moisture detrainment function
CLIMAParameters.Atmos.EDMF.χ(ps::EarthParameterSet) = ps.nt.χ # fraction of updraft air for buoyancy mixing in entrainment/detrainment (0≤χ≤1)
CLIMAParameters.Atmos.EDMF.c_γ(ps::EarthParameterSet) = ps.nt.c_γ # scaling factor for turbulent entrainment rate
CLIMAParameters.Atmos.EDMF.c_λ(ps::EarthParameterSet) = ps.nt.c_λ # scaling factor for TKE in entrainment scale calculations
CLIMAParameters.Atmos.EDMF.w_min(ps::EarthParameterSet) = ps.nt.w_min # minimum updraft velocity to avoid zero division in b/w²
CLIMAParameters.Atmos.EDMF.μ_0(ps::EarthParameterSet) = ps.nt.μ_0 # dimensional scale logistic function in the dry term in entrainment/detrainment
# mixing length parameters
CLIMAParameters.Atmos.EDMF.c_m(ps::EarthParameterSet) = ps.nt.c_m # tke diffusivity coefficient
CLIMAParameters.Atmos.EDMF.c_d(ps::EarthParameterSet) = ps.nt.c_d # tke dissipation coefficient
CLIMAParameters.Atmos.EDMF.c_b(ps::EarthParameterSet) = ps.nt.c_b # static stability coefficient
CLIMAParameters.Atmos.EDMF.κ_star²(ps::EarthParameterSet) = ps.nt.κ_star² # Ratio of TKE to squared friction velocity in surface layer
CLIMAParameters.Atmos.EDMF.Pr_n(ps::EarthParameterSet) = ps.nt.Pr_n # turbulent Prandtl number in neutral conditions
CLIMAParameters.Atmos.EDMF.ω_pr(ps::EarthParameterSet) = ps.nt.ω_pr # cospectral budget factor for turbulent Prandtl number
CLIMAParameters.Atmos.EDMF.Ri_c(ps::EarthParameterSet) = ps.nt.Ri_c # critical Richardson number
CLIMAParameters.Atmos.EDMF.smin_ub(ps::EarthParameterSet) = ps.nt.smin_ub #  lower limit for smin function
CLIMAParameters.Atmos.EDMF.smin_rm(ps::EarthParameterSet) = ps.nt.smin_rm #  upper ratio limit for smin function

#! format: off
function create_parameter_set(namelist)
    TC = TurbulenceConvection

    # this is needed to avoid slowdown when using large parameter vectors
    use_ran_features = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment") == "RF"
    use_nn = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment") == "NN"
    use_nn_nonlocal = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment") == "NN_nonlocal"
    use_fno = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment") == "FNO"
    use_linear = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment") == "Linear"

    entr_closure_kwargs = if use_ran_features
        (;
        c_rf_fix = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "rf_fix_ent_params"),
        c_rf_opt = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "rf_opt_ent_params"),
    )
    elseif use_nn
        (;
        c_nn_params = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_ent_params"),
        nn_arc = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_arc"),
    )
    elseif use_nn_nonlocal
        (;
        c_nn_params = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_ent_params"),
        nn_arc = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_arc"),
    )
    elseif use_fno
        (;
        w_fno = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fno_ent_width"),
        nm_fno = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fno_ent_n_modes"),
        c_fno = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fno_ent_params"),
    )
    elseif use_linear
        (;
        c_linear = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "linear_ent_params"),
    )
    else
        ()
    end

    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
        cp_d = 1004.0,
        cp_v = 1859.0,
        R_d = 287.1,
        R_v = 461.5,
        molmass_ratio = 461.5/287.1,
        τ_precip = TC.parse_namelist(namelist, "microphysics", "τ_precip"; default = 1000.0),
        τ_cond_evap = TC.parse_namelist(namelist, "microphysics", "τ_cond_evap"; default = 10.0),
        τ_sub_dep = TC.parse_namelist(namelist, "microphysics", "τ_sub_dep"; default = 10.0),
        τ_acnv_rai = TC.parse_namelist(namelist, "microphysics", "τ_acnv_rai"; default = 2500.0),
        τ_acnv_sno = TC.parse_namelist(namelist, "microphysics", "τ_acnv_sno"; default = 100.0),
        q_liq_threshold = TC.parse_namelist(namelist, "microphysics", "q_liq_threshold"; default = 0.5e-3),
        q_ice_threshold = TC.parse_namelist(namelist, "microphysics", "q_ice_threshold"; default = 1e-6),
        microph_scaling =         TC.parse_namelist(namelist, "microphysics", "microph_scaling"; default = 1.0),
        microph_scaling_dep_sub = TC.parse_namelist(namelist, "microphysics", "microph_scaling_dep_sub"; default = 1.0),
        microph_scaling_melt =    TC.parse_namelist(namelist, "microphysics", "microph_scaling_melt"; default = 1.0),
        E_liq_rai = TC.parse_namelist(namelist, "microphysics", "E_liq_rai"; default = 0.8),
        E_liq_sno = TC.parse_namelist(namelist, "microphysics", "E_liq_sno"; default = 0.1),
        E_ice_rai = TC.parse_namelist(namelist, "microphysics", "E_ice_rai"; default = 1.0),
        E_ice_sno = TC.parse_namelist(namelist, "microphysics", "E_ice_sno"; default = 0.1),
        E_rai_sno = TC.parse_namelist(namelist, "microphysics", "E_rai_sno"; default = 1.0),
        c_ε = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_factor"),
        c_div = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_massflux_div_factor"; default = 0.0),
        α_b = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_buoy_coeff1"),
        α_a = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_adv_coeff"),
        α_d = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_drag_coeff"),
        H_up_min = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_updraft_top"),
        ω_pr = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_scale"),
        c_δ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "detrainment_factor"),
        Π_norm = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pi_norm_consts"),
        β = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "sorting_power"),
        χ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "updraft_mixing_frac"),
        c_γ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "turbulent_entrainment_factor"),
        c_λ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_smin_tke_coeff"),
        w_min = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_upd_velocity"),
        μ_0 = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_scale"),
        γ_lim = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "area_limiter_scale"),
        β_lim = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "area_limiter_power"),
        c_m = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_ed_coeff"),
        c_d = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_diss_coeff"),
        c_b = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "static_stab_coeff"; default = 0.4), # this is here due to a value error in CliMAParmameters.jl
        κ_star² = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_surf_scale"),
        Pr_n = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_0"),
        Ri_c = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Ri_crit"),
        smin_ub = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_ub"),
        smin_rm = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_rm"),
        l_max = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "l_max"; default = 1.0e6),
        covar_lim = TC.parse_namelist(namelist, "thermodynamics", "diagnostic_covar_limiter"),
        prescribed_precip_frac_value = TC.parse_namelist(namelist, "microphysics", "prescribed_precip_frac_value"; default = 1.0),
        precip_fraction_limiter = TC.parse_namelist(namelist, "microphysics", "precip_fraction_limiter"; default = 0.3),
        entr_closure_kwargs...,
        ## Stochastic parameters
        c_gen_stoch = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "general_stochastic_ent_params"),
    )
    param_set = EarthParameterSet(nt)
    if !isbits(param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on
