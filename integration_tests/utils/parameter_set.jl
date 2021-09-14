import TurbulenceConvection
import CLIMAParameters
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters.Planet

struct EarthParameterSet{NT} <: AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CLIMAParameters.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
CLIMAParameters.Planet.cp_v(ps::EarthParameterSet) = ps.nt.cp_v
CLIMAParameters.Planet.R_d(ps::EarthParameterSet) = ps.nt.R_d
CLIMAParameters.Planet.R_v(ps::EarthParameterSet) = ps.nt.R_v
CLIMAParameters.Planet.molmass_ratio(ps::EarthParameterSet) = ps.nt.molmass_ratio
CLIMAParameters.Planet.T_freeze(ps::EarthParameterSet) = ps.nt.T_freeze
# microphysics 0-moment parameters
CLIMAParameters.Atmos.Microphysics_0M.τ_precip(ps::EarthParameterSet) = ps.nt.τ_precip
# microphysics 1-moment parameters
CLIMAParameters.Atmos.Microphysics.τ_acnv(ps::EarthParameterSet) = ps.nt.τ_acnv
# entrainment/detrainment parameters
CLIMAParameters.Atmos.EDMF.c_ε(ps::EarthParameterSet) = ps.nt.c_ε # factor multiplyer for dry term in entrainment/detrainment
CLIMAParameters.Atmos.EDMF.α_b(ps::EarthParameterSet) = ps.nt.α_b # factor multiplyer for pressure buoyancy terms (effective buoyancy is (1-α_b))
CLIMAParameters.Atmos.EDMF.α_a(ps::EarthParameterSet) = ps.nt.α_a # factor multiplyer for pressure advection
CLIMAParameters.Atmos.EDMF.α_d(ps::EarthParameterSet) = ps.nt.α_d # factor multiplyer for pressure drag
CLIMAParameters.Atmos.EDMF.H_up_min(ps::EarthParameterSet) = ps.nt.H_up_min # minimum updarft top to avoid zero devision in pressure drag and turb-entr
CLIMAParameters.Atmos.EDMF.ω_pr(ps::EarthParameterSet) = ps.nt.ω_pr # experiment scale factor for turbulent Prandtl number
CLIMAParameters.Atmos.EDMF.c_δ(ps::EarthParameterSet) = ps.nt.c_δ # factor multiplyer for moist term in entrainment/detrainment
CLIMAParameters.Atmos.EDMF.β(ps::EarthParameterSet) = ps.nt.β # sorting power for ad-hoc moisture detrainment function
CLIMAParameters.Atmos.EDMF.χ(ps::EarthParameterSet) = ps.nt.χ # fraction of updraft air for buoyancy mixing in entrainment/detrainment (0≤χ≤1)
CLIMAParameters.Atmos.EDMF.c_t(ps::EarthParameterSet) = ps.nt.c_t # factor multiplyer for turbulent term in entrainment/detrainment
CLIMAParameters.Atmos.EDMF.c_λ(ps::EarthParameterSet) = ps.nt.c_λ # scaling factor for TKE in entrainment scale calculations
CLIMAParameters.Atmos.EDMF.w_min(ps::EarthParameterSet) = ps.nt.w_min # minimum updraft velocity to aviod zero division in b/w²
CLIMAParameters.Atmos.EDMF.μ_0(ps::EarthParameterSet) = ps.nt.μ_0 # dimentional scale logistic function in the dry term in entrainment/detrainment
# mixing length parameters
CLIMAParameters.Atmos.EDMF.c_m(ps::EarthParameterSet) = ps.nt.c_m # tke diffusivity coefficient
CLIMAParameters.Atmos.EDMF.c_d(ps::EarthParameterSet) = ps.nt.c_d # tke dissipation coefficient
CLIMAParameters.Atmos.EDMF.Pr_n(ps::EarthParameterSet) = ps.nt.Pr_n # tke dissipation coefficient
CLIMAParameters.Atmos.EDMF.smin_ub(ps::EarthParameterSet) = ps.nt.smin_ub #  lower limit for smin function
CLIMAParameters.Atmos.EDMF.smin_rm(ps::EarthParameterSet) = ps.nt.smin_rm #  upper ratio limit for smin function

#! format: off
function create_parameter_set(namelist)
    TC = TurbulenceConvection
    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
        cp_d = 1004.0,
        cp_v = 1859.0,
        R_d = 287.1,
        R_v = 461.5,
        T_freeze = 0.0, # liquid only for now until we add ice consistently in the model
        molmass_ratio = 461.5/287.1,
        τ_precip = TC.parse_namelist(namelist, "microphysics", "τ_precip"; default = 1000.0),
        τ_acnv = TC.parse_namelist(namelist, "microphysics", "τ_acnv"; default = 1000.0),
        c_ε = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_factor"),
        c_div = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_massflux_div_factor"; default = 0.0),
        α_b = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_buoy_coeff1"),
        α_a = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_adv_coeff"),
        α_d = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_drag_coeff"),
        H_up_min = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_updraft_top"),
        ω_pr = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_scale"),
        c_δ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "detrainment_factor"),
        β = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "sorting_power"),
        χ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "updraft_mixing_frac"),
        c_t = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "turbulent_entrainment_factor"),
        c_λ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_smin_tke_coeff"),
        w_min = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_upd_velocity"),
        μ_0 = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_scale"),
        μ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_sigma"),
        c_m = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_ed_coeff"),
        c_d = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_diss_coeff"),
        c_b = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "static_stab_coeff"; default = 0.4), # this is here due to an value error in CliMAParmaeters.jl
        Pr_n = TC.parse_namelist(namelist, "turbulence", "prandtl_number_0"), # this is here due to an value error in CliMAParmaeters.jl
        smin_ub = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_ub"),
        smin_rm = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_rm"),
        l_max = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "l_max"; default = 1.0e6),

        ## Stochastic parameters
        # lognormal model
        stoch_ε_lognormal_var = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stochastic", "entr_lognormal_var"; default = 0.0),
        stoch_δ_lognormal_var = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stochastic", "detr_lognormal_var"; default = 0.0),
        # sde model
        sde_ϵ_θ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stochastic", "sde_entr_theta"; default = 1.0),
        sde_ϵ_σ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stochastic", "sde_entr_std"; default = 0.0),
        sde_δ_θ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stochastic", "sde_detr_theta"; default = 1.0),
        sde_δ_σ = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stochastic", "sde_detr_std"; default = 0.0),
    )
    return EarthParameterSet(nt)
end
#! format: on
