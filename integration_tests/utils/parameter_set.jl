import TurbulenceConvection
import CLIMAParameters
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters.Planet

struct EarthParameterSet{NT} <: AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CLIMAParameters.Atmos.Microphysics.τ_cond_evap(ps::EarthParameterSet) = ps.nt.τ_cond_evap
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

#! format: off
function create_parameter_set(namelist)
    TC = TurbulenceConvection
    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
        τ_cond_evap = TC.parse_namelist(namelist, "microphysics", "τ_cond_evap"; default = 10.0),
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
    )
    return EarthParameterSet(nt)
end
#! format: on
