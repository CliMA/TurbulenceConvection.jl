import TurbulenceConvection
import CLIMAParameters
const CP = CLIMAParameters

struct EarthParameterSet{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
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
CLIMAParameters.Atmos.EDMF.α_b(ps::EarthParameterSet) = ps.nt.α_b # factor multiplier for pressure buoyancy terms (effective buoyancy is (1-α_b))
CLIMAParameters.Atmos.EDMF.α_a(ps::EarthParameterSet) = ps.nt.α_a # factor multiplier for pressure advection
CLIMAParameters.Atmos.EDMF.α_d(ps::EarthParameterSet) = ps.nt.α_d # factor multiplier for pressure drag
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

    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
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
        α_b = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_buoy_coeff1"),
        α_a = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_adv_coeff"),
        α_d = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_drag_coeff"),
        ω_pr = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_scale"),
        c_m = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_ed_coeff"),
        c_d = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_diss_coeff"),
        c_b = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "static_stab_coeff"; default = 0.4), # this is here due to a value error in CliMAParmameters.jl
        κ_star² = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_surf_scale"),
        Pr_n = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_0"),
        Ri_c = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Ri_crit"),
        smin_ub = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_ub"),
        smin_rm = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_rm"),

        # Can be removed before CP overhaul:
        l_max = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "l_max"; default = 1.0e6),
    )
    param_set = EarthParameterSet(nt)
    if !isbits(param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on
