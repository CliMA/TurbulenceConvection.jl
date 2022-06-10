import TurbulenceConvection
import CLIMAParameters
const CP = CLIMAParameters

struct ThermodynamicsParameters{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end
CLIMAParameters.Planet.MSLP(ps::ThermodynamicsParameters) = ps.nt.MSLP

struct EarthParameterSet{TP, NT} <: CP.AbstractEarthParameterSet
    thermo_params::TP
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

#! format: off
function create_parameter_set(namelist)
    TC = TurbulenceConvection
    MSLP_kwarg = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
        )
    nt = (;
        MSLP_kwarg...,
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
    )
    param_set = EarthParameterSet(ThermodynamicsParameters((;MSLP_kwarg...)),nt)
    if !isbits(param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on
