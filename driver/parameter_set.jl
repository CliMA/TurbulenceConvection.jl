import TurbulenceConvection
import CLIMAParameters
const CP = CLIMAParameters

struct ThermodynamicsParameters{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end
struct MicrophysicsParameters{NT} <: CP.AbstractEarthParameterSet
    nt::NT
end
CLIMAParameters.Planet.MSLP(ps::ThermodynamicsParameters) = ps.nt.MSLP

struct EarthParameterSet{FT, TP, MP, NT} <: CP.AbstractEarthParameterSet
    thermo_params::TP
    microphys_params::MP
    nt::NT
end

EarthParameterSet{FT}(args...) where {FT} = EarthParameterSet{FT, typeof.(args)...}(args...)

Base.eltype(::EarthParameterSet{FT}) where {FT} = FT

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
# microphysics parameters
CP.Atmos.Microphysics_0M.τ_precip(ps::MicrophysicsParameters) = ps.nt.τ_precip
CP.Atmos.Microphysics.τ_cond_evap(ps::MicrophysicsParameters) = ps.nt.τ_cond_evap
CP.Atmos.Microphysics.τ_sub_dep(ps::MicrophysicsParameters) = ps.nt.τ_sub_dep
CP.Atmos.Microphysics.τ_acnv_rai(ps::MicrophysicsParameters) = ps.nt.τ_acnv_rai
CP.Atmos.Microphysics.τ_acnv_sno(ps::MicrophysicsParameters) = ps.nt.τ_acnv_sno
CP.Atmos.Microphysics.q_liq_threshold(ps::MicrophysicsParameters) = ps.nt.q_liq_threshold
CP.Atmos.Microphysics.q_ice_threshold(ps::MicrophysicsParameters) = ps.nt.q_ice_threshold
CP.Atmos.Microphysics.microph_scaling(ps::MicrophysicsParameters) = ps.nt.microph_scaling
CP.Atmos.Microphysics.microph_scaling_dep_sub(ps::MicrophysicsParameters) = ps.nt.microph_scaling_dep_sub
CP.Atmos.Microphysics.microph_scaling_melt(ps::MicrophysicsParameters) = ps.nt.microph_scaling_melt
CP.Atmos.Microphysics.E_liq_rai(ps::MicrophysicsParameters) = ps.nt.E_liq_rai
CP.Atmos.Microphysics.E_liq_sno(ps::MicrophysicsParameters) = ps.nt.E_liq_sno
CP.Atmos.Microphysics.E_ice_rai(ps::MicrophysicsParameters) = ps.nt.E_ice_rai
CP.Atmos.Microphysics.E_ice_sno(ps::MicrophysicsParameters) = ps.nt.E_ice_sno
CP.Atmos.Microphysics.E_rai_sno(ps::MicrophysicsParameters) = ps.nt.E_rai_sno

#! format: off
function create_parameter_set(::Type{FT}, namelist) where {FT}

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
    param_set = EarthParameterSet{FT}(
        ThermodynamicsParameters((;MSLP_kwarg...)),
        MicrophysicsParameters((;nt...)),
        nt
    )
    if !isbits(param_set)
        @warn "The parameter set SHOULD be isbits in order to be stack-allocated."
    end
    return param_set
end
#! format: on
