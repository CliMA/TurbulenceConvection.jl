import TurbulenceConvection
import CLIMAParameters
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters.Planet

struct EarthParameterSet{NT} <: AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CLIMAParameters.Atmos.Microphysics.τ_cond_evap(ps::EarthParameterSet) = ps.nt.τ_cond_evap
CLIMAParameters.Atmos.EDMF.c_ε(ps::EarthParameterSet) = ps.nt.c_ε

#! format: off
function create_parameter_set(namelist)
    TC = TurbulenceConvection
    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
        τ_cond_evap = TC.parse_namelist(namelist, "microphysics", "τ_cond_evap"; default = 10.0),
        c_ε = TC.parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_factor"),
    )
    return EarthParameterSet(nt)
end
#! format: on
