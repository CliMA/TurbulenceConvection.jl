import TurbulenceConvection
import CLIMAParameters
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters.Planet
using CLIMAParameters.Planet: grav

struct EarthParameterSet{NT} <: AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP
CLIMAParameters.Planet.R_v(ps::EarthParameterSet) = ps.nt.R_v
CLIMAParameters.Planet.R_d(ps::EarthParameterSet) = ps.nt.R_d
CLIMAParameters.Planet.cp_d(ps::EarthParameterSet) = ps.nt.cp_d
CLIMAParameters.Planet.cp_v(ps::EarthParameterSet) = ps.nt.cp_v
CLIMAParameters.Planet.grav(ps::EarthParameterSet) = ps.nt.grav
CLIMAParameters.Atmos.Microphysics.τ_cond_evap(ps::EarthParameterSet) = ps.nt.τ_cond_evap

#! format: off
function create_parameter_set(namelist)
    TC = TurbulenceConvection
    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
        grav = 9.80665,
        R_d = 287.1,
        R_v = 461.5,
        cp_d = 1004.0,
        cp_v = 1859.0,
        τ_cond_evap = TC.parse_namelist(namelist, "microphysics", "τ_cond_evap"; default = 10.0),
    )
    return EarthParameterSet(nt)
end
#! format: on
