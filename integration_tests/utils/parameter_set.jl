import CLIMAParameters
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters.Planet

struct EarthParameterSet{NT} <: AbstractEarthParameterSet
    nt::NT
end

CLIMAParameters.Planet.MSLP(ps::EarthParameterSet) = ps.nt.MSLP

function create_parameter_set(namelist)
    nt = (;
        MSLP = 100000.0, # or grab from, e.g., namelist[""][...]
    )
    return EarthParameterSet(nt)
end
