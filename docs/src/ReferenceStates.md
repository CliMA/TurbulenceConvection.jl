# Reference States

```@example
import TurbulenceConvection
import Plots
import NCDatasets
import CLIMAParameters
tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "integration_tests", "utils", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "Cases.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "parameter_set.jl"))
using .NameList
import .Cases
function export_ref_profile(case_name::String)
    namelist = default_namelist(case_name)
    param_set = create_parameter_set(namelist)
    namelist["meta"]["uuid"] = "01"
    grid = TurbulenceConvection.Grid(namelist)
    ref_state = TurbulenceConvection.ReferenceState(grid, param_set)
    Stats = TurbulenceConvection.NetCDFIO_Stats(namelist, grid)
    case = Cases.CasesFactory(namelist, grid, ref_state)
    Cases.initialize_reference(case, grid, ref_state, Stats)
    NCDatasets.Dataset(joinpath(Stats.path_plus_file), "r") do ds
        zc = ds.group["profiles"]["zc"][:]
        zf = ds.group["profiles"]["zf"][:]
        ρc_0 = ds.group["reference"]["rho0_half"][:]
        pc_0 = ds.group["reference"]["p0_half"][:]
        αc_0 = ds.group["reference"]["alpha0_half"][:]
        ρf_0 = ds.group["reference"]["rho0"][:]
        pf_0 = ds.group["reference"]["p0"][:]
        αf_0 = ds.group["reference"]["alpha0"][:]

        p1 = Plots.plot(ρc_0, zc ./ 1000;label="centers")
        Plots.plot!(ρf_0, zf ./ 1000;label="faces")
        Plots.plot!(size=(1000,400))
        Plots.plot!(margin=5Plots.mm)
        Plots.xlabel!("ρ_0")
        Plots.ylabel!("z (km)")
        Plots.title!("ρ_0")

        p2 = Plots.plot(pc_0 ./ 1000, zc ./ 1000;label="centers")
        Plots.plot!(pf_0 ./ 1000, zf ./ 1000;label="faces")
        Plots.plot!(size=(1000,400))
        Plots.plot!(margin=5Plots.mm)
        Plots.xlabel!("p_0 (kPa)")
        Plots.ylabel!("z (km)")
        Plots.title!("p_0 (kPa)")

        p3 = Plots.plot(αc_0, zc ./ 1000;label="centers")
        Plots.plot!(αf_0, zf ./ 1000;label="faces")
        Plots.plot!(size=(1000,400))
        Plots.plot!(margin=5Plots.mm)
        Plots.xlabel!("α_0")
        Plots.ylabel!("z (km)")
        Plots.title!("α_0")
        Plots.plot(p1,p2,p3; layout=(1,3))
        Plots.savefig("$case_name.svg")
    end
end
for case_name in (
    "Bomex",
    "life_cycle_Tan2018",
    "Soares",
    "Rico",
    "ARM_SGP",
    "DYCOMS_RF01",
    "GABLS",
    "SP",
    "DryBubble",
)
    export_ref_profile(case_name)
end;

# Note: temperatures in this case become extremely low.
CLIMAParameters.Planet.T_freeze(::EarthParameterSet) = 100.0
export_ref_profile("TRMM_LBA")
export_ref_profile("GATE_III")
```

## Bomex
![](Bomex.svg)

## life\_cycle\_Tan2018
![](life_cycle_Tan2018.svg)

## Soares
![](Soares.svg)

## Rico
![](Rico.svg)

## TRMM\_LBA
![](TRMM_LBA.svg)

## ARM\_SGP
![](ARM_SGP.svg)

## GATE\_III
![](GATE_III.svg)

## DYCOMS\_RF01
![](DYCOMS_RF01.svg)

## GABLS
![](GABLS.svg)

## SP
![](SP.svg)

## DryBubble
![](DryBubble.svg)


