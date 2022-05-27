import TurbulenceConvection
const TC = TurbulenceConvection
import Plots
import NCDatasets
import Logging
import CLIMAParameters
import ClimaCore
const CC = ClimaCore
const tc_dir = pkgdir(TurbulenceConvection)
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "driver", "Cases.jl"))
include(joinpath(tc_dir, "driver", "parameter_set.jl"))
include(joinpath(tc_dir, "driver", "dycore.jl"))
import .NameList
import .Cases
function export_ref_profile(case_name::String)
    TC = TurbulenceConvection
    namelist = NameList.default_namelist(case_name)
    param_set = create_parameter_set(namelist)
    namelist["meta"]["uuid"] = "01"

    FT = Float64
    grid = TC.Grid(FT(namelist["grid"]["dz"]), namelist["grid"]["nz"])
    case = Cases.get_case(namelist)
    surf_ref_state = Cases.surface_ref_state(case, param_set, namelist)
    ts_g = surf_ref_state

    aux_vars(FT) = (; ρ = FT(0), p = FT(0))
    cent = TC.FieldFromNamedTuple(TC.center_space(grid), aux_vars(FT))
    face = TC.FieldFromNamedTuple(TC.face_space(grid), aux_vars(FT))

    p_c = cent.p
    ρ_c = cent.ρ
    p_f = face.p
    ρ_f = face.ρ
    compute_ref_state!(p_c, ρ_c, p_f, ρ_f, grid, param_set; ts_g)

    zc = vec(grid.zc.z)
    zf = vec(grid.zf.z)

    ρ_c = vec(cent.ρ)
    p_c = vec(cent.p)
    ρ_f = vec(face.ρ)
    p_f = vec(face.p)

    p1 = Plots.plot(ρ_c, zc ./ 1000; label = "centers")
    Plots.plot!(ρ_f, zf ./ 1000; label = "faces")
    Plots.plot!(size = (1000, 400))
    Plots.plot!(margin = 5 * Plots.mm)
    Plots.xlabel!("ρ")
    Plots.ylabel!("z (km)")
    Plots.title!("ρ")

    p2 = Plots.plot(p_c ./ 1000, zc ./ 1000; label = "centers")
    Plots.plot!(p_f ./ 1000, zf ./ 1000; label = "faces")
    Plots.plot!(size = (1000, 400))
    Plots.plot!(margin = 5 * Plots.mm)
    Plots.xlabel!("p (kPa)")
    Plots.ylabel!("z (km)")
    Plots.title!("p (kPa)")
    Plots.savefig("$case_name.svg")

end

Logging.with_logger(Logging.NullLogger()) do # silence output
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
        "TRMM_LBA",
        "GATE_III",
    )
        export_ref_profile(case_name)
    end
end
