if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(path = ".")
end
import SnoopCompile

tinf = SnoopCompile.@snoopi_deep begin
    import TurbulenceConvection

    local_tc_dir = pkgdir(TurbulenceConvection)
    include(joinpath(local_tc_dir, "driver", "main.jl"))
    include(joinpath(local_tc_dir, "driver", "generate_namelist.jl"))
    import .NameList

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename, return_code = main(namelist)
end

import Profile

Profile.@profile begin
    import TurbulenceConvection

    local_tc_dir = pkgdir(TurbulenceConvection)
    include(joinpath(local_tc_dir, "driver", "main.jl"))
    include(joinpath(local_tc_dir, "driver", "generate_namelist.jl"))
    import .NameList


    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename, return_code = main(namelist)
end

import PyPlot # the GUI is dependent on PyPlot, must load it before the next line
mref, ax = SnoopCompile.pgdsgui(tinf);
# folder = "perf/pgd_output"
# mkpath(folder)
# PyPlot.savefig(joinpath(folder, "pgd.png"))
