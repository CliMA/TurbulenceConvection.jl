# From: https://timholy.github.io/SnoopCompile.jl/stable/snoopr/
using SnoopCompileCore
invalidations = @snoopr begin
    import Pkg
    import TurbulenceConvection

    const tc_dir = pkgdir(TurbulenceConvection)
    include(joinpath(tc_dir, "driver", "main.jl"))
    include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
    import .NameList

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["meta"]["uuid"] = "01_invalidations"
    ds_tc_filename, return_code = main(namelist)
end;

import ReportMetrics
ReportMetrics.report_invalidations(;
    job_name = "invalidations",
    invalidations,
    process_filename = x -> last(split(x, "packages/")),
)
