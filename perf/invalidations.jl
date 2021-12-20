# From: https://timholy.github.io/SnoopCompile.jl/stable/snoopr/
using SnoopCompileCore
invalidations = @snoopr begin
    if !("." in LOAD_PATH) # for easier local testing
        push!(LOAD_PATH, ".")
    end
    import TurbulenceConvection

    const tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
    include(joinpath(tc_dir, "driver", "main.jl"))
    include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
    import .NameList

    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["meta"]["uuid"] = "01"
    ds_tc_filename, return_code = main(namelist)
end;
import SnoopCompile
trees = SnoopCompile.invalidation_trees(invalidations);
@show length(SnoopCompile.uinvalidated(invalidations))

methinvs = trees[end]

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"
import Plots
n_invalidations = map(trees) do methinvs
    SnoopCompile.countchildren(methinvs)
end
p1 = Plots.plot(
    1:length(trees),
    n_invalidations;
    markershape = :circle,
    xlabel = "i-th method invalidation",
    label = "Number of children per method invalidations",
)
p2 = Plots.plot(
    1:length(trees),
    n_invalidations ./ sum(n_invalidations) .* 100;
    markershape = :circle,
    xlabel = "i-th method invalidation",
    label = "Percentage of invalidations",
)
Plots.plot(p1, p2; layout = (2, 1))
folder = "perf/invalidations_output"
mkpath(folder)
Plots.savefig(joinpath(folder, "invalidations.png"))

# More detail:
# root = methinvs.backedges[end]
# show(root; maxdepth=10)
# SnoopCompile.ascend(root)
