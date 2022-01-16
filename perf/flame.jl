if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import SnoopCompileCore
import TurbulenceConvection
tc_dir_glob = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir_glob, "perf", "common.jl"))

case_name = "Bomex"
println("Running $case_name...")
sim = init_sim(case_name)
sim.skip_io || open_files(sim.Stats) # #removeVarsHack
(prob, alg, kwargs) = solve_args(sim)

tinf = SnoopCompileCore.@snoopi_deep begin
    sol = ODE.solve(prob, alg; kwargs...)
    # ds_tc_filename, return_code = main(namelist)
end

sim.skip_io || close_files(sim.Stats) # #removeVarsHack

import ProfileView
import SnoopCompile # need SnoopCompile to iterate over InferenceTimingNode's
import FlameGraphs
fg = FlameGraphs.flamegraph(tinf)
ProfileView.view(fg) # looks good, even without initial compiled run

# It would have been nice to auto-generate these flame graphs
# as a part of CI, but they're really large and slow to load / navigate.
# ProfileView works much better.
# import ProfileSVG
# folder = "perf/flame_output"
# mkpath(folder)
# ProfileSVG.save(joinpath(folder, "flame.svg"), fg; maxframes = 40000, maxdepth = 100)
