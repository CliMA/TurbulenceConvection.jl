include(joinpath(@__DIR__, "common.jl"))
import SnoopCompileCore

case_name = "Bomex"
println("Running $case_name...")
sim = init_sim(case_name)
sim.skip_io || open_files(sim) # #removeVarsHack
(prob, alg, kwargs) = solve_args(sim)

tinf = SnoopCompileCore.@snoopi_deep begin
    sol = ODE.solve(prob, alg; kwargs...)
    # integrator, ds_tc_filenames, return_code = main(namelist)
end

sim.skip_io || close_files(sim) # #removeVarsHack

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
