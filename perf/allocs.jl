if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Coverage
import Plots

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

all_cases = [
    # "ARM_SGP",
    # "Bomex",
    # "DryBubble",
    # "DYCOMS_RF01",
    # "GABLS",
    # "GATE_III",
    # "life_cycle_Tan2018",
    # "Nieuwstadt",
    "Rico",
    # "Soares",
    # "SP",
    "TRMM_LBA",
    "LES_driven_SCM",
]

filter!(x -> x ≠ "GATE_III", all_cases) # no mse tables for GATE_III
filter!(x -> x ≠ "SP", all_cases) # not currently running SP
allocs = Dict()
for case in all_cases
    ENV["ALLOCATION_CASE_NAME"] = case
    run(`julia --project=test/ --track-allocation=user perf/alloc_per_case.jl`)
    allocs[case] = Coverage.analyze_malloc(".")

    # Clean up files
    all_files = [joinpath(root, f) for (root, dirs, files) in Base.Filesystem.walkdir(".") for f in files]
    all_mem_files = filter(x -> endswith(x, ".mem"), all_files)
    for f in all_mem_files
        rm(f)
    end
end

@info "Post-processing allocations"

function plot_allocs(case_name, allocs_per_case, n_unique_bytes)
    p = Plots.plot()
    @info "Allocations for $case_name"

    filename_only(fn) = first(split(fn, ".jl")) * ".jl"
    function compile_tc(fn, linenumber)
        c1 = endswith(filename_only(fn), "TurbulenceConvection.jl")
        c2 = linenumber == 1
        return c1 && c2
    end

    filter!(x -> x.bytes ≠ 0, allocs_per_case)
    filter!(x -> !compile_tc(x.filename, x.linenumber), allocs_per_case)

    for alloc in allocs_per_case
        println(alloc)
    end
    println("Number of allocating sites: $(length(allocs_per_case))")
    case_bytes = getproperty.(allocs_per_case, :bytes)[end:-1:1]
    case_filename = getproperty.(allocs_per_case, :filename)[end:-1:1]
    case_linenumber = getproperty.(allocs_per_case, :linenumber)[end:-1:1]
    all_bytes = Int[]
    filenames = String[]
    linenumbers = Int[]
    loc_ids = String[]
    for (bytes, filename, linenumber) in zip(case_bytes, case_filename, case_linenumber)
        compile_tc(filename, linenumber) && continue # Skip loading module
        loc_id = "$(filename_only(filename))" * "$linenumber"
        if !(bytes in all_bytes) && !(loc_id in loc_ids)
            push!(all_bytes, bytes)
            push!(filenames, filename)
            push!(linenumbers, linenumber)
            push!(loc_ids, loc_id)
            if length(all_bytes) ≥ n_unique_bytes
                break
            end
        end
    end

    all_bytes = all_bytes ./ 10^3
    max_bytes = maximum(all_bytes)
    @info "$case_name: $all_bytes"
    xtick_name(filename, linenumber) = "$filename, line number: $linenumber"
    markershape = (:square, :hexagon, :circle, :star, :utriangle, :dtriangle)
    for (bytes, filename, linenumber) in zip(all_bytes, filenames, linenumbers)
        Plots.plot!(
            [0],
            [bytes];
            seriestype = :scatter,
            label = xtick_name(filename_only(filename), linenumber),
            markershape = markershape[1],
            markersize = 1 + bytes / max_bytes * 10,
        )
        markershape = (markershape[end], markershape[1:(end - 1)]...)
    end
    p1 = Plots.plot!(ylabel = "Allocations (KB)", title = case_name)
    p2 = Plots.plot(
        1:length(allocs_per_case),
        getproperty.(allocs_per_case, :bytes)[end:-1:1] ./ 1000;
        xlabel = "i-th allocating line (sorted)",
        ylabel = "Allocations (KB)",
        markershape = :circle,
    )
    Plots.plot(p1, p2, layout = Plots.grid(2, 1))
    Plots.savefig(joinpath(folder, "allocations_$case_name.png"))
end

folder = "perf/allocations_output"
mkpath(folder)

@info "Allocated bytes for single tendency per case:"
for case in all_cases
    plot_allocs(case, allocs[case], 10)
end
