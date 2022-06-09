import OrderedCollections
import JSON

all_cases = [
    "ARM_SGP",
    "Bomex",
    "DryBubble",
    "DYCOMS_RF01",
    "DYCOMS_RF02",
    "GABLS",
    #"GATE_III",
    "life_cycle_Tan2018",
    "Nieuwstadt",
    "Rico",
    "Soares",
    "TRMM_LBA",
    "LES_driven_SCM",
]

include(joinpath("..", "post_processing", "mse_tables.jl"))

percent_reduction_mse = Dict()

computed_mse = OrderedCollections.OrderedDict()
files_skipped = OrderedCollections.OrderedDict()
for case in all_cases
    filename = "computed_mse_$case.json"
    if !isfile(filename)
        @warn "File $filename skipped"
        files_skipped[case] = true
        continue
    end
    jsonfile = JSON.parsefile(filename; dicttype = OrderedCollections.OrderedDict)
    files_skipped[case] = false
    computed_mse[case] = jsonfile
end

println("#################################")
println("################################# MSE tables")
println("#################################")
println("#")

println("all_best_mse = OrderedCollections.OrderedDict()\n#")
for case in keys(computed_mse)
    println("all_best_mse[\"$case\"] = OrderedCollections.OrderedDict()")
    for var in keys(computed_mse[case])
        if computed_mse[case][var] == "NA"
            println("all_best_mse[\"$case\"][\"$var\"] = \"$(computed_mse[case][var])\"")
        else
            println("all_best_mse[\"$case\"][\"$var\"] = $(computed_mse[case][var])")
        end
    end
    println("#")
end

println("#################################")
println("#################################")
println("#################################")

# Cleanup
for case in all_cases
    rm("computed_mse_$case.json"; force = true)
end

#####
##### min percentage reduction of mse across cases
#####

println("-- DO NOT COPY --")

for case in keys(computed_mse)
    percent_reduction_mse[case] = 0
    for var in keys(computed_mse[case])
        if haskey(all_best_mse[case], var)
            all_best_mse[case][var] isa Real || continue # skip if "NA"
            computed_mse[case][var] isa Real || continue # skip if "NA"
            percent_reduction_mse[case] = min(
                percent_reduction_mse[case],
                (all_best_mse[case][var] - computed_mse[case][var]) / all_best_mse[case][var] * 100,
            )
        else
            percent_reduction_mse[case] = "NA"
        end
    end
end

for case in keys(percent_reduction_mse)
    @info "percent_reduction_mse[$case] = $(percent_reduction_mse[case])"
end
@info "min mse reduction (%) over all cases = $(min(values(percent_reduction_mse)...))"

if any(values(files_skipped))
    @show files_skipped
    error("MSE printing skipped files skipped")
end
