using OrderedCollections
import JSON

all_cases = [
    "ARM_SGP",
    "Bomex",
    "DryBubble",
    "DYCOMS_RF01",
    "GABLS",
    "GATE_III",
    "life_cycle_Tan2018",
    "Nieuwstadt",
    "Rico",
    "Soares",
    "SP",
    "TRMM_LBA",
]

filter!(x -> x ≠ "GATE_III", all_cases) # no mse tables for GATE_III

dict = OrderedDict()
for case in all_cases
    dict[case] = JSON.parsefile("computed_mse_$case.json"; dicttype = OrderedCollections.OrderedDict)
end

println("#################################")
println("################################# MSE tables")
println("#################################")
println()

println("all_best_mse = OrderedDict()\n")
for case in keys(dict)
    println("all_best_mse[\"$case\"] = OrderedDict()")
    for var in keys(dict[case])
        println("all_best_mse[\"$case\"][\"$var\"] = $(dict[case][var])")
    end
    println()
end

println("#################################")
println("#################################")
println("#################################")

# Cleanup
for case in all_cases
    rm("computed_mse_$case.json"; force = true)
end
