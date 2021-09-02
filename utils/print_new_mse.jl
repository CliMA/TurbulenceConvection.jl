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
    "LES_driven_SCM",
]

filter!(x -> x ≠ "GATE_III", all_cases) # no mse tables for GATE_III

include(joinpath("..", "integration_tests", "utils", "mse_tables.jl"))

percent_reduction_mse = Dict()

computed_mse = OrderedDict()
for case in all_cases
    computed_mse[case] = JSON.parsefile("computed_mse_$case.json"; dicttype = OrderedCollections.OrderedDict)
end

println("#################################")
println("################################# MSE tables")
println("#################################")
println("#")

println("all_best_mse = OrderedDict()\n#")
for case in keys(computed_mse)
    println("all_best_mse[\"$case\"] = OrderedDict()")
    percent_reduction_mse[case] = 0
    for var in keys(computed_mse[case])
        println("all_best_mse[\"$case\"][\"$var\"] = $(computed_mse[case][var])")
        percent_reduction_mse[case] = min(
            percent_reduction_mse[case],
            (all_best_mse[case][var] - computed_mse[case][var]) / all_best_mse[case][var] * 100,
        )
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

println("-- DO NOT COPY --")
for case in keys(percent_reduction_mse)
    @info "percent_reduction_mse[$case] = $(percent_reduction_mse[case])"
end
@info "min mse reduction (%) over all cases = $(min(values(percent_reduction_mse)...))"
