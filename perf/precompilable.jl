include(joinpath(@__DIR__, "common.jl"))
using Test

case_name = "Bomex"

println("Running $case_name...")
namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "pc_01"
t_precompile = @elapsed main(namelist)

println("Running $case_name...")
namelist = NameList.default_namelist(case_name)
namelist["meta"]["uuid"] = "pc_02"
t_precompiled = @elapsed main(namelist)

@info "Precompiling run: $(t_precompile)"
@info "Precompiled  run: $(t_precompiled)"
@info "precompiled/precompiling: $(t_precompiled/t_precompile))"

@testset "Test runtime" begin
    @test t_precompiled / t_precompile < 0.14
end
