@testset "Post-run tests" begin
    isnan_or_inf(x) = isnan(x) || isinf(x)
    NC.Dataset(ds_tc_filename, "r") do ds
        profile = ds.group["profiles"]
        @test !any(isnan_or_inf.(Array(profile["qt_mean"])))
        @test !any(isnan_or_inf.(Array(profile["updraft_area"])))
        @test !any(isnan_or_inf.(Array(profile["updraft_w"])))
        @test !any(isnan_or_inf.(Array(profile["updraft_qt"])))
        @test !any(isnan_or_inf.(Array(profile["updraft_thetal"])))
        @test !any(isnan_or_inf.(Array(profile["u_mean"])))
        @test !any(isnan_or_inf.(Array(profile["tke_mean"])))
        @test !any(isnan_or_inf.(Array(profile["temperature_mean"])))
        @test !any(isnan_or_inf.(Array(profile["ql_mean"])))
        @test !any(isnan_or_inf.(Array(profile["qi_mean"])))
        @test !any(isnan_or_inf.(Array(profile["thetal_mean"])))
        @test !any(isnan_or_inf.(Array(profile["Hvar_mean"])))
        @test !any(isnan_or_inf.(Array(profile["QTvar_mean"])))
        @test !any(isnan_or_inf.(Array(profile["v_mean"])))
        @test !any(isnan_or_inf.(Array(profile["qr_mean"])))
    end
    nothing
end

@testset "Simulation completion" begin
    # Test that the simulation has actually finished,
    # and not aborted early.
    @test !(return_code == :simulation_aborted)
    @test return_code == :success
end
