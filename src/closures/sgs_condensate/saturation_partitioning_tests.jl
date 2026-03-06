using Test
using Distributions

# Assuming cond_qt_SD is defined here or imported

@testset "cond_qt_SD normalized shift tests (Pure Physics)" begin
    FT = Float64

    @testset "Case 1: Zero Variance (Delta Distribution)" begin
        qt, qt_var, q_sat, qc_real = FT(0.01), FT(0.0), FT(0.015), FT(0.0)

        qt_cond, cond_sd = cond_qt_SD(qt, qt_var, q_sat, qc_real)

        println("--- Case 1: Zero Variance ---")
        println("Call: cond_qt_SD(qt=$qt, qt_var=$qt_var, q_sat=$q_sat, qc_real=$qc_real)")
        println("Output: qt_cond = $qt_cond, cond_sd = $cond_sd\n")

        @test qt_cond == qt
        @test cond_sd == zero(FT)
    end

    @testset "Case 2: 50% Saturated (Mean sits exactly on q_sat)" begin
        qt, q_sat, qt_var, qc_real = FT(0.01), FT(0.01), FT(1e-6), FT(0.0)
        qt_std = sqrt(qt_var)

        expected_shift = qt_std * sqrt(FT(π) / 2)
        expected_sd_shift = sqrt(FT(π) / 2) # ≈ 1.253

        # Test A: Uncapped (pass σ_max = 5.0) to verify the pure math
        qt_cond, cond_sd_uncapped = cond_qt_SD(qt, qt_var, q_sat, qc_real, FT(5))

        # Test B: Capped (default σ_max = 1.0)
        _, cond_sd_capped = cond_qt_SD(qt, qt_var, q_sat, qc_real)

        println("--- Case 2: 50% Saturated ---")
        println("Call (Uncapped): cond_qt_SD(..., σ_max=5.0)")
        println("Output (Uncapped): qt_cond = $qt_cond, cond_sd = $cond_sd_uncapped")
        println("Expected (Uncapped): cond_sd = $expected_sd_shift")
        println("Output (Capped default): cond_sd = $cond_sd_capped\n")

        @test qt_cond ≈ (qt + expected_shift) rtol = 1e-5
        @test cond_sd_uncapped ≈ expected_sd_shift rtol = 1e-5
        @test cond_sd_capped == FT(1)
    end

    @testset "Case 3: Heavily Sub-saturated (Asymptotic limit)" begin
        # 10 standard deviations sub-saturated (z* ≈ 10.0)
        qt, q_sat, qt_var, qc_real = FT(0.005), FT(0.015), FT(1e-6), FT(1e-4)
        qt_std = sqrt(qt_var)

        # Test A: Uncapped (pass σ_max = 15.0 to see the pure z* limit)
        qt_cond_uncapped, cond_sd_uncapped = cond_qt_SD(qt, qt_var, q_sat, qc_real, FT(15))

        # Test B: Capped (default σ_max = 1.0)
        qt_cond_capped, cond_sd_capped = cond_qt_SD(qt, qt_var, q_sat, qc_real)

        println("--- Case 3: Heavily Sub-saturated ---")
        println("Call (Uncapped): cond_qt_SD(..., σ_max=15.0)")
        println("Output (Uncapped): cond_sd = $cond_sd_uncapped (Expect ≈ 10.0)")
        println("Call (Capped default): cond_qt_SD(..., σ_max=1.0)")
        println("Output (Capped): cond_sd = $cond_sd_capped (Expect exactly 1.0)\n")

        # Use ≈ (isapprox) for floating point comparisons!
        @test cond_sd_uncapped ≈ FT(10) rtol = 1e-5
        @test qt_cond_uncapped ≈ (qt + FT(10) * qt_std) rtol = 1e-5
        @test cond_sd_capped == FT(1) # We can use == here because min() outputs exact constants
        @test qt_cond_capped ≈ (qt + FT(1) * qt_std) rtol = 1e-5
    end

    @testset "Case 4: Highly Supersaturated (Grid is completely clouded)" begin
        qt, q_sat, qt_var, qc_real = FT(0.02), FT(0.01), FT(1e-6), FT(0.0)
        qt_std = sqrt(qt_var)

        qc_expected = qt - q_sat
        expected_shift = qt_var / qc_expected
        expected_sd_shift = expected_shift / qt_std

        qt_cond, cond_sd = cond_qt_SD(qt, qt_var, q_sat, qc_real)

        println("--- Case 4: Highly Supersaturated ---")
        println("Output: qt_cond = $qt_cond, cond_sd = $cond_sd")
        println("Expected: qt_cond = $(qt + expected_shift), cond_sd = $expected_sd_shift\n")

        @test qt_cond ≈ (qt + expected_shift) rtol = 1e-5
        @test cond_sd ≈ expected_sd_shift rtol = 1e-5
    end

    @testset "Case 5: Type Stability (Float32)" begin
        FT32 = Float32
        qt, qt_var, q_sat, qc_real = FT32(0.01), FT32(1e-6), FT32(0.01), FT32(1e-4)

        qt_cond, cond_sd = cond_qt_SD(qt, qt_var, q_sat, qc_real)

        println("--- Case 5: Float32 Type Stability ---")
        println("Output types: qt_cond is $(typeof(qt_cond)), cond_sd is $(typeof(cond_sd))\n")

        @test typeof(qt_cond) == FT32
        @test typeof(cond_sd) == FT32
    end
end
