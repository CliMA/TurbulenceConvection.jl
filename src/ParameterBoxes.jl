

## Precipitation Parameter Boxes
abstract type AbstractPrecipitationParameters end
struct NoPrecipitationParameters <: AbstractPrecipitationParameters end

struct CutOffPrecipitationParameters <: AbstractPrecipitationParameters end

struct Clima1MParameters{FT} <: AbstractPrecipitationParameters
    LH_s0::FT
    LH_v0::FT
    gas_constant::FT
    molmass_water::FT
    R_v::FT
    TPS::ThermodynamicsParameters
    MPS::Microphysics_1M_Parameters
end
function Clima1MParameters(param_set, TPS::ThermodynamicsParameters{FT} , MPS::Microphysics_1M_Parameters) where {FT}

    aliases = ["LH_s0", "LH_v0", "gas_constant", "molmass_water"]

    (LH_s0, LH_v0, gas_constant, molmass_water) = CLIMAParameters.get_parameter_values!(param_set, aliases, "Clima1M")

    #derived parameters
    R_v = gas_constant / molmass_water

    return Clima1MParameters{get_parametric_type(param_set)}(LH_s0, LH_v0, gas_constant, molmass_water, R_v, TPS, MPS)
end

## Entrainment Parameter Boxes

abstract type AbstractEntrainmentDimScaleParameters end

struct BuoyVelScaleParameters{FT} <: AbstractEntrainmentDimScaleParameters
    w_min::FT
    c_λ::FT
end
function BuoyVelScaleParameters(param_set)
    aliases = ["w_min", "c_λ"]
    (w_min, c_λ) = CLIMAParameters.get_parameter_values!(param_set, aliases, "BuoyVelScale")
    return BuoyVelScaleParameters{get_parametric_type(param_set)}(w_min, c_λ)
end

struct InvZScaleParameters{FT} <: AbstractEntrainmentDimScaleParameters end
struct InvMeterScaleParameters{FT} <: AbstractEntrainmentDimScaleParameters end

# the parameter structs mimic these abstract categories
abstract type AbstractEntrainmentClosureParameters end
abstract type AbstractEntrainmentClosureNonlocalParameters end

struct MDEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_ε::FT
    μ_0::FT
    β::FT
    χ::FT
    c_δ::FT
end
function MDEntrainmentParameters(param_set)

    aliases = ["w_min", "c_ε", "μ_0", "β", "χ", "c_δ"]
    (w_min, c_ε, μ_0, β, χ, c_δ) = CLIMAParameters.get_parameter_values!(param_set, aliases, "MDEntrainment")


    return MDEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(w_min, c_ε, μ_0, β, χ, c_δ)

end

struct FNOEntrainmentParameters{FT} <: AbstractEntrainmentClosureNonlocalParameters
    c_fno::AbstractVector{FT}
end
function FNOEntrainmentParameters(param_set)

    aliases = ["c_fno"]
    (c_fno) = CLIMAParameters.get_parameter_array_values!(param_set, aliases, "FNOEntrainment")


    return FNOEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(c_fno)

end

struct NNEntrainmentNonlocalParameters{FT} <: AbstractEntrainmentClosureNonlocalParameters
    c_gen::AbstractVector{FT}
end
function NNEntrainmentNonlocalParameters(param_set)

    array_aliases = ["c_gen"]
    (c_gen) = CLIMAParameters.get_parameter_array_values!(param_set, array_aliases, "NNEntrainmentNonlocal")


    return NNEntrainmentNonlocalParameters{CLIMAParameters.get_parametric_type(param_set)}(c_gen)

end

struct NNEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_gen::AbstractVector{FT}
end
function NNEntrainmentParameters(param_set)
    aliases = ["w_min"]
    (w_min) = CLIMAParameters.get_parameter_values!(param_set, aliases, "NNEntrainment")
    array_aliases = ["c_gen"]
    (c_gen) = CLIMAParameters.get_parameter_values!(param_set, array_aliases, "NNEntrainment")


    return NNEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(w_min, c_gen)

end

struct LinearEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_gen::AbstractVector{FT}
end
function LinearEntrainmentParameters(param_set)
    aliases = ["w_min"]
    (w_min) = CLIMAParameters.get_parameter_values!(param_set, aliases, "LinearEntrainment")
    array_aliases = ["c_gen"]
    (c_gen) = CLIMAParameters.get_parameter_values!(param_set, array_aliases, "LinearEntrainment")


    return LinearEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(w_min, c_gen)

end


struct RFEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_rf_fix::AbstractVector{FT}
    c_rf_opt::AbstractVector{FT}
end
function RFEntrainmentParameters(param_set)
    aliases = ["w_min"]
    (w_min) = CLIMAParameters.get_parameter_values!(param_set, aliases, "RFEntrainment")
    array_aliases = ["c_rf_fix", "c_rf_opt"]
    (c_rf_fix, c_rf_opt) = CLIMAParameters.get_parameter_values!(param_set, array_aliases, "RFEntrainment")


    return RFEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(w_min, c_rf_fix, c_rf_opt)

end

# parameters for the stochastic types
abstract type AbstractStochasticEntrainmentClosureParameters end

struct LogNormalScalingProcessParameters{FT, AECPS} <: AbstractStochasticEntrainmentClosureParameters
    c_gen_stoch::AbstractVector{FT}
    ECPS::AECPS
end
function LogNormalScalingProcessParameters(param_set, ECPS::AECPS) where {AECPS <: AbstractEntrainmentClosureParameters}
    array_aliases = ["c_gen_stoch"]
    (c_gen_stoch) = CLIMAParameters.get_parameter_values!(param_set, array_aliases, "LogNormalScalingProcess")

    return LogNormalScalingProcessParameters{CLIMAParameters.get_parametric_type(param_set), AECPS}(c_gen_stoch, ECPS)
end

struct NoisyRelaxationProcessParameters{FT, AECPS} <: AbstractStochasticEntrainmentClosureParameters
    c_gen_stoch::AbstractVector{FT}
    ECPS::AECPS
end
function NoisyRelaxationProcessParameters(param_set, ECPS::AECPS) where {AECPS <: AbstractEntrainmentClosureParameters}
    array_aliases = ["c_gen_stoch"]
    (c_gen_stoch) = CLIMAParameters.get_parameter_values!(param_set, array_aliases, "NoisyRelaxationProcess")

    return NoisyRelaxationProcessParameters{CLIMAParameters.get_parametric_type(param_set), AECPS}(c_gen_stoch, ECPS)
end

struct PrognosticNoisyRelaxationProcessParameters{FT, AECPS} <: AbstractStochasticEntrainmentClosureParameters
    ECPS::AECPS
end
function PrognosticNoisyRelaxationProcessParameters(
    param_set,
    ECPS::AECPS,
) where {AECPS <: AbstractEntrainmentClosureParameters}

    return PrognosticNoisyRelaxationProcessParameters{CLIMAParameters.get_parametric_type(param_set), AECPS}(ECPS)
end



struct EDMFParameters{FT, AECPS, AEDSPS, APPS}
    a_surf::FT
    a_min::FT
    a_max::FT
    Pr_n::FT
    α_d::FT
    α_a::FT
    α_b::FT
    ECPS::AECPS
    EDSPS::AEDSPS
    PPS::APPS
end


function EDMFParameters(
    param_set,
    ECPS::AECPS,
    EDSPS::AEDSPS,
    PPS::APPS,
) where {AECPS <: Union{AbstractEntrainmentClosureParameters,AbstractStochasticEntrainmentClosureParameters}, AEDSPS <: AbstractEntrainmentDimScaleParameters, APPs <: AbstractPrecipitationParameters}

    aliases = ["a_surf","a_min", "a_max", "Prandtl_air", "α_d", "α_a", "α_b"]

    (a_surf, a_min, a_max, Prandtl_air, α_d, α_a, α_b) = CLIMAParameters.get_parameter_values!(param_set,aliases,"EDMF")

    return EDMFParameters{get_parametric_type(param_set),AECPS,AEDSPS,APPS}(
        a_surf,
        a_min,
        a_max,
        Prandtl_air,
        α_d,
        α_a,
        α_b,
        ECPS,
        EDSPS,
        PPS,
    )

end

struct TurbulenceConvectionParameters{FT, APPS}
    EDMFPS::EDMFParameters
    TPS::ThermodynamicsParameters
    PPS::APPS
end

struct TurbulenceConvectionParameters(
    param_set,
    EDMFPS::EDMFParameters,
    TPS::ThermodynamicsParameters,
    PPS::APPS,
) where {APPs <: AbstractPrecipitationParameters}

    aliases = []

    () = CLIMAParameters.get_parameter_values!(param_set,aliases,"TurbulenceConvection")


    return EDMFParameters{get_parametric_type(param_set),APPS}(
        EDMFPS,
        TPS,
        PPS,
    )
end
