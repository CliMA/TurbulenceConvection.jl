# ======================================================================================================================================== #
#= Supported types are:
:Base
:exponential_T_scaling_ice
:exponential_T_scaling_ice_raw
:powerlaw_T_scaling_ice
:exponential_times_powerlaw_scaling_ice
:geometric_liq__geometric_ice
:geometric_liq__exponential_T_scaling_ice
:geometric_liq__powerlaw_T_scaling_ice
:geometric_liq__exponential_T_scaling_and_geometric_ice
:linear_combination
:linear_combination_with_w
:neural_network
:raymond_ice_test
=#

# ======================================================================================================================================== #

# see types.jl where these are used, and stored in NonEquilibriumMoisture() in EDMFModel
# Look into using Functors for this, see https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects

# ---------------------------------------------------------------------------------------------------------------------------------------- #
# testing for morrisonmilbrandt
abstract type AbstractSaturationRegime end

struct Supersaturated{HL, HI, BF} <: AbstractSaturationRegime
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool
end
Supersaturated(has_liquid::Bool, has_ice::Bool, below_freezing::Bool) = Supersaturated{has_liquid, has_ice, below_freezing}(has_liquid, has_ice, below_freezing)

struct WBF{HL, HI, BF} <: AbstractSaturationRegime
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool # allow this type in lieu of created another type for when we're in WBF regime but not below freezing
end
WBF(has_liquid::Bool, has_ice::Bool, below_freezing::Bool) = WBF{has_liquid, has_ice, below_freezing}(has_liquid, has_ice, below_freezing)
WBF(has_liquid::Bool, has_ice::Bool) = WBF{has_liquid, has_ice, true}(has_liquid, has_ice, true)

struct LowerSatLine{HL, HI, BF} <: AbstractSaturationRegime # you can end up here if you've reached ice sat but liq effective timescale is slower than ices so you're stuck.
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool # below freezing you're stuck on the ice sat line, above freezing you're stuck on the liq sat line but not really bc there's no ice formation
end


struct Subsaturated{HL, HI, BF}  <: AbstractSaturationRegime
    has_liquid::Bool
    has_ice::Bool
    below_freezing::Bool # does this really matter for subsaturated? (using it bc whether q_liq_sat or q_ice_sat is lower depends on it)
end
Subsaturated(has_liquid::Bool, has_ice::Bool, below_freezing::Bool) = Subsaturated{has_liquid, has_ice, below_freezing}(has_liquid, has_ice, below_freezing)



# ---------------------------------------------------------------------------------------------------------------------------------------- #
abstract type AbstractNonEquillibriumSourcesType end
abstract type AbstractRelaxationTimescaleType <: AbstractNonEquillibriumSourcesType end
struct KorolevMazin2007 <: AbstractNonEquillibriumSourcesType end  # placeholder to keep that code around
struct RelaxToEquilibrium <: AbstractNonEquillibriumSourcesType end  # placeholder to keep that code around

Base.broadcastable(x::AbstractNonEquillibriumSourcesType) = Ref(x) # permit broadcasting over these types, e.g. f.(::AbstractNonEquillibriumSourcesType, T, p, q,...)

const valid_relaxation_timescale_types::Set{Symbol} = Set([
    :Base,
    :exponential_T_scaling_ice,
    :exponential_T_scaling_ice_raw,
    :powerlaw_T_scaling_ice,
    :exponential_times_powerlaw_scaling_ice,
    :geometric_liq__geometric_ice,
    :geometric_liq__exponential_T_scaling_ice,
    :geometric_liq__powerlaw_T_scaling_ice,
    :geometric_liq__exponential_T_scaling_and_geometric_ice,
    :linear_combination,
    :linear_combination_with_w,
    :neural_network,
    :raymond_ice_test,]
)
# ---------------------------------------------------------------------------------------------------------------------------------------- #
struct RelaxationTimescaleArgs{FT}
    min_τ_liq::FT
    min_τ_ice::FT
    max_τ_liq::FT
    max_τ_ice::FT
end

"""
Trying to use this to reduce duplicated code
"""
function get_relaxation_timescale_args(param_set::APS) 
    FT = eltype(param_set)
    return RelaxationTimescaleArgs{FT}(
        TCP.get_isbits_nt(param_set.user_params, :min_τ_liq, FT(0)),
        TCP.get_isbits_nt(param_set.user_params, :min_τ_ice, FT(0)),
        TCP.get_isbits_nt(param_set.user_params, :max_τ_liq, FT(Inf)),
        TCP.get_isbits_nt(param_set.user_params, :max_τ_ice, FT(Inf)),
    )
end

# ---------------------------------------------------------------------------------------------------------------------------------------- #

# :Base
struct BaseRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    τ_liq::FT
    τ_ice::FT
    args::RelaxationTimescaleArgs{FT}
end

# :exponential_T_scaling_ice
struct ExponentialTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1i::FT
    c_2i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :exponential_T_scaling_ice_raw
struct ExponentialTScalingIceRawRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1i::FT
    c_2i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :powerlaw_T_scaling_ice
struct PowerlawTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1i::FT
    c_2i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# # :exponential_times_powerlaw_scaling_ice
struct ExponentialTimesPowerlawScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    # not implemented yet
end

# :geometric_liq__geometric_ice
struct GeometricLiqGeometricIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_1i::FT
    c_2i::FT
    c_3i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :geometric_liq__exponential_T_scaling_ice
struct GeometricLiqExponentialTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_1i::FT
    c_2i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :geometric_liq__powerlaw_T_scaling_ice
struct GeometricLiqPowerlawTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_1i::FT
    c_2i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
struct GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_1i::FT
    c_2i::FT
    c_3i::FT
    c_4i::FT
    c_5i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :linear_combination
struct LinearCombinationRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_1i::FT
    c_2i::FT
    c_3i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :linear_combination_with_w
struct LinearCombinationWithWRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_4l::FT
    c_1i::FT
    c_2i::FT
    c_3i::FT
    c_4i::FT
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end

# :neural_network
struct NeuralNetworkRelaxationTimescale{FT, MRELT,  Nx} <: AbstractRelaxationTimescaleType
    model_x_0_characteristic::NTuple{Nx, FT}
    # neural_network_params::NTuple{Nnn, FT} # this wouldn't be immutable and thus not `isbits`
    model_re_location::MRELT # should be Val{Symbol(model_re_location)}?
    args::RelaxationTimescaleArgs{FT}
end
# NeuralNetworkRelaxationTimescale(model_x_0_characteristic::NTuple{Nx, FT}, model_re_location::MRELT, args::RelaxationTimescaleArgs{FT}) where {FT, MRELT, Nx} = NeuralNetworkRelaxationTimescale{FT, MRELT, Nx}(model_x_0_characteristic, model_re_location, args)

# :raymond_ice_test
struct RaymondIceTestRelaxationTimescale{FT, N0T <: Union{FT, NamedTuple}} <: AbstractRelaxationTimescaleType
    N_r_closure::AbstractNRClosureType
    N_0::N0T # either a value or a named tuple (;z=Tuple, values=Tuple)
    args::RelaxationTimescaleArgs{FT}
end



# ======================================================================================================================================== #
# ======================================================================================================================================== #

get_relaxation_timescale_type(relaxation_timescale_type::Symbol, param_set::APS, microphys_params::ACMP) = get_relaxation_timescale_type(Val(relaxation_timescale_type), param_set, microphys_params)
get_relaxation_timescale_type(relaxation_timescale_type::Symbol, param_set::APS) = get_relaxation_timescale_type(Val(relaxation_timescale_type), param_set, TCP.microphysics_params(param_set))
get_relaxation_timescale_type(param_set::APS) = get_relaxation_timescale_type(param_set.user_args.supersaturation_type, param_set, TCP.microphysics_params(param_set))


# ---------------------------------------------------------------------------------------------------------------------------------------- #


get_adjust_ice_N(relaxation_timescale::AbstractRelaxationTimescaleType) = hasproperty(relaxation_timescale, :adjust_ice_N) ? relaxation_timescale.adjust_ice_N : false
# ---------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------- #

# :Base
function get_relaxation_timescale_type(::Val{:Base}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    # tau_weights = TCP.get_isbits_nt(param_set.user_params, :tau_weights, nothing) # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
    # if isnothing(tau_weights)


    if hasproperty(param_set.user_params, :tau_weights) # or hasfield(typeof(param_set.user_params), :tau_weights)
        liq_params, ice_params = param_set.user_params.tau_weights.liq.liq_params, param_set.user_params.tau_weights.ice.ice_params # gets used in eval below
        τ_liq, τ_ice = FT(10) .^ liq_params.log10_tau_liq, FT(10) .^ ice_params.log10_tau_ice # log fcn hand implementation...
    else
        τ_liq = CMNe.τ_relax(microphys_params, liq_type)
        τ_ice = CMNe.τ_relax(microphys_params, ice_type)
    end

    return BaseRelaxationTimescale(τ_liq, τ_ice, get_relaxation_timescale_args(param_set))
end

#: exponential_T_scaling_ice
function get_relaxation_timescale_type(::Val{:exponential_T_scaling_ice}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    c_1i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_ice_c_1, FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_2i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_ice_c_2, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)
    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)
    return ExponentialTScalingIceRelaxationTimescale(c_1i, c_2i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :exponential_T_scaling_ice_raw
function get_relaxation_timescale_type(::Val{:exponential_T_scaling_ice_raw}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    q_i0 = FT(1e-7)
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    D = FT(0.0000226)
    # derived from typical values and assume q_i = 4/3 * π * r^3 * ρ_i * N_i and N_i = c_1 e^(c_2 T)
    c_1i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_ice_c_1, FT(4 * π * D * (q_i0 / (4 / 3) * π * ρ_i)^(1 / 3) * (0.02)^(2 / 3)),) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_2i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_ice_c_2, FT(-0.6 * 2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022)
    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)
    return ExponentialTScalingIceRawRelaxationTimescale(c_1i, c_2i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :powerlaw_T_scaling_ice
function get_relaxation_timescale_type(::Val{:powerlaw_T_scaling_ice}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    c_1i = TCP.get_isbits_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_1, FT(-9)) # F23 (values taken from Frostenberg 2022)
    c_2i = TCP.get_isbits_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_2, FT(9)) # F23 (values taken from Frostenberg 2022)
    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)
    return PowerlawTScalingIceRelaxationTimescale(c_1i, c_2i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :exponential_times_powerlaw_scaling_ice
function get_relaxation_timescale_type(::Val{:exponential_times_powerlaw_scaling_ice}, param_set::APS, microphys_params::ACMP)
    error("NotImplmentedError: This relaxation_timescale_type functionality has not been implemented yet")
end

# :geometric_liq__geometric_ice
function get_relaxation_timescale_type(::Val{:geometric_liq__geometric_ice}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)   # only for prior

    c_1l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

    c_1i = TCP.get_isbits_nt(param_set.user_params, :geometric_ice_c_1, FT(1 / (4 / 3 * π * ρ_i * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2i = TCP.get_isbits_nt(param_set.user_params, :geometric_ice_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3i = TCP.get_isbits_nt(param_set.user_params, :geometric_ice_c_3, FT(N_i0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)
    return GeometricLiqGeometricIceRelaxationTimescale(c_1l, c_2l, c_3l, c_1i, c_2i, c_3i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :geometric_liq__exponential_T_scaling_ice
function get_relaxation_timescale_type(::Val{:geometric_liq__exponential_T_scaling_ice}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)

    c_1l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

    c_1i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_ice_c_1, FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_2i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_ice_c_2, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)
   
    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)
    return GeometricLiqExponentialTScalingIceRelaxationTimescale(c_1l, c_2l, c_3l, c_1i, c_2i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end


# :geometric_liq__powerlaw_T_scaling_ice
function get_relaxation_timescale_type(::Val{:geometric_liq__powerlaw_T_scaling_ice}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)

    c_1l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

    c_1i = TCP.get_isbits_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_1, FT(-9)) # F23 (values taken from Frostenberg 2022)
    c_2i = TCP.get_isbits_nt(param_set.user_params, :powerlaw_T_scaling_ice_c_2, FT(9)) # F23 (values taken from Frostenberg 2022)
    
    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)
    return GeometricLiqPowerlawTScalingIceRelaxationTimescale(c_1l, c_2l, c_3l, c_1i, c_2i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
function get_relaxation_timescale_type(::Val{:geometric_liq__exponential_T_scaling_and_geometric_ice}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
    D = FT(0.0000226)
    c_1l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since

    c_1i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_1, FT(1/(4/3 * π * ρ_i * r_r^2))) # Yeahhh.... idk for this one lol... just combined them serially from the homogenous case where c_3 is -1/3, and used .02 as the prefactor
    # c_2i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1 -- should this be the same as c_2g? It's the same mixing... 
    c_2i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_2, FT(1 / 2.0)) # Halfway between 0 and 1 -- should this be the same as c_2g? It's the same mixing... 
    c_3i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_3, FT(-7)) # just avalue
    c_4i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_4, FT(0.02), ) # Fletcher 1962 (values taken from Frostenberg 2022) and used .02 as the prefactor
    c_5i = TCP.get_isbits_nt(param_set.user_params, :exponential_T_scaling_and_geometric_ice_c_5, FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)

    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)

    return GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale(c_1l, c_2l, c_3l, c_1i, c_2i, c_3i, c_4i, c_5i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :linear_combination
function get_relaxation_timescale_type(::Val{:linear_combination}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

    c_1l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_1, FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_2, FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = TCP.get_isbits_nt(param_set.user_params, :geometric_liq_c_3, FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

    c_1i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_1, FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
    c_2i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
    c_3i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_3, FT(2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
   
    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)

    return LinearCombinationRelaxationTimescale(c_1l, c_2l, c_3l, c_1i, c_2i, c_3i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end

# :linear_combination_with_w
function get_relaxation_timescale_type(::Val{:linear_combination_with_w}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    T_fr = TCP.T_freeze(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    w_0 = FT(1e-3) # 1 mm/s

    c_1l = TCP.get_isbits_nt(param_set.user_params, :linear_combination_liq_c_1, FT(N_l0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
    c_2l = TCP.get_isbits_nt(param_set.user_params, :linear_combination_liq_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
    c_3l = TCP.get_isbits_nt(param_set.user_params, :linear_combination_liq_c_3, FT(-10)) # asssume nothing here? (keep 0 as upper bound?) 
    c_4l = TCP.get_isbits_nt(param_set.user_params, :linear_combination_liq_c_4, FT(0))

    c_1i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_1, FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
    c_2i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_2, FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
    c_3i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_3, FT(-10)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
    c_4i = TCP.get_isbits_nt(param_set.user_params, :linear_combination_ice_c_4, FT(0)) # start at 0

    adjust_ice_N = TCP.get_isbits_nt(param_set.user_params, :adjust_ice_N, false)

    return LinearCombinationWithWRelaxationTimescale(c_1l, c_2l, c_3l, c_4l, c_1i, c_2i, c_3i, c_4i, adjust_ice_N, get_relaxation_timescale_args(param_set))
end


# :neural_network
function get_relaxation_timescale_type(::Val{:neural_network}, param_set::APS, microphys_params::ACMP)
    neural_network_params = param_set.user_params.neural_microphysics_relaxation_network # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
    neural_network_params = Float32.(collect(neural_network_params)) # collect from ntuple, then convert to Float32 for NN since that's what it's supposed to be (save eltype in the jld2?)
    model_x_0_characteristic = TCP.get_isbits_nt(param_set.user_params, :model_x_0_characteristic, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)
    # model_re_location = TCP.get_isbits_nt(param_set.user_params, :model_re_location, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)
    model_re_location_isbits = param_set.user_params.model_re_location # use the isbits version in the return so it can go into edmf and ode solver... (do we really need to store this? we read it every time...)
    TC = TurbulenceConvection

    # Redefine no matter what bc if you're in an active session, you may need to redefine the NN [can't be const for same reason]
    # This still should only happen once during the run at the beginning since later get_τ will already have this object created.
    # This still beats defining it again at every use and also should beat reading from disk every time.


    model_re_location = string(TCP.unwrap_val(model_re_location_isbits))
    if !isnothing(model_re_location)
        model_re_location = string(model_re_location) # convert symbol to string...
    else
        error("model_re_location must be specified in the namelist")
    end

    re = model_destructure_re_from_file(model_re_location) # get the reconstruction function from the file ( this is probably slow, we could pass the string repr)
    neural_network = vec_to_NN(neural_network_params, re) # construct the NN from the parameters

    if isdefined(TC, :neural_network) # do this so we don't have to read from disk and recreate the network evry time bc that's super slow
        @info "neural_network was already defined, replacing..."
        @eval neural_network = $neural_network  # eval should run in globl scope? $ to interpolate the value of the variable. Can't store in the type because it's not isbits
    else
        @info "neural_network was not defined, defining..."
        @eval neural_network::typeof($neural_network) = $neural_network  # eval should run in globl scope? $ to interpolate the value of the variable. Can't store in the type because it's not isbits
    end

    # if isdefined(TC, :neural_network) # do this so we don't have to read from disk and recreate the network evry time bc that's super slow
    #     # TC.neural_network .= re(neural_network_params) # set the parameters to the ones we just read in      # These don't work
    #     # TC.neural_network = re(neural_network_params) # set the parameters to the ones we just read in  [not sure if changing the entire object is bad but idk how to just update all the parameters in the chain)   # These don't work
    #     @error("We're not using the NN definition so shouldn't be here...")                    
    # else
    #     if @isdefined neural_network
    #        @info "neural_network is defined, but may need to be overwritten (e.g. in repeated runs)"
    #     else
    #         @info "neural_network is not defined"
    #     end
    #     model_re_location = string(TCP.unwrap_val(model_re_location_isbits))
    #     if !isnothing(model_re_location)
    #         model_re_location = string(model_re_location) # convert symbol to string...
    #     else
    #         error("model_re_location must be specified in the namelist")
    #     end

    #     re = model_destructure_re_from_file(model_re_location) # get the reconstruction function from the file ( this is probably slow, we could pass the string repr)
    #     neural_network = vec_to_NN(neural_network_params, re) # construct the NN from the parameters

    #     # Either do this or tore re and params and rebuild the NN every time
    #     # isn't this bad for reuse? like wouldn't it break in the new version where we reuse the same job? (actually it should be changeable as long as it's the same type...)
    #     # but either way it would need to be changed in between runs....
    #     @eval const neural_network::typeof($neural_network) = $neural_network  # eval should run in globl scope? $ to interpolate the value of the variable. Can't store in the type because it's not isbits
    #     # global const neural_network::typeof(neural_network) = neural_network  # can't run bc not at top level
    #     if @isdefined neural_network
    #         @info "neural_network is now defined"
    #     else
    #         @error "neural_network is not defined but should be"
    #     end
    #     if isdefined(TC, :neural_network) # do this so we don't have to read from disk and recreate the network evry time bc that's super slow
    #         @info "TC.neural_network is now defined"
    #     else
    #         @info "TC.neural_network is still not defined"
    #     end
    # end
    return NeuralNetworkRelaxationTimescale(model_x_0_characteristic, model_re_location_isbits, get_relaxation_timescale_args(param_set))
end

# :raymond_ice_test
function get_relaxation_timescale_type(::Val{:raymond_ice_test}, param_set::APS, microphys_params::ACMP)
    FT = eltype(param_set)
    N_0 = TCP.get_isbits_nt(param_set.user_args, :N_0, FT(100 * 10e6)) # 100 cm^-3
    N_r_closure = TCP.get_isbits_nt(param_set.user_args, :N_r_closure, MonodisperseNRClosure())
    if N_r_closure === :monodisperse
        N_r_closure = MonodisperseNRClosure()
    elseif N_r_closure === :fixed_radius
        N_r_closure = FixedRadiusNRClosure()
    elseif N_r_closure === :inhomogeneous
        N_r_closure = InhomogeneousNRClosure()
    else
        error("NotImplementedError: N_r_closure $N_r_closure not implemented")
    end

    return RaymondIceTestRelaxationTimescale( N_r_closure, N_0, get_relaxation_timescale_args(param_set))
end

