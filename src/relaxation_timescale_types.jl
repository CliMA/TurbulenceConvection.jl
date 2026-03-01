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
:neural_network_extended
:raymond_ice_test
=#

# ======================================================================================================================================== #

# see types.jl where these are used, and stored in NonEquilibriumMoisture() in EDMFModel
# Look into using Functors for this, see https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects

# ---------------------------------------------------------------------------------------------------------------------------------------- #
# testing for morrisonmilbrandt

include("supersaturation_regimes.jl")

# ---------------------------------------------------------------------------------------------------------------------------------------- #
abstract type AbstractNonEquillibriumSourcesType end
abstract type AbstractRelaxationTimescaleType <: AbstractNonEquillibriumSourcesType end

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
    :neural_network_extended,
    :neural_network_no_weights,
    :neural_network_random_init,
    :neural_network_pca_noise,
    # :raymond_ice_test,
])
# ---------------------------------------------------------------------------------------------------------------------------------------- #
struct RelaxationTimescaleArgs{FT}
    min_τ_liq::FT
    min_τ_ice::FT
    max_τ_liq::FT
    max_τ_ice::FT
    #
    min_N_liq::FT
    min_N_ice::FT
    max_N_liq::FT
    max_N_ice::FT
end

"""
Trying to use this to reduce duplicated code
"""
function get_relaxation_timescale_args(namelist, FT)
    # FT = namelist["float_type"] == "Float32" ? Float32 : Float64 # how main() does it.
    return RelaxationTimescaleArgs{FT}(
        #
        get(namelist["user_params"], "min_τ_liq", FT(0)),
        get(namelist["user_params"], "min_τ_ice", FT(0)),
        get(namelist["user_params"], "max_τ_liq", FT(Inf)),
        get(namelist["user_params"], "max_τ_ice", FT(Inf)),
        #
        get(namelist["user_params"], "min_N_liq", FT(0)),
        get(namelist["user_params"], "min_N_ice", FT(0)),
        get(namelist["user_params"], "max_N_liq", FT(Inf)),
        get(namelist["user_params"], "max_N_ice", FT(Inf)),
    )
end

# ---------------------------------------------------------------------------------------------------------------------------------------- #


# struct KorolevMazin2007 <: AbstractNonEquillibriumSourcesType end  # placeholder to keep that code around
struct RelaxToEquilibrium{FT} <: AbstractNonEquillibriumSourcesType
    adjust_ice_N::Bool
    args::RelaxationTimescaleArgs{FT}
end  # placeholder to keep that code around

function RelaxToEquilibrium(param_set::APS, namelist)
    FT = eltype(param_set)
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set
    return RelaxToEquilibrium{FT}(adjust_ice_N, get_relaxation_timescale_args(namelist, FT))
end

# ---------------------------------------------------------------------------------------------------------------------------------------- #

# :Base
struct BaseRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    args::RelaxationTimescaleArgs{FT}
end

# :exponential_T_scaling_ice
struct ExponentialTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1i::FT
    c_2i::FT
    adjust_ice_N::Bool
    τ_sub_dep_scaling_factor::FT
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
    τ_sub_dep_scaling_factor::FT
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
    c_4l::FT
    c_1i::FT
    c_2i::FT
    c_3i::FT
    c_4i::FT
    adjust_liq_N::Bool
    adjust_ice_N::Bool
    τ_cond_evap_scaling_factor::FT
    τ_sub_dep_scaling_factor::FT
    args::RelaxationTimescaleArgs{FT}
end

# :geometric_liq__exponential_T_scaling_ice
struct GeometricLiqExponentialTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_4l::FT
    c_1i::FT
    c_2i::FT
    adjust_liq_N::Bool
    adjust_ice_N::Bool
    τ_cond_evap_scaling_factor::FT
    τ_sub_dep_scaling_factor::FT
    args::RelaxationTimescaleArgs{FT}
end

# :geometric_liq__powerlaw_T_scaling_ice
struct GeometricLiqPowerlawTScalingIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_4l::FT
    c_1i::FT
    c_2i::FT
    adjust_liq_N::Bool
    adjust_ice_N::Bool
    τ_cond_evap_scaling_factor::FT
    τ_sub_dep_scaling_factor::FT
    args::RelaxationTimescaleArgs{FT}
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
struct GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale{FT} <: AbstractRelaxationTimescaleType
    c_1l::FT
    c_2l::FT
    c_3l::FT
    c_4l::FT
    c_1i::FT
    c_2i::FT
    # c_3i::FT # testing deprecation for now
    c_4i::FT
    c_5i::FT
    adjust_liq_N::Bool
    adjust_ice_N::Bool
    τ_cond_evap_scaling_factor::FT
    τ_sub_dep_scaling_factor::FT
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
    adjust_liq_N::Bool
    adjust_ice_N::Bool
    τ_cond_evap_scaling_factor::FT
    τ_sub_dep_scaling_factor::FT
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
    adjust_liq_N::Bool
    adjust_ice_N::Bool
    τ_cond_evap_scaling_factor::FT
    τ_sub_dep_scaling_factor::FT
    args::RelaxationTimescaleArgs{FT}
end

# :neural_network
# struct NeuralNetworkRelaxationTimescale{FTNN, FT, NNT, NNPT, MRELT,  Nx} <: AbstractRelaxationTimescaleType
struct NeuralNetworkRelaxationTimescale{FTNN, FT, NNT, NNNP, Nx} <: AbstractRelaxationTimescaleType # I gave up on NNPType in favor of enforcing FTNN within the type itself (rather than just in a constructor)
    neural_network::NNT # Now that we have SimpleChain, it is isbits :)
    # neural_network_params::NNPT # use a SVector, is isbits [seems to massively slow down init though... I wonder if using the global was faster lol..]
    neural_network_params::SA.SVector{NNNP, FTNN} # use a SVector, is isbits [seems to massively slow down init though... I wonder if using the global was faster lol..]
    model_x_0_characteristic::NTuple{Nx, FT}
    # neural_network_params::NTuple{Nnn, FT} # this wouldn't be immutable and thus not `isbits`
    # model_re_location::MRELT # should be Val{Symbol(model_re_location)}? (now we have an isbits version I dont think we need to store it)
    adjust_liq_N::Bool
    adjust_ice_N::Bool # You can hope the NN stays reasonable but if extrapolating to unseen data who knows... but it will not if not well pretrained or if you extrapolate far enough
    args::RelaxationTimescaleArgs{FT}
end

const BaseNeuralNetworkRelaxationTimescale{FTNN, FT, NNT, NNNP} =
    NeuralNetworkRelaxationTimescale{FTNN, FT, NNT, NNNP, 5}
const ExtendedNeuralNetworkRelaxationTimescale{FTNN, FT, NNT, NNNP} =
    NeuralNetworkRelaxationTimescale{FTNN, FT, NNT, NNNP, 9}

# NeuralNetworkRelaxationTimescale(model_x_0_characteristic::NTuple{Nx, FT}, model_re_location::MRELT, args::RelaxationTimescaleArgs{FT}) where {FT, MRELT, Nx} = NeuralNetworkRelaxationTimescale{FT, MRELT, Nx}(model_x_0_characteristic, model_re_location, args)

# # :raymond_ice_test
# struct RaymondIceTestRelaxationTimescale{FT, N0T <: Union{FT, NamedTuple}} <: AbstractRelaxationTimescaleType
#     N_r_closure::AbstractNRClosureType
#     N_0::N0T # either a value or a named tuple (;z=Tuple, values=Tuple)
#     args::RelaxationTimescaleArgs{FT}
# end



# ======================================================================================================================================== #
# ======================================================================================================================================== #

get_relaxation_timescale_type(relaxation_timescale_type::Symbol, param_set::APS, microphys_params::ACMP, namelist) =
    get_relaxation_timescale_type(Val(relaxation_timescale_type), param_set, microphys_params, namelist)
get_relaxation_timescale_type(relaxation_timescale_type::Symbol, param_set::APS, namelist) =
    get_relaxation_timescale_type(
        Val(relaxation_timescale_type),
        param_set,
        TCP.microphysics_params(param_set),
        namelist,
    )
get_relaxation_timescale_type(param_set::APS, namelist) = get_relaxation_timescale_type(
    get(namelist["user_args"], "nonequilibrium_moisture_scheme", :Base),
    param_set,
    TCP.microphysics_params(param_set),
    namelist,
)

# ---------------------------------------------------------------------------------------------------------------------------------------- #
get_adjust_ice_N(relaxation_timescale::AbstractRelaxationTimescaleType) =
    hasproperty(relaxation_timescale, :adjust_ice_N) ? relaxation_timescale.adjust_ice_N : false
get_adjust_liq_N(relaxation_timescale::AbstractRelaxationTimescaleType) =
    hasproperty(relaxation_timescale, :adjust_liq_N) ? relaxation_timescale.adjust_liq_N : false

get_adjust_ice_N(relaxation_timescale::RelaxToEquilibrium) = false
get_adjust_liq_N(relaxation_timescale::RelaxToEquilibrium) = false
# ---------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------- #


# :Base
function get_relaxation_timescale_type(::Val{:Base}, param_set::APS, microphys_params::ACMP, namelist) # added namelist
    FT = eltype(param_set)
    return BaseRelaxationTimescale(get_relaxation_timescale_args(namelist, FT))
end

#: exponential_T_scaling_ice
function get_relaxation_timescale_type(
    ::Val{:exponential_T_scaling_ice},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    c_1i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_1", FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_2i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_2", FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set
    return ExponentialTScalingIceRelaxationTimescale(
        c_1i,
        c_2i,
        adjust_ice_N,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :exponential_T_scaling_ice_raw
function get_relaxation_timescale_type(
    ::Val{:exponential_T_scaling_ice_raw},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    q_i0 = FT(1e-7)
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    D = FT(0.0000226)
    # derived from typical values and assume q_i = 4/3 * π * r^3 * ρ_i * N_i and N_i = c_1 e^(c_2 T)
    c_1i = get(
        namelist["relaxation_timescale_params"],
        "exponential_T_scaling_ice_c_1",
        FT(4 * π * D * (q_i0 / (4 / 3) * π * ρ_i)^(1 / 3) * (0.02)^(2 / 3)),
    ) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_2i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_2", FT(-0.6 * 2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022)
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set
    return ExponentialTScalingIceRawRelaxationTimescale(
        c_1i,
        c_2i,
        adjust_ice_N,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :powerlaw_T_scaling_ice
function get_relaxation_timescale_type(::Val{:powerlaw_T_scaling_ice}, param_set::APS, microphys_params::ACMP, namelist)
    FT = eltype(param_set)
    c_1i = get(namelist["relaxation_timescale_params"], "powerlaw_T_scaling_ice_c_1", FT(-9)) # F23 (values taken from Frostenberg 2022)
    c_2i = get(namelist["relaxation_timescale_params"], "powerlaw_T_scaling_ice_c_2", FT(9)) # F23 (values taken from Frostenberg 2022)
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set
    return PowerlawTScalingIceRelaxationTimescale(
        c_1i,
        c_2i,
        adjust_ice_N,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :exponential_times_powerlaw_scaling_ice
function get_relaxation_timescale_type(
    ::Val{:exponential_times_powerlaw_scaling_ice},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    error("NotImplmentedError: This relaxation_timescale_type functionality has not been implemented yet")
end

# :geometric_liq__geometric_ice
function get_relaxation_timescale_type(
    ::Val{:geometric_liq__geometric_ice},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)   # only for prior

    c_1l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_1", FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_3", FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_4l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_4", FT(2e2 * 1e6)) # max N acceptable

    c_1i = get(namelist["relaxation_timescale_params"], "geometric_ice_c_1", FT(1 / (4 / 3 * π * ρ_i * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2i = get(namelist["relaxation_timescale_params"], "geometric_ice_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3i = get(namelist["relaxation_timescale_params"], "geometric_ice_c_3", FT(N_i0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_4i = get(namelist["relaxation_timescale_params"], "geometric_ice_c_4", FT(2e4 * 1e6)) # max N acceptable

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set

    τ_cond_evap_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_cond_evap_scaling_factor", FT(1.0)) # default to 1.0 if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set
    return GeometricLiqGeometricIceRelaxationTimescale(
        c_1l,
        c_2l,
        c_3l,
        c_4l,
        c_1i,
        c_2i,
        c_3i,
        c_4i,
        adjust_liq_N,
        adjust_ice_N,
        τ_cond_evap_scaling_factor,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :geometric_liq__exponential_T_scaling_ice
function get_relaxation_timescale_type(
    ::Val{:geometric_liq__exponential_T_scaling_ice},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)

    c_1l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_1", FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_3", FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_4l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_4", FT(2e2 * 1e6)) # max N acceptable

    c_1i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_1", FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_2i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_2", FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set

    τ_cond_evap_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_cond_evap_scaling_factor", FT(1.0)) # default to 1.0 if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set
    return GeometricLiqExponentialTScalingIceRelaxationTimescale(
        c_1l,
        c_2l,
        c_3l,
        c_4l,
        c_1i,
        c_2i,
        adjust_liq_N,
        adjust_ice_N,
        τ_cond_evap_scaling_factor,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end


# :geometric_liq__powerlaw_T_scaling_ice
function get_relaxation_timescale_type(
    ::Val{:geometric_liq__powerlaw_T_scaling_ice},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)


    c_1l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_1", FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_3", FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_4l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_4", FT(2e2 * 1e6)) # max N acceptable

    c_1i = get(namelist["relaxation_timescale_params"], "powerlaw_T_scaling_ice_c_1", FT(-9)) # F23 (values taken from Frostenberg 2022)
    c_2i = get(namelist["relaxation_timescale_params"], "powerlaw_T_scaling_ice_c_2", FT(9)) # F23 (values taken from Frostenberg 2022)

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set

    τ_cond_evap_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_cond_evap_scaling_factor", FT(1.0)) # default to 1.0 if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set
    return GeometricLiqPowerlawTScalingIceRelaxationTimescale(
        c_1l,
        c_2l,
        c_3l,
        c_4l,
        c_1i,
        c_2i,
        adjust_liq_N,
        adjust_ice_N,
        τ_cond_evap_scaling_factor,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
function get_relaxation_timescale_type(
    ::Val{:geometric_liq__exponential_T_scaling_and_geometric_ice},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
    D = FT(0.0000226)

    c_1l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_1", FT(1 / (4 / 3 * π * ρ_l * r_r^2))) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1
    c_3l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_3", FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_4l = get(namelist["relaxation_timescale_params"], "geometric_liq_c_4", FT(2e2 * 1e6)) # max N acceptable

    # [[ deprecated ]] - use the names from :exponential_T_scaling_ice and :geometric_ice
    # c_1i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_1", FT(1/(4/3 * π * ρ_i * r_r^2))) # Yeahhh.... idk for this one lol... just combined them serially from the homogenous case where c_3 is -1/3, and used .02 as the prefactor
    # # c_2i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_2", FT(1 / 2.0)) # Halfway between 0 and 1 -- should this be the same as c_2g? It's the same mixing...
    # c_2i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1 -- should this be the same as c_2g? It's the same mixing...
    # c_3i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_3", FT(-7)) # just a value   
    # c_4i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_4", FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022) and used .02 as the prefactor
    # c_5i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_5", FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)


    c_1i = get(namelist["relaxation_timescale_params"], "geometric_ice_c_1", FT(1 / (4 / 3 * π * ρ_i * r_r^2)))
    c_2i = get(namelist["relaxation_timescale_params"], "geometric_ice_c_2", FT(2 / 3.0)) # Halfway between 1/3 and 1
    # c_3i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_and_geometric_ice_c_3", FT(-7)) #  we use c_3 differently so it keeps its special form...

    c_4i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_1", FT(0.02)) # Fletcher 1962 (values taken from Frostenberg 2022)
    c_5i = get(namelist["relaxation_timescale_params"], "exponential_T_scaling_ice_c_2", FT(-0.6)) # Fletcher 1962 (values taken from Frostenberg 2022)

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set

    τ_cond_evap_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_cond_evap_scaling_factor", FT(1.0)) # default to 1.0 if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set

    # return GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale(c_1l, c_2l, c_3l, c_4l, c_1i, c_2i, c_3i, c_4i, c_5i, adjust_liq_N, adjust_ice_N, τ_cond_evap_scaling_factor, τ_sub_dep_scaling_factor, get_relaxation_timescale_args(namelist, FT))
    return GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale(
        c_1l,
        c_2l,
        c_3l,
        c_4l,
        c_1i,
        c_2i,
        c_4i,
        c_5i,
        adjust_liq_N,
        adjust_ice_N,
        τ_cond_evap_scaling_factor,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :linear_combination
function get_relaxation_timescale_type(::Val{:linear_combination}, param_set::APS, microphys_params::ACMP, namelist)
    FT = eltype(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)

    c_1l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_1", FT(N_l0 * r_0)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define
    c_2l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_2", FT(0)) # Halfway between 1/3 and 1
    c_3l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_3", FT(-10)) # Inhomogenous default assuming r_0 = 20 micron since `typical` N is harder to define

    c_1i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_1", FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
    c_2i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_2", FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
    c_3i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_3", FT(2 / 3)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set

    τ_cond_evap_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_cond_evap_scaling_factor", FT(1.0)) # default to 1.0 if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set

    return LinearCombinationRelaxationTimescale(
        c_1l,
        c_2l,
        c_3l,
        c_1i,
        c_2i,
        c_3i,
        adjust_liq_N,
        adjust_ice_N,
        τ_cond_evap_scaling_factor,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# :linear_combination_with_w
function get_relaxation_timescale_type(
    ::Val{:linear_combination_with_w},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    T_fr = TCP.T_freeze(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_r = FT(20 * 1e-6) # 20 microns
    r_0::FT = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    N_i0 = FT(1e-7 / (4 / 3 * π * r_r^3 * ρ_i)) # estimated total N assuming reasonable q.ice... (N = N_r + N_0)
    N_l0 = FT(1e-5 / (4 / 3 * π * r_r^3 * ρ_l)) # estimated total N assuming reasonable q.liq.. (N = N_r in homogenous)
    w_0 = FT(1e-3) # 1 mm/s

    c_1l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_1", FT(N_l0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
    c_2l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_2", FT(0)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
    c_3l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_3", FT(-10)) # asssume nothing here? (keep 0 as upper bound?)
    c_4l = get(namelist["relaxation_timescale_params"], "linear_combination_liq_c_4", FT(0))

    c_1i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_1", FT(N_i0 * r_0)) # I think at q=0, we need c_1 from linear = c_1 from geometric...
    c_2i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_2", FT(-0.6)) # Halfway between 1/3 and 1 (we know these can't be right?) but it has the same sign lmao so it still decays... (we would need to figure out how to match slopes at some arbitrary point near 0 that isn't 0 lmao)
    c_3i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_3", FT(-10)) # Fletcher 1962 (values taken from Frostenberg 2022), same sign again I suppose...
    c_4i = get(namelist["relaxation_timescale_params"], "linear_combination_ice_c_4", FT(0)) # start at 0

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set

    τ_cond_evap_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_cond_evap_scaling_factor", FT(1.0)) # default to 1.0 if not set
    τ_sub_dep_scaling_factor = get(namelist["relaxation_timescale_params"], "τ_sub_dep_scaling_factor", FT(1.0)) # default to 1.0 if not set


    return LinearCombinationWithWRelaxationTimescale(
        c_1l,
        c_2l,
        c_3l,
        c_4l,
        c_1i,
        c_2i,
        c_3i,
        c_4i,
        adjust_liq_N,
        adjust_ice_N,
        τ_cond_evap_scaling_factor,
        τ_sub_dep_scaling_factor,
        get_relaxation_timescale_args(namelist, FT),
    )
end


# :neural_network
# this version accepts the namelist, so we can ideally keep the neural network params out of the param_set to reduce bloat...
function get_relaxation_timescale_type(
    ::Union{
        Val{:neural_network},
        Val{:neural_network_no_weights},
        Val{:neural_network_random_init},
        Val{:neural_network_extended},
    },
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)
    # rn everything in user_params ends up in param_set... but that isn't really necessary given we're creating these objects... it was more of a convenience when we were testing random things...
    neural_network_params =
        create_svector(FTNN.(namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network"])) # convert to SVector for NN since that's what it's supposed to be (save eltype in the jld2?)
    model_x_0_characteristic = Tuple(FT.(namelist["relaxation_timescale_params"]["model_x_0_characteristic"])) # make this match the data, the type conversion will happen in prepare_for_NN()
    model_re_location = namelist["relaxation_timescale_params"]["model_re_location"]

    neural_network = simple_chain_model_from_file(model_re_location) # construct the NN from the parameters, and also get the simple chain params for the neural network
    # @warn "neural_network = $neural_network"

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set # we could add this... hopefully the NN stays reasonable but if extrapolating to unseen data who knows...

    return NeuralNetworkRelaxationTimescale(
        neural_network,
        neural_network_params,
        model_x_0_characteristic,
        adjust_liq_N,
        adjust_ice_N,
        get_relaxation_timescale_args(namelist, FT),
    )
end


# :neural_network_pca_noise
function get_relaxation_timescale_type(
    ::Val{:neural_network_pca_noise},
    param_set::APS,
    microphys_params::ACMP,
    namelist,
)
    FT = eltype(param_set)

    pca_mean_vec = namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_mean"]
    pca_components = namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_components"]
    pca_weights = namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_weights"] # convert to SVector for NN since that's what it's supposed to be (save eltype in the jld2?)
    pca_explained_variance =
        namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_pca_explained_variance"] # this is just for info, not used in the NN

    by_output = namelist["relaxation_timescale_params"]["pca_by_output"] # if true, then the pca_components are by output, otherwise they are by input
    reference_to_truth = namelist["relaxation_timescale_params"]["reference_to_truth"] # if true, then the pca_components are in reference to the truth, otherwise they are in reference to the model
    specific_toggle = namelist["relaxation_timescale_params"]["specific_toggle"] # if true, then the pca_components are specific to the model, otherwise they are general


    # pca_explained_variance = S^2/(n_samples - 1)
    # the real noise, σ = S / sqrt(n_samples - 1) = sqrt(pca_explained_variance) # this is the standard deviation of the noise, so we can scale the components by it

    # number_samples = 10000
    # S = sqrt.(pca_explained_variance .* (number_samples - 1)) # this is the standard deviation of the noise, so we can scale the components by it

    # we want the largest pca component to to move the needle maybe 1 unit.. meaning that its noise is about 1/σ_real[1].

    # if youre not in loss framework, you are less certain, if youre in real sapce (not reference to truth) then you want df to be within 1. so you choose dp accordingly...

    # [[ Note, this latte sum(A .* B')]] is the same as if we'd just done A * B ]]

    n_pca = namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_number_components"] # number of PCA components, this is just the length of the pca_components vector



    if !by_output
        (length(pca_explained_variance) == n_pca) || error(
            "The number of PCA explained variances must match the number of PCA components, got $(length(pca_explained_variance)) explained variances and $n_pca components.",
        )
        ((1 + n_pca) == length(pca_weights)) || error(
            "The number of PCA weights must match the 1 plus the number of PCA components, got $(length(pca_weights)) instead of $(1 + n_pca) weights.",
        )
        dp_weight = pca_weights[end] # this is the weight for the total motion, which is just the last weight in the PCA weights
        pca_weights = pca_weights[1:(end - 1)] # remove the last weight
        # pca_components is [nparams x ncomponents] 
        pca_weights = reshape(pca_weights, 1, length(pca_weights)) # make sure it's a row vector
        pca_explained_variance = reshape(pca_explained_variance, 1, length(pca_explained_variance)) # make sure it's a row vector to match the weights
        σ_real = sqrt.(pca_explained_variance) # this is the standard deviation of the noise, so we can scale the components by it

        dlossdp = FTNN.(pca_mean_vec .+ sum(((pca_weights .* σ_real) .* pca_components), dims = 2))[:]  # this is the PCA noise, so we add the mean and scale by the weights, and sum over all components
        if reference_to_truth
            learning_rate = FT(0.1) * inv.(maximum(abs.(dlossdp))) # this is the learning rate, which is just a hyperparameter to control how much we wiggle the parameters, so we can adjust it later if needed
        else
            learning_rate = FT(0.5)
        end
        dp = -learning_rate .* dlossdp # scale the parameters by the learning rate, this is just a hyperparameter to control how much we wiggle the parameters, so we can adjust it later if needed
        dp .*= dp_weight
    else
        n_outputs = length(pca_mean_vec) # number of neural network outputs
        (sum(length(pev) for pev in pca_explained_variance) == n_pca) || error(
            "The number of PCA explained variances must match the number of PCA components, got $(sum(length(pev) for pev in pca_explained_variance)) explained variances and $n_pca components.",
        )
        ((n_outputs + n_pca) == length(pca_weights)) || error(
            "The number of PCA weights must match the number of outputs plus the number of PCA components, got $(length(pca_weights)) weights instead of $(n_outputs + n_pca), given $n_outputs outputs.",
        )
        dp_weights = pca_weights[(end - (n_outputs - 1)):end] # get the weights for the outputs, which are not in the PCA components
        pca_weights = pca_weights[1:(end - n_outputs)] # remove the last n_outputs weights

        dlossdps = Vector{Vector{FTNN}}(undef, n_outputs) # this is the gradient of the loss with respect to the parameters, which is just the PCA noise scaled by the weights and summed over all components
        learning_rates = Vector{FT}(undef, n_outputs)
        dps = Vector{Vector{FT}}(undef, n_outputs) # this is the gradient of the loss with respect to the parameters, which is just the PCA noise scaled by the weights and summed over all components

        pca_explained_variance = map(pev -> reshape(pev, 1, length(pev)), pca_explained_variance) # make sure each explained variance is a row vector to match the weights

        starting_ind = 1
        for i_o in 1:n_outputs
            n_pca_each = length(pca_explained_variance[i_o]) # number of PCA components for this output

            pca_weights_here = pca_weights[starting_ind:(starting_ind + n_pca_each - 1)] # get the weights for this output
            pca_weights_here = reshape(pca_weights_here, 1, length(pca_weights_here)) # make sure it's a row vector to match the weights
            starting_ind += n_pca_each # move to the next output
            σ_real = sqrt.(pca_explained_variance[i_o]) # this is the standard deviation of the noise, so we can scale the
            # @warn "σ_real[$i_o] = $(σ_real)"
            dlossdps[i_o] =
                FTNN.(pca_mean_vec[i_o] .+ sum(((pca_weights_here .* σ_real) .* pca_components[i_o]), dims = 2))[:] # this is the PCA noise, so we add the mean and scale by the weights, and sum over all components

            if reference_to_truth # we want to take a principled step in the loss direction via gradient descent.
                learning_rates[i_o] = FT(0.1) * inv.(maximum(abs.(dlossdps[i_o]))) # this is the learning rate, which is just a hyperparameter to control how much we wiggle the parameters, so we can adjust it later if needed
                learning_rates[i_o] /= n_outputs
            else # we are just stirring things up, we want to take a step such that the maximum df is about 1. so df/dp * dp = 1, dp = 1/(df/dp)
                if !specific_toggle
                    learning_rates[i_o] = FT(0.2) * inv.(maximum(abs.(dlossdps[i_o]))) # dτ (dloss) is O(1) so dp = 1/max(dlossdp) # we go with 0.2 bc with the gaussian noise, that will go up to 0.5 ish either way and we don't wanna wiggle things too much as that de-trains the model. (then again we can go higher bc PCA narrows the effect to only a few neurons)
                # learning_rates[i_o] = FT(0.1) * (dlossdps[i_o] ./ (maximum(abs.(dlossdps[i_o])))) # dp is in gradient direction, max value O(1)
                # no /= n_outputs here, they can overlap so we just let cedmf figure it out, we want each contribution to reach peak on its own i think...
                else
                    learning_rates[i_o] = FT(1.0) # this is the learning rate (we're already in dp units, it's already in units of dp to move output by 1)
                end
            end
            dps[i_o] = -learning_rates[i_o] .* dlossdps[i_o] # scale the parameters by the learning rate, this is just a hyperparameter to control how much we wiggle the parameters, so we can adjust it later if needed, negative only matters really when we're doing descent, shouldn't matter too much othrewise since the weight will change it

            dps[i_o] .*= dp_weights[i_o] # scale the parameters by the output weight, this is just a hyperparameter to control how much we wiggle the parameters, so we can adjust it later if needed
            # @warn "pca_mean_vec[$i_o] = $(pca_mean_vec[i_o]); dlossdps[$i_o] = $(dlossdps[i_o]); pca_weights_here = $(pca_weights_here); pca_components[$i_o] = $(pca_components[i_o])"
        end

        # @warn "dlossdps = $(dlossdps); dps = $dps; " # pca_components = $(pca_components); pca_explained_variance = $(pca_explained_variance); learning_rates = $(learning_rates);" # @info doesn't seem to work on workers idk why...


        dp = reduce(+, dps) # sum the gradients over all outputs, this is the gradient of the loss with respect to the parameters, which is just the PCA noise scaled by the weights and summed over all components

    end



    neural_network_params = copy(namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network"]) # copy so no mutate this is the neural network params, which are just the PCA noise scaled by the weights and summed over all components

    # @warn "dp = $dp; neural_network_params = $neural_network_params"

    # essentially we're getting how much we'd have to wiggle each parameter to get a change in the otput, but we dont want to wiggle all the parameters.
    # the wiggle should look like 

    neural_network_params .+= dp
    neural_network_params = create_svector(neural_network_params[:]) # convert to SVector for NN since that's what it's supposed to be (save eltype in the jld2?)

    # neural_network_params = create_svector(FTNN.(pca_mean_vec .+ sum(pca_weights.* pca_components, dims=2) )[:]) # this is the PCA noise, so we add the mean and scale by the weights, and sum over all components

    # neural_network_params = create_svector(FTNN.(pca_mean_vec .+ sum( ((pca_weights .* inv.(pca_explained_variance)).* pca_components), dims=2) )[:]) # this is the PCA noise, so we add the mean and scale by the weights, and sum over all components
    # @warn "new Neural network params: $(neural_network_params)"
    flush(stdout)
    flush(stderr)

    # @warn "dp = $(dp); " # pca_weights = $(pca_weights); pca_components = $(pca_components); pca_explained_variance = $(pca_explained_variance)"
    flush(stdout)
    flush(stderr)
    # @error "dp = $(dp); old_neural_network_params = $(neural_network_params); pca_weights = $(pca_weights); pca_components = $(pca_components); pca_explained_variance = $(pca_explained_variance)"
    # flush(stdout); flush(stderr)
    # error("have i printed? if not why?")

    # rn everything in user_params ends up in param_set... but that isn't really necessary given we're creating these objects... it was more of a convenience when we were testing random things...
    model_x_0_characteristic = Tuple(FT.(namelist["relaxation_timescale_params"]["model_x_0_characteristic"])) # make this match the data, the type conversion will happen in prepare_for_NN()
    model_re_location = namelist["relaxation_timescale_params"]["model_re_location"]


    neural_network = simple_chain_model_from_file(model_re_location) # construct the NN from the parameters, and also get the simple chain params for the neural network

    adjust_liq_N = get(namelist["user_args"], "adjust_liq_N", false) # default to false if not set
    adjust_ice_N = get(namelist["user_args"], "adjust_ice_N", false) # default to false if not set # we could add this... hopefully the NN stays reasonable but if extrapolating to unseen data who knows...

    namelist["relaxation_timescale_params"]["neural_microphysics_relaxation_network_actual"] .= neural_network_params # update the namelist so the real params get written to namelist_SOCRATES.in

    return NeuralNetworkRelaxationTimescale(
        neural_network,
        neural_network_params,
        model_x_0_characteristic,
        adjust_liq_N,
        adjust_ice_N,
        get_relaxation_timescale_args(namelist, FT),
    )
end

# # :raymond_ice_test
# function get_relaxation_timescale_type(::Val{:raymond_ice_test}, param_set::APS, microphys_params::ACMP, namelist)
#     FT = eltype(param_set)
#     N_0 = get(namelist["user_args"], "N_0", FT(100 * 10e6)) # 100 cm^-3
#     N_r_closure = get(namelist["user_args"], "N_r_closure", :monodisperse) # default to monodisperse if not set
#     if N_r_closure === :monodisperse
#         N_r_closure = MonodisperseNRClosure()
#     elseif N_r_closure === :fixed_radius
#         N_r_closure = FixedRadiusNRClosure()
#     elseif N_r_closure === :inhomogeneous
#         N_r_closure = InhomogeneousNRClosure()
#     else
#         error("NotImplementedError: N_r_closure $N_r_closure not implemented")
#     end

#     return RaymondIceTestRelaxationTimescale( N_r_closure, N_0, get_relaxation_timescale_args(namelist, FT))
# end


# ======================================================================================================================================== #



# Timescales that predict INP and do not limit it based on q or anything else.
const INP_Aware_Timescale = Union{
    ExponentialTScalingIceRelaxationTimescale,
    ExponentialTScalingIceRawRelaxationTimescale,
    PowerlawTScalingIceRelaxationTimescale,
    GeometricLiqExponentialTScalingIceRelaxationTimescale,
    GeometricLiqPowerlawTScalingIceRelaxationTimescale,
    GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, # we made a special carveout fcn for this one.
    # The NN and Linear Combination take T into account but they tend to modify their N based on q, as per the training data, which doesn't include INP.
}
