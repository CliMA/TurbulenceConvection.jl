"""
    PrecipFormation

Storage for tendencies due to precipitation formation

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrecipFormation{FT}
    θ_liq_ice_tendency::FT
    qt_tendency::FT
    ql_tendency::FT
    qi_tendency::FT
    qr_tendency::FT
    qs_tendency::FT
    # added these below for use in storing for outputs
    ql_tendency_acnv::FT
    qi_tendency_acnv::FT
    ql_tendency_accr_liq_rai::FT
    ql_tendency_accr_liq_ice::FT
    ql_tendency_accr_liq_sno::FT
    qi_tendency_accr_ice_liq::FT
    qi_tendency_accr_ice_rai::FT
    qi_tendency_accr_ice_sno::FT
end

"""
    NoneqMoistureSources

Storage for tendencies due to nonequilibrium moisture formation

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct NoneqMoistureSources{FT}
    ql_tendency::FT
    qi_tendency::FT
end

"""
    NoneqMoistureSource

Storage for tendency in a condensate field

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct NoneqMoistureSource{FT}
    q_tendency::FT
end


"""
    Other Microphysics

Storage for tendencies due to other microphsyics moisture formation

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OtherMicrophysicsSources{FT}
    ql_tendency::FT
    qi_tendency::FT
    # added these for separate writing to outputs
    qi_tendency_homogeneous_freezing::FT
    qi_tendency_heterogeneous_freezing::FT
    qi_tendency_heterogeneous_icenuc::FT
    qi_tendency_melting::FT
end

"""
    EntrDetr

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct EntrDetr{FT}
    "Fractional dynamical entrainment [1/m]"
    ε_dyn::FT
    "Fractional dynamical detrainment [1/m]"
    δ_dyn::FT
    "Nondimensional fractional dynamical entrainment"
    ε_nondim::FT
    "Nondimensional fractional dynamical detrainment"
    δ_nondim::FT
end

Base.@kwdef struct εδModelParams{FT, AFT}
    w_min::FT # minimum updraft velocity to avoid zero division in b/w²
    c_ε::FT # factor multiplier for dry term in entrainment/detrainment
    μ_0::FT # dimensional scale logistic function in the dry term in entrainment/detrainment
    β::FT # sorting power for ad-hoc moisture detrainment function
    χ::FT # fraction of updraft air for buoyancy mixing in entrainment/detrainment (0≤χ≤1)
    c_λ::FT # scaling factor for TKE in entrainment scale calculations
    γ_lim::FT
    β_lim::FT
    limit_min_area::Bool
    min_area_limiter_scale::FT
    min_area_limiter_power::FT
    c_γ::FT # scaling factor for turbulent entrainment rate
    c_δ::FT # factor multiplier for moist term in entrainment/detrainment
    Π_norm::AFT
    entr_nondim_norm_factor::FT
    detr_nondim_norm_factor::FT

end

abstract type AbstractEntrDetrModel end
abstract type AbstractNoisyEntrDetrModel <: AbstractEntrDetrModel end
Base.@kwdef struct MDEntr{P} <: AbstractEntrDetrModel
    params::P
end  # existing model
Base.@kwdef struct EntrNone{P} <: AbstractEntrDetrModel
    params::P
end  # empthy model
Base.@kwdef struct NoisyRelaxationProcess{MT, T} <: AbstractNoisyEntrDetrModel
    mean_model::MT
    c_gen_stoch::T
end
Base.@kwdef struct LogNormalScalingProcess{MT, T} <: AbstractNoisyEntrDetrModel
    mean_model::MT
    c_gen_stoch::T
end

Base.@kwdef struct PrognosticNoisyRelaxationProcess{MT, T} <: AbstractNoisyEntrDetrModel
    mean_model::MT
    c_gen_stoch::T
end

abstract type AbstractMLEntrDetrModel end
abstract type AbstractMLNonLocalEntrDetrModel <: AbstractMLEntrDetrModel end
Base.@kwdef struct NNEntr{P, AFT, T} <: AbstractMLEntrDetrModel
    params::P
    c_nn_params::AFT
    nn_arc::T
    biases_bool::Bool
end
Base.@kwdef struct NNEntrNonlocal{P, AFT, T} <: AbstractMLNonLocalEntrDetrModel
    params::P
    c_nn_params::AFT
    nn_arc::T
    biases_bool::Bool
end
Base.@kwdef struct LinearEntr{P, T} <: AbstractMLEntrDetrModel
    params::P
    c_linear::T
    biases_bool::Bool
end
Base.@kwdef struct FNOEntr{P, T} <: AbstractMLNonLocalEntrDetrModel
    params::P
    w_fno::Int
    nm_fno::Int
    c_fno::T
end
struct RFEntr{d, m, FT, P, A, B} <: AbstractMLEntrDetrModel
    params::P
    c_rf_fix::A
    c_rf_opt::B
    function RFEntr(params::P, c_rf_fix::A, c_rf_opt::B, d::Int) where {P, A, B}
        c_rf_fix = reshape(c_rf_fix, 2, :, 1 + d) # 2 x m x (1 + d), fix
        c_rf_fix = SA.SArray{Tuple{size(c_rf_fix)...}}(c_rf_fix)
        FT = eltype(c_rf_fix)
        c_rf_fix = Array{FT}(c_rf_fix)
        m = size(c_rf_fix, 2)
        c_rf_opt = reshape(c_rf_opt, 2, m + 1 + d) # 2 x (m + 1 + d), learn
        c_rf_opt = SA.SArray{Tuple{size(c_rf_opt)...}}(c_rf_opt)
        c_rf_opt = Array{FT}(c_rf_opt)
        return new{d, m, FT, P, typeof(c_rf_fix), typeof(c_rf_opt)}(params, c_rf_fix, c_rf_opt)
    end
end

Base.eltype(::RFEntr{d, m, FT}) where {d, m, FT} = FT

εδ_params(m::AbstractEntrDetrModel) = m.params
εδ_params(m::AbstractMLEntrDetrModel) = m.params
εδ_params(m::AbstractNoisyEntrDetrModel) = m.mean_model.params

abstract type EntrModelFacTotalType end
struct FractionalEntrModel <: EntrModelFacTotalType end
struct TotalRateEntrModel <: EntrModelFacTotalType end

abstract type EntrDimScale end
struct BuoyVelEntrDimScale <: EntrDimScale end
struct InvScaleHeightEntrDimScale <: EntrDimScale end
struct InvZEntrDimScale <: EntrDimScale end
struct InvMeterEntrDimScale <: EntrDimScale end
struct PosMassFluxGradDimScale <: EntrDimScale end
struct NegMassFluxGradDimScale <: EntrDimScale end
struct AbsMassFluxGradDimScale <: EntrDimScale end
struct MassFluxGradDimScale <: EntrDimScale end
struct BOverWDimScale <: EntrDimScale end
struct WOverHeightDimScale <: EntrDimScale end
struct BOverSqrtTKEDimScale <: EntrDimScale end
struct SqrtBOverZDimScale <: EntrDimScale end
struct TKEBWDimScale <: EntrDimScale end
struct DwDzDimScale <: EntrDimScale end
struct DmDzOverRhoaDimScale <: EntrDimScale end


"""
    GradBuoy

Environmental buoyancy gradients.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct GradBuoy{FT}
    "environmental vertical buoyancy gradient"
    ∂b∂z::FT
    "vertical buoyancy gradient in the unsaturated part of the environment"
    ∂b∂z_unsat::FT
    "vertical buoyancy gradient in the saturated part of the environment"
    ∂b∂z_sat::FT
end

abstract type AbstractEnvBuoyGradClosure end
struct BuoyGradMean <: AbstractEnvBuoyGradClosure end
struct BuoyGradQuadratures <: AbstractEnvBuoyGradClosure end

"""
    EnvBuoyGrad

Variables used in the environmental buoyancy gradient computation.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct EnvBuoyGrad{FT, EBC <: AbstractEnvBuoyGradClosure}
    "temperature in the saturated part"
    t_sat::FT
    "vapor specific humidity  in the saturated part"
    qv_sat::FT
    "total specific humidity in the saturated part"
    qt_sat::FT
    "potential temperature in the saturated part"
    θ_sat::FT
    "liquid ice potential temperature in the saturated part"
    θ_liq_ice_sat::FT
    "virtual potential temperature gradient in the non saturated part"
    ∂θv∂z_unsat::FT
    "total specific humidity gradient in the saturated part"
    ∂qt∂z_sat::FT
    "liquid ice potential temperature gradient in the saturated part"
    ∂θl∂z_sat::FT
    "reference pressure"
    p::FT
    "cloud fraction"
    en_cld_frac::FT
    "density"
    ρ::FT
end
function EnvBuoyGrad(::EBG; t_sat::FT, bg_kwargs...) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}
    return EnvBuoyGrad{FT, EBG}(; t_sat, bg_kwargs...)
end

"""
    MixLen

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MixLen{FT}
    "minimum length number"
    min_len_ind::Int
    "mixing length"
    mixing_length::FT
    "length ratio"
    ml_ratio::FT
end

Base.@kwdef struct MixingLengthParams{FT}
    ω_pr::FT # cospectral budget factor for turbulent Prandtl number
    c_m::FT # tke diffusivity coefficient
    c_d::FT # tke dissipation coefficient
    c_b::FT # static stability coefficient
    κ_star²::FT # Ratio of TKE to squared friction velocity in surface layer
    Pr_n::FT # turbulent Prandtl number in neutral conditions
    Ri_c::FT # critical Richardson number
    Le::FT # Lewis number
    smin_ub::FT # lower limit for smin function
    l_max::FT
end

"""
    MinDisspLen

Minimum dissipation model

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MinDisspLen{FT}
    "height"
    z::FT
    "obukhov length"
    obukhov_length::FT
    "surface TKE values"
    tke_surf::FT
    "u star - surface velocity scale"
    ustar::FT
    "turbulent Prandtl number"
    Pr::FT
    "reference pressure"
    p::FT
    "vertical buoyancy gradient struct"
    ∇b::GradBuoy{FT}
    "env shear"
    Shear²::FT
    "environment turbulent kinetic energy"
    tke::FT
    "Updraft tke source"
    b_exch::FT
end

abstract type SurfaceAreaBC end
struct FixedSurfaceAreaBC <: SurfaceAreaBC end
struct PrognosticSurfaceAreaBC <: SurfaceAreaBC end
Base.@kwdef struct ClosureSurfaceAreaBC{P} <: SurfaceAreaBC
    params::P
end


Base.@kwdef struct PressureModelParams{FT}
    α_b::FT # factor multiplier for pressure buoyancy terms (effective buoyancy is (1-α_b))
    α_a::FT # factor multiplier for pressure advection
    α_d::FT # factor multiplier for pressure drag
end

abstract type AbstractMoistureModel end
struct EquilibriumMoisture <: AbstractMoistureModel end
# struct NonEquilibriumMoisture <: AbstractMoistureModel end

"""
NR Closures (assumptions about the size distribution) -- used in relaxation timescale types so make sure it's before
"""
abstract type AbstractNRClosureType end
struct MonodisperseNRClosure <: AbstractNRClosureType end
struct FixedRadiusNRClosure <: AbstractNRClosureType end
struct InhomogeneousNRClosure <: AbstractNRClosureType end


include("relaxation_timescale_types.jl")

struct NonEquilibriumMoisture{RTT <: AbstractNonEquillibriumSourcesType, FT <: AbstractFloat} <: AbstractMoistureModel
    scheme::RTT
    heterogeneous_ice_nucleation::Tuple{Bool, FT, FT}
end

function NonEquilibriumMoisture(param_set::APS)
    FT = eltype(param_set)
    nonequilibrium_moisture_scheme_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :nonequilibrium_moisture_scheme, :relax_to_equilibrium)

    nonequilibrium_moisture_scheme = if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
        RelaxToEquilibrium()
    elseif nonequilibrium_moisture_scheme_type === :KorolevMazin2007
        KorolevMazin2007()
    elseif nonequilibrium_moisture_scheme_type ∈ valid_relaxation_timescale_types
        get_relaxation_timescale_type(nonequilibrium_moisture_scheme_type, param_set)
    else
        error("Invalid nonequilibrium_moisture_scheme type: $nonequilibrium_moisture_scheme_type")
    end


    heterogeneous_ice_nucleation = TCP.parse_isbits_nt(param_set.user_args, :use_heterogeneous_ice_nucleation, false)
    heterogeneous_ice_nuclation_coefficient = TCP.parse_isbits_nt(param_set.user_params, :heterogeneous_ice_nuclation_coefficient, FT(1))
    heterogeneous_ice_nuclation_exponent = TCP.parse_isbits_nt(param_set.user_params, :heterogeneous_ice_nuclation_exponent, FT(1))


    RTT = typeof(nonequilibrium_moisture_scheme)
    return NonEquilibriumMoisture{RTT, FT}(nonequilibrium_moisture_scheme, (heterogeneous_ice_nucleation, heterogeneous_ice_nuclation_coefficient, heterogeneous_ice_nuclation_exponent), )
end

abstract type AbstractIntegrationScheme end
struct UpwindDifferencingScheme <: AbstractIntegrationScheme end
struct RightBiasedDifferencingScheme <: AbstractIntegrationScheme end
abstract type AbstractCloudSedimentationModel end
struct CloudSedimentationModel{FT <: AbstractFloat, LTVS <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, ITVS <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, SLNC <: Union{AbstractFloat, AbstractRelaxationTimescaleType}, SINC <: Union{AbstractFloat, AbstractRelaxationTimescaleType}, SIM <: AbstractIntegrationScheme} <: AbstractCloudSedimentationModel
    liq_terminal_velocity_scheme::LTVS
    ice_terminal_velocity_scheme::ITVS
    sedimentation_liq_number_concentration::SLNC
    sedimentation_ice_number_concentration::SINC
    liq_Dmax::FT
    ice_Dmax::FT
    liq_sedimentation_scaling_factor::FT
    ice_sedimentation_scaling_factor::FT
    sedimentation_differencing_scheme::SIM
    liq_ice_collision_efficiency::FT 
    liq_ice_collision_scaling_factor::FT # this should be deprecated one day but I'm using it for initial calibration just to check if we're even close
    grid_mean::Bool # whether or not sedimentation is applied to the grid mean only instead of in env/up separately
end

struct CloudNoSedimentationModel <: AbstractCloudSedimentationModel end

function CloudSedimentationModel(param_set)
    FT = eltype(param_set)
    liq_terminal_velocity_scheme = get_termvel_type(TCP.parse_isbits_nt(param_set.user_args, :liq_terminal_velocity_scheme, :Blk1MVel)) 
    ice_terminal_velocity_scheme = get_termvel_type(TCP.parse_isbits_nt(param_set.user_args, :ice_terminal_velocity_scheme, :Blk1MVel))

    sedimentation_liq_number_concentration = TCP.parse_isbits_nt(param_set.user_args, :sedimentation_liq_number_concentration, FT(NaN)) # testing NaN over nothing for type stability
    sedimentation_ice_number_concentration = TCP.parse_isbits_nt(param_set.user_args, :sedimentation_ice_number_concentration, FT(NaN))
    if sedimentation_liq_number_concentration isa Symbol
        sedimentation_liq_number_concentration = get_relaxation_timescale_type(sedimentation_liq_number_concentration, param_set) # preserve isbits
    end
    if sedimentation_ice_number_concentration isa Symbol
        sedimentation_ice_number_concentration = get_relaxation_timescale_type(sedimentation_ice_number_concentration, param_set) # preserve isbits
    end

    liq_Dmax = TCP.parse_isbits_nt(param_set.user_params, :liq_sedimentation_Dmax, FT(Inf))
    ice_Dmax = TCP.parse_isbits_nt(param_set.user_params, :ice_sedimentation_Dmax, FT(62.5e-6)) # maybe this should also be inf? for chen it's not clear...
   
    liq_sedimentation_scaling_factor = TCP.parse_isbits_nt(param_set.user_params, :liq_sedimentation_scaling_factor, FT(1.0))
    ice_sedimentation_scaling_factor = TCP.parse_isbits_nt(param_set.user_params, :ice_sedimentation_scaling_factor, FT(1.0))

    sedimentation_differencing_scheme = TCP.parse_isbits_nt(param_set.user_params, :sedimentation_differencing_scheme, :upwinding)
     if sedimentation_differencing_scheme === :upwinding
        sedimentation_differencing_scheme = UpwindDifferencingScheme()
    elseif sedimentation_differencing_scheme === :right_biased
        sedimentation_differencing_scheme = RightBiasedDifferencingScheme()
    else
        error("Invalid sedimentation_differencing_scheme: $sedimentation_differencing_scheme")
    end

    liq_ice_collision_efficiency = TCP.parse_isbits_nt(param_set.user_params, :liq_ice_collision_efficiency, FT(1.0))
    liq_ice_collision_scaling_factor = TCP.parse_isbits_nt(param_set.user_params, :liq_ice_collision_scaling_factor, FT(1.0))

    LTVS = typeof(liq_terminal_velocity_scheme)
    ITVS = typeof(ice_terminal_velocity_scheme)
    SLNC = typeof(sedimentation_liq_number_concentration)
    SINC = typeof(sedimentation_ice_number_concentration)
    SDS = typeof(sedimentation_differencing_scheme)

    grid_mean = TCP.parse_isbits_nt(param_set.user_args, :grid_mean_sedimentation, false) # only apply to grid mean values, never actually got it to work though


    return CloudSedimentationModel{FT, LTVS, ITVS, SLNC, SINC, SDS}(
        liq_terminal_velocity_scheme, 
        ice_terminal_velocity_scheme,
        sedimentation_liq_number_concentration,
        sedimentation_ice_number_concentration,
        liq_Dmax,
        ice_Dmax,
        liq_sedimentation_scaling_factor,
        ice_sedimentation_scaling_factor,
        sedimentation_differencing_scheme,
        liq_ice_collision_efficiency,
        liq_ice_collision_scaling_factor,
        grid_mean,
    )
end

abstract type AbstractCovarianceModel end
struct PrognosticThermoCovariances <: AbstractCovarianceModel end
struct DiagnosticThermoCovariances{FT} <: AbstractCovarianceModel
    covar_lim::FT
end

abstract type AbstractPrecipitationModel end
struct NoPrecipitation <: AbstractPrecipitationModel end
struct Clima0M <: AbstractPrecipitationModel end
struct Clima1M{FT, RTVST <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, STVST <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}} <: AbstractPrecipitationModel 
    rain_sedimentation_scaling_factor::FT
    snow_sedimentation_scaling_factor::FT
    rain_terminal_velocity_scheme::RTVST
    snow_terminal_velocity_scheme::STVST
end
function Clima1M(param_set)
    FT = eltype(param_set)
    rain_sedimentation_scaling_factor = TCP.parse_isbits_nt(param_set.user_params, :rain_sedimentation_scaling_factor, FT(1.0))
    snow_sedimentation_scaling_factor = TCP.parse_isbits_nt(param_set.user_params, :snow_sedimentation_scaling_factor, FT(1.0))
    rain_terminal_velocity_scheme = get_termvel_type(TCP.parse_isbits_nt(param_set.user_args, :rain_terminal_velocity_scheme, :Blk1MVel))
    snow_terminal_velocity_scheme = get_termvel_type(TCP.parse_isbits_nt(param_set.user_args, :snow_terminal_velocity_scheme, :Blk1MVel))
    RTVST = typeof(rain_terminal_velocity_scheme)
    STVST = typeof(snow_terminal_velocity_scheme)
    return Clima1M{FT, RTVST, STVST}(rain_sedimentation_scaling_factor, snow_sedimentation_scaling_factor, rain_terminal_velocity_scheme, snow_terminal_velocity_scheme)
end




"""
   AbstractRainFormationModel

A type to switch between different autoconversion and accretion options.
See the CloudMicrophysics.jl documentation for the Clima1M_default version
and the additional options with prescribed cloud droplet number concentration.
"""
abstract type AbstractRainFormationModel end
struct NoRainFormation <: AbstractRainFormationModel end
struct Clima1M_default <: AbstractRainFormationModel end
struct Clima2M{FT, PT} <: AbstractRainFormationModel
    prescribed_Nd::FT
    type::PT
end

abstract type AbstractSnowFormationModel end
struct DefaultSnowFormationModel <: AbstractSnowFormationModel end
struct NonEquilibriumSnowFormationModel <: AbstractSnowFormationModel
    # scheme::RTT # currently we're not using a relaxation_timescale_type for snow formation. If we do, we can bring this back but you should reuse the same one as for the moisture model and not recreate it.
end
function SnowFormationModel(param_set)

    nonequilibrium_moisture_scheme_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :nonequilibrium_moisture_scheme, :relax_to_equilibrium)

    snow_formation_model = if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
        DefaultSnowFormationModel()
    elseif nonequilibrium_moisture_scheme_type === :KorolevMazin2007
        DefaultSnowFormationModel() # should this be default? I guess maybe cause there's no N prediction?
    elseif nonequilibrium_moisture_scheme_type ∈ valid_relaxation_timescale_types
        NonEquilibriumSnowFormationModel()
    else
        error("Invalid nonequilibrium_moisture_scheme type: $nonequilibrium_moisture_scheme_type")
    end

    return snow_formation_model
end

abstract type AbstractPrecipFractionModel end
struct PrescribedPrecipFraction{FT} <: AbstractPrecipFractionModel
    prescribed_precip_frac_value::FT
end
struct DiagnosticPrecipFraction{FT} <: AbstractPrecipFractionModel
    precip_fraction_limiter::FT
end









"""
Making how we handle limiters consistent across tendencies
"""

# Default Tendency Limiter
abstract type AbstractLimiter end
Base.broadcastable(x::AbstractLimiter) = Ref(x) # permit broadcasting over these types, 
abstract type AbstractTendencyLimiter <: AbstractLimiter end
struct NoTendencyLimiter <: AbstractTendencyLimiter end
limit_tendency(::NoTendencyLimiter, x::FT, x_up::FT, Δt::FT) where {FT} = x
limit_tendency(::NoTendencyLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = x
struct BasicTendencyLimiter <: AbstractTendencyLimiter end # Truncate to the available value within a timestep
limit_tendency(::BasicTendencyLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = max(x, -x_avail/Δt) #
limit_tendency(::BasicTendencyLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = safe_clamp(x, -x_up/Δt, x_en/Δt)

# Maybe consider adding
# abstract type AbstractdtDependentTendencyLimiter <: AbstractTendencyLimiter end that can supertype all dt dependent tendencies... rn using !isa(NoTendencyLimiter, limiter) is a bit hacky
# could also try to find a way to unify NoTendency and NoMoistureSourcesLimiter ? idk.
# abstract type N

# Moisture Limiters
abstract type AbstractMoistureSourcesLimiter <: AbstractTendencyLimiter end
struct NoMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end # No limiting is done, the sources are just directly calculated
struct BasicMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end # Sources are truncated to their available value within a timestep
struct StandardSupersaturationMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end # Sources are truncated to their available value within a timestep but this is geared towards supersaturation relaxation. It does a lot of checks that are unnecessary if you're just relaxing to equilibrium e.g. WBF regimes etc.
struct MorrisonMilbrandt2015MoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end
struct MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end
# Fallbacks (in noneq_moisture sources we have our own setup, but for dispatch these may be needed elsewhere)
limit_tendency(::AbstractMoistureSourcesLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = max(x, -x_avail/Δt)
limit_tendency(::AbstractMoistureSourcesLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = safe_clamp(x, -x_up/Δt, x_en/Δt) 
limit_tendency(::NoMoistureSourcesLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = x
limit_tendency(::NoMoistureSourcesLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = x

# # Small Tendency Limiter
# abstract type AbstractSmallTendencyLimiter <: AbstractTendencyLimiter end
# struct NoSmallTendencyLimiter <: AbstractSmallTendencyLimiter end
# struct BasicSmallTendencyLimiter{FT} <: AbstractSmallTendencyLimiter
#     cutoff::FT
# end
# limit_tendency(::NoSmallTendencyLimiter, x::FT, cutoff::FT) where {FT} = x
# limit_tendency(::BasicSmallTendencyLimiter, x::FT, cutoff::FT) where {FT} = abs(x) < cutoff ? zero(FT) : x

# Entr/Detr Limiters -- use default types
# Precipitation Limiters (for precip formation/accretion + for evaporation/melt/etc) -- use default types

abstract type AbstractLimiterSet end
abstract type AbstractTendencyLimiterSet <: AbstractLimiterSet end

"""
Be cautious selecting your tendency limiters... if dt is adapive, you should try to not let it vary very much because the tendencies mix under limiting will be timestep dependent as tendencies take differing amounts of times to saturate.
e.g. the magnitude of an originally fast tendency relative to a slow tendency will increase with shorter timesteps as that tendency gets less and less limited.

For best results, maybe use NoTendencyLimiter for highly adaptive timesteps if the timestep is calculated from the tendencies (though this could slow to a crawl in theory)
For limiters like basic and the noneq ones, maybe only stick to consistent timesteps, e.g. within a factor of 2-5 or a smooth spinup to regular dt transition, but you may still note some jittering.


"""
struct TendencyLimiterSet{DTLT <: AbstractTendencyLimiter, MSLT <: AbstractMoistureSourcesLimiter, EDTLT <: AbstractTendencyLimiter, PTLT <: AbstractTendencyLimiter, FDTLT <: AbstractTendencyLimiter, FMSLT <: AbstractMoistureSourcesLimiter, FEDTLT <: AbstractTendencyLimiter, FPTLT <: AbstractTendencyLimiter} <: AbstractTendencyLimiterSet
    default_tendency_limiter::DTLT #  for all other tendencies (gm, up etc) and as a
    moisture_sources_limiter::MSLT # if we're not a nonequilibrium moisture model, there's no limiting...can just use NoMoistureSourcesLimiter
    entr_detr_tendency_limiter::EDTLT
    precipitation_tendency_limiter::PTLT # precipitation formation, accretion, evaporation, as well as melt, sublimation, etc.

    # should we add a ∑tendencies_limiter? Then everything could be NoTendencyLimiter and you could limit after the summation rather than using tendency limited dt... (this is kind like filtering and induces leaks/drifts)
    # rn we just use default_tendency_limiter for that but that could become confusing if we ever decide to use default_tendency_limiter for something else
    # up_en_tendency_limiter::UETL # this allows the sum of all the updraft tendencies to be limited separately from the gm that uses default_tendency_limiter. These are most important bc theyre composed. However, decomposing this change into changes that do and don't affect the grid-mean is impossible.
    # advection/sedimentation are self limiting by CFL enforcement

    # It would be nice to have a way to transition to a fallback (e.g. BasicLimiter) once we run out of N_dt_max_edmf_violate_dt_min_remaining for stability (or maybe just letting it crash is fine idk... but that's how it was originally written)
    # We could make the TendencyLimiterSet object mutable but that sounds harsh... ... just something to think about for now -- rn if you run out of N_dt_max_edmf_violate_dt_min_remaining using NoTendencyLimiters, you'll almost certainly have a crash which wasn't always true but is in principle more accurate
    fallback_default_tendency_limiter::FDTLT
    fallback_moisture_sources_limiter::FMSLT
    fallback_entr_detr_tendency_limiter::FEDTLT
    fallback_precipitation_tendency_limiter::FPTLT
end

get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:default}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_default_tendency_limiter : tendency_limiter_set.default_tendency_limiter
get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:moisture_sources}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_moisture_sources_limiter : tendency_limiter_set.moisture_sources_limiter
get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:entr_detr}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_entr_detr_tendency_limiter : tendency_limiter_set.entr_detr_tendency_limiter
get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:precipitation}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_precipitation_tendency_limiter : tendency_limiter_set.precipitation_tendency_limiter

function TendencyLimiterSet(param_set::APS, moisture_model_name::String)
    FT = eltype(param_set)

    # ------------------ Default Tendency Limiter ------------------ #
    default_tendency_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :default_tendency_limiter_type, :none)  # Default to none bc that's purest (basic is how anna wrote it though)

    default_tendency_limiter = if default_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif default_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    else
        error("Invalid default_tendency_limiter_type: $default_tendency_limiter_type, valid options are :none, :basic")
    end

    fallback_default_tendency_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :fallback_default_tendency_limiter_type, :none)  # Default to none bc that's purest (basic is how anna wrote it though)
    fallback_default_tendency_limiter = if fallback_default_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif fallback_default_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    else
        error("Invalid fallback_default_tendency_limiter_type: $fallback_default_tendency_limiter_type, valid options are :none, :basic")
    end

    # ------------------ Moisture Sources Limiter ------------------ #
    if moisture_model_name == "equilibrium" 
        # can we change this to create no object at all? Or not, idk... maybe it's not that wasteful and improves type inference/stability.
        moisture_sources_limiter = NoMoistureSourcesLimiter() # if we're not a nonequilibrium moisture model, there's no limiting just sat adjust (liq/ice are diagnosed)...can just use NoMoistureSourcesLimiter
        fallback_moisture_sources_limiter = NoMoistureSourcesLimiter()
    elseif moisture_model_name == "nonequilibrium"    

        nonequilibrium_moisture_scheme_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :nonequilibrium_moisture_scheme, :relax_to_equilibrium)
        nonequilibrium_moisture_sources_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :nonequilibrium_moisture_sources_limiter_type, :standard)  # this has w built in though, no way around it, maybe we should write one that's just the exponential decay part w/o anything else... ( # change the default to at least :morrison_milbrandt_2015_style_exponential_part_only soon )
        if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
            @assert nonequilibrium_moisture_sources_limiter_type ∈ [:None, :standard,] "Relaxation to equilibrium only supports no integrator limiter or basic integrator limiter"
        end
        

        moisture_sources_limiter = if nonequilibrium_moisture_sources_limiter_type === :none
            NoMoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :standard_supersaturation
            StandardSupersaturationMoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :basic
            BasicMoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style
            MorrisonMilbrandt2015MoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style_exponential_part_only
            MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter()
        else
            error("Invalid nonequilibrium_moisture_sources_limiter_type: $nonequilibrium_moisture_sources_limiter_type, valid options are :none, :standard, :morrison_milbrandt_2015_style, :morrison_milbrandt_2015_style_exponential_part_only")
        end

        fallback_nonequilibrium_moisture_sources_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :fallback_nonequilibrium_moisture_sources_limiter_type, :none)  # Default to none bc that's purest (basic is how anna wrote it though)
        if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
            @assert fallback_nonequilibrium_moisture_sources_limiter_type ∈ [:None, :standard,] "Relaxation to equilibrium only supports no integrator limiter or basic integrator limiter"
        end

        fallback_moisture_sources_limiter = if fallback_nonequilibrium_moisture_sources_limiter_type === :none
            NoMoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :basic
            BasicMoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :standard_supersaturation
            StandardSupersaturationMoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style
            MorrisonMilbrandt2015MoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style_exponential_part_only
            MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter()
        else
            error("Invalid fallback_nonequilibrium_moisture_sources_limiter_type: $fallback_nonequilibrium_moisture_sources_limiter_type, valid options are :none, :basic, :standard_supersaturation, :morrison_milbrandt_2015_style, :morrison_milbrandt_2015_style_exponential_part_only")
        end

    else
        error("Something went wrong. Invalid moisture model: '$moisture_model_name'")
    end
    

    # ------------------ Entrainment/Detrainment Limiter ------------------ #
    entr_detr_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :entr_detr_limiter_type, :none)  # note entr/detr as written are never depndent on env, so the values here can be truly uninged w/ no limiting. However, that's the purest form for composability with other tendencies.
    entr_detr_tendency_limiter = if entr_detr_limiter_type === :none
        NoTendencyLimiter()
    elseif entr_detr_limiter_type === :basic
        BasicTendencyLimiter()
    else
        error("Invalid entr_detr_limiter_type: $entr_detr_limiter_type, valid options are :none, :basic")
    end

    fallback_entr_detr_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :fallback_entr_detr_limiter_type, :none)  # Default to none bc that's purest (basic is how anna wrote it though)
    fallback_entr_detr_tendency_limiter = if fallback_entr_detr_limiter_type === :none
        NoTendencyLimiter()
    elseif fallback_entr_detr_limiter_type === :basic
        BasicTendencyLimiter()
    else
        error("Invalid fallback_entr_detr_limiter_type: $fallback_entr_detr_limiter_type, valid options are :none, :basic")
    end
    
    # ------------------ Precipitation Limiter ------------------ #
    precipitation_tendency_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :precipitation_tendency_limiter_type, :none)  # Default to none bc that's purest (basic is how anna wrote it though)
    precipitation_tendency_limiter = if precipitation_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif precipitation_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    else
        error("Invalid precipitation_tendency_limiter_type: $precipitation_tendency_limiter_type, valid options are :none, :basic")
    end

    fallback_precipitation_tendency_limiter_type::Symbol = TCP.parse_isbits_nt(param_set.user_args, :fallback_precipitation_tendency_limiter_type, :none)  # Default to none bc that's purest (basic is how anna wrote it though)
    fallback_precipitation_tendency_limiter = if fallback_precipitation_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif fallback_precipitation_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    else
        error("Invalid fallback_precipitation_tendency_limiter_type: $fallback_precipitation_tendency_limiter_type, valid options are :none, :basic")
    end

    # ------------------ Return ------------------ #
 
    DTLT = typeof(default_tendency_limiter)
    MSLT = typeof(moisture_sources_limiter)
    EDTLT = typeof(entr_detr_tendency_limiter)
    PTLT = typeof(precipitation_tendency_limiter)

    FDTLT = typeof(fallback_default_tendency_limiter)
    FMSLT = typeof(fallback_moisture_sources_limiter)
    FEDTLT = typeof(fallback_entr_detr_tendency_limiter)
    FPTLT = typeof(fallback_precipitation_tendency_limiter)

    return TendencyLimiterSet{DTLT, MSLT, EDTLT, PTLT, FDTLT, FMSLT, FEDTLT, FPTLT}(
        default_tendency_limiter,
        moisture_sources_limiter,
        entr_detr_tendency_limiter,
        precipitation_tendency_limiter,
        #
        fallback_default_tendency_limiter,
        fallback_moisture_sources_limiter,
        fallback_entr_detr_tendency_limiter,
        fallback_precipitation_tendency_limiter,
    )

end



# We don't put these in a set right now because they're use is so diverse, insteady let's deprecate and just write functions instead
# abstract type AbstractValueLimiter <: AbstractLimiter end
# struct NoValueLimiter <: AbstractValueLimiter end
# struct BasicValueLimiter{FT} <: AbstractValueLimiter 
#     cutoff::FT
# end
# struct ZeroOrPosLimiter <: AbstractValueLimiter end
# struct ZeroOrNegLimiter <: AbstractValueLimiter end
# struct MinValueLimiter{FT} <: AbstractValueLimiter
#     min_value::FT
# end
# struct MaxValueLimiter{FT} <: AbstractValueLimiter
#     max_value::FT
# end
# limit_value(::NoValueLimiter, x::FT, x_avail::FT) where {FT} = x
# limit_value(::BasicValueLimiter, x::FT, cutoff::FT) where {FT} = abs(x) < cutoff ? zero(FT) : x
# limit_value(::ZeroOrPosLimiter, x::FT) where {FT} = max(x, zero(FT))
# limit_value(::ZeroOrNegLimiter, x::FT) where {FT} = min(x, zero(FT))
# limit_value(::MinValueLimiter, x::FT, min_value::FT) where {FT} = max(x, min_value)
# limit_value(::MaxValueLimiter, x::FT, max_value::FT) where {FT} = min(x, max_value)
cutoff_small_values(x::FT, cutoff::FT) where {FT} = abs(x) < cutoff ? zero(FT) : x
cutoff_small_values_positive(x::FT, cutoff::FT) where {FT} =  x < cutoff ? zero(FT) : x
cutoff_small_values_negative(x::FT, cutoff::FT) where {FT} =  x > -cutoff ? zero(FT) : x



abstract type AbstractQuadratureType end
struct LogNormalQuad <: AbstractQuadratureType end
struct GaussianQuad <: AbstractQuadratureType end

abstract type AbstractEnvThermo end
struct SGSMean <: AbstractEnvThermo end
struct SGSQuadrature{N, QT, A, W} <: AbstractEnvThermo
    quadrature_type::QT
    a::A
    w::W
    function SGSQuadrature(::Type{FT}, namelist) where {FT}
        N = parse_namelist(namelist, "thermodynamics", "quadrature_order"; default = 3)
        quadrature_name = parse_namelist(namelist, "thermodynamics", "quadrature_type"; default = "log-normal")
        quadrature_type = if quadrature_name == "log-normal"
            LogNormalQuad()
        elseif quadrature_name == "gaussian"
            GaussianQuad()
        else
            error("Invalid thermodynamics quadrature $(quadrature_name)")
        end
        # TODO: double check this python-> julia translation
        # a, w = np.polynomial.hermite.hermgauss(N)
        a, w = FastGaussQuadrature.gausshermite(N)
        a, w = SA.SVector{N, FT}(a), SA.SVector{N, FT}(w)
        QT = typeof(quadrature_type)
        return new{N, QT, typeof(a), typeof(w)}(quadrature_type, a, w)
    end
end
quadrature_order(::SGSQuadrature{N}) where {N} = N
quad_type(::SGSQuadrature{N}) where {N} = N

abstract type FrictionVelocityType end
struct FixedFrictionVelocity <: FrictionVelocityType end
struct VariableFrictionVelocity <: FrictionVelocityType end

abstract type AbstractSurfaceParameters{FT <: Real} end

const FloatOrFunc{FT} = Union{FT, Function, Dierckx.Spline1D}

Base.@kwdef struct FixedSurfaceFlux{FT, FVT <: FrictionVelocityType, TS, QS, SHF, LHF} <: AbstractSurfaceParameters{FT}
    zrough::FT = FT(0)
    Tsurface::TS = FT(0)
    qsurface::QS = FT(0)
    shf::SHF = FT(0)
    lhf::LHF = FT(0)
    cq::FT = FT(0)
    Ri_bulk_crit::FT = FT(0)
    ustar::FT = FT(0)
    zero_uv_fluxes::Bool = false
end

function FixedSurfaceFlux(
    ::Type{FT},
    ::Type{FVT};
    Tsurface::FloatOrFunc{FT},
    qsurface::FloatOrFunc{FT},
    shf::FloatOrFunc{FT},
    lhf::FloatOrFunc{FT},
    kwargs...,
) where {FT, FVT}
    TS = typeof(Tsurface)
    QS = typeof(qsurface)
    SHF = typeof(shf)
    LHF = typeof(lhf)
    return FixedSurfaceFlux{FT, FVT, TS, QS, SHF, LHF}(; Tsurface, qsurface, shf, lhf, kwargs...)
end

Base.@kwdef struct FixedSurfaceCoeffs{FT, TS, QS, CH, CM} <: AbstractSurfaceParameters{FT}
    zrough::FT = FT(0)
    Tsurface::TS = FT(0)
    qsurface::QS = FT(0)
    ch::CH = FT(0)
    cm::CM = FT(0)
    Ri_bulk_crit::FT = FT(0)
end

function FixedSurfaceCoeffs(
    ::Type{FT};
    Tsurface::FloatOrFunc{FT},
    qsurface::FloatOrFunc{FT},
    ch::FloatOrFunc{FT},
    cm::FloatOrFunc{FT},
    kwargs...,
) where {FT}
    TS = typeof(Tsurface)
    QS = typeof(qsurface)
    CH = typeof(ch)
    CM = typeof(cm)
    return FixedSurfaceCoeffs{FT, TS, QS, CH, CM}(; Tsurface, qsurface, ch, cm, kwargs...)
end

Base.@kwdef struct MoninObukhovSurface{FT, TS, QS} <: AbstractSurfaceParameters{FT}
    zrough::FT = FT(0)
    Tsurface::TS = FT(0)
    qsurface::QS = FT(0)
    Ri_bulk_crit::FT = FT(0)
end

function MoninObukhovSurface(
    ::Type{FT};
    Tsurface::FloatOrFunc{FT},
    qsurface::FloatOrFunc{FT},
    kwargs...,
) where {FT}
    TS = typeof(Tsurface)
    QS = typeof(qsurface)
    return MoninObukhovSurface{FT, TS, QS}(; Tsurface, qsurface, kwargs...)
end

float_or_func(s::Function, t::Real) = s(t)
float_or_func(s::Dierckx.Spline1D, t::Real) = s(t)
float_or_func(s::Real, t::Real) = s

surface_temperature(s::AbstractSurfaceParameters, t::Real = 0) = float_or_func(s.Tsurface, t)
surface_q_tot(s::AbstractSurfaceParameters, t::Real = 0) = float_or_func(s.qsurface, t)
sensible_heat_flux(s::AbstractSurfaceParameters, t::Real = 0) = float_or_func(s.shf, t)
latent_heat_flux(s::AbstractSurfaceParameters, t::Real = 0) = float_or_func(s.lhf, t)

fixed_ustar(::FixedSurfaceFlux{FT, FixedFrictionVelocity}) where {FT} = true
fixed_ustar(::FixedSurfaceFlux{FT, VariableFrictionVelocity}) where {FT} = false

Base.@kwdef struct SurfaceBase{FT}
    shf::FT = 0
    lhf::FT = 0
    cm::FT = 0
    ch::FT = 0
    bflux::FT = 0
    ustar::FT = 0
    ρq_tot_flux::FT = 0
    ρq_liq_flux::FT = 0
    ρq_ice_flux::FT = 0
    ρθ_liq_ice_flux::FT = 0
    ρs_flux::FT = 0
    ρu_flux::FT = 0
    ρv_flux::FT = 0
    obukhov_length::FT = 0
    wstar::FT = 0
end

# struct EDMFModel{N_up, FT, SABC, MM, TCM, PM, RFM, PFM, ENT, EBGC, MLP, PMP, EC, MLEC, ET, EDS, DDS, EPG}
struct EDMFModel{N_up, FT, SFCA, SABC, MM, CSM, TCM, PM, RFM, SFM, PFM, ENT, EBGC, MLP, PMP, EC, MLEC, ET, EDS, DDS, EPG, TLT}
    # surface_area::FT
    surface_area::SFCA # trying to allow for non even split of surface area in updrafts... we'll read this in driver/initial_conditions.jl which will call src/EDMF_Functions.jl area_surface_bc()
    surface_area_bc::SABC
    max_area::FT
    minimum_area::FT
    moisture_model::MM
    cloud_sedimentation_model::CSM
    thermo_covariance_model::TCM
    precip_model::PM
    rain_formation_model::RFM
    snow_formation_model::SFM
    precip_fraction_model::PFM
    en_thermo::ENT
    bg_closure::EBGC
    mixing_length_params::MLP
    pressure_model_params::PMP
    entr_closure::EC
    ml_entr_closure::MLEC
    entrainment_type::ET
    entr_dim_scale::EDS
    detr_dim_scale::DDS
    entr_pi_subset::EPG
    set_src_seed::Bool
    H_up_min::FT # minimum updraft top to avoid zero division in pressure drag and turb-entr
    tendency_limiters::TLT
end
function EDMFModel(::Type{FT}, namelist, precip_model, rain_formation_model, param_set) where {FT}

    # Set the number of updrafts (1)
    n_updrafts = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "updraft_number"; default = 1)
    set_src_seed::Bool = namelist["set_src_seed"]

    pressure_func_drag_str = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "pressure_closure_drag";
        default = "normalmode",
        valid_options = ["normalmode", "normalmode_signdf"],
    )

    surface_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "surface_area"; default = 0.1)
    surface_area_bc = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "surface_area_bc";
        default = "Fixed",
        valid_options = ["Fixed", "Prognostic", "Closure"],
    )

    surface_area_bc = if surface_area_bc == "Fixed"
        FixedSurfaceAreaBC()
    elseif surface_area_bc == "Prognostic"
        PrognosticSurfaceAreaBC()
    elseif surface_area_bc == "Closure"
        surface_area_bc_params = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "surface_area_bc_params")
        ClosureSurfaceAreaBC(params = surface_area_bc_params)
    else
        error("Something went wrong. Invalid surface area boundary condition '$surface_area_bc'")
    end

    max_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "max_area"; default = 0.9)
    minimum_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_area"; default = 1e-5)

    moisture_model_name = parse_namelist(namelist, "thermodynamics", "moisture_model"; default = "equilibrium")

    moisture_model = if moisture_model_name == "equilibrium"
        EquilibriumMoisture()
    elseif moisture_model_name == "nonequilibrium"
        NonEquilibriumMoisture(param_set)
    else
        error("Something went wrong. Invalid moisture model: '$moisture_model_name'")
    end

    snow_formation_model = SnowFormationModel(param_set)

    use_sedimentation = TCP.parse_isbits_nt(param_set.user_args, :use_sedimentation, false)
    cloud_sedimentation_model = use_sedimentation ? CloudSedimentationModel(param_set) : CloudNoSedimentationModel()

    tendency_limiters = TendencyLimiterSet(param_set, moisture_model_name)

    thermo_covariance_model_name =
        parse_namelist(namelist, "thermodynamics", "thermo_covariance_model"; default = "prognostic")

    thermo_covariance_model = if thermo_covariance_model_name == "prognostic"
        PrognosticThermoCovariances()
    elseif thermo_covariance_model_name == "diagnostic"
        covar_lim = parse_namelist(namelist, "thermodynamics", "diagnostic_covar_limiter")
        DiagnosticThermoCovariances(covar_lim)
    else
        error("Something went wrong. Invalid thermo_covariance model: '$thermo_covariance_model_name'")
    end

    precip_fraction_model_name =
        parse_namelist(namelist, "microphysics", "precip_fraction_model"; default = "prescribed")

    precip_fraction_model = if precip_fraction_model_name == "prescribed"
        prescribed_precip_frac_value =
            parse_namelist(namelist, "microphysics", "prescribed_precip_frac_value"; default = 1.0)
        PrescribedPrecipFraction(FT(prescribed_precip_frac_value))
    elseif precip_fraction_model_name == "cloud_cover"
        precip_fraction_limiter = parse_namelist(namelist, "microphysics", "precip_fraction_limiter"; default = 0.3)
        DiagnosticPrecipFraction(FT(precip_fraction_limiter))
    else
        error("Something went wrong. Invalid `precip_fraction` model: `$precip_fraction_model_name`")
    end

    # Create the environment variable class (major diagnostic and prognostic variables)

    # Create the class for environment thermodynamics and buoyancy gradient computation
    en_sgs_name =
        parse_namelist(namelist, "thermodynamics", "sgs"; default = "mean", valid_options = ["mean", "quadrature"])
    en_thermo = if en_sgs_name == "mean"
        SGSMean()
    elseif en_sgs_name == "quadrature"
        SGSQuadrature(FT, namelist)
    else
        error("Something went wrong. Invalid environmental sgs type '$en_sgs_name'")
    end
    bg_closure = if en_sgs_name == "mean"
        BuoyGradMean()
    elseif en_sgs_name == "quadrature"
        BuoyGradQuadratures()
    else
        error("Something went wrong. Invalid environmental buoyancy gradient closure type '$en_sgs_name'")
    end
    if moisture_model_name == "nonequilibrium" && en_thermo == "quadrature"
        error("SGS quadratures are not yet implemented for non-equilibrium moisture. Please use the option: mean.")
    end

    # entr closure
    stoch_entr_type = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "stochastic_entrainment";
        default = "deterministic",
        valid_options = [
            "deterministic",
            "noisy_relaxation_process",
            "lognormal_scaling",
            "prognostic_noisy_relaxation_process",
        ],
    )
    entr_type = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "entrainment";
        default = "moisture_deficit",
        valid_options = ["moisture_deficit", "None"],
    )

    ml_entr_type = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "ml_entrainment";
        default = "None",
        valid_options = ["None", "NN", "NN_nonlocal", "FNO", "Linear", "RF"],
    )

    nn_biases = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "nn_ent_biases";
        default = true,
        valid_options = [true, false],
    )

    linear_biases = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "linear_ent_biases";
        default = true,
        valid_options = [true, false],
    )

    w_min = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_upd_velocity")
    c_ε = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_factor")
    μ_0 = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_scale")
    β = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "sorting_power")
    χ = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "updraft_mixing_frac")
    c_λ = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entrainment_smin_tke_coeff")
    γ_lim = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "area_limiter_scale")
    β_lim = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "area_limiter_power")
    limit_min_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "limit_min_area")
    min_area_limiter_scale = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_area_limiter_scale")
    min_area_limiter_power = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_area_limiter_power")
    c_γ = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "turbulent_entrainment_factor")
    c_δ = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "detrainment_factor")
    Π_norm = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pi_norm_consts")
    Π_norm = SA.SVector{length(Π_norm), FT}(Π_norm)
    entr_nondim_norm_factor =
        parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entr_nondim_norm_factor", default = 1.0)
    detr_nondim_norm_factor =
        parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "detr_nondim_norm_factor", default = 1.0)

    εδ_params = εδModelParams{FT, typeof(Π_norm)}(;
        w_min,
        c_ε,
        μ_0,
        β,
        χ,
        c_λ,
        γ_lim,
        β_lim,
        limit_min_area,
        min_area_limiter_scale,
        min_area_limiter_power,
        c_γ,
        c_δ,
        Π_norm,
        entr_nondim_norm_factor,
        detr_nondim_norm_factor,
    )

    entr_pi_subset = Tuple(parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entr_pi_subset"))

    mean_entr_closure = if entr_type == "moisture_deficit"
        MDEntr(; params = εδ_params)
    elseif entr_type == "None"
        EntrNone(; params = εδ_params)
    else
        error("Something went wrong. Invalid entrainment type '$entr_type'")
    end

    ml_entr_closure = if ml_entr_type == "None"
        EntrNone(; params = εδ_params)
    elseif ml_entr_type == "NN"
        c_nn_params = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_ent_params")
        nn_arc = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_arc")
        NNEntr(; params = εδ_params, biases_bool = nn_biases, c_nn_params, nn_arc)
    elseif ml_entr_type == "NN_nonlocal"
        c_nn_params = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_ent_params")
        nn_arc = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "nn_arc")
        NNEntrNonlocal(; params = εδ_params, biases_bool = nn_biases, c_nn_params, nn_arc)
    elseif ml_entr_type == "FNO"
        w_fno = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fno_ent_width")
        nm_fno = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fno_ent_n_modes")
        c_fno = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fno_ent_params")
        FNOEntr(; params = εδ_params, w_fno, nm_fno, c_fno)
    elseif ml_entr_type == "Linear"
        c_linear = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "linear_ent_params")
        LinearEntr(; params = εδ_params, biases_bool = linear_biases, c_linear)
    elseif ml_entr_type == "RF"
        c_rf_fix = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "rf_fix_ent_params")
        c_rf_opt = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "rf_opt_ent_params")
        RFEntr(εδ_params, c_rf_fix, c_rf_opt, length(entr_pi_subset))
    else
        error("Something went wrong. Invalid ML entrainment type '$ml_entr_type'")
    end

    # Overwrite `entr_closure` if a noisy relaxation process is used
    entr_closure = if stoch_entr_type == "noisy_relaxation_process"
        c_gen_stoch = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "general_stochastic_ent_params")
        NoisyRelaxationProcess(; mean_model = mean_entr_closure, c_gen_stoch)
    elseif stoch_entr_type == "lognormal_scaling"
        c_gen_stoch = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "general_stochastic_ent_params")
        LogNormalScalingProcess(; mean_model = mean_entr_closure, c_gen_stoch)
    elseif stoch_entr_type == "deterministic"
        mean_entr_closure
    elseif stoch_entr_type == "prognostic_noisy_relaxation_process"
        c_gen_stoch = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "general_stochastic_ent_params")
        PrognosticNoisyRelaxationProcess(; mean_model = mean_entr_closure, c_gen_stoch)
    else
        error("Something went wrong. Invalid stochastic entrainment type '$stoch_entr_type'")
    end

    # minimum updraft top to avoid zero division in pressure drag and turb-entr
    H_up_min = FT(parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_updraft_top"))


    entrainment_type = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "entrainment_type";
        default = "fractional",
        valid_options = ["fractional", "total_rate"],
    )

    if entrainment_type == "fractional"
        entrainment_type = FractionalEntrModel()
    elseif entrainment_type == "total_rate"
        entrainment_type = TotalRateEntrModel()
    else
        error("Something went wrong. Invalid entrainment type '$entrainment_type'")
    end

    entr_dim_scale = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "entr_dim_scale";
        default = "buoy_vel",
        valid_options = [
            "buoy_vel",
            "inv_scale_height",
            "inv_z",
            "pos_massflux",
            "neg_massflux",
            "abs_massflux",
            "mf_grad",
            "none",
            "b_w",
            "w_height",
            "b_sqrt_tke",
            "sqrt_b_z",
            "tke_b_w",
            "dw_dz",
            "mf_grad_rhoa",
        ],
    )

    pressure_model_params = PressureModelParams{FT}(;
        α_b = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_buoy_coeff1"),
        α_a = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_adv_coeff"),
        α_d = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "pressure_normalmode_drag_coeff"),
    )

    mixing_length_params = MixingLengthParams{FT}(;
        ω_pr = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_scale"),
        c_m = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_ed_coeff"),
        c_d = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_diss_coeff"),
        c_b = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "static_stab_coeff"; default = 0.4), # this is here due to a value error in CliMAParmameters.j,
        κ_star² = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_surf_scale"),
        Pr_n = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Prandtl_number_0"),
        Ri_c = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Ri_crit"),
        Le = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "Lewis_number"; default = 1.0),
        smin_ub = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "smin_ub"),
        l_max = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "l_max"; default = 1.0e6),
    )

    entr_dim_scale = if entr_dim_scale == "buoy_vel"
        BuoyVelEntrDimScale()
    elseif entr_dim_scale == "inv_scale_height"
        InvScaleHeightEntrDimScale()
    elseif entr_dim_scale == "inv_z"
        InvZEntrDimScale()
    elseif entr_dim_scale == "pos_massflux"
        PosMassFluxGradDimScale()
    elseif entr_dim_scale == "neg_massflux"
        NegMassFluxGradDimScale()
    elseif entr_dim_scale == "abs_massflux"
        AbsMassFluxGradDimScale()
    elseif entr_dim_scale == "mf_grad"
        MassFluxGradDimScale()
    elseif entr_dim_scale == "none"
        InvMeterEntrDimScale()
    elseif entr_dim_scale == "b_w" # b/w
        BOverWDimScale()
    elseif entr_dim_scale == "w_height" # w/z
        WOverHeightDimScale()
    elseif entr_dim_scale == "b_sqrt_tke" # b/sqrt(tke)
        BOverSqrtTKEDimScale()
    elseif entr_dim_scale == "sqrt_b_z" # sqrt(b/z)
        SqrtBOverZDimScale()
    elseif entr_dim_scale == "tke_b_w" # (tke*b)/w^3
        TKEBWDimScale()
    elseif entr_dim_scale == "dw_dz"
        DwDzDimScale()
    elseif entr_dim_scale == "mf_grad_rhoa"
        DmDzOverRhoaDimScale()
    else
        error("Something went wrong. Invalid entrainment dimension scale '$entr_dim_scale'")
    end

    detr_dim_scale = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "detr_dim_scale";
        default = "buoy_vel",
        valid_options = [
            "buoy_vel",
            "inv_scale_height",
            "inv_z",
            "pos_massflux",
            "neg_massflux",
            "abs_massflux",
            "mf_grad",
            "none",
            "b_w",
            "w_height",
            "b_sqrt_tke",
            "sqrt_b_z",
            "tke_b_w",
            "dw_dz",
            "mf_grad_rhoa",
        ],
    )

    detr_dim_scale = if detr_dim_scale == "buoy_vel"
        BuoyVelEntrDimScale()
    elseif detr_dim_scale == "inv_scale_height"
        InvScaleHeightEntrDimScale()
    elseif detr_dim_scale == "inv_z"
        InvZEntrDimScale()
    elseif detr_dim_scale == "pos_massflux"
        PosMassFluxGradDimScale()
    elseif detr_dim_scale == "neg_massflux"
        NegMassFluxGradDimScale()
    elseif detr_dim_scale == "abs_massflux"
        AbsMassFluxGradDimScale()
    elseif detr_dim_scale == "mf_grad"
        MassFluxGradDimScale()
    elseif detr_dim_scale == "none"
        InvMeterEntrDimScale()
    elseif detr_dim_scale == "b_w" # b/w
        BOverWDimScale()
    elseif detr_dim_scale == "w_height" # w/z
        WOverHeightDimScale()
    elseif detr_dim_scale == "b_sqrt_tke" # b/sqrt(tke)
        BOverSqrtTKEDimScale()
    elseif detr_dim_scale == "sqrt_b_z" # sqrt(b/z)
        SqrtBOverZDimScale()
    elseif detr_dim_scale == "tke_b_w" # (tke*b)/w^3
        TKEBWDimScale()
    elseif detr_dim_scale == "dw_dz"
        DwDzDimScale()
    elseif detr_dim_scale == "mf_grad_rhoa"
        DmDzOverRhoaDimScale()
    else
        error("Something went wrong. Invalid entrainment dimension scale '$detr_dim_scale'")
    end

    SABC = typeof(surface_area_bc)
    EDS = typeof(entr_dim_scale)
    ET = typeof(entrainment_type)
    DDS = typeof(detr_dim_scale)
    EC = typeof(entr_closure)
    MLEC = typeof(ml_entr_closure)
    MM = typeof(moisture_model)
    CSM = typeof(cloud_sedimentation_model)
    TCM = typeof(thermo_covariance_model)
    PM = typeof(precip_model)
    RFM = typeof(rain_formation_model)
    SFM = typeof(snow_formation_model)
    PFM = typeof(precip_fraction_model)
    EBGC = typeof(bg_closure)
    ENT = typeof(en_thermo)
    EPG = typeof(entr_pi_subset)
    MLP = typeof(mixing_length_params)
    PMP = typeof(pressure_model_params)
    TLT = typeof(tendency_limiters)

    SFCA = typeof(surface_area) # testing allowing this to be a vector for initializing updrafts with non equal areas
    # return EDMFModel{n_updrafts,FT, SABC, MM, TCM, PM, RFM, PFM, ENT, EBGC, MLP, PMP, EC, MLEC, ET, EDS, DDS, EPG, TLT}(
    return EDMFModel{
        n_updrafts,
        FT,
        SFCA,
        SABC,
        MM,
        CSM,
        TCM,
        PM,
        RFM,
        SFM,
        PFM,
        ENT,
        EBGC,
        MLP,
        PMP,
        EC,
        MLEC,
        ET,
        EDS,
        DDS,
        EPG,
        TLT
    }(
        surface_area,
        surface_area_bc,
        max_area,
        minimum_area,
        moisture_model,
        cloud_sedimentation_model,
        thermo_covariance_model,
        precip_model,
        rain_formation_model,
        snow_formation_model,
        precip_fraction_model,
        en_thermo,
        bg_closure,
        mixing_length_params,
        pressure_model_params,
        entr_closure,
        ml_entr_closure,
        entrainment_type,
        entr_dim_scale,
        detr_dim_scale,
        entr_pi_subset,
        set_src_seed,
        H_up_min,
        tendency_limiters,
    )
end

parameter_set(obj) = obj.param_set
n_updrafts(::EDMFModel{N_up}) where {N_up} = N_up
Base.eltype(::EDMFModel{N_up, FT}) where {N_up, FT} = FT
n_Π_groups(m::EDMFModel) = length(m.entr_pi_subset)
entrainment_Π_subset(m::EDMFModel) = m.entr_pi_subset
pressure_model_params(m::EDMFModel) = m.pressure_model_params
mixing_length_params(m::EDMFModel) = m.mixing_length_params

Base.broadcastable(edmf::EDMFModel) = Ref(edmf)

function Base.summary(io::IO, edmf::EDMFModel)
    pns = string.(propertynames(edmf))
    buf = maximum(length.(pns))
    keys = propertynames(edmf)
    vals = repeat.(" ", map(s -> buf - length(s) + 2, pns))
    bufs = (; zip(keys, vals)...)
    print(io, '\n')
    for pn in propertynames(edmf)
        prop = getproperty(edmf, pn)
        # Skip some data:
        prop isa Bool && continue
        prop isa NTuple && continue
        prop isa Int && continue
        prop isa Float64 && continue
        prop isa Float32 && continue
        s = string(
            "  ", # needed for some reason
            getproperty(bufs, pn),
            '`',
            string(pn),
            '`',
            "::",
            '`',
            typeof(prop),
            '`',
            '\n',
        )
        print(io, s)
    end
end


struct State{P, A, T, G}
    prog::P
    aux::A
    tendencies::T
    grid::G
end

function State(prog::P, aux::A, tendencies::T) where {P, A, T}
    grid = Grid(axes(prog.cent))
    return State{P, A, T, typeof(grid)}(prog, aux, tendencies, grid)
end

"""
    column_state(prog, aux, tendencies, colidx)

Create a columnar state given full 3D states
 - `prog` prognostic state
 - `aux` auxiliary state
 - `tendencies` tendencies state
 - `colidx` column index, from ClimaCore's `bycolumn` function

## Example
```julia
bycolumn(axes(prog.cent)) do colidx
    state = TC.column_state(prog, aux, tendencies, colidx)
    ...
end
"""
function column_state(prog, aux, tendencies, colidx)
    prog_cent_column = CC.column(prog.cent, colidx)
    prog_face_column = CC.column(prog.face, colidx)
    aux_cent_column = CC.column(aux.cent, colidx)
    aux_face_column = CC.column(aux.face, colidx)
    tends_cent_column = CC.column(tendencies.cent, colidx)
    tends_face_column = CC.column(tendencies.face, colidx)
    prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)
    tends_column = CC.Fields.FieldVector(cent = tends_cent_column, face = tends_face_column)

    return State(prog_column, aux_column, tends_column)
end

function column_prog_aux(prog, aux, colidx)
    prog_cent_column = CC.column(prog.cent, colidx)
    prog_face_column = CC.column(prog.face, colidx)
    aux_cent_column = CC.column(aux.cent, colidx)
    aux_face_column = CC.column(aux.face, colidx)
    prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)

    return State(prog_column, aux_column, nothing)
end

function column_diagnostics(diagnostics, colidx)
    diag_cent_column = CC.column(diagnostics.cent, colidx)
    diag_face_column = CC.column(diagnostics.face, colidx)
    diag_svpc_column = CC.column(diagnostics.svpc, colidx)
    return CC.Fields.FieldVector(cent = diag_cent_column, face = diag_face_column, svpc = diag_svpc_column)
end


Grid(state::State) = state.grid

float_type(state::State) = eltype(state.prog)
# float_type(field::CC.Fields.Field) = CC.Spaces.undertype(axes(field))
float_type(field::CC.Fields.Field) = eltype(parent(field))
