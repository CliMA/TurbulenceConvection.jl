
"""
Model Regions
"""
# --- Level 0: Root
abstract type AbstractDomain end

# --- Level 1: Categories
abstract type AbstractPrimitiveDomain <: AbstractDomain end
abstract type AbstractDerivedDomain   <: AbstractDomain end

# --- Level 2: Concrete domains
# Original EDMF domains
struct EnvDomain          <: AbstractPrimitiveDomain end
struct BulkDomain         <: AbstractPrimitiveDomain end
struct UpDomain           <: AbstractPrimitiveDomain end
# Regions for cloaking
struct CloakUpDomain      <: AbstractDerivedDomain end
struct CloakDownDomain    <: AbstractDerivedDomain end
struct EnvRemainingDomain <: AbstractDerivedDomain end

const EnvOrUpDomain = Union{EnvDomain, BulkDomain, UpDomain} # for functions that work on both env and updrafts, but this cant be instantiated
# const CloakRegions = Union{CloakUpDomain, CloakDownDomain}


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
    qi_tendency_acnv_dep::FT
    qi_tendency_acnv_dep_is::FT
    qi_tendency_acnv_dep_above::FT
    qi_tendency_acnv_agg_mix::FT
    qi_tendency_acnv_thresh::FT
    ql_tendency_accr_liq_rai::FT
    ql_tendency_accr_liq_ice::FT
    ql_tendency_accr_liq_sno::FT
    qi_tendency_accr_ice_liq::FT
    qi_tendency_accr_ice_rai::FT
    qi_tendency_accr_ice_sno::FT
    #
    qs_tendency_accr_rai_sno::FT # accretion QR by QS [QR -> QS]  (we store it this way, it's always to snow below freezing and [QS -> QR] above freezing.) - we calculate here becasue the outcome is temperature dependent, but we store in aux_tc bc precip is on the grid mean.
end
# null_PrecipitationSources(::Type{FT}; fill_value::FT = FT(NaN)) where {FT} = PrecipFormation{FT}((fill_value for f in fieldnames(PrecipFormation))...)

@generated function null_PrecipitationSources(::Type{FT}; fill_value::FT = FT(NaN)) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    # build expression: PrecipFormation{FT}(fill_value, fill_value, ...)
    ex = :(PrecipFormation{FT}($(map(f -> :(fill_value), fields)...)))
    return ex
end

# Base.:+(a::PrecipFormation{FT}, b::PrecipFormation{FT}) where {FT} = PrecipFormation{FT}((getfield(a, f) + getfield(b, f) for f in fieldnames(PrecipFormation))...)
# Base.:-(a::PrecipFormation{FT}, b::PrecipFormation{FT}) where {FT} = PrecipFormation{FT}((getfield(a, f) - getfield(b, f) for f in fieldnames(PrecipFormation))...)

@generated function Base.:+(a::PrecipFormation{FT}, b::PrecipFormation{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    ex = :(PrecipFormation{FT}($(map(f -> :(a.$f + b.$f), fields)...)))
    return ex
end

@generated function Base.:-(a::PrecipFormation{FT}, b::PrecipFormation{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    ex = :(PrecipFormation{FT}($(map(f -> :(a.$f - b.$f), fields)...)))
    return ex
end


# Base.:*(a::PrecipFormation{FT}, b::FT) where {FT} = PrecipFormation{FT}((getfield(a, f) * b for f in fieldnames(PrecipFormation))...)
# Base.:*(a::FT, b::PrecipFormation{FT}) where {FT} = PrecipFormation{FT}((a * getfield(b, f) for f in fieldnames(PrecipFormation))...)
# Base.:+(a::PrecipFormation{FT}, b::FT) where {FT} = PrecipFormation{FT}((getfield(a, f) + b for f in fieldnames(PrecipFormation))...)
# Base.:+(a::FT, b::PrecipFormation{FT}) where {FT} = PrecipFormation{FT}((a + getfield(b, f) for f in fieldnames(PrecipFormation))...)
# Base.:-(a::PrecipFormation{FT}, b::FT) where {FT} = PrecipFormation{FT}((getfield(a, f) - b for f in fieldnames(PrecipFormation))...)
# Base.:-(a::FT, b::PrecipFormation{FT}) where {FT} = PrecipFormation{FT}((a - getfield(b, f) for f in fieldnames(PrecipFormation))...)

@generated function Base.:*(a::PrecipFormation{FT}, b::FT) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    return :(PrecipFormation{FT}($(map(f -> :(a.$f * b), fields)...)))
end

@generated function Base.:*(a::FT, b::PrecipFormation{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    return :(PrecipFormation{FT}($(map(f -> :(a * b.$f), fields)...)))
end

@generated function Base.:+(a::PrecipFormation{FT}, b::FT) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    return :(PrecipFormation{FT}($(map(f -> :(a.$f + b), fields)...)))
end

@generated function Base.:+(a::FT, b::PrecipFormation{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    return :(PrecipFormation{FT}($(map(f -> :(a + b.$f), fields)...)))
end

@generated function Base.:-(a::PrecipFormation{FT}, b::FT) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    return :(PrecipFormation{FT}($(map(f -> :(a.$f - b), fields)...)))
end

@generated function Base.:-(a::FT, b::PrecipFormation{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(PrecipFormation)
    return :(PrecipFormation{FT}($(map(f -> :(a - b.$f), fields)...)))
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
# null constructor [separate name to not clobber kwdef method]
null_NoneqMoistureSources(::Type{FT}; fill_value::FT = FT(NaN)) where {FT} = NoneqMoistureSources{FT}(fill_value, fill_value)

Base.:+(a::NoneqMoistureSources{FT}, b::NoneqMoistureSources{FT}) where {FT} = NoneqMoistureSources{FT}(a.ql_tendency + b.ql_tendency, a.qi_tendency + b.qi_tendency)
Base.:-(a::NoneqMoistureSources{FT}, b::NoneqMoistureSources{FT}) where {FT} = NoneqMoistureSources{FT}(a.ql_tendency - b.ql_tendency, a.qi_tendency - b.qi_tendency)

Base.:*(a::NoneqMoistureSources{FT}, b::FT) where {FT} = NoneqMoistureSources{FT}(a.ql_tendency * b, a.qi_tendency * b)
Base.:*(a::FT, b::NoneqMoistureSources{FT}) where {FT} = NoneqMoistureSources{FT}(a * b.ql_tendency, a * b.qi_tendency)
Base.:+(a::NoneqMoistureSources{FT}, b::FT) where {FT} = NoneqMoistureSources{FT}(a.ql_tendency + b, a.qi_tendency + b)
Base.:+(a::FT, b::NoneqMoistureSources{FT}) where {FT} = NoneqMoistureSources{FT}(a + b.ql_tendency, a + b.qi_tendency)
Base.:-(a::NoneqMoistureSources{FT}, b::FT) where {FT} = NoneqMoistureSources{FT}(a.ql_tendency - b, a.qi_tendency - b)
Base.:-(a::FT, b::NoneqMoistureSources{FT}) where {FT} = NoneqMoistureSources{FT}(a - b.ql_tendency, a - b.qi_tendency)

"""
    NoneqMoistureSource

Storage for tendency in a condensate field

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct NoneqMoistureSource{FT}
    q_tendency::FT
end
# null constructor [separate name to not clobber kwdef method]
null_NoneqMoistureSource(::Type{FT}; fill_value::FT = FT(NaN)) where {FT} = NoneqMoistureSource{FT}(fill_value)

Base.:+(a::NoneqMoistureSource{FT}, b::NoneqMoistureSource{FT}) where {FT} = NoneqMoistureSource{FT}(a.q_tendency + b.q_tendency)
Base.:-(a::NoneqMoistureSource{FT}, b::NoneqMoistureSource{FT}) where {FT} = NoneqMoistureSource{FT}(a.q_tendency - b.q_tendency)

Base.:*(a::NoneqMoistureSource{FT}, b::FT) where {FT} = NoneqMoistureSource{FT}(a.q_tendency * b)
Base.:*(a::FT, b::NoneqMoistureSource{FT}) where {FT} = NoneqMoistureSource{FT}(a * b.q_tendency)
Base.:+(a::NoneqMoistureSource{FT}, b::FT) where {FT} = NoneqMoistureSource{FT}(a.q_tendency + b)
Base.:+(a::FT, b::NoneqMoistureSource{FT}) where {FT} = NoneqMoistureSource{FT}(a + b.q_tendency)
Base.:-(a::NoneqMoistureSource{FT}, b::FT) where {FT} = NoneqMoistureSource{FT}(a.q_tendency - b)
Base.:-(a::FT, b::NoneqMoistureSource{FT}) where {FT} = NoneqMoistureSource{FT}(a - b.q_tendency)


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
# null constructor [separate name to not clobber kwdef method]
# null_OtherMicrophysicsSources(::Type{FT}; fill_value::FT = FT(NaN)) where {FT} = OtherMicrophysicsSources{FT}(fill_value, fill_value, fill_value, fill_value, fill_value, fill_value)

@generated function null_OtherMicrophysicsSources(::Type{FT}; fill_value::FT = FT(NaN)) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    return :(OtherMicrophysicsSources{FT}($(map(f -> :(fill_value), fields)...)))
end

# Base.:+(a::OtherMicrophysicsSources{FT}, b::OtherMicrophysicsSources{FT}) where {FT} = OtherMicrophysicsSources{FT}((getfield(a, f) + getfield(b, f) for f in fieldnames(OtherMicrophysicsSources))...)
# Base.:-(a::OtherMicrophysicsSources{FT}, b::OtherMicrophysicsSources{FT}) where {FT} = OtherMicrophysicsSources{FT}((getfield(a, f) - getfield(b, f) for f in fieldnames(OtherMicrophysicsSources))...)

@generated function Base.:+(a::OtherMicrophysicsSources{FT}, b::OtherMicrophysicsSources{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}())
    # Build tuple of additions
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a.$f + b.$f), fields)...)))
    return ex
end
@generated function Base.:-(a::OtherMicrophysicsSources{FT}, b::OtherMicrophysicsSources{FT}) where {FT}
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a.$f - b.$f), fields)...)))
    return ex
end

# Base.:*(a::OtherMicrophysicsSources{FT}, b::FT) where {FT} = OtherMicrophysicsSources{FT}((getfield(a, f) * b for f in fieldnames(OtherMicrophysicsSources))...)
# Base.:*(a::FT, b::OtherMicrophysicsSources{FT}) where {FT} = OtherMicrophysicsSources{FT}((a * getfield(b, f) for f in fieldnames(OtherMicrophysicsSources))...)
# Base.:+(a::OtherMicrophysicsSources{FT}, b::FT) where {FT} = OtherMicrophysicsSources{FT}((getfield(a, f) + b for f in fieldnames(OtherMicrophysicsSources))...)
# Base.:+(a::FT, b::OtherMicrophysicsSources{FT}) where {FT} = OtherMicrophysicsSources{FT}((a + getfield(b, f) for f in fieldnames(OtherMicrophysicsSources))...)
# Base.:-(a::OtherMicrophysicsSources{FT}, b::FT) where {FT} = OtherMicrophysicsSources{FT}((getfield(a, f) - b for f in fieldnames(OtherMicrophysicsSources))...)
# Base.:-(a::FT, b::OtherMicrophysicsSources{FT}) where {FT} = OtherMicrophysicsSources{FT}((a - getfield(b, f) for f in fieldnames(OtherMicrophysicsSources))...)

@generated function Base.:*(a::OtherMicrophysicsSources{FT}, b::FT) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a.$f * b), fields)...)))
    return ex
end

@generated function Base.:*(a::FT, b::OtherMicrophysicsSources{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a * b.$f), fields)...)))
    return ex
end

@generated function Base.:+(a::OtherMicrophysicsSources{FT}, b::FT) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a.$f + b), fields)...)))
    return ex
end

@generated function Base.:+(a::FT, b::OtherMicrophysicsSources{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a + b.$f), fields)...)))
    return ex
end

@generated function Base.:-(a::OtherMicrophysicsSources{FT}, b::FT) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a.$f - b), fields)...)))
    return ex
end

@generated function Base.:-(a::FT, b::OtherMicrophysicsSources{FT}) where {FT} # generated runs at compile time to avoid runtime lookup
    fields = fieldnames(OtherMicrophysicsSources)
    ex = :(OtherMicrophysicsSources{FT}($(map(f -> :(a - b.$f), fields)...)))
    return ex
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

abstract type AbstractStalledUpdraftHandler end
struct StalledUpdraftDoNothing <: AbstractStalledUpdraftHandler end # don't do anything. this leads to just a completely stranded updraft. However, probably not a good idea because entr/detr turn off (outside of base_detrainment_rate_inv_s if available) and the remaining stranded updraft can do strange things...
struct StalledUpdraftMixToGridMean <: AbstractStalledUpdraftHandler end # mix the updraft tracers to grid mean but keep the area 
struct StalledUpdraftKill <: AbstractStalledUpdraftHandler end # the original default. set area to 0 and thus remove all tracers. This is most stable but the holes left in w=0 regions can make it hard for other updrafts to continue..., and possibly unnecessarily mix all that air into the env immediately, losing local buoyancy gains. if the updraft starts from the surface it's probably ok, but this becomes very bad for elevatd convection, which gets eroded from all sides by the strong area gradients. that develop.
struct StalledUpdraftDetrainDowndrafts <: AbstractStalledUpdraftHandler end # if the buoyancy is negative (θ_liq_ice_up < θ_liq_ice_gm), detrains enough area to bring it back into balance. this is becaue we don't have downdrafts so we can have runaway growth of negative buoyancy and temperature...



# ========== Convective TKE Handler ============================================================================== #

abstract type AbstractConvectiveTKEHandler end
struct NoConvectiveTKE <: AbstractConvectiveTKEHandler end
struct ConvectiveTKE{FT} <: AbstractConvectiveTKEHandler
    "buoyancy_coeff::FT"
    buoyancy_coeff::FT
    "advection_coeff::FT"
    advection_coeff::FT
    "dissipation_coeff::FT"
    dissipation_coeff::FT # smaller than generation
    "self_dissipation_coeff::FT"
    self_dissipation_coeff::FT # extra dissipation when tke is convectively generated
    "max_scaling_factor::FT"
    max_scaling_factor::FT
    # "use_separate_ed_coeff::Bool" # We'll just use if ed_scaling_factor != 1.0
    # use_separate_ed_coeff::Bool
    "ed_scaling_factor::FT"
    ed_scaling_factor::FT
    "transport_tke_by_advection::Bool"
    transport_tke_by_advection::Bool
    "transport_condensed_by_advection::Bool"
    transport_condensed_by_advection::Bool
    "transport_conserved_by_advection::Bool"
    transport_conserved_by_advection::Bool
    "entr_detr_rate_inv_s::FT"
    entr_detr_rate_inv_s::FT # Time rate for grafting onto updraft
end
function ConvectiveTKE(param_set::APS, namelist)
    FT = eltype(param_set)
    buoyancy_coeff = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_buoyancy_coeff"; default = FT(1.0))
    advection_coeff = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_advection_coeff"; default = FT(1.0))
    dissipation_coeff = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_dissipation_coeff"; default = FT(1.0))
    self_dissipation_coeff = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_self_dissipation_coeff"; default = FT(1.0))
    max_scaling_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_max_scaling_factor"; default = FT(1.0))
    # use_separate_ed_coeff = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_use_separate_ed_coeff"; default = false, valid_options = [true, false])
    ed_scaling_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_ed_scaling_factor"; default = FT(1.0))
    transport_tke_by_advection = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_transport_tke_by_advection"; default = true, valid_options = [true, false])
    transport_condensed_by_advection = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_transport_condensed_by_advection"; default = false, valid_options = [true, false]) # ql, qi, qr, qs
    transport_conserved_by_advection = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "convective_tke_transport_conserved_by_advection"; default = true, valid_options = [true, false]) # θ_liq_ice, qt
    entr_detr_rate_inv_s = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "tke_conv_entr_detr_rate_inv_s"; default = FT(0)) # default 0 means no grafting onto updraft (infinite time scale)
    return ConvectiveTKE{FT}(
        buoyancy_coeff,
        advection_coeff,
        dissipation_coeff,
        self_dissipation_coeff,
        max_scaling_factor,
        # use_separate_convective_tke_ed_coeff, # we'll just use if convective_tke_ed_scaling_factor != 1.0
        ed_scaling_factor,
        transport_tke_by_advection,
        transport_condensed_by_advection,
        transport_conserved_by_advection,
        entr_detr_rate_inv_s,
    )            
end

# ========== Entrainment/Detrainment Models ============================================================================== #

abstract type EntrModelFacTotalType end
# struct FractionalEntrModel <: EntrModelFacTotalType end
Base.@kwdef struct FractionalEntrModel{SUHT} <: EntrModelFacTotalType
    "stalled_updraft_handler::SUHT"
    stalled_updraft_handler::SUHT # Whether or not we mix the updraft to the grid mean when we get to w=0
end

# struct TotalRateEntrModel <: EntrModelFacTotalType end
Base.@kwdef struct TotalRateEntrModel{FT, SUHT} <: EntrModelFacTotalType
    "base detrainment rate inverse [s^-1]"
    base_detrainment_rate_inv_s::FT # testing having some base detrainment rate so that even when w goes to 0 we can detrain slowly.
    "base_entrainment_rate_inv_s::FT"
    base_entrainment_rate_inv_s::FT # testing having some base entrainment rate so [[ this gets applied if ∂θ_virt/∂z < 0 ]] So we can start updrafts from rest in unstable layers
    "stalled_updraft_handler::SUHT"
    stalled_updraft_handler::SUHT # Whether or not we mix the updraft to the grid mean when we get to w=0
end

# ========== Area Partition Models ============================================================================== #

abstract type AbstractAreaPartitionModel end
struct StandardAreaPartitionModel{FT} <: AbstractAreaPartitionModel
    "apply_second_order_flux_correction::Bool"
    apply_second_order_flux_correction::Bool
    "second_order_correction_limit_factor::FT"
    second_order_correction_limit_factor::FT # factor to limit the second order correction to avoid overshoots. 0 means no correction, 1 means full correction

    end;
function StandardAreaPartitionModel(param_set::APS, namelist)
    FT = eltype(param_set)
    apply_second_order_flux_correction = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "apply_second_order_flux_correction"; default = false)
    second_order_correction_limit_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "second_order_correction_limit_factor"; default = FT(Inf)) # default Inf (no limit)
    return StandardAreaPartitionModel{FT}(apply_second_order_flux_correction, second_order_correction_limit_factor)
end

#=
            Based on the idea that for a lot of convection, sure the core does most of the MF, but the cloak has much of the area and is needed for correctly diagnosing environmental advection fluxes

            Each updraft has an updraft cloak that has a_cloak_up = cloak_area_factor * core_area, 
            Then there's a downdraft cloak with area  a_cloak_dn = a_up + a_cloak_up.

            max_combined_updraft_area = edmf.max_area for simplicity... so in regions of max_area and strong detrainment, no cloak.

            The combined cloak area [[ a_cloak = a_cloak_up + a_cloak_dn]] is limited so that the total a_bulk + a_bulk_cloak <= max(a_bulk, max_combined_updraft_area), and a_cloak_up and a_cloak_dn are scaled accordingly

            Let f_cm = f_cloak_mix ==> The fractional combination of env and updraft that creates this updraft cloak. Default is 0.5

            Then
                q_cloak_up = f_cloak_mix * q_updraft + (1 - f_cloak_mix) * q_env
                w_updraft_cloak_up = f_cloak_mix * w_updraft #  + (1 - f_cloak_mix) * w_mean  (use w_mean to ensure positivity, w_mean should be 0...)

                # For the downdraft, we just have to close the budget. we need qi_mean to remain unchanged, and w_mean to remain unchanged
                q_cloak_dn = (qi_mean - (a_cloak_up * q_cloak_up) - (a_up * q_updraft)) / a_cloak_dn
                w_cloak_dn = (w_mean - (a_cloak_up * w_cloak_up) - (a_up * w_updraft)) / a_cloak_dn

                NOTE: We need to limit the above to ensure positivity for q... we are ok with negative w.
                This does pose some challenges since if q_env = 0, this precludes any vertical advection in the cloak when it might really be real...

            To calculate advective tendencies in the env then, we just do these cloak regions. We do need to limit losses so that we don't get negative ql or qi in the env
        
            We derive an upper limit for q_cloak_up such that q_cloak_dn >= 0.
                This is:  q_cloak_up <= (qi_mean - a_up * q_updraft) / a_cloak_up
        =#
struct CoreCloakAreaPartitionModel{FT} <: AbstractAreaPartitionModel
    # max_combined_updraft_area::FT # Might be hard to have this be separate from edmf.max_area...
    cloak_area_factor::FT # Scaling factor for how much larger the cloak area is than the core area
    cloak_dn_area_factor::FT # Scaling factor for the ratio of downdraft cloak to (updraft + updraft cloak) area. a_cloak_dn = cloak_dn_area_factor * (a_up + a_cloak_up) up to area limit
    cloak_mix_factor::FT # Scaling factor that decides how close the cloak is to the env vs the core, 1 means fully updraft, 0 means fully env
    confine_all_downdraft_to_cloak::Bool # If true, all downdraft area is put into the cloak, if false, downdraft area is split between cloak and env based on remaining area.
    apply_second_order_flux_correction::Bool # If true, the fluxes from core and cloak are combined before computing gradients, otherwise gradients are computed separately and then combined. This is NOT as necessary for StandardAreaPartitionMode because we have prognostic up and env, though env is backed out so it kind of matters
    second_order_correction_limit_factor::FT # factor to limit the second order correction to avoid overshoots. 0 means no correction, 1 means full correction
    apply_cloak_to_condensate_formation::Bool # Whether or not to apply the cloak area partitioning to condensate formation tendencies. Default true.
    fraction_of_area_above_max_area_allowed_in_cloak::FT
end

"""
    We could either
    - confine the downdraft to the cloak only (confine_all_downdraft_to_cloak = true)
    - construct the cloak downdraft such that w_en remains unchanged
    - have no downdraft cloak, and just let the remaining environment take the enhanced downdraft. in this case however, we'd have to change all the env properties... so for the math/code implementation it would probably be easier to make everything left in env be the cloak.

"""
function CoreCloakAreaPartitionModel(param_set::APS, namelist)
    FT = eltype(param_set)
    cloak_area_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "cloak_area_factor"; default = FT(4.0))
    cloak_dn_area_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "cloak_dn_area_factor"; default = FT(1.0))
    cloak_mix_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "cloak_mix_factor"; default = FT(0.5))
    confine_all_downdraft_to_cloak = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "confine_all_downdraft_to_cloak"; default = false)
    # combine_fluxes_before_gradients = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "combine_fluxes_before_gradients"; default = true) # I'm not sure this does anything w/o prognostic q in each region. Instead we opt for a second order correction....
    apply_second_order_flux_correction = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "apply_second_order_flux_correction"; default = false)
    second_order_correction_limit_factor = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "second_order_correction_limit_factor"; default = FT(Inf)) # no limit
    apply_cloak_to_condensate_formation = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "apply_cloak_to_condensate_formation"; default = true)
    fraction_of_area_above_max_area_allowed_in_cloak = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "fraction_of_area_above_max_area_allowed_in_cloak"; default = FT(0.5))
    return CoreCloakAreaPartitionModel{FT}(cloak_area_factor, cloak_dn_area_factor, cloak_mix_factor, confine_all_downdraft_to_cloak, apply_second_order_flux_correction, second_order_correction_limit_factor, apply_cloak_to_condensate_formation, fraction_of_area_above_max_area_allowed_in_cloak)
end
# ============================================================================================================= #


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
    "liquid specific humidity in the saturated part" # [[ my addition ]] # see https://github.com/CliMA/ClimaAtmos.jl/issues/2429 https://github.com/CliMA/ClimaAtmos.jl/pull/2401
    ql_sat::FT
    "ice specific humidity in the saturated part" # [[ my addition ]] # see https://github.com/CliMA/ClimaAtmos.jl/issues/2429 https://github.com/CliMA/ClimaAtmos.jl/pull/2401 
    qi_sat::FT
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

include("relaxation_timescale_types.jl")

abstract type AbstractMoistureModel end
# struct EquilibriumMoisture <: AbstractMoistureModel end
# struct NonEquilibriumMoisture <: AbstractMoistureModel end


struct EquilibriumMoisture{RTT <: AbstractNonEquillibriumSourcesType, FT <: AbstractFloat} <: AbstractMoistureModel
    scheme::RTT # this is to allow us choose how to diagnose N for example, even when we are not using the timescales. if it's a type that defaults to NaNs, that's ok too but at least it won't be 0. We can resolve the NaNs before saving.
    condensate_qt_SD::FT
end
function EquilibriumMoisture(param_set::APS, namelist)
    nonequilibrium_moisture_scheme_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "nonequilibrium_moisture_scheme"; default = :relax_to_equilibrium))

    nonequilibrium_moisture_scheme = if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
        RelaxToEquilibrium()
    elseif nonequilibrium_moisture_scheme_type === :KorolevMazin2007
        KorolevMazin2007()
    elseif (nonequilibrium_moisture_scheme_type ∈ valid_relaxation_timescale_types) && (nonequilibrium_moisture_scheme_type ∉ [:neural_network, :neural_network_no_weights])
        get_relaxation_timescale_type(nonequilibrium_moisture_scheme_type, param_set, namelist) # we really could just pass namelist here...
    elseif nonequilibrium_moisture_scheme_type ∈ [:neural_network, :neural_network_no_weights]
        # we decided to stop storing NN things in the param_set... it's redundant for all the objects but it's especially bad for the NN, as every call w/ param_set has to deal w/ massive tuple/svector
        get_relaxation_timescale_type(nonequilibrium_moisture_scheme_type, param_set, namelist)
    else
        error("Invalid nonequilibrium_moisture_scheme type: $nonequilibrium_moisture_scheme_type")
    end

    FT = eltype(param_set)
    condensate_qt_SD::FT = parse_namelist(namelist, "user_params", "condensate_qt_SD"; default = zero(FT))

    RTT = typeof(nonequilibrium_moisture_scheme)
    return EquilibriumMoisture{RTT, FT}(nonequilibrium_moisture_scheme, condensate_qt_SD)
end



struct NonEquilibriumMoisture{RTT <: AbstractNonEquillibriumSourcesType, FT <: AbstractFloat} <: AbstractMoistureModel
    scheme::RTT
    heterogeneous_ice_nucleation::NamedTuple{(:use_heterogeneous_ice_nucleation, :heterogeneous_ice_nucleation_coefficient, :heterogeneous_ice_nucleation_exponent, :use_ice_mult), Tuple{Bool, FT, FT, Bool}}
    # heterogeneous_ice_nucleation::Tuple{Bool, FT, FT, Bool}
    condensate_qt_SD::FT
end

function NonEquilibriumMoisture(param_set::APS, namelist)
    FT = eltype(param_set)
    nonequilibrium_moisture_scheme_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "nonequilibrium_moisture_scheme"; default = :relax_to_equilibrium))

    nonequilibrium_moisture_scheme = if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
        RelaxToEquilibrium()
    elseif nonequilibrium_moisture_scheme_type === :KorolevMazin2007
        KorolevMazin2007()
    elseif (nonequilibrium_moisture_scheme_type ∈ valid_relaxation_timescale_types) && (nonequilibrium_moisture_scheme_type ∉ [:neural_network, :neural_network_no_weights])
        get_relaxation_timescale_type(nonequilibrium_moisture_scheme_type, param_set, namelist) # we really could just pass namelist here...
    elseif nonequilibrium_moisture_scheme_type ∈ [:neural_network, :neural_network_no_weights]
        # we decided to stop storing NN things in the param_set... it's redundant for all the objects but it's especially bad for the NN, as every call w/ param_set has to deal w/ massive tuple/svector
        get_relaxation_timescale_type(nonequilibrium_moisture_scheme_type, param_set, namelist)
    else
        error("Invalid nonequilibrium_moisture_scheme type: $nonequilibrium_moisture_scheme_type")
    end

    use_heterogeneous_ice_nucleation = parse_namelist(namelist, "user_args", "use_heterogeneous_ice_nucleation"; default = false)
    heterogeneous_ice_nucleation_coefficient = parse_namelist(namelist, "user_params", "heterogeneous_ice_nucleation_coefficient"; default = FT(1))
    heterogeneous_ice_nucleation_exponent = parse_namelist(namelist, "user_params", "heterogeneous_ice_nucleation_exponent"; default = FT(1))
    use_ice_mult = parse_namelist(namelist, "user_args", "use_ice_mult"; default = false)

    # heterogeneous_ice_nucleation_named_tuple = (heterogeneous_ice_nucleation, heterogeneous_ice_nucleation_coefficient, heterogeneous_ice_nucleation_exponent, use_ice_mult) 
    heterogeneous_ice_nucleation_named_tuple = (; use_heterogeneous_ice_nucleation, heterogeneous_ice_nucleation_coefficient, heterogeneous_ice_nucleation_exponent, use_ice_mult) 

    condensate_qt_SD::FT = parse_namelist(namelist, "user_params", "condensate_qt_SD"; default = zero(FT))

    RTT = typeof(nonequilibrium_moisture_scheme)
    return NonEquilibriumMoisture{RTT, FT}(nonequilibrium_moisture_scheme, heterogeneous_ice_nucleation_named_tuple, condensate_qt_SD)
end

"""
NR Closures (assumptions about the size distribution) -- used in relaxation timescale types so make sure it's before
"""
# abstract type AbstractNRClosureType end
# struct GammaNRClosure <: AbstractNRClosureType end # default
# struct MonodisperseNRClosure <: AbstractNRClosureType end # redundant...
# struct FixedRadiusNRClosure <: AbstractNRClosureType end # r fixed
# struct InhomogeneousNRClosure <: AbstractNRClosureType end # N fixed


# Sedimentation
abstract type AbstractIntegrationScheme end
struct UpwindDifferencingScheme <: AbstractIntegrationScheme end
struct RightBiasedDifferencingScheme <: AbstractIntegrationScheme end
abstract type AbstractCloudSedimentationModel end
# struct CloudSedimentationModel{FT <: AbstractFloat, LTVS <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, ITVS <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, SLNC <: Union{AbstractFloat, AbstractRelaxationTimescaleType}, SINC <: Union{AbstractFloat, AbstractRelaxationTimescaleType}, SIM <: AbstractIntegrationScheme} <: AbstractCloudSedimentationModel
struct CloudSedimentationModel{FT <: AbstractFloat, LTVS <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, ITVS <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, SIM <: AbstractIntegrationScheme} <: AbstractCloudSedimentationModel
    liq_terminal_velocity_scheme::LTVS
    ice_terminal_velocity_scheme::ITVS
    # sedimentation_liq_number_concentration::SLNC # we now cache N_i, which derives directly from moisture_model, and term_vel_ice, derived from that N_i, so we don't need this.
    # sedimentation_ice_number_concentration::SINC # we now cache N_i, which derives directly from moisture_model, and term_vel_ice, derived from that N_i, so we don't need this.
    liq_Dmax::FT
    ice_Dmax::FT
    liq_sedimentation_scaling_factor::FT
    ice_sedimentation_scaling_factor::FT
    sedimentation_differencing_scheme::SIM
    E_liq_liq::FT # liq_liq_collision_efficiency (E_liq_liq to match other micorphysics efficiencies)
    E_liq_ice::FT # liq_ice_collision_efficiency (E_liq_ice to match other micorphysics efficiencies)
    E_ice_ice::FT # ice_ice_collision_efficiency (E_ice_ice to match other micorphysics efficiencies)
    E_ice_ice_mix::FT # testing
    # liq_ice_collision_scaling_factor::FT # this should be deprecated one day but I'm using it for initial calibration just to check if we're even close
    # liq_sno_accretion_scaling_factor::FT # could also go in snow_formation I guess [[ E_liq_sno already exists in microphysics so we'll just use that...]]
    grid_mean::Bool # whether or not sedimentation is applied to the grid mean only instead of in env/up separately
end

struct CloudNoSedimentationModel <: AbstractCloudSedimentationModel end

function CloudSedimentationModel(param_set::APS, moisture_model::AbstractMoistureModel, namelist)
    FT = eltype(param_set)

    liq_terminal_velocity_scheme = get_termvel_type(Symbol(parse_namelist(namelist, "user_args", "liq_terminal_velocity_scheme"; default = :Blk1MVel)))
    ice_terminal_velocity_scheme = get_termvel_type(Symbol(parse_namelist(namelist, "user_args", "ice_terminal_velocity_scheme"; default = :Blk1MVel)))


    liq_Dmax = parse_namelist(namelist, "user_params", "liq_sedimentation_Dmax"; default =  FT(Inf))
    ice_Dmax = parse_namelist(namelist, "user_params", "ice_sedimentation_Dmax"; default =  FT(62.5e-6 * 2)) # maybe this should also be inf? for chen it's not clear...

    # NaN's (e.g. from json) get converted to Inf here
    liq_Dmax = isnan(liq_Dmax) ? FT(Inf) : liq_Dmax
    ice_Dmax = isnan(ice_Dmax) ? FT(Inf) : ice_Dmax


    liq_sedimentation_scaling_factor = parse_namelist(namelist, "user_params", "liq_sedimentation_scaling_factor"; default = FT(1.0))
    ice_sedimentation_scaling_factor = parse_namelist(namelist, "user_params", "ice_sedimentation_scaling_factor"; default = FT(1.0))

    sedimentation_differencing_scheme_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "sedimentation_differencing_scheme"; default = :upwinding))
    sedimentation_differencing_scheme = if sedimentation_differencing_scheme_type === :upwinding
        UpwindDifferencingScheme()
    elseif sedimentation_differencing_scheme_type === :right_biased
        RightBiasedDifferencingScheme()
    else
        error("Invalid sedimentation_differencing_scheme: $sedimentation_differencing_scheme_type")
    end

    E_liq_liq = parse_namelist(namelist, "user_params", "E_liq_liq"; default = FT(1.0)) # this is the collision efficiency for liquid to liquid collisions, it should be 1.0 in most cases but we can use it to tune the model a bit
    E_liq_ice = parse_namelist(namelist, "user_params", "E_liq_ice"; default = FT(1.0)) # this is the collision efficiency for liquid to ice collisions, it should be 1.0 in most cases but we can use it to tune the model a bit
    E_ice_ice = parse_namelist(namelist, "user_params", "E_ice_ice"; default = FT(1.0)) # this is the collision efficiency for ice to ice collisions, it should be 1.0 in most cases but we can use it to tune the model a bit

    E_ice_ice_mix = parse_namelist(namelist, "user_params", "E_ice_ice_mix"; default = FT(1.0)) # this is the collision efficiency for ice to ice collisions in the mixed phase, it should be 1.0 in most cases but we can use it to tune the model a bit

    # liq_sno_accretion_scaling_factor = parse_namelist(namelist, "user_params", "liq_sno_accretion_scaling_factor"; default = FT(1.0)) # this is the scaling factor for the liquid to snow accretion term. It should be 1.0 in most cases but we can use it to tune the model a bit

    LTVS = typeof(liq_terminal_velocity_scheme)
    ITVS = typeof(ice_terminal_velocity_scheme)
    # SLNC = typeof(sedimentation_liq_number_concentration)
    # SINC = typeof(sedimentation_ice_number_concentration)
    SDS = typeof(sedimentation_differencing_scheme)

    grid_mean::Bool = parse_namelist(namelist, "user_args", "grid_mean_sedimentation"; default = false) # only apply to grid mean values, never actually got it to work though


    # return CloudSedimentationModel{FT, LTVS, ITVS, SLNC, SINC, SDS}(
    return CloudSedimentationModel{FT, LTVS, ITVS, SDS}(
        liq_terminal_velocity_scheme, 
        ice_terminal_velocity_scheme,
        # sedimentation_liq_number_concentration,
        # sedimentation_ice_number_concentration,
        liq_Dmax,
        ice_Dmax,
        liq_sedimentation_scaling_factor,
        ice_sedimentation_scaling_factor,
        sedimentation_differencing_scheme,
        E_liq_liq,
        E_liq_ice,
        E_ice_ice,
        E_ice_ice_mix, # testing
        # liq_ice_collision_scaling_factor,
        # liq_sno_accretion_scaling_factor,
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
# struct Clima1M{FT, RTVST <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, STVST <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}} <: AbstractPrecipitationModel 
struct Clima1M{RTVST <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}, STVST <: Union{CMT.Blk1MVelType, CMT.Chen2022Type}} <: AbstractPrecipitationModel 
    # rain_sedimentation_scaling_factor::FT # just use χv
    # snow_sedimentation_scaling_factor::FT # just use χv
    rain_terminal_velocity_scheme::RTVST
    snow_terminal_velocity_scheme::STVST
    # sedimentation_differencing_scheme::CMT.AbstractSedimentationDifferencingScheme # This is shared w/ condensate but there's only one option, so i feel like it should be in only one place
end
function Clima1M(param_set::APS, namelist)
    # FT = eltype(param_set)

    # rain_sedimentation_scaling_factor::FT = parse_namelist(namelist, "user_params", "rain_sedimentation_scaling_factor"; default =  FT(1.0))
    # snow_sedimentation_scaling_factor::FT = parse_namelist(namelist, "user_params", "snow_sedimentation_scaling_factor"; default =  FT(1.0))
    rain_terminal_velocity_scheme = get_termvel_type(Symbol(parse_namelist(namelist, "user_args", "rain_terminal_velocity_scheme"; default =  :Blk1MVel)))
    snow_terminal_velocity_scheme = get_termvel_type(Symbol(parse_namelist(namelist, "user_args", "snow_terminal_velocity_scheme"; default =  :Blk1MVel)))


    # sedimentation_differencing_scheme_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "sedimentation_differencing_scheme"; default = :upwinding))
    # sedimentation_differencing_scheme = if sedimentation_differencing_scheme_type === :upwinding
    #     UpwindDifferencingScheme()
    # elseif sedimentation_differencing_scheme_type === :right_biased
    #     RightBiasedDifferencingScheme()
    # else
    #     error("Invalid sedimentation_differencing_scheme: $sedimentation_differencing_scheme_type")
    # end

    RTVST = typeof(rain_terminal_velocity_scheme)
    STVST = typeof(snow_terminal_velocity_scheme)
    # SDST = typeof(sedimentation_differencing_scheme)
    # return Clima1M{FT, RTVST, STVST}(rain_sedimentation_scaling_factor, snow_sedimentation_scaling_factor, rain_terminal_velocity_scheme, snow_terminal_velocity_scheme) # , sedimentation_differencing_scheme)
    return Clima1M{RTVST, STVST}(rain_terminal_velocity_scheme, snow_terminal_velocity_scheme) # , sedimentation_differencing_scheme)
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
struct NonEquilibriumSnowFormationModel{FT} <: AbstractSnowFormationModel
    ice_dep_acnv_scaling_factor::FT # this is the scaling factor for the ice deposition/accretion term. It should be 1.0 in most cases but we can use it to tune the model a bit
    ice_dep_acnv_scaling_factor_above::FT
    ice_acnv_power::FT # this is the power for the ice accretion term, it should be 1.0 in most cases but we can use it to tune the model a bit
    r_ice_acnv_scaling_factor::FT # this is the scaling factor for the rain to ice accretion term. It should be 1.0 in most cases but we can use it to tune the model a bit
    r_ice_snow_threshold_scaling_factor::FT # factor for how much higher the forced <r> auto convert threshold is than r_ice_snow
end

function NonEquilibriumSnowFormationModel(param_set::APS, namelist)
    FT = eltype(param_set)
    ice_dep_acnv_scaling_factor = parse_namelist(namelist, "user_params", "ice_dep_acnv_scaling_factor"; default = 1.0)
    ice_dep_acnv_scaling_factor_above = parse_namelist(namelist, "user_params", "ice_dep_acnv_scaling_factor_above"; default = 1.0)
    ice_acnv_power = parse_namelist(namelist, "user_params", "ice_acnv_power"; default = 1.0)
    r_ice_acnv_scaling_factor = parse_namelist(namelist, "user_params", "r_ice_acnv_scaling_factor"; default = 1.0)
    r_ice_snow_threshold_scaling_factor = parse_namelist(namelist, "user_params", "r_ice_snow_threshold_scaling_factor"; default = 1.0)
    return NonEquilibriumSnowFormationModel{FT}(ice_dep_acnv_scaling_factor, ice_dep_acnv_scaling_factor_above, ice_acnv_power, r_ice_acnv_scaling_factor, r_ice_snow_threshold_scaling_factor)
end

function SnowFormationModel(param_set::APS, moisture_model::AbstractMoistureModel, namelist)

    if moisture_model isa NonEquilibriumMoisture
        snow_formation_model = if moisture_model.scheme isa RelaxToEquilibrium
            DefaultSnowFormationModel()
        elseif moisture_model.scheme isa KorolevMazin2007
            DefaultSnowFormationModel() # should this be default? I guess maybe cause there's no N prediction?
        elseif moisture_model.scheme isa AbstractRelaxationTimescaleType
            NonEquilibriumSnowFormationModel(param_set, namelist)
        else
            error("Invalid moisture_model type: $(typeof(moisture_model))")
        end
    elseif moisture_model isa EquilibriumMoisture
        snow_formation_model = DefaultSnowFormationModel()
    else
        error("Invalid moisture model type: $(typeof(moisture_model))")
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

struct TruncatedBasicTendencyLimiter{FT} <: AbstractTendencyLimiter 
    factor::FT
end # Same as BasicTendencyLimiter but scaled down so things take more than 1 timestep. For stability.

limit_tendency(limiter::TruncatedBasicTendencyLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = max(x, -x_avail/Δt * limiter.factor)
limit_tendency(limiter::TruncatedBasicTendencyLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = safe_clamp(x, -x_up/Δt * limiter.factor, x_en/Δt * limiter.factor)

# Maybe consider adding
# abstract type AbstractdtDependentTendencyLimiter <: AbstractTendencyLimiter end that can supertype all dt dependent tendencies... rn using !isa(NoTendencyLimiter, limiter) is a bit hacky
# could also try to find a way to unify NoTendency and NoMoistureSourcesLimiter ? idk.
# abstract type N

# Moisture Limiters
abstract type AbstractMoistureSourcesLimiter <: AbstractTendencyLimiter end
struct NoMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end # No limiting is done, the sources are just directly calculated
struct BasicMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end # Sources are truncated to their available value within a timestep
struct TruncatedBasicMoistureSourcesLimiter{FT} <: AbstractMoistureSourcesLimiter 
    factor::FT
end # Same as BasicMoistureSourcesLimiter but scaled down so things take more than 1 timestep. For stability.
struct StandardSupersaturationMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end # Sources are truncated to their available value within a timestep but this is geared towards supersaturation relaxation. It does a lot of checks that are unnecessary if you're just relaxing to equilibrium e.g. WBF regimes etc.
struct MorrisonMilbrandt2015MoistureSourcesLimiter <: AbstractMoistureSourcesLimiter end
struct MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter <: AbstractMoistureSourcesLimiter
    fallback_to_standard_supersaturation_limiter::Bool
end
# Fallbacks (in noneq_moisture sources we have our own setup, but for dispatch these may be needed elsewhere)
limit_tendency(::AbstractMoistureSourcesLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = max(x, -x_avail/Δt)
limit_tendency(::AbstractMoistureSourcesLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = safe_clamp(x, -x_up/Δt, x_en/Δt) 
limit_tendency(::NoMoistureSourcesLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = x
limit_tendency(::NoMoistureSourcesLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = x
limit_tendency(limiter::TruncatedBasicMoistureSourcesLimiter, x::FT, x_avail::FT, Δt::FT) where {FT} = max(x, -x_avail/Δt * limiter.factor)
limit_tendency(limiter::TruncatedBasicMoistureSourcesLimiter, x::FT, x_up::FT, x_en::FT, Δt::FT) where {FT} = safe_clamp(x, -x_up/Δt * limiter.factor, x_en/Δt * limiter.factor)

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
    use_tendency_resolver::Bool
    use_tendency_resolver_on_full_tendencies::Bool
end

get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:default}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_default_tendency_limiter : tendency_limiter_set.default_tendency_limiter
get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:moisture_sources}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_moisture_sources_limiter : tendency_limiter_set.moisture_sources_limiter
get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:entr_detr}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_entr_detr_tendency_limiter : tendency_limiter_set.entr_detr_tendency_limiter
get_tendency_limiter(tendency_limiter_set::TendencyLimiterSet, ::Val{:precipitation}, use_fallback::Bool) = use_fallback ? tendency_limiter_set.fallback_precipitation_tendency_limiter : tendency_limiter_set.precipitation_tendency_limiter

function TendencyLimiterSet(param_set::APS, moisture_model::AbstractMoistureModel, namelist)
    FT = eltype(param_set)

    factor::FT = parse_namelist(namelist, "user_args", "truncated_basic_limiter_factor"; default = FT(1.0)) # ideally in user_args so it's with the other limiter settings but i don't remember if scalars are ok there...

    # ------------------ Default Tendency Limiter ------------------ #
    default_tendency_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "default_tendency_limiter_type"; default = :none))  # Default to none bc that's purest (basic is how anna wrote it though)
    default_tendency_limiter = if default_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif default_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    elseif default_tendency_limiter_type === :truncated_basic
        TruncatedBasicTendencyLimiter(factor)
    else
        error("Invalid default_tendency_limiter_type: $default_tendency_limiter_type, valid options are :none, :basic")
    end

    fallback_default_tendency_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "fallback_default_tendency_limiter_type"; default = :none))  # Default to none bc that's purest (basic is how anna wrote it though)
    fallback_default_tendency_limiter = if fallback_default_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif fallback_default_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    elseif fallback_default_tendency_limiter_type === :truncated_basic
        TruncatedBasicTendencyLimiter(factor)
    else
        error("Invalid fallback_default_tendency_limiter_type: $fallback_default_tendency_limiter_type, valid options are :none, :basic")
    end

    # ------------------ Moisture Sources Limiter ------------------ #
    if moisture_model isa EquilibriumMoisture
        # can we change this to create no object at all? Or not, idk... maybe it's not that wasteful and improves type inference/stability.
        moisture_sources_limiter = NoMoistureSourcesLimiter() # if we're not a nonequilibrium moisture model, there's no limiting just sat adjust (liq/ice are diagnosed)...can just use NoMoistureSourcesLimiter
        fallback_moisture_sources_limiter = NoMoistureSourcesLimiter()
    elseif moisture_model isa NonEquilibriumMoisture

        nonequilibrium_moisture_scheme_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "nonequilibrium_moisture_scheme"; default = :relax_to_equilibrium))

        nonequilibrium_moisture_sources_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "nonequilibrium_moisture_sources_limiter_type"; default = :standard_supersaturation))

        local fallback_to_standard_supersaturation_limiter::Bool 
        if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
            @assert nonequilibrium_moisture_sources_limiter_type ∈ [:None, :standard_supersaturation,] "Relaxation to equilibrium only supports no integrator limiter or basic integrator limiter"
        end
        
        moisture_sources_limiter = if nonequilibrium_moisture_sources_limiter_type === :none
            NoMoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :standard_supersaturation
            StandardSupersaturationMoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :basic
            BasicMoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :truncated_basic
            TruncatedBasicMoistureSourcesLimiter(factor)
        elseif nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style
            MorrisonMilbrandt2015MoistureSourcesLimiter()
        elseif nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style_exponential_part_only
            fallback_to_standard_supersaturation_limiter = parse_namelist(namelist, "user_args", "fallback_to_standard_supersaturation_limiter"; default = false)
            MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter(fallback_to_standard_supersaturation_limiter)
        else
            error("Invalid nonequilibrium_moisture_sources_limiter_type: $nonequilibrium_moisture_sources_limiter_type, valid options are :none, :standard_supersaturation, :morrison_milbrandt_2015_style, :morrison_milbrandt_2015_style_exponential_part_only")
        end

        fallback_nonequilibrium_moisture_sources_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "fallback_nonequilibrium_moisture_sources_limiter_type"; default = :none))  # Default to none bc that's purest (basic is how anna wrote it though)
        if nonequilibrium_moisture_scheme_type === :relax_to_equilibrium
            @assert fallback_nonequilibrium_moisture_sources_limiter_type ∈ [:None, :standard_supersaturation,] "Relaxation to equilibrium only supports no integrator limiter or basic integrator limiter"
        end

        fallback_moisture_sources_limiter = if fallback_nonequilibrium_moisture_sources_limiter_type === :none
            NoMoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :basic
            BasicMoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :truncated_basic
            TruncatedBasicMoistureSourcesLimiter(factor)
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :standard_supersaturation
            StandardSupersaturationMoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style
            MorrisonMilbrandt2015MoistureSourcesLimiter()
        elseif fallback_nonequilibrium_moisture_sources_limiter_type === :morrison_milbrandt_2015_style_exponential_part_only
            fallback_to_standard_supersaturation_limiter = parse_namelist(namelist, "user_args", "fallback_to_standard_supersaturation_limiter"; default = false)
            MorrisonMilbrandt2015ExponentialPartOnlyMoistureSourcesLimiter(fallback_to_standard_supersaturation_limiter)
        else
            error("Invalid fallback_nonequilibrium_moisture_sources_limiter_type: $fallback_nonequilibrium_moisture_sources_limiter_type, valid options are :none, :basic, :standard_supersaturation, :morrison_milbrandt_2015_style, :morrison_milbrandt_2015_style_exponential_part_only")
        end

    else
        error("Something went wrong. Invalid moisture model: `$moisture_model`")
    end
    

    # ------------------ Entrainment/Detrainment Limiter ------------------ #
    entr_detr_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "entr_detr_limiter_type"; default = :none))  # Default to none bc that's purest (basic is how anna wrote it though)
    entr_detr_tendency_limiter = if entr_detr_limiter_type === :none
        NoTendencyLimiter()
    elseif entr_detr_limiter_type === :basic
        BasicTendencyLimiter()
    elseif entr_detr_limiter_type === :truncated_basic
        TruncatedBasicTendencyLimiter(factor)
    else
        error("Invalid entr_detr_limiter_type: $entr_detr_limiter_type, valid options are :none, :basic")
    end

    fallback_entr_detr_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "fallback_entr_detr_limiter_type"; default = :none))  # Default to none bc that's purest (basic is how anna wrote it though)
    fallback_entr_detr_tendency_limiter = if fallback_entr_detr_limiter_type === :none
        NoTendencyLimiter()
    elseif fallback_entr_detr_limiter_type === :basic
        BasicTendencyLimiter()
    elseif fallback_entr_detr_limiter_type === :truncated_basic
        TruncatedBasicTendencyLimiter(factor)
    else
        error("Invalid fallback_entr_detr_limiter_type: $fallback_entr_detr_limiter_type, valid options are :none, :basic")
    end
    
    # ------------------ Precipitation Limiter ------------------ #
    precipitation_tendency_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "precipitation_tendency_limiter_type"; default = :none))  # Default to none bc that's purest (basic is how anna wrote it though)
    precipitation_tendency_limiter = if precipitation_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif precipitation_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    elseif precipitation_tendency_limiter_type === :truncated_basic
        TruncatedBasicTendencyLimiter(factor)
    else
        error("Invalid precipitation_tendency_limiter_type: $precipitation_tendency_limiter_type, valid options are :none, :basic")
    end

    fallback_precipitation_tendency_limiter_type::Symbol = Symbol(parse_namelist(namelist, "user_args", "fallback_precipitation_tendency_limiter_type"; default = :none))
    fallback_precipitation_tendency_limiter = if fallback_precipitation_tendency_limiter_type === :none
        NoTendencyLimiter()
    elseif fallback_precipitation_tendency_limiter_type === :basic
        BasicTendencyLimiter()
    elseif fallback_precipitation_tendency_limiter_type === :truncated_basic
        TruncatedBasicTendencyLimiter(factor)
    else
        error("Invalid fallback_precipitation_tendency_limiter_type: $fallback_precipitation_tendency_limiter_type, valid options are :none, :basic")
    end


    # ------------------ Tendency Resolver ------------------ #
    tendency_resolver_setup::Symbol = Symbol(parse_namelist(namelist, "user_args", "tendency_resolver_setup"; default = :none)) # if you do it just on entr/det's it's funamentally unstable....
    if tendency_resolver_setup == :none
        use_tendency_resolver = false
        use_tendency_resolver_on_full_tendencies = false
    elseif tendency_resolver_setup == :full_tendencies
        use_tendency_resolver = true
        use_tendency_resolver_on_full_tendencies = true
    elseif tendency_resolver_setup == :normal
        use_tendency_resolver = true
        use_tendency_resolver_on_full_tendencies = false
    else
        error("tendency_resolver_setup must be :none, :full_tendencies, or :normal, not $tendency_resolver_setup")
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
        #
        use_tendency_resolver,
        use_tendency_resolver_on_full_tendencies,
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


# =============== #
abstract type AbstractQuadratureType end
struct LogNormalQuad <: AbstractQuadratureType end
struct GaussianQuad <: AbstractQuadratureType end

abstract type AbstractEnvThermo end
struct SGSMean <: AbstractEnvThermo end

abstract type AbstractSGSQuadratureType end
struct SGSQuadrature{N, QT, A, W} <: AbstractSGSQuadratureType
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


struct SGSMeanWQuadratureAdjustedNoneqMoistureSources{N, QT, A, W} <: AbstractSGSQuadratureType
    quadrature_type::QT
    a::A
    w::W
    function SGSMeanWQuadratureAdjustedNoneqMoistureSources(::Type{FT}, namelist) where {FT}
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
quadrature_order(::SGSMeanWQuadratureAdjustedNoneqMoistureSources{N}) where {N} = N
quad_type(::SGSMeanWQuadratureAdjustedNoneqMoistureSources{N}) where {N} = N
# =============== #

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
struct EDMFModel{N_up, FT, SFCA, SABC, MM <: AbstractMoistureModel, CSM <: AbstractCloudSedimentationModel, TCM, PM, RFM, SFM, PFM, ENT, EBGC, MLP, PMP, EC, MLEC, ET, EDS, DDS, EPG, TLT <: AbstractTendencyLimiterSet, APM <: AbstractAreaPartitionModel, CTKE <: AbstractConvectiveTKEHandler}
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
    area_partition_model::APM
    convective_tke_handler::CTKE
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
        EquilibriumMoisture(param_set, namelist)
    elseif moisture_model_name == "nonequilibrium"
        NonEquilibriumMoisture(param_set, namelist)
    else
        error("Something went wrong. Invalid moisture model: '$moisture_model_name'")
    end

    snow_formation_model = SnowFormationModel(param_set, moisture_model, namelist) # rn this has to go after moisture_model, if we need it to go before later we can go back to dispatching on symbols etc...

    use_sedimentation = parse_namelist(namelist, "user_args", "use_sedimentation"; default = false)
    cloud_sedimentation_model = use_sedimentation ? CloudSedimentationModel(param_set, moisture_model, namelist) : CloudNoSedimentationModel()

    tendency_limiters = TendencyLimiterSet(param_set, moisture_model, namelist)

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
        parse_namelist(namelist, "thermodynamics", "sgs"; default = "mean", valid_options = ["mean", "quadrature", "mean_w_quadrature_adjusted_noneq_moisture_sources"])
    en_thermo = if en_sgs_name == "mean"
        SGSMean()
    elseif en_sgs_name == "quadrature"
        SGSQuadrature(FT, namelist)
    elseif en_sgs_name == "mean_w_quadrature_adjusted_noneq_moisture_sources"
        SGSMeanWQuadratureAdjustedNoneqMoistureSources(FT, namelist) # special case for nonequilibrium moisture with quadrature thermo, we do mean thermo but quadrature adjusted non-eq sources
    else
        error("Something went wrong. Invalid environmental sgs type '$en_sgs_name'")
    end
    bg_closure = if en_sgs_name == "mean"
        BuoyGradMean()
    elseif en_sgs_name == "quadrature"
        BuoyGradQuadratures()
    elseif en_sgs_name == "mean_w_quadrature_adjusted_noneq_moisture_sources"
        BuoyGradMean() # same as quadrature, we just do mean thermo but quadrature adjusted non-eq sources
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


    # -- Stalled updraft handler ----------------------------------------------------- #

    stalled_updraft_handler = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "stalled_updraft_handler", default = "kill", valid_options = ["none", "mix_to_grid_mean", "kill", "detrain_downdrafts"]) # kill bc that was original default
    stalled_updraft_handler = if stalled_updraft_handler == "none"
        StalledUpdraftDoNothing()
    elseif stalled_updraft_handler == "mix_to_grid_mean"
        StalledUpdraftMixToGridMean()
    elseif stalled_updraft_handler == "kill"
        StalledUpdraftKill()
    elseif stalled_updraft_handler == "detrain_downdrafts"
        StalledUpdraftDetrainDowndrafts()
    else
        error("Something went wrong. Invalid stalled updraft handler '$stalled_updraft_handler'")
    end
    # ---------------------------------------------------------------------------------- #

    # -- Convective TKE ---------------------------------------------------------------- #

    use_convective_tke = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "use_convective_tke"; default = false, valid_options = [true, false])
    convective_tke = use_convective_tke ? ConvectiveTKE(param_set, namelist) : NoConvectiveTKE()


    # -- Entrainment Type --------------------------------------------------------------- #


    entrainment_type = parse_namelist(
        namelist,
        "turbulence",
        "EDMF_PrognosticTKE",
        "entrainment_type";
        default = "fractional",
        valid_options = ["fractional", "total_rate"],
    )


    if entrainment_type == "fractional"
        # entrainment_type = FractionalEntrModel()
        entrainment_type = FractionalEntrModel(;
        stalled_updraft_handler = stalled_updraft_handler,
        )
    elseif entrainment_type == "total_rate"
        # entrainment_type = TotalRateEntrModel()
        entrainment_type = TotalRateEntrModel(;
        base_detrainment_rate_inv_s = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "base_detrainment_rate_inv_s", default = 0.0), # assume infinite timescale so 0 inv timescale (consider choosing a better default like 6 hours or something...)
        base_entrainment_rate_inv_s = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "base_entrainment_rate_inv_s", default = 0.0), # assume infinite timescale so 0 inv timescale (consider choosing a better default like 6 hours or something...)
        stalled_updraft_handler = stalled_updraft_handler, 
        )
    else
        error("Something went wrong. Invalid entrainment type '$entrainment_type'")
    end
    @info "entrainment_type: $entrainment_type"


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

    # -- Area partition model ----------------------------------------------------- #

    area_partition_model = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "area_partition_model", default = "standard", valid_options = ["standard", "core_cloak"],) # standard bc that was original default
    area_partition_model = if area_partition_model == "standard"
        StandardAreaPartitionModel(param_set, namelist)
    elseif area_partition_model == "core_cloak"
        CoreCloakAreaPartitionModel(param_set, namelist)
    else
        error("Something went wrong. Invalid area partition model '$area_partition_model'")
    end

    # ================================================================================= #


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
    APM = typeof(area_partition_model)
    CTKE = typeof(convective_tke)

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
        TLT,
        APM,
        CTKE,
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
        area_partition_model,
        convective_tke,
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


struct State{P, A, T, G <: Grid}
    prog::P
    aux::A
    tendencies::T
    grid::G
    calibrate_io::Bool # put it in here so we don't have to carry it around...
end

function State(prog::P, aux::A, tendencies::T, calibrate_io::Bool) where {P, A, T}
    grid = Grid(axes(prog.cent))
    return State{P, A, T, typeof(grid)}(prog, aux, tendencies, grid, calibrate_io)
end

"""
    column_state(prog, aux, tendencies, colidx)

Create a columnar state given full 3D states
 - `prog` prognostic state
 - `aux` auxiliary state
 - `tendencies` tendencies state
 - `colidx` column index, from ClimaCore's `bycolumn` function

## Example
```
julia
bycolumn(axes(prog.cent)) do colidx
    state = TC.column_state(prog, aux, tendencies, colidx, calibrate_io)
    ...
end
```
"""
# function column_state(prog, aux, tendencies, colidx, calibrate_io::Bool)
#     prog_cent_column = CC.column(prog.cent, colidx)
#     prog_face_column = CC.column(prog.face, colidx)
#     aux_cent_column = CC.column(aux.cent, colidx)
#     aux_face_column = CC.column(aux.face, colidx)
#     tends_cent_column = CC.column(tendencies.cent, colidx)
#     tends_face_column = CC.column(tendencies.face, colidx)
#     prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
#     aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)
#     tends_column = CC.Fields.FieldVector(cent = tends_cent_column, face = tends_face_column)

#     return State(prog_column, aux_column, tends_column, calibrate_io)
# end

function column_state(prog, aux, tendencies, colidx, calibrate_io::Bool)
    # USE BRACKETS [colidx]
    prog_cent_column = prog.cent[colidx]
    prog_face_column = prog.face[colidx]
    aux_cent_column  = aux.cent[colidx]
    aux_face_column  = aux.face[colidx]
    tends_cent_column = tendencies.cent[colidx]
    tends_face_column = tendencies.face[colidx]
    
    prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)
    tends_column = CC.Fields.FieldVector(cent = tends_cent_column, face = tends_face_column)

    return State(prog_column, aux_column, tends_column, calibrate_io)
end

# function column_prog_aux(prog, aux, colidx, calibrate_io::Bool)
#     prog_cent_column = CC.column(prog.cent, colidx)
#     prog_face_column = CC.column(prog.face, colidx)
#     aux_cent_column = CC.column(aux.cent, colidx)
#     aux_face_column = CC.column(aux.face, colidx)
#     prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
#     aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)

#     return State(prog_column, aux_column, nothing, calibrate_io)
# end


"""
    column_prog_aux(prog, aux, colidx, calibrate_io)

Creates a local state vector for a single column.

# ClimaCore 0.14+ Update
In ClimaCore 0.14, using `CC.column(field, colidx)` is discouraged for 1D physics broadcasting because it returns a 3D `VIJFH` layout view, which causes `MethodError: no method matching copyto!(::VIJFH, ...)` during broadcasting.

Instead, use `field[colidx]`. This invokes the specific `getindex` method that projects the DataLayout from 3D (`VIJFH`) down to 1D (`VF`), enabling standard broadcasting and kernel optimization.

# Example
# OLD (ClimaCore 0.13):
# prog_cent_column = CC.column(prog.cent, colidx)

# NEW (ClimaCore 0.14):
prog_cent_column = prog.cent[colidx]
"""
function column_prog_aux(prog, aux, colidx, calibrate_io::Bool)
    # USE BRACKETS [colidx] to get a 1D VF (Vertical Field) view
    # CC.column(...) returns a 3D VIJFH view, which crashes the broadcaster
    
    prog_cent_column = prog.cent[colidx]
    prog_face_column = prog.face[colidx]
    aux_cent_column  = aux.cent[colidx]
    aux_face_column  = aux.face[colidx]

    # Reconstruct FieldVectors using these 1D views
    prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)

    return State(prog_column, aux_column, nothing, calibrate_io)
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


function full_raw_string(f::CC.Fields.Field)
    # No truncation
    # Use parent to get the raw array and join elements horizontally
    raw_array = parent(f)

    out = "$(eltype(f))-valued Field:  ["

    out *= string(Base.join((raw_array), " ")) * "]"

    return out
end


function fv_field_summary_helper(f::CC.Fields.Field)
    pns = propertynames(f)

    if iszero(length(pns)) # no properties
        # Base.show(io, f)
        # s = repr(f)
        s = full_raw_string(f) # no truncation
        # s = replace(s, "\n" => " \n     ") # otherwise @info smushes onto same line, not sure why"
        # print(io, s)
        return s
    else
        s = map(pn -> (pn, fv_field_summary_helper(getproperty(f, pn))), pns)
        # s = map(pn -> string("`", pn, "`::", fv_field_summary_helper(getproperty(f, pn))), pns)
        # s = join(s, ",\n ")
        return s
    end
end


function do_print(io, fv_strings; indent::Int=0, return_str::Bool=false)
    function action(io, xargs..., ; s, return_str::Bool=return_str)
        if return_str
            return s * string(xargs...)
        else
            print(io, xargs...)
            return string() # for type stability, don't return `nothing`
        end
    end
    
    big_space_and_newline = "\n              " # for some reason, this triggers nicer output w/ @info
    s::String = "" # you have to start with this to get the nice formatting w/ @info, otherwise it just prints in a long line
    if fv_strings isa String
        s = action(io, fv_strings, big_space_and_newline; s=s)
    elseif fv_strings isa Tuple || s isa Base.Generator
        for (pn, psummary) in fv_strings
            s = action(io, " "^indent; s=s)
            s = action(io, "`", pn, "`::",; s=s)
            if !(psummary isa String)
                s = action(io, big_space_and_newline; s=s)
            end

            if return_str
                s = s * do_print(io, psummary; indent = indent + 2, return_str = return_str)
            else
                do_print(io, psummary; indent = indent + 2, return_str = return_str)
            end
        end
    else
        error("Unexpected type in fv_field_summary: $(typeof(s))")
    end
    if return_str
        return s
    else
        return string() # for type stability, don't return `nothing`
    end

end

function fv_field_summary(io::IO, f::CC.Fields.Field)
    s = fv_field_summary_helper(f)
    do_print(io, s)
end

function fv_field_summary(f::CC.Fields.Field)
    s = fv_field_summary_helper(f);
    big_space_and_newline = "\n              " # for some reason, this triggers nicer output w/ @info
    s =  big_space_and_newline  * do_print(stdout, s; return_str = true) # add leading big space to outermost call
    return s
end
    

function Base.summary(io::IO, f::CC.Fields.Field)
    fv_field_summary(io, f)
end

function Base.summary(f::CC.Fields.Field) 
    fv_field_summary(f)
end

function Base.summary(io::IO, fv::CC.Fields.FieldVector)

    pns = string.(propertynames(fv))

    if length(pns) == 0
        print(io, fv)
        return
    end
    
    for pn in propertynames(fv)[1:2]
        prop = getproperty(fv, pn)

        if prop isa CC.Fields.Field
           fv_field_summary(io, prop)
        elseif prop isa CC.Fields.FieldVector
            Base.summary(io, prop)
        else
            error("Unexpected type in FieldVector summary: $(typeof(prop))")
        end

    end
end


function Base.summary(fv::CC.Fields.FieldVector)
    s = ""
    pns = string.(propertynames(fv))

    if length(pns) == 0
        s = string(fv)
        return s
    end

    for pn in propertynames(fv)[1:2]
        prop = getproperty(fv, pn)

        if prop isa CC.Fields.Field
            s *= fv_field_summary(prop)
        elseif prop isa CC.Fields.FieldVector
            s *= Base.summary(prop)
        else
            error("Unexpected type in FieldVector summary: $(typeof(prop))")
        end

    end
    return s
end