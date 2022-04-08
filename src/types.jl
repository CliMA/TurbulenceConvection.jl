"""
    PrecipFormation

Storage for tendencies due to precipitation formation

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrecipFormation{FT}
    θ_liq_ice_tendency::FT
    qt_tendency::FT
    qr_tendency::FT
    qs_tendency::FT
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
    "Turbulent entrainment"
    ε_turb::FT
    "Nondimensional fractional dynamical entrainment"
    ε_nondim::FT
    "Nondimensional fractional dynamical detrainment"
    δ_nondim::FT
end

"""
    GeneralizedEntr

A general set of variables entrainment might depend on.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct GeneralizedEntr{FT}
    "updraft condensate (liquid water + ice)"
    q_cond_up::FT
    "environment condensate (liquid water + ice)"
    q_cond_en::FT
    "updraft vertical velocity"
    w_up::FT
    "environment vertical velocity"
    w_en::FT
    "updraft buoyancy"
    b_up::FT
    "environment buoyancy"
    b_en::FT
    "grid mean tke"
    tke_gm::FT
    "environment tke"
    tke_en::FT
    "updraft momentum divergence"
    dMdz::FT
    "updraft momentum"
    M::FT
    "updraft area fraction"
    a_up::FT
    "environment area fraction"
    a_en::FT
    "plume scale height"
    H_up::FT
    "updraft relative humidity"
    RH_up::FT
    "environment relative humidity"
    RH_en::FT
    "maximum updraft area"
    max_area::FT
    "vertical coordinate"
    zc_i::FT
    "Model time step"
    Δt::FT
    # For stochastic modelling
    "nondimensional fractional dynamical entrainment"
    ε_nondim::FT
    "nondimensional fractional dynamical detrainment"
    δ_nondim::FT
    "convective velocity"
    wstar::FT
end

abstract type AbstractEntrainmentDimScaleParameters end

struct BuoyVelScaleParameters{FT} <: AbstractEntrainmentDimScaleParameters
    w_min::FT
    c_λ::FT
end
function BuoyVelScaleParameters(
    param_set)
    aliases = ["w_min", "c_λ",]
    (w_min, c_λ) = CLIMAParameters.get_parameter_values!(
        param_set,
        aliases,
        "BuoyVelScale"
    )
    return BuoyVelScaleParameters{get_parametric_type(param_set)}(
        w_min,
        c_λ
    )
end

struct InvZScaleParameters{FT} <: AbstractEntrainmentDimScaleParameters end


Base.eltype(::GeneralizedEntr{FT}) where {FT} = FT

abstract type AbstractEntrDetrModel end
abstract type AbstractNonLocalEntrDetrModel end
struct MDEntr <: AbstractEntrDetrModel end  # existing model
struct NNEntr <: AbstractEntrDetrModel end
struct NNEntrNonlocal <: AbstractNonLocalEntrDetrModel end
struct LinearEntr <: AbstractEntrDetrModel end
struct FNOEntr <: AbstractNonLocalEntrDetrModel end
struct RFEntr <: AbstractEntrDetrModel end

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
function MDEntrainmentParameters(
    param_set)

    aliases = ["w_min","c_ε","μ_0","β","χ","c_δ"]
    (w_min,c_ε,μ_0,β,χ,c_δ) = CLIMAParameters.get_parameter_values!(
        param_set,
        aliases,
        "MDEntrainment"
    )
    

    return MDEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(
        w_min,
        c_ε,
        μ_0,
        β,
        χ,
        c_δ,
    )

end

struct FNOEntrainmentParameters{FT} <: AbstractEntrainmentClosureNonlocalParameters
        c_fno::AbstractVector{FT}
end
function FNOEntrainmentParameters(
    param_set)

    aliases = ["c_fno"]
    (c_fno) = CLIMAParameters.get_parameter_array_values!(
        param_set,
        aliases,
        "FNOEntrainment"
    )
    

    return FNOEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(
        c_fno
    )

end

struct NNEntrainmentNonlocalParameters{FT} <: AbstractEntrainmentClosureNonlocalParameters
        c_gen::AbstractVector{FT}
end
function NNEntrainmentNonlocalParameters(
    param_set)

    array_aliases = ["c_gen"]
    (c_gen) = CLIMAParameters.get_parameter_array_values!(
        param_set,
        array_aliases,
        "NNEntrainmentNonlocal"
    )
    

    return NNEntrainmentNonlocalParameters{CLIMAParameters.get_parametric_type(param_set)}(
        c_gen,
    )

end

struct NNEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_gen::AbstractVector{FT}
end
function NNEntrainmentParameters(
    param_set)
    aliases = ["w_min"]
    (w_min) = CLIMAParameters.get_parameter_values!(
        param_set,
        aliases,
        "NNEntrainment"
    )
    array_aliases = ["c_gen"]
    (c_gen) = CLIMAParameters.get_parameter_values!(
        param_set,
        array_aliases,
        "NNEntrainment"
    )
    

    return NNEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(
        w_min,
        c_gen,
    )

end

struct LinearEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_gen::AbstractVector{FT}
end
function LinearEntrainmentParameters(
    param_set)
    aliases = ["w_min"]
    (w_min) = CLIMAParameters.get_parameter_values!(
        param_set,
        aliases,
        "LinearEntrainment"
    )
    array_aliases = ["c_gen"]
    (c_gen) = CLIMAParameters.get_parameter_values!(
        param_set,
        array_aliases,
        "LinearEntrainment"
    )
    

    return LinearEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(
        w_min,
        c_gen,
    )

end


struct RFEntrainmentParameters{FT} <: AbstractEntrainmentClosureParameters
    w_min::FT
    c_rf_fix::AbstractVector{FT}
    c_rf_opt::AbstractVector{FT}
end
function RFEntrainmentParameters(
    param_set)
    aliases = ["w_min"]
    (w_min) = CLIMAParameters.get_parameter_values!(
        param_set,
        aliases,
        "RFEntrainment"
    )
    array_aliases = ["c_rf_fix","c_rf_opt"]
    (c_rf_fix,c_rf_opt) = CLIMAParameters.get_parameter_values!(
        param_set,
        array_aliases,
        "RFEntrainment"
    )
    

    return RFEntrainmentParameters{CLIMAParameters.get_parametric_type(param_set)}(
        w_min,
        c_rf_fix,
        c_rf_opt,
    )

end



Base.@kwdef struct NoisyRelaxationProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end
Base.@kwdef struct LogNormalScalingProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end

Base.@kwdef struct PrognosticNoisyRelaxationProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end

#parameters for the stochastic types
abstract type AbstractStochasticEntrainmentClosureParameters end

struct LogNormalScalingProcessParameters{FT, AECPS} <: AbstractStochasticEntrainmentClosureParameters
    c_gen_stoch::AbstractVector{FT}
    ECPS::AECPS
end
function LogNormalScalingProcessParameters(
    param_set,
    ECPS::AECPS) where {AECPS <: AbstractEntrainmentClosureParameters}
    array_aliases = ["c_gen_stoch"]
    (c_gen_stoch) = CLIMAParameters.get_parameter_values!(
        param_set,
        array_aliases,
        "LogNormalScalingProcess"
    )

    return LogNormalScalingProcessParameters{CLIMAParameters.get_parametric_type(param_set), AECPS}(
        c_gen_stoch,
        ECPS,
    )

end

struct NoisyRelaxationProcessParameters{FT, AECPS} <: AbstractStochasticEntrainmentClosureParameters
    c_gen_stoch::AbstractVector{FT}
    ECPS::AECPS
end
function NoisyRelaxationProcessParameters(
    param_set,
    ECPS::AECPS) where {AECPS <: AbstractEntrainmentClosureParameters}
    array_aliases = ["c_gen_stoch"]
    (c_gen_stoch) = CLIMAParameters.get_parameter_values!(
        param_set,
        array_aliases,
        "NoisyRelaxationProcess"
    )

    return NoisyRelaxationProcessParameters{CLIMAParameters.get_parametric_type(param_set), AECPS}(
        c_gen_stoch,
        ECPS,
    )
end

struct PrognosticNoisyRelaxationProcessParameters{FT, AECPS} <: AbstractStochasticEntrainmentClosureParameters
    ECPS::AECPS
end
function PrognosticNoisyRelaxationProcessParameters(
    param_set,
    ECPS::AECPS) where {AECPS <: AbstractEntrainmentClosureParameters}

    return PrognosticNoisyRelaxationProcessParameters{CLIMAParameters.get_parametric_type(param_set), AECPS}(
        ECPS,
    )
end



abstract type EntrDimScale end
struct BuoyVelEntrDimScale <: EntrDimScale end
struct InvZEntrDimScale <: EntrDimScale end

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
    p0::FT
    "cloud fraction"
    en_cld_frac::FT
    "specific volume"
    alpha0::FT
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
    p0::FT
    "vertical buoyancy gradient struct"
    ∇b::GradBuoy{FT}
    "env shear"
    Shear²::FT
    "environment turbulent kinetic energy"
    tke::FT
    "Updraft tke source"
    b_exch::FT
end

abstract type AbstractPrecipitationModel end
struct NoPrecipitation <: AbstractPrecipitationModel end
struct CutoffPrecipitation <: AbstractPrecipitationModel end
struct Clima1M <: AbstractPrecipitationModel end





abstract type AbstractQuadratureType end
struct LogNormalQuad <: AbstractQuadratureType end

abstract type AbstractEnvThermo end
struct SGSMean <: AbstractEnvThermo end
struct SGSQuadrature{N, QT, A, W} <: AbstractEnvThermo
    quadrature_type::QT
    a::A
    w::W
    function SGSQuadrature(namelist)
        N = parse_namelist(namelist, "thermodynamics", "quadrature_order"; default = 3)
        quadrature_name = parse_namelist(namelist, "thermodynamics", "quadrature_type"; default = "gaussian")
        quadrature_type = if quadrature_name == "log-normal"
            LogNormalQuad()
        else
            error("Invalid thermodynamics quadrature $(quadrature_name)")
        end
        # TODO: double check this python-> julia translation
        # a, w = np.polynomial.hermite.hermgauss(N)
        a, w = FastGaussQuadrature.gausshermite(N)
        a, w = SA.SVector{N}(a), SA.SVector{N}(w)
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
) where {FT, FVT}
    TS = typeof(Tsurface)
    QS = typeof(qsurface)
    CH = typeof(ch)
    CM = typeof(cm)
    return FixedSurfaceCoeffs{FT, TS, QS, CH, CM}(; Tsurface, qsurface, ch, cm, kwargs...)
end

Base.@kwdef struct MoninObukhovSurface{FT, TS, QS, SHF, LHF} <: AbstractSurfaceParameters{FT}
    zrough::FT = FT(0)
    Tsurface::TS = FT(0)
    qsurface::QS = FT(0)
    shf::SHF = FT(0)
    lhf::LHF = FT(0)
    Ri_bulk_crit::FT = FT(0)
end

function MoninObukhovSurface(
    ::Type{FT};
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
    return MoninObukhovSurface{FT, TS, QS, SHF, LHF}(; Tsurface, qsurface, shf, lhf, kwargs...)
end

Base.@kwdef struct SullivanPattonSurface{FT, TS, QS, SHF, LHF} <: AbstractSurfaceParameters{FT}
    zrough::FT = FT(0)
    Tsurface::TS = FT(0)
    qsurface::QS = FT(0)
    shf::SHF = FT(0)
    lhf::LHF = FT(0)
    cq::FT = FT(0)
    Ri_bulk_crit::FT = FT(0)
    ustar::FT = FT(0)
end

function SullivanPattonSurface(
    ::Type{FT};
    Tsurface::FloatOrFunc{FT},
    qsurface::FloatOrFunc{FT},
    shf::FloatOrFunc{FT},
    lhf::FloatOrFunc{FT},
    kwargs...,
) where {FT}
    TS = typeof(Tsurface)
    QS = typeof(qsurface)
    SHF = typeof(shf)
    LHF = typeof(lhf)
    return SullivanPattonSurface{FT, TS, QS, SHF, LHF}(; Tsurface, qsurface, shf, lhf, kwargs...)
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
    ρθ_liq_ice_flux::FT = 0
    ρu_flux::FT = 0
    ρv_flux::FT = 0
    obukhov_length::FT = 0
    wstar::FT = 0
end

struct ForcingNone end
struct ForcingStandard end
struct ForcingDYCOMS_RF01 end
struct ForcingLES end

struct RadiationNone end
struct RadiationDYCOMS_RF01 end
struct RadiationLES end

Base.@kwdef struct LESData
    "Start time index of LES"
    imin::Int = 0
    "End time index of LES"
    imax::Int = 0
    "Path to LES stats file used to drive SCM"
    les_filename::String = nothing
    "Drive SCM with LES data from t = [end - t_interval_from_end_s, end]"
    t_interval_from_end_s::Float64 = 6 * 3600.0
    "Length of time to average over for SCM initialization"
    initial_condition_averaging_window_s::Float64 = 3600.0
end

"""
    ForcingBase

LES-driven forcing

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ForcingBase{T, R}
    "Boolean specifying whether Coriolis forcing is applied"
    apply_coriolis::Bool = false
    "Boolean specifying whether subsidence forcing is applied"
    apply_subsidence::Bool = false
    "Coriolis parameter"
    coriolis_param::Float64 = 0
    "Momentum relaxation timescale"
    nudge_tau::Float64 = 0.0
    "Radiative forcing"
    rad::R
end

ForcingBase(::Type{T}; rad = nothing, kwargs...) where {T} = ForcingBase{T, typeof(rad)}(; rad, kwargs...)

force_type(::ForcingBase{T}) where {T} = T

Base.@kwdef struct RadiationBase{T}
    divergence::Float64 = 0
    alpha_z::Float64 = 0
    kappa::Float64 = 0
    F0::Float64 = 0
    F1::Float64 = 0
end

rad_type(::RadiationBase{T}) where {T} = T

abstract type AbstractInversion end
struct CriticalRiInversion <: AbstractInversion end
struct max∇θInversion <: AbstractInversion end
struct θρInversion <: AbstractInversion end

Base.@kwdef struct CasesBase{T, IT, SURFP, F, R, LESDataT}
    case::T
    casename::String
    inversion_type::IT
    surf_params::SURFP
    Fo::F
    Rad::R
    LESDat::LESDataT
end

function CasesBase(case::T; inversion_type, surf_params, Fo, Rad, LESDat = nothing, kwargs...) where {T}
    F = typeof(Fo)
    R = typeof(Rad)
    IT = typeof(inversion_type)
    SURFP = typeof(surf_params)
    LESDataT = typeof(LESDat)
    CasesBase{T, IT, SURFP, F, R, LESDataT}(;
        case = case,
        casename = string(nameof(T)),
        inversion_type,
        surf_params,
        Fo,
        Rad,
        LESDat,
        kwargs...,
    )
end

struct EDMFParameters{FT,AECPS,AEDSPS}
    ECPS::AECPS
    EDSPS::AEDSPS
end
function EDMFParameters(
    param_set,
    ECPS::AECPS,
    EDSPS::AEDSPS
) where {AECPS <: Union{AbstractEntrainmentClosureParameters,AbstractStochasticEntrainmentClosureParameters}, AEDSPS <: AbstractEntrainmentDimScaleParameters}

    aliases = []

    () = CLIMAParameters.get_parameter_values!(param_set,aliases,"EDMF")
    
    return EDMFParameters{get_parametric_type(param_set),AECPS,AEDSPS}(
        ECPS,
        EDSPS,
    )
    
end



struct EDMFModel{N_up, FT, PM, ENT, EBGC, EC, ECP, EDS, EDSP}
    surface_area::FT
    max_area::FT
    minimum_area::FT
    precip_model::PM
    en_thermo::ENT
    prandtl_number::FT
    bg_closure::EBGC
    entr_closure::EC
    entr_closure_parameters::ECP
    entr_dim_scale::EDS
    entr_dim_scale_parameters::EDSP
    function EDMFModel(namelist, precip_model) where {PS}
        # TODO: move this into arg list
        FT = Float64
        # get values from namelist
        prandtl_number = namelist["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"]

        # Set the number of updrafts (1)
        n_updrafts = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "updraft_number"; default = 1)

        pressure_func_drag_str = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "pressure_closure_drag";
            default = "normalmode",
            valid_options = ["normalmode", "normalmode_signdf"],
        )

        surface_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "surface_area"; default = 0.1)
        max_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "max_area"; default = 0.9)
        minimum_area = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "min_area"; default = 1e-5)

        precip_model = precip_model
        # Create the environment variable class (major diagnostic and prognostic variables)

        # Create the class for environment thermodynamics
        en_thermo_name = parse_namelist(namelist, "thermodynamics", "sgs"; default = "mean")
        en_thermo = if en_thermo_name == "mean"
            SGSMean()
        elseif en_thermo_name == "quadrature"
            SGSQuadrature(namelist)
        end

        bg_type = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "env_buoy_grad";
            default = "mean",
            valid_options = ["mean", "quadratures"],
        )
        bg_closure = if bg_type == "mean"
            BuoyGradMean()
        elseif bg_type == "quadratures"
            BuoyGradQuadratures()
        else
            error("Something went wrong. Invalid environmental buoyancy gradient closure type '$bg_type'")
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
            valid_options = ["moisture_deficit", "NN", "NN_nonlocal", "FNO", "Linear", "RF"],
        )

        (mean_entr_closure, mean_entr_parameters) = if entr_type == "moisture_deficit"
            (MDEntr(), MDEntrainmentParameters(param_set))
        elseif entr_type == "NN"
            (NNEntr(), NNEntrainmentParameters(param_set))
        elseif entr_type == "NN_nonlocal"
            (NNEntrNonlocal(), NNEntrainmentNonlocalParameters(param_set))
        elseif entr_type == "FNO"
            (FNOEntr(), FNOEntrainmentParameters(param_set))
        elseif entr_type == "Linear"
            (LinearEntr(), LinearEntrainmentParameters(param_set))
        elseif entr_type == "RF"
            (RFEntr(), RFEntrainmentParameters(param_set))
        else
            error("Something went wrong. Invalid entrainment type '$entr_type'")
        end

        # Overwrite `entr_closure` if a noisy relaxation process is used
        (entr_closure, entr_closure_parameters) = if stoch_entr_type == "noisy_relaxation_process"
            (NoisyRelaxationProcess(mean_model = mean_entr_closure),
             NoisyRelaxationProcessParameters(param_set, mean_entr_parameters))
        elseif stoch_entr_type == "lognormal_scaling"
            (LogNormalScalingProcess(mean_model = mean_entr_closure),
             LogNormalScalingProcessParameters(param_set, mean_entr_parameters))
        elseif stoch_entr_type == "prognostic_noisy_relaxation_process"
            (PrognosticNoisyRelaxationProcess(mean_model = mean_entr_closure),
             PrognosticNoisyRelaxationProcessParameters(param_set, mean_entr_parameters))
        elseif stoch_entr_type == "deterministic"
            (mean_entr_closure, mean_entr_parameters)
        else
            error("Something went wrong. Invalid stochastic entrainment type '$stoch_entr_type'")
        end

        entr_dim_scale = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "entr_dim_scale";
            default = "buoy_vel",
            valid_options = ["buoy_vel", "inv_z"],
        )

        (entr_dim_scale, entr_dim_scale_parameters) = if entr_dim_scale == "buoy_vel"
            (BuoyVelEntrDimScale(), BuoyVelScaleParameters(param_set))
        elseif entr_dim_scale == "inv_z"
            (InvZEntrDimScale(),  InvZScaleParameters())
        else
            error("Something went wrong. Invalid entrainment dimension scale '$entr_dim_scale'")
        end

        EDSP = typeof(entr_dim_scale_parameters)
        EDS = typeof(entr_dim_scale)
        ECP = typeof(entr_closure_parameters)
        EC = typeof(entr_closure)
        PM = typeof(precip_model)
        EBGC = typeof(bg_closure)
        ENT = typeof(en_thermo)
        return new{n_updrafts, FT, PM, ENT, EBGC, EC, ECP, EDS, EDSP}(
            surface_area,
            max_area,
            minimum_area,
            precip_model,
            en_thermo,
            prandtl_number,
            bg_closure,
            entr_closure,
            entr_closure_parameters,
            entr_dim_scale,
            entr_dim_scale_parameters
        )
    end
end
parameter_set(obj) = obj.param_set
n_updrafts(::EDMFModel{N_up}) where {N_up} = N_up
Base.eltype(::EDMFModel{N_up, FT}) where {N_up, FT} = FT

struct State{P, A, T}
    prog::P
    aux::A
    tendencies::T
end
