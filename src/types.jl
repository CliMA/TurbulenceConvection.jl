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
    "Dynamical entrainment"
    ε_dyn::FT
    "Dynamical detrainment"
    δ_dyn::FT
    "Turbulent entrainment"
    ε_turb::FT
    "nondimensional dynamical entrainment"
    ε_dyn_nondim::FT
    "nondimensional dynamical detrainment"
    δ_dyn_nondim::FT
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
    "environment tke"
    tke::FT
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
    "nondimensional entrainment"
    nondim_entr_sc::FT
    "nondimensional detrainment"
    nondim_detr_sc::FT
end

Base.eltype(::GeneralizedEntr{FT}) where {FT} = FT

abstract type AbstractEntrDetrModel end
abstract type AbstractNonLocalEntrDetrModel end
struct MDEntr <: AbstractEntrDetrModel end  # existing model
struct NNEntr <: AbstractEntrDetrModel end
struct NNEntrNonlocal <: AbstractNonLocalEntrDetrModel end
struct LinearEntr <: AbstractEntrDetrModel end
struct FNOEntr <: AbstractNonLocalEntrDetrModel end

Base.@kwdef struct NoisyRelaxationProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end
Base.@kwdef struct LogNormalScalingProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end

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

struct GridMeanVariables{PS}
    param_set::PS
end

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

abstract type AbstractTurbConvModel end

struct EDMF_PrognosticTKE{N_up, PM, ENT, EBGC, EC} <: AbstractTurbConvModel
    surface_area::Float64
    max_area::Float64
    minimum_area::Float64
    precip_model::PM
    en_thermo::ENT
    prandtl_number::Float64
    bg_closure::EBGC
    entr_closure::EC
    function EDMF_PrognosticTKE(namelist, grid::Grid, param_set::PS) where {PS}
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

        # Get values from namelist
        # set defaults at some point?
        surface_area = namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"]
        max_area = namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"]
        # entrainment parameters

        # Need to code up as namelist option?
        minimum_area = 1e-5

        # Create the class for precipitation

        precip_name = parse_namelist(
            namelist,
            "microphysics",
            "precipitation_model";
            default = "None",
            valid_options = ["None", "cutoff", "clima_1m"],
        )

        precip_model = if precip_name == "None"
            NoPrecipitation()
        elseif precip_name == "cutoff"
            CutoffPrecipitation()
        elseif precip_name == "clima_1m"
            Clima1M()
        else
            error("Invalid precip_name $(precip_name)")
        end

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
            valid_options = ["deterministic", "noisy_relaxation_process", "lognormal_scaling"],
        )
        entr_type = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "entrainment";
            default = "moisture_deficit",
            valid_options = ["moisture_deficit", "NN", "NN_nonlocal", "FNO", "Linear"],
        )
        mean_entr_closure = if entr_type == "moisture_deficit"
            MDEntr()
        elseif entr_type == "NN"
            NNEntr()
        elseif entr_type == "NN_nonlocal"
            NNEntrNonlocal()
        elseif entr_type == "FNO"
            FNOEntr()
        elseif entr_type == "Linear"
            LinearEntr()
        else
            error("Something went wrong. Invalid entrainment type '$entr_type'")
        end

        # Overwrite `entr_closure` if a noisy relaxation process is used
        entr_closure = if stoch_entr_type == "noisy_relaxation_process"
            NoisyRelaxationProcess(mean_model = mean_entr_closure)
        elseif stoch_entr_type == "lognormal_scaling"
            LogNormalScalingProcess(mean_model = mean_entr_closure)
        elseif stoch_entr_type == "deterministic"
            mean_entr_closure
        else
            error("Something went wrong. Invalid stochastic entrainment type '$stoch_entr_type'")
        end

        EC = typeof(entr_closure)
        PM = typeof(precip_model)
        EBGC = typeof(bg_closure)
        ENT = typeof(en_thermo)
        return new{n_updrafts, PM, ENT, EBGC, EC}(
            surface_area,
            max_area,
            minimum_area,
            precip_model,
            en_thermo,
            prandtl_number,
            bg_closure,
            entr_closure,
        )
    end
end
parameter_set(obj) = obj.param_set
n_updrafts(::EDMF_PrognosticTKE{N_up}) where {N_up} = N_up

struct State{P, A, T}
    prog::P
    aux::A
    tendencies::T
end
