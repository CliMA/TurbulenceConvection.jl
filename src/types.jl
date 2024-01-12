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
struct NonEquilibriumMoisture <: AbstractMoistureModel end

abstract type AbstractCovarianceModel end
struct PrognosticThermoCovariances <: AbstractCovarianceModel end
struct DiagnosticThermoCovariances{FT} <: AbstractCovarianceModel
    covar_lim::FT
end

abstract type AbstractPrecipitationModel end
struct NoPrecipitation <: AbstractPrecipitationModel end
struct Clima0M <: AbstractPrecipitationModel end
struct Clima1M <: AbstractPrecipitationModel end

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

abstract type AbstractPrecipFractionModel end
struct PrescribedPrecipFraction{FT} <: AbstractPrecipFractionModel
    prescribed_precip_frac_value::FT
end
struct DiagnosticPrecipFraction{FT} <: AbstractPrecipFractionModel
    precip_fraction_limiter::FT
end

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
) where {FT, FVT}
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
) where {FT, FVT}
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

struct EDMFModel{N_up, FT, SABC, MM, TCM, PM, RFM, PFM, ENT, EBGC, MLP, PMP, EC, MLEC, ET, EDS, DDS, EPG}
    surface_area::FT
    surface_area_bc::SABC
    max_area::FT
    minimum_area::FT
    moisture_model::MM
    thermo_covariance_model::TCM
    precip_model::PM
    rain_formation_model::RFM
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
end
function EDMFModel(::Type{FT}, namelist, precip_model, rain_formation_model) where {FT}

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
        NonEquilibriumMoisture()
    else
        error("Something went wrong. Invalid moisture model: '$moisture_model_name'")
    end

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
        PrescribedPrecipFraction(prescribed_precip_frac_value)
    elseif precip_fraction_model_name == "cloud_cover"
        precip_fraction_limiter = parse_namelist(namelist, "microphysics", "precip_fraction_limiter"; default = 0.3)
        DiagnosticPrecipFraction(precip_fraction_limiter)
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
    TCM = typeof(thermo_covariance_model)
    PM = typeof(precip_model)
    RFM = typeof(rain_formation_model)
    PFM = typeof(precip_fraction_model)
    EBGC = typeof(bg_closure)
    ENT = typeof(en_thermo)
    EPG = typeof(entr_pi_subset)
    MLP = typeof(mixing_length_params)
    PMP = typeof(pressure_model_params)
    return EDMFModel{n_updrafts, FT, SABC, MM, TCM, PM, RFM, PFM, ENT, EBGC, MLP, PMP, EC, MLEC, ET, EDS, DDS, EPG}(
        surface_area,
        surface_area_bc,
        max_area,
        minimum_area,
        moisture_model,
        thermo_covariance_model,
        precip_model,
        rain_formation_model,
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
