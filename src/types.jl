"""
    PrecipFormation

Storage for tendencies due to precipitation formation

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PrecipFormation{FT}
    e_tot_tendency::FT
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
    "Turbulent entrainment"
    ε_turb::FT
    "Nondimensional fractional dynamical entrainment"
    ε_nondim::FT
    "Nondimensional fractional dynamical detrainment"
    δ_nondim::FT
end

abstract type AbstractEntrDetrModel end
abstract type AbstractNonLocalEntrDetrModel end
struct MDEntr <: AbstractEntrDetrModel end  # existing model

Base.@kwdef struct NNEntr <: AbstractEntrDetrModel
    biases_bool::Bool
end
Base.@kwdef struct NNEntrNonlocal <: AbstractNonLocalEntrDetrModel
    biases_bool::Bool
end

struct LinearEntr <: AbstractEntrDetrModel end
struct FNOEntr <: AbstractNonLocalEntrDetrModel end
struct RFEntr <: AbstractEntrDetrModel end

Base.@kwdef struct NoisyRelaxationProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end
Base.@kwdef struct LogNormalScalingProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end

Base.@kwdef struct PrognosticNoisyRelaxationProcess{MT} <: AbstractEntrDetrModel
    mean_model::MT
end

abstract type EntrDimScale end
struct BuoyVelEntrDimScale <: EntrDimScale end
struct InvZEntrDimScale <: EntrDimScale end
struct InvMeterEntrDimScale <: EntrDimScale end

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

abstract type AbstractMoistureModel end
struct EquilibriumMoisture <: AbstractMoistureModel end
struct NonEquilibriumMoisture <: AbstractMoistureModel end

abstract type AbstractPrecipitationModel end
struct NoPrecipitation <: AbstractPrecipitationModel end
struct CutoffPrecipitation <: AbstractPrecipitationModel end
struct Clima1M <: AbstractPrecipitationModel end

abstract type AbstractQuadratureType end
struct LogNormalQuad <: AbstractQuadratureType end
struct GaussianQuad <: AbstractQuadratureType end

abstract type AbstractEnvThermo end
struct SGSMean <: AbstractEnvThermo end
struct SGSQuadrature{N, QT, A, W} <: AbstractEnvThermo
    quadrature_type::QT
    a::A
    w::W
    function SGSQuadrature(namelist)
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
    ρq_liq_flux::FT = 0
    ρq_ice_flux::FT = 0
    ρe_tot_flux::FT = 0
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

struct EDMFModel{N_up, FT, MM, PM, ENT, EBGC, EC, EDS, EPG}
    surface_area::FT
    max_area::FT
    minimum_area::FT
    moisture_model::MM
    precip_model::PM
    en_thermo::ENT
    prandtl_number::FT
    bg_closure::EBGC
    entr_closure::EC
    entr_dim_scale::EDS
    entr_pi_subset::EPG
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

        moisture_model_name = parse_namelist(namelist, "thermodynamics", "moisture_model"; default = "equilibrium")

        moisture_model = if moisture_model_name == "equilibrium"
            EquilibriumMoisture()
        elseif moisture_model_name == "nonequilibrium"
            NonEquilibriumMoisture()
        else
            error("Something went wrong. Invalid moisture model: '$moisture_model_name'")
        end

        precip_model = precip_model
        # Create the environment variable class (major diagnostic and prognostic variables)

        # Create the class for environment thermodynamics and buoyancy gradient computation
        en_sgs_name =
            parse_namelist(namelist, "thermodynamics", "sgs"; default = "mean", valid_options = ["mean", "quadrature"])
        en_thermo = if en_sgs_name == "mean"
            SGSMean()
        elseif en_sgs_name == "quadrature"
            SGSQuadrature(namelist)
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
        if moisture_model_name == "nonequilibrium" && en_thermo_name == "quadrature"
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
            valid_options = ["moisture_deficit", "NN", "NN_nonlocal", "FNO", "Linear", "RF"],
        )

        nn_biases = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "nn_ent_biases";
            default = false,
            valid_options = [true, false],
        )

        mean_entr_closure = if entr_type == "moisture_deficit"
            MDEntr()
        elseif entr_type == "NN"
            NNEntr(biases_bool = nn_biases)
        elseif entr_type == "NN_nonlocal"
            NNEntrNonlocal(biases_bool = nn_biases)
        elseif entr_type == "FNO"
            FNOEntr()
        elseif entr_type == "Linear"
            LinearEntr()
        elseif entr_type == "RF"
            RFEntr()
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
        elseif stoch_entr_type == "prognostic_noisy_relaxation_process"
            PrognosticNoisyRelaxationProcess(mean_model = mean_entr_closure)
        else
            error("Something went wrong. Invalid stochastic entrainment type '$stoch_entr_type'")
        end

        entr_dim_scale = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "entr_dim_scale";
            default = "buoy_vel",
            valid_options = ["buoy_vel", "inv_z", "none"],
        )

        entr_dim_scale = if entr_dim_scale == "buoy_vel"
            BuoyVelEntrDimScale()
        elseif entr_dim_scale == "inv_z"
            InvZEntrDimScale()
        elseif entr_dim_scale == "none"
            InvMeterEntrDimScale()
        else
            error("Something went wrong. Invalid entrainment dimension scale '$entr_dim_scale'")
        end

        entr_pi_subset = parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "entr_pi_subset")

        EDS = typeof(entr_dim_scale)
        EC = typeof(entr_closure)
        MM = typeof(moisture_model)
        PM = typeof(precip_model)
        EBGC = typeof(bg_closure)
        ENT = typeof(en_thermo)
        EPG = typeof(entr_pi_subset)
        return new{n_updrafts, FT, MM, PM, ENT, EBGC, EC, EDS, EPG}(
            surface_area,
            max_area,
            minimum_area,
            moisture_model,
            precip_model,
            en_thermo,
            prandtl_number,
            bg_closure,
            entr_closure,
            entr_dim_scale,
            entr_pi_subset,
        )
    end
end
parameter_set(obj) = obj.param_set
n_updrafts(::EDMFModel{N_up}) where {N_up} = N_up
Base.eltype(::EDMFModel{N_up, FT}) where {N_up, FT} = FT
n_Π_groups(m::EDMFModel) = length(m.entr_pi_subset)
entrainment_Π_subset(m::EDMFModel) = m.entr_pi_subset
Base.broadcastable(edmf::EDMFModel) = Ref(edmf)

struct State{P, A, T}
    prog::P
    aux::A
    tendencies::T
end

"""
    column_state(prog, aux, tendencies, inds...)

Create a columnar state given full 3D states
 - `prog` prognostic state
 - `aux` auxiliary state
 - `tendencies` tendencies state
 - `inds` `i`, `j`, `h` indices

## Example
```julia
local_geom = ClimaCore.Spaces.local_geometry_data(space)
Ni, Nj, _, _, Nh = size(local_geom)
for h in 1:Nh, j in 1:Nj, i in 1:Ni
    inds = (i, j, h)
    state = TC.column_state(prog, aux, tendencies, inds...)
    ...
end
"""
function column_state(prog, aux, tendencies, inds...)
    prog_cent_column = CC.column(prog.cent, inds...)
    prog_face_column = CC.column(prog.face, inds...)
    aux_cent_column = CC.column(aux.cent, inds...)
    aux_face_column = CC.column(aux.face, inds...)
    tends_cent_column = CC.column(tendencies.cent, inds...)
    tends_face_column = CC.column(tendencies.face, inds...)
    prog_column = CC.Fields.FieldVector(cent = prog_cent_column, face = prog_face_column)
    aux_column = CC.Fields.FieldVector(cent = aux_cent_column, face = aux_face_column)
    tends_column = CC.Fields.FieldVector(cent = tends_cent_column, face = tends_face_column)

    return State(prog_column, aux_column, tends_column)
end
