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
    "Horizontal eddy-diffusivity"
    K_ε::FT
end

"""
    MoistureDeficitEntr

Entrainment detrainment model from Cohen et al (2020)

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MoistureDeficitEntr{FT}
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
    "pressure plume spacing"
    R_up::FT
    "updraft relative humidity"
    RH_up::FT
    "environment relative humidity"
    RH_en::FT
    "maximum updraft area"
    max_area::FT
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

abstract type EnvBuoyGradClosure end
struct BuoyGradMean <: EnvBuoyGradClosure end
struct BuoyGradQuadratures <: EnvBuoyGradClosure end

"""
    EnvBuoyGrad

Variables used in the environmental buoyancy gradient computation.

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct EnvBuoyGrad{FT, EBC <: EnvBuoyGradClosure}
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
function EnvBuoyGrad(::EBG; t_sat::FT, bg_kwargs...) where {FT <: Real, EBG <: EnvBuoyGradClosure}
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
Base.@kwdef struct MinDisspLen{FT, T}
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
    "environment area"
    a_en::FT
    "environment velocity at cell center"
    wc_en::FT
    "updraft velocity at cell center"
    wc_up::T
    "updraft area"
    a_up::T
    "up draft turbulent entrainment"
    ε_turb::T
    "updraft dynamic detrainment"
    δ_dyn::T
    "number of updraft"
    N_up::Int
end

Base.@kwdef mutable struct PrecipVariables
    precipitation_model::String = "default_precipitation_model"
    mean_rwp::Float64 = 0
    mean_swp::Float64 = 0
    cutoff_precipitation_rate::Float64 = 0
end
function PrecipVariables(namelist, grid::Grid)

    precipitation_model = parse_namelist(
        namelist,
        "microphysics",
        "precipitation_model";
        default = "None",
        valid_options = ["None", "cutoff", "clima_1m"],
    )

    if !(precipitation_model in ["None", "cutoff", "clima_1m"])
        error("precipitation model not recognized")
    end

    return PrecipVariables(; precipitation_model)
end

mutable struct UpdraftVariables{A1}
    n_updrafts::Int
    cloud_base::A1
    cloud_top::A1
    cloud_cover::A1
    updraft_top::A1
    lwp::Float64
    iwp::Float64
    function UpdraftVariables(nu, namelist, grid::Grid)
        n_updrafts = nu

        # cloud and precipitation diagnostics for output
        cloud_base = zeros(nu)
        cloud_top = zeros(nu)
        cloud_cover = zeros(nu)
        updraft_top = zeros(nu)

        lwp = 0.0
        iwp = 0.0

        A1 = typeof(cloud_base)
        return new{A1}(n_updrafts, cloud_base, cloud_top, cloud_cover, updraft_top, lwp, iwp)
    end
end

Base.@kwdef mutable struct GridMeanVariables{PS}
    param_set::PS
    lwp::Float64
    iwp::Float64
    cloud_base::Float64
    cloud_top::Float64
    cloud_cover::Float64
    EnvThermo_scheme::String
end
function GridMeanVariables(namelist, grid::Grid, param_set::PS) where {PS}
    lwp = 0.0
    iwp = 0.0

    cloud_base = 0.0
    cloud_top = 0.0
    cloud_cover = 0.0

    EnvThermo_scheme = parse_namelist(namelist, "thermodynamics", "sgs"; default = "mean")

    return GridMeanVariables(; param_set, lwp, iwp, cloud_base, cloud_top, cloud_cover, EnvThermo_scheme)
end

Base.@kwdef mutable struct EnvironmentVariables
    cloud_base::Float64 = 0
    cloud_top::Float64 = 0
    cloud_cover::Float64 = 0
    lwp::Float64 = 0
    iwp::Float64 = 0
    EnvThermo_scheme::String = "default_EnvThermo_scheme"
end
function EnvironmentVariables(namelist, grid::Grid)
    # TODO: EnvThermo_scheme is repeated in GridMeanVariables
    EnvThermo_scheme = parse_namelist(namelist, "thermodynamics", "sgs"; default = "mean")
    return EnvironmentVariables(; EnvThermo_scheme)
end

struct EnvironmentThermodynamics
    quadrature_order::Int
    quadrature_type::String
    function EnvironmentThermodynamics(namelist, grid::Grid)
        quadrature_order = parse_namelist(namelist, "thermodynamics", "quadrature_order"; default = 3)
        quadrature_type = parse_namelist(namelist, "thermodynamics", "quadrature_type"; default = "gaussian")
        return new(quadrature_order, quadrature_type)
    end
end

# Stochastic entrainment/detrainment closures:
struct NoneClosureType end
struct LogNormalClosureType end
struct SDEClosureType end

# Stochastic differential equation memory
Base.@kwdef mutable struct sde_struct{T}
    u0::Float64
    dt::Float64
end

# SurfaceMoninObukhovDry:
#     Needed for dry cases (qt=0). They have to be initialized with nonzero qtg for the
#     reference profiles. This surface subroutine sets the latent heat flux to zero
#     to prevent errors due to nonzero qtg in vars such as the obukhov_length.
# SurfaceSullivanPatton
#     Not fully implemented yet. Maybe not needed - Ignacio
struct SurfaceNone end
struct SurfaceFixedFlux end
struct SurfaceFixedCoeffs end
struct SurfaceMoninObukhov end
struct SurfaceMoninObukhovDry end
struct SurfaceSullivanPatton end

Base.@kwdef mutable struct SurfaceBase{T}
    zrough::Float64 = 0
    interactive_zrough::Bool = false
    Tsurface::Float64 = 0
    qsurface::Float64 = 0
    shf::Float64 = 0
    lhf::Float64 = 0
    cm::Float64 = 0
    ch::Float64 = 0
    cq::Float64 = 0
    bflux::Float64 = 0
    windspeed::Float64 = 0
    ustar::Float64 = 0
    rho_qtflux::Float64 = 0
    rho_hflux::Float64 = 0
    rho_uflux::Float64 = 0
    rho_vflux::Float64 = 0
    obukhov_length::Float64 = 0
    Ri_bulk_crit::Float64 = 0
    ustar_fixed::Bool = false
    ref_params::NamedTuple = NamedTuple()
end

function SurfaceBase(::Type{T}; namelist::Dict, ref_params) where {T}
    Ri_bulk_crit = namelist["turbulence"]["Ri_bulk_crit"]
    return SurfaceBase{T}(; Ri_bulk_crit, ref_params)
end

struct ForcingBaseType end
struct ForcingNone end
struct ForcingStandard end
struct ForcingDYCOMS_RF01 end
struct ForcingLES end

struct RadiationBaseType end
struct RadiationNone end
struct RadiationStandard end
struct RadiationDYCOMS_RF01 end
struct RadiationLES end

Base.@kwdef mutable struct LESData
    imin::Int = 0
    imax::Int = 0
    les_filename::String = nothing
end

"""
    ForcingBase

LES-driven forcing

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct ForcingBase{T}
    "Boolean specifying whether Coriolis forcing is applied"
    apply_coriolis::Bool = false
    "Boolean specifying whether subsidence forcing is applied"
    apply_subsidence::Bool = false
    "Coriolis parameter"
    coriolis_param::Float64 = 0
    "Momentum relaxation timescale"
    nudge_tau::Float64 = 0.0
    "Conversion function from forcing to prognostic"
    convert_forcing_prog_fp::Function = x -> x
end

force_type(::ForcingBase{T}) where {T} = T

Base.@kwdef mutable struct RadiationBase{T}
    divergence::Float64 = 0
    alpha_z::Float64 = 0
    kappa::Float64 = 0
    F0::Float64 = 0
    F1::Float64 = 0
end

rad_type(::RadiationBase{T}) where {T} = T

Base.@kwdef mutable struct CasesBase{T}
    casename::String = "default_casename"
    inversion_option::String = "default_inversion_option"
    les_filename::String = "None"
    Sur::SurfaceBase
    Fo::ForcingBase
    Rad::RadiationBase
    rad_time::StepRangeLen = range(10, 360; length = 36) .* 60
    rad::AbstractMatrix{Float64} = zeros(1, 1)
    lhf0::Float64 = 0
    shf0::Float64 = 0
    LESDat::Union{LESData, Nothing} = nothing
end

CasesBase(case::T; kwargs...) where {T} = CasesBase{T}(; casename = string(nameof(T)), kwargs...)

mutable struct EDMF_PrognosticTKE{A1}
    Ri_bulk_crit::Float64
    zi::Float64
    n_updrafts::Int
    asp_label::String
    extrapolate_buoyancy::Bool
    surface_area::Float64
    max_area::Float64
    aspect_ratio::Float64
    tke_ed_coeff::Float64
    tke_diss_coeff::Float64
    static_stab_coeff::Float64
    lambda_stab::Float64
    minimum_area::Float64
    Precip::PrecipVariables
    UpdVar::UpdraftVariables
    EnvVar::EnvironmentVariables
    EnvThermo::EnvironmentThermodynamics
    area_surface_bc::A1
    w_surface_bc::A1
    h_surface_bc::A1
    qt_surface_bc::A1
    pressure_plume_spacing::A1
    prandtl_number::Float64
    wstar::Float64
    entr_surface_bc::Float64
    detr_surface_bc::Float64
    sde_model::sde_struct
    bg_closure::EnvBuoyGradClosure
    function EDMF_PrognosticTKE(namelist, grid::Grid, param_set::PS) where {PS}
        # get values from namelist
        prandtl_number = namelist["turbulence"]["prandtl_number_0"]
        Ri_bulk_crit = namelist["turbulence"]["Ri_bulk_crit"]
        zi = 0.0

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

        asp_label = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "pressure_closure_asp_label";
            default = "const",
        )

        extrapolate_buoyancy =
            parse_namelist(namelist, "turbulence", "EDMF_PrognosticTKE", "extrapolate_buoyancy"; default = true)

        # Get values from namelist
        # set defaults at some point?
        surface_area = namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"]
        max_area = namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"]
        # entrainment parameters
        # pressure parameters
        aspect_ratio = namelist["turbulence"]["EDMF_PrognosticTKE"]["aspect_ratio"]

        # mixing length parameters
        tke_ed_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"]
        tke_diss_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"]
        static_stab_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"]
        lambda_stab = namelist["turbulence"]["EDMF_PrognosticTKE"]["lambda_stab"] # Latent heat stability effect
        # Need to code up as namelist option?
        minimum_area = 1e-5

        # Create the class for precipitation
        Precip = PrecipVariables(namelist, grid)

        # Create the updraft variable class (major diagnostic and prognostic variables)
        UpdVar = UpdraftVariables(n_updrafts, namelist, grid)

        # Create the environment variable class (major diagnostic and prognostic variables)
        EnvVar = EnvironmentVariables(namelist, grid)
        # Create the class for environment thermodynamics
        EnvThermo = EnvironmentThermodynamics(namelist, grid)

        # Near-surface BC of updraft area fraction
        area_surface_bc = zeros(n_updrafts)
        w_surface_bc = zeros(n_updrafts)
        h_surface_bc = zeros(n_updrafts)
        qt_surface_bc = zeros(n_updrafts)
        pressure_plume_spacing = zeros(n_updrafts)

        # Initialize SDE parameters
        dt = parse_namelist(namelist, "time_stepping", "dt"; default = 1.0)
        closure = parse_namelist(
            namelist,
            "turbulence",
            "EDMF_PrognosticTKE",
            "stochastic",
            "closure";
            default = "none",
            valid_options = ["none", "lognormal", "sde"],
        )
        closure_type = if closure == "none"
            NoneClosureType
        elseif closure == "lognormal"
            LogNormalClosureType
        elseif closure == "sde"
            SDEClosureType
        else
            error("Something went wrong. Invalid stochastic closure type '$closure'")
        end
        sde_model = sde_struct{closure_type}(u0 = 1, dt = dt)

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

        wstar = 0
        entr_surface_bc = 0
        detr_surface_bc = 0
        A1 = typeof(area_surface_bc)
        return new{A1}(
            Ri_bulk_crit,
            zi,
            n_updrafts,
            asp_label,
            extrapolate_buoyancy,
            surface_area,
            max_area,
            aspect_ratio,
            tke_ed_coeff,
            tke_diss_coeff,
            static_stab_coeff,
            lambda_stab,
            minimum_area,
            Precip,
            UpdVar,
            EnvVar,
            EnvThermo,
            area_surface_bc,
            w_surface_bc,
            h_surface_bc,
            qt_surface_bc,
            pressure_plume_spacing,
            prandtl_number,
            wstar,
            entr_surface_bc,
            detr_surface_bc,
            sde_model,
            bg_closure,
        )
    end
end
diffusivity_m(edmf::EDMF_PrognosticTKE) = edmf.KM
diffusivity_h(edmf::EDMF_PrognosticTKE) = edmf.KH
Ri_bulk_crit(edmf::EDMF_PrognosticTKE) = edmf.Ri_bulk_crit
parameter_set(obj) = obj.param_set
