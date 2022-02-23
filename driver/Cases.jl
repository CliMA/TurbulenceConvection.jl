module Cases

import NCDatasets
const NC = NCDatasets

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

import ClimaCore
const CC = ClimaCore
const CCO = CC.Operators

import AtmosphericProfilesLibrary
const APL = AtmosphericProfilesLibrary

import Dierckx
import Statistics
import Random
import UnPack

import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet
const APS = CP.AbstractEarthParameterSet

import Thermodynamics
const TD = Thermodynamics

import ..TurbulenceConvection
const TC = TurbulenceConvection

import ..TurbulenceConvection: ForcingBase
import ..TurbulenceConvection: RadiationBase

using ..TurbulenceConvection: CasesBase
using ..TurbulenceConvection: pyinterp
using ..TurbulenceConvection: Grid
using ..TurbulenceConvection: real_center_indices
using ..TurbulenceConvection: real_face_indices
using ..TurbulenceConvection: get_inversion

# Support interpolating with z-points
(prof::Dierckx.Spline1D)(z::CC.Geometry.ZPoint{FT}) where {FT} = prof(z.z)
(prof::Dierckx.Spline2D)(t, z::CC.Geometry.ZPoint{FT}) where {FT} = prof(t, z.z)

#=
    arr_type(x)

We're keeping this around in case we need
to move some initialized data to the GPU.
In this case, we may be able to modify this
function to do this.
=#
arr_type(x) = x

#####
##### Case types
#####

# For dispatching to inherited class
struct BaseCase end

abstract type AbstractCaseType end

""" [Soares2004](@cite) """
struct Soares <: AbstractCaseType end

""" [Nieuwstadt1993](@cite) """
struct Nieuwstadt <: AbstractCaseType end

struct Bomex <: AbstractCaseType end

""" [Tan2018](@cite) """
struct life_cycle_Tan2018 <: AbstractCaseType end

struct Rico <: AbstractCaseType end

""" [Grabowski2006](@cite) """
struct TRMM_LBA <: AbstractCaseType end

""" [Brown2002](@cite) """
struct ARM_SGP <: AbstractCaseType end

""" [Khairoutdinov2009](@cite) """
struct GATE_III <: AbstractCaseType end

""" [Stevens2005](@cite) """
struct DYCOMS_RF01 <: AbstractCaseType end

""" [Ackerman2009](@cite) """
struct DYCOMS_RF02 <: AbstractCaseType end

struct GABLS <: AbstractCaseType end

struct SP <: AbstractCaseType end

struct DryBubble <: AbstractCaseType end

struct LES_driven_SCM <: AbstractCaseType end

# Abstract parameter struct for case-based parameters
abstract type AbstractCaseParameters end

#####
##### Case methods
#####

include("Radiation.jl")
include("Forcing.jl")

get_case(namelist::Dict) = get_case(namelist["meta"]["casename"])
get_case(casename::String) = get_case(Val(Symbol(casename)))
get_case(::Val{:Soares}) = Soares()
get_case(::Val{:Nieuwstadt}) = Nieuwstadt()
get_case(::Val{:Bomex}) = Bomex()
get_case(::Val{:life_cycle_Tan2018}) = life_cycle_Tan2018()
get_case(::Val{:Rico}) = Rico()
get_case(::Val{:TRMM_LBA}) = TRMM_LBA()
get_case(::Val{:ARM_SGP}) = ARM_SGP()
get_case(::Val{:GATE_III}) = GATE_III()
get_case(::Val{:DYCOMS_RF01}) = DYCOMS_RF01()
get_case(::Val{:DYCOMS_RF02}) = DYCOMS_RF02()
get_case(::Val{:GABLS}) = GABLS()
get_case(::Val{:SP}) = SP()
get_case(::Val{:DryBubble}) = DryBubble()
get_case(::Val{:LES_driven_SCM}) = LES_driven_SCM()

get_case_name(case_type::AbstractCaseType) = string(case_type)

#####
##### Case configurations
#####

inversion_type(::AbstractCaseType) = TC.CriticalRiInversion()
inversion_type(::TRMM_LBA) = TC.max∇θInversion()
inversion_type(::ARM_SGP) = TC.max∇θInversion()
inversion_type(::GATE_III) = TC.max∇θInversion()
inversion_type(::DYCOMS_RF01) = TC.max∇θInversion()
inversion_type(::DYCOMS_RF02) = TC.max∇θInversion()
inversion_type(::DryBubble) = TC.θρInversion()

get_forcing_type(::AbstractCaseType) = TC.ForcingStandard # default
get_forcing_type(::Soares) = TC.ForcingNone
get_forcing_type(::Nieuwstadt) = TC.ForcingNone
get_forcing_type(::DYCOMS_RF01) = TC.ForcingDYCOMS_RF01
get_forcing_type(::DYCOMS_RF02) = TC.ForcingDYCOMS_RF01
get_forcing_type(::DryBubble) = TC.ForcingNone
get_forcing_type(::LES_driven_SCM) = TC.ForcingLES

get_radiation_type(::AbstractCaseType) = TC.RadiationNone # default
get_radiation_type(::DYCOMS_RF01) = TC.RadiationDYCOMS_RF01
get_radiation_type(::DYCOMS_RF02) = TC.RadiationDYCOMS_RF01
get_radiation_type(::LES_driven_SCM) = TC.RadiationLES

RadiationBase(case::AbstractCaseType, param_set::ACPS) where {ACPS <: AbstractCaseParameters} = RadiationBase{Cases.get_radiation_type(case)}()

forcing_kwargs(::AbstractCaseType, namelist) = ()
surface_param_kwargs(::AbstractCaseType, namelist) = ()

ForcingBase(case::AbstractCaseType, param_set::APS; kwargs...) = ForcingBase(get_forcing_type(case); kwargs...)

#####
##### Default CasesBase behavior:
#####

initialize_radiation(self::CasesBase, grid, state, param_set) = nothing

update_forcing(self::CasesBase, grid, state, t::Real, param_set) = nothing
initialize_forcing(self::CasesBase, grid::Grid, state, param_set) = initialize(self.Fo, grid, state)


#####
##### Soares
#####

struct SoaresParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    zrough::FT
    Tsurface::FT
    qsurface::FT
    θ_flux::FT
    qt_flux::FT
    ustar::FT    
    TPS::ThermodynamicsParameters{FT}
end

function SoaresParameters(
    param_set
)

    aliases =["Pg","qtg","Tg","zrough", "Tsurface","qsurface", "θ_flux", "qt_flux", "ustar"]
    (Pg, qtg, Tg, zrough, Tsurface, qsurface, θ_flux, qt_flux, ustar) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "Soares")
    
    return SoaresParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        zrough,
        Tsurface,
        qsurface,
        θ_flux,
        qt_flux,
        ustar,
        TPS,
    )
    
end


ForcingBase(case::Soares, param_set::SoaresParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::Soares, param_set::SoaresParameters, namelist)
    Pg = param_set.Pg
    qtg = param_set.qtg
    Tg = param_set.Tg
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{Soares}, grid::Grid, param_set::SoaresParameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    FT = eltype(grid)
    prof_q_tot = APL.Soares_q_tot(FT)
    prof_θ_liq_ice = APL.Soares_θ_liq_ice(FT)
    prof_u = APL.Soares_u(FT)
    prof_tke = APL.Soares_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::Soares, grid::TC.Grid, surf_ref_state, param_set::SoaresParameters; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set.TPS,surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set.TPS,surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = param_set.zrough #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    Tsurface = param_set.Tsurface
    qsurface = param_set.qsurface
    θ_flux = param_set.θ_flux
    qt_flux = param_set.qt_flux 
    ts = TD.PhaseEquil_pTq(param_set.TPS, p0_f_surf, Tsurface, qsurface)
    lhf = qt_flux * ρ0_f_surf * TD.latent_heat_vapor(param_set.TPS, ts)
    shf = θ_flux * TD.cp_m(param_set.TPS, ts) * ρ0_f_surf
    ustar = param_set.ustar # just to initilize grid mean covariances
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

#####
##### Nieuwstadt
#####

struct NieuwstadtParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    zrough::FT
    Tsurface::FT
    qsurface::FT
    θ_flux::FT
    lhf::FT
    ustar::FT
    TPS::ThermodynamicsParameters{FT}
end

function NieuwstadtParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases =["Pg","qtg","Tg","zrough", "Tsurface","qsurface", "θ_flux",  "lhf", "ustar"]
    (Pg, qtg, Tg, zrough, Tsurface, qsurface, θ_flux, lhf, ustar) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "Nieuwstadt")
    
    return NieuwstadtParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        zrough,
        Tsurface,
        qsurface,
        θ_flux,
        lhf,
        ustar,
        TPS,
    )
    
end


ForcingBase(case::Nieuwstadt, param_set::NieuwstadtParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::Nieuwstadt, param_set::NieuwstadtParameters, namelist)
    Pg = param_set.Pg
    Tg = param_set.Tg
    qtg = param_set.qtg # Total water mixing ratio
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{Nieuwstadt}, grid::Grid, param_set::NieuwstadtParameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)

    FT = eltype(grid)
    prof_θ_liq_ice = APL.Nieuwstadt_θ_liq_ice(FT)
    prof_u = APL.Nieuwstadt_u(FT)
    prof_tke = APL.Nieuwstadt_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::Nieuwstadt, grid::TC.Grid, surf_ref_state, param_set::NieuwstadtParameters; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set.TPS, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set.TPS, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = param_set.zrough #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    Tsurface = param_set.Tsurface
    qsurface = param_set.qsurface
    θ_flux = param_set.θ_flux
    lhf = param_set.lhf # It would be 0.0 if we follow Nieuwstadt.
    ts = TD.PhaseEquil_pTq(param_set.TPS, p0_f_surf, Tsurface, qsurface)
    shf = θ_flux * TD.cp_m(param_set.TPS, ts) * ρ0_f_surf
    ustar = param_set.ustar # just to initilize grid mean covariances
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

#####
##### Bomex
#####


struct BomexParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    zrough::FT
    qsurface::FT
    θ_surface::FT
    θ_flux::FT
    qtflux::FT
    ustar::FT
    coriolis_param::FT
    TPS::ThermodynamicsParameters{FT}
end

function BomexParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases =["Pg","qtg","Tg","zrough", "qsurface","θ_surface", "θ_flux",  "qt_flux", "ustar","coriolis_param"]
    (Pg, qtg, Tg, zrough, qsurface, θ_surface, θ_flux, qt_flux, ustar, coriolis_param) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "Bomex")
    
    return BomexParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        zrough,
        qsurface,
        θ_surface,
        θ_flux,
        qt_flux,
        ustar,
        coriolis_param,
        TPS,
    )
    
end

# replace coriolis_param
ForcingBase(case::Bomex, param_set::BomexParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = true, apply_subsidence = true, coriolis_param = param_set.coriolis_param) #= s^{-1} =#

function surface_ref_state(::Bomex, param_set::BomexParameters, namelist)
    Pg = param_set.Pg #Pressure at ground
    Tg = param_set.Tg #Temperature at ground
    qtg = param_set.qtg#Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{Bomex}, grid::Grid, param_set::BomexParameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)

    FT = eltype(grid)
    prof_q_tot = APL.Bomex_q_tot(FT)
    prof_θ_liq_ice = APL.Bomex_θ_liq_ice(FT)
    prof_u = APL.Bomex_u(FT)
    prof_tke = APL.Bomex_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::Bomex, grid::TC.Grid, surf_ref_state, param_set::BomexParameters; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set.TPS, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set.TPS, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = param_set.zrough
    qsurface = param_set.qsurface # kg/kg
    θ_surface = param_set.θ_surface 
    θ_flux = param_set.θ_flux
    qt_flux = param_set.qt_flux
    ts = TD.PhaseEquil_pθq(param_set.TPS, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set.TPS, ts)
    lhf = qt_flux * ρ0_f_surf * TD.latent_heat_vapor(param_set.TPS, ts)
    shf = θ_flux * TD.cp_m(param_set.TPS, ts) * ρ0_f_surf
    ustar = param_set.ustar # m/s
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{Bomex}, grid::Grid, state, param_set::BomexParameters)
    initialize(self.Fo, grid, state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)

    FT = eltype(grid)
    prof_ug = APL.Bomex_geostrophic_u(FT)
    prof_dTdt = APL.Bomex_dTdt(FT)
    prof_dqtdt = APL.Bomex_dqtdt(FT)
    prof_subsidence = APL.Bomex_subsidence(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        # Geostrophic velocity profiles. vg = 0
        aux_gm.ug[k] = prof_ug(z)
        ts = TD.PhaseEquil_pθq(param_set.TPS,  p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(param_set.TPS, ts)
        # Set large-scale cooling
        aux_gm.dTdt[k] = prof_dTdt(Π, z)
        # Set large-scale drying
        aux_gm.dqtdt[k] = prof_dqtdt(z)
        #Set large scale subsidence
        aux_gm.subsidence[k] = prof_subsidence(z)
    end
    return nothing
end

#####
##### life_cycle_Tan2018
#####


struct life_cycle_Tan2018Parameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    zrough::FT
    qsurface::FT
    θ_surface::FT
    θ_flux::FT
    qtflux::FT
    ustar::FT
    coriolis_param::FT
    grav::FT
    molmass_dryair::FT
    molmass_water::FT
    molmass_ratio::FT
    TPS::ThermodynamicsParameters{FT}
end

function life_cycle_Tan2018Parameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases =["Pg","qtg","Tg","zrough", "qsurface","θ_surface", "θ_flux",  "qt_flux", "ustar", "coriolis_param",
              "grav", "molmass_dryair", "molmass_water"]
    (Pg, qtg, Tg, zrough, qsurface, θ_surface, θ_flux, qt_flux, ustar, coriolis_param, grav, molmass_dryair, molmass_water) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "life_cycle_Tan2018")
    #derived parameters
    molmass_ratio = molmass_dryair/molmass_water
      
    return life_cycle_Tan2018Parameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        zrough,
        qsurface,
        θ_surface,
        θ_flux,
        qt_flux,
        ustar,
        coriolis_param,
        grav,
        molmass_dryair,
        molmass_water,
        molmass_ratio,
        TPS,
    )
    
end



ForcingBase(case::life_cycle_Tan2018, param_set::life_cycle_Tan2018Parameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = true, apply_subsidence = true, coriolis_param = param_set.coriolis_param) #= s^{-1} =#

function surface_ref_state(::life_cycle_Tan2018, param_set::life_cycle_Tan2018Parameters, namelist)
    Pg = param_set.Pg #Pressure at ground
    Tg = param_set.Tg  #Temperature at ground
    qtg = param_set.qtg   #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{life_cycle_Tan2018}, grid::Grid, param_set::life_cycle_Tan2018Parameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)

    FT = eltype(grid)
    prof_q_tot = APL.LifeCycleTan2018_q_tot(FT)
    prof_θ_liq_ice = APL.LifeCycleTan2018_θ_liq_ice(FT)
    prof_u = APL.LifeCycleTan2018_u(FT)
    prof_tke = APL.LifeCycleTan2018_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function life_cycle_buoyancy_flux(param_set::life_cycle_Tan2018Parameters, weight = 1)
    g = param_set.grav
    molmass_ratio = param_set.molmass_ratio
    return g * (
        (8.0e-3 * weight + (molmass_ratio - 1) * (299.1 * 5.2e-5 * weight + 22.45e-3 * 8.0e-3 * weight)) /
        (299.1 * (1 + (molmass_ratio - 1) * 22.45e-3))
    )
end

function surface_params(case::life_cycle_Tan2018, grid::TC.Grid, surf_ref_state, param_set::life_cycle_Tan2018Parameters; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set.TPS, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set.TPS, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = param_set.zrough # not actually used, but initialized to reasonable value
    qsurface = param_set.qsurface # kg/kg
    θ_surface = param_set.θ_surface
    θ_flux = param_set.θ_flux
    qt_flux = param_set.qt_flux
    ts = TD.PhaseEquil_pθq(param_set.TPS, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set.TPS, ts)
    lhf0 = qt_flux * ρ0_f_surf * TD.latent_heat_vapor(param_set.TPS, ts)
    shf0 = θ_flux * TD.cp_m(param_set.TPS, ts) * ρ0_f_surf

    weight_factor(t) = 0.01 + 0.99 * (cos(2.0 * π * t / 3600.0) + 1.0) / 2.0
    weight = 1.0
    lhf = t -> lhf0 * (weight * weight_factor(t))
    shf = t -> shf0 * (weight * weight_factor(t))

    ustar = param_set.ustar # m/s
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{life_cycle_Tan2018}, grid::Grid, state, param_set::life_cycle_Tan2018Parameters)
    initialize(self.Fo, grid, state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)

    FT = eltype(grid)
    prof_ug = APL.LifeCycleTan2018_geostrophic_u(FT)
    prof_dTdt = APL.LifeCycleTan2018_dTdt(FT)
    prof_dqtdt = APL.LifeCycleTan2018_dqtdt(FT)
    prof_subsidence = APL.LifeCycleTan2018_subsidence(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        # Geostrophic velocity profiles. vg = 0
        aux_gm.ug[k] = prof_ug(z)
        ts = TD.PhaseEquil_pθq(param_set.TPS, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(param_set.TPS, ts)
        # Set large-scale cooling
        aux_gm.dTdt[k] = prof_dTdt(Π, z)
        # Set large-scale drying
        aux_gm.dqtdt[k] = prof_dqtdt(z)
        #Set large scale subsidence
        aux_gm.subsidence[k] = prof_subsidence(z)
    end
    return nothing
end

#####
##### Rico
#####

struct RicoParameters{FT} <: AbstractCaseParameters
    Pg::FT
    Tg::FT
    zrough::FT
    cm::FT
    ch::FT
    cq::FT
    Tsurface::FT
    molmass_dryair::FT
    molmass_water::FT
    molmass_ratio::FT
    Omega::FT
    latitude::FT
    TPS::ThermodynamicsParameters{FT}
end

function RicoParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases =["Pg","Tg","zrough", "cm", "ch", "cq", "Tsurface","molmass_dryair", "molmass_water", "Omega","latitude"]
    (Pg, Tg, zrough, qsurface, cm, ch, cq, Tsurface, molmass_dryair, molmass_water, Omega, latitude) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "Rico")
    #derived parameters
    molmass_ratio = molmass_dryair/molmass_water
    
    return RicoParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        Tg,
        zrough,
        cm,
        ch,
        cq,
        Tsurface,
        molmass_dryair,
        molmass_water,
        molmass_ratio,
        Omega,
        latitude,
        TPS,
    )
    
end


function ForcingBase(case::Rico, param_set::RicoParameters; kwargs...)
    latitude = param_set.latitude
    Omega = param_set.Omega
    return ForcingBase(
        get_forcing_type(case);
        apply_coriolis = true,
        apply_subsidence = true,
        coriolis_param = 2.0 * Omega * sin(latitude * π / 180.0),
    ) #= s^{-1} =#
end

function surface_ref_state(::Rico, param_set::RicoParameters, namelist)
    molmass_ratio = param_set.molmass_ratio
    Pg = param_set.Pg #Pressure at ground
    Tg = param_set.Tg  #Temperature at ground
    pvg = TD.saturation_vapor_pressure(param_set.TPS, Tg, TD.Liquid())
    qtg = (1 / molmass_ratio) * pvg / (Pg - pvg)   #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{Rico}, grid::Grid, param_set::RicoParameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_tc = TC.center_aux_turbconv(state)
    p0 = TC.center_ref_state(state).p0

    FT = eltype(grid)
    prof_u = APL.Rico_u(FT)
    prof_v = APL.Rico_v(FT)
    prof_q_tot = APL.Rico_q_tot(FT)
    prof_θ_liq_ice = APL.Rico_θ_liq_ice(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.u[k] = prof_u(z)
        prog_gm.v[k] = prof_v(z)
        prog_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.q_tot[k] = prof_q_tot(z)
    end

    # Need to get θ_virt
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set.TPS, p0[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        aux_gm.θ_virt[k] = TD.virtual_pottemp(param_set.TPS, ts)
    end
    zi = 0.6 * get_inversion(grid, state, param_set, 0.2)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.tke[k] = if z <= zi
            1 - z / zi
        else
            FT(0)
        end
    end
end

function surface_params(case::Rico, grid::TC.Grid, surf_ref_state, param_set::RicoParameters; kwargs...)
    p0_f_surf = TD.air_pressure(param_set.TPS, surf_ref_state)
    FT = eltype(p0_f_surf)

    zrough = param_set.zrough
    cm = param_set.cm
    ch = param_set.ch
    cq = param_set.cq
    # Adjust for non-IC grid spacing
    grid_adjust = (log(20.0 / zrough) / log(TC.zc_surface(grid) / zrough))^2
    cm = cm * grid_adjust
    ch = ch * grid_adjust
    cq = cq * grid_adjust # TODO: not yet used..
    Tsurface = param_set.Tsurface

    # For Rico we provide values of transfer coefficients
    ts = TD.PhaseEquil_pTq(param_set.TPS, p0_f_surf, Tsurface, FT(0)) # TODO: is this correct?
    qsurface = TD.q_vap_saturation(param_set.TPS, ts)
    kwargs = (; zrough, Tsurface, qsurface, cm, ch)
    return TC.FixedSurfaceCoeffs(FT; kwargs...)
end

function initialize_forcing(self::CasesBase{Rico}, grid::Grid, state, param_set::RicoParameters)
    initialize(self.Fo, grid, state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)

    FT = eltype(grid)
    prof_ug = APL.Rico_geostrophic_ug(FT)
    prof_vg = APL.Rico_geostrophic_vg(FT)
    prof_dTdt = APL.Rico_dTdt(FT)
    prof_dqtdt = APL.Rico_dqtdt(FT)
    prof_subsidence = APL.Rico_subsidence(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        ts = TD.PhaseEquil_pθq(param_set.TPS, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(param_set.TPS, ts)
        # Geostrophic velocity profiles
        aux_gm.ug[k] = prof_ug(z)
        aux_gm.vg[k] = prof_vg(z)
        aux_gm.dTdt[k] = prof_dTdt(Π, z) # Set large-scale cooling
        aux_gm.dqtdt[k] = prof_dqtdt(z) # Set large-scale moistening
        aux_gm.subsidence[k] = prof_subsidence(z) #Set large scale subsidence
    end
    return nothing
end

#####
##### TRMM_LBA
#####

struct TRMM_LBAParameters{FT} <: AbstractCaseParameters
    Pg::FT
    Tg::FT
    qsurface::FT
    θ_surface::FT
    ustar::FT
    molmass_dryair::FT
    molmass_water::FT
    molmass_ratio::FT
    TPS::ThermodynamicsParameters{FT}
end

function TRMM_LBAParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases =["Pg", "Tg", "qsurface","θ_surface", "ustar", "molmass_dryair", "molmass_water"]
    (Pg,  Tg, qsurface, θ_surface, θ_flux, qt_flux, ustar,  molmass_dryair, molmass_water) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "TRMM_LBA")

    #derived parameters
    molmass_ratio = molmass_dryair/molmass_water
      
    return TRMM_LBAParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        Tg,  
        qsurface,
        θ_surface,
        ustar,
        molmass_dryair,
        molmass_water,
        molmass_ratio,
        TPS,
    )
    
end


ForcingBase(case::TRMM_LBA, param_set::TRMM_LBAParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false, kwargs...)

function surface_ref_state(::TRMM_LBA, param_set::TRMM_LBAParameters, namelist)
    molmass_ratio = param_set.molmass_ratio
    Pg = param_set.Pg  #Pressure at ground
    Tg = param_set.Tg   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    pvg = TD.saturation_vapor_pressure(param_set.TPS, Tg, TD.Liquid())
    qtg = (1 / molmass_ratio) * pvg / (Pg - pvg) #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{TRMM_LBA}, grid::Grid, param_set::TRMM_LBAParameters, state)
    p0 = TC.center_ref_state(state).p0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)

    FT = eltype(grid)
    # Get profiles from AtmosphericProfilesLibrary.jl
    prof_p = APL.TRMM_LBA_p(FT)
    prof_T = APL.TRMM_LBA_T(FT)
    prof_RH = APL.TRMM_LBA_RH(FT)
    prof_u = APL.TRMM_LBA_u(FT)
    prof_v = APL.TRMM_LBA_v(FT)
    prof_tke = APL.TRMM_LBA_tke(FT)

    zc = grid.zc
    molmass_ratio = param_set.molmass_ratio
    prog_gm = TC.center_prog_grid_mean(state)

    prog_gm.u .= prof_u.(zc)
    prog_gm.v .= prof_v.(zc)
    aux_gm.T .= prof_T.(zc)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        pv_star = TD.saturation_vapor_pressure(param_set.TPS, aux_gm.T[k], TD.Liquid())
        # eq. 37 in pressel et al and the def of RH
        RH = prof_RH(z)
        denom = (prof_p(z) - pv_star + (1 / molmass_ratio) * pv_star * RH / 100.0)
        qv_star = pv_star * (1 / molmass_ratio) / denom
        prog_gm.q_tot[k] = qv_star * RH / 100
        phase_part = TD.PhasePartition(prog_gm.q_tot[k], 0.0, 0.0) # initial state is not saturated
        prog_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp_given_pressure(param_set.TPS, aux_gm.T[k], p0[k], phase_part)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::TRMM_LBA, grid::TC.Grid, surf_ref_state, param_set::TRMM_LBAParameters; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set.TPS,surf_ref_state)
    FT = eltype(p0_f_surf)
    qsurface = param_set.qsurface # kg/kg
    θ_surface = param_set.θ_surface
    ts = TD.PhaseEquil_pθq(param_set.TPS, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set.TPS, ts)
    ustar = param_set.ustar # this is taken from Bomex -- better option is to approximate from LES tke above the surface
    lhf = t -> 554.0 * max(0, cos(π / 2 * ((5.25 * 3600.0 - t) / 5.25 / 3600.0)))^1.3
    shf = t -> 270.0 * max(0, cos(π / 2 * ((5.25 * 3600.0 - t) / 5.25 / 3600.0)))^1.5
    kwargs = (; Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit, zero_uv_fluxes = true)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function forcing_kwargs(::TRMM_LBA, namelist)
    return (; rad = APL.TRMM_LBA_radiation(Float64))
end

function initialize_forcing(self::CasesBase{TRMM_LBA}, grid::Grid, state, param_set::TRMM_LBAParameters)
    aux_gm = TC.center_aux_grid_mean(state)
    rad = self.Fo.rad
    @inbounds for k in real_center_indices(grid)
        aux_gm.dTdt[k] = rad(0, grid.zc[k].z)
    end
    return nothing
end

function update_forcing(self::CasesBase{TRMM_LBA}, grid, state, t::Real, param_set::TRMM_LBAParameters)
    aux_gm = TC.center_aux_grid_mean(state)
    rad = self.Fo.rad
    @inbounds for k in real_center_indices(grid)
        aux_gm.dTdt[k] = rad(t, grid.zc[k].z)
    end
end

#####
##### ARM_SGP
#####

#TODO: Array functionality for SH,LH etc.
struct ARM_SGPParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    qsurface::FT
    θ_surface::FT
    ustar::FT
    coriolis_param::FT
    TPS::ThermodynamicsParameters{FT}
end

function ARM_SGPParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg", "Tg", "qsurface","θ_surface", "ustar", "coriolis_param"]
    (Pg, qtg, Tg, qsurface, θ_surface, θ_flux, qt_flux, ustar,  coriolis_param) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "ARM_SGP")

    # array_aliases = ["t_Sur_in","SH","LH"]
    #(t_Sur_in,SH,LH) = CLiMAParameters_get_parameter_value_arrays!(param_set,array_aliases,"ARM_SGP")
    
      
    return ARM_SGPParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,  
        qsurface,
        θ_surface,
        ustar,
        coriolis_param,
        TPS,
    )
    
end

ForcingBase(case::ARM_SGP, param_set::ARM_SGPParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = true, apply_subsidence = false, coriolis_param = 8.5e-5)

function surface_ref_state(::ARM_SGP, param_set::ARM_SGPParameters, namelist)
    Pg = param_set.Pg  #Pressure at ground
    Tg = param_set.Tg   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    qtg = param_set.qtg #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{ARM_SGP}, grid::Grid, param_set::ARM_SGPParameters, state)
    # ARM_SGP inputs
    p0 = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)

    FT = eltype(grid)
    prof_u = APL.ARM_SGP_u(FT)
    prof_q_tot = APL.ARM_SGP_q_tot(FT)
    prof_θ_liq_ice = APL.ARM_SGP_θ_liq_ice(FT)
    prof_tke = APL.ARM_SGP_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        # TODO figure out how to use ts here
        phase_part = TD.PhasePartition(prog_gm.q_tot[k], aux_gm.q_liq[k], 0.0)
        Π = TD.exner_given_pressure(param_set.TPS, p0[k], phase_part)
        prog_gm.u[k] = prof_u(z)
        prog_gm.q_tot[k] = prof_q_tot(z)
        aux_gm.T[k] = prof_θ_liq_ice(z) * Π
        prog_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp_given_pressure(param_set.TPS, aux_gm.T[k], p0[k], phase_part)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::ARM_SGP, grid::TC.Grid, surf_ref_state, param_set::ARM_SGPParameters; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set.TPS, surf_ref_state)
    FT = eltype(p0_f_surf)
    qsurface = param_set.qsurface # kg/kg
    θ_surface =  param_set.θ_surface
    ts = TD.PhaseEquil_pθq(param_set.TPS, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set.TPS, ts)
    ustar =  param_set.ustar # this is taken from Bomex -- better option is to approximate from LES tke above the surface

    # TODO ORAD: Array type parameters
    t_Sur_in = arr_type([0.0, 4.0, 6.5, 7.5, 10.0, 12.5, 14.5]) .* 3600 #LES time is in sec
    SH = arr_type([-30.0, 90.0, 140.0, 140.0, 100.0, -10, -10]) # W/m^2
    LH = arr_type([5.0, 250.0, 450.0, 500.0, 420.0, 180.0, 0.0]) # W/m^2
    shf = Dierckx.Spline1D(t_Sur_in, SH; k = 1)
    lhf = Dierckx.Spline1D(t_Sur_in, LH; k = 1)

    kwargs = (; Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit, zero_uv_fluxes = true)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{ARM_SGP}, grid::Grid, state, param_set::ARM_SGPParameters)
    aux_gm = TC.center_aux_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        aux_gm.ug[k] = 10.0
        aux_gm.vg[k] = 0.0
    end
    return nothing
end

function update_forcing(self::CasesBase{ARM_SGP}, grid, state, t::Real, param_set::ARM_SGPParameters)
    aux_gm = TC.center_aux_grid_mean(state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    FT = eltype(grid)
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set.TPS, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        Π = TD.exner(param_set.TPS, ts)
        z = grid.zc[k]
        aux_gm.dTdt[k] = APL.ARM_SGP_dTdt(FT)(t, z)
        aux_gm.dqtdt[k] = APL.ARM_SGP_dqtdt(FT)(Π, t, z)

    end
end

#####
##### GATE_III
#####

#TODO: Bug in original code (no zrough defined)
struct GATE_IIIParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    zrough::FT
    qsurface::FT
    cm::FT
    ch::FT
    cq::FT
    Tsurface::FT
    TPS::ThermodynamicsParameters{FT}
end

function GATE_IIIParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg", "Tg", "zrough", "qsurface","cm", "ch", "cq", "Tsurface"]
    (Pg, qtg, Tg, zrough, qsurface, cm, ch, cq, Tsurface) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "GATE_III") 
      
    return GATE_IIIParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,  
        zrough,
        qsurface,
        cm,
        ch,
        cq,
        Tsurface,
        TPS,
    )
    
end

ForcingBase(case::GATE_III, param_set::GATE_IIIParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::GATE_III, param_set::GATE_IIIParameters, namelist)
    Pg = param_set.Pg   #Pressure at ground
    Tg = param_set.Tg   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    qtg = param_set.qtg  #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{GATE_III}, grid::Grid, param_set::GATE_IIIParameters, state)
    p0 = TC.center_ref_state(state).p0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.q_tot[k] = APL.GATE_III_q_tot(FT)(z)
        aux_gm.T[k] = APL.GATE_III_T(FT)(z)
        prog_gm.u[k] = APL.GATE_III_u(FT)(z)
        ts = TD.PhaseEquil_pTq(param_set.TPS, p0[k], aux_gm.T[k], prog_gm.q_tot[k])
        prog_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp(param_set.TPS, ts)
        aux_gm.tke[k] = APL.GATE_III_tke(FT)(z)
    end
end

function surface_params(case::GATE_III, grid::TC.Grid, surf_ref_state, param_set::GATE_IIIParameters; kwargs...)
    p0_f_surf = TD.air_pressure(param_set.TPS, surf_ref_state)
    FT = eltype(p0_f_surf)

    qsurface = param_set.qsurface # kg/kg
    cm = param_set.cm
    ch = param_set.ch
    cq = param_set.cq
    Tsurface =  param_set.Tsurface

    # For GATE_III we provide values of transfer coefficients
    ts = TD.PhaseEquil_pθq(param_set.TPS, param_set, p0_f_surf, Tsurface, qsurface)
    qsurface = TD.q_vap_saturation(param_set.TPS, ts)
    kwargs = (; zrough, Tsurface, qsurface, cm, ch)
    return TC.FixedSurfaceCoeffs(FT; kwargs...)
end

function initialize_forcing(self::CasesBase{GATE_III}, grid::Grid, state, param_set::GATE_IIIParameters)
    aux_gm = TC.center_aux_grid_mean(state)
    for k in TC.real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.dqtdt[k] = APL.GATE_III_dqtdt(FT)(z)
        aux_gm.dTdt[k] = APL.GATE_III_dTdt(FT)(z)
    end
end

#####
##### DYCOMS_RF01
#####

struct DYCOMS_RF01Parameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    θ_surf::FT
    zrough::FT
    ustar::FT
    shf::FT
    lhf::FT
    Tsurface::FT
    qsurface::FT
    divergence::FT
    alpha_z::FT
    kappa::FT
    F0::FT
    F1::FT
    gas_constant::FT
    molmass_dryair::FT
    kappa_d::FT
    R_d::FT
    cp_d::FT
    TPS::ThermodynamicsParameters{FT}
end

function DYCOMS_RF01Parameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg", "θ_surf", "zrough", "ustar", "shf", "lhf",
               "Tsurface", "divergence", "alpha_z", "kappa", "F0", "F1",
               "gas_constant", "molmass_dryair", "kappa_d"
               ]
    (Pg, qtg, θ_surf, zrough, ustar, shf, lhf, Tsurface, qsurface, divergence, alpha_z, kappa, F0, F1,
     gas_constant, molmass_dryair, kappa_d) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "DYCOMS_RF01") 

    #derived parameters
    R_d = gas_constant / molmass_dryair
    cp_d = R_d / kappa_d
    
    return DYCOMS_RF01Parameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        θ_surf,  
        zrough,
        ustar,
        shf,
        lhf,
        Tsurface,
        qsurface,
        divergence,
        alpha_z,
        kappa,
        F0,
        F1,
        gas_constant,
        molmass_dryair,
        kappa_d,
        R_d,
        cp_d,
        TPS,
    )
    
end


function surface_ref_state(::DYCOMS_RF01, param_set::DYCOMS_RF01Parameters, namelist)
    Pg = param_set.Pg
    qtg = param_set.qtg
    θ_surf = param_set.θ_surf
    ts = TD.PhaseEquil_pθq(param_set.TPS, Pg, θ_surf, qtg)
    Tg = TD.air_temperature(param_set.TPS, ts)
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{DYCOMS_RF01}, grid::Grid, param_set::DYCOMS_RF01Parameters, state)
    FT = eltype(grid)
    p0 = TC.center_ref_state(state).p0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        # thetal profile as defined in DYCOMS
        z = grid.zc[k]
        prog_gm.q_tot[k] = APL.Dycoms_RF01_q_tot(FT)(z)
        prog_gm.θ_liq_ice[k] = APL.Dycoms_RF01_θ_liq_ice(FT)(z)
        # velocity profile (geostrophic)
        prog_gm.u[k] = APL.Dycoms_RF01_u0(FT)(z)
        prog_gm.v[k] = APL.Dycoms_RF01_v0(FT)(z)
        aux_gm.tke[k] = APL.Dycoms_RF01_tke(FT)(z)
    end
end

function surface_params(case::DYCOMS_RF01, grid::TC.Grid, surf_ref_state, param_set::DYCOMS_RF01Parameters; Ri_bulk_crit)
    FT = eltype(grid)
    zrough = param_set.zrough
    ustar = param_set.ustar # just to initilize grid mean covariances
    shf = param_set.shf # sensible heat flux
    lhf = param_set.lhf # latent heat flux
    Tsurface = param_set.Tsurface    # K      # i.e. the SST from DYCOMS setup
    qsurface = param_set.qsurface  # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?
    #density_surface  = 1.22     # kg/m^3

    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{DYCOMS_RF01}, grid::Grid, state, param_set::DYCOMS_RF01Parameters)
    aux_gm = TC.center_aux_grid_mean(state)

    # geostrophic velocity profiles
    parent(aux_gm.ug) .= 7.0
    parent(aux_gm.vg) .= -5.5

    # large scale subsidence
    divergence = self.Rad.divergence
    @inbounds for k in real_center_indices(grid)
        aux_gm.subsidence[k] = -grid.zc[k] * divergence
    end

    # no large-scale drying
    parent(aux_gm.dqtdt) .= 0 #kg/(kg * s)
end

function RadiationBase(case::DYCOMS_RF01, param_set::DYCOMS_RF01Parameters)
    return RadiationBase{Cases.get_radiation_type(case)}(;
        divergence = param_set.divergence,
        alpha_z = param_set.alpha_z,
        kappa = param_set.kappa,
        F0 = param_set.F0,
        F1 = param_set.F1,
    )
end

function initialize_radiation(self::CasesBase{DYCOMS_RF01}, grid::Grid, state, param_set::DYCOMS_RF01Parameters)
    aux_gm = TC.center_aux_grid_mean(state)

    # no large-scale drying
    parent(aux_gm.dqtdt_rad) .= 0 #kg/(kg * s)

    # Radiation based on eq. 3 in Stevens et. al., (2005)
    # cloud-top cooling + cloud-base warming + cooling in free troposphere
    update_radiation(self.Rad, grid, state, param_set)
end

#####
##### DYCOMS_RF02
#####
struct DYCOMS_RF02Parameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    θ_surf::FT
    zrough::FT
    ustar::FT
    shf::FT
    lhf::FT
    Tsurface::FT
    qsurface::FT
    divergence::FT
    alpha_z::FT
    kappa::FT
    F0::FT
    F1::FT
    gas_constant::FT
    molmass_dryair::FT
    kappa_d::FT
    R_d::FT
    cp_d::FT
    TPS::ThermodynamicsParameters{FT}
end

function DYCOMS_RF02Parameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg", "θ_surf", "zrough", "ustar", "shf", "lhf",
               "Tsurface", "divergence", "alpha_z", "kappa", "F0", "F1",
               "gas_constant", "molmass_dryair", "kappa_d"]
    (Pg, qtg, θ_surf, zrough, ustar, shf, lhf, Tsurface, qsurface, divergence, alpha_z, kappa, F0, F1,
     gas_constant, molmass_dryair, kappa_d) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "DYCOMS_RF02") 

    #derived parameters
    R_d = gas_constant / molmass_dryair
    cp_d = R_d / kappa_d
    
    return DYCOMS_RF02Parameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        θ_surf,  
        zrough,
        ustar,
        shf,
        lhf,
        Tsurface,
        qsurface,
        divergence,
        alpha_z,
        kappa,
        F0,
        F1,
        gas_constant,
        molmass_dryair,
        kappa_d,
        R_d,
        cp_d,
        TPS,
    )
    
end

function surface_ref_state(::DYCOMS_RF02, param_set::DYCOMS_RF02Parameters, namelist)
    Pg = param_set.Pg
    qtg = param_set.qtg
    θ_surf = param_set.θ_surf
    ts = TD.PhaseEquil_pθq(param_set.TPS, Pg, θ_surf, qtg)
    Tg = TD.air_temperature(param_set.TPS, ts)
    return TD.PhaseEquil_pTq(param_set.TPS,  Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{DYCOMS_RF02}, grid::Grid, param_set::DYCOMS_RF02Parameters, state)
    FT = eltype(grid)
    p0 = TC.center_ref_state(state).p0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        # θ_liq_ice profile as defined in DYCOM RF02
        z = grid.zc[k]
        prog_gm.q_tot[k] = APL.Dycoms_RF02_q_tot(FT)(z)
        prog_gm.θ_liq_ice[k] = APL.Dycoms_RF02_θ_liq_ice(FT)(z)
        # velocity profile
        prog_gm.v[k] = APL.Dycoms_RF02_v(FT)(z)
        prog_gm.u[k] = APL.Dycoms_RF02_u(FT)(z)
        aux_gm.tke[k] = APL.Dycoms_RF02_tke(FT)(z)
    end
end

function surface_params(case::DYCOMS_RF02, grid::TC.Grid, surf_ref_state, param_set::DYCOMS_RF02Parameters; Ri_bulk_crit)
    FT = eltype(grid)
    zrough = param_set.zrough  #TODO - not needed?
    ustar = param_set.ustar
    shf = param_set.shf# sensible heat flux
    lhf = param_set.lhf # latent heat flux
    Tsurface = param_set.Tsurface    # K      # i.e. the SST from DYCOMS setup
    qsurface = param_set.qsurface  # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?

    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{DYCOMS_RF02}, grid::Grid, state, param_set::DYCOMS_RF02Parameters)
    # the same as in DYCOMS_RF01
    aux_gm = TC.center_aux_grid_mean(state)

    # geostrophic velocity profiles
    parent(aux_gm.ug) .= 5.0
    parent(aux_gm.vg) .= -5.5

    # large scale subsidence
    divergence = self.Rad.divergence
    @inbounds for k in real_center_indices(grid)
        aux_gm.subsidence[k] = -grid.zc[k] * divergence
    end

    # no large-scale drying
    parent(aux_gm.dqtdt) .= 0 #kg/(kg * s)
end

function RadiationBase(case::DYCOMS_RF02, param_set::DYCOMS_RF02Parameters)
    return RadiationBase{Cases.get_radiation_type(case)}(;
        divergence = param_set.divergence,
        alpha_z = param_set.alpha_z,
        kappa = param_set.kappa,
        F0 = param_set.F0,
        F1 = param_set.F1,
    )
end

function initialize_radiation(self::CasesBase{DYCOMS_RF02}, grid::Grid, state, param_set::DYCOMS_RF02Parameters)
    # the same as in DYCOMS_RF01
    aux_gm = TC.center_aux_grid_mean(state)

    # no large-scale drying
    parent(aux_gm.dqtdt_rad) .= 0 #kg/(kg * s)

    # Radiation based on eq. 3 in Stevens et. al., (2005)
    # cloud-top cooling + cloud-base warming + cooling in free troposphere
    update_radiation(self.Rad, grid, state, param_set)
end

#####
##### GABLS
#####

struct GABLSParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    zrough::FT
    shf::FT
    lhf::FT
    qsurface::FT
    coriolis_param::FT
    TPS::ThermodynamicsParameters{FT}
end

function GABLSParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg","Tg",  "zrough", "shf", "lhf", "qsurface", "coriolis_param"]
    (Pg, qtg, Tg,  zrough, shf, lhf, qsurface, coriolis_param) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "GABLS") 
    
    return GABLSParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        zrough,
        shf,
        lhf,
        qsurface,
        coriolis_param,
        TPS,
    )
    
end


function ForcingBase(case::GABLS, param_set::GABLSParameters; kwargs...)
    # Omega = CPP.Omega(param_set)
    # coriolis_param = 2 * Omega * sin(latitude * π / 180 ) # s^{-1}
    ForcingBase(
        get_forcing_type(case);
        apply_coriolis = true,
        apply_subsidence = false,
        coriolis_param = param_set.coriolis_param,
    )
end

function surface_ref_state(::GABLS, param_set::GABLSParameters, namelist)
    Pg = param_set.Pg  #Pressure at ground,
    Tg = param_set.Tg  #Temperature at ground,
    qtg = param_set.qtg
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{GABLS}, grid::Grid, param_set::GABLSParameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    FT = eltype(grid)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        #Set wind velocity profile
        prog_gm.u[k] = APL.GABLS_u(FT)(z)
        prog_gm.v[k] = APL.GABLS_v(FT)(z)
        prog_gm.θ_liq_ice[k] = APL.GABLS_θ_liq_ice(FT)(z)
        prog_gm.q_tot[k] = APL.GABLS_q_tot(FT)(z)
        aux_gm.tke[k] = APL.GABLS_tke(FT)(z)
        aux_gm.Hvar[k] = aux_gm.tke[k]
    end
end

function surface_params(case::GABLS, grid::TC.Grid, surf_ref_state, param_set::GABLSParameters; kwargs...)
    FT = eltype(grid)
    Tsurface = t -> 265.0 - (0.25 / 3600.0) * t
    qsurface = param_set.qsurface
    shf = param_set.shf  # only prevent zero division in SF.jl lmo
    lhf = param_set.lhf # only prevent zero division in SF.jl lmo
    zrough = param_set.zrough

    kwargs = (; Tsurface, qsurface, shf, lhf, zrough)
    return TC.MoninObukhovSurface(FT; kwargs...)
end

function initialize_forcing(self::CasesBase{GABLS}, grid::Grid, state, param_set::GABLSParameters)
    FT = eltype(grid)
    initialize(self.Fo, grid, state)
    aux_gm = TC.center_aux_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        # Geostrophic velocity profiles.
        z = grid.zc[k]
        aux_gm.ug[k] = APL.GABLS_geostrophic_ug(FT)(z)
        aux_gm.vg[k] = APL.GABLS_geostrophic_vg(FT)(z)
    end
    return nothing
end

#####
##### SP
#####

# Not fully implemented yet - Ignacio
struct SPParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    coriolis_param::FT
    TPS::ThermodynamicsParameters{FT}
end

function SPParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg","Tg", "coriolis_param"]
    (Pg, qtg, Tg, coriolis_param) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "SP") 
    
    return SPParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        coriolis_param,
        TPS,
    )
    
end

function ForcingBase(case::SP, param_set::SPParameters; kwargs...)
    # Omega = CPP.Omega(param_set)
    # coriolis_param = 2 * Omega * sin(latitude * π / 180 ) # s^{-1}
    ForcingBase(
        get_forcing_type(case);
        apply_coriolis = true,
        apply_subsidence = false,
        coriolis_param = param_set.coriolis_param, #s^{-1}
    )
end

function surface_ref_state(::SP, param_set::SPParameters, namelist)
    Pg = param_set.Pg  #Pressure at ground
    Tg = param_set.Tg  #Temperature at ground
    qtg = param_set.qtg   #Total water mixing ratio at TC. if set to 0, alpha0, rho0, p0 are NaN.
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{SP}, grid::Grid, param_set::SPParameters, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    FT = eltype(grid)
    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.u[k] = APL.SP_u(FT)(z)
        prog_gm.v[k] = APL.SP_v(FT)(z)
        prog_gm.θ_liq_ice[k] = APL.SP_θ_liq_ice(FT)(z)
        prog_gm.q_tot[k] = APL.SP_q_tot(FT)(z)
        aux_gm.tke[k] = APL.SP_tke(FT)(z)
    end
end

function initialize_forcing(self::CasesBase{SP}, grid::Grid, state, param_set::SPParameters)
    aux_gm = TC.center_aux_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.ug[k] = APL.SP_geostrophic_u(FT)(z)
        aux_gm.vg[k] = APL.SP_geostrophic_v(FT)(z)
    end
end

#####
##### DryBubble
#####

struct DryBubbleParameters{FT} <: AbstractCaseParameters
    Pg::FT
    qtg::FT
    Tg::FT
    Tsurface::FT
    qsurface::FT
    shf::FT
    lhf::FT
    ustar::FT
    TPS::ThermodynamicsParameters{FT}
end

function DryBubbleParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["Pg","qtg","Tg",  "Tsurface", "qsurface", "shf", "lhf", "ustar"]
    (Pg, qtg, Tg, Tsurface, qsurface, shf, lhf, ustar) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "DryBubble") 
    
    return DryBubbleParameters{CLIMAParameters.get_parametric_type(param_set)}(
        Pg,
        qtg,
        Tg,
        Tsurface,
        qsurface,
        shf,
        lhf,
        ustar,
        TPS,
    )
    
end



ForcingBase(case::DryBubble, param_set::DryBubbleParameters; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::DryBubble, param_set::DryBubbleParameters, namelist)
    Pg = param_set.Pg  #Pressure at ground
    Tg = param_set.Tg
    qtg = param_set.qtg
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{DryBubble}, grid::Grid, param_set::DryBubbleParameters, state)

    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    # initialize Grid Mean Profiles of thetali and qt
    zc_in = grid.zc
    FT = eltype(grid)
    prof_θ_liq_ice = APL.DryBubble_θ_liq_ice(FT)
    prog_gm.θ_liq_ice .= prof_θ_liq_ice.(zc_in)
    parent(prog_gm.u) .= 0.01
    parent(prog_gm.q_tot) .= 0
    parent(aux_gm.tke) .= 0
    parent(aux_gm.Hvar) .= 0
    parent(aux_gm.QTvar) .= 0
    parent(aux_gm.HQTcov) .= 0
end

function surface_params(case::DryBubble, grid::TC.Grid, surf_ref_state, param_set::DryBubbleParameters; Ri_bulk_crit)
    FT = eltype(grid)
    Tsurface = param_set.Tsurface
    qsurface = param_set.qsurface
    shf = param_set.shf # only prevent zero division in SF.jl lmo
    lhf = param_set.lhf # only prevent zero division in SF.jl lmo
    ustar = param_set.ustar

    kwargs = (; Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

#####
##### LES_driven_SCM
#####

struct LES_driven_SCMParameters{FT} <: AbstractCaseParameters
    zrough::FT
    ustar::FT
    coriolis_param::FT
    TPS::ThermodynamicsParameters{FT}
end

function LES_driven_SCMParameters(
    param_set,
    ThermodynamicsParameters{FT}
) where {FT}
    aliases = ["zrough","ustar","coriolis_param"]
    (zrough,ustar,coriolis_param) =
        CLIMAParameters.get_parameter_values!(param_set, aliases, "LES_driven_SCM") 
    
    return LES_driven_SCMParameters{CLIMAParameters.get_parametric_type(param_set)}(
        zrough,
        ustar,
        coriolis_param,
        TPS,
    )
    
end




forcing_kwargs(::LES_driven_SCM, namelist) = (; nudge_tau = namelist["forcing"]["nudging_timescale"])

function surface_param_kwargs(::LES_driven_SCM, namelist)
    les_filename = namelist["meta"]["lesfile"]
    # load data here
    LESDat = NC.Dataset(les_filename, "r") do data
        t = data.group["profiles"]["t"][:]
        t_interval_from_end_s = namelist["t_interval_from_end_s"]
        t_from_end_s = t .- t[end]
        # find inds within time interval
        time_interval_bool = findall(>(-t_interval_from_end_s), t_from_end_s)
        imin = time_interval_bool[1]
        imax = time_interval_bool[end]

        TC.LESData(imin, imax, les_filename, t_interval_from_end_s, namelist["initial_condition_averaging_window_s"])
    end
    return (; LESDat)
end

function ForcingBase(case::LES_driven_SCM, param_set::LES_driven_SCMParameters; nudge_tau)
    ForcingBase(
        get_forcing_type(case);
        apply_coriolis = false,
        apply_subsidence = true,
        coriolis_param = param_set.coriolis_param,
        nudge_tau = nudge_tau,
    )
end

function surface_ref_state(::LES_driven_SCM, param_set::LES_driven_SCMParameters, namelist)
    les_filename = namelist["meta"]["lesfile"]

    Pg, Tg, qtg = NC.Dataset(les_filename, "r") do data
        Pg = data.group["reference"]["p0_full"][1] #Pressure at ground
        Tg = data.group["reference"]["temperature0"][1] #Temperature at ground
        ql_ground = data.group["reference"]["ql0"][1]
        qv_ground = data.group["reference"]["qv0"][1]
        qi_ground = data.group["reference"]["qi0"][1]
        qtg = ql_ground + qv_ground + qi_ground #Total water mixing ratio at surface
        (Pg, Tg, qtg)
    end
    return TD.PhaseEquil_pTq(param_set.TPS, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{LES_driven_SCM}, grid::Grid, param_set::LES_driven_SCMParameters, state)

    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    NC.Dataset(self.LESDat.les_filename, "r") do data
        t = data.group["profiles"]["t"][:]
        # define time interval
        half_window = self.LESDat.initial_condition_averaging_window_s / 2
        t_window_start = (self.LESDat.t_interval_from_end_s + half_window)
        t_window_end = (self.LESDat.t_interval_from_end_s - half_window)
        t_from_end_s = Array(t) .- t[end]
        # find inds within time interval
        in_window(x) = -t_window_start <= x <= -t_window_end
        time_interval_bool = findall(in_window, t_from_end_s)
        imin = time_interval_bool[1]
        imax = time_interval_bool[end]
        zc_les = Array(TC.get_nc_data(data, "zc"))
        parent(prog_gm.θ_liq_ice) .=
            pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "thetali_mean", imin, imax))
        parent(prog_gm.q_tot) .= pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "qt_mean", imin, imax))
        parent(prog_gm.u) .= pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "u_mean", imin, imax))
        parent(prog_gm.v) .= pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "v_mean", imin, imax))
    end
    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.tke[k] = if z <= 2500.0
            1.0 - z / 3000.0
        else
            0.0
        end
    end
end

function surface_params(case::LES_driven_SCM, grid::TC.Grid, surf_ref_state, param_set::LES_driven_SCMParameters; Ri_bulk_crit, LESDat)
    FT = eltype(grid)
    nt = NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zrough = param_set.zrough
        Tsurface = Statistics.mean(data.group["timeseries"]["surface_temperature"][:][imin:imax], dims = 1)[1]
        # get surface value of q
        mean_qt_prof = Statistics.mean(data.group["profiles"]["qt_mean"][:][:, imin:imax], dims = 2)[:]
        field = TC.FieldFromNamedTuple(TC.face_space(grid), (; q_tot = FT(0)))
        Ic = CCO.InterpolateF2C()
        q_tot_c = Ic.(field.q_tot)
        qsurface = q_tot_c[TC.kc_surface(grid)]
        lhf = Statistics.mean(data.group["timeseries"]["lhf_surface_mean"][:][imin:imax], dims = 1)[1]
        shf = Statistics.mean(data.group["timeseries"]["shf_surface_mean"][:][imin:imax], dims = 1)[1]
        (; zrough, Tsurface, qsurface, lhf, shf)
    end
    UnPack.@unpack zrough, Tsurface, qsurface, lhf, shf = nt

    ustar = param_set.ustar # TODO: why is initialization missing?
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

initialize_forcing(self::CasesBase{LES_driven_SCM}, grid::Grid, state, param_set::LES_driven_SCMParameters) =
    initialize(self.Fo, grid, state, self.LESDat)

initialize_radiation(self::CasesBase{LES_driven_SCM}, grid::Grid, state, param_set::LES_driven_SCMParameters) =
    initialize(self.Rad, grid, state, self.LESDat)

end # module Cases
