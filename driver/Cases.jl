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

RadiationBase(case::AbstractCaseType) = RadiationBase{Cases.get_radiation_type(case)}()

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

ForcingBase(case::Soares, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::Soares, param_set::APS, namelist)
    Pg = 1000.0 * 100.0
    qtg = 5.0e-3
    Tg = 300.0
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{Soares}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0
    FT = eltype(grid)
    prof_q_tot = APL.Soares_q_tot(FT)
    prof_θ_liq_ice = APL.Soares_θ_liq_ice(FT)
    prof_u = APL.Soares_u(FT)
    prof_tke = APL.Soares_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::Soares, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = 0.16 #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    Tsurface = 300.0
    qsurface = 5.0e-3
    θ_flux = 6.0e-2
    qt_flux = 2.5e-5
    ts = TD.PhaseEquil_pTq(param_set, p0_f_surf, Tsurface, qsurface)
    lhf = qt_flux * ρ0_f_surf * TD.latent_heat_vapor(param_set, ts)
    shf = θ_flux * TD.cp_m(param_set, ts) * ρ0_f_surf
    ustar = 0.28 # just to initilize grid mean covariances
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

#####
##### Nieuwstadt
#####

ForcingBase(case::Nieuwstadt, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::Nieuwstadt, param_set::APS, namelist)
    Pg = 1000.0 * 100.0
    Tg = 300.0
    qtg = 1.0e-12 # Total water mixing ratio
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{Nieuwstadt}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0

    FT = eltype(grid)
    prof_θ_liq_ice = APL.Nieuwstadt_θ_liq_ice(FT)
    prof_u = APL.Nieuwstadt_u(FT)
    prof_tke = APL.Nieuwstadt_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::Nieuwstadt, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = 0.16 #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    Tsurface = 300.0
    qsurface = 0.0
    θ_flux = 6.0e-2
    lhf = 0.0 # It would be 0.0 if we follow Nieuwstadt.
    ts = TD.PhaseEquil_pTq(param_set, p0_f_surf, Tsurface, qsurface)
    shf = θ_flux * TD.cp_m(param_set, ts) * ρ0_f_surf
    ustar = 0.28 # just to initilize grid mean covariances
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

#####
##### Bomex
#####

ForcingBase(case::Bomex, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = true, apply_subsidence = true, coriolis_param = 0.376e-4) #= s^{-1} =#

function surface_ref_state(::Bomex, param_set::APS, namelist)
    Pg = 1.015e5 #Pressure at ground
    Tg = 300.4 #Temperature at ground
    qtg = 0.02245#Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{Bomex}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0

    FT = eltype(grid)
    prof_q_tot = APL.Bomex_q_tot(FT)
    prof_θ_liq_ice = APL.Bomex_θ_liq_ice(FT)
    prof_u = APL.Bomex_u(FT)
    prof_tke = APL.Bomex_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::Bomex, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = 1.0e-4
    qsurface = 22.45e-3 # kg/kg
    θ_surface = 299.1
    θ_flux = 8.0e-3
    qt_flux = 5.2e-5
    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set, ts)
    lhf = qt_flux * ρ0_f_surf * TD.latent_heat_vapor(param_set, ts)
    shf = θ_flux * TD.cp_m(param_set, ts) * ρ0_f_surf
    ustar = 0.28 # m/s
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{Bomex}, grid::Grid, state, param_set)
    initialize(self.Fo, grid, state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    ts_gm = aux_gm.ts

    FT = eltype(grid)
    prof_ug = APL.Bomex_geostrophic_u(FT)
    prof_dTdt = APL.Bomex_dTdt(FT)
    prof_dqtdt = APL.Bomex_dqtdt(FT)
    prof_subsidence = APL.Bomex_subsidence(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        # Geostrophic velocity profiles. vg = 0
        aux_gm.ug[k] = prof_ug(z)
        Π = TD.exner(param_set, ts_gm[k])
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

ForcingBase(case::life_cycle_Tan2018, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = true, apply_subsidence = true, coriolis_param = 0.376e-4) #= s^{-1} =#

function surface_ref_state(::life_cycle_Tan2018, param_set::APS, namelist)
    Pg = 1.015e5  #Pressure at ground
    Tg = 300.4  #Temperature at ground
    qtg = 0.02245   #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{life_cycle_Tan2018}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0

    FT = eltype(grid)
    prof_q_tot = APL.LifeCycleTan2018_q_tot(FT)
    prof_θ_liq_ice = APL.LifeCycleTan2018_θ_liq_ice(FT)
    prof_u = APL.LifeCycleTan2018_u(FT)
    prof_tke = APL.LifeCycleTan2018_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        prog_gm.u[k] = prof_u(z)
        aux_gm.tke[k] = prof_tke(z)
    end
end

function life_cycle_buoyancy_flux(param_set, weight = 1)
    g = CPP.grav(param_set)
    molmass_ratio = CPP.molmass_ratio(param_set)
    return g * (
        (8.0e-3 * weight + (molmass_ratio - 1) * (299.1 * 5.2e-5 * weight + 22.45e-3 * 8.0e-3 * weight)) /
        (299.1 * (1 + (molmass_ratio - 1) * 22.45e-3))
    )
end

function surface_params(case::life_cycle_Tan2018, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    ρ0_f_surf = TD.air_density(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)
    zrough = 1.0e-4 # not actually used, but initialized to reasonable value
    qsurface = 22.45e-3 # kg/kg
    θ_surface = 299.1
    θ_flux = 8.0e-3
    qt_flux = 5.2e-5
    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set, ts)
    lhf0 = qt_flux * ρ0_f_surf * TD.latent_heat_vapor(param_set, ts)
    shf0 = θ_flux * TD.cp_m(param_set, ts) * ρ0_f_surf

    weight_factor(t) = 0.01 + 0.99 * (cos(2.0 * π * t / 3600.0) + 1.0) / 2.0
    weight = 1.0
    lhf = t -> lhf0 * (weight * weight_factor(t))
    shf = t -> shf0 * (weight * weight_factor(t))

    ustar = 0.28 # m/s
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{life_cycle_Tan2018}, grid::Grid, state, param_set)
    initialize(self.Fo, grid, state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    ts_gm = aux_gm.ts

    FT = eltype(grid)
    prof_ug = APL.LifeCycleTan2018_geostrophic_u(FT)
    prof_dTdt = APL.LifeCycleTan2018_dTdt(FT)
    prof_dqtdt = APL.LifeCycleTan2018_dqtdt(FT)
    prof_subsidence = APL.LifeCycleTan2018_subsidence(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        # Geostrophic velocity profiles. vg = 0
        aux_gm.ug[k] = prof_ug(z)
        Π = TD.exner(param_set, ts_gm[k])
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

function ForcingBase(case::Rico, param_set::APS; kwargs...)
    latitude = 18.0
    Omega = CPP.Omega(param_set)
    return ForcingBase(
        get_forcing_type(case);
        apply_coriolis = true,
        apply_subsidence = true,
        coriolis_param = 2.0 * Omega * sin(latitude * π / 180.0),
    ) #= s^{-1} =#
end

function surface_ref_state(::Rico, param_set::APS, namelist)
    molmass_ratio = CPP.molmass_ratio(param_set)
    Pg = 1.0154e5  #Pressure at ground
    Tg = 299.8  #Temperature at ground
    pvg = TD.saturation_vapor_pressure(param_set, Tg, TD.Liquid())
    qtg = (1 / molmass_ratio) * pvg / (Pg - pvg)   #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{Rico}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_tc = TC.center_aux_turbconv(state)
    p0 = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0

    FT = eltype(grid)
    prof_u = APL.Rico_u(FT)
    prof_v = APL.Rico_v(FT)
    prof_q_tot = APL.Rico_q_tot(FT)
    prof_θ_liq_ice = APL.Rico_θ_liq_ice(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.u[k] = prof_u(z)
        prog_gm.v[k] = prof_v(z)
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
    end

    # Need to get θ_virt
    @inbounds for k in real_center_indices(grid)
        # Thermo state field cache is not yet
        # defined, so we can't use it yet.
        ts = TD.PhaseEquil_pθq(param_set, p0[k], aux_gm.θ_liq_ice[k], aux_gm.q_tot[k])
        aux_gm.θ_virt[k] = TD.virtual_pottemp(param_set, ts)
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

function surface_params(case::Rico, zc_surf, surf_ref_state, param_set; kwargs...)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)

    zrough = 0.00015
    cm = 0.001229
    ch = 0.001094
    cq = 0.001133
    # Adjust for non-IC grid spacing
    grid_adjust = (log(20.0 / zrough) / log(zc_surf / zrough))^2
    cm = cm * grid_adjust
    ch = ch * grid_adjust
    cq = cq * grid_adjust # TODO: not yet used..
    Tsurface = 299.8

    # For Rico we provide values of transfer coefficients
    ts = TD.PhaseEquil_pTq(param_set, p0_f_surf, Tsurface, FT(0)) # TODO: is this correct?
    qsurface = TD.q_vap_saturation(param_set, ts)
    kwargs = (; zrough, Tsurface, qsurface, cm, ch)
    return TC.FixedSurfaceCoeffs(FT; kwargs...)
end

function initialize_forcing(self::CasesBase{Rico}, grid::Grid, state, param_set)
    initialize(self.Fo, grid, state)
    p0_c = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    ts_gm = aux_gm.ts

    FT = eltype(grid)
    prof_ug = APL.Rico_geostrophic_ug(FT)
    prof_vg = APL.Rico_geostrophic_vg(FT)
    prof_dTdt = APL.Rico_dTdt(FT)
    prof_dqtdt = APL.Rico_dqtdt(FT)
    prof_subsidence = APL.Rico_subsidence(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        Π = TD.exner(param_set, ts_gm[k])
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

ForcingBase(case::TRMM_LBA, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false, kwargs...)

function surface_ref_state(::TRMM_LBA, param_set::APS, namelist)
    molmass_ratio = CPP.molmass_ratio(param_set)
    Pg = 991.3 * 100  #Pressure at ground
    Tg = 296.85   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    pvg = TD.saturation_vapor_pressure(param_set, Tg, TD.Liquid())
    qtg = (1 / molmass_ratio) * pvg / (Pg - pvg) #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{TRMM_LBA}, grid::Grid, param_set, state)
    p0 = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
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
    molmass_ratio = CPP.molmass_ratio(param_set)
    prog_gm = TC.center_prog_grid_mean(state)

    prog_gm.u .= prof_u.(zc)
    prog_gm.v .= prof_v.(zc)
    aux_gm.T .= prof_T.(zc)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        pv_star = TD.saturation_vapor_pressure(param_set, aux_gm.T[k], TD.Liquid())
        # eq. 37 in pressel et al and the def of RH
        RH = prof_RH(z)
        denom = (prof_p(z) - pv_star + (1 / molmass_ratio) * pv_star * RH / 100.0)
        qv_star = pv_star * (1 / molmass_ratio) / denom
        aux_gm.q_tot[k] = qv_star * RH / 100
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        phase_part = TD.PhasePartition(aux_gm.q_tot[k], 0.0, 0.0) # initial state is not saturated
        aux_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp_given_pressure(param_set, aux_gm.T[k], p0[k], phase_part)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::TRMM_LBA, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)
    # zrough = 1.0e-4 # not actually used, but initialized to reasonable value
    qsurface = 22.45e-3 # kg/kg
    θ_surface = (273.15 + 23)
    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set, ts)
    ustar = 0.28 # this is taken from Bomex -- better option is to approximate from LES tke above the surface
    lhf = t -> 554.0 * max(0, cos(π / 2 * ((5.25 * 3600.0 - t) / 5.25 / 3600.0)))^1.3
    shf = t -> 270.0 * max(0, cos(π / 2 * ((5.25 * 3600.0 - t) / 5.25 / 3600.0)))^1.5
    kwargs = (; Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit, zero_uv_fluxes = true)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function forcing_kwargs(::TRMM_LBA, namelist)
    return (; rad = APL.TRMM_LBA_radiation(Float64))
end

function initialize_forcing(self::CasesBase{TRMM_LBA}, grid::Grid, state, param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    rad = self.Fo.rad
    @inbounds for k in real_center_indices(grid)
        aux_gm.dTdt[k] = rad(0, grid.zc[k].z)
    end
    return nothing
end

function update_forcing(self::CasesBase{TRMM_LBA}, grid, state, t::Real, param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    rad = self.Fo.rad
    @inbounds for k in real_center_indices(grid)
        aux_gm.dTdt[k] = rad(t, grid.zc[k].z)
    end
end

#####
##### ARM_SGP
#####

ForcingBase(case::ARM_SGP, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = true, apply_subsidence = false, coriolis_param = 8.5e-5)

function surface_ref_state(::ARM_SGP, param_set::APS, namelist)
    Pg = 970.0 * 100 #Pressure at ground
    Tg = 299.0   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    qtg = 15.2 / 1000 #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{ARM_SGP}, grid::Grid, param_set, state)
    # ARM_SGP inputs
    p0 = TC.center_ref_state(state).p0
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0

    FT = eltype(grid)
    prof_u = APL.ARM_SGP_u(FT)
    prof_q_tot = APL.ARM_SGP_q_tot(FT)
    prof_θ_liq_ice = APL.ARM_SGP_θ_liq_ice(FT)
    prof_tke = APL.ARM_SGP_tke(FT)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        # TODO figure out how to use ts here
        phase_part = TD.PhasePartition(aux_gm.q_tot[k], aux_gm.q_liq[k], 0.0)
        Π = TD.exner_given_pressure(param_set, p0[k], phase_part)
        prog_gm.u[k] = prof_u(z)
        aux_gm.q_tot[k] = prof_q_tot(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.T[k] = prof_θ_liq_ice(z) * Π
        aux_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp_given_pressure(param_set, aux_gm.T[k], p0[k], phase_part)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.tke[k] = prof_tke(z)
    end
end

function surface_params(case::ARM_SGP, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)
    qsurface = 15.2e-3 # kg/kg
    θ_surface = 299.0
    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, θ_surface, qsurface)
    Tsurface = TD.air_temperature(param_set, ts)
    ustar = 0.28 # this is taken from Bomex -- better option is to approximate from LES tke above the surface

    t_Sur_in = arr_type([0.0, 4.0, 6.5, 7.5, 10.0, 12.5, 14.5]) .* 3600 #LES time is in sec
    SH = arr_type([-30.0, 90.0, 140.0, 140.0, 100.0, -10, -10]) # W/m^2
    LH = arr_type([5.0, 250.0, 450.0, 500.0, 420.0, 180.0, 0.0]) # W/m^2
    shf = Dierckx.Spline1D(t_Sur_in, SH; k = 1)
    lhf = Dierckx.Spline1D(t_Sur_in, LH; k = 1)

    kwargs = (; Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit, zero_uv_fluxes = true)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{ARM_SGP}, grid::Grid, state, param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        aux_gm.ug[k] = 10.0
        aux_gm.vg[k] = 0.0
    end
    return nothing
end

function update_forcing(self::CasesBase{ARM_SGP}, grid, state, t::Real, param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    p0_c = TC.center_ref_state(state).p0
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    FT = eltype(grid)
    @inbounds for k in real_center_indices(grid)
        Π = TD.exner(param_set, ts_gm[k])
        z = grid.zc[k]
        aux_gm.dTdt[k] = APL.ARM_SGP_dTdt(FT)(t, z)
        aux_gm.dqtdt[k] = APL.ARM_SGP_dqtdt(FT)(Π, t, z)

    end
end

#####
##### GATE_III
#####

ForcingBase(case::GATE_III, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::GATE_III, param_set::APS, namelist)
    Pg = 1013.0 * 100  #Pressure at ground
    Tg = 299.184   # surface values for reference state (RS) which outputs p0 rho0 alpha0
    qtg = 16.5 / 1000 #Total water mixing ratio at surface
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{GATE_III}, grid::Grid, param_set, state)
    FT = eltype(grid)
    p0 = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        aux_gm.q_tot[k] = APL.GATE_III_q_tot(FT)(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.T[k] = APL.GATE_III_T(FT)(z)
        prog_gm.u[k] = APL.GATE_III_u(FT)(z)
        ts = TD.PhaseEquil_pTq(param_set, p0[k], aux_gm.T[k], aux_gm.q_tot[k])
        aux_gm.θ_liq_ice[k] = TD.liquid_ice_pottemp(param_set, ts)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.tke[k] = APL.GATE_III_tke(FT)(z)
    end
end

function surface_params(case::GATE_III, zc_surf, surf_ref_state, param_set; kwargs...)
    p0_f_surf = TD.air_pressure(param_set, surf_ref_state)
    FT = eltype(p0_f_surf)

    qsurface = 16.5 / 1000.0 # kg/kg
    cm = 0.0012
    ch = 0.0034337
    cq = 0.0034337
    Tsurface = 299.184

    # For GATE_III we provide values of transfer coefficients
    ts = TD.PhaseEquil_pθq(param_set, p0_f_surf, Tsurface, qsurface)
    qsurface = TD.q_vap_saturation(param_set, ts)
    kwargs = (; Tsurface, qsurface, cm, ch)
    return TC.FixedSurfaceCoeffs(FT; kwargs...)
end

function initialize_forcing(self::CasesBase{GATE_III}, grid::Grid, state, param_set)
    FT = eltype(grid)
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

function surface_ref_state(::DYCOMS_RF01, param_set::APS, namelist)
    Pg = 1017.8 * 100.0
    qtg = 9.0 / 1000.0
    θ_surf = 289.0
    ts = TD.PhaseEquil_pθq(param_set, Pg, θ_surf, qtg)
    Tg = TD.air_temperature(param_set, ts)
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{DYCOMS_RF01}, grid::Grid, param_set, state)
    FT = eltype(grid)
    p0 = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        # thetal profile as defined in DYCOMS
        z = grid.zc[k]
        aux_gm.q_tot[k] = APL.Dycoms_RF01_q_tot(FT)(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.θ_liq_ice[k] = APL.Dycoms_RF01_θ_liq_ice(FT)(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        # velocity profile (geostrophic)
        prog_gm.u[k] = APL.Dycoms_RF01_u0(FT)(z)
        prog_gm.v[k] = APL.Dycoms_RF01_v0(FT)(z)
        aux_gm.tke[k] = APL.Dycoms_RF01_tke(FT)(z)
    end
end

function surface_params(case::DYCOMS_RF01, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    FT = eltype(surf_ref_state)
    zrough = 1.0e-4
    ustar = 0.28 # just to initilize grid mean covariances
    shf = 15.0 # sensible heat flux
    lhf = 115.0 # latent heat flux
    Tsurface = 292.5    # K      # i.e. the SST from DYCOMS setup
    qsurface = 13.84e-3 # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?
    #density_surface  = 1.22     # kg/m^3

    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{DYCOMS_RF01}, grid::Grid, state, param_set)
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

function RadiationBase(case::DYCOMS_RF01)
    return RadiationBase{Cases.get_radiation_type(case)}(;
        divergence = 3.75e-6,
        alpha_z = 1.0,
        kappa = 85.0,
        F0 = 70.0,
        F1 = 22.0,
    )
end

function initialize_radiation(self::CasesBase{DYCOMS_RF01}, grid::Grid, state, param_set)
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

function surface_ref_state(::DYCOMS_RF02, param_set::APS, namelist)
    Pg = 1017.8 * 100.0
    qtg = 9.0 / 1000.0
    θ_surf = 288.3
    ts = TD.PhaseEquil_pθq(param_set, Pg, θ_surf, qtg)
    Tg = TD.air_temperature(param_set, ts)
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{DYCOMS_RF02}, grid::Grid, param_set, state)
    FT = eltype(grid)
    p0 = TC.center_ref_state(state).p0
    ρ0_c = TC.center_ref_state(state).ρ0
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    @inbounds for k in real_center_indices(grid)
        # θ_liq_ice profile as defined in DYCOM RF02
        z = grid.zc[k]
        aux_gm.q_tot[k] = APL.Dycoms_RF02_q_tot(FT)(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.θ_liq_ice[k] = APL.Dycoms_RF02_θ_liq_ice(FT)(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        # velocity profile
        prog_gm.v[k] = APL.Dycoms_RF02_v(FT)(z)
        prog_gm.u[k] = APL.Dycoms_RF02_u(FT)(z)
        aux_gm.tke[k] = APL.Dycoms_RF02_tke(FT)(z)
    end
end

function surface_params(case::DYCOMS_RF02, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    FT = eltype(surf_ref_state)
    zrough = 1.0e-4  #TODO - not needed?
    ustar = 0.25
    shf = 16.0 # sensible heat flux
    lhf = 93.0 # latent heat flux
    Tsurface = 292.5    # K      # i.e. the SST from DYCOMS setup
    qsurface = 13.84e-3 # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?

    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

function initialize_forcing(self::CasesBase{DYCOMS_RF02}, grid::Grid, state, param_set)
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

function RadiationBase(case::DYCOMS_RF02)
    return RadiationBase{Cases.get_radiation_type(case)}(;
        divergence = 3.75e-6,
        alpha_z = 1.0,
        kappa = 85.0,
        F0 = 70.0,
        F1 = 22.0,
    )
end

function initialize_radiation(self::CasesBase{DYCOMS_RF02}, grid::Grid, state, param_set)
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

function ForcingBase(case::GABLS, param_set::APS; kwargs...)
    coriolis_param = 1.39e-4 # s^{-1}
    # Omega = CPP.Omega(param_set)
    # coriolis_param = 2 * Omega * sin(latitude * π / 180 ) # s^{-1}
    ForcingBase(
        get_forcing_type(case);
        apply_coriolis = true,
        apply_subsidence = false,
        coriolis_param = coriolis_param,
    )
end

function surface_ref_state(::GABLS, param_set::APS, namelist)
    Pg = 1.0e5  #Pressure at ground,
    Tg = 265.0  #Temperature at ground,
    qtg = 0.0
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end
function initialize_profiles(self::CasesBase{GABLS}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0
    FT = eltype(grid)

    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        #Set wind velocity profile
        prog_gm.u[k] = APL.GABLS_u(FT)(z)
        prog_gm.v[k] = APL.GABLS_v(FT)(z)
        aux_gm.θ_liq_ice[k] = APL.GABLS_θ_liq_ice(FT)(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.q_tot[k] = APL.GABLS_q_tot(FT)(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.tke[k] = APL.GABLS_tke(FT)(z)
        aux_gm.Hvar[k] = aux_gm.tke[k]
    end
end

function surface_params(case::GABLS, zc_surf, surf_ref_state, param_set; kwargs...)
    FT = eltype(surf_ref_state)
    Tsurface = t -> 265.0 - (0.25 / 3600.0) * t
    qsurface = 0.0
    shf = 0.0001 # only prevent zero division in SF.jl lmo
    lhf = 0.0001 # only prevent zero division in SF.jl lmo
    # ustar = 0.1 # TODO: remove, this isn't actually used
    zrough = 0.1

    kwargs = (; Tsurface, qsurface, shf, lhf, zrough)
    return TC.MoninObukhovSurface(FT; kwargs...)
end

function initialize_forcing(self::CasesBase{GABLS}, grid::Grid, state, param_set)
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
function ForcingBase(case::SP, param_set::APS; kwargs...)
    coriolis_param = 1.0e-4 # s^{-1}
    # Omega = CPP.Omega(param_set)
    # coriolis_param = 2 * Omega * sin(latitude * π / 180 ) # s^{-1}
    ForcingBase(
        get_forcing_type(case);
        apply_coriolis = true,
        apply_subsidence = false,
        coriolis_param = coriolis_param,
    )
end

function surface_ref_state(::SP, param_set::APS, namelist)
    Pg = 1.0e5  #Pressure at ground
    Tg = 300.0  #Temperature at ground
    qtg = 1.0e-4   #Total water mixing ratio at TC. if set to 0, alpha0, rho0, p0 are NaN.
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{SP}, grid::Grid, param_set, state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0
    FT = eltype(grid)
    @inbounds for k in real_center_indices(grid)
        z = grid.zc[k]
        prog_gm.u[k] = APL.SP_u(FT)(z)
        prog_gm.v[k] = APL.SP_v(FT)(z)
        aux_gm.θ_liq_ice[k] = APL.SP_θ_liq_ice(FT)(z)
        prog_gm.ρθ_liq_ice[k] = ρ0_c[k] * aux_gm.θ_liq_ice[k]
        aux_gm.q_tot[k] = APL.SP_q_tot(FT)(z)
        prog_gm.ρq_tot[k] = ρ0_c[k] * aux_gm.q_tot[k]
        aux_gm.tke[k] = APL.SP_tke(FT)(z)
    end
end

function initialize_forcing(self::CasesBase{SP}, grid::Grid, state, param_set)
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

ForcingBase(case::DryBubble, param_set::APS; kwargs...) =
    ForcingBase(get_forcing_type(case); apply_coriolis = false, apply_subsidence = false)

function surface_ref_state(::DryBubble, param_set::APS, namelist)
    Pg = 1.0e5  #Pressure at ground
    Tg = 296.0
    qtg = 1.0e-5
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{DryBubble}, grid::Grid, param_set, state)

    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0
    # initialize Grid Mean Profiles of thetali and qt
    zc_in = grid.zc
    FT = eltype(grid)
    prof_θ_liq_ice = APL.DryBubble_θ_liq_ice(FT)
    aux_gm.θ_liq_ice .= prof_θ_liq_ice.(zc_in)
    @. prog_gm.ρθ_liq_ice = ρ0_c * aux_gm.θ_liq_ice
    parent(prog_gm.u) .= 0.01
    parent(prog_gm.ρq_tot) .= 0
    parent(aux_gm.q_tot) .= 0
    parent(aux_gm.tke) .= 0
    parent(aux_gm.Hvar) .= 0
    parent(aux_gm.QTvar) .= 0
    parent(aux_gm.HQTcov) .= 0
end

function surface_params(case::DryBubble, zc_surf, surf_ref_state, param_set; Ri_bulk_crit)
    FT = eltype(surf_ref_state)
    Tsurface = 300.0
    qsurface = 0.0
    shf = 0.0001 # only prevent zero division in SF.jl lmo
    lhf = 0.0001 # only prevent zero division in SF.jl lmo
    ustar = 0.1

    kwargs = (; Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.FixedFrictionVelocity; kwargs...)
end

#####
##### LES_driven_SCM
#####

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

function ForcingBase(case::LES_driven_SCM, param_set::APS; nudge_tau)
    ForcingBase(
        get_forcing_type(case);
        apply_coriolis = false,
        apply_subsidence = true,
        coriolis_param = 0.376e-4,
        nudge_tau = nudge_tau,
    )
end

function surface_ref_state(::LES_driven_SCM, param_set::APS, namelist)
    les_filename = namelist["meta"]["lesfile"]

    Pg, Tg, qtg = NC.Dataset(les_filename, "r") do data
        pg_str = haskey(data.group["reference"], "p0_full") ? "p0_full" : "p0"
        Pg = data.group["reference"][pg_str][1] #Pressure at ground
        Tg = data.group["reference"]["temperature0"][1] #Temperature at ground
        ql_ground = data.group["reference"]["ql0"][1]
        qv_ground = data.group["reference"]["qv0"][1]
        qi_ground = data.group["reference"]["qi0"][1]
        qtg = ql_ground + qv_ground + qi_ground #Total water mixing ratio at surface
        (Pg, Tg, qtg)
    end
    return TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
end

function initialize_profiles(self::CasesBase{LES_driven_SCM}, grid::Grid, param_set, state)

    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ0_c = TC.center_ref_state(state).ρ0
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
        parent(aux_gm.θ_liq_ice) .=
            pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "thetali_mean", imin, imax))
        parent(aux_gm.q_tot) .= pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "qt_mean", imin, imax))
        parent(prog_gm.u) .= pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "u_mean", imin, imax))
        parent(prog_gm.v) .= pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "v_mean", imin, imax))
        @. prog_gm.ρθ_liq_ice = ρ0_c * aux_gm.θ_liq_ice
        @. prog_gm.ρq_tot .= ρ0_c * aux_gm.q_tot
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

function surface_params(case::LES_driven_SCM, zc_surf, surf_ref_state, param_set; Ri_bulk_crit, LESDat)
    FT = eltype(surf_ref_state)
    nt = NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zrough = 1.0e-4
        Tsurface = Statistics.mean(data.group["timeseries"]["surface_temperature"][:][imin:imax], dims = 1)[1]
        # get surface value of q
        mean_qt_prof = Statistics.mean(data.group["profiles"]["qt_mean"][:][:, imin:imax], dims = 2)[:]
        # TODO: this will need to be changed if/when we don't prescribe surface fluxes
        qsurface = FT(0)
        lhf = Statistics.mean(data.group["timeseries"]["lhf_surface_mean"][:][imin:imax], dims = 1)[1]
        shf = Statistics.mean(data.group["timeseries"]["shf_surface_mean"][:][imin:imax], dims = 1)[1]
        (; zrough, Tsurface, qsurface, lhf, shf)
    end
    UnPack.@unpack zrough, Tsurface, qsurface, lhf, shf = nt

    ustar = FT(0) # TODO: why is initialization missing?
    kwargs = (; zrough, Tsurface, qsurface, shf, lhf, ustar, Ri_bulk_crit)
    return TC.FixedSurfaceFlux(FT, TC.VariableFrictionVelocity; kwargs...)
end

initialize_forcing(self::CasesBase{LES_driven_SCM}, grid::Grid, state, param_set) =
    initialize(self.Fo, grid, state, self.LESDat)

initialize_radiation(self::CasesBase{LES_driven_SCM}, grid::Grid, state, param_set) =
    initialize(self.Rad, grid, state, self.LESDat)

end # module Cases
