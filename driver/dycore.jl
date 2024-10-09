import UnPack
import LinearAlgebra as LA
import LinearAlgebra: ×

import TurbulenceConvection as TC
import TurbulenceConvection.Parameters as TCP
const APS = TCP.AbstractTurbulenceConvectionParameters
import Thermodynamics as TD
import ClimaCore as CC
import ClimaCore.Geometry as CCG
import OrdinaryDiffEq as ODE

import CLIMAParameters as CP

include(joinpath(@__DIR__, "dycore_variables.jl"))

#####
##### Methods
#####

####
#### Reference state
####

"""
    compute_ref_state!(
        state,
        grid::Grid,
        param_set::PS;
        ts_g,
    ) where {PS}

TODO: add better docs once the API converges

The reference profiles, given
 - `grid` the grid
 - `param_set` the parameter set
 - `ts_g` the surface reference state (a thermodynamic state)
"""
function compute_ref_state!(state, grid::TC.Grid, param_set::PS; ts_g) where {PS}
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    p_f = aux_gm_f.p
    ρ_f = aux_gm_f.ρ
    compute_ref_state!(p_c, ρ_c, p_f, ρ_f, grid, param_set; ts_g)
end

function compute_ref_state!(
    p_c::CC.Fields.Field,
    ρ_c::CC.Fields.Field,
    p_f::CC.Fields.Field,
    ρ_f::CC.Fields.Field,
    grid::TC.Grid,
    param_set::PS;
    ts_g,
) where {PS}
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = TC.float_type(p_c)
    kf_surf = TC.kf_surface(grid)
    qtg = TD.total_specific_humidity(thermo_params, ts_g)
    Φ = TC.geopotential(param_set, grid.zf[kf_surf].z)
    mse_g = TD.moist_static_energy(thermo_params, ts_g, Φ)
    Pg = TD.air_pressure(thermo_params, ts_g)

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    logp = log(Pg)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure
    function minus_inv_scale_height(logp, u, z)
        p_ = exp(logp)
        grav = FT(TCP.grav(param_set))
        Φ = TC.geopotential(param_set, z)
        h = TC.enthalpy(mse_g, Φ)
        ts = TD.PhaseEquil_phq(thermo_params, p_, h, qtg)
        R_m = TD.gas_constant_air(thermo_params, ts)
        T = TD.air_temperature(thermo_params, ts)
        return -FT(TCP.grav(param_set)) / (T * R_m)
    end

    # Perform the integration
    z_span = (grid.zmin, grid.zmax)
    prob = ODE.ODEProblem(minus_inv_scale_height, logp, z_span)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    parent(p_f) .= sol.(vec(grid.zf.z))
    parent(p_c) .= sol.(vec(grid.zc.z))

    p_f .= exp.(p_f)
    p_c .= exp.(p_c)

    # Compute reference state thermodynamic profiles
    @inbounds for k in TC.real_center_indices(grid)
        Φ = TC.geopotential(param_set, grid.zc[k].z)
        h = TC.enthalpy(mse_g, Φ)
        ts = TD.PhaseEquil_phq(thermo_params, p_c[k], h, qtg)
        ρ_c[k] = TD.air_density(thermo_params, ts)
    end

    @inbounds for k in TC.real_face_indices(grid)
        Φ = TC.geopotential(param_set, grid.zf[k].z)
        h = TC.enthalpy(mse_g, Φ)
        ts = TD.PhaseEquil_phq(thermo_params, p_f[k], h, qtg)
        ρ_f[k] = TD.air_density(thermo_params, ts)
    end
    return nothing
end


function set_thermo_state_from_prog!(state, grid, moisture_model, param_set)
    Ic = CCO.InterpolateF2C()
    thermo_params = TCP.thermodynamics_params(param_set)
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ

    @inbounds for k in TC.real_center_indices(grid)
        thermo_args = if moisture_model isa TC.EquilibriumMoisture
            ()
        elseif moisture_model isa TC.NonEquilibriumMoisture
            (prog_gm.q_liq[k], prog_gm.q_ice[k])
        else
            error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
        end
        ts_gm[k] = TC.thermo_state_pθq(
            param_set,
            p_c[k],
            prog_gm.ρθ_liq_ice[k] / ρ_c[k],
            prog_gm.ρq_tot[k] / ρ_c[k],
            thermo_args...,
        )
    end
    return nothing
end

function set_thermo_state_from_aux!(state, grid, moisture_model, param_set)
    Ic = CCO.InterpolateF2C()
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    p_c = aux_gm.p
    @inbounds for k in TC.real_center_indices(grid)
        thermo_args = if moisture_model isa TC.EquilibriumMoisture
            ()
        elseif moisture_model isa TC.NonEquilibriumMoisture
            (prog_gm.q_liq[k], prog_gm.q_ice[k])
        else
            error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
        end
        ts_gm[k] = TC.thermo_state_pθq(param_set, p_c[k], aux_gm.θ_liq_ice[k], aux_gm.q_tot[k], thermo_args...)
    end
    return nothing
end

function set_grid_mean_from_thermo_state!(param_set, state, grid)
    thermo_params = TCP.thermodynamics_params(param_set)
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    ρ_c = prog_gm.ρ
    @. prog_gm.ρθ_liq_ice = ρ_c * TD.liquid_ice_pottemp(thermo_params, ts_gm)
    @. prog_gm.ρq_tot = ρ_c * TD.total_specific_humidity(thermo_params, ts_gm)
    return nothing
end

function assign_thermo_aux!(state, grid, moisture_model, param_set)
    thermo_params = TCP.thermodynamics_params(param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    ρ_c = prog_gm.ρ
    @inbounds for k in TC.real_center_indices(grid)
        ts = ts_gm[k]
        aux_gm.q_tot[k] = prog_gm.ρq_tot[k] / ρ_c[k]
        aux_gm.θ_liq_ice[k] = prog_gm.ρθ_liq_ice[k] / ρ_c[k]
        aux_gm.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
        aux_gm.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)
        aux_gm.T[k] = TD.air_temperature(thermo_params, ts)
        aux_gm.RH[k] = TD.relative_humidity(thermo_params, ts)
    end
    return
end

function ∑stoch_tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    UnPack.@unpack edmf, param_set, surf_params, aux = params

    # set all tendencies to zero
    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0

    CC.Fields.bycolumn(axes(prog.cent)) do colidx

        state = TC.column_state(prog, aux, tendencies, colidx)
        grid = TC.Grid(state)
        surf = get_surface(surf_params, grid, state, t, param_set)

        # compute updraft stochastic tendencies
        TC.compute_up_stoch_tendencies!(edmf, grid, state, param_set, surf)
    end
end

# Compute the sum of tendencies for the scheme
function ∑tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    CC.Fields.bycolumn(axes(prog.cent)) do colidx
        UnPack.@unpack edmf, precip_model, param_set, case = params
        UnPack.@unpack surf_params, radiation, forcing, aux, TS = params

        thermo_params = TCP.thermodynamics_params(param_set)

        parent(tendencies.face[colidx]) .= 0
        parent(tendencies.cent[colidx]) .= 0
        state = TC.column_state(prog, aux, tendencies, colidx)
        grid = TC.Grid(state)

        set_thermo_state_from_prog!(state, grid, edmf.moisture_model, param_set)
        assign_thermo_aux!(state, grid, edmf.moisture_model, param_set)

        aux_gm = TC.center_aux_grid_mean(state)

        @. aux_gm.θ_virt = TD.virtual_pottemp(thermo_params, aux_gm.ts)

        # reduce dt during spinup period for stability... (edit in side TS just in case it's called anywhere else (don't think it it but can't be too careful... original setting is in callbacks if using adaptive, but we need to be able to set it generally otherwise...)
        Δt = TS.dt # we need to extract Δt first and then edit it... we cannot directly edit TS.dt bc repeated calls would decimate TS.dt... hopefully this causes no dissonance w/ callbacks setting in adapt_dt, or just dt=dt_min in regular case...
        if !TS.spinup_adapt_dt
            Δt = TS.dt_min # enforce no adapt_dt if not allowed
        end
        if t < TS.spinup_half_t_max
            Δt *= TS.spinup_dt_factor # can we fade this out gradually between t=spinup_half_t_max and t=2*spinup_half_t_max?, and ramp up from Δt*spinup_dt_factor to Δt?
        elseif t < 2 * TS.spinup_half_t_max # if spinup_half_t_max is 0, we'll bypass this
            Δt *= TS.spinup_dt_factor + (1 - TS.spinup_dt_factor) * (t - TS.spinup_half_t_max) / TS.spinup_half_t_max # a smooth ramp from TS.spinup_dt_factor to 1 between t=spinup_half_t_max and t=2*spinup_half_t_max
        end
  
        surf = get_surface(surf_params, grid, state, t, param_set)

        TC.affect_filter!(edmf, grid, state, param_set, surf, t)

        # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
        # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
        # treat them as auxiliary variables for now, until we disentangle the tendency computations.
        if isa(case, Cases.SOCRATES)
            Cases.update_forcing(case, grid, state, t, param_set, forcing) # should dispatch to socrates
        else
            Cases.update_forcing(case, grid, state, t, param_set)
        end
        Cases.update_radiation(radiation, grid, state, t, param_set)

        TC.update_aux!(edmf, grid, state, surf, param_set, t, Δt)

        # compute tendencies
        # causes division error in dry bubble first time step
        TC.compute_precipitation_sink_tendencies(precip_model, edmf, grid, state, param_set, Δt)
        TC.compute_precipitation_advection_tendencies(precip_model, edmf, grid, state, param_set)

        TC.compute_turbconv_tendencies!(edmf, grid, state, param_set, surf, Δt)
        compute_gm_tendencies!(edmf, grid, state, surf, radiation, forcing, param_set)
    end

    return nothing
end

function compute_les_Γᵣ(z::FT, τᵣ::FT = 24.0 * 3600.0, zᵢ::FT = 3000.0, zᵣ::FT = 3500.0) where {FT <: Real}
    # returns height-dependent relaxation timescale from eqn. 9 in `Shen et al. 2021`
    if z < zᵢ
        return FT(0)
    elseif zᵢ <= z <= zᵣ
        cos_arg = pi * ((z - zᵢ) / (zᵣ - zᵢ))
        return (FT(0.5) / τᵣ) * (1 - cos(cos_arg))
    elseif z > zᵣ
        return (1 / τᵣ)
    end
end

function compute_gm_tendencies!(
    edmf::TC.EDMFModel,
    grid::TC.Grid,
    state::TC.State,
    surf::TC.SurfaceBase,
    radiation::Cases.RadiationBase,
    force::Cases.ForcingBase,
    param_set::APS,
)
    thermo_params = TCP.thermodynamics_params(param_set)
    Ic = CCO.InterpolateF2C()
    tendencies_gm = TC.center_tendencies_grid_mean(state)
    kc_toa = TC.kc_top_of_atmos(grid)
    kf_surf = TC.kf_surface(grid)
    FT = TC.float_type(state)
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    ∇θ_liq_ice_gm = TC.center_aux_grid_mean(state).∇θ_liq_ice_gm
    ∇q_tot_gm = TC.center_aux_grid_mean(state).∇q_tot_gm
    aux_en = TC.center_aux_environment(state)
    aux_en_f = TC.face_aux_environment(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_bulk = TC.center_aux_bulk(state)
    ρ_f = aux_gm_f.ρ
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_tc = TC.center_aux_turbconv(state)
    ts_gm = TC.center_aux_grid_mean(state).ts

    θ_liq_ice_gm_toa = prog_gm.ρθ_liq_ice[kc_toa] / ρ_c[kc_toa]
    q_tot_gm_toa = prog_gm.ρq_tot[kc_toa] / ρ_c[kc_toa]
    # RBθ = CCO.RightBiasedC2F(; top = CCO.SetValue(θ_liq_ice_gm_toa)) # right biased
    # RBq = CCO.RightBiasedC2F(; top = CCO.SetValue(q_tot_gm_toa)) # right biased
    wvec = CC.Geometry.WVector
    ∇c = CCO.DivergenceF2C() # F2C to come back from C2F
    # @. ∇θ_liq_ice_gm = ∇c(wvec(RBθ(prog_gm.ρθ_liq_ice / ρ_c)))
    # @. ∇q_tot_gm = ∇c(wvec(RBq(prog_gm.ρq_tot / ρ_c)))

    # if edmf.moisture_model isa TC.NonEquilibriumMoisture
    #     ∇q_liq_gm = TC.center_aux_grid_mean(state).∇q_liq_gm
    #     ∇q_ice_gm = TC.center_aux_grid_mean(state).∇q_ice_gm
    #     q_liq_gm_toa = prog_gm.q_liq[kc_toa]
    #     q_ice_gm_toa = prog_gm.q_ice[kc_toa]
    #     RBq_liq = CCO.RightBiasedC2F(; top = CCO.SetValue(q_liq_gm_toa))
    #     RBq_ice = CCO.RightBiasedC2F(; top = CCO.SetValue(q_ice_gm_toa))
    #     @. ∇q_liq_gm = ∇c(wvec(RBq_liq(prog_gm.q_liq)))
    #     @. ∇q_ice_gm = ∇c(wvec(RBq_ice(prog_gm.q_ice)))
    # end

    # Apply forcing and radiation
    prog_gm_uₕ = TC.grid_mean_uₕ(state)
    aux_gm_uₕ_g = TC.grid_mean_uₕ_g(state)
    tendencies_gm_uₕ = TC.tendencies_grid_mean_uₕ(state)
    prog_gm_u = TC.physical_grid_mean_u(state)
    prog_gm_v = TC.physical_grid_mean_v(state)

    # Coriolis
    coriolis_param = force.coriolis_param
    # TODO: union split over sphere or box
    # lat = CC.Fields.coordinate_field(axes(ρ_c)).lat
    coords = CC.Fields.coordinate_field(axes(ρ_c))
    coriolis_fn(coord) = CCG.WVector(coriolis_param)
    f = @. CCG.Contravariant3Vector(coriolis_fn(coords))

    C123 = CCG.Covariant123Vector
    C12 = CCG.Contravariant12Vector
    lg = CC.Fields.local_geometry_field(axes(ρ_c))
    @. tendencies_gm_uₕ -= f × (C12(C123(prog_gm_uₕ)) - C12(C123(aux_gm_uₕ_g)))


    # ========================================================================================================================================== #
    # trying subsidence here outside the loop w/ upwinding...
    # 1) it seems we before just used a straight product, subsidence * dq/dz
    # 2) This could be bad because q is right biased C2F (towards ground), so dq/dz F2C encodes information shifted towards the ground (i.e. from the q above).
    #    Then, if we have descent, w x dqdz on the center is w times q_above-q and the flux is good from upstream :)
    #    But ascent is propagating q-q_above into the current node, when the information should be q-q_below. 
    #    This is bad, if the point above is smaller, for example, this point will continue to grow -- even if the point below is smaller and w is upwards, as it has no way to access that information -- and blowup :(
    # 3) here we'll try upwinding d/dz (wq) instead. This hopefully will make ascent stable :)
    # 4) unclear if things have to be rightbiased in converting c to f, instead of just interpolatec2f (with what BCs?), so we'll see if it's stable...
    # 5) assume no penetration BCs on subsidence w

    UBsub = CCO.UpwindBiasedProductC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # upwinding, extrapolate bc we don't know the boa/toa derivatives

    kc_surf = TC.kc_surface(grid)
    θ_liq_ice_gm_boa = prog_gm.ρθ_liq_ice[kc_surf] / ρ_c[kc_surf]
    q_tot_gm_boa = prog_gm.ρq_tot[kc_surf] / ρ_c[kc_surf]

    # subsidence should be face valued so we'll need a C2F call
    C2Fsub = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0))) # not sure but maybe we should use this for stability as it's not biased for your upwind algorithm... no penetration. # I think this hsould be the more stable option? Though according to cgarkue github issue, it doesn't always work (see https://github.com/CliMA/TurbulenceConvection.jl/pull/1146)
    C2Fθ = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(θ_liq_ice_gm_boa)), top = CCO.SetValue(FT(θ_liq_ice_gm_toa))) # not sure if this should be right biased or not
    C2Fq = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(q_tot_gm_boa)), top = CCO.SetValue(FT(q_tot_gm_toa))) # not sure if this should be right biased or not

    # we to C2F, then ∇c goes F2C, then UBsub got C2F, but we wanna end on C [ there's no UpwindBiasedProductF2C and we had UpwindBiasedProductF2C(u[F], x[C])
    F2Csub = CCO.InterpolateF2C(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # We might have put 0 for no penetration on boundary, but that's not exactly true in our dataset...
    CV32FT = x -> x[1] # convert Contravariant3Vector to Float64

    # regular upwinding 
    ∇θ_liq_ice_gm = TC.center_aux_grid_mean(state).∇θ_liq_ice_gm
    @. ∇θ_liq_ice_gm = ∇c(wvec(C2Fθ(prog_gm.ρθ_liq_ice / ρ_c)))
    w∇θ_liq_ice_gm = @. UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇θ_liq_ice_gm) # works, but output type is different than tendencies type
    w∇θ_liq_ice_gm = @. F2Csub(w∇θ_liq_ice_gm) # convert back to C
    w∇θ_liq_ice_gm = @. CV32FT(w∇θ_liq_ice_gm) # convert back to Float64 from Contravariant3Vector

    ∇q_tot_gm = TC.center_aux_grid_mean(state).∇q_tot_gm
    @. ∇q_tot_gm = ∇c(wvec(C2Fq(prog_gm.ρq_tot / ρ_c)))
    w∇q_tot_gm = @. UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_tot_gm) # works, but output type is different than tendencies type
    w∇q_tot_gm = @. F2Csub(w∇q_tot_gm) # convert back to C
    w∇q_tot_gm = @. CV32FT(w∇q_tot_gm) # convert back to Float64 from Contravariant3Vector


    if edmf.moisture_model isa TC.NonEquilibriumMoisture
        q_liq_gm_toa = prog_gm.q_liq[kc_toa]
        q_ice_gm_toa = prog_gm.q_ice[kc_toa]
        q_liq_boa = prog_gm.q_liq[kc_surf]
        q_ice_boa = prog_gm.q_ice[kc_surf]
        C2Fq_liq = CCO.InterpolateC2F(bottom = CCO.SetValue(FT(q_liq_boa)), top = CCO.SetValue(FT(q_liq_gm_toa))) # not sure if this should be right biased or not
        C2Fq_ice = CCO.InterpolateC2F(bottom = CCO.SetValue(FT(q_ice_boa)), top = CCO.SetValue(FT(q_ice_gm_toa))) # not sure if this should be right biased or not

        ∇q_liq_gm = TC.center_aux_grid_mean(state).∇q_liq_gm
        @. ∇q_liq_gm = ∇c(wvec(C2Fq_liq(prog_gm.q_liq)))
        w∇q_liq_gm = @. UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_liq_gm) # works, but output type is different than tendencies type
        w∇q_liq_gm = @. F2Csub(w∇q_liq_gm) # convert back to C
        w∇q_liq_gm = @. CV32FT(w∇q_liq_gm) # convert back to Float64 from Contravariant3Vector

        ∇q_ice_gm = TC.center_aux_grid_mean(state).∇q_ice_gm
        @. ∇q_ice_gm = ∇c(wvec(C2Fq_ice(prog_gm.q_ice)))
        w∇q_ice_gm = @. UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_ice_gm) # works, but output type is different than tendencies type
        w∇q_ice_gm = @. F2Csub(w∇q_ice_gm) # convert back to C
        w∇q_ice_gm = @. CV32FT(w∇q_ice_gm) # convert back to Float64 from Contravariant3Vector
    end


    @inbounds for k in TC.real_center_indices(grid)
        Π = TD.exner(thermo_params, ts_gm[k])

        # LS Subsidence (any positive subsidence (ascent) leads to numerical instability) -- why? fast descent is unstable too...)
        # tendencies_gm.ρθ_liq_ice[k] -= ρ_c[k] * aux_gm.subsidence[k] * ∇θ_liq_ice_gm[k] # trying turning this off
        tendencies_gm.ρθ_liq_ice[k] -= ρ_c[k] * w∇θ_liq_ice_gm[k] # trying turning this off
        # tendencies_gm.ρq_tot[k] -= ρ_c[k] * aux_gm.subsidence[k] * ∇q_tot_gm[k] # trying turning this off
        tendencies_gm.ρq_tot[k] -= ρ_c[k] * w∇q_tot_gm[k] # trying turning this off
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            # tendencies_gm.q_liq[k] -= ∇q_liq_gm[k] * aux_gm.subsidence[k]
            tendencies_gm.q_liq[k] -= w∇q_liq_gm[k]
            # tendencies_gm.q_ice[k] -= ∇q_ice_gm[k] * aux_gm.subsidence[k]
            tendencies_gm.q_ice[k] -= w∇q_ice_gm[k]
        end

        # Radiation
        if Cases.rad_type(radiation) <:
           Union{Cases.RadiationDYCOMS_RF01, Cases.RadiationLES, Cases.RadiationTRMM_LBA, Cases.RadiationSOCRATES}
            tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * aux_gm.dTdt_rad[k] / Π .* 1  # testing maybe it's an order of operations thing and the radiation needs to happen first before precip then nudging? idk
        end
        # LS advection
        tendencies_gm.ρq_tot[k] += ρ_c[k] * aux_gm.dqtdt_hadv[k] .* 1
        if !(Cases.force_type(force) <: Cases.ForcingDYCOMS_RF01)
            tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * aux_gm.dTdt_hadv[k] / Π .* 1 # testing maybe it's an order of operations thing and the radiation needs to happen first before precip then nudging? idk
        end
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tendencies_gm.q_liq[k] += aux_gm.dqldt[k] # we never set these, e.g. in socrates... presumably they're zero?
            tendencies_gm.q_ice[k] += aux_gm.dqidt[k]
        end

        # LES specific forcings
        if Cases.force_type(force) <: Cases.ForcingLES
            H_fluc = aux_gm.dTdt_fluc[k] / Π

            gm_U_nudge_k = (aux_gm.u_nudge[k] - prog_gm_u[k]) / force.wind_nudge_τᵣ
            gm_V_nudge_k = (aux_gm.v_nudge[k] - prog_gm_v[k]) / force.wind_nudge_τᵣ

            Γᵣ = compute_les_Γᵣ(grid.zc[k].z, force.scalar_nudge_τᵣ, force.scalar_nudge_zᵢ, force.scalar_nudge_zᵣ)
            gm_H_nudge_k = Γᵣ * (aux_gm.H_nudge[k] - aux_gm.θ_liq_ice[k])
            gm_q_tot_nudge_k = Γᵣ * (aux_gm.qt_nudge[k] - aux_gm.q_tot[k])
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                gm_q_liq_nudge_k = Γᵣ * (aux_gm.ql_nudge[k] - prog_gm.q_liq[k])
                gm_q_ice_nudge_k = Γᵣ * (aux_gm.qi_nudge[k] - prog_gm.q_ice[k])
            end

            tendencies_gm.ρq_tot[k] += ρ_c[k] * (gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k])
            tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * (gm_H_nudge_k + H_fluc)
            tendencies_gm_uₕ[k] += CCG.Covariant12Vector(CCG.UVVector(gm_U_nudge_k, gm_V_nudge_k), lg[k])
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                tendencies_gm.q_liq[k] += aux_gm.dqldt_hadv[k] + gm_q_liq_nudge_k + aux_gm.dqldt_fluc[k]
                tendencies_gm.q_ice[k] += aux_gm.dqidt_hadv[k] + gm_q_ice_nudge_k + aux_gm.dqidt_fluc[k]
            end
        end

        if Cases.force_type(force) <: Cases.ForcingSOCRATES
            # temperature horizontal divergence is handeled above
            # geostrophic handled elsewhere (only affects u,v tendency and maybe surface fluxes?)

            # nudge back to era 5 horizosntal
            gm_U_nudge_k = (aux_gm.u_nudge[k] - prog_gm_u[k]) / force.wind_nudge_τᵣ
            gm_V_nudge_k = (aux_gm.v_nudge[k] - prog_gm_v[k]) / force.wind_nudge_τᵣ
            tendencies_gm_uₕ[k] += CCG.Covariant12Vector(CCG.UVVector(gm_U_nudge_k, gm_V_nudge_k), lg[k]) # add nudge (again we have dqtdt_hadv so what re dqldt and dqidt derived from?)

            # nudge liq-ice pottemp back towards our forcing profile
            H_fluc = aux_gm.dTdt_fluc[k] / Π # should be 0
            gm_H_nudge_k = (aux_gm.H_nudge[k] - aux_gm.θ_liq_ice[k]) / force.scalar_nudge_τᵣ # go back to regular tau relaxation (entire column now unlike zhaoyi's func )
            tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * (gm_H_nudge_k + H_fluc) # we have no fluc so get that out (but we have the other terms so brought it back maybe itll lfix intsbailitis?)

            # nudge total water back towards our forcing profile
            gm_q_tot_nudge_k = (aux_gm.qt_nudge[k] - aux_gm.q_tot[k]) / force.scalar_nudge_τᵣ
            tendencies_gm.ρq_tot[k] += ρ_c[k] * (gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k])

            # nudge liquid and ice -- but we don't have a ql_nudge qi_nudge from external forcing in socrates? only qt_nudge -- what are ql_nudge and qi_nudge derived from
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                # gm_q_liq_nudge_k = (aux_gm.ql_nudge[k] - prog_gm.q_liq[k]) / force.scalar_nudge_τᵣ
                # gm_q_ice_nudge_k = (aux_gm.qi_nudge[k] - prog_gm.q_ice[k]) / force.scalar_nudge_τᵣ
                # tendencies_gm.q_liq[k] += aux_gm.dqldt_hadv[k] + gm_q_liq_nudge_k + aux_gm.dqldt_fluc[k] # we dont have ql_nudge or dqldt_hadv,fluc in socrates
                # tendencies_gm.q_ice[k] += aux_gm.dqidt_hadv[k] + gm_q_ice_nudge_k + aux_gm.dqidt_fluc[k] # we dont have qi_nudge or dqidt_hadv,fluc in socrates
            end
        end

        # Apply precipitation tendencies
        tendencies_gm.ρq_tot[k] +=
            ρ_c[k] * (
                aux_bulk.qt_tendency_precip_formation[k] +
                aux_en.qt_tendency_precip_formation[k] +
                aux_tc.qt_tendency_precip_sinks[k]
            )

        tendencies_gm.ρθ_liq_ice[k] +=
            ρ_c[k] * (
                aux_bulk.θ_liq_ice_tendency_precip_formation[k] +
                aux_en.θ_liq_ice_tendency_precip_formation[k] +
                aux_tc.θ_liq_ice_tendency_precip_sinks[k]
            )

        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tendencies_gm.q_liq[k] += aux_bulk.ql_tendency_precip_formation[k] + aux_en.ql_tendency_precip_formation[k]
            tendencies_gm.q_ice[k] += aux_bulk.qi_tendency_precip_formation[k] + aux_en.qi_tendency_precip_formation[k]

            # Additionally apply cloud liquid and ice formation tendencies ( no effect on qtot, θ_li)
            tendencies_gm.q_liq[k] += aux_bulk.ql_tendency_noneq[k] + aux_en.ql_tendency_noneq[k]
            tendencies_gm.q_ice[k] += aux_bulk.qi_tendency_noneq[k] + aux_en.qi_tendency_noneq[k]
        end
    end


    # ====== sedimentation grid mean only (stability fix...)
    if TC.get_isbits_nt(param_set.user_args, :use_sedimentation, false)
        if TC.get_isbits_nt(param_set.user_args, :grid_mean_sedimentation, false)
            # @info "grid mean sedimenting"
            ts_gm = TC.center_aux_grid_mean(state).ts
            mph, _ = TC.calculate_sedimentation_sources(param_set, grid, ρ_c, ts_gm; grid_mean=true)

            L_v0 = TCP.LH_v0(param_set)
            L_s0 = TCP.LH_s0(param_set)
            @inbounds for k in TC.real_center_indices(grid)
                ql_sedimentation_tendency = mph.ql_tendency[k]
                qi_sedimentation_tendency = mph.qi_tendency[k]
                qt_sedimentation_tendency = ql_sedimentation_tendency + qi_sedimentation_tendency

                # How does updrafts/bulk know about these? does it just all end up in the environment?
                # We need to split somehow so that the updrafts/bulk can also know about these tendencies...
                # How do the large-scale tendencies work? do they just work on the environment?
                # if so, that leads to instabilities here...

                if edmf.moisture_model isa TC.NonEquilibriumMoisture
                    tendencies_gm.q_liq[k] += ql_sedimentation_tendency
                    tendencies_gm.q_ice[k] += qi_sedimentation_tendency
                end
                tendencies_gm.ρq_tot[k] += ρ_c[k] * qt_sedimentation_tendency

                Π_m = TD.exner(thermo_params, ts_gm[k])
                c_pm = TD.cp_m(thermo_params, ts_gm[k])
                θ_liq_ice_sedimentation_tendency = 1 / Π_m / c_pm * ( L_v0 * ql_sedimentation_tendency + L_s0 * qi_sedimentation_tendency )
                tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * θ_liq_ice_sedimentation_tendency
            end
        else # TC.get_isbits_nt(param_set.user_args, :grid_mean_sedimentation, false) # separate env/updraft sedimentation calcs

            # kmax = argmax([abs(aux_en.qt_tendency_sedimentation[k]) for k in TC.real_center_indices(grid)])
            # kmax = TC.real_center_indices(grid)[kmax] # convert to real index
            # @info "z", grid.zc[kmax].z,
            # " ql_tendency_sedimentation[kmax]: ", (aux_en.ql_tendency_sedimentation[kmax], aux_bulk.ql_tendency_sedimentation[kmax]),
            # " qi_tendency_sedimentation[kmax]: ", (aux_en.qi_tendency_sedimentation[kmax], aux_bulk.qi_tendency_sedimentation[kmax]),
            # " qt_tendency_sedimentation[kmax]: ", (aux_en.qt_tendency_sedimentation[kmax], aux_bulk.qt_tendency_sedimentation[kmax])

            @inbounds for k in TC.real_center_indices(grid)
                # @info aux_en.ql_tendency_sedimentation[k] + aux_bulk.ql_tendency_sedimentation[k]
                # @info aux_en.qi_tendency_sedimentation[k] + aux_bulk.qi_tendency_sedimentation[k]
                # @info aux_en.qt_tendency_sedimentation[k] + aux_bulk.qt_tendency_sedimentation[k]

                if edmf.moisture_model isa TC.NonEquilibriumMoisture
                    tendencies_gm.q_liq[k] += aux_en.ql_tendency_sedimentation[k] + aux_bulk.ql_tendency_sedimentation[k]
                    tendencies_gm.q_ice[k] += aux_en.qi_tendency_sedimentation[k] + aux_bulk.qi_tendency_sedimentation[k]
                end

                tendencies_gm.ρq_tot[k] +=
                ρ_c[k] * (
                    aux_bulk.qt_tendency_sedimentation[k] +
                    aux_en.qt_tendency_sedimentation[k]
                )

                tendencies_gm.ρθ_liq_ice[k] +=
                ρ_c[k] * (
                    aux_bulk.θ_liq_ice_tendency_sedimentation[k] +
                    aux_en.θ_liq_ice_tendency_sedimentation[k]
                )
            end
        end
    end
    # =================================================== #


    TC.compute_sgs_flux!(edmf, grid, state, surf)
    sgs_flux_θ_liq_ice = aux_gm_f.sgs_flux_θ_liq_ice
    sgs_flux_q_tot = aux_gm_f.sgs_flux_q_tot
    sgs_flux_uₕ = aux_gm_f.sgs_flux_uₕ
    tends_ρθ_liq_ice = tendencies_gm.ρθ_liq_ice
    tends_ρq_tot = tendencies_gm.ρq_tot
    tends_uₕ = TC.tendencies_grid_mean_uₕ(state)

    ∇sgs = CCO.DivergenceF2C()
    @. tends_ρθ_liq_ice += -∇sgs(wvec(sgs_flux_θ_liq_ice))
    @. tends_ρq_tot += -∇sgs(wvec(sgs_flux_q_tot))
    @. tends_uₕ += -∇sgs(sgs_flux_uₕ) / ρ_c

    if edmf.moisture_model isa TC.NonEquilibriumMoisture
        sgs_flux_q_liq = aux_gm_f.sgs_flux_q_liq
        sgs_flux_q_ice = aux_gm_f.sgs_flux_q_ice

        tends_q_liq = tendencies_gm.q_liq
        tends_q_ice = tendencies_gm.q_ice
        @. tends_q_liq += -∇sgs(wvec(sgs_flux_q_liq)) / ρ_c
        @. tends_q_ice += -∇sgs(wvec(sgs_flux_q_ice)) / ρ_c
    end

    return nothing
end
