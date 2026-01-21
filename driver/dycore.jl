import UnPack
import LinearAlgebra as LA
import LinearAlgebra: ×

import TurbulenceConvection as TC
import TurbulenceConvection.Parameters as TCP
const APS = TCP.AbstractTurbulenceConvectionParameters
import Thermodynamics as TD
import ClimaCore as CC
import ClimaCore.Geometry as CCG
# import OrdinaryDiffEq as ODE
import OrdinaryDiffEqTsit5

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
        param_set::PS;
        ts_g,
    ) where {PS}

TODO: add better docs once the API converges

The reference profiles, given
 - `grid` the grid
 - `param_set` the parameter set
 - `ts_g` the surface reference state (a thermodynamic state)
"""
function compute_ref_state!(state, param_set::PS; ts_g) where {PS}
    grid = TC.Grid(state)
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
    prob = SciMLBase.ODEProblem{false, SciMLBase.FullSpecialize}(minus_inv_scale_height, logp, z_span) # false means not in place
    sol = SciMLBase.solve(prob, OrdinaryDiffEqTsit5.Tsit5(), reltol = 1e-12, abstol = 1e-12)
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


function set_thermo_state_from_prog!(state::TC.State, moisture_model::TC.AbstractMoistureModel, param_set::APS)
    grid = TC.Grid(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
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
        ts_gm[k] = TC.thermo_state_pθq(param_set, p_c[k], prog_gm.ρθ_liq_ice[k] / ρ_c[k], prog_gm.ρq_tot[k] / ρ_c[k], thermo_args..., )
    end
    return nothing
end

"""
This function is only used for initializatoin
"""
function set_thermo_state_from_aux!(state, moisture_model, param_set)
    grid = TC.Grid(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    prog_gm = TC.center_prog_grid_mean(state)
    aux_gm = TC.center_aux_grid_mean(state)
    p_c = aux_gm.p
    @inbounds for k in TC.real_center_indices(grid)
        thermo_args = if moisture_model isa TC.EquilibriumMoisture
            ()
        elseif moisture_model isa TC.NonEquilibriumMoisture
            # (prog_gm.q_liq[k], prog_gm.q_ice[k])
            prog_gm.q_liq[k] = aux_gm.q_liq[k] # make sure they agree
            prog_gm.q_ice[k] = aux_gm.q_ice[k] # make sure they agree
            (aux_gm.q_liq[k], aux_gm.q_ice[k])
        else
            error("Something went wrong. The moisture_model options are equilibrium or nonequilibrium")
        end
        ts_gm[k] = TC.thermo_state_pθq(param_set, p_c[k], aux_gm.θ_liq_ice[k], aux_gm.q_tot[k], thermo_args...)
    end
    return nothing
end

function set_grid_mean_from_thermo_state!(param_set, state)
    grid = TC.Grid(state)
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

function assign_thermo_aux!(state, param_set)
    grid = TC.Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    ρ_c = prog_gm.ρ
    p_c = aux_gm.p
    @inbounds for k in TC.real_center_indices(grid)
        ts = ts_gm[k]
        aux_gm.q_tot[k] = prog_gm.ρq_tot[k] / ρ_c[k]
        aux_gm.θ_liq_ice[k] = prog_gm.ρθ_liq_ice[k] / ρ_c[k]
        aux_gm.q_liq[k] = TD.liquid_specific_humidity(thermo_params, ts)
        aux_gm.q_ice[k] = TD.ice_specific_humidity(thermo_params, ts)
        aux_gm.T[k] = TD.air_temperature(thermo_params, ts)
        
        # @assert [for search optimization, but don't use assert, use an explicit error]
        # (0 < aux_gm.T[k] < Inf) || error("Negative/zero or nonfinite temperature will cause issues here, status is z = $(grid.zc[k].z), T = $(aux_gm.T[k]), q = $(ts_gm[k].q), θ_liq_ice = $(aux_gm.θ_liq_ice[k]), ts = $ts") # debugging, remove later -- will slow down code

        # if !(0 < aux_gm.T[k] < Inf)
            # @info "status: T = $(TC.full_print(aux_gm.T)); θ_liq_ice = $(TC.full_print(aux_gm.θ_liq_ice)); q_tot = $(TC.full_print(aux_gm.q_tot)); q_liq = $(TC.full_print(aux_gm.q_liq)); q_ice = $(TC.full_print(aux_gm.q_ice))"
            # aux_up = TC.center_aux_updrafts(state)
            # aux_up_f = TC.face_aux_updrafts(state)
            # aux_en = TC.center_aux_environment(state)
            # aux_en_f = TC.face_aux_environment(state)
            # prog_up = TC.center_prog_updrafts(state)
            # prog_up_f = TC.face_prog_updrafts(state)
            # prog_gm = TC.center_prog_grid_mean(state)
            # aux_gm = TC.center_aux_grid_mean(state)

            # println("---------------------------------------------------")
            # # summary(stdout, state.prog)
            # @info summary(state.prog)
            # println("---------------------------------------------------")
            # # summary(stdout, state.aux)
            # @info summary(state.aux)
            # println("---------------------------------------------------")
            # flush(stdout); flush(stderr)
            # error("Negative or nonfinite temperature will cause issues here, status is z = $(grid.zc[k].z), T = $(aux_gm.T[k]), q = $(ts_gm[k].q), θ_liq_ice = $(aux_gm.θ_liq_ice[k]), ts = $ts")
        # end
        
        aux_gm.RH[k] = TD.relative_humidity(thermo_params, ts)
        if !state.calibrate_io
            aux_gm.RH_liq[k] = TC.relative_humidity_over_liquid(thermo_params, ts)
            aux_gm.RH_ice[k] = TC.relative_humidity_over_ice(thermo_params, ts)
        end
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

        state = TC.column_state(prog, aux, tendencies, colidx, params.calibrate_io)
        grid = TC.Grid(state)
        surf = get_surface(surf_params, state, t, param_set)

        # compute updraft stochastic tendencies
        TC.compute_up_stoch_tendencies!(edmf, state)
    end
end

"""
"""
function my_unstable_check_test(dt::Real, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    TS = params.TS
    TS.isoutofdomain = false # reset to false at the beginning of the function

    @debug "Checking for nan in prognostic variables"
    
    CC.Fields.bycolumn(axes(prog.cent)) do colidx # have to do this bc it's a do block which is technically an anonymous function...
        UnPack.@unpack edmf, precip_model, param_set, case = params
        UnPack.@unpack surf_params, radiation, forcing, aux, TS = params

        state = TC.column_prog_aux(prog, aux, colidx, params.calibrate_io)
        FT = TC.float_type(state)
        # check for a nan in every variable
        error("not implemented yet")
    end
end

"""
    ``If those don't work, call out the big guns. One of them is isoutofdomain, where you can define a boolean function which will cause step rejections whenever it is not satisfied. For example, isoutofdomain = (u,p,t)->any(x->x<0,u) will make the solver reject any step which cases any variable u to go negative. Now, using any pure-Julia solver with this option, it's impossible to get a negative in the result.``
    Originally, I wanted to use du in a callback to manually find the largest dt
    
"""
neg_or_nan(x::FT) where {FT} = x < FT(0) || isnan(x)
function isoutofdomain(prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
# need to figure something out
    # bycolumn() returns nothing, and also we'd like to short circuit if we find a negative value
    # basically it's passing do_block(colidx) = /block/ as first argument to bycolumn

    
    # not sure if we need this before the by-column call (or does it get automatically rolled back?)
    TS = params.TS
    TS.isoutofdomain = false # reset to false at the beginning of the function

    @debug "Checking for negative values in prognostic variables in isoutofdomain() at t = $t, TS.isoutofdomain = $(params.TS.isoutofdomain)"
    
    CC.Fields.bycolumn(axes(prog.cent)) do colidx # have to do this bc it's a do block which is technically an anonymous function...
        UnPack.@unpack edmf, precip_model, param_set, case = params
        UnPack.@unpack surf_params, radiation, forcing, aux, TS = params

        state = TC.column_prog_aux(prog, aux, colidx, params.calibrate_io)

        FT = TC.float_type(state)

        grid = TC.Grid(state)
        N_up = TC.n_updrafts(edmf)
        kc_surf = TC.kc_surface(grid)
        kf_surf = TC.kf_surface(grid)
    
        aux_tc = TC.center_aux_turbconv(state)
        aux_up = TC.center_aux_updrafts(state)
        aux_en = TC.center_aux_environment(state)
        aux_en_f = TC.face_aux_environment(state)
        aux_up_f = TC.face_aux_updrafts(state)
        aux_bulk = TC.center_aux_bulk(state)
        aux_bulk_f = TC.face_aux_bulk(state) # same as using aux_tc_f = face_aux_turbconv(state) and then using aux_tc_f.bulk.w for example
        prog_up = TC.center_prog_updrafts(state)
        prog_up_f = TC.face_prog_updrafts(state)
        prog_gm = TC.center_prog_grid_mean(state)
        aux_gm = TC.center_aux_grid_mean(state)
        aux_gm_f = TC.face_aux_grid_mean(state)

        ρ_c = prog_gm.ρ
        ρ_f = aux_gm_f.ρ

        # check for negative values in updraft prognostic variables
        @inbounds for i in 1:N_up
            area = aux_up[i].area
            θ_liq_ice = aux_up[i].θ_liq_ice
            q_tot = aux_up[i].q_tot
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                q_liq = aux_up[i].q_liq
                q_ice = aux_up[i].q_ice
            end
            ρaw = prog_up_f[i].ρaw

            @inbounds for k in TC.real_center_indices(grid)
                if (area[k] < 0) || (θ_liq_ice[k] < 0) || (q_tot[k] < 0)
                    @debug "Negative value in updraft area, θ_liq_ice, or q_tot ($(area[k]), $(θ_liq_ice[k]), $(q_tot[k]))"
                    TS.isoutofdomain = true
                    return true
                end
                if edmf.moisture_model isa TC.NonEquilibriumMoisture
                    if (q_liq[k] < 0) || (q_ice[k] < 0)
                        @debug "Negative value in updraft q_liq or q_ice ($(q_liq[k]), $(q_ice[k]))"
                        TS.isoutofdomain = true
                        return true
                    end
                end
            end
            @inbounds for k in TC.real_face_indices(grid)
                if ρaw[k] < 0 # arguably this should be allowed but lol
                    @debug "Negative value in updraft ρaw ($(ρaw[k]))"
                    TS.isoutofdomain = true
                    return true
                end
            end
        end

        # check for negative values in grid mean prognostic variables
        @inbounds for k in TC.real_center_indices(grid)
            if (prog_gm.ρθ_liq_ice[k] < 0) || (prog_gm.ρq_tot[k] < 0)
                @debug "Negative value in grid mean ρθ_liq_ice or ρq_tot ($(prog_gm.ρθ_liq_ice[k]), $(prog_gm.ρq_tot[k]))"
                TS.isoutofdomain = true
                return true
            end
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                if (prog_gm.q_liq[k] < 0) || (prog_gm.q_ice[k] < 0)
                    @debug "Negative value in grid mean q_liq or q_ice ($(prog_gm.q_liq[k]), $(prog_gm.q_ice[k]))"
                    TS.isoutofdomain = true
                    return true
                end
            end
        end

        # check for negative values in grid mean auxiliary variables (some of these are set before tendencies are calculated... see set_thermo_state_from_prog!() etc)
        # I'm not sure how these stay bad when everything else stays good? maybe it's just an unrealistic state idk...
        @inbounds for k in TC.real_center_indices(grid)
            if (aux_gm.θ_liq_ice[k] < 0) || (aux_gm.q_tot[k] < 0)
                @debug "Negative value in grid mean θ_liq_ice or q_tot ($(aux_gm.θ_liq_ice[k]), $(aux_gm.q_tot[k]))"
                TS.isoutofdomain = true
                return true
            end
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                if (aux_gm.q_liq[k] < 0) || (aux_gm.q_ice[k] < 0)
                    @debug "Negative value in grid mean q_liq or q_ice ($(aux_gm.q_liq[k]), $(aux_gm.q_ice[k]))"
                    TS.isoutofdomain = true
                    return true
                end
            end
            if (aux_gm.T[k] < 0) || (aux_gm.RH[k] < 0)
                @debug "Negative value in grid mean RH or T ($(aux_gm.RH[k]), $(aux_gm.T[k]))"
                TS.isoutofdomain = true
                return true
            end
        end

        # check for negative values in environment auxiliary variables
        @inbounds for k in TC.real_center_indices(grid)
            if (aux_en.θ_liq_ice[k] < 0) || (aux_en.q_tot[k] < 0)
                @debug "Negative value in environment θ_liq_ice or q_tot ($(aux_en.θ_liq_ice[k]), $(aux_en.q_tot[k]))"
                TS.isoutofdomain = true
                return true
            end
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                if (aux_en.q_liq[k] < 0) || (aux_en.q_ice[k] < 0)
                    @debug "Negative value in environment q_liq or q_ice ($(aux_en.q_liq[k]), $(aux_en.q_ice[k]))"
                    TS.isoutofdomain = true
                    return true
                end
            end
        end

        @debug "No negative values found in prognostic variables in isoutofdomain() at t = $t"
        return false # default to returning false from the do block
    end

    @debug "TS.isoutofdomain = $(params.TS.isoutofdomain)"
    return params.TS.isoutofdomain
end

"""
    While well intentioned, this function is not truly good bc it doesn't feed back on gm...
    You would need to assess the change in these tendencies and apply them to the grid mean variables. 
    A lot of that work is done in compute_gm_tendencies!()

    Problematically, it's not clear what limiting the gm tendencies have... so unless that's NoTendencyLimiter(), not recalculating could be an error...
    The safest (and easiest to maintain and keep bug-free) solution, despite being slower, is to recalculate all the tendencies.
"""
# function update_noneq_moisture_sources_tendencies!(tendencies::FV, prog::FV, params::NT, t::Real, use_fallback_tendency_limiters::Bool) where {NT, FV <: CC.Fields.FieldVector}
#     error("this function has been deprecated")
#     return nothing
# end

function ∑tendencies_null!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    # do nothing
    return nothing
end

function ∑tendencies_robust!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    # fails gracefully by setting all tendencies to Inf, for use w/ adaptive timestepping that rejects steps

    if isoutofdomain(prog, params, t)
        @debug "Negative values in prognostic variables detected. Setting all tendencies to Inf so timestep can be gracefully rejected."
        CC.Fields.bycolumn(axes(prog.cent)) do colidx
            # UnPack.@unpack edmf, precip_model, param_set, case = params
            # UnPack.@unpack surf_params, radiation, forcing, aux, TS = params
            # thermo_params = TCP.thermodynamics_params(param_set)

            parent(tendencies.face[colidx]) .= -Inf # reset tendencies to 0?
            parent(tendencies.cent[colidx]) .= -Inf # reset tendencies to 0?
        end
    else
        ∑tendencies!(tendencies, prog, params, t)
    end
end

# Compute the sum of tendencies for the scheme
function ∑tendencies!(tendencies::FV, prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    CC.Fields.bycolumn(axes(prog.cent)) do colidx
        UnPack.@unpack edmf, precip_model, param_set, case = params
        UnPack.@unpack surf_params, radiation, forcing, aux, TS = params

        # thermo_params = TCP.thermodynamics_params(param_set)

        parent(tendencies.face[colidx]) .= 0 # reset tendencies to 0?
        parent(tendencies.cent[colidx]) .= 0 # reset tendencies to 0?
        state = TC.column_state(prog, aux, tendencies, colidx, params.calibrate_io) # also creates face and cent

        # TC.filter_gm_vars(edmf, state)  # filter gm vars before setting thermo state so that we don't have to do it again later
        # TC.affect_filter!(edmf, state, param_set, surf, cfl_limit, Δt) # This should alos be done in the filter callback...

        # set_thermo_state_from_prog!(state, edmf.moisture_model, param_set) # moved to update_aux_caller!()
        # assign_thermo_aux!(state, param_set) # moved to update_aux_caller!()
        # aux_gm = TC.center_aux_grid_mean(state) # moved to update_aux_caller!()
        # @. aux_gm.θ_virt = TD.virtual_pottemp(thermo_params, aux_gm.ts) # moved to update_aux_caller!()

        # Δt = TS.dt 
        Δt = TS.dt_limit_tendencies_factor * (TS.limit_tendencies_by_dt_min ? TS.dt_min : TS.dt)
        @debug "calc'ing tends: t = $t, Δt = $Δt, use_fallback_tendency_limiters = $(TS.use_fallback_tendency_limiters)"
        cfl_limit = TS.cfl_limit

        surf = get_surface(surf_params, state, t, param_set)

        # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
        # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
        # treat them as auxiliary variables for now, until we disentangle the tendency computations.
        if isa(case, Cases.SOCRATES)
            Cases.update_forcing(case, state, t, param_set, forcing) # should dispatch to socrates
        else
            Cases.update_forcing(case, state, t, param_set)
        end
        Cases.update_radiation(radiation, state, t, param_set)

        # TC.update_aux!(edmf, state, surf, param_set, t, Δt, cfl_limit, TS.use_fallback_tendency_limiters) # update aux vars, reset tendencies... call to microphsyics!() also resets some aux vars and calculates en tendencies
        TC.update_aux_tendencies!(edmf, state, param_set, Δt, TS.use_fallback_tendency_limiters) # move update aux to after step, in a callback I guess
        # compute tendencies
        # causes division error in dry bubble first time step
        TC.compute_precipitation_sink_tendencies(precip_model, edmf, state, param_set, Δt, TS.use_fallback_tendency_limiters)
        TC.compute_precipitation_advection_tendencies(precip_model, edmf, state, param_set, TS.use_fallback_tendency_limiters)

        TC.compute_turbconv_tendencies!(edmf, state, param_set, surf, Δt, cfl_limit, TS.use_fallback_tendency_limiters)
        compute_gm_tendencies!(edmf, state, surf, radiation, forcing, param_set, TS.use_fallback_tendency_limiters, Δt)

        # Now that we have updraft, gm from updraft, we can calculate any post tendencies that we need for things that are backed out.

        # I haven't come up with a limiter that's better than filtering... (we can stop the model from crashing but you're essentially getting random results, just use tendency-adaptive dt)
    end

    return nothing
end


# Compute the sum of tendencies for the scheme
function update_aux_caller!(prog::FV, params::NT, t::Real) where {NT, FV <: CC.Fields.FieldVector}
    CC.Fields.bycolumn(axes(prog.cent)) do colidx
        UnPack.@unpack edmf, precip_model, param_set, case = params
        UnPack.@unpack surf_params, radiation, forcing, aux, TS = params

        thermo_params = TCP.thermodynamics_params(param_set)
        # state = TC.column_state(prog, aux, tendencies, colidx, params.calibrate_io) # also creates face and cent
        state = TC.column_prog_aux(prog, aux, colidx, params.calibrate_io)

        # filter gm_vars should have been done in the filter callback

        set_thermo_state_from_prog!(state, edmf.moisture_model, param_set) # sets ts_gm using prog variables
        assign_thermo_aux!(state, param_set) # sets aux_gm thermo vars
        aux_gm = TC.center_aux_grid_mean(state)
        @. aux_gm.θ_virt = TD.virtual_pottemp(thermo_params, aux_gm.ts)

        Δt = TS.dt_limit_tendencies_factor * (TS.limit_tendencies_by_dt_min ? TS.dt_min : TS.dt)
        @debug "calc'ing tends: t = $t, Δt = $Δt, use_fallback_tendency_limiters = $(TS.use_fallback_tendency_limiters)"
        cfl_limit = TS.cfl_limit

        surf = get_surface(surf_params, state, t, param_set)

        # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
        # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
        # treat them as auxiliary variables for now, until we disentangle the tendency computations.

        TC.update_aux!(edmf, state, surf, param_set, Δt, cfl_limit) # update aux vars, reset tendencies... call to microphsyics!() also resets some aux vars and calculates en tendencies

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


# This function is the "Firewall".
# The compiler analyzes this ONE time. 
# It generates code that accepts ANY val_bottom/val_top.
# This function handles the heavy lifting.
# The compiler compiles this ONCE for any Float64 inputs.
# ------------------------------------------------------------------
#  HELPER FUNCTION (Strict Types Enforced)
# ------------------------------------------------------------------
function apply_subsidence_barrier!(
    tendency::F, subsidence::F, scalar::F, scalar_F::F2, scalar_F2::F2, val_bottom::FT, val_top::FT, prefactor::ForFT
) where {FT, F <: CC.Fields.Field, F2 <: CC.Fields.Field, ForFT <: Union{F, FT}} # Strictly enforces ClimaCore Fields
    InterpScalar = CCO.InterpolateC2F(bottom = CCO.SetValue(val_bottom), top = CCO.SetValue(val_top))
    InterpW = CCO.InterpolateC2F(bottom = CCO.SetValue(zero(FT)), top = CCO.Extrapolate())
    Div = CCO.DivergenceF2C()
    Upwind = CCO.UpwindBiasedProductC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    Reconstruct = TC.Ic # Global const
    W = CCG.WVector

    # @. scalar_F = InterpScalar(scalar)
    # @. scalar = Div(W(scalar_F))
    # @. scalar_F2 = InterpW(subsidence)  

    # scalar_F = @. TC.lazy(InterpScalar(scalar))
    # scalar = @. TC.lazy(Div(W(scalar_F)))
    # scalar_F2 = @. TC.lazy(InterpW(subsidence))

    # @. tendency -= prefactor * (Reconstruct(W(Upwind(W(scalar_F2), scalar))).components.data.:1) # second W seems to be important for units
    # @. tendency -= prefactor * TC.toscalar((Reconstruct(W(Upwind(W(scalar_F2), scalar))))) # second W seems to be important for units

    # L = TC.lazy
    # @lazy = TC.@lazy
    # out = @. L(prefactor * L(Reconstruct(L(TC.toscalar(L(W(L(Upwind(L(W(L(InterpW(subsidence)), L(Div(L(W(L(InterpScalar(scalar)))))))))))))))))

    # out = @lazy @. (prefactor * (Reconstruct((TC.toscalar((W((Upwind((W((InterpW(subsidence)), (Div((W((InterpScalar(scalar)))))))))))))))))
    # @. tendency -= Base.materialize(out)


    
    @. scalar_F = InterpScalar(scalar) # 1. Fill Face Scratch (Reads scalar, Writes scalar_F)  ( MUST be done first before we overwrite scalar.)
    @. scalar = Div(W(scalar_F)) # 2. Fill Center Scratch (Reads scalar_F, Overwrites scalar) ( We use 'scalar' itself to store the Divergence result, as you did.)
    @. scalar_F2 = InterpW(subsidence) # 3. Fill Face Scratch (Reads subsidence, Writes scalar_F2)
    out = @. TC.lazy(prefactor * Reconstruct(TC.toscalar(Upwind(W(scalar_F2), scalar)))) # 4. Final Fused Tendency Update [[ We use TC.lazy ONLY here.  This fuses Reconstruct + toscalar + prefactor into the Upwind loop. Upwind reads the materialized 'scalar_F2' and 'scalar' from above.
    @. tendency -= Base.materialize(out)
        

    return nothing
end
function save_subsidence_barrier!(dest::F, subsidence::F, scalar::F, scalar_F::F2, scalar_F2::F2, val_bottom::FT, val_top::FT, prefactor::ForFT
) where {FT, F <: CC.Fields.Field, F2 <: CC.Fields.Field, ForFT <: Union{F, FT}}
    InterpScalar = CCO.InterpolateC2F(bottom = CCO.SetValue(val_bottom), top = CCO.SetValue(val_top))
    InterpW = CCO.InterpolateC2F(bottom = CCO.SetValue(zero(FT)), top = CCO.Extrapolate())
    Div = CCO.DivergenceF2C()
    Upwind = CCO.UpwindBiasedProductC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    Reconstruct = TC.Ic # Global const
    W = CCG.WVector


    @. scalar_F = InterpScalar(scalar)
    @. scalar = Div(W(scalar_F))
    @. scalar_F2 = InterpW(subsidence)

    # @. dest = -prefactor * (Reconstruct(W(Upwind(W(scalar_F2), scalar))).components.data.:1) # second W seems to be important for units
    @. dest = -prefactor * TC.toscalar((Reconstruct(W(Upwind(W(scalar_F2), scalar))))) # second W seems to be important for units

    return nothing
end

function compute_gm_tendencies!(
    edmf::TC.EDMFModel,
    state::TC.State,
    surf::TC.SurfaceBase,
    radiation::Cases.RadiationBase,
    force::Cases.ForcingBase,
    param_set::APS,
    use_fallback_tendency_limiters::Bool,
    Δt::Real,
)
    grid = TC.Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    tendencies_gm = TC.center_tendencies_grid_mean(state)
    tendencies_pr = TC.center_tendencies_precipitation(state)
    kc_toa = TC.kc_top_of_atmos(grid)
    kf_surf = TC.kf_surface(grid)
    FT = TC.float_type(state)
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    prog_pr = TC.center_prog_precipitation(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    aux_en = TC.center_aux_environment(state)
    aux_en_f = TC.face_aux_environment(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_bulk = TC.center_aux_bulk(state)
    ρ_f = aux_gm_f.ρ
    p_c = aux_gm.p
    ρ_c = prog_gm.ρ
    aux_tc = TC.center_aux_turbconv(state)
    aux_tc_f = TC.face_aux_turbconv(state)
    ts_gm = TC.center_aux_grid_mean(state).ts
    N_up = TC.n_updrafts(edmf)

    θ_liq_ice_gm_toa = prog_gm.ρθ_liq_ice[kc_toa] / ρ_c[kc_toa]
    q_tot_gm_toa = prog_gm.ρq_tot[kc_toa] / ρ_c[kc_toa]



    wvec = CC.Geometry.WVector
    ∇c = TC.∇c

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

    f = TC.center_aux_turbconv(state).k̂
    @. f = CCG.Contravariant3Vector(coriolis_fn(coords))


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

    # UBsub = CCO.UpwindBiasedProductC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # upwinding, extrapolate bc we don't know the boa/toa derivatives

    kc_surf = TC.kc_surface(grid)
    θ_liq_ice_gm_boa = prog_gm.ρθ_liq_ice[kc_surf] / ρ_c[kc_surf]
    q_tot_gm_boa = prog_gm.ρq_tot[kc_surf] / ρ_c[kc_surf]

    # subsidence should be face valued so we'll need a C2F call
    # C2Fsub = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.Extrapolate()) # not sure but maybe we should use this for stability as it's not biased for your upwind algorithm... no penetration. #

    # C2Fθ = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(θ_liq_ice_gm_boa)), top = CCO.SetValue(FT(θ_liq_ice_gm_toa))) # not sure if this should be right biased or not
    # C2Fq = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(q_tot_gm_boa)), top = CCO.SetValue(FT(q_tot_gm_toa))) # not sure if this should be right biased or not

    # we to C2F, then ∇c goes F2C, then UBsub got C2F, but we wanna end on C [ there's no UpwindBiasedProductF2C and we had UpwindBiasedProductF2C(u[F], x[C])
    # F2Csub = CCO.InterpolateF2C(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # We might have put 0 for no penetration on boundary, but that's not exactly true in our dataset...
    # CV32FT = x -> x[1] # convert Contravariant3Vector to Float64 


    TC.zero_field!(aux_tc.dqvdt) # reset to zero, do here instead of in update_aux bc we need only some of these values...
    TC.zero_field!(aux_tc.dTdt) # reset to zero, do here instead of in update_aux bc we need only some of these values...
    

    # We have access to three center temp vars: ϕ, ψ, and Φ, and one face temp var ϕ [[ These are okay because we do not make any nested calls here... ]]
    # ∇Tr = aux_tc.temporary_1 # reusing temporary storage for w∇ stuff :: Tracer
    # w∇Tr = aux_tc.temporary_2 # reusing temporary storage for w∇ stuff :: w∇Tracer
    # w∇TrF = aux_tc_f.temporary_f1 # reusing temporary storage for w∇ stuff on faces :: w∇Tracer on faces
    # w∇CV3 = similar(f)



    # considering deprecation pending Dennis's response (or dividing by dz)

    #=
        I don't know for sure why but this version works much better
        I think the problem is that convergence in the other one isn't reflected in density so you get weird temperature and other changes
        Because we ignore density changes, we really do need to stick to this one and just the w∇q part, and leave out the q∇w part...

        :: We also always divide out the ρ_c to allow for adiabatic compression. Otherwise, advecting ρq downwards for example could imply a reduction in q, when really it will just compress. It's not like true advection e.g. in updrafts, since grid mean ρ is fixed
    =#

    # regular upwinding 
    # ∇θ_liq_ice_gm = ∇Tr # alias
    # w∇θ_liq_ice_gm_F = w∇TrF # alias
    # w∇θ_liq_ice_gm_CV3 = w∇CV3 # alias
    # w∇θ_liq_ice_gm = w∇Tr # alias
    # @. ∇θ_liq_ice_gm = ∇c(wvec(C2Fθ(prog_gm.ρθ_liq_ice / ρ_c)))
    # @. w∇θ_liq_ice_gm_F = wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇θ_liq_ice_gm)) # works. outer wvec important for getting physical units [[ # w∇θ_liq_ice_gm = @. UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇θ_liq_ice_gm) # works, but output type is different than tendencies type ]]
    # @. w∇θ_liq_ice_gm_CV3 = F2Csub(w∇θ_liq_ice_gm_F) # convert back to C
    # @. w∇θ_liq_ice_gm = w∇θ_liq_ice_gm_CV3.components.data.:1 # from Dennis [[ x = @. CV32FT(x) # convert back to Float64 from Contravariant3Vector I think maybe also works?]]
    # @. w∇θ_liq_ice_gm = (F2Csub(wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇c(wvec(C2Fθ(prog_gm.ρθ_liq_ice / ρ_c))))))).components.data.:1
    # @. tendencies_gm.ρθ_liq_ice -= ρ_c * w∇θ_liq_ice_gm

    # 1. Reuse one scratch field for division
    w∇θ_liq_ice_gm = aux_tc.temporary_2
    @. w∇θ_liq_ice_gm = prog_gm.ρθ_liq_ice / ρ_c
    apply_subsidence_barrier!(tendencies_gm.ρθ_liq_ice, aux_gm.subsidence, w∇θ_liq_ice_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, θ_liq_ice_gm_boa, θ_liq_ice_gm_toa, ρ_c)


    # ∇q_tot_gm = ∇Tr # alias
    # w∇q_tot_gm_F = w∇TrF # alias
    # w∇q_tot_gm_CV3 = w∇CV3 # alias
    # w∇q_tot_gm = w∇Tr # alias
    # @. ∇q_tot_gm = ∇c(wvec(C2Fq(prog_gm.ρq_tot / ρ_c)))
    # @. w∇q_tot_gm_F = wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_tot_gm)) # works, outer wvec important for getting physical units
    # @. w∇q_tot_gm_CV3 = F2Csub(w∇q_tot_gm_F) # convert back to C
    # @. w∇q_tot_gm = w∇q_tot_gm_CV3.components.data.:1 # from Dennis
    # @. w∇q_tot_gm = (F2Csub(wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇c(wvec(C2Fq(prog_gm.ρq_tot / ρ_c))))))).components.data.:1
    # @. tendencies_gm.ρq_tot -= ρ_c * w∇q_tot_gm
    w∇q_tot_gm = aux_tc.temporary_2
    @. w∇q_tot_gm = prog_gm.ρq_tot / ρ_c
    if state.calibrate_io
        apply_subsidence_barrier!(tendencies_gm.ρq_tot, aux_gm.subsidence, w∇q_tot_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_tot_gm_boa, q_tot_gm_toa, ρ_c)
    else
        save_subsidence_barrier!(aux_tc.qt_tendency_ls_vert_adv, aux_gm.subsidence, w∇q_tot_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_tot_gm_boa, q_tot_gm_toa, one(FT))
        @. tendencies_gm.ρq_tot -= ρ_c * aux_tc.qt_tendency_ls_vert_adv
    end
    # if !state.calibrate_io
    #     @. aux_tc.qt_tendency_ls_vert_adv = -w∇q_tot_gm # FOR STORAGE
    # end

    @. aux_tc.dqvdt -= w∇q_tot_gm # grid mean contribution to dqt/dt

    if edmf.moisture_model isa TC.NonEquilibriumMoisture
        q_liq_gm_toa = prog_gm.q_liq[kc_toa]
        q_ice_gm_toa = prog_gm.q_ice[kc_toa]
        q_liq_gm_boa = prog_gm.q_liq[kc_surf]
        q_ice_gm_boa = prog_gm.q_ice[kc_surf]
        # C2Fq_liq = CCO.InterpolateC2F(bottom = CCO.SetValue(FT(q_liq_gm_boa)), top = CCO.SetValue(FT(q_liq_gm_toa))) # not sure if this should be right biased or not
        # C2Fq_ice = CCO.InterpolateC2F(bottom = CCO.SetValue(FT(q_ice_gm_boa)), top = CCO.SetValue(FT(q_ice_gm_toa))) # not sure if this should be right biased or not

        # ∇q_liq_gm = ∇Tr # alias
        # w∇q_liq_gm_F = w∇TrF # alias
        # w∇q_liq_gm_CV3 = w∇CV3 # alias
        # w∇q_liq_gm = w∇Tr # alias
        # @. ∇q_liq_gm = ∇c(wvec(C2Fq_liq(prog_gm.q_liq)))
        # @. w∇q_liq_gm_F = wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_liq_gm)) # works, outer wvec important for getting physical units
        # @. w∇q_liq_gm_CV3 = F2Csub(w∇q_liq_gm_F) # convert back to C
        # @. w∇q_liq_gm = w∇q_liq_gm_CV3.components.data.:1 # from Dennis
        # @. w∇q_liq_gm = (F2Csub(wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇c(wvec(C2Fq_liq(prog_gm.q_liq))))))).components.data.:1
        # @. tendencies_gm.q_liq -= w∇q_liq_gm
        w∇q_liq_gm = aux_tc.temporary_2
        @.w∇q_liq_gm = prog_gm.q_liq
        if state.calibrate_io
            apply_subsidence_barrier!(tendencies_gm.q_liq, aux_gm.subsidence, w∇q_liq_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_liq_gm_boa, q_liq_gm_toa, one(FT))
        else
            save_subsidence_barrier!(aux_tc.ql_tendency_ls_vert_adv, aux_gm.subsidence, w∇q_liq_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_liq_gm_boa, q_liq_gm_toa, one(FT))
            @. tendencies_gm.q_liq -= aux_tc.ql_tendency_ls_vert_adv
        end
        # if !state.calibrate_io
        #     @. aux_tc.ql_tendency_ls_vert_adv = -w∇q_liq_gm # FOR STORAGE
        # end

        # ∇q_ice_gm = ∇Tr # alias
        # w∇q_ice_gm_F = w∇TrF # alias
        # w∇q_ice_gm_CV3 = w∇CV3 # alias
        # w∇q_ice_gm = w∇Tr # alias
        # @. ∇q_ice_gm = ∇c(wvec(C2Fq_ice(prog_gm.q_ice)))
        # @. w∇q_ice_gm_F = wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_ice_gm)) # works, outer wvec important for getting physical units
        # @. w∇q_ice_gm_CV3 = F2Csub(w∇q_ice_gm_F) # convert back to C
        # @. w∇q_ice_gm = w∇q_ice_gm_CV3.components.data.:1 # from Dennis
        # @. w∇q_ice_gm = (F2Csub(wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇c(wvec(C2Fq_ice(prog_gm.q_ice))))))).components.data.:1
        # @. tendencies_gm.q_ice -= w∇q_ice_gm
        w∇q_ice_gm = aux_tc.temporary_2
        @. w∇q_ice_gm = prog_gm.q_ice
        
        if state.calibrate_io
            apply_subsidence_barrier!(tendencies_gm.q_ice, aux_gm.subsidence, w∇q_ice_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_ice_gm_boa, q_ice_gm_toa, one(FT))
        else
            save_subsidence_barrier!(aux_tc.qi_tendency_ls_vert_adv, aux_gm.subsidence, w∇q_ice_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_ice_gm_boa, q_ice_gm_toa, one(FT))
            @. tendencies_gm.q_ice -= aux_tc.qi_tendency_ls_vert_adv
        end
        # if !state.calibrate_io
        #     @. aux_tc.qi_tendency_ls_vert_adv = -w∇q_ice_gm # FOR STORAGE
        # end
    end

    # LS Vert for precipitation
    q_rai_gm_toa = prog_pr.q_rai[kc_toa]
    q_rai_gm_boa = prog_pr.q_rai[kc_surf]
    q_sno_gm_toa = prog_pr.q_sno[kc_toa]
    q_sno_gm_boa = prog_pr.q_sno[kc_surf]

    # ∇q_rai_gm = ∇Tr # alias
    # w∇q_rai_gm_F = w∇TrF # alias
    # w∇q_rai_gm_CV3 = w∇CV3 # alias
    # w∇q_rai_gm = w∇Tr # alias
    # @. ∇q_rai_gm = ∇c(wvec(C2Fsub(prog_pr.q_rai * ρ_c))) / ρ_c
    # @. w∇q_rai_gm_F = wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_rai_gm))
    # @. w∇q_rai_gm_CV3 = F2Csub(w∇q_rai_gm_F) # convert back to C
    # @. w∇q_rai_gm = w∇q_rai_gm_CV3.components.data.:1 # from Dennis
    # @. w∇q_rai_gm = (F2Csub(wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇c(wvec(C2Fsub(prog_pr.q_rai))))))).components.data.:1
    # @. tendencies_pr.q_rai += -w∇q_rai_gm
    w∇q_rai_gm = aux_tc.temporary_2
    @. w∇q_rai_gm = prog_pr.q_rai
    if state.calibrate_io
        apply_subsidence_barrier!(tendencies_pr.q_rai, aux_gm.subsidence, w∇q_rai_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_rai_gm_boa, q_rai_gm_toa, one(FT))
    else
        save_subsidence_barrier!(aux_tc.qr_tendency_ls_vert_adv, aux_gm.subsidence, w∇q_rai_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_rai_gm_boa, q_rai_gm_toa, one(FT))
        @. tendencies_pr.q_rai -= aux_tc.qr_tendency_ls_vert_adv
    end
    # if !state.calibrate_io
    #     @. aux_tc.qr_tendency_ls_vert_adv = -w∇q_rai_gm # FOR STORAGE
    # end

    # ∇q_sno_gm = ∇Tr # alias
    # w∇q_sno_gm_F = w∇TrF # alias
    # w∇q_sno_gm_CV3 = w∇CV3 # alias
    # w∇q_sno_gm = w∇Tr # alias
    # @. ∇q_sno_gm = ∇c(wvec(C2Fsub(prog_pr.q_sno * ρ_c))) / ρ_c
    # @. w∇q_sno_gm_F = wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇q_sno_gm))
    # @. w∇q_sno_gm_CV3 = F2Csub(w∇q_sno_gm_F) # convert back to C
    # @. w∇q_sno_gm = w∇q_sno_gm_CV3.components.data.:1 # from Dennis
    # @. w∇q_sno_gm = (F2Csub(wvec(UBsub(wvec(C2Fsub(aux_gm.subsidence)), ∇c(wvec(C2Fsub(prog_pr.q_sno))))))).components.data.:1
    w∇q_sno_gm = aux_tc.temporary_2
    @. w∇q_sno_gm = prog_pr.q_sno
    if state.calibrate_io
        apply_subsidence_barrier!(tendencies_pr.q_sno, aux_gm.subsidence, w∇q_sno_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_sno_gm_boa, q_sno_gm_toa, one(FT))
    else
        save_subsidence_barrier!(aux_tc.qs_tendency_ls_vert_adv, aux_gm.subsidence, w∇q_sno_gm, aux_tc_f.temporary_f1, aux_tc_f.temporary_f2, q_sno_gm_boa, q_sno_gm_toa, one(FT))
        @. tendencies_pr.q_sno -= aux_tc.qs_tendency_ls_vert_adv
    end
    # apply_subsidence_barrier!(tendencies_pr.q_sno, w∇q_sno_gm, aux_gm.subsidence, q_sno_gm_boa, q_sno_gm_toa, -FT(1))
    # @. tendencies_pr.q_sno += -w∇q_sno_gm
    # if !state.calibrate_io
    #     @. aux_tc.qs_tendency_ls_vert_adv = -w∇q_sno_gm # FOR STORAGE
    # end



    @inbounds for k in TC.real_center_indices(grid)
        Π = TD.exner(thermo_params, ts_gm[k])

        # Radiation
        if Cases.rad_type(radiation) <: Union{Cases.RadiationDYCOMS_RF01, Cases.RadiationLES, Cases.RadiationTRMM_LBA, Cases.RadiationSOCRATES}
            tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * aux_gm.dTdt_rad[k] / Π 
            aux_tc.dTdt[k] += aux_gm.dTdt_rad[k] # [store things that matter for our dTdt] [ dTdt priority ] [ idk if this totally helps, it's typically strong cloud top cooling so maybe it drives liquid generation?]
        end

        # LS advection
        tendencies_gm.ρq_tot[k] += ρ_c[k] * aux_gm.dqtdt_hadv[k]
        if !(Cases.force_type(force) <: Cases.ForcingDYCOMS_RF01)
            tendencies_gm.ρθ_liq_ice[k] += ρ_c[k] * aux_gm.dTdt_hadv[k] / Π 
            aux_tc.dTdt[k] += aux_gm.dTdt_hadv[k] # [store things that matter for our dTdt] [ dTdt priority ]

        end
        if edmf.moisture_model isa TC.NonEquilibriumMoisture
            tendencies_gm.q_liq[k] += aux_gm.dqldt[k] # we never set these, e.g. in socrates... presumably they're zero?
            tendencies_gm.q_ice[k] += aux_gm.dqidt[k]
        end

        aux_tc.dqvdt[k] += aux_gm.dqtdt_hadv[k] # [store things that matter for our dqvdt] [ dqvdt priority ]

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

            aux_tc.dqvdt[k] += (gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k]) # [store things that matter for our dqvdt] [  dqvdt priority ]

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

            aux_gm.dqtdt_nudge[k] = gm_q_tot_nudge_k

            aux_tc.dqvdt[k] += (gm_q_tot_nudge_k + aux_gm.dqtdt_fluc[k]) # [store things that matter for our dqvdt] [dqvdt priority ]
            aux_tc.dTdt[k] += (gm_H_nudge_k* Π + aux_gm.dTdt_fluc[k]) #  [store things that matter for our aux_tc.dTdt] [dTdt priority ] [ go opposite way with exner] # DO NOT DO THIS, it's not at all the same thing as temperature!!! because of condensate. it's really not the same thing as θ_liq_ice, so we don't do this.

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



        # ====== sedimentation grid mean only (stability fix...)
        if edmf.cloud_sedimentation_model isa TC.CloudSedimentationModel 
            if edmf.cloud_sedimentation_model.grid_mean
                error("Grid mean sedimentation not implemented yet")
            else

            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                tendencies_gm.q_liq[k] +=
                    aux_en.ql_tendency_sedimentation[k] + aux_bulk.ql_tendency_sedimentation[k]
                tendencies_gm.q_ice[k] +=
                    aux_en.qi_tendency_sedimentation[k] + aux_bulk.qi_tendency_sedimentation[k]
            end

            tendencies_gm.ρq_tot[k] +=
                ρ_c[k] * (aux_bulk.qt_tendency_sedimentation[k] + aux_en.qt_tendency_sedimentation[k])

            tendencies_gm.ρθ_liq_ice[k] +=
                ρ_c[k] * (aux_bulk.θ_liq_ice_tendency_sedimentation[k] + aux_en.θ_liq_ice_tendency_sedimentation[k])
            end
        end
    end


    # =================================================== #
    # note, these `sgs` fluxes also include the mass fluxes from vertical advection

    TC.compute_sgs_flux!(edmf, state, surf, param_set, Δt)
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

        # sgs_tendency_q_liq = aux_gm.sgs_tendency_q_liq
        # sgs_tendency_q_ice = aux_gm.sgs_tendency_q_ice
        # @. sgs_tendency_q_liq = -∇sgs(wvec(sgs_flux_q_liq)) / ρ_c
        # @. sgs_tendency_q_ice = -∇sgs(wvec(sgs_flux_q_ice)) / ρ_c
        # @. tends_q_liq += sgs_tendency_q_liq
        # @. tends_q_ice += sgs_tendency_q_ice

        @. tends_q_liq += -∇sgs(wvec(sgs_flux_q_liq)) / ρ_c
        @. tends_q_ice += -∇sgs(wvec(sgs_flux_q_ice)) / ρ_c

    end

    # my addition -- note this is the sgs flux only, doesn't include advection
    sgs_flux_q_rai = aux_gm_f.sgs_flux_q_rai
    sgs_flux_q_sno = aux_gm_f.sgs_flux_q_sno
    # sgs_tendency_q_rai = aux_gm.sgs_tendency_q_rai
    # sgs_tendency_q_sno = aux_gm.sgs_tendency_q_sno

    # @. sgs_tendency_q_rai = -∇sgs(wvec(sgs_flux_q_rai)) / ρ_c
    # @. sgs_tendency_q_sno = -∇sgs(wvec(sgs_flux_q_sno)) / ρ_c
    # @. tendencies_pr.q_rai += sgs_tendency_q_rai
    # @. tendencies_pr.q_sno += sgs_tendency_q_sno

    @. tendencies_pr.q_rai += -∇sgs(wvec(sgs_flux_q_rai)) / ρ_c
    @. tendencies_pr.q_sno += -∇sgs(wvec(sgs_flux_q_sno)) / ρ_c



    return nothing
end

