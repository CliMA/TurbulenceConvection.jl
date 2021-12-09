function condition_io(u, t, integrator)
    UnPack.@unpack TS, Stats = integrator.p
    TS.dt_io += TS.dt
    io_flag = false
    if TS.dt_io > Stats.frequency
        TS.dt_io = 0
        io_flag = true
    end
    return io_flag || t ≈ 0 || t ≈ TS.t_max
end

condition_every_iter(u, t, integrator) = true

function affect_io!(integrator)
    UnPack.@unpack edmf, aux, grid, io_nt, diagnostics, case, gm, TS, Stats, skip_io = integrator.p
    skip_io && return nothing

    state = TC.State(integrator.u, aux, integrator.du)

    param_set = TC.parameter_set(gm)
    # TODO: is this the best location to call diagnostics?
    compute_diagnostics!(edmf, gm, grid, state, diagnostics, case, TS)

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack
    # TurbulenceConvection.io(sim) # #removeVarsHack
    TC.write_simulation_time(Stats, TS.t) # #removeVarsHack

    TC.io(io_nt.aux, Stats, state)
    TC.io(io_nt.diagnostics, Stats, diagnostics)

    TC.io(gm, grid, state, Stats) # #removeVarsHack
    TC.io(case, grid, state, Stats) # #removeVarsHack
    TC.io(edmf, grid, state, Stats, TS, param_set) # #removeVarsHack

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end

function affect_filter!(integrator)
    UnPack.@unpack edmf, grid, gm, aux, case, TS = integrator.p
    state = TC.State(integrator.u, aux, integrator.du)
    TC.affect_filter!(edmf, grid, state, gm, case, TS)

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end

function adaptive_dt!(integrator)
    UnPack.@unpack edmf, TS, dt_min = integrator.p
    TS.dt = min(TS.dt_max, max(edmf.dt_max, dt_min))
    SciMLBase.set_proposed_dt!(integrator, TS.dt)
    ODE.u_modified!(integrator, false)
end

function dt_max!(integrator)
    UnPack.@unpack gm, grid, edmf, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, integrator.du)
    Δz = grid.Δz
    CFL_limit = TS.cfl_limit
    N_up = TC.n_updrafts(edmf)

    dt_max = TS.dt_max # initialize dt_max

    aux_tc = TC.center_aux_turbconv(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_en_f = TC.face_aux_environment(state)
    KM = aux_tc.KM
    KH = aux_tc.KH

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_face_indices(grid)
        TC.is_surface_face(grid, k) && continue
        @inbounds for i in 1:N_up
            dt_max = min(dt_max, CFL_limit * Δz / (abs(aux_up_f[i].w[k]) + eps(Float32)))
        end
        dt_max = min(dt_max, CFL_limit * Δz / (abs(aux_en_f.w[k]) + eps(Float32)))
    end
    @inbounds for k in TC.real_center_indices(grid)
        vel_max = max(term_vel_rain[k], term_vel_snow[k])
        # Check terminal rain/snow velocity CFL
        dt_max = min(dt_max, CFL_limit * Δz / (vel_max + eps(Float32)))
        # Check diffusion CFL (i.e., Fourier number)
        dt_max = min(dt_max, CFL_limit * Δz^2 / (max(KH[k], KM[k]) + eps(Float32)))
    end
    edmf.dt_max = dt_max

    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(integrator)
    UnPack.@unpack gm, grid, edmf, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, integrator.du)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = TS.cfl_limit

    aux_tc = TC.center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_center_indices(grid)
        # check stability criterion
        CFL_out_rain = Δt / Δz * term_vel_rain[k]
        CFL_out_snow = Δt / Δz * term_vel_snow[k]
        if TC.is_toa_center(grid, k)
            CFL_in_rain = 0.0
            CFL_in_snow = 0.0
        else
            CFL_in_rain = Δt / Δz * term_vel_rain[k + 1]
            CFL_in_snow = Δt / Δz * term_vel_snow[k + 1]
        end
        if max(CFL_in_rain, CFL_in_snow, CFL_out_rain, CFL_out_snow) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end
    end

    ODE.u_modified!(integrator, false)
end
