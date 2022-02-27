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
    UnPack.@unpack turbconv, calibrate_io, precip_model, aux, grid, io_nt, diagnostics, case, param_set, Stats, skip_io =
        integrator.p
    skip_io && return nothing
    t = integrator.t

    state = TC.State(integrator.u, aux, ODE.get_du(integrator))

    # TODO: is this the best location to call diagnostics?
    if !calibrate_io
        compute_diagnostics!(turbconv, precip_model, param_set, grid, state, diagnostics, Stats, case, t)
    end

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack
    # TurbulenceConvection.io(sim) # #removeVarsHack
    write_simulation_time(Stats, t) # #removeVarsHack

    io(io_nt.aux, Stats, state)
    io(io_nt.diagnostics, Stats, diagnostics)

    if !calibrate_io
        surf = get_surface(case.surf_params, grid, state, t, param_set)
        io(surf, case.surf_params, grid, state, Stats, t)
    end

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end

function affect_filter!(integrator)
    UnPack.@unpack turbconv, grid, param_set, aux, case = integrator.p
    t = integrator.t
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    surf = get_surface(case.surf_params, grid, state, t, param_set)
    if turbconv isa TC.EDMFModel
        TC.affect_filter!(turbconv, grid, state, param_set, surf, case.casename, t)
    end

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end

function adaptive_dt!(integrator)
    UnPack.@unpack TS, dt_min = integrator.p
    TS.dt = min(TS.dt_max, max(TS.dt_max_edmf, dt_min))
    SciMLBase.set_proposed_dt!(integrator, TS.dt)
    ODE.u_modified!(integrator, false)
end

function edmf_dt_max!(integrator)
    UnPack.@unpack grid, turbconv, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    Δzc = TC.get_Δz(prog_gm.u)
    Δzf = TC.get_Δz(prog_gm_f.w)
    CFL_limit = TS.cfl_limit
    N_up = TC.n_updrafts(turbconv)

    dt_max = TS.dt_max # initialize dt_max

    aux_tc = TC.center_aux_turbconv(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_en_f = TC.face_aux_environment(state)
    KM = aux_tc.KM
    KH = aux_tc.KH

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_face_indices(grid)
        TC.is_surface_face(grid, k) && continue
        @inbounds for i in 1:N_up
            dt_max = min(dt_max, CFL_limit * Δzf[k] / (abs(aux_up_f[i].w[k]) + eps(Float32)))
        end
        dt_max = min(dt_max, CFL_limit * Δzf[k] / (abs(aux_en_f.w[k]) + eps(Float32)))
    end
    @inbounds for k in TC.real_center_indices(grid)
        vel_max = max(term_vel_rain[k], term_vel_snow[k])
        # Check terminal rain/snow velocity CFL
        dt_max = min(dt_max, CFL_limit * Δzc[k] / (vel_max + eps(Float32)))
        # Check diffusion CFL (i.e., Fourier number)
        dt_max = min(dt_max, CFL_limit * Δzc[k]^2 / (max(KH[k], KM[k]) + eps(Float32)))
    end
    TS.dt_max_edmf = dt_max

    ODE.u_modified!(integrator, false)
end

function diffusivity_dt_max!(integrator)
    UnPack.@unpack grid, turbconv, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    Δzc = TC.get_Δz(prog_gm.u)
    Δzf = TC.get_Δz(prog_gm_f.w)
    CFL_limit = TS.cfl_limit

    dt_max = TS.dt_max # initialize dt_max

    aux_tc = TC.center_aux_turbconv(state)
    aux_gm_f = TC.face_aux_grid_mean(state)

    # # helper to calculate the rain velocity
    # # TODO: assuming gm.W = 0
    # # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_center_indices(grid)
        vel_max = max(term_vel_rain[k], term_vel_snow[k])
        # Check terminal rain/snow velocity CFL
        dt_max = min(dt_max, CFL_limit * Δzc[k] / (vel_max + eps(Float32)))
        # Check diffusion CFL (i.e., Fourier number)
    end
    @inbounds for k in TC.real_face_indices(grid)
        dt_max = min(dt_max, CFL_limit * Δzc[k]^2 / (aux_gm_f.ν[k] + eps(Float32)))
    end
    TS.dt_max_edmf = dt_max

    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(integrator)
    UnPack.@unpack gm, grid, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    prog_gm = TC.center_prog_grid_mean(state)
    Δz = TC.get_Δz(prog_gm.u)
    Δt = TS.dt
    CFL_limit = TS.cfl_limit

    aux_tc = TC.center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in TC.real_center_indices(grid)
        # check stability criterion
        CFL_out_rain = Δt / Δz[k] * term_vel_rain[k]
        CFL_out_snow = Δt / Δz[k] * term_vel_snow[k]
        if TC.is_toa_center(grid, k)
            CFL_in_rain = 0.0
            CFL_in_snow = 0.0
        else
            CFL_in_rain = Δt / Δz[k] * term_vel_rain[k + 1]
            CFL_in_snow = Δt / Δz[k] * term_vel_snow[k + 1]
        end
        if max(CFL_in_rain, CFL_in_snow, CFL_out_rain, CFL_out_snow) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end
    end

    ODE.u_modified!(integrator, false)
end
