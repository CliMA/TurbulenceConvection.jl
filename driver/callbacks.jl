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
    UnPack.@unpack edmf, calibrate_io, precip_model, aux, grid, io_nt, diagnostics, case, param_set, Stats, skip_io =
        integrator.p
    skip_io && return nothing
    t = integrator.t

    # TODO: remove `vars` hack that avoids
    # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    # opening/closing files every step should be okay. #removeVarsHack
    # TurbulenceConvection.io(sim) # #removeVarsHack
    write_simulation_time(Stats, t) # #removeVarsHack

    state = TC.State(integrator.u, aux, ODE.get_du(integrator))

    # TODO: is this the best location to call diagnostics?
    compute_diagnostics!(edmf, precip_model, param_set, grid, state, diagnostics, Stats, case, t, calibrate_io)

    cent = TC.Cent(1)
    diag_svpc = svpc_diagnostics_grid_mean(diagnostics)
    diag_tc_svpc = svpc_diagnostics_turbconv(diagnostics)
    write_ts(Stats, "lwp_mean", diag_svpc.lwp_mean[cent])
    write_ts(Stats, "iwp_mean", diag_svpc.iwp_mean[cent])
    write_ts(Stats, "rwp_mean", diag_svpc.rwp_mean[cent])
    write_ts(Stats, "swp_mean", diag_svpc.swp_mean[cent])

    if !calibrate_io
        write_ts(Stats, "updraft_cloud_cover", diag_tc_svpc.updraft_cloud_cover[cent])
        write_ts(Stats, "updraft_cloud_base", diag_tc_svpc.updraft_cloud_base[cent])
        write_ts(Stats, "updraft_cloud_top", diag_tc_svpc.updraft_cloud_top[cent])
        write_ts(Stats, "env_cloud_cover", diag_tc_svpc.env_cloud_cover[cent])
        write_ts(Stats, "env_cloud_base", diag_tc_svpc.env_cloud_base[cent])
        write_ts(Stats, "env_cloud_top", diag_tc_svpc.env_cloud_top[cent])
        write_ts(Stats, "env_lwp", diag_tc_svpc.env_lwp[cent])
        write_ts(Stats, "env_iwp", diag_tc_svpc.env_iwp[cent])
        write_ts(Stats, "Hd", diag_tc_svpc.Hd[cent])
        write_ts(Stats, "updraft_lwp", diag_tc_svpc.updraft_lwp[cent])
        write_ts(Stats, "updraft_iwp", diag_tc_svpc.updraft_iwp[cent])

        write_ts(Stats, "cutoff_precipitation_rate", diag_svpc.cutoff_precipitation_rate[cent])
        write_ts(Stats, "cloud_cover_mean", diag_svpc.cloud_cover_mean[cent])
        write_ts(Stats, "cloud_base_mean", diag_svpc.cloud_base_mean[cent])
        write_ts(Stats, "cloud_top_mean", diag_svpc.cloud_top_mean[cent])

        write_ts(Stats, "integ_total_flux_qt", diag_svpc.integ_total_flux_qt[cent])
        write_ts(Stats, "integ_total_flux_s", diag_svpc.integ_total_flux_s[cent])
    end

    io(io_nt.aux, Stats, state)
    io(io_nt.diagnostics, Stats, diagnostics)

    surf = get_surface(case.surf_params, grid, state, t, param_set)
    io(surf, case.surf_params, grid, state, Stats, t)

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end

function affect_filter!(integrator)
    UnPack.@unpack edmf, grid, param_set, aux, case = integrator.p
    t = integrator.t
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    surf = get_surface(case.surf_params, grid, state, t, param_set)
    TC.affect_filter!(edmf, grid, state, param_set, surf, case.casename, t)

    # We're lying to OrdinaryDiffEq.jl, in order to avoid
    # paying for an additional `∑tendencies!` call, which is required
    # to support supplying a continuous representation of the
    # solution.
    ODE.u_modified!(integrator, false)
end

function adaptive_dt!(integrator)
    UnPack.@unpack edmf, TS, dt_min = integrator.p
    TS.dt = min(TS.dt_max, max(TS.dt_max_edmf, dt_min))
    SciMLBase.set_proposed_dt!(integrator, TS.dt)
    ODE.u_modified!(integrator, false)
end

function dt_max!(integrator)
    UnPack.@unpack grid, edmf, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    Δzc = TC.get_Δz(prog_gm.ρ)
    Δzf = TC.get_Δz(prog_gm_f.w)
    CFL_limit = TS.cfl_limit
    N_up = TC.n_updrafts(edmf)

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

function monitor_cfl!(integrator)
    UnPack.@unpack grid, edmf, aux, TS = integrator.p
    state = TC.State(integrator.u, aux, ODE.get_du(integrator))
    prog_gm = TC.center_prog_grid_mean(state)
    Δz = TC.get_Δz(prog_gm.ρ)
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
