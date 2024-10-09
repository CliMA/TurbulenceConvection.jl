function condition_io(u, t, integrator)
    UnPack.@unpack TS, Stats = integrator.p
    TS.dt_io += TS.dt
    io_flag = false
    if TS.dt_io > Stats[CC.Fields.ColumnIndex((1, 1), 1)].frequency
        TS.dt_io = 0
        io_flag = true
    end
    return io_flag || t ≈ 0 || t ≈ TS.t_max
end

condition_every_iter(u, t, integrator) = true

function affect_io!(integrator)
    integrator.p.skip_io && return nothing

    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, calibrate_io, precip_model, aux, io_nt, diagnostics, surf_params, param_set, Stats) = integrator.p
        t = integrator.t
        prog = integrator.u
        stats = Stats[colidx]
        # TODO: remove `vars` hack that avoids
        # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
        # opening/closing files every step should be okay. #removeVarsHack
        # TurbulenceConvection.io(sim) # #removeVarsHack
        write_simulation_time(stats, t) # #removeVarsHack

        state = TC.column_prog_aux(prog, aux, colidx)
        grid = TC.Grid(state)
        diag_col = TC.column_diagnostics(diagnostics, colidx)

        # TODO: is this the best location to call diagnostics?
        compute_diagnostics!(edmf, precip_model, param_set, grid, state, diag_col, stats, surf_params, t, calibrate_io)

        cent = TC.Cent(1)
        diag_svpc = svpc_diagnostics_grid_mean(diag_col)
        diag_tc_svpc = svpc_diagnostics_turbconv(diag_col)
        write_ts(stats, "lwp_mean", diag_svpc.lwp_mean[cent])
        write_ts(stats, "iwp_mean", diag_svpc.iwp_mean[cent])
        write_ts(stats, "rwp_mean", diag_svpc.rwp_mean[cent])
        write_ts(stats, "swp_mean", diag_svpc.swp_mean[cent])

        write_ts(stats, "integ_total_flux_qt", diag_svpc.integ_total_flux_qt[cent])
        write_ts(stats, "integ_total_flux_s", diag_svpc.integ_total_flux_s[cent])

        if !calibrate_io
            write_ts(stats, "updraft_cloud_cover", diag_tc_svpc.updraft_cloud_cover[cent])
            write_ts(stats, "updraft_cloud_base", diag_tc_svpc.updraft_cloud_base[cent])
            write_ts(stats, "updraft_cloud_top", diag_tc_svpc.updraft_cloud_top[cent])
            write_ts(stats, "env_cloud_cover", diag_tc_svpc.env_cloud_cover[cent])
            write_ts(stats, "env_cloud_base", diag_tc_svpc.env_cloud_base[cent])
            write_ts(stats, "env_cloud_top", diag_tc_svpc.env_cloud_top[cent])
            write_ts(stats, "env_lwp", diag_tc_svpc.env_lwp[cent])
            write_ts(stats, "env_iwp", diag_tc_svpc.env_iwp[cent])
            write_ts(stats, "Hd", diag_tc_svpc.Hd[cent])
            write_ts(stats, "updraft_lwp", diag_tc_svpc.updraft_lwp[cent])
            write_ts(stats, "updraft_iwp", diag_tc_svpc.updraft_iwp[cent])

            write_ts(stats, "cutoff_precipitation_rate", diag_svpc.cutoff_precipitation_rate[cent])
            write_ts(stats, "cloud_cover_mean", diag_svpc.cloud_cover_mean[cent])
            write_ts(stats, "cloud_base_mean", diag_svpc.cloud_base_mean[cent])
            write_ts(stats, "cloud_top_mean", diag_svpc.cloud_top_mean[cent])
        end

        io(io_nt.aux, stats, state)
        io(io_nt.diagnostics, stats, diag_col)

        surf = get_surface(surf_params, grid, state, t, param_set)
        io(surf, surf_params, grid, state, stats, t)
        nothing
    end

    ODE.u_modified!(integrator, false) # We're legitamately not mutating `u` (the state vector)
end

function affect_filter!(integrator)
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, param_set, aux, case, surf_params) = integrator.p
        t = integrator.t
        prog = integrator.u

        state = TC.column_prog_aux(prog, aux, colidx)
        grid = TC.Grid(state)
        surf = get_surface(surf_params, grid, state, t, param_set)
        TC.affect_filter!(edmf, grid, state, param_set, surf, t)
        nothing
    end

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

function compute_dt_max(state::TC.State, edmf::TC.EDMFModel, dt_max::FT, CFL_limit::FT) where {FT <: Real}
    grid = TC.Grid(state)

    prog_gm = TC.center_prog_grid_mean(state)
    prog_gm_f = TC.face_prog_grid_mean(state)
    Δzc = TC.get_Δz(prog_gm.ρ)
    Δzf = TC.get_Δz(prog_gm_f.w)
    N_up = TC.n_updrafts(edmf)

    aux_tc = TC.center_aux_turbconv(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_en_f = TC.face_aux_environment(state)
    KM = aux_tc.KM
    KH = aux_tc.KH
    KQ = aux_tc.KQ

    # helper to calculate the rain velocity
    # TODO: assuming w_gm = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow
    ε = FT(eps(FT))

    @inbounds for k in TC.real_face_indices(grid)
        TC.is_surface_face(grid, k) && continue
        @inbounds for i in 1:N_up
            dt_max = min(dt_max, CFL_limit * Δzf[k] / (abs(aux_up_f[i].w[k]) + ε))
        end
        dt_max = min(dt_max, CFL_limit * Δzf[k] / (abs(aux_en_f.w[k]) + ε))
    end
    @inbounds for k in TC.real_center_indices(grid)
        vel_max = max(term_vel_rain[k], term_vel_snow[k])
        # Check terminal rain/snow velocity CFL
        dt_max = min(dt_max, CFL_limit * Δzc[k] / (vel_max + ε))
        # Check diffusion CFL (i.e., Fourier number)
        dt_max = min(dt_max, CFL_limit * Δzc[k]^2 / (max(KH[k], KM[k], KQ[k]) + ε))
    end
    return dt_max
end

dt_max!(integrator) = dt_max!(integrator, integrator.p.TS.dt_max)
function dt_max!(integrator, dt_max)
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; edmf, aux, TS) = integrator.p
        prog = integrator.u
        state = TC.column_prog_aux(prog, aux, colidx)
        dt_max = compute_dt_max(state, edmf, dt_max, TS.cfl_limit)
    end
    to_float(f) = f isa ForwardDiff.Dual ? ForwardDiff.value(f) : f
    (; TS) = integrator.p
    TS.dt_max_edmf = to_float(dt_max)

    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(integrator)
    CC.Fields.bycolumn(axes(integrator.u.cent)) do colidx
        (; aux, TS) = integrator.p
        prog = integrator.u
        state = TC.column_prog_aux(prog, aux, colidx)
        monitor_cfl!(state, TS.dt, TS.cfl_limit)
    end
    ODE.u_modified!(integrator, false)
end

function monitor_cfl!(state, Δt, CFL_limit)
    grid = TC.Grid(state)
    prog_gm = TC.center_prog_grid_mean(state)
    Δz = TC.get_Δz(prog_gm.ρ)
    FT = TC.float_type(state)
    frac_tol = FT(0.05)

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
            CFL_in_rain = FT(0)
            CFL_in_snow = FT(0)
        else
            CFL_in_rain = Δt / Δz[k] * term_vel_rain[k + 1]
            CFL_in_snow = Δt / Δz[k] * term_vel_snow[k + 1]
        end
        CFL_meteor, i_CFL_meteor = findmax((CFL_in_rain, CFL_in_snow, CFL_out_rain, CFL_out_snow))
        if CFL_meteor > CFL_limit * (1 + frac_tol)
            prog_pr = TC.center_prog_precipitation(state)
            prog_gm = TC.center_prog_grid_mean(state)
            aux_gm = TC.center_aux_grid_mean(state)
            aux_bulk = TC.center_aux_bulk(state)
            aux_en = TC.center_aux_environment(state)

            if i_CFL_meteor == 1
                term_vel = term_vel_rain[k + 1]
                species = "rain"
                q = prog_pr.q_rai[k + 1]
                q_liq, q_ice = prog_gm.q_liq[k + 1], prog_gm.q_ice[k + 1]
                q_tot = aux_gm.q_tot[k + 1]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k + 1], aux_bulk.q_liq[k + 1], aux_bulk.q_ice[k + 1]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k + 1], aux_en.q_liq[k + 1], aux_en.q_ice[k + 1]
                updraft_area = aux_bulk.area[k + 1]
                T = aux_gm.T[k + 1]
            elseif i_CFL_meteor == 2
                term_vel = term_vel_snow[k + 1]
                species = "snow"
                q = prog_pr.q_sno[k + 1]
                q_liq, q_ice = prog_gm.q_liq[k + 1], prog_gm.q_ice[k + 1]
                q_tot = aux_gm.q_tot[k + 1]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k + 1], aux_bulk.q_liq[k + 1], aux_bulk.q_ice[k + 1]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k + 1], aux_en.q_liq[k + 1], aux_en.q_ice[k + 1]
                updraft_area = aux_bulk.area[k + 1]
                T = aux_gm.T[k + 1]
            elseif i_CFL_meteor == 3
                term_vel = term_vel_rain[k]
                species = "rain"
                q = prog_pr.q_rai[k]
                q_liq, q_ice = prog_gm.q_liq[k], prog_gm.q_ice[k]
                q_tot = aux_gm.q_tot[k]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k], aux_bulk.q_liq[k], aux_bulk.q_ice[k]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k]
                updraft_area = aux_bulk.area[k]
                T = aux_gm.T[k]
            else
                term_vel = term_vel_snow[k]
                species = "snow"
                q = prog_pr.q_sno[k]
                q_liq, q_ice = prog_gm.q_liq[k], prog_gm.q_ice[k]
                q_tot = aux_gm.q_tot[k]
                q_tot_up, q_liq_up, q_ice_up = aux_bulk.q_tot[k], aux_bulk.q_liq[k], aux_bulk.q_ice[k]
                q_tot_env, q_liq_env, q_ice_env = aux_en.q_tot[k], aux_en.q_liq[k], aux_en.q_ice[k]
                updraft_area = aux_bulk.area[k]
                T = aux_gm.T[k]
            end
            q_tot
            error(
                "Time step is too large for hydrometeor fall velocity! Hydrometeor CFL is $(CFL_meteor) > $(CFL_limit), Δt = $(Δt), Δz = $(Δz[k]), term_vel = $(term_vel), species = $(species), q_$(species) = $(q), q_liq = $(q_liq), q_ice = $(q_ice), q_tot = $(q_tot) | q_tot_up = $(q_tot_up), q_liq_up = $(q_liq_up), q_ice_up = $(q_ice_up) | q_tot_env = $(q_tot_env), q_liq_env = $(q_liq_env), q_ice_env = $(q_ice_env) | updraft_area = $(updraft_area) | T = $(T)",
            )
        end
    end
    return nothing
end
