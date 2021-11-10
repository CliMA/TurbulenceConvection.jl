# TODO: should this live in its own module?

""" Purely diagnostic fields for the host model """
diagnostics(state, fl) = getproperty(state, TC.field_loc(fl))

center_diagnostics_grid_mean(state) = diagnostics(state, TC.CentField())
center_diagnostics_turbconv(state) = diagnostics(state, TC.CentField()).turbconv
face_diagnostics_turbconv(state) = diagnostics(state, TC.FaceField()).turbconv

#=
    io_dictionary_diagnostics()

These functions return a dictionary whose
 - `keys` are the nc variable names
 - `values` are NamedTuples corresponding to
    - `dims` (`("z")`  or `("z", "t")`) and
    - `group` (`"reference"` or `"profiles"`)

This dictionary is for purely diagnostic quantities--which
are not required to compute in order to run a simulation.
=#

#! format: off
function io_dictionary_diagnostics()
    DT = NamedTuple{(:dims, :group, :field), Tuple{Tuple{String, String}, String, Any}}
    io_dict = Dict{String, DT}(
        "nh_pressure" => (; dims = ("zf", "t"), group = "profiles", field = state -> face_diagnostics_turbconv(state).nh_pressure),
        "nh_pressure_adv" => (; dims = ("zf", "t"), group = "profiles", field = state -> face_diagnostics_turbconv(state).nh_pressure_adv,),
        "nh_pressure_drag" => (; dims = ("zf", "t"), group = "profiles", field = state -> face_diagnostics_turbconv(state).nh_pressure_drag,),
        "nh_pressure_b" => (; dims = ("zf", "t"), group = "profiles", field = state -> face_diagnostics_turbconv(state).nh_pressure_b,),
        "turbulent_entrainment" => (; dims = ("zc", "t"), group = "profiles", field = state -> center_diagnostics_turbconv(state).frac_turb_entr,),
        "horiz_K_eddy" => (; dims = ("zc", "t"), group = "profiles", field = state -> center_diagnostics_turbconv(state).horiz_K_eddy,),
        "entrainment_sc" => (; dims = ("zc", "t"), group = "profiles", field = state -> center_diagnostics_turbconv(state).entr_sc),
        "detrainment_sc" => (; dims = ("zc", "t"), group = "profiles", field = state -> center_diagnostics_turbconv(state).detr_sc),
        "asp_ratio" => (; dims = ("zc", "t"), group = "profiles", field = state -> center_diagnostics_turbconv(state).asp_ratio),
        "massflux" => (; dims = ("zc", "t"), group = "profiles", field = state -> center_diagnostics_turbconv(state).massflux),
    )
    return io_dict
end
#! format: on

#=
    compute_diagnostics!

Computes diagnostic quantities. The state _should not_ depend
on any quantities here. I.e., we should be able to shut down
diagnostics and still run, at which point we'll be able to export
the state, auxiliary fields (which the state does depend on), and
tendencies.
=#
function compute_diagnostics!(edmf, gm, grid, state, diagnostics, Case, TS)
    gm.lwp = 0.0
    gm.iwp = 0.0
    ρ0_c = TC.center_ref_state(state).ρ0
    p0_c = TC.center_ref_state(state).p0
    aux_gm = TC.center_aux_grid_mean(state)
    aux_en = TC.center_aux_environment(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_tc_f = TC.face_aux_turbconv(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_pr = TC.center_prog_precipitation(state)
    aux_bulk = TC.center_aux_bulk(state)
    a_up_bulk = aux_bulk.area
    kc_toa = TC.kc_top_of_atmos(grid)
    gm.cloud_base = grid.zc[kc_toa]
    gm.cloud_top = 0.0
    param_set = TC.parameter_set(gm)
    prog_gm = TC.center_prog_grid_mean(state)
    up = edmf.UpdVar
    en = edmf.EnvVar
    precip = edmf.Precip
    en_thermo = edmf.EnvThermo
    n_updrafts = up.n_updrafts
    diag_tc = center_diagnostics_turbconv(diagnostics)
    diag_tc_f = face_diagnostics_turbconv(diagnostics)

    @inbounds for k in TC.real_center_indices(grid)
        gm.lwp += ρ0_c[k] * aux_gm.q_liq[k] * grid.Δz
        gm.iwp += ρ0_c[k] * aux_gm.q_ice[k] * grid.Δz
        if TD.has_condensate(aux_gm.q_liq[k] + aux_gm.q_ice[k])
            gm.cloud_base = min(gm.cloud_base, grid.zc[k])
            gm.cloud_top = max(gm.cloud_top, grid.zc[k])
        end
        ts = TC.thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        aux_gm.s[k] = TD.specific_entropy(ts)
        ts_en = TC.thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])
        aux_en.s[k] = TD.specific_entropy(ts_en)
        @inbounds for i in 1:n_updrafts
            if aux_up[i].area[k] > 0.0
                ts_up = TC.thermo_state_pθq(param_set, p0_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k])
                aux_up[i].s[k] = TD.specific_entropy(ts_up)
            end
        end
    end
    m_bcs = (; bottom = TC.SetValue(0), top = TC.SetValue(0))
    @inbounds for k in TC.real_face_indices(grid)
        s_dual = TC.dual_centers(aux_en.s, grid, k)
        # TODO(ilopezgp): Fix bottom gradient
        ∇s_f = TC.∇c2f(s_dual, grid, k; bottom = TC.SetGradient(0), top = TC.SetGradient(0))
        aux_gm_f.diffusive_flux_s[k] = -aux_tc_f.ρ_ae_KH[k] * ∇s_f
        s_en_f = TC.interpc2f(aux_en.s, grid, k; m_bcs...)
        @inbounds for i in 1:n_updrafts
            s_up_f = TC.interpc2f(aux_up[i].s, grid, k; m_bcs...)
            aux_gm_f.massflux_s[k] += aux_up_f[i].massflux[k] * (s_up_f - s_en_f)
        end
    end

    up.lwp = 0.0
    up.iwp = 0.0

    @inbounds for i in 1:(up.n_updrafts)
        up.cloud_base[i] = TC.zc_toa(grid)
        up.cloud_top[i] = 0.0
        up.updraft_top[i] = 0.0
        up.cloud_cover[i] = 0.0

        @inbounds for k in TC.real_center_indices(grid)
            if aux_up[i].area[k] > 1e-3
                up.updraft_top[i] = max(up.updraft_top[i], grid.zc[k])
                up.lwp += ρ0_c[k] * aux_up[i].q_liq[k] * aux_up[i].area[k] * grid.Δz
                up.iwp += ρ0_c[k] * aux_up[i].q_ice[k] * aux_up[i].area[k] * grid.Δz

                if TD.has_condensate(aux_up[i].q_liq[k] + aux_up[i].q_ice[k])
                    up.cloud_base[i] = min(up.cloud_base[i], grid.zc[k])
                    up.cloud_top[i] = max(up.cloud_top[i], grid.zc[k])
                    up.cloud_cover[i] = max(up.cloud_cover[i], aux_up[i].area[k])
                end
            end
        end
    end

    en.cloud_top = 0.0
    en.cloud_base = TC.zc_toa(grid)
    en.cloud_cover = 0.0
    en.lwp = 0.0
    en.iwp = 0.0

    @inbounds for k in TC.real_center_indices(grid)
        en.lwp += ρ0_c[k] * aux_en.q_liq[k] * aux_en.area[k] * grid.Δz
        en.iwp += ρ0_c[k] * aux_en.q_ice[k] * aux_en.area[k] * grid.Δz

        if TD.has_condensate(aux_en.q_liq[k] + aux_en.q_ice[k]) && aux_en.area[k] > 1e-6
            en.cloud_base = min(en.cloud_base, grid.zc[k])
            en.cloud_top = max(en.cloud_top, grid.zc[k])
            en.cloud_cover = max(en.cloud_cover, aux_en.area[k] * aux_en.cloud_fraction[k])
        end
    end

    precip.mean_rwp = 0.0
    precip.mean_swp = 0.0
    precip.cutoff_precipitation_rate = 0.0

    @inbounds for k in TC.real_center_indices(grid)
        precip.mean_rwp += ρ0_c[k] * prog_pr.q_rai[k] * grid.Δz
        precip.mean_swp += ρ0_c[k] * prog_pr.q_sno[k] * grid.Δz

        # precipitation rate from cutoff microphysics scheme defined as a total amount of removed water
        # per timestep per EDMF surface area [mm/h]
        if (precip.precipitation_model == "cutoff")
            precip.cutoff_precipitation_rate -=
                (aux_en.qt_tendency_precip_formation[k] + aux_bulk.qt_tendency_precip_formation[k]) *
                ρ0_c[k] *
                grid.Δz / TC.rho_cloud_liq *
                3.6 *
                1e6
        end
    end

    @inbounds for k in TC.real_center_indices(grid)
        a_up_bulk_k = a_up_bulk[k]
        if a_up_bulk_k > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                aux_up_i = aux_up[i]
                diag_tc.massflux[k] += TC.interpf2c(aux_up_f[i].massflux, grid, k)
                diag_tc.entr_sc[k] += aux_up_i.area[k] * aux_up_i.entr_sc[k] / a_up_bulk_k
                diag_tc.detr_sc[k] += aux_up_i.area[k] * aux_up_i.detr_sc[k] / a_up_bulk_k
                diag_tc.asp_ratio[k] += aux_up_i.area[k] * aux_up_i.asp_ratio[k] / a_up_bulk_k
                diag_tc.frac_turb_entr[k] += aux_up_i.area[k] * aux_up_i.frac_turb_entr[k] / a_up_bulk_k
                diag_tc.horiz_K_eddy[k] += aux_up_i.area[k] * aux_up_i.horiz_K_eddy[k] / a_up_bulk_k
            end
        end
    end

    @inbounds for k in TC.real_face_indices(grid)
        a_up_bulk_f = TC.interpc2f(
            a_up_bulk,
            grid,
            k;
            bottom = TC.SetValue(sum(edmf.area_surface_bc)),
            top = TC.SetZeroGradient(),
        )
        if a_up_bulk_f > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                a_up_f = TC.interpc2f(
                    aux_up[i].area,
                    grid,
                    k;
                    bottom = TC.SetValue(edmf.area_surface_bc[i]),
                    top = TC.SetZeroGradient(),
                )
                diag_tc_f.nh_pressure[k] += a_up_f * aux_up_f[i].nh_pressure[k] / a_up_bulk_f
                diag_tc_f.nh_pressure_b[k] += a_up_f * aux_up_f[i].nh_pressure_b[k] / a_up_bulk_f
                diag_tc_f.nh_pressure_adv[k] += a_up_f * aux_up_f[i].nh_pressure_adv[k] / a_up_bulk_f
                diag_tc_f.nh_pressure_drag[k] += a_up_f * aux_up_f[i].nh_pressure_drag[k] / a_up_bulk_f
            end
        end
    end


    return
end
