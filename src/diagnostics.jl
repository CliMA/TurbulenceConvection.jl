#####
##### Diagnostics
#####

#=
    io_dictionary_ref_state()
    io_dictionary_aux()
    io_dictionary_diagnostics()
    io_dictionary_state()
    io_dictionary_tendencies()

All of these functions return a dictionary whose
 - `keys` are the nc variable names
 - `values` are NamedTuples corresponding
    - `dims` (`("z")`  or `("z", "t")`) and
    - `group` (`"reference"` or `"profiles"`)
=#

function io_dictionary_ref_state(state)
    DT = NamedTuple{(:dims, :group, :field), Tuple{Tuple{String}, String, Any}}
    cent_ref_state = center_ref_state # so that things nicely align :)
    io_dict = Dict{String, DT}(
        "ρ0_f" => (; dims = ("zf",), group = "reference", field = face_ref_state(state).ρ0),
        "ρ0_c" => (; dims = ("zc",), group = "reference", field = cent_ref_state(state).ρ0),
        "p0_f" => (; dims = ("zf",), group = "reference", field = face_ref_state(state).p0),
        "p0_c" => (; dims = ("zc",), group = "reference", field = cent_ref_state(state).p0),
        "α0_f" => (; dims = ("zf",), group = "reference", field = face_ref_state(state).α0),
        "α0_c" => (; dims = ("zc",), group = "reference", field = cent_ref_state(state).α0),
    )
    return io_dict
end

#! format: off
# TODO: We probably don't need to split the aux/prog dictionaries. Only static vs dynamic.
function io_dictionary_aux(state)
    DT = NamedTuple{(:dims, :group, :field), Tuple{Tuple{String, String}, String, Any}}
    io_dict = Dict{String, DT}(
        "updraft_area" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.area),
        "updraft_ql" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.q_liq),
        "updraft_RH" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.RH),
        "updraft_qt" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.q_tot),
        "updraft_w" => (; dims = ("zf", "t"), group = "profiles", field = face_aux_tc(state).bulk.w),
        "updraft_temperature" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.T),
        "updraft_thetal" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.θ_liq_ice),
        "updraft_buoyancy" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.buoy),
        "H_third_m" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).H_third_m),
        "W_third_m" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).W_third_m),
        "QT_third_m" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).H_third_m),
        "cloud_fraction" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).cloud_fraction), # was this "cloud_fraction_mean"?
        "buoyancy_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).buoy),
        "temperature_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).T),
        "RH_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).RH),
        "s_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).s),
        "ql_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).q_liq),
        "qi_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).q_ice),
        "tke_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).tke),
        "Hvar_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).Hvar),
        "QTvar_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).QTvar),
        "HQTcov_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).HQTcov),
        "u_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).u),
        "v_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).v),
        "w_mean" => (; dims = ("zf", "t"), group = "profiles", field = face_prog_grid_mean(state).w),
        "qt_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).q_tot),
        "thetal_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).θ_liq_ice),
        "eddy_viscosity" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).KM),
        "eddy_diffusivity" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).KH),
        "env_tke" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_environment(state).tke),
        "env_Hvar" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_environment(state).Hvar),
        "env_QTvar" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_environment(state).QTvar),

        "env_HQTcov" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_environment(state).HQTcov),

        "tke_buoy" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.buoy),
        "tke_pressure" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.press),
        "tke_dissipation" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.dissipation),
        "tke_entr_gain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.entr_gain),
        "tke_detr_loss" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.detr_loss),
        "tke_shear" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.shear),
        "tke_interdomain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).tke.interdomain),

        "Hvar_dissipation" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).Hvar.dissipation),
        "Hvar_entr_gain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).Hvar.entr_gain),
        "Hvar_detr_loss" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).Hvar.detr_loss),
        "Hvar_interdomain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).Hvar.interdomain),
        "Hvar_rain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).Hvar.rain_src),
        "Hvar_shear" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).Hvar.shear),

        "QTvar_dissipation" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).QTvar.dissipation),
        "QTvar_entr_gain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).QTvar.entr_gain),
        "QTvar_detr_loss" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).QTvar.detr_loss),
        "QTvar_shear" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).QTvar.shear),
        "QTvar_rain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).QTvar.rain_src),
        "QTvar_interdomain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).QTvar.interdomain),

        "HQTcov_rain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).HQTcov.rain_src),
        "HQTcov_dissipation" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).HQTcov.dissipation),
        "HQTcov_entr_gain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).HQTcov.entr_gain),
        "HQTcov_detr_loss" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).HQTcov.detr_loss),
        "HQTcov_shear" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).HQTcov.shear),
        "HQTcov_interdomain" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment_2m(state).HQTcov.interdomain),

        "env_w" => (; dims = ("zf", "t"), group = "profiles", field = face_aux_environment(state).w),
        "env_qt" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).q_tot),
        "env_ql" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).q_liq),
        "env_area" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).area),
        "env_temperature" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).T),
        "env_RH" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).RH),
        "env_thetal" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).θ_liq_ice),
        "env_cloud_fraction" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_environment(state).cloud_fraction),
        "massflux_s" => (; dims = ("zf", "t"), group = "profiles", field = face_aux_grid_mean(state).massflux_s),
        "diffusive_flux_s" => (; dims = ("zf", "t"), group = "profiles", field = face_aux_grid_mean(state).diffusive_flux_s),
        "total_flux_s" => (; dims = ("zf", "t"), group = "profiles", field = face_aux_grid_mean(state).massflux_s .+ face_aux_grid_mean(state).diffusive_flux_s),

        "qr_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_rain(state).qr),

    )
    return io_dict
end
#! format: on
io_dictionary_diagnostics(state) = Dict()
io_dictionary_state(state) = Dict()
io_dictionary_tendencies(state) = Dict()

function initialize_io(io_dict::Dict, Stats::NetCDFIO_Stats)
    for var_name in keys(io_dict)
        add_field(Stats, var_name; dims = io_dict[var_name].dims, group = io_dict[var_name].group)
    end
end

function io(io_dict::Dict, Stats::NetCDFIO_Stats)
    for var in keys(io_dict)
        write_field(Stats, var, io_dict[var].field; group = io_dict[var].group)
    end
end

#=
    compute_diagnostics!

Computes diagnostic quantities. The state _should not_ depend
on any quantities here. I.e., we should be able to shut down
diagnostics and still run, at which point we'll be able to export
the state, auxiliary fields (which the state does depend on), and
tendencies.
=#
function compute_diagnostics!(edmf, gm, grid, state, Case, TS)
    gm.lwp = 0.0
    gm.iwp = 0.0
    ρ0_c = center_ref_state(state).ρ0
    p0_c = center_ref_state(state).p0
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_tc_f = face_aux_tc(state)
    aux_gm_f = face_aux_grid_mean(state)
    kc_toa = kc_top_of_atmos(grid)
    gm.cloud_base = grid.zc[kc_toa]
    gm.cloud_top = 0.0
    param_set = parameter_set(gm)
    prog_gm = center_prog_grid_mean(state)
    n_updrafts = edmf.UpdVar.n_updrafts

    @inbounds for k in real_center_indices(grid)
        gm.lwp += ρ0_c[k] * aux_gm.q_liq[k] * grid.Δz
        gm.iwp += ρ0_c[k] * aux_gm.q_ice[k] * grid.Δz
        if TD.has_condensate(aux_gm.q_liq[k] + aux_gm.q_ice[k])
            gm.cloud_base = min(gm.cloud_base, grid.zc[k])
            gm.cloud_top = max(gm.cloud_top, grid.zc[k])
        end
        ts = thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        aux_gm.s[k] = TD.specific_entropy(ts)
        ts_en = thermo_state_pθq(param_set, p0_c[k], aux_en.θ_liq_ice[k], aux_en.q_tot[k])
        aux_en.s[k] = TD.specific_entropy(ts_en)
        @inbounds for i in 1:n_updrafts
            if aux_up[i].area[k] > 0.0
                ts_up = thermo_state_pθq(param_set, p0_c[k], aux_up[i].θ_liq_ice[k], aux_up[i].q_tot[k])
                aux_up[i].s[k] = TD.specific_entropy(ts_up)
            end
        end
    end
    m_bcs = (; bottom = SetValue(0), top = SetValue(0))
    @inbounds for k in real_face_indices(grid)
        s_dual = dual_centers(aux_en.s, grid, k)
        # TODO(ilopezgp): Fix bottom gradient
        ∇s_f = ∇c2f(s_dual, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        aux_gm_f.diffusive_flux_s[k] = -aux_tc_f.ρ_ae_KH[k] * ∇s_f
        s_en_f = interpc2f(aux_en.s, grid, k; m_bcs...)
        @inbounds for i in 1:n_updrafts
            s_up_f = interpc2f(aux_up[i].s, grid, k; m_bcs...)
            aux_gm_f.massflux_s[k] += edmf.m[i, k] * (s_up_f - s_en_f)
        end
    end
    return
end
