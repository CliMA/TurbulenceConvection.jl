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
        "ql_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).q_liq),
        "tke_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).tke),
        "Hvar_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).Hvar),
        "QTvar_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).QTvar),
        "HQTcov_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_grid_mean(state).HQTcov),
        "u_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).u),
        "v_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).v),
        "w_mean" => (; dims = ("zf", "t"), group = "profiles", field = face_prog_grid_mean(state).w),
        "qt_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).q_tot),
        "thetal_mean" => (; dims = ("zc", "t"), group = "profiles", field = center_prog_grid_mean(state).θ_liq_ice),
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
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    kc_toa = kc_top_of_atmos(grid)
    gm.cloud_base = grid.zc[kc_toa]
    gm.cloud_top = 0.0

    @inbounds for k in real_center_indices(grid)
        gm.lwp += ρ0_c[k] * aux_gm.q_liq[k] * grid.Δz
        if aux_gm.q_liq[k] > 1e-8
            gm.cloud_base = min(gm.cloud_base, grid.zc[k])
            gm.cloud_top = max(gm.cloud_top, grid.zc[k])
        end
    end
    return
end
