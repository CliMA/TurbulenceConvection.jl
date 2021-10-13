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
function io_dictionary_aux(state)
    DT = NamedTuple{(:dims, :group, :field), Tuple{Tuple{String, String}, String, Any}}
    io_dict = Dict{String, DT}(
        "updraft_area" => (; dims = ("zc", "t"), group = "profiles", field = center_aux_tc(state).bulk.area),
    )
    return io_dict
end
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
    kc_toa = kc_top_of_atmos(grid)
    gm.cloud_base = grid.zc[kc_toa]
    gm.cloud_top = 0.0

    @inbounds for k in real_center_indices(grid)
        gm.lwp += ρ0_c[k] * gm.QL.values[k] * grid.Δz
        if gm.QL.values[k] > 1e-8
            gm.cloud_base = min(gm.cloud_base, grid.zc[k])
            gm.cloud_top = max(gm.cloud_top, grid.zc[k])
        end
    end
    return
end
