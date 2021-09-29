#####
##### Diagnostics
#####

#=
    io_dictionary_ref_state()
    io_dictionary_aux()
    io_dictionary_state()
    io_dictionary_tendencies()

All of these functions return a dictionary whose
 - `keys` are the nc variable names
 - `values` are NamedTuples corresponding
    - `dims` (`("z")`  or `("z", "t")`) and
    - `group` (`"reference"` or `"profiles"`)
=#

function io_dictionary_ref_state(aux)
    DT = NamedTuple{(:dims, :group, :field), Tuple{Tuple{String}, String, Any}}
    cent_ref_state = center_ref_state # so that things nicely align :)
    io_dict = Dict{String, DT}(
        "ρ0_f" => (; dims = ("zf",), group = "reference", field = face_ref_state(aux).ρ0),
        "ρ0_c" => (; dims = ("zc",), group = "reference", field = cent_ref_state(aux).ρ0),
        "p0_f" => (; dims = ("zf",), group = "reference", field = face_ref_state(aux).p0),
        "p0_c" => (; dims = ("zc",), group = "reference", field = cent_ref_state(aux).p0),
        "α0_f" => (; dims = ("zf",), group = "reference", field = face_ref_state(aux).α0),
        "α0_c" => (; dims = ("zc",), group = "reference", field = cent_ref_state(aux).α0),
    )
    return io_dict
end
io_dictionary_aux(aux) = Dict()
io_dictionary_state(state) = Dict()
io_dictionary_tendencies(tendencies) = Dict()

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
