import NCDatasets
const NC = NCDatasets
import JSON
import TurbulenceConvection
const TC = TurbulenceConvection

# TODO: remove `vars` hack that avoids https://github.com/Alexander-Barth/NCDatasets.jl/issues/135

mutable struct NetCDFIO_Stats
    root_grp::NC.NCDataset{Nothing}
    profiles_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    ts_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    uuid::String
    frequency::Float64
    stats_path::String
    path_plus_file::String
    vars::Dict{String, Any} # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    function NetCDFIO_Stats(namelist, z_centers::AbstractVector, z_faces::AbstractVector)

        # Initialize properties with valid type:
        uuid = string(namelist["meta"]["uuid"])
        frequency = namelist["stats_io"]["frequency"]
        # Setup the statistics output path
        simname = namelist["meta"]["simname"]
        casename = namelist["meta"]["casename"]
        outpath = joinpath(namelist["output"]["output_root"], "Output.$simname.$uuid")
        stats_path = joinpath(outpath, namelist["stats_io"]["stats_dir"])

        mkpath(outpath)
        mkpath(stats_path)

        path_plus_file = joinpath(stats_path, "Stats.$simname.nc")

        # Write namelist file to output directory
        open(joinpath(outpath, "namelist_$casename.in"), "w") do io
            JSON.print(io, namelist, 4)
        end

        # Remove the NC file if it exists, in case it accidentally wasn't closed
        isfile(path_plus_file) && rm(path_plus_file; force = true)

        root_grp = NC.Dataset(path_plus_file, "c") do root_grp
            # Set profile dimensions
            profile_grp = NC.defGroup(root_grp, "profiles")
            NC.defDim(profile_grp, "zf", length(z_faces))
            NC.defDim(profile_grp, "zc", length(z_centers))
            NC.defDim(profile_grp, "t", Inf)
            NC.defVar(profile_grp, "zf", z_faces, ("zf",))
            NC.defVar(profile_grp, "zc", z_centers, ("zc",))
            NC.defVar(profile_grp, "t", Float64, ("t",))

            reference_grp = NC.defGroup(root_grp, "reference")
            NC.defDim(reference_grp, "zf", length(z_faces))
            NC.defDim(reference_grp, "zc", length(z_centers))
            NC.defVar(reference_grp, "zf", z_faces, ("zf",))
            NC.defVar(reference_grp, "zc", z_centers, ("zc",))

            ts_grp = NC.defGroup(root_grp, "timeseries")
            NC.defDim(ts_grp, "t", Inf)
            NC.defVar(ts_grp, "t", Float64, ("t",))
            root_grp
        end
        vars = Dict{String, Any}()
        return new(root_grp, profiles_grp, ts_grp, uuid, frequency, stats_path, path_plus_file, vars)
    end
end

NetCDFIO_Stats(namelist, grid::TC.Grid) = NetCDFIO_Stats(namelist, vec(grid.zc), vec(grid.zf))

function open_files(self::NetCDFIO_Stats)
    self.root_grp = NC.Dataset(self.path_plus_file, "a")
    self.profiles_grp = self.root_grp.group["profiles"]
    self.ts_grp = self.root_grp.group["timeseries"]
    vars = self.vars

    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    vars["profiles"] = Dict{String, Any}()
    for k in keys(self.profiles_grp)
        vars["profiles"][k] = self.profiles_grp[k]
    end
    vars["timeseries"] = Dict{String, Any}()
    for k in keys(self.ts_grp)
        vars["timeseries"][k] = self.ts_grp[k]
    end
end

function close_files(self::NetCDFIO_Stats)
    close(self.root_grp)
end

#####
##### Generic field
#####

function add_field(self::NetCDFIO_Stats, var_name::String; dims, group)
    NC.Dataset(self.path_plus_file, "a") do root_grp
        profile_grp = root_grp.group[group]
        new_var = NC.defVar(profile_grp, var_name, Float64, dims)
    end
end

#####
##### Time-series data
#####

function add_ts(self::NetCDFIO_Stats, var_name::String)
    NC.Dataset(self.path_plus_file, "a") do root_grp
        ts_grp = root_grp.group["timeseries"]
        new_var = NC.defVar(ts_grp, var_name, Float64, ("t",))
    end
end

#####
##### Performance critical IO
#####

# Field wrapper
write_field(self::NetCDFIO_Stats, var_name::String, data; group) = write_field(self, var_name, vec(data); group = group)

function write_field(self::NetCDFIO_Stats, var_name::String, data::T; group) where {T <: AbstractArray{Float64, 1}}
    if group == "profiles"
        # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
        @inbounds self.vars[group][var_name][:, end] = data
        # Ideally, we remove self.vars and use:
        # var = self.profiles_grp[var_name]
        # Not sure why `end` instead of `end+1`, but `end+1` produces garbage output
        # @inbounds var[end, :] = data :: T
    elseif group == "reference"
        NC.Dataset(self.path_plus_file, "a") do root_grp
            reference_grp = root_grp.group[group]
            var = reference_grp[var_name]
            var .= data::T
        end
    else
        error("Bad group given")
    end
end

function write_ts(self::NetCDFIO_Stats, var_name::String, data::Float64)
    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    @inbounds self.vars["timeseries"][var_name][end] = data::Float64
    # Ideally, we remove self.vars and use:
    # var = self.ts_grp[var_name]
    # @inbounds var[end+1] = data :: Float64
end

function write_simulation_time(self::NetCDFIO_Stats, t::Float64)
    # # Write to profiles group
    profile_t = self.profiles_grp["t"]
    @inbounds profile_t[end + 1] = t::Float64

    # # Write to timeseries group
    ts_t = self.ts_grp["t"]
    @inbounds ts_t[end + 1] = t::Float64
end
