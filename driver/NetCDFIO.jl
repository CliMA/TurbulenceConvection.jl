import NCDatasets
const NC = NCDatasets
import JSON
import TurbulenceConvection
const TC = TurbulenceConvection

# TODO: remove `vars` hack that avoids https://github.com/Alexander-Barth/NCDatasets.jl/issues/135

function nc_fileinfo(namelist)

    uuid = string(namelist["meta"]["uuid"])
    simname = namelist["meta"]["simname"]
    outpath = joinpath(namelist["output"]["output_root"], "Output.$simname.$uuid")
    @info "Output folder: `$outpath`"
    mkpath(outpath)

    nc_filename = joinpath(outpath, namelist["stats_io"]["stats_dir"])
    mkpath(nc_filename)
    @info "NC filename path: `$nc_filename`"

    nc_filename = joinpath(nc_filename, "Stats.$simname.nc")
    @info "NC filename: `$nc_filename`"
    return nc_filename, outpath
end

mutable struct NetCDFIO_Stats
    root_grp::NC.NCDataset{Nothing}
    profiles_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    ts_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    frequency::Float64
    nc_filename::String
    vars::Dict{String, Any} # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
end

# Convenience backward compatible outer constructor
function NetCDFIO_Stats(namelist, grid::TC.Grid)
    frequency = namelist["stats_io"]["frequency"]
    # Setup the statistics output path
    casename = namelist["meta"]["casename"]
    nc_filename, outpath = nc_fileinfo(namelist)
    # Write namelist file to output directory
    open(joinpath(outpath, "namelist_$casename.in"), "w") do io
        JSON.print(io, namelist, 4)
    end
    NetCDFIO_Stats(; nc_filename, frequency, z_faces = vec(grid.zf), z_centers = vec(grid.zc))
end

function NetCDFIO_Stats(; nc_filename, frequency, z_faces, z_centers)
    # Initialize properties with valid type:
    tmp = tempname()
    root_grp = NC.Dataset(tmp, "c")
    NC.defGroup(root_grp, "profiles")
    NC.defGroup(root_grp, "timeseries")
    profiles_grp = root_grp.group["profiles"]
    ts_grp = root_grp.group["timeseries"]
    close(root_grp)

    # Remove the NC file if it exists, in case it accidentally wasn't closed
    isfile(nc_filename) && rm(nc_filename; force = true)

    NC.Dataset(nc_filename, "c") do root_grp
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
    end
    vars = Dict{String, Any}()
    return NetCDFIO_Stats(root_grp, profiles_grp, ts_grp, frequency, nc_filename, vars)
end


function open_files(self::NetCDFIO_Stats)
    self.root_grp = NC.Dataset(self.nc_filename, "a")
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

function add_field(ds, var_name::String, dims, group)
    profile_grp = ds.group[group]
    new_var = NC.defVar(profile_grp, var_name, Float64, dims)
    return nothing
end

#####
##### Time-series data
#####

function add_ts(ds, var_name::String)
    ts_grp = ds.group["timeseries"]
    new_var = NC.defVar(ts_grp, var_name, Float64, ("t",))
    return nothing
end

#####
##### Performance critical IO
#####

function write_field(self::NetCDFIO_Stats, var_name::String, data::T, group) where {T <: AbstractArray{Float64, 1}}
    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    @inbounds self.vars[group][var_name][:, end] = data
    # Ideally, we remove self.vars and use:
    # var = self.profiles_grp[var_name]
    # Not sure why `end` instead of `end+1`, but `end+1` produces garbage output
    # @inbounds var[end, :] = data :: T
end

function add_write_field(ds, var_name::String, data::T, group, dims) where {T <: AbstractArray{Float64, 1}}
    grp = ds.group[group]
    NC.defVar(grp, var_name, Float64, dims)
    var = grp[var_name]
    var .= data::T
    return nothing
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
