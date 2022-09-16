import NCDatasets as NC
import JSON
import TurbulenceConvection as TC

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

struct NetCDFIO_Stats{FT}
    frequency::FT
    nc_filename::String
end

# Convenience backward compatible outer constructor
function NetCDFIO_Stats(nc_filename, frequency, grid::TC.Grid)
    NetCDFIO_Stats(; nc_filename, frequency, z_faces = vec(grid.zf.z), z_centers = vec(grid.zc.z))
end

# Avoid quadratic lookup due to CFConventions:
# https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
get_var(group, var_name) = NC.cfvariable(group, var_name, _parentname = nothing)
# get_var(group, var_name) = group[var_name] # slow due to NCDatasets.jl/issues/135

function NetCDFIO_Stats(; nc_filename, frequency, z_faces, z_centers)
    FT = eltype(z_faces)
    # Remove the NC file if it exists, in case it accidentally wasn't closed
    isfile(nc_filename) && rm(nc_filename; force = true)

    NC.Dataset(nc_filename, "c") do ds
        # Set profile dimensions
        profile_grp = NC.defGroup(ds, "profiles")
        NC.defDim(profile_grp, "zf", length(z_faces))
        NC.defDim(profile_grp, "zc", length(z_centers))
        NC.defDim(profile_grp, "t", Inf)
        NC.defVar(profile_grp, "zf", z_faces, ("zf",))
        NC.defVar(profile_grp, "zc", z_centers, ("zc",))
        NC.defVar(profile_grp, "t", FT, ("t",))

        reference_grp = NC.defGroup(ds, "reference")
        NC.defDim(reference_grp, "zf", length(z_faces))
        NC.defDim(reference_grp, "zc", length(z_centers))
        NC.defVar(reference_grp, "zf", z_faces, ("zf",))
        NC.defVar(reference_grp, "zc", z_centers, ("zc",))

        ts_grp = NC.defGroup(ds, "timeseries")
        NC.defDim(ts_grp, "t", Inf)
        NC.defVar(ts_grp, "t", FT, ("t",))
    end
    return NetCDFIO_Stats{FT}(frequency, nc_filename)
end
