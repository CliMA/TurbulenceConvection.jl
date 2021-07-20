using NCDatasets
using JSON

# TODO: remove `vars` hack that avoids https://github.com/Alexander-Barth/NCDatasets.jl/issues/135

mutable struct NetCDFIO_Stats
    root_grp::NCDatasets.NCDataset{Nothing}
    profiles_grp::NCDatasets.NCDataset{NCDatasets.NCDataset{Nothing}}
    ts_grp::NCDatasets.NCDataset{NCDatasets.NCDataset{Nothing}}
    Gr::Grid
    last_output_time::Float64
    uuid::String
    frequency::Float64
    stats_path::String
    path_plus_file::String
    vars::Dict{String, Any} # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    function NetCDFIO_Stats(namelist, paramlist, Gr::Grid)

        # Initialize properties with valid type:
        tmp = tempname()
        root_grp = Dataset(tmp, "c")
        defGroup(root_grp, "profiles")
        defGroup(root_grp, "timeseries")
        profiles_grp = root_grp.group["profiles"]
        ts_grp = root_grp.group["timeseries"]
        close(root_grp)

        last_output_time = 0.0
        uuid = string(namelist["meta"]["uuid"])

        frequency = namelist["stats_io"]["frequency"]

        # Setup the statistics output path
        simname = namelist["meta"]["simname"]
        casename = paramlist["meta"]["casename"]
        outpath = joinpath(namelist["output"]["output_root"], "Output.$simname.$uuid")
        mkpath(outpath)

        stats_path = joinpath(outpath, namelist["stats_io"]["stats_dir"])
        mkpath(stats_path)

        path_plus_file = joinpath(stats_path, "Stats.$simname.nc")

        # TODO: uncomment restart
        # if isfile(path_plus_file)
        #   @inbounds for i in 1:100
        #         res_name = "Restart_$i"
        #         if isfile(path_plus_file)
        #             path_plus_file = stats_path * "Stats.$simname.$res_name.nc"
        #         else
        #             break
        #         end
        #     end
        # end

        # Write namefile and paramfile to output directory
        open(joinpath(outpath, "$simname.in"), "w") do io
            JSON.print(io, namelist, 4)
        end
        open(joinpath(outpath, "paramlist_$casename.in"), "w") do io
            JSON.print(io, namelist, 4)
        end

        # TODO: make cell centers and cell faces different sizes
        cinterior = Gr.cinterior
        finterior = Gr.finterior

        Dataset(path_plus_file, "c") do root_grp

            zf = Gr.z[finterior]
            zc = Gr.z_half[cinterior]

            # Set profile dimensions
            profile_grp = defGroup(root_grp, "profiles")
            defDim(profile_grp, "z", Gr.nz)
            defDim(profile_grp, "t", Inf)
            defVar(profile_grp, "z", zf, ("z",))
            defVar(profile_grp, "z_half", zc, ("z_half",))
            defVar(profile_grp, "t", Float64, ("t",))

            reference_grp = defGroup(root_grp, "reference")
            defDim(reference_grp, "z", Gr.nz)
            defVar(reference_grp, "z", zf, ("z",))
            defVar(reference_grp, "z_half", zc, ("z_half",))

            ts_grp = defGroup(root_grp, "timeseries")
            defDim(ts_grp, "t", Inf)
            defVar(ts_grp, "t", Float64, ("t",))
        end
        vars = Dict{String, Any}()
        return new(
            root_grp,
            profiles_grp,
            ts_grp,
            Gr,
            last_output_time,
            uuid,
            frequency,
            stats_path,
            path_plus_file,
            vars,
        )
    end
end


function open_files(self)
    self.root_grp = Dataset(self.path_plus_file, "a")
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

function add_profile(self::NetCDFIO_Stats, var_name::String)
    Dataset(self.path_plus_file, "a") do root_grp
        profile_grp = root_grp.group["profiles"]
        new_var = defVar(profile_grp, var_name, Float64, ("z", "t"))
    end
end

function add_reference_profile(self::NetCDFIO_Stats, var_name::String)
    Dataset(self.path_plus_file, "a") do root_grp
        reference_grp = root_grp.group["reference"]
        new_var = defVar(reference_grp, var_name, Float64, ("z",))
    end
end

function add_ts(self::NetCDFIO_Stats, var_name::String)
    Dataset(self.path_plus_file, "a") do root_grp
        ts_grp = root_grp.group["timeseries"]
        new_var = defVar(ts_grp, var_name, Float64, ("t",))
    end
end

""" Writes a profile to the reference group NetCDF Stats file.
The variable must have already been added to the NetCDF file
using `add_reference_profile`.

Parameters
----------
var_name :: name of variables
data :: data to be written to file
"""
function write_reference_profile(self::NetCDFIO_Stats, var_name, data::T) where {T <: AbstractArray{Float64, 1}}

    Dataset(self.path_plus_file, "a") do root_grp
        reference_grp = root_grp.group["reference"]
        var = reference_grp[var_name]
        var .= data::T
    end
end

#####
##### Performance critical IO
#####

function write_profile(self::NetCDFIO_Stats, var_name::String, data::T) where {T <: AbstractArray{Float64, 1}}
    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    @inbounds self.vars["profiles"][var_name][:, end] = data
    # Ideally, we remove self.vars and use:
    # var = self.profiles_grp[var_name]
    # Not sure why `end` instead of `end+1`, but `end+1` produces garbage output
    # @inbounds var[end, :] = data :: T
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
