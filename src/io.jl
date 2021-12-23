#= This file is being deprecated in favor of using diagnostics.jl =#

initialize_io(self::ForcingBase, Stats) = nothing

initialize_io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing
io(self::RadiationBase, grid, state, Stats::NetCDFIO_Stats) = nothing

function initialize_io(en::EnvironmentVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "env_cloud_base")
    add_ts(Stats, "env_cloud_top")
    add_ts(Stats, "env_cloud_cover")
    add_ts(Stats, "env_lwp")
    add_ts(Stats, "env_iwp")
    return
end

function io(en::EnvironmentVariables, grid, state, Stats::NetCDFIO_Stats)
    # Assuming amximum overlap in environmental clouds
    write_ts(Stats, "env_cloud_cover", en.cloud_cover)
    write_ts(Stats, "env_cloud_base", en.cloud_base)
    write_ts(Stats, "env_cloud_top", en.cloud_top)
    return
end

function initialize_io(precip::PrecipVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "rwp_mean")
    add_ts(Stats, "swp_mean")
    add_ts(Stats, "cutoff_precipitation_rate")
    return
end

io(precip::PrecipVariables, grid, state, Stats::NetCDFIO_Stats) = nothing

function initialize_io(up::UpdraftVariables, Stats::NetCDFIO_Stats)

    add_ts(Stats, "updraft_cloud_cover")
    add_ts(Stats, "updraft_cloud_base")
    add_ts(Stats, "updraft_cloud_top")
    add_ts(Stats, "updraft_lwp")
    add_ts(Stats, "updraft_iwp")
    return
end

function io(up::UpdraftVariables, grid, state, Stats::NetCDFIO_Stats)
    # Note definition of cloud cover : each updraft is associated with a cloud cover equal to the maximum
    # area fraction of the updraft where ql > 0. Each updraft is assumed to have maximum overlap with respect to
    # itup (i.e. no consideration of tilting due to shear) while the updraft classes are assumed to have no overlap
    # at all. Thus total updraft cover is the sum of each updraft"s cover
    write_ts(Stats, "updraft_cloud_cover", sum(up.cloud_cover))
    write_ts(Stats, "updraft_cloud_base", minimum(abs.(up.cloud_base)))
    write_ts(Stats, "updraft_cloud_top", maximum(abs.(up.cloud_top)))
    return
end

function initialize_io(gm::GridMeanVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "lwp_mean")
    add_ts(Stats, "iwp_mean")
    add_ts(Stats, "cloud_base_mean")
    add_ts(Stats, "cloud_top_mean")
    add_ts(Stats, "cloud_cover_mean")
    return
end

function io(gm::GridMeanVariables, grid, state, Stats::NetCDFIO_Stats)
    write_ts(Stats, "cloud_cover_mean", gm.cloud_cover)
    write_ts(Stats, "cloud_base_mean", gm.cloud_base)
    write_ts(Stats, "cloud_top_mean", gm.cloud_top)
    return
end

# Initialize the IO pertaining to this class
function initialize_io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)

    initialize_io(edmf.UpdVar, Stats)
    initialize_io(edmf.EnvVar, Stats)
    initialize_io(edmf.Precip, Stats)

    add_ts(Stats, "rd")
    return
end

function io(edmf::EDMF_PrognosticTKE, grid, state, Stats::NetCDFIO_Stats)
    io(edmf.UpdVar, grid, state, Stats)
    io(edmf.EnvVar, grid, state, Stats)
    io(edmf.Precip, grid, state, Stats)
    write_ts(Stats, "rd", StatsBase.mean(edmf.pressure_plume_spacing))
    return
end
