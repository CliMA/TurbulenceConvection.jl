#= This file is being deprecated in favor of using diagnostics.jl =#

initialize_io(self::ForcingBase, Stats) = nothing

initialize_io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing

function initialize_io(precip::AbstractPrecipitationModel, Stats::NetCDFIO_Stats)
    add_ts(Stats, "rwp_mean")
    add_ts(Stats, "swp_mean")
    add_ts(Stats, "cutoff_precipitation_rate")
    return nothing
end

function initialize_io(gm::GridMeanVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "lwp_mean")
    add_ts(Stats, "iwp_mean")
    add_ts(Stats, "cloud_base_mean")
    add_ts(Stats, "cloud_top_mean")
    add_ts(Stats, "cloud_cover_mean")
    return nothing
end

# Initialize the IO pertaining to this class
function initialize_io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)
    add_ts(Stats, "env_cloud_base")
    add_ts(Stats, "env_cloud_top")
    add_ts(Stats, "env_cloud_cover")
    add_ts(Stats, "env_lwp")
    add_ts(Stats, "env_iwp")
    add_ts(Stats, "updraft_cloud_cover")
    add_ts(Stats, "updraft_cloud_base")
    add_ts(Stats, "updraft_cloud_top")
    add_ts(Stats, "updraft_lwp")
    add_ts(Stats, "updraft_iwp")
    initialize_io(edmf.precip_model, Stats)
    add_ts(Stats, "rd")
    return nothing
end
