mutable struct TimeStepping{FT}
    dt::FT
    dt_min::FT
    t_max::FT
    t::FT
    nstep::Int
    cfl_limit::FT
    dt_max::FT
    dt_max_edmf::FT
    dt_io::FT
    spinup_half_t_max::FT # 1/2 of the spinup period -- we use Δt * spinup_dt_factor where spinup_dt_factor = spinup_dt_factor from 0 to spinup_half_t_max and then ramps up from spinup_dt_factor to 1 from t=spinup_half_t_max to t=2*spinup_half_t_max
    spinup_dt_factor::FT
    spinup_adapt_dt::Bool
end

function TimeStepping(::Type{FT}, namelist) where {FT}
    dt = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = FT(1.0))
    dt_min = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = FT(1.0)) # placeholder to keep storing it in case spinup needs it since it gets edited in callbacks
    t_max = TC.parse_namelist(namelist, "time_stepping", "t_max"; default = FT(7200.0))
    cfl_limit = TC.parse_namelist(namelist, "time_stepping", "cfl_limit"; default = FT(0.5))
    dt_max = TC.parse_namelist(namelist, "time_stepping", "dt_max"; default = FT(10.0))
    dt_max_edmf = FT(0)

    spinup_half_t_max =  TC.parse_namelist(namelist, "time_stepping", "spinup_half_t_max"; default = FT(0.0))
    spinup_dt_factor =  TC.parse_namelist(namelist, "time_stepping", "spinup_dt_factor"; default = FT(1.0))
    adapt_dt = TC.parse_namelist(namelist, "time_stepping", "adapt_dt"; default = false)
    spinup_adapt_dt = TC.parse_namelist(namelist, "time_stepping", "spinup_adapt_dt"; default = adapt_dt)

    @assert (FT(0) <= spinup_half_t_max) "spinup_half_t_max must be greater than or equal to 0, given value $(spinup_half_t_max) is invalid"
    @assert (FT(0) < spinup_dt_factor ) "spinup_dt_factor must be greater than 0, given value $(spinup_dt_factor) is invalid"

    # set time
    t = FT(0)
    dt_io = FT(0)
    nstep = 0

    return TimeStepping{FT}(dt, dt_min, t_max, t, nstep, cfl_limit, dt_max, dt_max_edmf, dt_io, spinup_half_t_max, spinup_dt_factor, spinup_adapt_dt)
end
