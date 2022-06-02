mutable struct TimeStepping
    dt::Float64
    t_max::Float64
    t::Float64
    nstep::Int
    cfl_limit::Float64
    dt_max::Float64
    dt_max_edmf::Float64
    dt_io::Float64
    initialized::Bool
end

function TimeStepping(namelist)
    dt = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = 1.0)
    t_max = TC.parse_namelist(namelist, "time_stepping", "t_max"; default = 7200.0)
    cfl_limit = TC.parse_namelist(namelist, "time_stepping", "cfl_limit"; default = 0.5)
    dt_max = TC.parse_namelist(namelist, "time_stepping", "dt_max"; default = 10.0)
    dt_max_edmf = 0.0
    initialized = false

    # set time
    t = 0.0
    dt_io = 0.0
    nstep = 0

    return TimeStepping(dt, t_max, t, nstep, cfl_limit, dt_max, dt_max_edmf, dt_io, initialized)
end
