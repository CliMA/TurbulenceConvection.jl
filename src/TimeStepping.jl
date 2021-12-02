mutable struct TimeStepping
    dt::Float64
    t_max::Float64
    t::Float64
    nstep::Int
    cfl_limit::Float64
    dt_max::Float64
    dt_io::Float64
    adapt_dt::Bool
end

function TimeStepping(namelist)
    dt = parse_namelist(namelist, "time_stepping", "dt_min"; default = 1.0)
    t_max = parse_namelist(namelist, "time_stepping", "t_max"; default = 7200.0)
    cfl_limit = parse_namelist(namelist, "time_stepping", "cfl_limit"; default = 0.5)
    dt_max = parse_namelist(namelist, "time_stepping", "dt_max"; default = 10.0)

    adapt_dt = parse_namelist(namelist, "time_stepping", "adapt_dt"; default = true)

    # set time
    t = 0.0
    dt_io = 0.0
    nstep = 0

    return TimeStepping(dt, t_max, t, nstep, cfl_limit, dt_max, dt_io, adapt_dt)
end

function update(self)
    self.t += self.dt
    self.nstep += 1
    return
end
