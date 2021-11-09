mutable struct TimeStepping
    dt::Float64
    dti::Float64
    t_max::Float64
    dt_min::Float64
    t::Float64
    nstep::Int
end

function TimeStepping(namelist)
    dt = parse_namelist(namelist, "time_stepping", "dt"; default = 1.0)
    dti = 1.0 / dt
    t_max = parse_namelist(namelist, "time_stepping", "t_max"; default = 7200.0)
    dt_min = parse_namelist(namelist, "time_stepping", "dt_min"; default = 7200.0)

    # set time
    t = 0.0
    nstep = 0

    return TimeStepping(dt, dti, t_max, dt_min, t, nstep)
end

function update(self)
    self.t += self.dt
    self.nstep += 1
    return
end
