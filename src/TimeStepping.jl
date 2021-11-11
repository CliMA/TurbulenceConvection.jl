mutable struct TimeStepping
    dt::Float64
    t_max::Float64
    t::Float64
    nstep::Int
end

function TimeStepping(namelist)
    dt = parse_namelist(namelist, "time_stepping", "dt_min"; default = 1.0)
    t_max = parse_namelist(namelist, "time_stepping", "t_max"; default = 7200.0)

    # set time
    t = 0.0
    nstep = 0

    return TimeStepping(dt, t_max, t, nstep)
end

function update(self)
    self.t += self.dt
    self.nstep += 1
    return
end
