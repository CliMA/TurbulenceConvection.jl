mutable struct TimeStepping
    dt::Float64
    dti::Float64
    t_max::Float64
    t::Float64
    nstep::Int
    i_iter::Int
end

function TimeStepping(namelist)
    dt = try
        namelist["time_stepping"]["dt"]
    catch
        1.0
    end

    dti = 1.0/dt

    t_max = try
        namelist["time_stepping"]["t_max"]
    catch
        7200.0
    end

    # set time
    t = 0.0
    nstep = 0
    i_iter = 0

    return TimeStepping(dt,dti,t_max,t,nstep, i_iter)
end

function update(self)
    self.t += self.dt
    self.nstep += 1
    return
end

function update_iter(self)
    self.i_iter += 1
    return
end
