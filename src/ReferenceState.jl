#### RefState

using OrdinaryDiffEq
using UnPack

mutable struct ReferenceState{PS, A1}
    param_set::PS
    p0::A1
    p0_half::A1
    alpha0::A1
    alpha0_half::A1
    rho0::A1
    rho0_half::A1
    Pg::Float64
    Tg::Float64
    qtg::Float64
    sg::Float64
    function ReferenceState(grid::Grid, param_set::PS) where {PS}

        p0 = face_field(grid)
        p0_half = center_field(grid)
        alpha0 = face_field(grid)
        alpha0_half = center_field(grid)
        rho0 = face_field(grid)
        rho0_half = center_field(grid)
        Pg::Float64 = 0
        Tg::Float64 = 0
        qtg::Float64 = 0
        sg::Float64 = 0
        return new{PS, typeof(p0)}(param_set, p0, p0_half, alpha0, alpha0_half, rho0, rho0_half, Pg, Tg, qtg, sg)
    end
end


"""
Initialize the reference profiles. The function is typically called from the case
specific initialization function defined in Initialization.pyx

:param grid: Grid class
:param Thermodynamics: Thermodynamics class
:param NS: StatsIO class
:param Pa:  ParallelMPI class
:return
"""
function initialize(self::ReferenceState, grid::Grid, Stats::NetCDFIO_Stats)

    self.sg = t_to_entropy_c(self.Pg, self.Tg, self.qtg, 0.0, 0.0)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure

    function rhs(p, u, z)
        ret = eos(exp(p), self.qtg, self.sg; t_to_prog = t_to_entropy_c, prog_to_t = eos_first_guess_entropy)
        q_i = 0.0
        q_l = ret.ql
        T = ret.T
        return -g / (Rd * T * (1.0 - self.qtg + eps_vi * (self.qtg - q_l - q_i)))
    end

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    p0 = log(self.Pg)
    logp = p0

    p = face_field(grid)
    p_half = center_field(grid)

    # Perform the integration
    # TODO: replace with OrdinaryDiffEq

    z_span = (grid.zmin, grid.zmax)
    @show z_span
    prob = ODEProblem(rhs, logp, z_span)
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
    cinterior = kc_surface(grid):kc_top_of_atmos(grid)
    finterior = kf_surface(grid):(kf_top_of_atmos(grid) + 1) # TODO: this should not have +1
    p_0 = [sol(grid.z[k]) for k in finterior]
    p_0_half = [sol(grid.z_half[k]) for k in cinterior]

    p[finterior] .= p_0
    p_half[cinterior] .= p_0_half

    # Set boundary conditions (in log-space) by mirroring log-pressure
    # TODO: re-generalize, is setting the BCs like this correct?
    p[1] = p[2]
    p[end] = p[end - 1]

    p_half[2] = p_half[3]
    p_half[end - 1] = p_half[end - 2]
    p_half[1] = p_half[4]
    p_half[end] = p_half[end - 3]

    p .= exp.(p)
    p_half .= exp.(p_half)

    p_ = deepcopy(p)
    p_half_ = deepcopy(p_half)
    temperature = face_field(grid)
    temperature_half = center_field(grid)
    alpha = face_field(grid)
    alpha_half = center_field(grid)

    ql = face_field(grid)
    qi = face_field(grid)
    qv = face_field(grid)

    ql_half = center_field(grid)
    qi_half = center_field(grid)
    qv_half = center_field(grid)

    # Compute reference state thermodynamic profiles

    @inbounds for k in center_indicies(grid)
        ret = eos(p_half_[k], self.qtg, self.sg; t_to_prog = t_to_entropy_c, prog_to_t = eos_first_guess_entropy)
        temperature_half[k] = ret.T
        ql_half[k] = ret.ql
        qv_half[k] = self.qtg - (ql_half[k] + qi_half[k])
        alpha_half[k] = alpha_c(p_half_[k], temperature_half[k], self.qtg, qv_half[k])
    end

    @inbounds for k in face_indicies(grid)
        ret = eos(p_[k], self.qtg, self.sg; t_to_prog = t_to_entropy_c, prog_to_t = eos_first_guess_entropy)
        temperature[k] = ret.T
        ql[k] = ret.ql
        qv[k] = self.qtg - (ql[k] + qi[k])
        alpha[k] = alpha_c(p_[k], temperature[k], self.qtg, qv[k])
    end

    # Now do a sanity check to make sure that the Reference State entropy profile is uniform following
    # saturation adjustment
    local s
    @inbounds for k in center_indicies(grid)
        s = t_to_entropy_c(p_half[k], temperature_half[k], self.qtg, ql_half[k], qi_half[k])
        if abs(s - self.sg) / self.sg > 0.01
            println("Error in reference profiles entropy not constant !")
            println("Likely error in saturation adjustment")
        end
    end


    self.alpha0_half = alpha_half
    self.alpha0 = alpha
    self.p0 = p_
    self.p0_half = p_half
    self.rho0 = 1.0 ./ self.alpha0
    self.rho0_half = 1.0 ./ self.alpha0_half

    # TODO: centers and faces are sliced with equal sizes,
    # they should be unequal.
    cinterior = grid.cinterior
    finterior = grid.finterior
    add_reference_profile(Stats, "alpha0")
    write_reference_profile(Stats, "alpha0", alpha[finterior])
    add_reference_profile(Stats, "alpha0_half")
    write_reference_profile(Stats, "alpha0_half", alpha_half[cinterior])


    add_reference_profile(Stats, "p0")
    write_reference_profile(Stats, "p0", p_[finterior])
    add_reference_profile(Stats, "p0_half")
    write_reference_profile(Stats, "p0_half", p_half[cinterior])

    add_reference_profile(Stats, "rho0")
    write_reference_profile(Stats, "rho0", 1.0 ./ alpha[finterior])
    add_reference_profile(Stats, "rho0_half")
    write_reference_profile(Stats, "rho0_half", 1.0 ./ alpha_half[cinterior])

    return
end
