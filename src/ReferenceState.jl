#### RefState

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

    FT = eltype(grid)
    P_g = self.Pg
    T_g = self.Tg
    q_tot_g = self.qtg
    param_set = parameter_set(self)
    q_pt_g = TD.PhasePartition(q_tot_g)
    ts_g = TD.PhaseEquil_pTq(param_set, P_g, T_g, q_tot_g)
    θ_liq_ice_g = TD.liquid_ice_pottemp(ts_g)

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    logp = log(P_g)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure
    function rhs(logp, u, z)
        p_ = exp(logp)
        ts = TD.PhaseEquil_pθq(param_set, p_, θ_liq_ice_g, q_tot_g)
        R_m = TD.gas_constant_air(ts)
        T = TD.air_temperature(ts)
        return -FT(CPP.grav(param_set)) / (T * R_m)
    end

    p = face_field(grid)
    p_half = center_field(grid)

    # Perform the integration
    z_span = (grid.zmin, grid.zmax)
    @show z_span
    prob = ODE.ODEProblem(rhs, logp, z_span)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)

    p .= sol.(vec(grid.zf))
    p_half .= sol.(vec(grid.zc))

    p .= exp.(p)
    p_half .= exp.(p_half)

    alpha = face_field(grid)
    alpha_half = center_field(grid)

    # Compute reference state thermodynamic profiles
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p_half[k], θ_liq_ice_g, q_tot_g)
        alpha_half[k] = TD.specific_volume(ts)
    end

    @inbounds for k in real_face_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p[k], θ_liq_ice_g, q_tot_g)
        alpha[k] = TD.specific_volume(ts)
    end

    self.alpha0_half = alpha_half
    self.alpha0 = alpha
    self.p0 = p
    self.p0_half = p_half
    self.rho0 = 1.0 ./ self.alpha0
    self.rho0_half = 1.0 ./ self.alpha0_half

    add_reference_profile(Stats, "alpha0")
    write_reference_profile(Stats, "alpha0", alpha)
    add_reference_profile(Stats, "alpha0_half")
    write_reference_profile(Stats, "alpha0_half", alpha_half)


    add_reference_profile(Stats, "p0")
    write_reference_profile(Stats, "p0", p)
    add_reference_profile(Stats, "p0_half")
    write_reference_profile(Stats, "p0_half", p_half)

    add_reference_profile(Stats, "rho0")
    write_reference_profile(Stats, "rho0", 1.0 ./ alpha)
    add_reference_profile(Stats, "rho0_half")
    write_reference_profile(Stats, "rho0_half", 1.0 ./ alpha_half)

    return
end
