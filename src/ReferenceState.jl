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
    prob = ODEProblem(rhs, logp, z_span)
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
    cinterior = kc_surface(grid):kc_top_of_atmos(grid)
    finterior = kf_surface(grid):kf_top_of_atmos(grid)
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
    @inbounds for k in center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p_half_[k], θ_liq_ice_g, q_tot_g)
        temperature_half[k] = TD.air_temperature(ts)
        ql_half[k] = TD.liquid_specific_humidity(ts)
        qv_half[k] = TD.vapor_specific_humidity(ts)
        alpha_half[k] = TD.specific_volume(ts)
    end

    @inbounds for k in face_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p_[k], θ_liq_ice_g, q_tot_g)
        temperature[k] = TD.air_temperature(ts)
        ql[k] = TD.liquid_specific_humidity(ts)
        qv[k] = TD.vapor_specific_humidity(ts)
        alpha[k] = TD.specific_volume(ts)
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
