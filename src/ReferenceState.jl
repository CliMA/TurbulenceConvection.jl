####
#### Reference state
####

"""
    ReferenceState(
        grid::Grid,
        param_set::PS,
        Stats::NetCDFIO_Stats;
    ) where {FT}

The reference profiles, given
 - `grid` the grid
 - `param_set` the parameter set
 - `Stats` the NC file handler struct
"""
struct ReferenceState{PS, A1}
    param_set::PS
    p0::A1
    p0_half::A1
    alpha0::A1
    alpha0_half::A1
    rho0::A1
    rho0_half::A1
end

function ReferenceState(grid::Grid, param_set::PS, Stats::NetCDFIO_Stats; Pg::FT, Tg::FT, qtg::FT) where {PS, FT}

    q_pt_g = TD.PhasePartition(qtg)
    ts_g = TD.PhaseEquil_pTq(param_set, Pg, Tg, qtg)
    θ_liq_ice_g = TD.liquid_ice_pottemp(ts_g)

    # We are integrating the log pressure so need to take the log of the
    # surface pressure
    logp = log(Pg)

    # Form a right hand side for integrating the hydrostatic equation to
    # determine the reference pressure
    function rhs(logp, u, z)
        p_ = exp(logp)
        ts = TD.PhaseEquil_pθq(param_set, p_, θ_liq_ice_g, qtg)
        R_m = TD.gas_constant_air(ts)
        T = TD.air_temperature(ts)
        return -FT(CPP.grav(param_set)) / (T * R_m)
    end

    # Perform the integration
    z_span = (grid.zmin, grid.zmax)
    @show z_span
    prob = ODE.ODEProblem(rhs, logp, z_span)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)

    p0 = face_field(grid)
    p0_half = center_field(grid)
    alpha0 = face_field(grid)
    alpha0_half = center_field(grid)

    p0 .= sol.(vec(grid.zf))
    p0_half .= sol.(vec(grid.zc))

    p0 .= exp.(p0)
    p0_half .= exp.(p0_half)

    # Compute reference state thermodynamic profiles
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_half[k], θ_liq_ice_g, qtg)
        alpha0_half[k] = TD.specific_volume(ts)
    end

    @inbounds for k in real_face_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0[k], θ_liq_ice_g, qtg)
        alpha0[k] = TD.specific_volume(ts)
    end

    rho0 = 1 ./ alpha0
    rho0_half = 1 ./ alpha0_half

    args = (param_set, p0, p0_half, alpha0, alpha0_half, rho0, rho0_half)
    return ReferenceState{PS, typeof(p0)}(args...)
end
