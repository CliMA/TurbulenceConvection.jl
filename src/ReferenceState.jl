####
#### Reference state
####

"""
    compute_ref_state!(
        state,
        grid::Grid,
        param_set::PS,
        Stats::NetCDFIO_Stats;
        Pg::FT,
        Tg::FT,
        qtg::FT,
    ) where {PS, FT}

TODO: add better docs once the API converges

The reference profiles, given
 - `grid` the grid
 - `param_set` the parameter set
 - `Stats` the NC file handler struct
"""
function compute_ref_state!(
    state,
    grid::Grid,
    param_set::PS,
    Stats::NetCDFIO_Stats;
    Pg::FT,
    Tg::FT,
    qtg::FT,
) where {PS, FT}

    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    α0_c = center_ref_state(state).α0
    p0_f = face_ref_state(state).p0
    ρ0_f = face_ref_state(state).ρ0
    α0_f = face_ref_state(state).α0

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

    parent(p0_f) .= sol.(vec(grid.zf))
    parent(p0_c) .= sol.(vec(grid.zc))

    p0_f .= exp.(p0_f)
    p0_c .= exp.(p0_c)

    # Compute reference state thermodynamic profiles
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], θ_liq_ice_g, qtg)
        α0_c[k] = TD.specific_volume(ts)
    end

    @inbounds for k in real_face_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_f[k], θ_liq_ice_g, qtg)
        α0_f[k] = TD.specific_volume(ts)
    end

    ρ0_f .= 1 ./ α0_f
    ρ0_c .= 1 ./ α0_c
end
