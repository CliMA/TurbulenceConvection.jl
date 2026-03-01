mutable struct TimeStepping{FT, AT <: Val}
    const algorithm::AT # Val{:AlgorithmType} # does this really need to be in the type? we don't use it anywhere in the code really. ( I feel like having this is smaller than putting the actual algorithm in the type, which is a ODE.OrdinaryDiffEqAlgorithm (and I think not SciMLBase.AbstractODEAlgorithm)
    const fixed_step_solver::Bool
    const explicit_solver::Bool
    #
    adapt_dt::Bool
    dt::FT
    const dt_min::FT
    const dt_max::FT
    const t_max::FT
    t::FT
    const abstol::FT
    const reltol::FT
    #
    # nstep::Int # does this do anything?
    const cfl_limit::FT
    dt_max_edmf::FT
    dt_io::FT
    cfl_dt_max::FT
    #
    const spinup_half_t_max::FT # 1/2 of the spinup period -- we use Δt * spinup_dt_factor where spinup_dt_factor = spinup_dt_factor from 0 to spinup_half_t_max and then ramps up from spinup_dt_factor to 1 from t=spinup_half_t_max to t=2*spinup_half_t_max
    const spinup_dt_factor::FT
    const allow_spinup_adapt_dt::Bool # whether to allow dt to adapt during spinup. if adapt_dt is false, this should have no effect
    #
    #
    const use_tendency_timestep_limiter::Bool # whether or not to use tendencies to limit timestepping. [this is distinct (the reverse of) from using timstep to limit tendencies using tendency_limiters...]
    const use_isoutofdomain_limiter::Bool # whether or not to use isoutofdomain to limit timestepping.
    N_dt_max_edmf_violate_dt_min_remaining::Int # The number of times we allow dt_max_edmf to violate dt_min before we stop allowing violations
    const adaptive_depth_limit::Int # how many recursive calls we allow to recalculate tendencies for new proposed dt
    const allow_cfl_dt_max_violate_dt_min::Bool # whether to let falling precip violate the cfl (we calculated w/ a max of 5 m/s but maybe that's arbitrary -- rain can in principle go up to 10 or so).
    use_fallback_tendency_limiters::Bool # check if this could be in a callback or something w/ remake to just remake the ode problem (remake edmf w/ limiters hopefully...)
    const use_fallback_during_spinup::Bool # whether to use fallback limiters during spinup
    const dt_limit_tendencies_factor::FT # whether to calculate tendencies using dt_min or the current dt -- I think we should default to dt_min bc (though there's some risk they dont compose nicely? hopefully fine in steady state... could do dt_max but dt_max could be arbitrarily large and then cold never take shorter steps... could also pick some fixed number...)
    const limit_tendencies_by_dt_min::Bool # whether to limit by dt_min or current dt.
    #
    isoutofdomain::Bool # whether or not we are out of domain
end

# To handle the dichotomy between starting w/ N_dt_max_edmf_violate_dt_min_remaining > 0, using them all up and transitioning to Fallbacks once you get to 0
# vs
# starting w/ N_dt_max_edmf_violate_dt_min_remaining = 0 and never switching to Fallbacks
# when you reach N_dt_max_edmf_violate_dt_min_remaining by decline, we set it to -1 when we set use_fallback_tendency_limiters to true so we can disambiguate states...
# This also asllows use_fallback_during_spinup to set use_fallback_tendency_limiters to true, and then when coming out of spinup, we know which state to revert to, true if -1, false if >= 0 


function TimeStepping(::Type{FT}, namelist) where {FT}

    local explicit_solver::Bool
    local fixed_step_solver::Bool

    algorithm = TC.parse_namelist(namelist, "time_stepping", "algorithm"; default = "Euler") # I believe strings are fine and won't go to param_set
    algorithm = Val(Symbol(algorithm)) # create a Val type for the algorithm

    # always use their adaptive dt if it's an adaptive solver (adaptive = true I believe is the default)
    if algorithm isa Val{:Euler}
        explicit_solver = true
        fixed_step_solver = true
    elseif algorithm isa Val{:Heun}
        explicit_solver = true
        fixed_step_solver = false
    elseif algorithm isa Val{:RK4}
        explicit_solver = true
        fixed_step_solver = false
    elseif algorithm isa Val{:Tsit5}
        explicit_solver = true
        fixed_step_solver = false
    elseif algorithm isa Val{:ImplicitEuler}
        explicit_solver = false
        fixed_step_solver = false
    elseif algorithm isa Val{:KenCarp47}
        explicit_solver = false
        fixed_step_solver = false
    else
        error("Unsupported algorithm: $(algorithm)")
    end


    if fixed_step_solver
        adapt_dt = TC.parse_namelist(namelist, "time_stepping", "adapt_dt"; default = false)
        dt = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = FT(1.0))
        dt_min = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = FT(1.0)) # placeholder to keep storing it in case spinup needs it since it gets edited in callbacks
        dt_max = TC.parse_namelist(namelist, "time_stepping", "dt_max"; default = FT(10.0))

        adaptive_depth_limit = TC.parse_namelist(namelist, "time_stepping", "adaptive_depth_limit"; default = 2) # how many recursive calls we allow to recalculate tendencies for new proposed dt
        @assert adaptive_depth_limit >= 1 "adaptive_depth_limit must be greater than or equal to 0, given value $(adaptive_depth_limit) is invalid"
        use_tendency_timestep_limiter =
            explicit_solver ?
            TC.parse_namelist(namelist, "time_stepping", "use_tendency_timestep_limiter"; default = false) : false
        use_isoutofdomain_limiter = false # isoutofdomain is only used for adaptive timestepping to accept or reject timesteps

        if use_tendency_timestep_limiter
            @assert (adapt_dt == true) "Tendency limiting timesteps is adaptive. If you want to use the tendency limiter, you must also set adapt_dt to true."
        end

        spinup_half_t_max = TC.parse_namelist(namelist, "time_stepping", "spinup_half_t_max"; default = FT(0.0))
        spinup_dt_factor = TC.parse_namelist(namelist, "time_stepping", "spinup_dt_factor"; default = FT(1.0))
        allow_spinup_adapt_dt =
            adapt_dt ?
            TC.parse_namelist(
                namelist,
                "time_stepping",
                "allow_spinup_adapt_dt";
                default = (adapt_dt || use_tendency_timestep_limiter),
            ) : false
        N_dt_max_edmf_violate_dt_min_remaining =
            TC.parse_namelist(namelist, "time_stepping", "N_dt_max_edmf_violate_dt_min"; default = 0)

        # -------------------------------- #
        reltol = FT(NaN) # not used here
        abstol = FT(NaN) # not used here

    else
        adapt_dt = true
        dt = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = FT(NaN)) # use their defaults
        dt_min = TC.parse_namelist(namelist, "time_stepping", "dt_min"; default = FT(NaN)) # use their defaults
        dt_max = TC.parse_namelist(namelist, "time_stepping", "dt_max"; default = FT(Inf)) # use their defaults
        adaptive_depth_limit = TC.parse_namelist(namelist, "time_stepping", "adaptive_depth_limit"; default = -1) # default to some negative number so we can leave their defaults in place

        use_tendency_timestep_limiter = false
        use_isoutofdomain_limiter = true # use this to accept/reject steps

        reltol = TC.parse_namelist(namelist, "time_stepping", "reltol"; default = FT(NaN)) # use their defaults
        abstol = TC.parse_namelist(namelist, "time_stepping", "abstol"; default = FT(NaN)) # use their defaults

        # -------------------------------- #
        spinup_half_t_max = FT(0.0) # not used w/ adaptive
        spinup_dt_factor = FT(1.0) # not used w/ adaptive
        allow_spinup_adapt_dt = false # not used w/ adaptive
        N_dt_max_edmf_violate_dt_min_remaining = 0 # Not used here

    end


    if !isnan(dt_min) && !isnan(dt_max)
        @assert dt_min <= dt_max "dt_min must be less than or equal to dt_max, given values dt_min = $(dt_min) and dt_max = $(dt_max) are invalid"
    end

    t_max = TC.parse_namelist(namelist, "time_stepping", "t_max"; default = FT(7200.0))
    cfl_limit = TC.parse_namelist(namelist, "time_stepping", "cfl_limit"; default = FT(0.5))
    cfl_dt_max = TC.parse_namelist(namelist, "time_stepping", "cfl_dt_max"; default = dt_max)
    allow_cfl_dt_max_violate_dt_min =
        TC.parse_namelist(namelist, "time_stepping", "allow_cfl_dt_max_violate_dt_min"; default = false)

    use_fallback_tendency_limiters = false # initialize as false for everything... it'll get set to true later if needed.
    use_fallback_during_spinup =
        TC.parse_namelist(namelist, "time_stepping", "use_fallback_during_spinup"; default = false)

    # [[ deprecated ]] [ this was bad because, say if dt was longer than dt_limit_tendencies, then you're actually undoing your limiter. it's the real TS.dt that matters for most things, not the Δt you pass to ∑tendencies!()
    # so passing in a custom dt to calculate the tendencies, but then using the real one in ODE solver did bad things. We now still offer a factor but base it off the real TS.dt ]
    # limit_tendencies_using_dt_min_factor = TC.parse_namelist(namelist, "time_stepping", "limit_tendencies_using_dt_min_factor"; default = missing) # nothing triggers a special path in parse_namelist lmao
    # if ismissing(limit_tendencies_using_dt_min_factor)
    #     dt_limit_tendencies = TC.parse_namelist(namelist, "time_stepping", "dt_limit_tendencies"; default = dt_min) # we dont need this, you can just calculate it yourself, deprecate this whole section...
    # else
    #     dt_limit_tendencies = dt_min * limit_tendencies_using_dt_min_factor
    # end

    dt_limit_tendencies_factor =
        TC.parse_namelist(namelist, "time_stepping", "dt_limit_tendencies_factor"; default = FT(1)) # for consistency, tendencies are limited by this factor times dt_min (so adapt_dt doesn't oscillate)
    limit_tendencies_by_dt_min =
        TC.parse_namelist(namelist, "time_stepping", "limit_tendencies_by_dt_min"; default = false) # whether to limit by dt_min instead of current dt. The former is at least consistent in magnitude, but if dt is changing bad things can happen. Limiting by dt is prone to wild swings in the limited tendencies though.

    dt_max_edmf = FT(0)  # initialize


    @assert (FT(0) <= spinup_half_t_max) "spinup_half_t_max must be greater than or equal to 0, given value $(spinup_half_t_max) is invalid"
    @assert (FT(0) < spinup_dt_factor) "spinup_dt_factor must be greater than 0, given value $(spinup_dt_factor) is invalid"

    if spinup_half_t_max > FT(0)
        dt = dt * spinup_dt_factor # enforce this spinup factor for the first timestep.... otherwise it doesnt seem to happen until after the first timestep when callbacks are called.
    end


    isoutofdomain = false # start false

    # set time
    t = FT(0)
    dt_io = FT(0)
    # nstep = 0

    AT = typeof(algorithm)

    return TimeStepping{FT, AT}(
        algorithm,
        fixed_step_solver,
        explicit_solver,
        #
        adapt_dt,
        dt,
        dt_min,
        dt_max,
        t_max,
        t,
        abstol,
        reltol,
        #
        # nstep,
        cfl_limit,
        dt_max_edmf,
        dt_io,
        cfl_dt_max,
        #
        spinup_half_t_max,
        spinup_dt_factor,
        allow_spinup_adapt_dt,
        #
        use_tendency_timestep_limiter,
        use_isoutofdomain_limiter,
        N_dt_max_edmf_violate_dt_min_remaining,
        adaptive_depth_limit,
        allow_cfl_dt_max_violate_dt_min,
        use_fallback_tendency_limiters,
        use_fallback_during_spinup,
        # dt_limit_tendencies, # for consistency, tendencies are limited by this factor times dt_min (so adapt_dt doesn't oscillate) [ this is bad because as the timestep changes, you can induce instability, e.g., consider you have a dt longer than dt_limit_tendencies, then you've underdone your limiter! Essentially the relative magnitudes of things were constantly changing if the real dt was. if you allowed more adaptivity, this would be even worse.] [ essentially, this gave nice looking consistently smooth tendencies, but the integrated values suffered. the reverse is actually more important.]
        dt_limit_tendencies_factor, # just stick to a fixed factor of dt. this won't be as stationary in the calculated tendencies, but the integrated tendencies, i.e. the real values, will be stable and properly limited, which is what actually matters.
        limit_tendencies_by_dt_min,
        #
        isoutofdomain,
    )
end
