initialize(::ForcingBase, state) = nothing

function initialize(::ForcingBase{ForcingLES}, state, LESDat::LESData)
    grid = TC.Grid(state)
    aux_gm = TC.center_aux_grid_mean(state)
    nt = NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zc_les = Array(TC.get_nc_data(data, "zc"))
        getvar(var) = TC.pyinterp(vec(grid.zc.z), zc_les, TC.mean_nc_data(data, "profiles", var, imin, imax))
        get_nudgevar(var) = TC.pyinterp(vec(grid.zc.z), zc_les, TC.init_nc_data(data, "profiles", var))

        dTdt_hadv = getvar("dtdt_hadv")
        H_nudge = get_nudgevar("thetali_mean")
        dTdt_fluc = getvar("dtdt_fluc")
        dqtdt_hadv = getvar("dqtdt_hadv")
        qt_nudge = get_nudgevar("qt_mean")
        dqtdt_fluc = getvar("dqtdt_fluc")
        subsidence = getvar("ls_subsidence")
        u_nudge = get_nudgevar("u_mean")
        v_nudge = get_nudgevar("v_mean")
        (; dTdt_hadv, H_nudge, dTdt_fluc, dqtdt_hadv, qt_nudge, dqtdt_fluc, subsidence, u_nudge, v_nudge)
    end
    @inbounds for k in TC.real_center_indices(grid)
        aux_gm.dTdt_hadv[k] = nt.dTdt_hadv[k]
        aux_gm.H_nudge[k] = nt.H_nudge[k]
        aux_gm.dTdt_fluc[k] = nt.dTdt_fluc[k]
        aux_gm.dqtdt_hadv[k] = nt.dqtdt_hadv[k]
        aux_gm.qt_nudge[k] = nt.qt_nudge[k]
        aux_gm.dqtdt_fluc[k] = nt.dqtdt_fluc[k]
        aux_gm.subsidence[k] = nt.subsidence[k]
        aux_gm.u_nudge[k] = nt.u_nudge[k]
        aux_gm.v_nudge[k] = nt.v_nudge[k]
    end
end

function initialize(forcing::ForcingBase{ForcingSOCRATES}, state) #where {T <: ForcingSOCRATES} # added param_set so we can calculate stuff....
    grid = TC.Grid(state)
    FT = TC.float_type(state)
    aux_gm = TC.center_aux_grid_mean(state)

    forcing_keys = (:dTdt_hadv, :H_nudge, :dqtdt_hadv, :qt_nudge, :subsidence, :u_nudge, :v_nudge) # all socrates forcings except those handled separately
    # ug_keys = (:ug_nudge, :vg_nudge)
    # rad_keys = (:dTdt_rad,)
    forcing_funcs = forcing.forcing_funcs[]

    # set the geostrophic velocity profiles -- need to check if we actually have a fcn that can take in t=0 and return a profile... ( i think ours is one func for each z so whoops... maybe need to reconstruct?)
    # g_func = (f) -> f([FT(0)])[1]
    # g_func = (f) -> f(FT(0)) # take fcn and evluate at time 0
    # prof_ug = g_func.(forcing_funcs[:ug_nudge]) # map over the forcing funcs to get the profile at t=0
    # prof_vg = g_func.(forcing_funcs[:vg_nudge])
    # # it wants a fcn out, could edit src/Fields.jl I guess to add another method but maybe it needs to face/center points idk...
    # prof_ug = Dierckx.Spline1D(vec(grid.zc.z), vec(prof_ug); k = 1)
    # prof_vg = Dierckx.Spline1D(vec(grid.zc.z), vec(prof_vg); k = 1)

    aux_gm_uₕ_g = TC.grid_mean_uₕ_g(state)
    # TC.set_z!(aux_gm_uₕ_g, prof_ug, prof_vg)

    @inbounds for k in TC.real_center_indices(grid)
        TC.set_z!(k, aux_gm_uₕ_g, forcing_funcs[:ug_nudge][k](FT(0)), forcing_funcs[:vg_nudge][k](FT(0))) # We do the interp in SSCF so we dont need to redo it
    end

    forcing_funcs = forcing_funcs[forcing_keys] # keys we don't need.
    for (name, funcs) in zip(keys(forcing_funcs), forcing_funcs) # iterate over the named tuple of our forcings...
        @inbounds for k in TC.real_center_indices(grid)
            func = funcs[k]
            # getproperty(aux_gm, name)[k] = func([FT(0)])[1] # apply to time = 0 and apply to aux_gm, turn to vec cause needs to be cast as in https://github.com/CliMA/TurbulenceConvection.jl/blob/a9ebce1f5f15f049fc3719a013ddbc4a9662943a/src/utility_functions.jl#L48
            getproperty(aux_gm, name)[k] = func(FT(0)) # apply to time = 0 and apply to aux_gm, turn to vec cause needs to be cast as in
        end
    end
    return nothing
end
