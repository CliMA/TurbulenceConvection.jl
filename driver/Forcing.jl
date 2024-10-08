initialize(::ForcingBase, grid, state) = nothing

function initialize(::ForcingBase{ForcingLES}, grid, state, LESDat::LESData)
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
    for k in TC.real_center_indices(grid)
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
