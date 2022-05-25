initialize(self::TC.ForcingBase, grid, state) = nothing

function initialize(self::TC.ForcingBase{TC.ForcingLES}, grid, state, LESDat::TC.LESData)
    aux_gm = TC.center_aux_grid_mean(state)
    nt = NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zc_les = Array(TC.get_nc_data(data, "zc"))

        dTdt_hadv = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "dtdt_hadv", imin, imax))
        T_nudge = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "temperature_mean", imin, imax))
        dTdt_fluc = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "dtdt_fluc", imin, imax))
        dqtdt_hadv = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "dqtdt_hadv", imin, imax))
        qt_nudge = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "qt_mean", imin, imax))
        dqtdt_fluc = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "dqtdt_fluc", imin, imax))
        subsidence = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "ls_subsidence", imin, imax))
        u_nudge = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "u_mean", imin, imax))
        v_nudge = TC.pyinterp(grid.zc, zc_les, TC.mean_nc_data(data, "profiles", "v_mean", imin, imax))
        (; dTdt_hadv, T_nudge, dTdt_fluc, dqtdt_hadv, qt_nudge, dqtdt_fluc, subsidence, u_nudge, v_nudge)
    end
    for k in TC.real_center_indices(grid)
        aux_gm.dTdt_hadv[k] = nt.dTdt_hadv[k]
        aux_gm.T_nudge[k] = nt.T_nudge[k]
        aux_gm.dTdt_fluc[k] = nt.dTdt_fluc[k]
        aux_gm.dqtdt_hadv[k] = nt.dqtdt_hadv[k]
        aux_gm.qt_nudge[k] = nt.qt_nudge[k]
        aux_gm.dqtdt_fluc[k] = nt.dqtdt_fluc[k]
        aux_gm.subsidence[k] = nt.subsidence[k]
        aux_gm.u_nudge[k] = nt.u_nudge[k]
        aux_gm.v_nudge[k] = nt.v_nudge[k]
    end
end
