function initialize(self::ForcingBase, grid, ::ForcingBaseType)
    self.subsidence = center_field(grid)
    self.dTdt = center_field(grid)
    self.dqtdt = center_field(grid)
    self.ug = center_field(grid)
    self.vg = center_field(grid)
    return
end

initialize(self::ForcingBase{ForcingNone}, grid) = initialize(self, grid, ForcingBaseType())
initialize(self::ForcingBase{ForcingStandard}, grid) = initialize(self, grid, ForcingBaseType())
initialize(self::ForcingBase{ForcingDYCOMS_RF01}, grid) = initialize(self, grid, ForcingBaseType())

function initialize(self::ForcingBase{ForcingLES}, grid, LESDat::LESData)
    initialize(self, grid, ForcingBaseType())
    NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zc_les = Array(get_nc_data(data, "zc"))

        self.dtdt_hadv = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "dtdt_hadv", imin, imax))
        self.dtdt_nudge = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "dtdt_nudge", imin, imax))
        self.dtdt_fluc = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "dtdt_fluc", imin, imax))
        self.dqtdt_hadv = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "dqtdt_hadv", imin, imax))
        self.dqtdt_nudge = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "dqtdt_nudge", imin, imax))
        self.dqtdt_fluc = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "dqtdt_fluc", imin, imax))
        self.subsidence = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "ls_subsidence", imin, imax))
        self.u_nudge = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "u_mean", imin, imax))
        self.v_nudge = pyinterp(grid.zc, zc_les, mean_nc_data(data, "profiles", "v_mean", imin, imax))
    end
end
