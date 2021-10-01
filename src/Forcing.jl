function initialize(self::ForcingBase, GMV::GridMeanVariables, ::ForcingBaseType)
    self.subsidence = center_field(self.grid)
    self.dTdt = center_field(self.grid)
    self.dqtdt = center_field(self.grid)
    self.ug = center_field(self.grid)
    self.vg = center_field(self.grid)
    return
end

initialize_io(self::ForcingBase, Stats) = nothing
io(self::ForcingBase, Stats) = nothing
initialize(self::ForcingBase{ForcingNone}, GMV::GridMeanVariables) = initialize(self, GMV, ForcingBaseType())
initialize(self::ForcingBase{ForcingStandard}, GMV::GridMeanVariables) = initialize(self, GMV, ForcingBaseType())
initialize(self::ForcingBase{ForcingStandard}, GMV) = initialize(self, GMV, ForcingBaseType())
initialize(self::ForcingBase{ForcingDYCOMS_RF01}, GMV::GridMeanVariables) = initialize(self, GMV, ForcingBaseType())

function initialize(self::ForcingBase{ForcingLES}, GMV::GridMeanVariables, grid::Grid, LESDat::LESData)
    initialize(self, GMV, ForcingBaseType())
    NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zc_les = get_nc_data(data, "zc")

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
