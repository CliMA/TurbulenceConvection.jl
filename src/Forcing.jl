function initialize(self::ForcingBase, GMV::GridMeanVariables, ::ForcingBaseType)
    self.subsidence = center_field(self.Gr)
    self.dTdt = center_field(self.Gr)
    self.dqtdt = center_field(self.Gr)
    self.ug = center_field(self.Gr)
    self.vg = center_field(self.Gr)
    return
end

update(self::ForcingBase, GMV::GridMeanVariables) = nothing
initialize_io(self::ForcingBase, Stats) = nothing
io(self::ForcingBase, Stats) = nothing
initialize(self::ForcingBase{ForcingNone}, GMV::GridMeanVariables) = initialize(self, GMV, ForcingBaseType())
initialize(self::ForcingBase{ForcingStandard}, GMV::GridMeanVariables) = initialize(self, GMV, ForcingBaseType())
initialize(self::ForcingBase{ForcingStandard}, GMV) = initialize(self, GMV, ForcingBaseType())
initialize(self::ForcingBase{ForcingDYCOMS_RF01}, GMV::GridMeanVariables) = initialize(self, GMV, ForcingBaseType())

function initialize(self::ForcingBase{ForcingLES}, GMV::GridMeanVariables, Gr::Grid, LESDat::LESData)
    initialize(self, GMV, ForcingBaseType())
    Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        z_les_half = data.group["profiles"]["z_half"][:, 1]

        self.dtdt_hadv = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "dtdt_hadv", imin, imax))
        self.dtdt_nudge = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "dtdt_nudge", imin, imax))
        self.dtdt_fluc = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "dtdt_fluc", imin, imax))
        self.dqtdt_hadv = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "dqtdt_hadv", imin, imax))
        self.dqtdt_nudge = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "dqtdt_nudge", imin, imax))
        self.dqtdt_fluc = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "dqtdt_fluc", imin, imax))
        self.subsidence = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "ls_subsidence", imin, imax))
        self.u_nudge = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "u_mean", imin, imax))
        self.v_nudge = pyinterp(Gr.z_half, z_les_half, get_nc_data(data, "profiles", "v_mean", imin, imax))
    end
end
