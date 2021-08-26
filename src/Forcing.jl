
function initialize(self::ForcingBase, GMV::GridMeanVariables, ::ForcingBaseType)
    self.subsidence = center_field(self.Gr)
    self.dTdt = center_field(self.Gr)
    self.dqtdt = center_field(self.Gr)
    self.ug = center_field(self.Gr)
    self.vg = center_field(self.Gr)
    return
end

update(self::ForcingBase, GMV::GridMeanVariables) = nothing

function coriolis_force(self::ForcingBase, U::VariablePrognostic, V::VariablePrognostic, ::ForcingBaseType)
    @inbounds for k in real_center_indicies(self.Gr)
        U.tendencies[k] -= self.coriolis_param * (self.vg[k] - V.values[k])
        V.tendencies[k] += self.coriolis_param * (self.ug[k] - U.values[k])
    end
    return
end

initialize_io(self::ForcingBase, Stats) = nothing
io(self::ForcingBase, Stats) = nothing


function initialize(self::ForcingBase{ForcingNone}, GMV::GridMeanVariables)
    initialize(self, GMV, ForcingBaseType())
end

update(self::ForcingBase{ForcingNone}, GMV::GridMeanVariables) = nothing
coriolis_force(self::ForcingBase{ForcingNone}, U, V) = nothing
initialize_io(self::ForcingBase{ForcingNone}, Stats) = nothing
io(self::ForcingBase{ForcingNone}, Stats) = nothing


function initialize(self::ForcingBase{ForcingStandard}, GMV::GridMeanVariables)
    initialize(self, GMV, ForcingBaseType())
end

initialize(self::ForcingBase{ForcingStandard}, GMV) = initialize(self, GMV, ForcingBaseType())

function update(self::ForcingBase{ForcingStandard}, GMV::GridMeanVariables)
    grid = self.Gr
    @inbounds for k in real_center_indicies(grid)
        # Apply large-scale horizontal advection tendencies
        GMV.H.tendencies[k] += self.dTdt[k] / exner_c(self.Ref.p0_half[k])
        GMV.QT.tendencies[k] += self.dqtdt[k]
    end
    if self.apply_subsidence
        @inbounds for k in real_center_indicies(grid)
            # Apply large-scale subsidence tendencies
            H_cut = cut_onesided(GMV.H.values, grid, k)
            q_tot_cut = cut_onesided(GMV.QT.values, grid, k)
            ∇H = ∇_onesided(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
            ∇q_tot = ∇_onesided(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
            GMV.H.tendencies[k] -= ∇H * self.subsidence[k]
            GMV.QT.tendencies[k] -= ∇q_tot * self.subsidence[k]
        end
    end
    if self.apply_coriolis
        coriolis_force(self, GMV.U, GMV.V, ForcingBaseType())
    end
    return
end
initialize_io(self::ForcingBase{ForcingStandard}, Stats) = nothing
io(self::ForcingBase{ForcingStandard}, Stats) = nothing

function initialize(self::ForcingBase{ForcingDYCOMS_RF01}, GMV::GridMeanVariables)
    initialize(self, GMV, ForcingBaseType())
    return
end

function coriolis_force(self::ForcingBase{ForcingDYCOMS_RF01}, U::VariablePrognostic, V::VariablePrognostic)
    coriolis_force(self, U, V, ForcingBaseType())
    return
end

function update(self::ForcingBase{ForcingDYCOMS_RF01}, GMV::GridMeanVariables)
    grid = self.Gr
    @inbounds for k in real_center_indicies(grid)
        H_cut = cut_onesided(GMV.H.values, grid, k)
        q_tot_cut = cut_onesided(GMV.QT.values, grid, k)
        ∇H = ∇_onesided(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
        ∇q_tot = ∇_onesided(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))

        GMV.QT.tendencies[k] += self.dqtdt[k]
        # Apply large-scale subsidence tendencies
        GMV.H.tendencies[k] -= ∇H * self.subsidence[k]
        GMV.QT.tendencies[k] -= ∇q_tot * self.subsidence[k]
    end

    if self.apply_coriolis
        coriolis_force(self, GMV.U, GMV.V)
    end

    return
end

function initialize_io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    return
end

function io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    return
end

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

function update(self::ForcingBase{ForcingLES}, GMV::GridMeanVariables)

    grid = self.Gr
    @inbounds for k in real_center_indicies(grid)
        H_horz_adv = self.dtdt_hadv[k] / exner_c(self.Ref.p0_half[k])
        H_nudge = self.dtdt_nudge[k] / exner_c(self.Ref.p0_half[k])
        H_fluc = self.dtdt_fluc[k] / exner_c(self.Ref.p0_half[k])

        GMV_U_nudge_k = (self.u_nudge[k] - GMV.U.values[k]) / self.nudge_tau
        GMV_V_nudge_k = (self.v_nudge[k] - GMV.V.values[k]) / self.nudge_tau
        if self.apply_subsidence
            # Apply large-scale subsidence tendencies

            H_cut = cut_onesided(GMV.H.values, grid, k)
            q_tot_cut = cut_onesided(GMV.QT.values, grid, k)
            ∇H = ∇_onesided(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
            ∇q_tot = ∇_onesided(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
            GMV_H_subsidence_k = -∇H * self.subsidence[k]
            GMV_QT_subsidence_k = -∇q_tot * self.subsidence[k]
        else
            GMV_H_subsidence_k = 0.0
            GMV_QT_subsidence_k = 0.0
        end

        GMV.H.tendencies[k] += H_horz_adv + H_nudge + H_fluc + GMV_H_subsidence_k
        GMV.QT.tendencies[k] += self.dqtdt_hadv[k] + self.dqtdt_nudge[k] + GMV_QT_subsidence_k + self.dqtdt_fluc[k]

        GMV.U.tendencies[k] += GMV_U_nudge_k
        GMV.V.tendencies[k] += GMV_V_nudge_k
    end

    if self.apply_coriolis
        coriolis_force(self, GMV.U, GMV.V, ForcingBaseType())
    end
    return
end
initialize_io(self::ForcingBase{ForcingLES}, Stats) = nothing
io(self::ForcingBase{ForcingLES}, Stats) = nothing
