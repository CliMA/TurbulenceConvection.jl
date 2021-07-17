
function initialize(self::ForcingBase, GMV::GridMeanVariables, ::ForcingBaseType)
    self.subsidence = pyzeros(self.Gr.nzg)
    self.dTdt = pyzeros(self.Gr.nzg)
    self.dqtdt = pyzeros(self.Gr.nzg)
    self.ug = pyzeros(self.Gr.nzg)
    self.vg = pyzeros(self.Gr.nzg)
    self.convert_forcing_prog_fp = convert_forcing_thetal
    return
end

update(self::ForcingBase, GMV::GridMeanVariables) = nothing

function coriolis_force(
    self::ForcingBase,
    U::VariablePrognostic,
    V::VariablePrognostic,
    ::ForcingBaseType)
    gw = self.Gr.gw
    @inbounds for k in xrange(gw, self.Gr.nzg-gw)
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

    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        # Apply large-scale horizontal advection tendencies
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        GMV.H.tendencies[k] += self.convert_forcing_prog_fp(self.Ref.p0_half[k],GMV.QT.values[k],
                                                            qv, GMV.T.values[k], self.dqtdt[k], self.dTdt[k])
        GMV.QT.tendencies[k] += self.dqtdt[k]
    end
    if self.apply_subsidence
        @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
            # Apply large-scale subsidence tendencies
            GMV.H.tendencies[k] -= (GMV.H.values[k+1]-GMV.H.values[k]) * self.Gr.dzi * self.subsidence[k]
            GMV.QT.tendencies[k] -= (GMV.QT.values[k+1]-GMV.QT.values[k]) * self.Gr.dzi * self.subsidence[k]
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

function coriolis_force(
        self::ForcingBase{ForcingDYCOMS_RF01},
        U::VariablePrognostic,
        V::VariablePrognostic
    )
    coriolis_force(self, U, V, ForcingBaseType())
    return
end

function update(self::ForcingBase{ForcingDYCOMS_RF01}, GMV::GridMeanVariables)

    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        GMV.QT.tendencies[k] += self.dqtdt[k]
        # Apply large-scale subsidence tendencies
        GMV.H.tendencies[k]  -= (GMV.H.values[k+1]-GMV.H.values[k]) * self.Gr.dzi * self.subsidence[k]
        GMV.QT.tendencies[k] -= (GMV.QT.values[k+1]-GMV.QT.values[k]) * self.Gr.dzi * self.subsidence[k]
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


function initialize(self::ForcingBase{ForcingLES}, GMV::GridMeanVariables)
    initialize(self, GMV, ForcingBaseType())

    les_data = Dataset(self.les_filename, "r")
    # find the indexes for a 6h interval at the end of the simualtion
    t = data.group["profiles"]["t"][:];
    imin = Int(3600.0*6.0/(t[2]-t[1]))

    self.les_dtdt_hadv  = mean(data.group["profiles"]["dtdt_hadv"][:][:,imin:end],dims = 2)
    self.les_dtdt_nudge  = mean(data.group["profiles"]["dtdt_nudge"][:][:,imin:end],dims = 2)
    self.les_dtdt_fluc = mean(data.group["profiles"]["dtdt_fluc"][:][:,imin:end],dims = 2)
    self.les_dqtdt_hadv  = mean(data.group["profiles"]["dqtdt_hadv"][:][:,imin:end],dims = 2)
    self.les_dqtdt_nudge = mean(data.group["profiles"]["dqtdt_nudge"][:][:,imin:end],dims = 2)
    self.les_dqtdt_fluc  = mean(data.group["profiles"]["dqtdt_fluc"][:][:,imin:end],dims = 2)
    self.les_subsidence  = mean(data.group["profiles"]["ls_subsidence"][:][:,imin:end],dims = 2)
    self.les_u_nudge  = mean(data.group["profiles"]["u_mean"][:][:,imin:end],dims = 2)
    self.les_v_nudge  = mean(data.group["profiles"]["v_mean"][:][:,imin:end],dims = 2)
end

initialize(self::ForcingBase{ForcingLES}, GMV) = initialize(self, GMV, ForcingBaseType())

function update(self::ForcingBase{ForcingLES}, GMV::GridMeanVariables)

    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        GMV.H.horz_adv[k] = self.convert_forcing_prog_fp(self.Ref.p0_half[k],GMV.QT.values[k], qv, GMV.T.values[k], qv,self.dtdt_hadv[k])
        GMV.H.nudge[k] = self.convert_forcing_prog_fp(self.Ref.p0_half[k],GMV.QT.values[k], qv, GMV.T.values[k], qv,self.dtdt_nudge[k])
        GMV.H.fluc[k] = self.convert_forcing_prog_fp(self.Ref.p0_half[k],GMV.QT.values[k], qv, GMV.T.values[k], qv, self.dtdt_fluc[k])
        GMV.QT.horz_adv[k] = self.dqtdt_hadv[k]
        GMV.QT.nudge[k] = self.dqtdt_nudge[k]
        GMV.QT.fluc[k] = self.dqtdt_fluc[k]
        GMV.U.nudge[k] = (self.u_nudge[k] - GMV.U.values[k])/self.nudge_tau
        GMV.V.nudge[k] = (self.v_nudge[k] - GMV.V.values[k])/self.nudge_tau
        if self.apply_subsidence:
            # Apply large-scale subsidence tendencies
            GMV.H.subsidence[k] =  -(GMV.H.values[k+1]-GMV.H.values[k]) * self.Gr.dzi * self.scm_subsidence[k]
            GMV.QT.subsidence[k] =  -(GMV.QT.values[k+1]-GMV.QT.values[k]) * self.Gr.dzi * self.scm_subsidence[k]
        else:
            GMV.H.subsidence[k] =  0.0
            GMV.QT.subsidence[k] = 0.0

        GMV.H.tendencies[k] += GMV.H.horz_adv[k] + GMV.H.nudge[k] + GMV.H.subsidence[k] + GMV.H.fluc[k]
        GMV.QT.tendencies[k] += GMV.QT.horz_adv[k] + GMV.QT.nudge[k] + GMV.QT.subsidence[k] + GMV.QT.fluc[k]

        GMV.U.tendencies[k] += GMV.U.nudge[k]
        GMV.V.tendencies[k] += GMV.V.nudge[k]

        # Apply large-scale horizontal advection tendencies
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        GMV.H.tendencies[k] += self.convert_forcing_prog_fp(self.Ref.p0_half[k],GMV.QT.values[k],
                                                            qv, GMV.T.values[k], self.dqtdt[k], self.dTdt[k])
        GMV.QT.tendencies[k] += self.dqtdt[k]
    end
    if self.apply_subsidence
        @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
            # Apply large-scale subsidence tendencies
            GMV.H.tendencies[k] -= (GMV.H.values[k+1]-GMV.H.values[k]) * self.Gr.dzi * self.subsidence[k]
            GMV.QT.tendencies[k] -= (GMV.QT.values[k+1]-GMV.QT.values[k]) * self.Gr.dzi * self.subsidence[k]
        end
    end
    if self.apply_coriolis
        coriolis_force(self, GMV.U, GMV.V, ForcingBaseType())
    end
    return
end
initialize_io(self::ForcingBase{ForcingLES}, Stats) = nothing
io(self::ForcingBase{ForcingLES}, Stats) = nothing
