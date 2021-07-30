
function initialize(self::ForcingBase, GMV::GridMeanVariables, ::ForcingBaseType)
    self.subsidence = center_field(self.Gr)
    self.dTdt = center_field(self.Gr)
    self.dqtdt = center_field(self.Gr)
    self.ug = center_field(self.Gr)
    self.vg = center_field(self.Gr)
    self.convert_forcing_prog_fp = convert_forcing_thetal
    return
end

update(self::ForcingBase, GMV::GridMeanVariables) = nothing

function coriolis_force(self::ForcingBase, U::VariablePrognostic, V::VariablePrognostic, ::ForcingBaseType)
    gw = self.Gr.gw
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

    @inbounds for k in real_center_indicies(self.Gr)
        # Apply large-scale horizontal advection tendencies
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        GMV.H.tendencies[k] += self.convert_forcing_prog_fp(
            self.Ref.p0_half[k],
            GMV.QT.values[k],
            qv,
            GMV.T.values[k],
            self.dqtdt[k],
            self.dTdt[k],
        )
        GMV.QT.tendencies[k] += self.dqtdt[k]
    end
    if self.apply_subsidence
        @inbounds for k in real_center_indicies(self.Gr)
            # Apply large-scale subsidence tendencies
            GMV.H.tendencies[k] -= ∇_upwind(GMV.H.values, self.Gr, k) * self.subsidence[k]
            GMV.QT.tendencies[k] -= ∇_upwind(GMV.QT.values, self.Gr, k) * self.subsidence[k]
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

    @inbounds for k in real_center_indicies(self.Gr)
        GMV.QT.tendencies[k] += self.dqtdt[k]
        # Apply large-scale subsidence tendencies
        GMV.H.tendencies[k] -= ∇_upwind(GMV.H.values, self.Gr, k) * self.subsidence[k]
        GMV.QT.tendencies[k] -= ∇_upwind(GMV.QT.values, self.Gr, k) * self.subsidence[k]
    end

    if self.apply_coriolis
        coriolis_force(self, GMV.U, GMV.V)
    end

    return
end

initialize_io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats) = nothing
io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats) = nothing
