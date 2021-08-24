
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

initialize_io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats) = nothing
io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats) = nothing
