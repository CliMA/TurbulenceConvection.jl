
function initialize(self::ForcingBase, GMV::GridMeanVariables, ::ForcingBaseType)
    self.subsidence = pyzeros(self.Gr.nzg)
    self.dTdt = pyzeros(self.Gr.nzg)
    self.dqtdt = pyzeros(self.Gr.nzg)
    self.ug = pyzeros(self.Gr.nzg)
    self.vg = pyzeros(self.Gr.nzg)

    if GMV.H.name == "s"
        self.convert_forcing_prog_fp = convert_forcing_entropy
    elseif GMV.H.name == "thetal"
        self.convert_forcing_prog_fp = convert_forcing_thetal
    end
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

    self.alpha_z    = 1.
    self.kappa      = 85.
    self.F0         = 70.
    self.F1         = 22.
    self.divergence = 3.75e-6  # divergence is defined twice: here and in initialize_forcing method of DYCOMS_RF01 case class
                               # where it is used to initialize large scale subsidence

    self.f_rad = pyzeros(self.Gr.nzg + 1) # radiative flux at cell edges
    return
end

"""
see eq. 3 in Stevens et. al. 2005 DYCOMS paper
"""
function calculate_radiation(self::ForcingBase{ForcingDYCOMS_RF01}, GMV::GridMeanVariables)
    # find zi (level of 8.0 g/kg isoline of qt)
    # TODO: report bug: zi and rhoi are not initialized
    zi = 0
    rhoi = 0
    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg - self.Gr.gw)
        if (GMV.QT.values[k] < 8.0 / 1000)
            idx_zi = k
            # will be used at cell edges
            zi     = self.Gr.z[idx_zi]
            rhoi   = self.Ref.rho0[idx_zi]
            break
        end
    end
    # cloud-top cooling
    q_0 = 0.0

    self.f_rad = pyzeros(self.Gr.nzg + 1)
    self.f_rad[self.Gr.nzg] = self.F0 * exp(-q_0)
    @inbounds for k in revxrange(self.Gr.nzg-1, 0, -1)
        q_0           += self.kappa * self.Ref.rho0_half[k] * GMV.QL.values[k] * self.Gr.dz
        self.f_rad[k]  = self.F0 * exp(-q_0)
    end

    # cloud-base warming
    q_1 = 0.0
    self.f_rad[0] += self.F1 * exp(-q_1)
    @inbounds for k in xrange(1, self.Gr.nzg + 1)
        q_1           += self.kappa * self.Ref.rho0_half[k - 1] * GMV.QL.values[k - 1] * self.Gr.dz
        self.f_rad[k] += self.F1 * exp(-q_1)
    end

    # cooling in free troposphere
    @inbounds for k in xrange(0, self.Gr.nzg)
        if self.Gr.z[k] > zi
            cbrt_z         = cbrt(self.Gr.z[k] - zi)
            self.f_rad[k] += rhoi * dycoms_cp * self.divergence * self.alpha_z * (power(cbrt_z, 4) / 4.0 + zi * cbrt_z)
        end
    end
    # TODO: report bug: k is not initialized
    # in the following expression, check that this
    # is equivalent.
    k_last = last(xrange(0, self.Gr.nzg))
    # condition at the top
    cbrt_z                   = cbrt(self.Gr.z[k_last] + self.Gr.dz - zi)
    self.f_rad[self.Gr.nzg] += rhoi * dycoms_cp * self.divergence * self.alpha_z * (power(cbrt_z, 4) / 4.0 + zi * cbrt_z)
    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg - self.Gr.gw)
        self.dTdt[k] = - (self.f_rad[k + 1] - self.f_rad[k]) / self.Gr.dz / self.Ref.rho0_half[k] / dycoms_cp
    end

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

    calculate_radiation(self, GMV)

    @inbounds for k in xrange(self.Gr.gw, self.Gr.nzg-self.Gr.gw)
        # Apply large-scale horizontal advection tendencies
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        GMV.H.tendencies[k]  += self.convert_forcing_prog_fp(self.Ref.p0_half[k],GMV.QT.values[k], qv, GMV.T.values[k], self.dqtdt[k], self.dTdt[k])
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
    add_profile(Stats, "rad_dTdt")
    add_profile(Stats, "subsidence")
    add_profile(Stats, "rad_flux")
    return
end

function io(self::ForcingBase{ForcingDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    cinterior = self.Gr.cinterior
    finterior = self.Gr.finterior
    write_profile(Stats, "rad_dTdt", self.dTdt[cinterior])
    write_profile(Stats, "subsidence", self.subsidence[cinterior])
    write_profile(Stats, "rad_flux", self.f_rad[finterior])
    return
end
