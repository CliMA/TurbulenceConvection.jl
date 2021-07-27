
function initialize(self::RadiationBase, GMV::GridMeanVariables, ::RadiationBaseType)
    self.dTdt = center_field(self.Gr)
    self.dqtdt = center_field(self.Gr)

    self.convert_forcing_prog_fp = convert_forcing_thetal
    return
end

update(self::RadiationBase, GMV::GridMeanVariables) = nothing

initialize_io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing
io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing

function initialize(self::RadiationBase{RadiationNone}, GMV::GridMeanVariables)
    initialize(self, GMV, RadiationBaseType())
end

update(self::RadiationBase{RadiationNone}, GMV::GridMeanVariables) = nothing
initialize_io(self::RadiationBase{RadiationNone}, Stats::NetCDFIO_Stats) = nothing
io(self::RadiationBase{RadiationNone}, Stats::NetCDFIO_Stats) = nothing


function initialize(self::RadiationBase{RadiationDYCOMS_RF01}, GMV::GridMeanVariables)
    initialize(self, GMV, RadiationBaseType())


    self.divergence = 3.75e-6  # divergence is defined twice: here and in initialize_forcing method of DYCOMS_RF01 case class
    # where it is used to initialize large scale subsidence
    self.alpha_z = 1.0
    self.kappa = 85.0
    self.F0 = 70.0
    self.F1 = 22.0
    self.divergence = 3.75e-6
    self.f_rad = pyzeros(self.Gr.nzg + 1) # radiative flux at cell faces
    return
end

"""
see eq. 3 in Stevens et. al. 2005 DYCOMS paper
"""
function calculate_radiation(self::RadiationBase{RadiationDYCOMS_RF01}, GMV::GridMeanVariables)
    # find zi (level of 8.0 g/kg isoline of qt)
    # TODO: report bug: zi and rhoi are not initialized
    zi = 0
    rhoi = 0
    @inbounds for k in real_center_indicies(self.Gr)
        if (GMV.QT.values[k] < 8.0 / 1000)
            idx_zi = k
            # will be used at cell faces
            zi = self.Gr.z[idx_zi]
            rhoi = self.Ref.rho0[idx_zi]
            break
        end
    end

    # cloud-top cooling
    q_0 = 0.0

    self.f_rad = pyzeros(self.Gr.nzg + 1)
    self.f_rad[self.Gr.nzg] = self.F0 * exp(-q_0)
    @inbounds for k in revxrange(self.Gr.nzg - 1, 0, -1)
        q_0 += self.kappa * self.Ref.rho0_half[k] * GMV.QL.values[k] * self.Gr.dz
        self.f_rad[k] = self.F0 * exp(-q_0)
    end

    # cloud-base warming
    q_1 = 0.0
    self.f_rad[0] += self.F1 * exp(-q_1)
    @inbounds for k in xrange(1, self.Gr.nzg + 1)
        q_1 += self.kappa * self.Ref.rho0_half[k - 1] * GMV.QL.values[k - 1] * self.Gr.dz
        self.f_rad[k] += self.F1 * exp(-q_1)
    end

    # cooling in free troposphere
    @inbounds for k in face_indicies(self.Gr)
        if self.Gr.z[k] > zi
            cbrt_z = cbrt(self.Gr.z[k] - zi)
            self.f_rad[k] += rhoi * dycoms_cp * self.divergence * self.alpha_z * (power(cbrt_z, 4) / 4.0 + zi * cbrt_z)
        end
    end
    # TODO: report bug: k is not initialized
    # in the following expression, check that this
    # is equivalent.
    k_last = last(xrange(0, self.Gr.nzg))
    # condition at the top
    cbrt_z = cbrt(self.Gr.z[k_last] + self.Gr.dz - zi)
    self.f_rad[self.Gr.nzg] +=
        rhoi * dycoms_cp * self.divergence * self.alpha_z * (power(cbrt_z, 4) / 4.0 + zi * cbrt_z)

    @inbounds for k in real_center_indicies(self.Gr)
        self.dTdt[k] = -(self.f_rad[k + 1] - self.f_rad[k]) / self.Gr.dz / self.Ref.rho0_half[k] / dycoms_cp
    end

    return
end

function update(self::RadiationBase{RadiationDYCOMS_RF01}, GMV::GridMeanVariables)

    calculate_radiation(self, GMV)

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
    end

    return
end

function initialize_io(self::RadiationBase{RadiationDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    add_profile(Stats, "rad_dTdt")
    add_profile(Stats, "rad_flux")
    return
end

function io(self::RadiationBase{RadiationDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    cinterior = self.Gr.cinterior
    finterior = self.Gr.finterior
    write_profile(Stats, "rad_dTdt", self.dTdt[cinterior])
    write_profile(Stats, "rad_flux", self.f_rad[finterior])
    return
end
