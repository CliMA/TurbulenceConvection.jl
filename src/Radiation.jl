
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
    self.f_rad = face_field(self.Gr)
    return
end

"""
see eq. 3 in Stevens et. al. 2005 DYCOMS paper
"""
function calculate_radiation(self::RadiationBase{RadiationDYCOMS_RF01}, GMV::GridMeanVariables)
    # find zi (level of 8.0 g/kg isoline of qt)
    # TODO: report bug: zi and ρ_i are not initialized
    zi = 0
    ρ_i = 0
    @inbounds for k in real_center_indicies(self.Gr)
        if (GMV.QT.values[k] < 8.0 / 1000)
            idx_zi = k
            # will be used at cell faces
            zi = self.Gr.z[idx_zi]
            ρ_i = self.Ref.rho0[idx_zi]
            break
        end
    end

    ρ_z = Dierckx.Spline1D([self.Gr.z_half...], [self.Ref.rho0_half...]; k = 1)
    q_liq_z = Dierckx.Spline1D([self.Gr.z_half...], [GMV.QL.values...]; k = 1)

    integrand(ρq_l, params, z) = params.κ * ρ_z(z) * q_liq_z(z)
    rintegrand(ρq_l, params, z) = -integrand(ρq_l, params, z)

    z_span = (self.Gr.zmin, self.Gr.zmax)
    rz_span = (self.Gr.zmax, self.Gr.zmin)
    params = (; κ = self.kappa)

    rprob = ODEProblem(rintegrand, 0.0, rz_span, params; dt = self.Gr.dz)
    rsol = solve(rprob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_0 = [rsol(self.Gr.z_half[k]) for k in face_indicies(self.Gr)]

    prob = ODEProblem(integrand, 0.0, z_span, params; dt = self.Gr.dz)
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_1 = [sol(self.Gr.z_half[k]) for k in face_indicies(self.Gr)]
    self.f_rad .= self.F0 .* exp.(-q_0)
    self.f_rad .+= self.F1 .* exp.(-q_1)

    # cooling in free troposphere
    @inbounds for k in face_indicies(self.Gr)
        if self.Gr.z[k] > zi
            cbrt_z = cbrt(self.Gr.z[k] - zi)
            self.f_rad[k] += ρ_i * dycoms_cp * self.divergence * self.alpha_z * (cbrt_z^4 / 4 + zi * cbrt_z)
        end
    end

    @inbounds for k in real_center_indicies(self.Gr)
        self.dTdt[k] = -(self.f_rad[k + 1] - self.f_rad[k]) / self.Gr.dz / self.Ref.rho0_half[k] / dycoms_cp
    end

    return
end

function update(self::RadiationBase{RadiationDYCOMS_RF01}, GMV::GridMeanVariables)

    calculate_radiation(self, GMV)

    @inbounds for k in real_center_indicies(self.Gr)
        # Apply large-scale horizontal advection tendencies
        pp = TD.PhasePartition(GMV.QT.values[k], GMV.QL.values[k], 0.0)
        qv = TD.vapor_specific_humidity(pp)
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
