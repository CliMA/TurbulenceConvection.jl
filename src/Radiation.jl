
function initialize(self::RadiationBase, grid, ::RadiationBaseType)
    self.dTdt = center_field(grid)
    self.dqtdt = center_field(grid)
    return
end

update(self::RadiationBase, grid, state, GMV::GridMeanVariables, param_set) = nothing

initialize_io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing
io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing

initialize(self::RadiationBase{RadiationNone}, grid) = initialize(self, grid, RadiationBaseType())

function initialize(self::RadiationBase{RadiationDYCOMS_RF01}, grid)
    initialize(self, grid, RadiationBaseType())


    self.divergence = 3.75e-6  # divergence is defined twice: here and in initialize_forcing method of DYCOMS_RF01 case class
    # where it is used to initialize large scale subsidence
    self.alpha_z = 1.0
    self.kappa = 85.0
    self.F0 = 70.0
    self.F1 = 22.0
    self.divergence = 3.75e-6
    self.f_rad = face_field(grid)
    return
end

"""
see eq. 3 in Stevens et. al. 2005 DYCOMS paper
"""
function calculate_radiation(self::RadiationBase{RadiationDYCOMS_RF01}, grid, state, GMV::GridMeanVariables, param_set)
    cp_d = CPP.cp_d(param_set)
    ρ0_f = face_ref_state(state).ρ0
    ρ0_c = center_ref_state(state).ρ0
    # find zi (level of 8.0 g/kg isoline of qt)
    # TODO: report bug: zi and ρ_i are not initialized
    zi = 0
    ρ_i = 0
    kc_surf = kc_surface(grid)
    q_tot_surf = GMV.QT.values[kc_surf]
    @inbounds for k in real_face_indices(grid)
        q_tot_f = interpc2f(GMV.QT.values, grid, k; bottom = SetValue(q_tot_surf), top = SetGradient(0))
        if (q_tot_f < 8.0 / 1000)
            idx_zi = k
            # will be used at cell faces
            zi = grid.zf[k]
            ρ_i = ρ0_f[k]
            break
        end
    end

    ρ_z = Dierckx.Spline1D(vec(grid.zc), vec(ρ0_c); k = 1)
    q_liq_z = Dierckx.Spline1D(vec(grid.zc), vec(GMV.QL.values); k = 1)

    integrand(ρq_l, params, z) = params.κ * ρ_z(z) * q_liq_z(z)
    rintegrand(ρq_l, params, z) = -integrand(ρq_l, params, z)

    z_span = (grid.zmin, grid.zmax)
    rz_span = (grid.zmax, grid.zmin)
    params = (; κ = self.kappa)

    rprob = ODE.ODEProblem(rintegrand, 0.0, rz_span, params; dt = grid.Δz)
    rsol = ODE.solve(rprob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_0 = rsol.(vec(grid.zf))

    prob = ODE.ODEProblem(integrand, 0.0, z_span, params; dt = grid.Δz)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_1 = sol.(vec(grid.zf))
    self.f_rad .= self.F0 .* exp.(-q_0)
    self.f_rad .+= self.F1 .* exp.(-q_1)

    # cooling in free troposphere
    @inbounds for k in real_face_indices(grid)
        if grid.zf[k] > zi
            cbrt_z = cbrt(grid.zf[k] - zi)
            self.f_rad[k] += ρ_i * cp_d * self.divergence * self.alpha_z * (cbrt_z^4 / 4 + zi * cbrt_z)
        end
    end

    @inbounds for k in real_center_indices(grid)
        f_rad_dual = dual_faces(self.f_rad, grid, k)
        ∇f_rad = ∇_staggered(f_rad_dual, grid)
        self.dTdt[k] = -∇f_rad / ρ0_c[k] / cp_d
    end

    return
end

update(self::RadiationBase{RadiationDYCOMS_RF01}, grid, state, GMV::GridMeanVariables, param_set) =
    calculate_radiation(self, grid, state, GMV, param_set)

function initialize_io(self::RadiationBase{RadiationDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    add_profile(Stats, "rad_dTdt")
    add_profile(Stats, "rad_flux")
    return
end

function io(self::RadiationBase{RadiationDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    write_profile(Stats, "rad_dTdt", self.dTdt)
    write_profile(Stats, "rad_flux", self.f_rad)
    return
end

function initialize(self::RadiationBase{RadiationLES}, grid, LESDat::LESData)
    initialize(self, grid, RadiationBaseType())
    # load from LES
    NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        # interpolate here
        zc_les = get_nc_data(data, "zc")
        meandata = mean_nc_data(data, "profiles", "dtdt_rad", imin, imax)
        self.dTdt = pyinterp(grid.zc, zc_les, meandata)
    end
    return
end

function initialize_io(self::RadiationBase{RadiationLES}, Stats::NetCDFIO_Stats)
    add_profile(Stats, "rad_dTdt")
    add_profile(Stats, "rad_flux")
    return
end

function io(self::RadiationBase{RadiationLES}, Stats::NetCDFIO_Stats)
    write_profile(Stats, "rad_dTdt", self.dTdt)
    write_profile(Stats, "rad_flux", self.f_rad)
    return
end
