
function initialize(self::EDMF_PrognosticTKE, Case::CasesBase, GMV::GridMeanVariables, Ref::ReferenceState)
    if Case.casename == "DryBubble"
        initialize_DryBubble(self.UpdVar, GMV, Ref)
    else
        initialize(self.UpdVar, GMV)
    end
    return
end

# Initialize the IO pertaining to this class
function initialize_io(self::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)

    initialize_io(self.UpdVar, Stats)
    initialize_io(self.EnvVar, Stats)
    initialize_io(self.Rain, Stats)

    add_profile(Stats, "eddy_viscosity")
    add_profile(Stats, "eddy_diffusivity")
    add_profile(Stats, "entrainment_sc")
    add_profile(Stats, "detrainment_sc")
    add_profile(Stats, "nh_pressure")
    add_profile(Stats, "nh_pressure_adv")
    add_profile(Stats, "nh_pressure_drag")
    add_profile(Stats, "nh_pressure_b")
    add_profile(Stats, "asp_ratio")
    add_profile(Stats, "b_coeff")

    add_profile(Stats, "horizontal_KM")
    add_profile(Stats, "horizontal_KH")
    add_profile(Stats, "sorting_function")
    add_profile(Stats, "b_mix")
    add_ts(Stats, "rd")
    add_profile(Stats, "turbulent_entrainment")
    add_profile(Stats, "turbulent_entrainment_full")
    add_profile(Stats, "turbulent_entrainment_W")
    add_profile(Stats, "turbulent_entrainment_H")
    add_profile(Stats, "turbulent_entrainment_QT")
    add_profile(Stats, "massflux")
    add_profile(Stats, "massflux_h")
    add_profile(Stats, "massflux_qt")
    add_profile(Stats, "massflux_tendency_h")
    add_profile(Stats, "massflux_tendency_qt")
    add_profile(Stats, "diffusive_flux_h")
    add_profile(Stats, "diffusive_flux_u")
    add_profile(Stats, "diffusive_flux_v")
    add_profile(Stats, "diffusive_flux_qt")
    add_profile(Stats, "diffusive_tendency_h")
    add_profile(Stats, "diffusive_tendency_qt")
    add_profile(Stats, "total_flux_h")
    add_profile(Stats, "total_flux_qt")
    add_profile(Stats, "mixing_length")
    add_profile(Stats, "updraft_qt_precip")
    add_profile(Stats, "updraft_thetal_precip")
    # Diff mixing lengths: Ignacio
    add_profile(Stats, "ed_length_scheme")
    add_profile(Stats, "mixing_length_ratio")
    add_profile(Stats, "entdet_balance_length")
    add_profile(Stats, "interdomain_tke_t")
    if self.calc_tke
        add_profile(Stats, "tke_buoy")
        add_profile(Stats, "tke_dissipation")
        add_profile(Stats, "tke_entr_gain")
        add_profile(Stats, "tke_detr_loss")
        add_profile(Stats, "tke_shear")
        add_profile(Stats, "tke_pressure")
        add_profile(Stats, "tke_interdomain")
        add_profile(Stats, "tke_transport")
        add_profile(Stats, "tke_advection")
    end

    if self.calc_scalar_var
        add_profile(Stats, "Hvar_dissipation")
        add_profile(Stats, "QTvar_dissipation")
        add_profile(Stats, "HQTcov_dissipation")
        add_profile(Stats, "Hvar_entr_gain")
        add_profile(Stats, "QTvar_entr_gain")
        add_profile(Stats, "Hvar_detr_loss")
        add_profile(Stats, "QTvar_detr_loss")
        add_profile(Stats, "HQTcov_detr_loss")
        add_profile(Stats, "HQTcov_entr_gain")
        add_profile(Stats, "Hvar_shear")
        add_profile(Stats, "QTvar_shear")
        add_profile(Stats, "HQTcov_shear")
        add_profile(Stats, "Hvar_rain")
        add_profile(Stats, "QTvar_rain")
        add_profile(Stats, "HQTcov_rain")
        add_profile(Stats, "Hvar_interdomain")
        add_profile(Stats, "QTvar_interdomain")
        add_profile(Stats, "HQTcov_interdomain")
    end
    return
end

function io(self::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats, TS::TimeStepping)

    grid = get_grid(self)
    ref_state = reference_state(self)
    cinterior = grid.cinterior
    finterior = grid.finterior

    mean_entr_sc = center_field(grid)
    mean_nh_pressure = center_field(grid)
    mean_nh_pressure_adv = center_field(grid)
    mean_nh_pressure_drag = center_field(grid)
    mean_nh_pressure_b = center_field(grid)
    mean_asp_ratio = center_field(grid)
    mean_b_coeff = center_field(grid)

    mean_detr_sc = center_field(grid)
    massflux = face_field(grid)
    mf_h = face_field(grid)
    mf_qt = face_field(grid)
    mean_frac_turb_entr = center_field(grid)
    mean_frac_turb_entr_full = center_field(grid)
    mean_turb_entr_W = center_field(grid)
    mean_turb_entr_H = center_field(grid)
    mean_turb_entr_QT = center_field(grid)
    mean_horizontal_KM = center_field(grid)
    mean_horizontal_KH = center_field(grid)
    mean_sorting_function = center_field(grid)
    mean_b_mix = center_field(grid)

    io(self.UpdVar, Stats, ref_state)
    io(self.EnvVar, Stats, ref_state)
    io(self.Rain, Stats, ref_state, self.UpdThermo, self.EnvThermo, TS)

    write_profile(Stats, "eddy_viscosity", diffusivity_m(self).values[cinterior])
    write_profile(Stats, "eddy_diffusivity", diffusivity_h(self).values[cinterior])
    write_ts(Stats, "rd", mean(self.pressure_plume_spacing))
    @inbounds for k in real_center_indicies(grid)
        mf_h[k] = interp2pt(self.massflux_h[k], self.massflux_h[k - 1])
        mf_qt[k] = interp2pt(self.massflux_qt[k], self.massflux_qt[k - 1])
        if self.UpdVar.Area.bulkvalues[k] > 0.0
            @inbounds for i in xrange(self.n_updrafts)
                massflux[k] += interp2pt(self.m[i, k], self.m[i, k - 1])
                mean_entr_sc[k] += self.UpdVar.Area.values[i, k] * self.entr_sc[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_detr_sc[k] += self.UpdVar.Area.values[i, k] * self.detr_sc[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_nh_pressure[k] +=
                    self.UpdVar.Area.values[i, k] * self.nh_pressure[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_nh_pressure_b[k] +=
                    self.UpdVar.Area.values[i, k] * self.nh_pressure_b[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_nh_pressure_adv[k] +=
                    self.UpdVar.Area.values[i, k] * self.nh_pressure_adv[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_nh_pressure_drag[k] +=
                    self.UpdVar.Area.values[i, k] * self.nh_pressure_drag[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_asp_ratio[k] +=
                    self.UpdVar.Area.values[i, k] * self.asp_ratio[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_b_coeff[k] += self.UpdVar.Area.values[i, k] * self.b_coeff[i, k] / self.UpdVar.Area.bulkvalues[k]

                mean_frac_turb_entr_full[k] +=
                    self.UpdVar.Area.values[i, k] * self.frac_turb_entr_full[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_frac_turb_entr[k] +=
                    self.UpdVar.Area.values[i, k] * self.frac_turb_entr[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_turb_entr_W[k] +=
                    self.UpdVar.Area.values[i, k] * self.turb_entr_W[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_turb_entr_H[k] +=
                    self.UpdVar.Area.values[i, k] * self.turb_entr_H[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_turb_entr_QT[k] +=
                    self.UpdVar.Area.values[i, k] * self.turb_entr_QT[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_horizontal_KM[k] +=
                    self.UpdVar.Area.values[i, k] * self.horizontal_KM[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_horizontal_KH[k] +=
                    self.UpdVar.Area.values[i, k] * self.horizontal_KH[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_sorting_function[k] +=
                    self.UpdVar.Area.values[i, k] * self.sorting_function[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_b_mix[k] += self.UpdVar.Area.values[i, k] * self.b_mix[i, k] / self.UpdVar.Area.bulkvalues[k]
            end
        end
    end

    write_profile(Stats, "turbulent_entrainment", mean_frac_turb_entr[cinterior])
    write_profile(Stats, "turbulent_entrainment_full", mean_frac_turb_entr_full[cinterior])
    write_profile(Stats, "turbulent_entrainment_W", mean_turb_entr_W[cinterior])
    write_profile(Stats, "turbulent_entrainment_H", mean_turb_entr_H[cinterior])
    write_profile(Stats, "turbulent_entrainment_QT", mean_turb_entr_QT[cinterior])
    write_profile(Stats, "horizontal_KM", mean_horizontal_KM[cinterior])
    write_profile(Stats, "horizontal_KH", mean_horizontal_KH[cinterior])
    write_profile(Stats, "entrainment_sc", mean_entr_sc[cinterior])
    write_profile(Stats, "detrainment_sc", mean_detr_sc[cinterior])
    write_profile(Stats, "sorting_function", mean_sorting_function[cinterior])
    write_profile(Stats, "b_mix", mean_b_mix[cinterior])
    write_profile(Stats, "nh_pressure", mean_nh_pressure[cinterior])
    write_profile(Stats, "nh_pressure_adv", mean_nh_pressure_adv[cinterior])
    write_profile(Stats, "nh_pressure_drag", mean_nh_pressure_drag[cinterior])
    write_profile(Stats, "nh_pressure_b", mean_nh_pressure_b[cinterior])
    write_profile(Stats, "asp_ratio", mean_asp_ratio[cinterior])
    write_profile(Stats, "b_coeff", mean_b_coeff[cinterior])

    write_profile(Stats, "massflux", massflux[cinterior])
    write_profile(Stats, "massflux_h", mf_h[cinterior])
    write_profile(Stats, "massflux_qt", mf_qt[cinterior])
    write_profile(Stats, "massflux_tendency_h", self.massflux_tendency_h[cinterior])
    write_profile(Stats, "massflux_tendency_qt", self.massflux_tendency_qt[cinterior])
    write_profile(Stats, "diffusive_flux_h", self.diffusive_flux_h[cinterior])
    write_profile(Stats, "diffusive_flux_qt", self.diffusive_flux_qt[cinterior])
    write_profile(Stats, "diffusive_flux_u", self.diffusive_flux_u[cinterior])
    write_profile(Stats, "diffusive_flux_v", self.diffusive_flux_v[cinterior])
    write_profile(Stats, "diffusive_tendency_h", self.diffusive_tendency_h[cinterior])
    write_profile(Stats, "diffusive_tendency_qt", self.diffusive_tendency_qt[cinterior])
    write_profile(Stats, "total_flux_h", mf_h[cinterior] .+ self.diffusive_flux_h[cinterior])
    write_profile(Stats, "total_flux_qt", mf_qt[cinterior] .+ self.diffusive_flux_qt[cinterior])
    write_profile(Stats, "mixing_length", self.mixing_length[cinterior])
    write_profile(Stats, "updraft_qt_precip", self.UpdThermo.prec_source_qt_tot[cinterior])
    write_profile(Stats, "updraft_thetal_precip", self.UpdThermo.prec_source_h_tot[cinterior])

    #Different mixing lengths : Ignacio
    write_profile(Stats, "ed_length_scheme", self.mls[cinterior])
    write_profile(Stats, "mixing_length_ratio", self.ml_ratio[cinterior])
    write_profile(Stats, "entdet_balance_length", self.l_entdet[cinterior])
    write_profile(Stats, "interdomain_tke_t", self.b[cinterior])
    if self.calc_tke
        compute_covariance_dissipation(self, self.EnvVar.TKE)
        write_profile(Stats, "tke_dissipation", self.EnvVar.TKE.dissipation[cinterior])
        write_profile(Stats, "tke_entr_gain", self.EnvVar.TKE.entr_gain[cinterior])
        compute_covariance_detr(self, self.EnvVar.TKE)
        write_profile(Stats, "tke_detr_loss", self.EnvVar.TKE.detr_loss[cinterior])
        write_profile(Stats, "tke_shear", self.EnvVar.TKE.shear[cinterior])
        write_profile(Stats, "tke_buoy", self.EnvVar.TKE.buoy[cinterior])
        write_profile(Stats, "tke_pressure", self.EnvVar.TKE.press[cinterior])
        write_profile(Stats, "tke_interdomain", self.EnvVar.TKE.interdomain[cinterior])
        compute_tke_transport(self)
        write_profile(Stats, "tke_transport", self.tke_transport[cinterior])
        compute_tke_advection(self)
        write_profile(Stats, "tke_advection", self.tke_advection[cinterior])
    end

    if self.calc_scalar_var
        compute_covariance_dissipation(self, self.EnvVar.Hvar)
        write_profile(Stats, "Hvar_dissipation", self.EnvVar.Hvar.dissipation[cinterior])
        compute_covariance_dissipation(self, self.EnvVar.QTvar)
        write_profile(Stats, "QTvar_dissipation", self.EnvVar.QTvar.dissipation[cinterior])
        compute_covariance_dissipation(self, self.EnvVar.HQTcov)
        write_profile(Stats, "HQTcov_dissipation", self.EnvVar.HQTcov.dissipation[cinterior])
        write_profile(Stats, "Hvar_entr_gain", self.EnvVar.Hvar.entr_gain[cinterior])
        write_profile(Stats, "QTvar_entr_gain", self.EnvVar.QTvar.entr_gain[cinterior])
        write_profile(Stats, "HQTcov_entr_gain", self.EnvVar.HQTcov.entr_gain[cinterior])
        compute_covariance_detr(self, self.EnvVar.Hvar)
        compute_covariance_detr(self, self.EnvVar.QTvar)
        compute_covariance_detr(self, self.EnvVar.HQTcov)
        write_profile(Stats, "Hvar_detr_loss", self.EnvVar.Hvar.detr_loss[cinterior])
        write_profile(Stats, "QTvar_detr_loss", self.EnvVar.QTvar.detr_loss[cinterior])
        write_profile(Stats, "HQTcov_detr_loss", self.EnvVar.HQTcov.detr_loss[cinterior])
        write_profile(Stats, "Hvar_shear", self.EnvVar.Hvar.shear[cinterior])
        write_profile(Stats, "QTvar_shear", self.EnvVar.QTvar.shear[cinterior])
        write_profile(Stats, "HQTcov_shear", self.EnvVar.HQTcov.shear[cinterior])
        write_profile(Stats, "Hvar_rain", self.EnvVar.Hvar.rain_src[cinterior])
        write_profile(Stats, "QTvar_rain", self.EnvVar.QTvar.rain_src[cinterior])
        write_profile(Stats, "HQTcov_rain", self.EnvVar.HQTcov.rain_src[cinterior])
        write_profile(Stats, "Hvar_interdomain", self.EnvVar.Hvar.interdomain[cinterior])
        write_profile(Stats, "QTvar_interdomain", self.EnvVar.QTvar.interdomain[cinterior])
        write_profile(Stats, "HQTcov_interdomain", self.EnvVar.HQTcov.interdomain[cinterior])
    end
    return
end

# Perform the update of the scheme
function update(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    grid = get_grid(self)
    update_inversion(self, GMV, Case.inversion_option)
    compute_pressure_plume_spacing(self, GMV, Case)
    self.wstar = get_wstar(Case.Sur.bflux, self.base.zi)
    if TS.nstep == 0
        decompose_environment(self, GMV, "values")

        if Case.casename == "DryBubble"
            saturation_adjustment(self.EnvThermo, self.EnvVar)
            buoyancy(self.UpdThermo, self.UpdVar, self.EnvVar, GMV, self.extrapolate_buoyancy)
        end

        microphysics(self.EnvThermo, self.EnvVar, self.Rain, TS.dt)
        initialize_covariance(self, GMV, Case)

        @inbounds for k in center_indicies(grid)
            if self.calc_tke
                self.EnvVar.TKE.values[k] = GMV.TKE.values[k]
            end
            if self.calc_scalar_var
                self.EnvVar.Hvar.values[k] = GMV.Hvar.values[k]
                self.EnvVar.QTvar.values[k] = GMV.QTvar.values[k]
                self.EnvVar.HQTcov.values[k] = GMV.HQTcov.values[k]
            end
        end
    end

    decompose_environment(self, GMV, "values")
    compute_prognostic_updrafts(self, GMV, Case, TS)

    # TODO -maybe not needed? - both diagnostic and prognostic updrafts end with decompose_environment
    # But in general ok here without thermodynamics because MF doesnt depend directly on buoyancy
    decompose_environment(self, GMV, "values")
    update_GMV_MF(self, GMV, TS)
    # (###)
    # decompose_environment +  EnvThermo.saturation_adjustment + UpdThermo.buoyancy should always be used together
    # This ensures that
    #   - the buoyancy of updrafts and environment is up to date with the most recent decomposition,
    #   - the buoyancy of updrafts and environment is updated such that
    #     the mean buoyancy with repect to reference state alpha_0 is zero.

    decompose_environment(self, GMV, "mf_update")
    microphysics(self.EnvThermo, self.EnvVar, self.Rain, TS.dt) # saturation adjustment + rain creation
    # Sink of environmental QT and H due to rain creation is applied in tridiagonal solver
    buoyancy(self.UpdThermo, self.UpdVar, self.EnvVar, GMV, self.extrapolate_buoyancy)
    compute_eddy_diffusivities_tke(self, GMV, Case)
    update_GMV_ED(self, GMV, Case, TS)
    compute_covariance(self, GMV, Case, TS)

    if self.Rain.rain_model == "clima_1m"
        # sum updraft and environment rain into bulk rain
        sum_subdomains_rain(self.Rain, self.UpdThermo, self.EnvThermo)

        # rain fall (all three categories are assumed to be falling though "grid-mean" conditions
        solve_rain_fall(self.rainphysics, self.Rain, GMV, TS, self.Rain.QR, self.Rain.RainArea)
        solve_rain_fall(self.rainphysics, self.Rain, GMV, TS, self.Rain.Upd_QR, self.Rain.Upd_RainArea)
        solve_rain_fall(self.rainphysics, self.Rain, GMV, TS, self.Rain.Env_QR, self.Rain.Env_RainArea)

        # rain evaporation (all three categories are assumed to be evaporating in "grid-mean" conditions
        solve_rain_evap(self.rainphysics, self.Rain, GMV, TS, self.Rain.QR, self.Rain.RainArea)
        solve_rain_evap(self.rainphysics, self.Rain, GMV, TS, self.Rain.Upd_QR, self.Rain.Upd_RainArea)
        solve_rain_evap(self.rainphysics, self.Rain, GMV, TS, self.Rain.Env_QR, self.Rain.Env_RainArea)
    end
    # update grid-mean cloud fraction and cloud cover
    @inbounds for k in center_indicies(grid)
        self.EnvVar.Area.values[k] = 1.0 - self.UpdVar.Area.bulkvalues[k]
        GMV.cloud_fraction.values[k] =
            self.EnvVar.Area.values[k] * self.EnvVar.cloud_fraction.values[k] +
            self.UpdVar.Area.bulkvalues[k] * self.UpdVar.cloud_fraction[k]
    end
    GMV.cloud_cover = min(self.EnvVar.cloud_cover + sum(self.UpdVar.cloud_cover), 1)
    # Back out the tendencies of the grid mean variables for the whole timestep
    # by differencing GMV.new and GMV.values
    update(self.base, GMV, Case, TS)
    return
end

function compute_prognostic_updrafts(
    self::EDMF_PrognosticTKE,
    GMV::GridMeanVariables,
    Case::CasesBase,
    TS::TimeStepping,
)

    time_elapsed = 0.0

    set_subdomain_bcs(self)
    set_new_with_values(self.UpdVar)
    set_old_with_values(self.UpdVar)

    set_updraft_surface_bc(self, GMV, Case)
    self.dt_upd = min(TS.dt, 0.5 * get_grid(self).dz / fmax(maximum(self.UpdVar.W.values), 1e-10))

    clear_precip_sources(self.UpdThermo)

    while time_elapsed < TS.dt
        compute_entrainment_detrainment(self, GMV, Case)
        if self.turbulent_entrainment_factor > 1.0e-6
            compute_horizontal_eddy_diffusivities(self, GMV)
            compute_turbulent_entrainment(self, GMV, Case)
        end
        compute_nh_pressure(self)
        solve_updraft_velocity_area(self)
        solve_updraft_scalars(self, GMV)
        microphysics(self.UpdThermo, self.UpdVar, self.Rain, TS.dt) # causes division error in dry bubble first time step

        set_values_with_new(self.UpdVar)
        zero_area_fraction_cleanup(self, GMV)
        time_elapsed += self.dt_upd
        self.dt_upd = min(TS.dt - time_elapsed, 0.5 * get_grid(self).dz / fmax(maximum(self.UpdVar.W.values), 1e-10))
        # (####)
        # TODO - see comment (###)
        # It would be better to have a simple linear rule for updating environment here
        # instead of calling EnvThermo saturation adjustment scheme for every updraft.
        decompose_environment(self, GMV, "values")
        saturation_adjustment(self.EnvThermo, self.EnvVar)
        buoyancy(self.UpdThermo, self.UpdVar, self.EnvVar, GMV, self.extrapolate_buoyancy)
        set_subdomain_bcs(self)
    end

    update_total_precip_sources(self.UpdThermo)
    return
end

function update_inversion(self, GMV::GridMeanVariables, option)
    update_inversion(self.base, GMV, option)
    return
end

function compute_mixing_length(self, obukhov_length, ustar, GMV::GridMeanVariables)

    grid = get_grid(self)
    ref_state = reference_state(self)
    gw = grid.gw
    dzi = grid.dzi
    tau = get_mixing_tau(self.base.zi, self.wstar)
    l = pyzeros(3)
    m_eps = 1.0e-9 # Epsilon to avoid zero
    @inbounds for k in real_center_indicies(grid)
        z_ = grid.z_half[k]
        # kz scale (surface layer)
        if obukhov_length < 0.0 #unstable
            l2 =
                vkb * z_ / (sqrt(self.EnvVar.TKE.values[grid.gw] / ustar / ustar) * self.tke_ed_coeff) *
                fmin((1.0 - 100.0 * z_ / obukhov_length)^0.2, 1.0 / vkb)
        else # neutral or stable
            l2 = vkb * z_ / (sqrt(self.EnvVar.TKE.values[grid.gw] / ustar / ustar) * self.tke_ed_coeff)
        end

        # Buoyancy-shear-subdomain exchange-dissipation TKE equilibrium scale
        U_cut = cut(GMV.U.values, k)
        V_cut = cut(GMV.V.values, k)
        shear2 =
            pow(∇_collocated(U_cut, grid), 2) +
            pow(∇_collocated(V_cut, grid), 2) +
            pow(∇f2c(self.EnvVar.W.values, grid, k), 2)
        qt_dry = self.EnvThermo.qt_dry[k]
        th_dry = self.EnvThermo.th_dry[k]
        t_cloudy = self.EnvThermo.t_cloudy[k]
        qv_cloudy = self.EnvThermo.qv_cloudy[k]
        qt_cloudy = self.EnvThermo.qt_cloudy[k]
        th_cloudy = self.EnvThermo.th_cloudy[k]
        lh = latent_heat(t_cloudy)
        cpm = cpm_c(qt_cloudy)

        QT_cut = cut(self.EnvVar.QT.values, k)
        grad_qt = c∇(QT_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        THL_cut = cut(self.EnvVar.THL.values, k)
        grad_thl = c∇(THL_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))

        # g/theta_ref
        prefactor = g * (Rd / ref_state.alpha0_half[k] / ref_state.p0_half[k]) * exner_c(ref_state.p0_half[k])
        d_buoy_thetal_dry = prefactor * (1.0 + (eps_vi - 1.0) * qt_dry)
        d_buoy_qt_dry = prefactor * th_dry * (eps_vi - 1.0)

        if self.EnvVar.cloud_fraction.values[k] > 0.0
            d_buoy_thetal_cloudy = (
                prefactor * (1.0 + eps_vi * (1.0 + lh / Rv / t_cloudy) * qv_cloudy - qt_cloudy) /
                (1.0 + lh * lh / cpm / Rv / t_cloudy / t_cloudy * qv_cloudy)
            )
            d_buoy_qt_cloudy = (lh / cpm / t_cloudy * d_buoy_thetal_cloudy - prefactor) * th_cloudy
        else
            d_buoy_thetal_cloudy = 0.0
            d_buoy_qt_cloudy = 0.0
        end
        d_buoy_thetal_total = (
            self.EnvVar.cloud_fraction.values[k] * d_buoy_thetal_cloudy +
            (1.0 - self.EnvVar.cloud_fraction.values[k]) * d_buoy_thetal_dry
        )
        d_buoy_qt_total = (
            self.EnvVar.cloud_fraction.values[k] * d_buoy_qt_cloudy +
            (1.0 - self.EnvVar.cloud_fraction.values[k]) * d_buoy_qt_dry
        )

        # Partial buoyancy gradients
        grad_b_thl = grad_thl * d_buoy_thetal_total
        grad_b_qt = grad_qt * d_buoy_qt_total
        ri_grad = fmin(grad_b_thl / fmax(shear2, m_eps) + grad_b_qt / fmax(shear2, m_eps), 0.25)

        # Turbulent Prandtl number
        if obukhov_length > 0.0 && ri_grad > 0.0 #stable
            # CSB (Dan Li, 2019), with Pr_neutral=0.74 and w1=40.0/13.0
            self.prandtl_nvec[k] =
                prandtl_number(self) * (
                    2.0 * ri_grad /
                    (1.0 + (53.0 / 13.0) * ri_grad - sqrt((1.0 + (53.0 / 13.0) * ri_grad)^2.0 - 4.0 * ri_grad))
                )
        else
            self.prandtl_nvec[k] = prandtl_number(self)
        end

        # Production/destruction terms
        a =
            self.tke_ed_coeff *
            (shear2 - grad_b_thl / self.prandtl_nvec[k] - grad_b_qt / self.prandtl_nvec[k]) *
            sqrt(self.EnvVar.TKE.values[k])
        # Dissipation term
        c_neg = self.tke_diss_coeff * self.EnvVar.TKE.values[k] * sqrt(self.EnvVar.TKE.values[k])
        # Subdomain exchange term
        self.b[k] = 0.0
        for nn in xrange(self.n_updrafts)
            wc_upd_nn = (self.UpdVar.W.values[nn, k] + self.UpdVar.W.values[nn, k - 1]) / 2.0
            wc_env = (self.EnvVar.W.values[k] + self.EnvVar.W.values[k - 1]) / 2.0
            frac_turb_entr_half = interp2pt(self.frac_turb_entr_full[nn, k], self.frac_turb_entr_full[nn, k - 1])
            self.b[k] +=
                self.UpdVar.Area.values[nn, k] * wc_upd_nn * self.detr_sc[nn, k] /
                (1.0 - self.UpdVar.Area.bulkvalues[k]) *
                ((wc_upd_nn - wc_env) * (wc_upd_nn - wc_env) / 2.0 - self.EnvVar.TKE.values[k]) -
                self.UpdVar.Area.values[nn, k] * wc_upd_nn * (wc_upd_nn - wc_env) * frac_turb_entr_half * wc_env /
                (1.0 - self.UpdVar.Area.bulkvalues[k])
        end

        if abs(a) > m_eps && 4.0 * a * c_neg > -self.b[k] * self.b[k]
            self.l_entdet[k] = fmax(-self.b[k] / 2.0 / a + sqrt(self.b[k] * self.b[k] + 4.0 * a * c_neg) / 2.0 / a, 0.0)
        elseif abs(a) < m_eps && abs(self.b[k]) > m_eps
            self.l_entdet[k] = c_neg / self.b[k]
        end
        l3 = self.l_entdet[k]

        # Limiting stratification scale (Deardorff, 1976)
        p0_cut = cut(ref_state.p0_half, k)
        T_cut = cut(self.EnvVar.T.values, k)
        QT_cut = cut(self.EnvVar.QT.values, k)
        QL_cut = cut(self.EnvVar.QL.values, k)
        thv_cut = theta_virt_c.(p0_cut, T_cut, QT_cut, QL_cut)
        thv = thv_cut[2]
        grad_thv = c∇(thv_cut, grid, k; bottom = SetGradient(0), top = Extrapolate())

        # Effective static stability using environmental mean.
        # Set lambda for now to environmental cloud_fraction (TBD: Rain)
        grad_th_eff =
            (1.0 - self.EnvVar.cloud_fraction.values[k]) * grad_thv +
            self.EnvVar.cloud_fraction.values[k] * (
                1.0 / exp(
                    -latent_heat(self.EnvVar.T.values[k]) * self.EnvVar.QL.values[k] / cpm_c(self.EnvVar.QT.values[k]) /
                    self.EnvVar.T.values[k],
                ) * (
                    (1.0 + (eps_vi - 1.0) * self.EnvVar.QT.values[k]) * grad_thl +
                    (eps_vi - 1.0) * self.EnvVar.THL.values[k] * grad_qt
                )
            )
        N = sqrt(fmax(g / thv * grad_th_eff, 0.0))
        if N > 0.0
            l1 = fmin(sqrt(fmax(self.static_stab_coeff * self.EnvVar.TKE.values[k], 0.0)) / N, 1.0e6)
        else
            l1 = 1.0e6
        end

        l[0] = l1
        l[1] = l3
        l[2] = l2

        j = 0
        while (j < length(l))
            if l[j] < m_eps || l[j] > 1.0e6
                l[j] = 1.0e6
            end
            j += 1
        end

        self.mls[k] = argmin(l)
        self.mixing_length[k] = lamb_smooth_minimum(l, 0.1, 1.5)
        self.ml_ratio[k] = self.mixing_length[k] / l[Int(self.mls[k])]
    end
    return
end


function compute_eddy_diffusivities_tke(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    grid = get_grid(self)
    if self.similarity_diffusivity
        compute_eddy_diffusivities_similarity(self.base, GMV, Case)
    else
        compute_mixing_length(self, Case.Sur.obukhov_length, Case.Sur.ustar, GMV)
        KM = diffusivity_m(self)
        KH = diffusivity_h(self)
        @inbounds for k in real_center_indicies(grid)
            lm = self.mixing_length[k]
            pr = self.prandtl_nvec[k]
            KM.values[k] = self.tke_ed_coeff * lm * sqrt(fmax(self.EnvVar.TKE.values[k], 0.0))
            KH.values[k] = KM.values[k] / pr
        end
    end
    return
end

function compute_horizontal_eddy_diffusivities(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)
    grid = get_grid(self)
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues
    l = pyzeros(2)

    @inbounds for k in real_center_indicies(grid)
        @inbounds for i in xrange(self.n_updrafts)
            if self.UpdVar.Area.values[i, k] > 0.0
                self.horizontal_KM[i, k] =
                    self.UpdVar.Area.values[i, k] *
                    self.turbulent_entrainment_factor *
                    sqrt(fmax(self.EnvVar.TKE.values[k], 0.0)) *
                    self.pressure_plume_spacing[i]
                self.horizontal_KH[i, k] = self.horizontal_KM[i, k] / self.prandtl_nvec[k]
            else
                self.horizontal_KM[i, k] = 0.0
                self.horizontal_KH[i, k] = 0.0
            end
        end
    end

    return
end

function set_updraft_surface_bc(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)

    update_inversion(self, GMV, Case.inversion_option)
    self.wstar = get_wstar(Case.Sur.bflux, self.base.zi)

    gw = get_grid(self).gw
    dzi = get_grid(self).dzi
    zLL = get_grid(self).z_half[gw]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    alpha0LL = reference_state(self).alpha0_half[gw]
    qt_var = get_surface_variance(Case.Sur.rho_qtflux * alpha0LL, Case.Sur.rho_qtflux * alpha0LL, ustar, zLL, oblength)
    h_var = get_surface_variance(Case.Sur.rho_hflux * alpha0LL, Case.Sur.rho_hflux * alpha0LL, ustar, zLL, oblength)

    if Case.Sur.bflux > 0.0
        a_total = self.surface_area
        self.entr_surface_bc = 2.0 * dzi
        self.detr_surface_bc = 0.0
    else
        # a_total = self.surface_area
        a_total = self.minimum_area * 0.9
        self.entr_surface_bc = 0.0
        self.detr_surface_bc = 2.0 * dzi
    end

    a_ = a_total / self.n_updrafts
    @inbounds for i in xrange(self.n_updrafts)
        surface_scalar_coeff = percentile_bounds_mean_norm(1.0 - a_total + i * a_, 1.0 - a_total + (i + 1) * a_, 1000)
        self.area_surface_bc[i] = a_
        self.w_surface_bc[i] = 0.0
        self.h_surface_bc[i] = (GMV.H.values[gw] + surface_scalar_coeff * sqrt(h_var))
        self.qt_surface_bc[i] = (GMV.QT.values[gw] + surface_scalar_coeff * sqrt(qt_var))
    end
    return
end

function reset_surface_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    flux1 = Case.Sur.rho_hflux
    flux2 = Case.Sur.rho_qtflux
    grid = get_grid(self)
    ref_state = reference_state(self)
    zLL = get_grid(self).z_half[grid.gw]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    alpha0LL = ref_state.alpha0_half[grid.gw]

    if self.calc_tke
        self.EnvVar.TKE.values[grid.gw] =
            get_surface_tke(Case.Sur.ustar, self.wstar, grid.z_half[grid.gw], Case.Sur.obukhov_length)
        get_GMV_CoVar(
            self,
            self.UpdVar.Area,
            self.UpdVar.W,
            self.UpdVar.W,
            self.EnvVar.W,
            self.EnvVar.W,
            self.EnvVar.TKE,
            GMV.W.values,
            GMV.W.values,
            GMV.TKE.values,
        )
    end

    if self.calc_scalar_var
        self.EnvVar.Hvar.values[grid.gw] =
            get_surface_variance(flux1 * alpha0LL, flux1 * alpha0LL, ustar, zLL, oblength)
        self.EnvVar.QTvar.values[grid.gw] =
            get_surface_variance(flux2 * alpha0LL, flux2 * alpha0LL, ustar, zLL, oblength)
        self.EnvVar.HQTcov.values[grid.gw] =
            get_surface_variance(flux1 * alpha0LL, flux2 * alpha0LL, ustar, zLL, oblength)
        get_GMV_CoVar(
            self,
            self.UpdVar.Area,
            self.UpdVar.H,
            self.UpdVar.H,
            self.EnvVar.H,
            self.EnvVar.H,
            self.EnvVar.Hvar,
            GMV.H.values,
            GMV.H.values,
            GMV.Hvar.values,
        )
        get_GMV_CoVar(
            self,
            self.UpdVar.Area,
            self.UpdVar.QT,
            self.UpdVar.QT,
            self.EnvVar.QT,
            self.EnvVar.QT,
            self.EnvVar.QTvar,
            GMV.QT.values,
            GMV.QT.values,
            GMV.QTvar.values,
        )
        get_GMV_CoVar(
            self,
            self.UpdVar.Area,
            self.UpdVar.H,
            self.UpdVar.QT,
            self.EnvVar.H,
            self.EnvVar.QT,
            self.EnvVar.HQTcov,
            GMV.H.values,
            GMV.QT.values,
            GMV.HQTcov.values,
        )
    end
    return
end

# Find values of environmental variables by subtracting updraft values from grid mean values
# whichvals used to check which substep we are on--correspondingly use "GMV.SomeVar.value" (last timestep value)
# or GMV.SomeVar.mf_update (GMV value following massflux substep)
function decompose_environment(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, whichvals)

    # first make sure the "bulkvalues" of the updraft variables are updated
    set_means(self.UpdVar, GMV)
    grid = get_grid(self)

    gw = grid.gw

    if whichvals == "values"
        @inbounds for k in real_face_indicies(grid)
            val1 = 1.0 / (1.0 - self.UpdVar.Area.bulkvalues[k])
            val2 = self.UpdVar.Area.bulkvalues[k] * val1

            self.EnvVar.Area.values[k] = 1.0 - self.UpdVar.Area.bulkvalues[k]
            self.EnvVar.QT.values[k] = fmax(val1 * GMV.QT.values[k] - val2 * self.UpdVar.QT.bulkvalues[k], 0.0) #Yair - this is here to prevent negative QT
            self.EnvVar.H.values[k] = val1 * GMV.H.values[k] - val2 * self.UpdVar.H.bulkvalues[k]
            # Have to account for staggering of W--interpolate area fraction to the "full" grid points
            # Assuming GMV.W = 0!
            au_full = 0.5 * (self.UpdVar.Area.bulkvalues[k + 1] + self.UpdVar.Area.bulkvalues[k])
            self.EnvVar.W.values[k] = -au_full / (1.0 - au_full) * self.UpdVar.W.bulkvalues[k]
        end

        if self.calc_tke
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.W,
                self.UpdVar.W,
                self.EnvVar.W,
                self.EnvVar.W,
                self.EnvVar.TKE,
                GMV.W.values,
                GMV.W.values,
                GMV.TKE.values,
            )
        end
        if self.calc_scalar_var
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.H,
                self.UpdVar.H,
                self.EnvVar.H,
                self.EnvVar.H,
                self.EnvVar.Hvar,
                GMV.H.values,
                GMV.H.values,
                GMV.Hvar.values,
            )
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.QT,
                self.UpdVar.QT,
                self.EnvVar.QT,
                self.EnvVar.QT,
                self.EnvVar.QTvar,
                GMV.QT.values,
                GMV.QT.values,
                GMV.QTvar.values,
            )
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.H,
                self.UpdVar.QT,
                self.EnvVar.H,
                self.EnvVar.QT,
                self.EnvVar.HQTcov,
                GMV.H.values,
                GMV.QT.values,
                GMV.HQTcov.values,
            )
        end


    elseif whichvals == "mf_update"
        # same as above but replace GMV.SomeVar.values with GMV.SomeVar.mf_update

        @inbounds for k in real_face_indicies(grid)
            val1 = 1.0 / (1.0 - self.UpdVar.Area.bulkvalues[k])
            val2 = self.UpdVar.Area.bulkvalues[k] * val1

            self.EnvVar.QT.values[k] = fmax(val1 * GMV.QT.mf_update[k] - val2 * self.UpdVar.QT.bulkvalues[k], 0.0)#Yair - this is here to prevent negative QT
            self.EnvVar.H.values[k] = val1 * GMV.H.mf_update[k] - val2 * self.UpdVar.H.bulkvalues[k]
            # Have to account for staggering of W
            # Assuming GMV.W = 0!
            au_full = 0.5 * (self.UpdVar.Area.bulkvalues[k + 1] + self.UpdVar.Area.bulkvalues[k])
            self.EnvVar.W.values[k] = -au_full / (1.0 - au_full) * self.UpdVar.W.bulkvalues[k]
        end

        if self.calc_tke
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.W,
                self.UpdVar.W,
                self.EnvVar.W,
                self.EnvVar.W,
                self.EnvVar.TKE,
                GMV.W.values,
                GMV.W.values,
                GMV.TKE.values,
            )
        end

        if self.calc_scalar_var
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.H,
                self.UpdVar.H,
                self.EnvVar.H,
                self.EnvVar.H,
                self.EnvVar.Hvar,
                GMV.H.values,
                GMV.H.values,
                GMV.Hvar.values,
            )
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.QT,
                self.UpdVar.QT,
                self.EnvVar.QT,
                self.EnvVar.QT,
                self.EnvVar.QTvar,
                GMV.QT.values,
                GMV.QT.values,
                GMV.QTvar.values,
            )
            get_GMV_CoVar(
                self,
                self.UpdVar.Area,
                self.UpdVar.H,
                self.UpdVar.QT,
                self.EnvVar.H,
                self.EnvVar.QT,
                self.EnvVar.HQTcov,
                GMV.H.values,
                GMV.QT.values,
                GMV.HQTcov.values,
            )
        end
    end

    return
end

# Note: this assumes all variables are defined on half levels not full levels (i.e. phi, psi are not w)
# if covar_e.name is not "tke".
function get_GMV_CoVar(
    self::EDMF_PrognosticTKE,
    au::UpdraftVariable,
    phi_u::UpdraftVariable,
    psi_u::UpdraftVariable,
    phi_e::EnvironmentVariable,
    psi_e::EnvironmentVariable,
    covar_e::EnvironmentVariable_2m,
    gmv_phi,
    gmv_psi,
    gmv_covar,
)

    grid = get_grid(self)
    ae = pyones(grid.nzg) .- au.bulkvalues
    tke_factor = 1.0
    is_tke = covar_e.name == "tke"

    if is_tke
        @inbounds for k in face_indicies(grid)
            tke_factor = 0.5
            # TODO: report bug: k-1 for k = 0 yields
            # -1, indexing phi_e.values[-1] yields the
            # _last_ value in the array. This is certainly
            # not intended
            if k ≠ 0
                phi_diff = interp2pt(phi_e.values[k - 1] - gmv_phi[k - 1], phi_e.values[k] - gmv_phi[k])
                psi_diff = interp2pt(psi_e.values[k - 1] - gmv_psi[k - 1], psi_e.values[k] - gmv_psi[k])
            else # just use 0th order approximation
                phi_diff = phi_e.values[k] - gmv_phi[k]
                psi_diff = psi_e.values[k] - gmv_psi[k]
            end

            gmv_covar[k] = tke_factor * ae[k] * phi_diff * psi_diff + ae[k] * covar_e.values[k]
            @inbounds for i in xrange(self.n_updrafts)
                # TODO: report bug: k-1 for k = 0 yields
                # -1, indexing phi_e.values[-1] yields the
                # _last_ value in the array. This is certainly
                # not intended
                if k ≠ 0
                    phi_diff = interp2pt(phi_u.values[i, k - 1] - gmv_phi[k - 1], phi_u.values[i, k] - gmv_phi[k])
                    psi_diff = interp2pt(psi_u.values[i, k - 1] - gmv_psi[k - 1], psi_u.values[i, k] - gmv_psi[k])
                else # just use 0th order approximation
                    phi_diff = phi_u.values[i, k] - gmv_phi[k]
                    psi_diff = psi_u.values[i, k] - gmv_psi[k]
                end
                gmv_covar[k] += tke_factor * au.values[i, k] * phi_diff * psi_diff
            end
        end
    else

        @inbounds for k in center_indicies(grid)
            tke_factor = 1.0
            phi_diff = phi_e.values[k] - gmv_phi[k]
            psi_diff = psi_e.values[k] - gmv_psi[k]

            gmv_covar[k] = tke_factor * ae[k] * phi_diff * psi_diff + ae[k] * covar_e.values[k]
            @inbounds for i in xrange(self.n_updrafts)
                phi_diff = phi_u.values[i, k] - gmv_phi[k]
                psi_diff = psi_u.values[i, k] - gmv_psi[k]
                gmv_covar[k] += tke_factor * au.values[i, k] * phi_diff * psi_diff
            end
        end
    end
    return
end


function get_env_covar_from_GMV(
    self::EDMF_PrognosticTKE,
    au::UpdraftVariable,
    phi_u::UpdraftVariable,
    psi_u::UpdraftVariable,
    phi_e::EnvironmentVariable,
    psi_e::EnvironmentVariable,
    covar_e::EnvironmentVariable_2m,
    gmv_phi,
    gmv_psi,
    gmv_covar,
)

    grid = get_grid(self)
    ae = pyones(grid.nzg) .- au.bulkvalues
    tke_factor = 1.0
    is_tke = covar_e.name == "tke"
    if is_tke
        tke_factor = 0.5
    end

    if is_tke
        @inbounds for k in face_indicies(grid)
            if ae[k] > 0.0
                phi_diff = interp2pt(phi_e.values[k - 1] - gmv_phi[k - 1], phi_e.values[k] - gmv_phi[k])
                psi_diff = interp2pt(psi_e.values[k - 1] - gmv_psi[k - 1], psi_e.values[k] - gmv_psi[k])

                covar_e.values[k] = gmv_covar[k] - tke_factor * ae[k] * phi_diff * psi_diff
                @inbounds for i in xrange(self.n_updrafts)
                    phi_diff = interp2pt(phi_u.values[i, k - 1] - gmv_phi[k - 1], phi_u.values[i, k] - gmv_phi[k])
                    psi_diff = interp2pt(psi_u.values[i, k - 1] - gmv_psi[k - 1], psi_u.values[i, k] - gmv_psi[k])
                    covar_e.values[k] -= tke_factor * au.values[i, k] * phi_diff * psi_diff
                end
                covar_e.values[k] = covar_e.values[k] / ae[k]
            else
                covar_e.values[k] = 0.0
            end
        end
    else
        @inbounds for k in center_indicies(grid)
            if ae[k] > 0.0
                phi_diff = phi_e.values[k] - gmv_phi[k]
                psi_diff = psi_e.values[k] - gmv_psi[k]

                covar_e.values[k] = gmv_covar[k] - tke_factor * ae[k] * phi_diff * psi_diff
                @inbounds for i in xrange(self.n_updrafts)
                    phi_diff = phi_u.values[i, k] - gmv_phi[k]
                    psi_diff = psi_u.values[i, k] - gmv_psi[k]

                    covar_e.values[k] -= tke_factor * au.values[i, k] * phi_diff * psi_diff
                end
                covar_e.values[k] = covar_e.values[k] / ae[k]
            else
                covar_e.values[k] = 0.0
            end
        end
    end

    return
end

function compute_turbulent_entrainment(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)

    grid = get_grid(self)
    ref_state = reference_state(self)
    tau = get_mixing_tau(self.base.zi, self.wstar)

    @inbounds for i in xrange(self.n_updrafts)
        @inbounds for k in real_center_indicies(grid)
            a = self.UpdVar.Area.values[i, k]
            a_full = interp2pt(self.UpdVar.Area.values[i, k], self.UpdVar.Area.values[i, k + 1])
            R_up = self.pressure_plume_spacing[i]
            R_up_full = self.pressure_plume_spacing[i]
            wu_half = interp2pt(self.UpdVar.W.values[i, k], self.UpdVar.W.values[i, k - 1])
            we_half = interp2pt(self.EnvVar.W.values[k], self.EnvVar.W.values[k - 1])
            if a * wu_half > 0.0
                self.turb_entr_H[i, k] =
                    (2.0 / R_up^2.0) *
                    ref_state.rho0_half[k] *
                    a *
                    self.horizontal_KH[i, k] *
                    (self.EnvVar.H.values[k] - self.UpdVar.H.values[i, k])
                self.turb_entr_QT[i, k] =
                    (2.0 / R_up^2.0) *
                    ref_state.rho0_half[k] *
                    a *
                    self.horizontal_KH[i, k] *
                    (self.EnvVar.QT.values[k] - self.UpdVar.QT.values[i, k])
                self.frac_turb_entr[i, k] = (2.0 / R_up^2.0) * self.horizontal_KH[i, k] / wu_half / a
            else
                self.turb_entr_H[i, k] = 0.0
                self.turb_entr_QT[i, k] = 0.0
            end

            if a_full * self.UpdVar.W.values[i, k] > 0.0
                K_full = interp2pt(self.horizontal_KM[i, k], self.horizontal_KM[i, k - 1])
                b_upd_full = interp2pt(self.UpdVar.B.values[i, k], self.UpdVar.B.values[i, k - 1])
                b_env_full = interp2pt(self.EnvVar.B.values[k], self.EnvVar.B.values[k - 1])
                env_tke_full = interp2pt(self.EnvVar.TKE.buoy[k], self.EnvVar.TKE.buoy[k - 1])

                self.turb_entr_W[i, k] =
                    (2.0 / R_up_full^2.0) *
                    ref_state.rho0[k] *
                    a_full *
                    K_full *
                    (self.EnvVar.W.values[k] - self.UpdVar.W.values[i, k])
                self.frac_turb_entr_full[i, k] = (2.0 / R_up_full^2.0) * K_full / self.UpdVar.W.values[i, k] / a_full
            else
                self.turb_entr_W[i, k] = 0.0
            end
        end
    end
    return
end

function compute_entrainment_detrainment(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    ret = entr_struct()
    input = entr_in_struct()
    sa = eos_struct()
    quadrature_order = 3
    grid = get_grid(self)

    upd_cloud_diagnostics(self.UpdVar, reference_state(self))

    input.wstar = self.wstar

    input.dz = get_grid(self).dz
    input.zbl = compute_zbl_qt_grad(self, GMV)
    input.sort_pow = self.sorting_power
    input.c_ent = self.entrainment_factor
    input.c_det = self.detrainment_factor
    input.c_div = self.entrainment_Mdiv_factor
    input.c_mu = self.entrainment_sigma
    input.c_mu0 = self.entrainment_scale
    input.c_ed_mf = self.entrainment_ed_mf_sigma
    input.chi_upd = self.updraft_mixing_frac
    input.tke_coef = self.entrainment_smin_tke_coeff
    @inbounds for i in xrange(self.n_updrafts)
        input.zi = self.UpdVar.cloud_base[i]
        @inbounds for k in real_center_indicies(grid)
            input.buoy_ed_flux = self.EnvVar.TKE.buoy[k]
            if self.UpdVar.Area.values[i, k] > 0.0
                input.b_upd = self.UpdVar.B.values[i, k]
                input.w_upd = interp2pt(self.UpdVar.W.values[i, k], self.UpdVar.W.values[i, k - 1])
                input.z = grid.z_half[k]
                input.a_upd = self.UpdVar.Area.values[i, k]
                input.a_env = 1.0 - self.UpdVar.Area.bulkvalues[k]
                input.tke = self.EnvVar.TKE.values[k]
                input.ql_env = self.EnvVar.QL.values[k]
                input.b_env = self.EnvVar.B.values[k]
                input.w_env = interp2pt(self.EnvVar.W.values[k], self.EnvVar.W.values[k - 1])
                input.ql_up = self.UpdVar.QL.values[i, k]
                input.ql_env = self.EnvVar.QL.values[k]
                input.nh_pressure = self.nh_pressure[i, k]
                input.RH_upd = self.UpdVar.RH.values[i, k]
                input.RH_env = self.EnvVar.RH.values[k]

                # compute dMdz at half levels
                gmv_w_k = interp2pt(GMV.W.values[k], GMV.W.values[k - 1])
                gmv_w_km = interp2pt(GMV.W.values[k - 1], GMV.W.values[k - 2])
                upd_w_km = interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k - 2])
                input.M = input.a_upd * (input.w_upd - gmv_w_k)
                Mm = self.UpdVar.Area.values[i, k - 1] * (upd_w_km - gmv_w_km)
                input.dMdz = (input.M - Mm) * grid.dzi

                if self.calc_tke
                    input.tke = self.EnvVar.TKE.values[k]
                end

                ret = self.entr_detr_fp(input)
                self.entr_sc[i, k] = ret.entr_sc
                self.detr_sc[i, k] = ret.detr_sc
                self.sorting_function[i, k] = ret.sorting_function
                self.b_mix[i, k] = ret.b_mix

            else
                self.entr_sc[i, k] = 0.0
                self.detr_sc[i, k] = 0.0
                self.sorting_function[i, k] = 0.0
                self.b_mix[i, k] = self.EnvVar.B.values[k]
            end
        end
    end
    return
end

function compute_zbl_qt_grad(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)
    # computes inversion height as z with max gradient of qt
    zbl_qt = 0.0
    qt_grad = 0.0
    grid = get_grid(self)
    dzi = grid.dzi

    @inbounds for k in real_center_indicies(grid)
        z_ = grid.z_half[k]
        qt_up = GMV.QT.values[k + 1]
        qt_ = GMV.QT.values[k]

        if fabs(qt_up - qt_) * dzi > qt_grad
            qt_grad = fabs(qt_up - qt_) * dzi
            zbl_qt = z_
        end
    end
    return zbl_qt
end

function compute_pressure_plume_spacing(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)

    @inbounds for i in xrange(self.n_updrafts)
        if self.use_const_plume_spacing
            self.pressure_plume_spacing[i] = self.constant_plume_spacing
        else
            self.pressure_plume_spacing[i] =
                fmax(self.aspect_ratio * self.UpdVar.updraft_top[i], 500.0 * self.aspect_ratio)
        end
    end
    return
end

function compute_nh_pressure(self::EDMF_PrognosticTKE)
    ret_b = pressure_buoy_struct()
    ret_w = pressure_drag_struct()
    input = pressure_in_struct()
    grid = get_grid(self)

    # TODO: this should eventually match the grid
    cinterior = (get_grid(self).gw):(get_grid(self).nzg - get_grid(self).gw)

    @inbounds for i in xrange(self.n_updrafts)
        input.updraft_top = self.UpdVar.updraft_top[i]
        alen = max(length(argwhere(self.UpdVar.Area.values[i, cinterior])) - 1, 0)
        avals = off_arr(self.UpdVar.Area.values[i, cinterior])
        input.a_med = Statistics.median(avals[0:alen])
        @inbounds for k in real_center_indicies(grid)
            input.a_kfull = interp2pt(self.UpdVar.Area.values[i, k], self.UpdVar.Area.values[i, k + 1])
            if input.a_kfull >= self.minimum_area
                input.dzi = get_grid(self).dzi

                input.b_kfull = interp2pt(self.UpdVar.B.values[i, k], self.UpdVar.B.values[i, k + 1])
                input.rho0_kfull = reference_state(self).rho0[k]
                input.alpha1 = self.pressure_normalmode_buoy_coeff1
                input.alpha2 = self.pressure_normalmode_buoy_coeff2
                input.beta1 = self.pressure_normalmode_adv_coeff
                input.beta2 = self.pressure_normalmode_drag_coeff
                input.w_kfull = self.UpdVar.W.values[i, k]
                input.w_kmfull = self.UpdVar.W.values[i, k - 1]
                input.w_kenv = self.EnvVar.W.values[k]

                if self.asp_label == "z_dependent"
                    input.asp_ratio = input.updraft_top / 2.0 / sqrt(input.a_kfull) / input.rd
                elseif self.asp_label == "median"
                    input.asp_ratio = input.updraft_top / 2.0 / sqrt(input.a_med) / input.rd
                elseif self.asp_label == "const"
                    input.asp_ratio = 1.0
                end

                if input.a_kfull > 0.0
                    ret_b = self.pressure_func_buoy(input)
                    ret_w = self.pressure_func_drag(input)
                    self.nh_pressure_b[i, k] = ret_b.nh_pressure_b
                    self.nh_pressure_adv[i, k] = ret_w.nh_pressure_adv
                    self.nh_pressure_drag[i, k] = ret_w.nh_pressure_drag

                    self.b_coeff[i, k] = ret_b.b_coeff
                    self.asp_ratio[i, k] = input.asp_ratio

                else
                    self.nh_pressure_b[i, k] = 0.0
                    self.nh_pressure_adv[i, k] = 0.0
                    self.nh_pressure_drag[i, k] = 0.0

                    self.b_coeff[i, k] = 0.0
                    self.asp_ratio[i, k] = 0.0
                end
            else
                self.nh_pressure_b[i, k] = 0.0
                self.nh_pressure_adv[i, k] = 0.0
                self.nh_pressure_drag[i, k] = 0.0

                self.b_coeff[i, k] = 0.0
                self.asp_ratio[i, k] = 0.0
            end

            self.nh_pressure[i, k] = self.nh_pressure_b[i, k] + self.nh_pressure_adv[i, k] + self.nh_pressure_drag[i, k]
        end
    end

    return
end

function zero_area_fraction_cleanup(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)

    grid = get_grid(self)
    @inbounds for k in real_center_indicies(grid)
        @inbounds for i in xrange(self.n_updrafts)
            if self.UpdVar.Area.values[i, k] < self.minimum_area
                self.UpdVar.Area.values[i, k] = 0.0
                self.UpdVar.W.values[i, k] = GMV.W.values[k]
                self.UpdVar.B.values[i, k] = GMV.B.values[k]
                self.UpdVar.H.values[i, k] = GMV.H.values[k]
                self.UpdVar.QT.values[i, k] = GMV.QT.values[k]
                self.UpdVar.T.values[i, k] = GMV.T.values[k]
                self.UpdVar.QL.values[i, k] = GMV.QL.values[k]
                self.UpdVar.THL.values[i, k] = GMV.THL.values[k]
            end
        end

        if sum(self.UpdVar.Area.values[:, k]) == 0.0
            self.EnvVar.W.values[k] = GMV.W.values[k]
            self.EnvVar.B.values[k] = GMV.B.values[k]
            self.EnvVar.H.values[k] = GMV.H.values[k]
            self.EnvVar.QT.values[k] = GMV.QT.values[k]
            self.EnvVar.T.values[k] = GMV.T.values[k]
            self.EnvVar.QL.values[k] = GMV.QL.values[k]
            self.EnvVar.THL.values[k] = GMV.THL.values[k]
        end
    end

    return
end


function set_subdomain_bcs(self::EDMF_PrognosticTKE)

    grid = get_grid(self)
    set_bcs(self.UpdVar.W, grid)
    set_bcs(self.UpdVar.Area, grid)
    set_bcs(self.UpdVar.H, grid)
    set_bcs(self.UpdVar.THL, grid)
    set_bcs(self.UpdVar.QT, grid)
    set_bcs(self.UpdVar.T, grid)
    set_bcs(self.UpdVar.B, grid)

    set_bcs(self.EnvVar.W, grid)
    set_bcs(self.EnvVar.H, grid)
    set_bcs(self.EnvVar.THL, grid)
    set_bcs(self.EnvVar.T, grid)
    set_bcs(self.EnvVar.QL, grid)
    set_bcs(self.EnvVar.QT, grid)

    return
end

function solve_updraft_velocity_area(self::EDMF_PrognosticTKE)
    grid = get_grid(self)
    ref_state = reference_state(self)
    gw = grid.gw
    dzi = grid.dzi
    dti_ = 1.0 / self.dt_upd
    dt_ = 1.0 / dti_

    @inbounds for i in xrange(self.n_updrafts)
        self.entr_sc[i, gw] = self.entr_surface_bc
        self.detr_sc[i, gw] = self.detr_surface_bc
        self.UpdVar.W.new[i, gw - 1] = self.w_surface_bc[i]
        self.UpdVar.Area.new[i, gw] = self.area_surface_bc[i]
        au_lim = self.max_area

        @inbounds for k in real_center_indicies(grid)

            # First solve for updated area fraction at k+1
            whalf_kp = interp2pt(self.UpdVar.W.values[i, k], self.UpdVar.W.values[i, k + 1])
            whalf_k = interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k])
            adv =
                -ref_state.alpha0_half[k + 1] *
                dzi *
                (
                    ref_state.rho0_half[k + 1] * self.UpdVar.Area.values[i, k + 1] * whalf_kp -
                    ref_state.rho0_half[k] * self.UpdVar.Area.values[i, k] * whalf_k
                )
            entr_term = self.UpdVar.Area.values[i, k + 1] * whalf_kp * (self.entr_sc[i, k + 1])
            detr_term = self.UpdVar.Area.values[i, k + 1] * whalf_kp * (-self.detr_sc[i, k + 1])

            self.UpdVar.Area.new[i, k + 1] =
                fmax(dt_ * (adv + entr_term + detr_term) + self.UpdVar.Area.values[i, k + 1], 0.0)

            if self.UpdVar.Area.new[i, k + 1] > au_lim
                self.UpdVar.Area.new[i, k + 1] = au_lim
                if self.UpdVar.Area.values[i, k + 1] > 0.0
                    self.detr_sc[i, k + 1] = (
                        ((au_lim - self.UpdVar.Area.values[i, k + 1]) * dti_ - adv - entr_term) /
                        (-self.UpdVar.Area.values[i, k + 1] * whalf_kp)
                    )
                else
                    # this detrainment rate won"t affect scalars but would affect velocity
                    self.detr_sc[i, k + 1] =
                        (((au_lim - self.UpdVar.Area.values[i, k + 1]) * dti_ - adv - entr_term) / (-au_lim * whalf_kp))
                end
            end

            # Now solve for updraft velocity at k
            rho_ratio = ref_state.rho0[k - 1] / ref_state.rho0[k]
            anew_k = interp2pt(self.UpdVar.Area.new[i, k], self.UpdVar.Area.new[i, k + 1])

            if anew_k >= self.minimum_area
                a_k = interp2pt(self.UpdVar.Area.values[i, k], self.UpdVar.Area.values[i, k + 1])
                a_km = interp2pt(self.UpdVar.Area.values[i, k - 1], self.UpdVar.Area.values[i, k])
                entr_w = interp2pt(self.entr_sc[i, k], self.entr_sc[i, k + 1])
                detr_w = interp2pt(self.detr_sc[i, k], self.detr_sc[i, k + 1])
                B_k = interp2pt(self.UpdVar.B.values[i, k], self.UpdVar.B.values[i, k + 1])

                adv = (
                    ref_state.rho0[k] * a_k * self.UpdVar.W.values[i, k] * self.UpdVar.W.values[i, k] * dzi -
                    ref_state.rho0[k - 1] *
                    a_km *
                    self.UpdVar.W.values[i, k - 1] *
                    self.UpdVar.W.values[i, k - 1] *
                    dzi
                )
                exch = (
                    ref_state.rho0[k] *
                    a_k *
                    self.UpdVar.W.values[i, k] *
                    (entr_w * self.EnvVar.W.values[k] - detr_w * self.UpdVar.W.values[i, k]) +
                    self.turb_entr_W[i, k]
                )
                buoy = ref_state.rho0[k] * a_k * B_k
                self.UpdVar.W.new[i, k] =
                    (
                        ref_state.rho0[k] * a_k * self.UpdVar.W.values[i, k] * dti_ - adv +
                        exch +
                        buoy +
                        self.nh_pressure[i, k]
                    ) / (ref_state.rho0[k] * anew_k * dti_)

                if self.UpdVar.W.new[i, k] <= 0.0
                    self.UpdVar.W.new[i, k] = 0.0
                    self.UpdVar.Area.new[i, k + 1] = 0.0
                    #break
                end
            else
                self.UpdVar.W.new[i, k] = 0.0
                self.UpdVar.Area.new[i, k + 1] = 0.0
                # keep this in mind if we modify updraft top treatment!
            end
        end
    end
    return
end

function solve_updraft_scalars(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)
    grid = get_grid(self)
    ref_state = reference_state(self)
    dzi = grid.dzi
    dti_ = 1.0 / self.dt_upd
    gw = grid.gw
    sa = eos_struct()

    @inbounds for i in xrange(self.n_updrafts)

        # starting from the bottom do entrainment at each level
        @inbounds for k in real_center_indicies(grid)
            if k == gw
                # at the surface
                if self.UpdVar.Area.new[i, k] >= self.minimum_area
                    self.UpdVar.H.new[i, k] = self.h_surface_bc[i]
                    self.UpdVar.QT.new[i, k] = self.qt_surface_bc[i]
                else
                    self.UpdVar.H.new[i, k] = GMV.H.values[k]
                    self.UpdVar.QT.new[i, k] = GMV.QT.values[k]
                end

                # saturation adjustment
                sa = eos(
                    self.UpdThermo.t_to_prog_fp,
                    self.UpdThermo.prog_to_t_fp,
                    ref_state.p0_half[k],
                    self.UpdVar.QT.new[i, k],
                    self.UpdVar.H.new[i, k],
                )
                self.UpdVar.QL.new[i, k] = sa.ql
                self.UpdVar.T.new[i, k] = sa.T
                continue
            end


            H_entr = self.EnvVar.H.values[k]
            QT_entr = self.EnvVar.QT.values[k]

            # write the discrete equations in form
            # c1 * phi_new[k] = c2 * phi[k] + c3 * phi[k-1] + c4 * phi_entr
            if self.UpdVar.Area.new[i, k] >= self.minimum_area
                m_k = (
                    ref_state.rho0_half[k] *
                    self.UpdVar.Area.values[i, k] *
                    interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k])
                )
                m_km = (
                    ref_state.rho0_half[k - 1] *
                    self.UpdVar.Area.values[i, k - 1] *
                    interp2pt(self.UpdVar.W.values[i, k - 2], self.UpdVar.W.values[i, k - 1])
                )
                c1 = ref_state.rho0_half[k] * self.UpdVar.Area.new[i, k] * dti_
                c2 = (ref_state.rho0_half[k] * self.UpdVar.Area.values[i, k] * dti_ - m_k * (dzi + self.detr_sc[i, k]))
                c3 = m_km * dzi
                c4 = m_k * self.entr_sc[i, k]

                self.UpdVar.H.new[i, k] =
                    (
                        c2 * self.UpdVar.H.values[i, k] +
                        c3 * self.UpdVar.H.values[i, k - 1] +
                        c4 * H_entr +
                        self.turb_entr_H[i, k]
                    ) / c1
                self.UpdVar.QT.new[i, k] =
                    (
                        c2 * self.UpdVar.QT.values[i, k] +
                        c3 * self.UpdVar.QT.values[i, k - 1] +
                        c4 * QT_entr +
                        self.turb_entr_QT[i, k]
                    ) / c1

            else
                self.UpdVar.H.new[i, k] = GMV.H.values[k]
                self.UpdVar.QT.new[i, k] = GMV.QT.values[k]
            end

            # saturation adjustment
            sa = eos(
                self.UpdThermo.t_to_prog_fp,
                self.UpdThermo.prog_to_t_fp,
                ref_state.p0_half[k],
                self.UpdVar.QT.new[i, k],
                self.UpdVar.H.new[i, k],
            )
            self.UpdVar.QL.new[i, k] = sa.ql
            self.UpdVar.T.new[i, k] = sa.T
        end
    end

    return
end

# After updating the updraft variables themselves
# 1. compute the mass fluxes (currently not stored as class members, probably will want to do this
# for output purposes)
# 2. Apply mass flux tendencies and updraft microphysical tendencies to GMV.SomeVar.Values (old time step values)
# thereby updating to GMV.SomeVar.mf_update
# mass flux tendency is computed as 1st order upwind

function update_GMV_MF(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, TS::TimeStepping)
    grid = get_grid(self)
    ref_state = reference_state(self)
    rho0 = ref_state.rho0
    gw = grid.gw
    mf_tend_h = 0.0
    mf_tend_qt = 0.0
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues # area of environment

    self.massflux_h .= 0.0
    self.massflux_qt .= 0.0

    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in xrange(self.n_updrafts)
        self.m[i, gw - 1] = 0.0
        @inbounds for k in real_face_indicies(grid)
            a = interp2pt(self.UpdVar.Area.values[i, k], self.UpdVar.Area.values[i, k + 1])
            self.m[i, k] =
                rho0[k] * a * interp2pt(ae[k], ae[k + 1]) * (self.UpdVar.W.values[i, k] - self.EnvVar.W.values[k])
        end
    end


    self.massflux_h[gw - 1] = 0.0
    self.massflux_qt[gw - 1] = 0.0

    @inbounds for k in real_face_indicies(grid)
        self.massflux_h[k] = 0.0
        self.massflux_qt[k] = 0.0
        env_h_interp = interp2pt(self.EnvVar.H.values[k], self.EnvVar.H.values[k + 1])
        env_qt_interp = interp2pt(self.EnvVar.QT.values[k], self.EnvVar.QT.values[k + 1])
        @inbounds for i in xrange(self.n_updrafts)
            self.massflux_h[k] +=
                self.m[i, k] * (interp2pt(self.UpdVar.H.values[i, k], self.UpdVar.H.values[i, k + 1]) - env_h_interp)
            self.massflux_qt[k] +=
                self.m[i, k] * (interp2pt(self.UpdVar.QT.values[i, k], self.UpdVar.QT.values[i, k + 1]) - env_qt_interp)
        end
    end

    # Compute the  mass flux tendencies
    # Adjust the values of the grid mean variables
    @inbounds for k in real_center_indicies(grid)
        mf_tend_h = -(self.massflux_h[k] - self.massflux_h[k - 1]) * (ref_state.alpha0_half[k] * grid.dzi)
        mf_tend_qt = -(self.massflux_qt[k] - self.massflux_qt[k - 1]) * (ref_state.alpha0_half[k] * grid.dzi)

        GMV.H.mf_update[k] = GMV.H.values[k] + TS.dt * mf_tend_h + self.UpdThermo.prec_source_h_tot[k]
        GMV.QT.mf_update[k] = GMV.QT.values[k] + TS.dt * mf_tend_qt + self.UpdThermo.prec_source_qt_tot[k]

        #No mass flux tendency for U, V
        GMV.U.mf_update[k] = GMV.U.values[k]
        GMV.V.mf_update[k] = GMV.V.values[k]
        # Prepare the output
        self.massflux_tendency_h[k] = mf_tend_h
        self.massflux_tendency_qt[k] = mf_tend_qt
    end
    return
end

# Update the grid mean variables with the tendency due to eddy diffusion
# Km and Kh have already been updated
# 2nd order finite differences plus implicit time step allows solution with tridiagonal matrix solver
# Update from GMV.SomeVar.mf_update to GMV.SomeVar.new
function update_GMV_ED(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)
    grid = get_grid(self)
    gw = grid.gw
    nzg = grid.nzg
    nz = grid.nz
    dzi = grid.dzi
    ref_state = reference_state(self)
    a = pyzeros(nz) # for tridiag solver
    b = pyzeros(nz) # for tridiag solver
    c = pyzeros(nz) # for tridiag solver
    x = pyzeros(nz) # for tridiag solver
    ae = pyones(nzg) .- self.UpdVar.Area.bulkvalues # area of environment
    rho_ae_K = pyzeros(nzg)
    KM = diffusivity_m(self)
    KH = diffusivity_h(self)
    @inbounds for k in xrange(nzg - 1)
        rho_ae_K[k] = 0.5 * (ae[k] * KH.values[k] + ae[k + 1] * KH.values[k + 1]) * ref_state.rho0[k]
    end

    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    construct_tridiag_diffusion(grid, TS.dt, rho_ae_K, ref_state.rho0_half, ae, a, b, c)
    # Solve QT
    @inbounds for k in xrange(nz)
        x[k] = self.EnvVar.QT.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_qtflux * dzi * ref_state.alpha0_half[gw] / ae[gw]
    x .= tridiag_solve(x, a, b, c)

    @inbounds for k in xrange(nz)
        GMV.QT.new[k + gw] = fmax(
            GMV.QT.mf_update[k + gw] +
            ae[k + gw] * (x[k] - self.EnvVar.QT.values[k + gw]) +
            self.EnvThermo.prec_source_qt[k + gw] +
            self.rainphysics.rain_evap_source_qt[k + gw],
            0.0,
        )
        self.diffusive_tendency_qt[k + gw] = (GMV.QT.new[k + gw] - GMV.QT.mf_update[k + gw]) * TS.dti
    end
    # get the diffusive flux
    self.diffusive_flux_qt[gw] = interp2pt(
        Case.Sur.rho_qtflux,
        -rho_ae_K[gw] * dzi * (self.EnvVar.QT.values[gw + 1] - self.EnvVar.QT.values[gw]),
    )
    @inbounds for k in real_center_indicies(grid)
        k == gw && continue # skip bc
        self.diffusive_flux_qt[k] =
            -0.5 *
            ref_state.rho0_half[k] *
            ae[k] *
            KH.values[k] *
            dzi *
            (self.EnvVar.QT.values[k + 1] - self.EnvVar.QT.values[k - 1])
    end

    # Solve H
    @inbounds for k in xrange(nz)
        x[k] = self.EnvVar.H.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_hflux * dzi * ref_state.alpha0_half[gw] / ae[gw]
    x .= tridiag_solve(x, a, b, c)
    @inbounds for k in xrange(nz)
        GMV.H.new[k + gw] =
            GMV.H.mf_update[k + gw] +
            ae[k + gw] * (x[k] - self.EnvVar.H.values[k + gw]) +
            self.EnvThermo.prec_source_h[k + gw] +
            self.rainphysics.rain_evap_source_h[k + gw]
        self.diffusive_tendency_h[k + gw] = (GMV.H.new[k + gw] - GMV.H.mf_update[k + gw]) * TS.dti
    end
    # get the diffusive flux
    self.diffusive_flux_h[gw] =
        interp2pt(Case.Sur.rho_hflux, -rho_ae_K[gw] * dzi * (self.EnvVar.H.values[gw + 1] - self.EnvVar.H.values[gw]))
    @inbounds for k in real_center_indicies(grid)
        k == gw && continue # skip bc
        self.diffusive_flux_h[k] =
            -0.5 *
            ref_state.rho0_half[k] *
            ae[k] *
            KH.values[k] *
            dzi *
            (self.EnvVar.H.values[k + 1] - self.EnvVar.H.values[k - 1])
    end

    # Solve U
    @inbounds for k in xrange(nzg - 1)
        rho_ae_K[k] = 0.5 * (ae[k] * KM.values[k] + ae[k + 1] * KM.values[k + 1]) * ref_state.rho0[k]
    end

    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    construct_tridiag_diffusion(grid, TS.dt, rho_ae_K, ref_state.rho0_half, ae, a, b, c)
    @inbounds for k in xrange(nz)
        x[k] = GMV.U.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_uflux * dzi * ref_state.alpha0_half[gw] / ae[gw]
    x .= tridiag_solve(x, a, b, c)

    @inbounds for k in xrange(nz)
        GMV.U.new[k + gw] = x[k]
    end
    self.diffusive_flux_u[gw] =
        interp2pt(Case.Sur.rho_uflux, -rho_ae_K[gw] * dzi * (GMV.U.values[gw + 1] - GMV.U.values[gw]))
    @inbounds for k in real_center_indicies(grid)
        k == gw && continue # skip bc
        self.diffusive_flux_u[k] =
            -0.5 * ref_state.rho0_half[k] * ae[k] * KM.values[k] * dzi * (GMV.U.values[k + 1] - GMV.U.values[k - 1])
    end

    # Solve V
    @inbounds for k in xrange(nz)
        x[k] = GMV.V.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_vflux * dzi * ref_state.alpha0_half[gw] / ae[gw]
    x .= tridiag_solve(x, a, b, c)
    @inbounds for k in xrange(nz)
        GMV.V.new[k + gw] = x[k]
    end
    self.diffusive_flux_v[gw] =
        interp2pt(Case.Sur.rho_vflux, -rho_ae_K[gw] * dzi * (GMV.V.values[gw + 1] - GMV.V.values[gw]))
    @inbounds for k in real_center_indicies(grid)
        k == gw && continue # skip bc
        self.diffusive_flux_v[k] =
            -0.5 * ref_state.rho0_half[k] * ae[k] * KM.values[k] * dzi * (GMV.V.values[k + 1] - GMV.V.values[k - 1])
    end
    set_bcs(GMV.QT, grid)
    set_bcs(GMV.THL, grid)
    set_bcs(GMV.H, grid)
    set_bcs(GMV.U, grid)
    set_bcs(GMV.V, grid)

    return
end

function compute_tke_buoy(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)
    grid = get_grid(self)
    ref_state = reference_state(self)
    grad_thl_minus = 0.0
    grad_qt_minus = 0.0
    grad_thl_plus = 0
    grad_qt_plus = 0
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues
    KH = diffusivity_h(self).values

    # Note that source terms at the gw grid point are not really used because that is where tke boundary condition is
    # enforced (according to MO similarity). Thus here I am being sloppy about lowest grid point
    @inbounds for k in real_center_indicies(grid)
        qt_dry = self.EnvThermo.qt_dry[k]
        th_dry = self.EnvThermo.th_dry[k]
        t_cloudy = self.EnvThermo.t_cloudy[k]
        qv_cloudy = self.EnvThermo.qv_cloudy[k]
        qt_cloudy = self.EnvThermo.qt_cloudy[k]
        th_cloudy = self.EnvThermo.th_cloudy[k]
        lh = latent_heat(t_cloudy)
        cpm = cpm_c(qt_cloudy)
        grad_thl_minus = grad_thl_plus
        grad_qt_minus = grad_qt_plus
        grad_thl_plus = (self.EnvVar.THL.values[k + 1] - self.EnvVar.THL.values[k]) * grid.dzi
        grad_qt_plus = (self.EnvVar.QT.values[k + 1] - self.EnvVar.QT.values[k]) * grid.dzi
        prefactor = Rd * exner_c(ref_state.p0_half[k]) / ref_state.p0_half[k]
        d_alpha_thetal_dry = prefactor * (1.0 + (eps_vi - 1.0) * qt_dry)
        d_alpha_qt_dry = prefactor * th_dry * (eps_vi - 1.0)

        if self.EnvVar.cloud_fraction.values[k] > 0.0
            d_alpha_thetal_cloudy = (
                prefactor * (1.0 + eps_vi * (1.0 + lh / Rv / t_cloudy) * qv_cloudy - qt_cloudy) /
                (1.0 + lh * lh / cpm / Rv / t_cloudy / t_cloudy * qv_cloudy)
            )
            d_alpha_qt_cloudy = (lh / cpm / t_cloudy * d_alpha_thetal_cloudy - prefactor) * th_cloudy
        else
            d_alpha_thetal_cloudy = 0.0
            d_alpha_qt_cloudy = 0.0
        end

        d_alpha_thetal_total = (
            self.EnvVar.cloud_fraction.values[k] * d_alpha_thetal_cloudy +
            (1.0 - self.EnvVar.cloud_fraction.values[k]) * d_alpha_thetal_dry
        )
        d_alpha_qt_total = (
            self.EnvVar.cloud_fraction.values[k] * d_alpha_qt_cloudy +
            (1.0 - self.EnvVar.cloud_fraction.values[k]) * d_alpha_qt_dry
        )

        # TODO - check
        self.EnvVar.TKE.buoy[k] =
            g / ref_state.alpha0_half[k] *
            ae[k] *
            ref_state.rho0_half[k] *
            (
                -KH[k] * interp2pt(grad_thl_plus, grad_thl_minus) * d_alpha_thetal_total -
                KH[k] * interp2pt(grad_qt_plus, grad_qt_minus) * d_alpha_qt_total
            )
    end
    return
end

function compute_tke_pressure(self::EDMF_PrognosticTKE)
    grid = get_grid(self)

    @inbounds for k in real_center_indicies(grid)
        self.EnvVar.TKE.press[k] = 0.0
        @inbounds for i in xrange(self.n_updrafts)
            wu_half = interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k])
            we_half = interp2pt(self.EnvVar.W.values[k - 1], self.EnvVar.W.values[k])
            press_half = interp2pt(self.nh_pressure[i, k - 1], self.nh_pressure[i, k])
            self.EnvVar.TKE.press[k] += (we_half - wu_half) * press_half
        end
    end
    return
end

function update_GMV_diagnostics(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)

    grid = get_grid(self)
    p0_half = reference_state(self).p0_half
    @inbounds for k in real_center_indicies(grid)
        GMV.QL.values[k] = (
            self.UpdVar.Area.bulkvalues[k] * self.UpdVar.QL.bulkvalues[k] +
            (1.0 - self.UpdVar.Area.bulkvalues[k]) * self.EnvVar.QL.values[k]
        )

        GMV.T.values[k] = (
            self.UpdVar.Area.bulkvalues[k] * self.UpdVar.T.bulkvalues[k] +
            (1.0 - self.UpdVar.Area.bulkvalues[k]) * self.EnvVar.T.values[k]
        )

        qv = GMV.QT.values[k] - GMV.QL.values[k]

        GMV.THL.values[k] = t_to_thetali_c(p0_half[k], GMV.T.values[k], GMV.QT.values[k], GMV.QL.values[k], 0.0)

        GMV.B.values[k] = (
            self.UpdVar.Area.bulkvalues[k] * self.UpdVar.B.bulkvalues[k] +
            (1.0 - self.UpdVar.Area.bulkvalues[k]) * self.EnvVar.B.values[k]
        )
    end

    return
end

function compute_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    if self.similarity_diffusivity # otherwise, we computed mixing length when we computed
        compute_mixing_length(self, Case.Sur.obukhov_length, Case.Sur.ustar, GMV)
    end
    if self.calc_tke
        compute_tke_buoy(self, GMV)
        compute_covariance_entr(
            self,
            self.EnvVar.TKE,
            self.UpdVar.W,
            self.UpdVar.W,
            self.EnvVar.W,
            self.EnvVar.W,
            GMV.W,
            GMV.W,
        )
        compute_covariance_shear(
            self,
            GMV,
            self.EnvVar.TKE,
            self.UpdVar.W.values,
            self.UpdVar.W.values,
            self.EnvVar.W.values,
            self.EnvVar.W.values,
        )
        compute_covariance_interdomain_src(
            self,
            self.UpdVar.Area,
            self.UpdVar.W,
            self.UpdVar.W,
            self.EnvVar.W,
            self.EnvVar.W,
            self.EnvVar.TKE,
        )
        compute_tke_pressure(self)
    end
    if self.calc_scalar_var
        compute_covariance_entr(
            self,
            self.EnvVar.Hvar,
            self.UpdVar.H,
            self.UpdVar.H,
            self.EnvVar.H,
            self.EnvVar.H,
            GMV.H,
            GMV.H,
        )
        compute_covariance_entr(
            self,
            self.EnvVar.QTvar,
            self.UpdVar.QT,
            self.UpdVar.QT,
            self.EnvVar.QT,
            self.EnvVar.QT,
            GMV.QT,
            GMV.QT,
        )
        compute_covariance_entr(
            self,
            self.EnvVar.HQTcov,
            self.UpdVar.H,
            self.UpdVar.QT,
            self.EnvVar.H,
            self.EnvVar.QT,
            GMV.H,
            GMV.QT,
        )
        compute_covariance_shear(
            self,
            GMV,
            self.EnvVar.Hvar,
            self.UpdVar.H.values,
            self.UpdVar.H.values,
            self.EnvVar.H.values,
            self.EnvVar.H.values,
        )
        compute_covariance_shear(
            self,
            GMV,
            self.EnvVar.QTvar,
            self.UpdVar.QT.values,
            self.UpdVar.QT.values,
            self.EnvVar.QT.values,
            self.EnvVar.QT.values,
        )
        compute_covariance_shear(
            self,
            GMV,
            self.EnvVar.HQTcov,
            self.UpdVar.H.values,
            self.UpdVar.QT.values,
            self.EnvVar.H.values,
            self.EnvVar.QT.values,
        )
        compute_covariance_interdomain_src(
            self,
            self.UpdVar.Area,
            self.UpdVar.H,
            self.UpdVar.H,
            self.EnvVar.H,
            self.EnvVar.H,
            self.EnvVar.Hvar,
        )
        compute_covariance_interdomain_src(
            self,
            self.UpdVar.Area,
            self.UpdVar.QT,
            self.UpdVar.QT,
            self.EnvVar.QT,
            self.EnvVar.QT,
            self.EnvVar.QTvar,
        )
        compute_covariance_interdomain_src(
            self,
            self.UpdVar.Area,
            self.UpdVar.H,
            self.UpdVar.QT,
            self.EnvVar.H,
            self.EnvVar.QT,
            self.EnvVar.HQTcov,
        )
        compute_covariance_rain(self, TS, GMV) # need to update this one

        GMV_third_m(self, GMV.H_third_m, self.EnvVar.Hvar, self.EnvVar.H, self.UpdVar.H)
        GMV_third_m(self, GMV.QT_third_m, self.EnvVar.QTvar, self.EnvVar.QT, self.UpdVar.QT)
        GMV_third_m(self, GMV.W_third_m, self.EnvVar.TKE, self.EnvVar.W, self.UpdVar.W)
    end

    reset_surface_covariance(self, GMV, Case)
    if self.calc_tke
        update_covariance_ED(
            self,
            GMV,
            Case,
            TS,
            GMV.W,
            GMV.W,
            GMV.TKE,
            self.EnvVar.TKE,
            self.EnvVar.W,
            self.EnvVar.W,
            self.UpdVar.W,
            self.UpdVar.W,
        )
    end
    if self.calc_scalar_var
        update_covariance_ED(
            self,
            GMV,
            Case,
            TS,
            GMV.H,
            GMV.H,
            GMV.Hvar,
            self.EnvVar.Hvar,
            self.EnvVar.H,
            self.EnvVar.H,
            self.UpdVar.H,
            self.UpdVar.H,
        )
        update_covariance_ED(
            self,
            GMV,
            Case,
            TS,
            GMV.QT,
            GMV.QT,
            GMV.QTvar,
            self.EnvVar.QTvar,
            self.EnvVar.QT,
            self.EnvVar.QT,
            self.UpdVar.QT,
            self.UpdVar.QT,
        )
        update_covariance_ED(
            self,
            GMV,
            Case,
            TS,
            GMV.H,
            GMV.QT,
            GMV.HQTcov,
            self.EnvVar.HQTcov,
            self.EnvVar.H,
            self.EnvVar.QT,
            self.UpdVar.H,
            self.UpdVar.QT,
        )
        cleanup_covariance(self, GMV)
    end
    return
end

function initialize_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)

    ws = self.wstar
    us = Case.Sur.ustar
    zs = self.base.zi
    grid = get_grid(self)

    reset_surface_covariance(self, GMV, Case)

    if self.calc_tke
        if ws > 0.0
            @inbounds for k in center_indicies(grid)
                z = grid.z_half[k]
                GMV.TKE.values[k] =
                    ws * 1.3 * cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) * sqrt(fmax(1.0 - z / zs, 0.0))
            end
        end
        # TKE initialization from Beare et al, 2006
        if Case.casename == "GABLS"
            @inbounds for k in center_indicies(grid)
                z = grid.z_half[k]
                if (z <= 250.0)
                    GMV.TKE.values[k] = 0.4 * (1.0 - z / 250.0) * (1.0 - z / 250.0) * (1.0 - z / 250.0)
                end
            end

        elseif Case.casename == "Bomex"
            @inbounds for k in center_indicies(grid)
                z = grid.z_half[k]
                if (z <= 2500.0)
                    GMV.TKE.values[k] = 1.0 - z / 3000.0
                end
            end

        elseif Case.casename == "Soares" || Case.casename == "Nieuwstadt"
            @inbounds for k in center_indicies(grid)
                z = grid.z_half[k]
                if (z <= 1600.0)
                    GMV.TKE.values[k] = 0.1 * 1.46 * 1.46 * (1.0 - z / 1600.0)
                end
            end
        end
    end

    if self.calc_scalar_var
        if ws > 0.0
            @inbounds for k in center_indicies(grid)
                z = grid.z_half[k]
                # need to rethink of how to initilize the covarinace profiles - for nowmI took the TKE profile
                GMV.Hvar.values[k] =
                    GMV.Hvar.values[grid.gw] *
                    ws *
                    1.3 *
                    cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) *
                    sqrt(fmax(1.0 - z / zs, 0.0))
                GMV.QTvar.values[k] =
                    GMV.QTvar.values[grid.gw] *
                    ws *
                    1.3 *
                    cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) *
                    sqrt(fmax(1.0 - z / zs, 0.0))
                GMV.HQTcov.values[k] =
                    GMV.HQTcov.values[grid.gw] *
                    ws *
                    1.3 *
                    cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) *
                    sqrt(fmax(1.0 - z / zs, 0.0))
            end
        end

        # TKE initialization from Beare et al, 2006
        if Case.casename == "GABLS"
            @inbounds for k in center_indicies(grid)
                z = grid.z_half[k]
                if (z <= 250.0)
                    GMV.Hvar.values[k] = 0.4 * (1.0 - z / 250.0) * (1.0 - z / 250.0) * (1.0 - z / 250.0)
                end
                GMV.QTvar.values[k] = 0.0
                GMV.HQTcov.values[k] = 0.0
            end
        end
    end
    return
end

function cleanup_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)

    tmp_eps = 1e-18
    grid = get_grid(self)

    @inbounds for k in real_center_indicies(grid)
        if GMV.TKE.values[k] < tmp_eps
            GMV.TKE.values[k] = 0.0
        end
        if GMV.Hvar.values[k] < tmp_eps
            GMV.Hvar.values[k] = 0.0
        end
        if GMV.QTvar.values[k] < tmp_eps
            GMV.QTvar.values[k] = 0.0
        end
        if fabs(GMV.HQTcov.values[k]) < tmp_eps
            GMV.HQTcov.values[k] = 0.0
        end
        if fabs(GMV.W_third_m.values[k]) < tmp_eps
            GMV.W_third_m.values[k] = 0.0
        end
        if fabs(GMV.H_third_m.values[k]) < tmp_eps
            GMV.H_third_m.values[k] = 0.0
        end
        if fabs(GMV.QT_third_m.values[k]) < tmp_eps
            GMV.QT_third_m.values[k] = 0.0
        end
        if self.EnvVar.Hvar.values[k] < tmp_eps
            self.EnvVar.Hvar.values[k] = 0.0
        end
        if self.EnvVar.TKE.values[k] < tmp_eps
            self.EnvVar.TKE.values[k] = 0.0
        end
        if self.EnvVar.QTvar.values[k] < tmp_eps
            self.EnvVar.QTvar.values[k] = 0.0
        end
        if fabs(self.EnvVar.HQTcov.values[k]) < tmp_eps
            self.EnvVar.HQTcov.values[k] = 0.0
        end
    end
end


function compute_covariance_shear(
    self::EDMF_PrognosticTKE,
    GMV::GridMeanVariables,
    Covar::EnvironmentVariable_2m,
    UpdVar1,
    UpdVar2,
    EnvVar1,
    EnvVar2,
)

    grid = get_grid(self)
    dzi = grid.dzi
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues
    diff_var1 = 0.0
    diff_var2 = 0.0
    du = 0.0
    dv = 0.0
    tke_factor = 1.0
    du_high = 0.0
    dv_high = 0.0
    KH = diffusivity_h(self).values
    rho0_half = reference_state(self).rho0_half

    @inbounds for k in real_center_indicies(grid)
        if Covar.name == "tke"
            du_low = du_high
            dv_low = dv_high
            du_high = (GMV.U.values[k + 1] - GMV.U.values[k]) * dzi
            dv_high = (GMV.V.values[k + 1] - GMV.V.values[k]) * dzi
            diff_var2 = (EnvVar2[k] - EnvVar2[k - 1]) * dzi
            diff_var1 = (EnvVar1[k] - EnvVar1[k - 1]) * dzi
            tke_factor = 0.5
            k_eddy = diffusivity_m(self).values[k]
        else
            # Defined correctly only for covariance between half-level variables.
            du_low = 0.0
            dv_low = 0.0
            du_high = 0.0
            dv_high = 0.0
            diff_var2 = interp2pt((EnvVar2[k + 1] - EnvVar2[k]), (EnvVar2[k] - EnvVar2[k - 1])) * dzi
            diff_var1 = interp2pt((EnvVar1[k + 1] - EnvVar1[k]), (EnvVar1[k] - EnvVar1[k - 1])) * dzi
            tke_factor = 1.0
            k_eddy = KH[k]
        end
        Covar.shear[k] =
            tke_factor *
            2.0 *
            (
                rho0_half[k] *
                ae[k] *
                k_eddy *
                (diff_var1 * diff_var2 + pow(interp2pt(du_low, du_high), 2.0) + pow(interp2pt(dv_low, dv_high), 2.0))
            )
    end
    return
end

function compute_covariance_interdomain_src(
    self::EDMF_PrognosticTKE,
    au::UpdraftVariable,
    phi_u::UpdraftVariable,
    psi_u::UpdraftVariable,
    phi_e::EnvironmentVariable,
    psi_e::EnvironmentVariable,
    Covar::EnvironmentVariable_2m,
)

    grid = get_grid(self)
    @inbounds for k in xrange(grid.nzg)
        Covar.interdomain[k] = 0.0
        @inbounds for i in xrange(self.n_updrafts)
            if Covar.name == "tke"
                tke_factor = 0.5
                # TODO: report bug: k-1 for k = 0 yields
                # -1, indexing phi_e.values[-1] yields the
                # _last_ value in the array. This is certainly
                # not intended
                if k ≠ 0
                    phi_diff =
                        interp2pt(phi_u.values[i, k - 1], phi_u.values[i, k]) -
                        interp2pt(phi_e.values[k - 1], phi_e.values[k])
                    psi_diff =
                        interp2pt(psi_u.values[i, k - 1], psi_u.values[i, k]) -
                        interp2pt(psi_e.values[k - 1], psi_e.values[k])
                else # use 0th order approximation
                    phi_diff = phi_u.values[i, k] - phi_e.values[k]
                    psi_diff = psi_u.values[i, k] - psi_e.values[k]
                end
            else
                tke_factor = 1.0
                phi_diff = phi_u.values[i, k] - phi_e.values[k]
                psi_diff = psi_u.values[i, k] - psi_e.values[k]
            end

            Covar.interdomain[k] += tke_factor * au.values[i, k] * (1.0 - au.values[i, k]) * phi_diff * psi_diff
        end
    end
    return
end

function compute_covariance_entr(
    self::EDMF_PrognosticTKE,
    Covar::EnvironmentVariable_2m,
    UpdVar1::UpdraftVariable,
    UpdVar2::UpdraftVariable,
    EnvVar1::EnvironmentVariable,
    EnvVar2::EnvironmentVariable,
    GmvVar1::VariablePrognostic,
    GmvVar2::VariablePrognostic,
)

    grid = get_grid(self)
    ref_state = reference_state(self)
    rho0_half = ref_state.rho0_half
    @inbounds for k in real_center_indicies(grid)
        Covar.entr_gain[k] = 0.0
        Covar.detr_loss[k] = 0.0
        @inbounds for i in xrange(self.n_updrafts)
            if self.UpdVar.Area.values[i, k] > self.minimum_area
                R_up = self.pressure_plume_spacing[i]
                if Covar.name == "tke" # (interp needed)
                    tke_factor = 0.5
                    updvar1 = interp2pt(UpdVar1.values[i, k], UpdVar1.values[i, k - 1])
                    updvar2 = interp2pt(UpdVar2.values[i, k], UpdVar2.values[i, k - 1])
                    envvar1 = interp2pt(EnvVar1.values[k], EnvVar1.values[k - 1])
                    envvar2 = interp2pt(EnvVar2.values[k], EnvVar2.values[k - 1])
                    gmvvar1 = interp2pt(GmvVar1.values[k], GmvVar1.values[k - 1])
                    gmvvar2 = interp2pt(GmvVar2.values[k], GmvVar2.values[k - 1])
                    eps_turb = interp2pt(self.frac_turb_entr_full[i, k], self.frac_turb_entr_full[i, k - 1])
                else
                    tke_factor = 1.0
                    updvar1 = UpdVar1.values[i, k]
                    updvar2 = UpdVar2.values[i, k]
                    envvar1 = EnvVar1.values[k]
                    envvar2 = EnvVar2.values[k]
                    gmvvar1 = GmvVar1.values[k]
                    gmvvar2 = GmvVar2.values[k]
                    eps_turb = self.frac_turb_entr[i, k]
                end

                w_u = interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k])
                dynamic_entr =
                    tke_factor *
                    rho0_half[k] *
                    self.UpdVar.Area.values[i, k] *
                    fabs(w_u) *
                    self.detr_sc[i, k] *
                    (updvar1 - envvar1) *
                    (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    rho0_half[k] *
                    self.UpdVar.Area.values[i, k] *
                    fabs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                Covar.entr_gain[k] += dynamic_entr + turbulent_entr
                Covar.detr_loss[k] +=
                    tke_factor *
                    rho0_half[k] *
                    self.UpdVar.Area.values[i, k] *
                    fabs(w_u) *
                    (self.entr_sc[i, k] + eps_turb) *
                    Covar.values[k]
            end
        end
    end

    return
end

function compute_covariance_detr(self::EDMF_PrognosticTKE, Covar::EnvironmentVariable_2m)
    grid = get_grid(self)
    rho0_half = reference_state(self).rho0_half
    @inbounds for k in real_center_indicies(grid)
        Covar.detr_loss[k] = 0.0
        @inbounds for i in xrange(self.n_updrafts)
            w_u = interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k])
            Covar.detr_loss[k] += self.UpdVar.Area.values[i, k] * fabs(w_u) * self.entr_sc[i, k]
        end
        Covar.detr_loss[k] *= rho0_half[k] * Covar.values[k]
    end
    return
end

function compute_covariance_rain(self::EDMF_PrognosticTKE, TS::TimeStepping, GMV::GridMeanVariables)
    # TODO defined again in compute_covariance_shear and compute_covaraince
    grid = get_grid(self)
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues # area of environment
    rho0_half = reference_state(self).rho0_half

    @inbounds for k in real_center_indicies(grid)
        self.EnvVar.TKE.rain_src[k] = 0.0
        self.EnvVar.Hvar.rain_src[k] = rho0_half[k] * ae[k] * 2.0 * self.EnvThermo.Hvar_rain_dt[k] * TS.dti
        self.EnvVar.QTvar.rain_src[k] = rho0_half[k] * ae[k] * 2.0 * self.EnvThermo.QTvar_rain_dt[k] * TS.dti
        self.EnvVar.HQTcov.rain_src[k] = rho0_half[k] * ae[k] * self.EnvThermo.HQTcov_rain_dt[k] * TS.dti
    end

    return
end


function compute_covariance_dissipation(self::EDMF_PrognosticTKE, Covar::EnvironmentVariable_2m)
    grid = get_grid(self)
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues
    rho0_half = reference_state(self).rho0_half

    @inbounds for k in real_center_indicies(grid)
        Covar.dissipation[k] = (
            rho0_half[k] * ae[k] * Covar.values[k] * pow(fmax(self.EnvVar.TKE.values[k], 0), 0.5) /
            fmax(self.mixing_length[k], 1.0e-3) * self.tke_diss_coeff
        )
    end
    return
end


function compute_tke_advection(self::EDMF_PrognosticTKE)

    grid = get_grid(self)
    ref_state = reference_state(self)
    rho0_half = ref_state.rho0_half
    gw = get_grid(self).gw
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues # area of environment
    drho_ae_we_e_plus = 0.0

    @inbounds for k in real_face_indicies(grid)
        drho_ae_we_e_minus = drho_ae_we_e_plus
        drho_ae_we_e_plus =
            (
                rho0_half[k + 1] *
                ae[k + 1] *
                self.EnvVar.TKE.values[k + 1] *
                (self.EnvVar.W.values[k + 1] + self.EnvVar.W.values[k]) / 2.0 -
                rho0_half[k] *
                ae[k] *
                self.EnvVar.TKE.values[k] *
                (self.EnvVar.W.values[k] + self.EnvVar.W.values[k - 1]) / 2.0
            ) * grid.dzi
        self.tke_advection[k] = interp2pt(drho_ae_we_e_minus, drho_ae_we_e_plus)
    end
    return
end

function compute_tke_transport(self::EDMF_PrognosticTKE)

    grid = get_grid(self)
    dzi = grid.dzi
    gw = grid.gw
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues # area of environment
    dtke_high = 0.0
    drho_ae_K_m_de_plus = 0.0
    KM = diffusivity_m(self).values
    rho0_half = reference_state(self).rho0_half

    @inbounds for k in real_face_indicies(grid)
        drho_ae_K_m_de_low = drho_ae_K_m_de_plus
        drho_ae_K_m_de_plus =
            (
                rho0_half[k + 1] *
                ae[k + 1] *
                KM[k + 1] *
                (self.EnvVar.TKE.values[k + 2] - self.EnvVar.TKE.values[k]) *
                0.5 *
                dzi -
                rho0_half[k] *
                ae[k] *
                KM[k] *
                (self.EnvVar.TKE.values[k + 1] - self.EnvVar.TKE.values[k - 1]) *
                0.5 *
                dzi
            ) * dzi
        self.tke_transport[k] = interp2pt(drho_ae_K_m_de_low, drho_ae_K_m_de_plus)
    end
    return
end

function update_covariance_ED(
    self::EDMF_PrognosticTKE,
    GMV::GridMeanVariables,
    Case::CasesBase,
    TS::TimeStepping,
    GmvVar1::VariablePrognostic,
    GmvVar2::VariablePrognostic,
    GmvCovar::VariableDiagnostic,
    Covar::EnvironmentVariable_2m,
    EnvVar1::EnvironmentVariable,
    EnvVar2::EnvironmentVariable,
    UpdVar1::UpdraftVariable,
    UpdVar2::UpdraftVariable,
)

    grid = get_grid(self)
    gw = grid.gw
    nzg = grid.nzg
    nz = grid.nz
    dzi = grid.dzi
    dti = TS.dti
    Ref = reference_state(self)
    alpha0LL = Ref.alpha0_half[grid.gw]
    zLL = grid.z_half[grid.gw]
    a = pyzeros(nz)
    b = pyzeros(nz)
    c = pyzeros(nz)
    x = pyzeros(nz)
    ae = pyones(nzg) .- self.UpdVar.Area.bulkvalues
    ae_old = pyones(nzg) .- up_sum(self.UpdVar.Area.old)
    rho_ae_K_m = pyzeros(nzg)
    whalf = pyzeros(nzg)
    D_env = 0.0
    KM = diffusivity_m(self)
    KH = diffusivity_h(self)

    @inbounds for k in real_face_indicies(grid)
        if Covar.name == "tke"
            K = KM.values[k]
            Kp = KM.values[k + 1]
        else
            K = KH.values[k]
            Kp = KH.values[k + 1]
        end
        rho_ae_K_m[k] = 0.5 * (ae[k] * K + ae[k + 1] * Kp) * Ref.rho0[k]
        whalf[k] = interp2pt(self.EnvVar.W.values[k - 1], self.EnvVar.W.values[k])
    end
    wu_half = interp2pt(self.UpdVar.W.bulkvalues[gw - 1], self.UpdVar.W.bulkvalues[gw])

    # Not necessary if BCs for variances are applied to environment.
    # if GmvCovar.name=="tke"
    #     GmvCovar.values[gw] =get_surface_tke(Case.Sur.ustar, self.wstar, get_grid(self).z_half[gw], Case.Sur.obukhov_length)
    # elseif GmvCovar.name=="thetal_var"
    #     GmvCovar.values[gw] = get_surface_variance(Case.Sur.rho_hflux * alpha0LL, Case.Sur.rho_hflux * alpha0LL, Case.Sur.ustar, zLL, Case.Sur.obukhov_length)
    # elseif GmvCovar.name=="qt_var"
    #     GmvCovar.values[gw] = get_surface_variance(Case.Sur.rho_qtflux * alpha0LL, Case.Sur.rho_qtflux * alpha0LL, Case.Sur.ustar, zLL, Case.Sur.obukhov_length)
    # elseif GmvCovar.name=="thetal_qt_covar"
    #     GmvCovar.values[gw] = get_surface_variance(Case.Sur.rho_hflux * alpha0LL, Case.Sur.rho_qtflux * alpha0LL, Case.Sur.ustar, zLL, Case.Sur.obukhov_length)
    # self.get_env_covar_from_GMV(self.UpdVar.Area, UpdVar1, UpdVar2, EnvVar1, EnvVar2, Covar, GmvVar1.values, GmvVar2.values, GmvCovar.values)

    Covar_surf = Covar.values[gw]

    @inbounds for kk in xrange(nz)
        k = kk + gw
        D_env = 0.0

        @inbounds for i in xrange(self.n_updrafts)
            if self.UpdVar.Area.values[i, k] > self.minimum_area
                if Covar.name == "tke"
                    turb_entr = interp2pt(self.frac_turb_entr_full[i, k - 1], self.frac_turb_entr_full[i, k])
                else
                    turb_entr = self.frac_turb_entr[i, k]
                end

                R_up = self.pressure_plume_spacing[i]
                wu_half = interp2pt(self.UpdVar.W.values[i, k - 1], self.UpdVar.W.values[i, k])
                D_env += Ref.rho0_half[k] * self.UpdVar.Area.values[i, k] * wu_half * (self.entr_sc[i, k] + turb_entr)
            else
                D_env = 0.0
            end
        end

        a[kk] = (-rho_ae_K_m[k - 1] * dzi * dzi)
        b[kk] = (
            Ref.rho0_half[k] * ae[k] * dti - Ref.rho0_half[k] * ae[k] * whalf[k] * dzi +
            rho_ae_K_m[k] * dzi * dzi +
            rho_ae_K_m[k - 1] * dzi * dzi +
            D_env +
            Ref.rho0_half[k] * ae[k] * self.tke_diss_coeff * sqrt(fmax(self.EnvVar.TKE.values[k], 0)) /
            fmax(self.mixing_length[k], 1.0)
        )
        c[kk] = (Ref.rho0_half[k + 1] * ae[k + 1] * whalf[k + 1] * dzi - rho_ae_K_m[k] * dzi * dzi)
        x[kk] = (
            Ref.rho0_half[k] * ae_old[k] * Covar.values[k] * dti +
            Covar.press[k] +
            Covar.buoy[k] +
            Covar.shear[k] +
            Covar.entr_gain[k] +
            Covar.rain_src[k]
        ) #

        a[0] = 0.0
        b[0] = 1.0
        c[0] = 0.0
        x[0] = Covar_surf

        b[nz - 1] += c[nz - 1]
        c[nz - 1] = 0.0
    end
    x .= tridiag_solve(x, a, b, c)

    @inbounds for kk in xrange(nz)
        k = kk + gw
        if Covar.name == "thetal_qt_covar"
            Covar.values[k] = fmax(x[kk], -sqrt(self.EnvVar.Hvar.values[k] * self.EnvVar.QTvar.values[k]))
            Covar.values[k] = fmin(x[kk], sqrt(self.EnvVar.Hvar.values[k] * self.EnvVar.QTvar.values[k]))
        else
            Covar.values[k] = fmax(x[kk], 0.0)
        end
    end
    set_bcs(Covar, grid)

    get_GMV_CoVar(
        self,
        self.UpdVar.Area,
        UpdVar1,
        UpdVar2,
        EnvVar1,
        EnvVar2,
        Covar,
        GmvVar1.values,
        GmvVar2.values,
        GmvCovar.values,
    )

    return
end

function GMV_third_m(
    self::EDMF_PrognosticTKE,
    Gmv_third_m::VariableDiagnostic,
    env_covar::EnvironmentVariable_2m,
    env_mean::EnvironmentVariable,
    upd_mean::UpdraftVariable,
)

    grid = get_grid(self)
    ae = pyones(grid.nzg) .- self.UpdVar.Area.bulkvalues
    au = self.UpdVar.Area.values

    @inbounds for k in real_center_indicies(grid)
        GMVv_ = ae[k] * env_mean.values[k]
        @inbounds for i in xrange(self.n_updrafts)
            GMVv_ += au[i, k] * upd_mean.values[i, k]
        end

        # TODO: report bug: i used outside of scope.
        # This is only valid (assuming correct) for 1
        # updraft.
        i_last = last(xrange(self.n_updrafts))
        if env_covar.name == "tke"
            Envcov_ =
                -self.horizontal_KM[i_last, k] * (self.EnvVar.W.values[k + 1] - self.EnvVar.W.values[k - 1]) /
                (2.0 * grid.dz)
        else
            Envcov_ = env_covar.values[k]
        end

        Upd_cubed = 0.0
        GMVcov_ = ae[k] * (Envcov_ + (env_mean.values[k] - GMVv_)^2.0)
        @inbounds for i in xrange(self.n_updrafts)
            GMVcov_ += au[i, k] * (upd_mean.values[i, k] - GMVv_)^2.0
            Upd_cubed += au[i, k] * upd_mean.values[i, k]^3
        end

        Gmv_third_m.values[k] =
            Upd_cubed + ae[k] * (env_mean.values[k]^3 + 3.0 * env_mean.values[k] * Envcov_) - GMVv_^3.0 -
            3.0 * GMVcov_ * GMVv_
    end
    Gmv_third_m.values[grid.gw] = 0.0 # this is here as first value is biased with BC area fraction
    return
end
