
function initialize(
    self::EDMF_PrognosticTKE,
    Case::CasesBase,
    GMV::GridMeanVariables,
    Ref::ReferenceState,
    TS::TimeStepping,
)
    param_set = parameter_set(GMV)

    if Case.casename == "DryBubble"
        initialize_DryBubble(self.UpdVar, GMV, Ref)
    else
        initialize(self.UpdVar, GMV)
    end
    decompose_environment(self, GMV, "values")
    saturation_adjustment(self.EnvThermo, self.EnvVar)
    buoyancy(self.UpdThermo, self.UpdVar, self.EnvVar, GMV, self.extrapolate_buoyancy)
    update_inversion(self, GMV, Case.inversion_option)
    self.wstar = get_wstar(Case.Sur.bflux, self.zi)
    microphysics(param_set, self.EnvThermo, self.EnvVar, self.Rain, TS.dt)
    initialize_covariance(self, GMV, Case)
    set_subdomain_bcs(self)
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
    add_profile(Stats, "horiz_K_eddy")
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
    add_profile(Stats, "tke_buoy")
    add_profile(Stats, "tke_dissipation")
    add_profile(Stats, "tke_entr_gain")
    add_profile(Stats, "tke_detr_loss")
    add_profile(Stats, "tke_shear")
    add_profile(Stats, "tke_pressure")
    add_profile(Stats, "tke_interdomain")
    add_profile(Stats, "tke_transport")
    add_profile(Stats, "tke_advection")
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
    mean_detr_sc = center_field(grid)
    massflux = face_field(grid)
    mf_h = face_field(grid)
    mf_qt = face_field(grid)
    mean_frac_turb_entr = center_field(grid)
    mean_horiz_K_eddy = center_field(grid)
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
                mean_frac_turb_entr[k] +=
                    self.UpdVar.Area.values[i, k] * self.frac_turb_entr[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_horiz_K_eddy[k] +=
                    self.UpdVar.Area.values[i, k] * self.horiz_K_eddy[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_sorting_function[k] +=
                    self.UpdVar.Area.values[i, k] * self.sorting_function[i, k] / self.UpdVar.Area.bulkvalues[k]
                mean_b_mix[k] += self.UpdVar.Area.values[i, k] * self.b_mix[i, k] / self.UpdVar.Area.bulkvalues[k]
            end
        end
    end

    write_profile(Stats, "turbulent_entrainment", mean_frac_turb_entr[cinterior])
    write_profile(Stats, "horiz_K_eddy", mean_horiz_K_eddy[cinterior])
    write_profile(Stats, "entrainment_sc", mean_entr_sc[cinterior])
    write_profile(Stats, "detrainment_sc", mean_detr_sc[cinterior])
    write_profile(Stats, "nh_pressure", mean_nh_pressure[cinterior])
    write_profile(Stats, "nh_pressure_adv", mean_nh_pressure_adv[cinterior])
    write_profile(Stats, "nh_pressure_drag", mean_nh_pressure_drag[cinterior])
    write_profile(Stats, "nh_pressure_b", mean_nh_pressure_b[cinterior])
    write_profile(Stats, "asp_ratio", mean_asp_ratio[cinterior])
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
    return
end

# Perform the update of the scheme
function update(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    param_set = parameter_set(self)

    grid = get_grid(self)
    set_old_with_values(self.UpdVar)
    set_updraft_surface_bc(self, GMV, Case)

    # compute RHS
    compute_pressure_plume_spacing(self, GMV, Case)
    compute_updraft_closures(self, GMV, Case)
    compute_eddy_diffusivities_tke(self, GMV, Case)
    compute_GMV_MF(self, GMV, TS)
    compute_covariance_rhs(self, GMV, Case, TS)

    # update
    solve_updraft_velocity_area(self, TS)
    solve_updraft_scalars(self, GMV, TS)
    update_GMV_ED(self, GMV, Case, TS)
    update_covariance(self, GMV, Case, TS)
    update_GMV_turbulence(self, GMV, Case, TS)

    # update microphysics
    microphysics(param_set, self.UpdThermo, self.UpdVar, self.Rain, TS.dt) # causes division error in dry bubble first time step
    microphysics(param_set, self.EnvThermo, self.EnvVar, self.Rain, TS.dt) # saturation adjustment + rain creation
    update_total_precip_sources(self.UpdThermo)

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

    # set values
    set_values_with_new(self.UpdVar)
    zero_area_fraction_cleanup(self, GMV)
    decompose_environment(self, GMV, "values")
    saturation_adjustment(self.EnvThermo, self.EnvVar)
    buoyancy(self.UpdThermo, self.UpdVar, self.EnvVar, GMV, self.extrapolate_buoyancy)
    set_subdomain_bcs(self)
    clear_precip_sources(self.UpdThermo)
    return
end

function update_GMV_turbulence(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    @inbounds for k in real_center_indicies(self.Gr)
        GMV.H.tendencies[k] += (GMV.H.new[k] - GMV.H.values[k]) * TS.dti
        GMV.QT.tendencies[k] += (GMV.QT.new[k] - GMV.QT.values[k]) * TS.dti
        GMV.U.tendencies[k] += (GMV.U.new[k] - GMV.U.values[k]) * TS.dti
        GMV.V.tendencies[k] += (GMV.V.new[k] - GMV.V.values[k]) * TS.dti
    end

    return
end

function compute_mixing_length(self, obukhov_length, ustar, GMV::GridMeanVariables)

    param_set = parameter_set(self)
    grid = get_grid(self)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    ref_state = reference_state(self)
    kc_surf = kc_surface(grid)
    tau = get_mixing_tau(self.zi, self.wstar)
    l = zeros(3)
    ϵ = 1.0e-9 # Epsilon to avoid zero
    ω_pr = CPEDMF.ω_pr(param_set)

    @inbounds for k in real_center_indicies(grid)
        z_ = grid.z_half[k]
        # kz scale (surface layer)
        if obukhov_length < 0.0 #unstable
            l2 =
                vkb * z_ / (sqrt(self.EnvVar.TKE.values[kc_surf] / ustar / ustar) * self.tke_ed_coeff) *
                min((1.0 - 100.0 * z_ / obukhov_length)^0.2, 1.0 / vkb)
        else # neutral or stable
            l2 = vkb * z_ / (sqrt(self.EnvVar.TKE.values[kc_surf] / ustar / ustar) * self.tke_ed_coeff)
        end

        # Buoyancy-shear-subdomain exchange-dissipation TKE equilibrium scale
        U_cut = cut(GMV.U.values, grid, k)
        V_cut = cut(GMV.V.values, grid, k)
        w_dual = dual_faces(self.EnvVar.W.values, grid, k)
        ∇U = c∇(U_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇V = c∇(V_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇w = ∇f2c(w_dual, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        Shear² = ∇U^2 + ∇V^2 + ∇w^2
        qt_dry = self.EnvThermo.qt_dry[k]
        th_dry = self.EnvThermo.th_dry[k]
        t_cloudy = self.EnvThermo.t_cloudy[k]
        qv_cloudy = self.EnvThermo.qv_cloudy[k]
        qt_cloudy = self.EnvThermo.qt_cloudy[k]
        th_cloudy = self.EnvThermo.th_cloudy[k]
        lh = TD.latent_heat_vapor(param_set, t_cloudy)
        #TODO - add ql and qi to cloudy output to better compute cpm here
        cpm = TD.cp_m(param_set, TD.PhasePartition(qt_cloudy, qt_cloudy - qv_cloudy, 0.0))

        QT_cut = cut(self.EnvVar.QT.values, grid, k)
        grad_qt = c∇(QT_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        THL_cut = cut(self.EnvVar.H.values, grid, k)
        grad_thl = c∇(THL_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))

        # g/theta_ref
        prefactor = g * (Rd / ref_state.alpha0_half[k] / ref_state.p0_half[k]) * TD.exner_given_pressure(param_set, ref_state.p0_half[k])
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
        ∂b∂θ_l = grad_thl * d_buoy_thetal_total
        ∂b∂q_tot = grad_qt * d_buoy_qt_total

        ∇Ri = gradient_Richardson_number(∂b∂θ_l, ∂b∂q_tot, Shear², ϵ)
        self.prandtl_nvec[k] = turbulent_Prandtl_number(obukhov_length, ∇Ri, prandtl_number(self), ω_pr)

        # Production/destruction terms
        a =
            self.tke_ed_coeff *
            (Shear² - ∂b∂θ_l / self.prandtl_nvec[k] - ∂b∂q_tot / self.prandtl_nvec[k]) *
            sqrt(self.EnvVar.TKE.values[k])
        # Dissipation term
        c_neg = self.tke_diss_coeff * self.EnvVar.TKE.values[k] * sqrt(self.EnvVar.TKE.values[k])
        # Subdomain exchange term
        self.b[k] = 0.0
        for nn in xrange(self.n_updrafts)
            wc_upd_nn = (self.UpdVar.W.values[nn, k] + self.UpdVar.W.values[nn, k - 1]) / 2.0
            wc_env = (self.EnvVar.W.values[k] + self.EnvVar.W.values[k - 1]) / 2.0
            self.b[k] +=
                self.UpdVar.Area.values[nn, k] * wc_upd_nn * self.detr_sc[nn, k] /
                (1.0 - self.UpdVar.Area.bulkvalues[k]) *
                ((wc_upd_nn - wc_env) * (wc_upd_nn - wc_env) / 2.0 - self.EnvVar.TKE.values[k]) -
                self.UpdVar.Area.values[nn, k] *
                wc_upd_nn *
                (wc_upd_nn - wc_env) *
                self.frac_turb_entr[nn, k] *
                wc_env / (1.0 - self.UpdVar.Area.bulkvalues[k])
        end

        if abs(a) > ϵ && 4.0 * a * c_neg > -self.b[k] * self.b[k]
            self.l_entdet[k] = max(-self.b[k] / 2.0 / a + sqrt(self.b[k] * self.b[k] + 4.0 * a * c_neg) / 2.0 / a, 0.0)
        elseif abs(a) < ϵ && abs(self.b[k]) > ϵ
            self.l_entdet[k] = c_neg / self.b[k]
        end
        l3 = self.l_entdet[k]

        # Limiting stratification scale (Deardorff, 1976)
        p0_cut = cut(ref_state.p0_half, grid, k)
        T_cut = cut(self.EnvVar.T.values, grid, k)
        QT_cut = cut(self.EnvVar.QT.values, grid, k)
        QL_cut = cut(self.EnvVar.QL.values, grid, k)
        thv_cut = theta_virt_c.(p0_cut, T_cut, QT_cut, QL_cut)
        grad_thv = c∇(thv_cut, grid, k; bottom = SetGradient(0), top = Extrapolate())

        p0_k = ref_state.p0_half[k]
        T_k = self.EnvVar.T.values[k]
        QT_k = self.EnvVar.QT.values[k]
        QL_k = self.EnvVar.QL.values[k]
        thv = theta_virt_c(p0_k, T_k, QT_k, QL_k)

        # Effective static stability using environmental mean.
        # Set lambda for now to environmental cloud_fraction (TBD: Rain)
        # also TODO - assumes no ice
        grad_th_eff =
            (1.0 - self.EnvVar.cloud_fraction.values[k]) * grad_thv +
            self.EnvVar.cloud_fraction.values[k] * (
                1.0 / exp(
                    -TD.latent_heat_vapor(param_set, self.EnvVar.T.values[k]) *
                     self.EnvVar.QL.values[k] /
                     TD.cp_m(param_set, TD.PhasePartition(self.EnvVar.QT.values[k], self.EnvVar.QL.values[k], 0.0)) /
                     self.EnvVar.T.values[k],
                ) * (
                    (1.0 + (eps_vi - 1.0) * self.EnvVar.QT.values[k]) * grad_thl +
                    (eps_vi - 1.0) * self.EnvVar.H.values[k] * grad_qt
                )
            )
        N = sqrt(max(g / thv * grad_th_eff, 0.0))
        if N > 0.0
            l1 = min(sqrt(max(self.static_stab_coeff * self.EnvVar.TKE.values[k], 0.0)) / N, 1.0e6)
        else
            l1 = 1.0e6
        end

        l[1] = l1
        l[2] = l3
        l[3] = l2

        j = 1
        while (j <= length(l))
            if l[j] < ϵ || l[j] > 1.0e6
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
    compute_mixing_length(self, Case.Sur.obukhov_length, Case.Sur.ustar, GMV)
    KM = diffusivity_m(self)
    KH = diffusivity_h(self)
    @inbounds for k in real_center_indicies(grid)
        lm = self.mixing_length[k]
        pr = self.prandtl_nvec[k]
        KM.values[k] = self.tke_ed_coeff * lm * sqrt(max(self.EnvVar.TKE.values[k], 0.0))
        KH.values[k] = KM.values[k] / pr
    end
    return
end

function set_updraft_surface_bc(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    grid = get_grid(self)
    kc_surf = kc_surface(grid)
    update_inversion(self, GMV, Case.inversion_option)
    self.wstar = get_wstar(Case.Sur.bflux, self.zi)

    dzi = grid.dzi
    zLL = grid.z_half[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    alpha0LL = reference_state(self).alpha0_half[kc_surf]
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
        self.h_surface_bc[i] = (GMV.H.values[kc_surf] + surface_scalar_coeff * sqrt(h_var))
        self.qt_surface_bc[i] = (GMV.QT.values[kc_surf] + surface_scalar_coeff * sqrt(qt_var))
    end
    return
end

function reset_surface_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    flux1 = Case.Sur.rho_hflux
    flux2 = Case.Sur.rho_qtflux
    grid = get_grid(self)
    kc_surf = kc_surface(grid)
    ref_state = reference_state(self)
    zLL = grid.z_half[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    alpha0LL = ref_state.alpha0_half[kc_surf]
    self.EnvVar.TKE.values[kc_surf] =
        get_surface_tke(Case.Sur.ustar, self.wstar, grid.z_half[kc_surf], Case.Sur.obukhov_length)
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

    self.EnvVar.Hvar.values[kc_surf] = get_surface_variance(flux1 * alpha0LL, flux1 * alpha0LL, ustar, zLL, oblength)
    self.EnvVar.QTvar.values[kc_surf] = get_surface_variance(flux2 * alpha0LL, flux2 * alpha0LL, ustar, zLL, oblength)
    self.EnvVar.HQTcov.values[kc_surf] = get_surface_variance(flux1 * alpha0LL, flux2 * alpha0LL, ustar, zLL, oblength)
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
    return
end

# Find values of environmental variables by subtracting updraft values from grid mean values
# whichvals used to check which substep we are on--correspondingly use "GMV.SomeVar.value" (last timestep value)
# or GMV.SomeVar.mf_update (GMV value following massflux substep)
function decompose_environment(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, whichvals)

    # first make sure the "bulkvalues" of the updraft variables are updated
    set_means(self.UpdVar, GMV)
    grid = get_grid(self)

    if whichvals == "values"
        @inbounds for k in real_face_indicies(grid)
            val1 = 1.0 / (1.0 - self.UpdVar.Area.bulkvalues[k])
            val2 = self.UpdVar.Area.bulkvalues[k] * val1

            self.EnvVar.Area.values[k] = 1.0 - self.UpdVar.Area.bulkvalues[k]
            self.EnvVar.QT.values[k] = max(val1 * GMV.QT.values[k] - val2 * self.UpdVar.QT.bulkvalues[k], 0.0) #Yair - this is here to prevent negative QT
            self.EnvVar.H.values[k] = val1 * GMV.H.values[k] - val2 * self.UpdVar.H.bulkvalues[k]
            # Have to account for staggering of W--interpolate area fraction to the "full" grid points
            # Assuming GMV.W = 0!
            au_full = 0.5 * (self.UpdVar.Area.bulkvalues[k + 1] + self.UpdVar.Area.bulkvalues[k])
            self.EnvVar.W.values[k] = -au_full / (1.0 - au_full) * self.UpdVar.W.bulkvalues[k]
        end
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


    elseif whichvals == "mf_update"
        # same as above but replace GMV.SomeVar.values with GMV.SomeVar.mf_update

        @inbounds for k in real_face_indicies(grid)
            val1 = 1.0 / (1.0 - self.UpdVar.Area.bulkvalues[k])
            val2 = self.UpdVar.Area.bulkvalues[k] * val1

            self.EnvVar.QT.values[k] = max(val1 * GMV.QT.mf_update[k] - val2 * self.UpdVar.QT.bulkvalues[k], 0.0)#Yair - this is here to prevent negative QT
            self.EnvVar.H.values[k] = val1 * GMV.H.mf_update[k] - val2 * self.UpdVar.H.bulkvalues[k]
            # Have to account for staggering of W
            # Assuming GMV.W = 0!
            au_full = 0.5 * (self.UpdVar.Area.bulkvalues[k + 1] + self.UpdVar.Area.bulkvalues[k])
            self.EnvVar.W.values[k] = -au_full / (1.0 - au_full) * self.UpdVar.W.bulkvalues[k]
        end

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
    ae = 1 .- au.bulkvalues
    tke_factor = 1.0
    is_tke = covar_e.name == "tke"

    if is_tke
        @inbounds for k in face_indicies(grid)
            tke_factor = 0.5
            # TODO: report bug: k-1 for k = 0 yields
            # -1, indexing phi_e.values[-1] yields the
            # _last_ value in the array. This is certainly
            # not intended
            if k ≠ 1
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
                if k ≠ 1
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
    ae = 1 .- au.bulkvalues
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

function compute_updraft_closures(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)
    grid = get_grid(self)
    param_set = parameter_set(self)
    ref_state = reference_state(self)

    upd_cloud_diagnostics(self.UpdVar, ref_state)

    @inbounds for k in real_center_indicies(grid)
        @inbounds for i in xrange(self.n_updrafts)
            # entrainment
            if self.UpdVar.Area.values[i, k] > 0.0
                w_upd = interpf2c(self.UpdVar.W.values, grid, k, i)
                # compute dMdz at half levels
                gmv_w_k = interpf2c(GMV.W.values, grid, k)
                gmv_w_km = interpf2c(GMV.W.values, grid, k - 1)
                upd_w_km = interpf2c(self.UpdVar.W.values, grid, k - 1, i)
                M = self.UpdVar.Area.values[i, k] * (w_upd - gmv_w_k)
                Mm = self.UpdVar.Area.values[i, k - 1] * (upd_w_km - gmv_w_km)
                dMdz = (M - Mm) * grid.dzi
                w_min = 0.001

                εδ_model = MoistureDeficitEntr(;
                    q_liq_up = self.UpdVar.QL.values[i, k],
                    q_liq_en = self.EnvVar.QL.values[k],
                    w_up = interpf2c(self.UpdVar.W.values, grid, k, i),
                    w_en = interpf2c(self.EnvVar.W.values, grid, k),
                    b_up = self.UpdVar.B.values[i, k],
                    b_en = self.EnvVar.B.values[k],
                    tke = self.EnvVar.TKE.values[k],
                    dMdz = dMdz,
                    M = M,
                    a_up = self.UpdVar.Area.values[i, k],
                    a_en = self.EnvVar.Area.values[k],
                    R_up = self.pressure_plume_spacing[i],
                    RH_up = self.UpdVar.RH.values[i, k],
                    RH_en = self.EnvVar.RH.values[k],
                )

                er = entr_detr(param_set, εδ_model)
                self.entr_sc[i, k] = er.ε_dyn
                self.detr_sc[i, k] = er.δ_dyn
                self.frac_turb_entr[i, k] = er.ε_turb
                self.horiz_K_eddy[i, k] = er.K_ε
            else
                self.entr_sc[i, k] = 0.0
                self.detr_sc[i, k] = 0.0
                self.frac_turb_entr[i, k] = 0.0
                self.horiz_K_eddy[i, k] = 0.0
            end

            # pressure
            a_kfull = interpc2f(self.UpdVar.Area.values, grid, k, i)
            if a_kfull > 0.0
                b_kfull = interpc2f(self.UpdVar.B.values, grid, k, i)
                ∇w_up = ∇f2c(self.UpdVar.W.values, grid, k, i)
                asp_ratio = 1.0
                self.nh_pressure_b[i, k], self.nh_pressure_adv[i, k], self.nh_pressure_drag[i, k] =
                    perturbation_pressure(
                        param_set,
                        self.UpdVar.updraft_top[i],
                        a_kfull,
                        b_kfull,
                        ref_state.rho0[k],
                        self.UpdVar.W.values[i, k],
                        ∇w_up,
                        self.EnvVar.W.values[k],
                        asp_ratio,
                    )
            else
                self.nh_pressure_b[i, k] = 0.0
                self.nh_pressure_adv[i, k] = 0.0
                self.nh_pressure_drag[i, k] = 0.0
            end
            self.nh_pressure[i, k] = self.nh_pressure_b[i, k] + self.nh_pressure_adv[i, k] + self.nh_pressure_drag[i, k]
        end
    end
    return
end

function compute_pressure_plume_spacing(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)

    param_set = parameter_set(self)
    H_up_min = CPEDMF.H_up_min(param_set)
    @inbounds for i in xrange(self.n_updrafts)
        self.pressure_plume_spacing[i] =
            max(self.aspect_ratio * self.UpdVar.updraft_top[i], H_up_min * self.aspect_ratio)
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
            end
        end

        if sum(self.UpdVar.Area.values[:, k]) == 0.0
            self.EnvVar.W.values[k] = GMV.W.values[k]
            self.EnvVar.B.values[k] = GMV.B.values[k]
            self.EnvVar.H.values[k] = GMV.H.values[k]
            self.EnvVar.QT.values[k] = GMV.QT.values[k]
            self.EnvVar.T.values[k] = GMV.T.values[k]
            self.EnvVar.QL.values[k] = GMV.QL.values[k]
        end
    end

    return
end


function set_subdomain_bcs(self::EDMF_PrognosticTKE)

    grid = get_grid(self)
    set_bcs(self.UpdVar.W, grid)
    set_bcs(self.UpdVar.Area, grid)
    set_bcs(self.UpdVar.H, grid)
    set_bcs(self.UpdVar.QT, grid)
    set_bcs(self.UpdVar.T, grid)
    set_bcs(self.UpdVar.B, grid)

    set_bcs(self.EnvVar.W, grid)
    set_bcs(self.EnvVar.H, grid)
    set_bcs(self.EnvVar.T, grid)
    set_bcs(self.EnvVar.QL, grid)
    set_bcs(self.EnvVar.QT, grid)

    return
end

function solve_updraft_velocity_area(self::EDMF_PrognosticTKE, TS::TimeStepping)
    grid = get_grid(self)
    ref_state = reference_state(self)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dzi = grid.dzi
    dti_ = 1.0 / TS.dt
    dt_ = 1.0 / dti_

    @inbounds for i in xrange(self.n_updrafts)
        self.entr_sc[i, kc_surf] = self.entr_surface_bc
        self.detr_sc[i, kc_surf] = self.detr_surface_bc
        self.UpdVar.W.new[i, kf_surf] = self.w_surface_bc[i]
        self.UpdVar.Area.new[i, kc_surf] = self.area_surface_bc[i]
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
                max(dt_ * (adv + entr_term + detr_term) + self.UpdVar.Area.values[i, k + 1], 0.0)

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
                entr_w = interp2pt(
                    self.entr_sc[i, k] + self.frac_turb_entr[i, k],
                    self.entr_sc[i, k + 1] + self.frac_turb_entr[i, k + 1],
                )
                detr_w = interp2pt(
                    self.detr_sc[i, k] + self.frac_turb_entr[i, k],
                    self.detr_sc[i, k + 1] + self.frac_turb_entr[i, k + 1],
                )
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
                    (entr_w * self.EnvVar.W.values[k] - detr_w * self.UpdVar.W.values[i, k])
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

function solve_updraft_scalars(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, TS::TimeStepping)

    param_set = parameter_set(GMV)

    grid = get_grid(self)
    ref_state = reference_state(self)
    dzi = grid.dzi
    dti_ = 1.0 / TS.dt

    @inbounds for i in xrange(self.n_updrafts)
        # starting from the bottom do entrainment at each level
        @inbounds for k in real_center_indicies(grid)
            if is_surface_center(grid, k)
                # at the surface
                if self.UpdVar.Area.new[i, k] >= self.minimum_area
                    self.UpdVar.H.new[i, k] = self.h_surface_bc[i]
                    self.UpdVar.QT.new[i, k] = self.qt_surface_bc[i]
                else
                    self.UpdVar.H.new[i, k] = GMV.H.values[k]
                    self.UpdVar.QT.new[i, k] = GMV.QT.values[k]
                end

                # saturation adjustment
                ts = TD.PhaseEquil_pθq(
                    param_set,
                    ref_state.p0_half[k],
                    self.UpdVar.H.new[i, k],
                    self.UpdVar.QT.new[i, k],
                )
                ql = TD.liquid_specific_humidity(ts)
                T = TD.air_temperature(ts)
                self.UpdVar.QL.new[i, k] = ql
                self.UpdVar.T.new[i, k] = T
                continue
            end

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

                adv = (m_k * self.UpdVar.H.values[i, k] - m_km * self.UpdVar.H.values[i, k - 1]) * dzi
                entr = (self.entr_sc[i, k] + self.frac_turb_entr[i, k]) * self.EnvVar.H.values[k]
                detr = (self.detr_sc[i, k] + self.frac_turb_entr[i, k]) * self.UpdVar.H.values[i, k]
                self.UpdVar.H.new[i, k] =
                    (
                        ref_state.rho0_half[k] * self.UpdVar.Area.values[i, k] * dti_ * self.UpdVar.H.values[i, k] -
                        adv + m_k * (entr - detr)
                    ) / (ref_state.rho0_half[k] * self.UpdVar.Area.new[i, k] * dti_)

                adv = (m_k * self.UpdVar.QT.values[i, k] - m_km * self.UpdVar.QT.values[i, k - 1]) * dzi
                entr = (self.entr_sc[i, k] + self.frac_turb_entr[i, k]) * self.EnvVar.QT.values[k]
                detr = (self.detr_sc[i, k] + self.frac_turb_entr[i, k]) * self.UpdVar.QT.values[i, k]
                self.UpdVar.QT.new[i, k] = max(
                    (
                        ref_state.rho0_half[k] * self.UpdVar.Area.values[i, k] * dti_ * self.UpdVar.QT.values[i, k] - adv + m_k * (entr - detr)
                    ) / (ref_state.rho0_half[k] * self.UpdVar.Area.new[i, k] * dti_),
                    0.0,
                )

            else
                self.UpdVar.H.new[i, k] = GMV.H.values[k]
                self.UpdVar.QT.new[i, k] = GMV.QT.values[k]
            end

            # saturation adjustment
            ts = TD.PhaseEquil_pθq(param_set, ref_state.p0_half[k], self.UpdVar.H.new[i, k], self.UpdVar.QT.new[i, k])
            ql = TD.liquid_specific_humidity(ts)
            T = TD.air_temperature(ts)
            self.UpdVar.QL.new[i, k] = ql
            self.UpdVar.T.new[i, k] = T
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

function compute_GMV_MF(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, TS::TimeStepping)
    grid = get_grid(self)
    ref_state = reference_state(self)
    rho0 = ref_state.rho0
    kf_surf = kf_surface(grid)
    mf_tend_h = 0.0
    mf_tend_qt = 0.0
    ae = 1 .- self.UpdVar.Area.bulkvalues # area of environment

    self.massflux_h .= 0.0
    self.massflux_qt .= 0.0

    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in xrange(self.n_updrafts)
        self.m[i, kf_surf] = 0.0
        @inbounds for k in real_face_indicies(grid)
            a = interp2pt(self.UpdVar.Area.values[i, k], self.UpdVar.Area.values[i, k + 1])
            self.m[i, k] =
                rho0[k] * a * interp2pt(ae[k], ae[k + 1]) * (self.UpdVar.W.values[i, k] - self.EnvVar.W.values[k])
        end
    end


    self.massflux_h[kf_surf] = 0.0
    self.massflux_qt[kf_surf] = 0.0

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
    nzg = grid.nzg
    nz = grid.nz
    dzi = grid.dzi
    ref_state = reference_state(self)
    a = center_field(grid) # for tridiag solver
    b = center_field(grid) # for tridiag solver
    c = center_field(grid) # for tridiag solver
    x = center_field(grid) # for tridiag solver
    ae = 1 .- self.UpdVar.Area.bulkvalues # area of environment
    rho_ae_K = face_field(grid)
    KM = diffusivity_m(self)
    KH = diffusivity_h(self)
    @inbounds for k in real_face_indicies(grid)
        rho_ae_K[k] = 0.5 * (ae[k] * KH.values[k] + ae[k + 1] * KH.values[k + 1]) * ref_state.rho0[k]
    end

    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    construct_tridiag_diffusion(grid, TS.dt, rho_ae_K, ref_state.rho0_half, ae, a, b, c)
    # Solve QT
    @inbounds for k in real_center_indicies(grid)
        x[k] = self.EnvVar.QT.values[k]
        if is_surface_center(grid, k)
            x[k] = x[k] + TS.dt * Case.Sur.rho_qtflux * dzi * ref_state.alpha0_half[k] / ae[k]
        end
    end
    cinterior = real_center_indicies(grid)
    x[cinterior] .= tridiag_solve(x[cinterior], a[cinterior], b[cinterior], c[cinterior])

    @inbounds for k in real_center_indicies(grid)
        GMV.QT.new[k] = max(
            GMV.QT.mf_update[k] +
            ae[k] * (x[k] - self.EnvVar.QT.values[k]) +
            self.EnvThermo.prec_source_qt[k] +
            self.rainphysics.rain_evap_source_qt[k],
            0.0,
        )
        self.diffusive_tendency_qt[k] = (GMV.QT.new[k] - GMV.QT.mf_update[k]) * TS.dti
    end
    # get the diffusive flux
    @inbounds for k in real_center_indicies(grid)
        if is_surface_center(grid, k)
            self.diffusive_flux_qt[k] = interp2pt(
                Case.Sur.rho_qtflux,
                -rho_ae_K[k] * dzi * (self.EnvVar.QT.values[k + 1] - self.EnvVar.QT.values[k]),
            )
        else
            self.diffusive_flux_qt[k] =
                -0.5 *
                ref_state.rho0_half[k] *
                ae[k] *
                KH.values[k] *
                dzi *
                (self.EnvVar.QT.values[k + 1] - self.EnvVar.QT.values[k - 1])
        end
    end

    # Solve H
    @inbounds for k in real_center_indicies(grid)
        x[k] = self.EnvVar.H.values[k]
        if is_surface_center(grid, k)
            x[k] = x[k] + TS.dt * Case.Sur.rho_hflux * dzi * ref_state.alpha0_half[k] / ae[k]
        end
    end
    x[cinterior] .= tridiag_solve(x[cinterior], a[cinterior], b[cinterior], c[cinterior])
    @inbounds for k in real_center_indicies(grid)
        GMV.H.new[k] =
            GMV.H.mf_update[k] +
            ae[k] * (x[k] - self.EnvVar.H.values[k]) +
            self.EnvThermo.prec_source_h[k] +
            self.rainphysics.rain_evap_source_h[k]
        self.diffusive_tendency_h[k] = (GMV.H.new[k] - GMV.H.mf_update[k]) * TS.dti
    end
    # get the diffusive flux
    @inbounds for k in real_center_indicies(grid)
        if is_surface_center(grid, k)
            self.diffusive_flux_h[k] = interp2pt(
                Case.Sur.rho_hflux,
                -rho_ae_K[k] * dzi * (self.EnvVar.H.values[k + 1] - self.EnvVar.H.values[k]),
            )
        else
            self.diffusive_flux_h[k] =
                -0.5 *
                ref_state.rho0_half[k] *
                ae[k] *
                KH.values[k] *
                dzi *
                (self.EnvVar.H.values[k + 1] - self.EnvVar.H.values[k - 1])
        end
    end

    # Solve U
    @inbounds for k in real_face_indicies(grid)
        rho_ae_K[k] = 0.5 * (ae[k] * KM.values[k] + ae[k + 1] * KM.values[k + 1]) * ref_state.rho0[k]
    end

    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    construct_tridiag_diffusion(grid, TS.dt, rho_ae_K, ref_state.rho0_half, ae, a, b, c)
    @inbounds for k in real_center_indicies(grid)
        x[k] = GMV.U.values[k]
        if is_surface_center(grid, k)
            x[k] = x[k] + TS.dt * Case.Sur.rho_uflux * dzi * ref_state.alpha0_half[k] / ae[k]
        end
    end
    x[cinterior] .= tridiag_solve(x[cinterior], a[cinterior], b[cinterior], c[cinterior])

    @inbounds for k in real_center_indicies(grid)
        GMV.U.new[k] = x[k]
    end
    @inbounds for k in real_center_indicies(grid)
        if is_surface_center(grid, k)
            self.diffusive_flux_u[k] =
                interp2pt(Case.Sur.rho_uflux, -rho_ae_K[k] * dzi * (GMV.U.values[k + 1] - GMV.U.values[k]))
        else
            self.diffusive_flux_u[k] =
                -0.5 * ref_state.rho0_half[k] * ae[k] * KM.values[k] * dzi * (GMV.U.values[k + 1] - GMV.U.values[k - 1])
        end
    end

    # Solve V
    @inbounds for k in real_center_indicies(grid)
        x[k] = GMV.V.values[k]
        if is_surface_center(grid, k)
            x[k] = x[k] + TS.dt * Case.Sur.rho_vflux * dzi * ref_state.alpha0_half[k] / ae[k]
        end
    end
    x[cinterior] .= tridiag_solve(x[cinterior], a[cinterior], b[cinterior], c[cinterior])
    @inbounds for k in real_center_indicies(grid)
        GMV.V.new[k] = x[k]
    end
    @inbounds for k in real_center_indicies(grid)
        if is_surface_center(grid, k)
            self.diffusive_flux_v[k] =
                interp2pt(Case.Sur.rho_vflux, -rho_ae_K[k] * dzi * (GMV.V.values[k + 1] - GMV.V.values[k]))
        else
            self.diffusive_flux_v[k] =
                -0.5 * ref_state.rho0_half[k] * ae[k] * KM.values[k] * dzi * (GMV.V.values[k + 1] - GMV.V.values[k - 1])
        end
    end
    set_bcs(GMV.QT, grid)
    set_bcs(GMV.H, grid)
    set_bcs(GMV.U, grid)
    set_bcs(GMV.V, grid)

    return
end

function compute_tke_buoy(self::EDMF_PrognosticTKE, GMV::GridMeanVariables)

    param_set = parameter_set(self)
    grid = get_grid(self)
    ref_state = reference_state(self)
    param_set = parameter_set(GMV)
    g = CPP.grav(param_set)
    grad_thl_minus = 0.0
    grad_qt_minus = 0.0
    grad_thl_plus = 0
    grad_qt_plus = 0
    ae = 1 .- self.UpdVar.Area.bulkvalues
    KH = diffusivity_h(self).values

    # Note that source terms at the first center interior grid point are not really used because
    # that is where tke boundary condition is enforced (according to MO similarity). Thus here I
    # am being sloppy about lowest grid point
    @inbounds for k in real_center_indicies(grid)
        qt_dry = self.EnvThermo.qt_dry[k]
        th_dry = self.EnvThermo.th_dry[k]
        t_cloudy = self.EnvThermo.t_cloudy[k]
        qv_cloudy = self.EnvThermo.qv_cloudy[k]
        qt_cloudy = self.EnvThermo.qt_cloudy[k]
        th_cloudy = self.EnvThermo.th_cloudy[k]

        lh = TD.latent_heat_vapor(param_set, t_cloudy)        # TODO pass in phase partition for cloudy conditions
        cpm = TD.cp_m(param_set, TD.PhasePartition(qt_cloudy, qt_cloudy-qv_cloudy, 0.0))

        grad_thl_minus = grad_thl_plus
        grad_qt_minus = grad_qt_plus
        grad_thl_plus = (self.EnvVar.H.values[k + 1] - self.EnvVar.H.values[k]) * grid.dzi
        grad_qt_plus = (self.EnvVar.QT.values[k + 1] - self.EnvVar.QT.values[k]) * grid.dzi
        prefactor = Rd * TD.exner_given_pressure(param_set, ref_state.p0_half[k]) / ref_state.p0_half[k]
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

        pp = TD.PhasePartition(GMV.QT.values[k], GMV.QL.values[k], 0.0)
        qv = TD.vapor_specific_humidity(pp)

        GMV.B.values[k] = (
            self.UpdVar.Area.bulkvalues[k] * self.UpdVar.B.bulkvalues[k] +
            (1.0 - self.UpdVar.Area.bulkvalues[k]) * self.EnvVar.B.values[k]
        )
    end

    return
end

function compute_covariance_rhs(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

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

    reset_surface_covariance(self, GMV, Case)
    return
end

function update_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)
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
    return
end

function initialize_covariance(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, Case::CasesBase)

    ws = self.wstar
    us = Case.Sur.ustar
    zs = self.zi
    grid = get_grid(self)
    kc_surf = kc_surface(grid)

    reset_surface_covariance(self, GMV, Case)


    if ws > 0.0
        @inbounds for k in center_indicies(grid)
            z = grid.z_half[k]
            GMV.TKE.values[k] =
                ws * 1.3 * cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) * sqrt(max(1.0 - z / zs, 0.0))
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

    if ws > 0.0
        @inbounds for k in center_indicies(grid)
            z = grid.z_half[k]
            # need to rethink of how to initilize the covarinace profiles - for now took the TKE profile
            GMV.Hvar.values[k] =
                GMV.Hvar.values[kc_surf] *
                ws *
                1.3 *
                cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) *
                sqrt(max(1.0 - z / zs, 0.0))
            GMV.QTvar.values[k] =
                GMV.QTvar.values[kc_surf] *
                ws *
                1.3 *
                cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) *
                sqrt(max(1.0 - z / zs, 0.0))
            GMV.HQTcov.values[k] =
                GMV.HQTcov.values[kc_surf] *
                ws *
                1.3 *
                cbrt((us * us * us) / (ws * ws * ws) + 0.6 * z / zs) *
                sqrt(max(1.0 - z / zs, 0.0))
        end
    else
        @inbounds for k in center_indicies(grid)
            GMV.Hvar.values[k] = 0.0
            GMV.QTvar.values[k] = 0.0
            GMV.HQTcov.values[k] = 0.0
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

    @inbounds for k in center_indicies(grid)
        self.EnvVar.TKE.values[k] = GMV.TKE.values[k]
        self.EnvVar.Hvar.values[k] = GMV.Hvar.values[k]
        self.EnvVar.QTvar.values[k] = GMV.QTvar.values[k]
        self.EnvVar.HQTcov.values[k] = GMV.HQTcov.values[k]
    end
    return
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
    ae = 1 .- self.UpdVar.Area.bulkvalues
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
                (diff_var1 * diff_var2 + interp2pt(du_low, du_high)^2 + interp2pt(dv_low, dv_high)^2)
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
                if k ≠ 1
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
                    eps_turb = interp2pt(self.frac_turb_entr[i, k], self.frac_turb_entr[i, k - 1])
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
                    abs(w_u) *
                    self.detr_sc[i, k] *
                    (updvar1 - envvar1) *
                    (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    rho0_half[k] *
                    self.UpdVar.Area.values[i, k] *
                    abs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                Covar.entr_gain[k] += dynamic_entr + turbulent_entr
                Covar.detr_loss[k] +=
                    tke_factor *
                    rho0_half[k] *
                    self.UpdVar.Area.values[i, k] *
                    abs(w_u) *
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
            Covar.detr_loss[k] += self.UpdVar.Area.values[i, k] * abs(w_u) * self.entr_sc[i, k]
        end
        Covar.detr_loss[k] *= rho0_half[k] * Covar.values[k]
    end
    return
end

function compute_covariance_rain(self::EDMF_PrognosticTKE, TS::TimeStepping, GMV::GridMeanVariables)
    # TODO defined again in compute_covariance_shear and compute_covaraince
    grid = get_grid(self)
    ae = 1 .- self.UpdVar.Area.bulkvalues # area of environment
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
    ae = 1 .- self.UpdVar.Area.bulkvalues
    rho0_half = reference_state(self).rho0_half

    @inbounds for k in real_center_indicies(grid)
        Covar.dissipation[k] = (
            rho0_half[k] * ae[k] * Covar.values[k] * max(self.EnvVar.TKE.values[k], 0)^0.5 /
            max(self.mixing_length[k], 1.0e-3) * self.tke_diss_coeff
        )
    end
    return
end


function compute_tke_advection(self::EDMF_PrognosticTKE)

    grid = get_grid(self)
    ref_state = reference_state(self)
    rho0_half = ref_state.rho0_half
    ae = 1 .- self.UpdVar.Area.bulkvalues # area of environment
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
    ae = 1 .- self.UpdVar.Area.bulkvalues # area of environment
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
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dzi = grid.dzi
    dti = TS.dti
    Ref = reference_state(self)
    alpha0LL = Ref.alpha0_half[kc_surf]
    zLL = grid.z_half[kc_surf]
    a = center_field(grid)
    b = center_field(grid)
    c = center_field(grid)
    x = center_field(grid)
    ae = 1 .- self.UpdVar.Area.bulkvalues
    ae_old = 1 .- up_sum(self.UpdVar.Area.old)
    rho_ae_K_m = face_field(grid)
    whalf = center_field(grid)
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
    wu_half = interp2pt(self.UpdVar.W.bulkvalues[kf_surf], self.UpdVar.W.bulkvalues[kf_surf + 1])

    # Not necessary if BCs for variances are applied to environment.
    # if GmvCovar.name=="tke"
    #     GmvCovar.values[kc_surf] =get_surface_tke(Case.Sur.ustar, self.wstar, get_grid(self).z_half[kc_surf], Case.Sur.obukhov_length)
    # elseif GmvCovar.name=="thetal_var"
    #     GmvCovar.values[kc_surf] = get_surface_variance(Case.Sur.rho_hflux * alpha0LL, Case.Sur.rho_hflux * alpha0LL, Case.Sur.ustar, zLL, Case.Sur.obukhov_length)
    # elseif GmvCovar.name=="qt_var"
    #     GmvCovar.values[kc_surf] = get_surface_variance(Case.Sur.rho_qtflux * alpha0LL, Case.Sur.rho_qtflux * alpha0LL, Case.Sur.ustar, zLL, Case.Sur.obukhov_length)
    # elseif GmvCovar.name=="thetal_qt_covar"
    #     GmvCovar.values[kc_surf] = get_surface_variance(Case.Sur.rho_hflux * alpha0LL, Case.Sur.rho_qtflux * alpha0LL, Case.Sur.ustar, zLL, Case.Sur.obukhov_length)
    # self.get_env_covar_from_GMV(self.UpdVar.Area, UpdVar1, UpdVar2, EnvVar1, EnvVar2, Covar, GmvVar1.values, GmvVar2.values, GmvCovar.values)

    Covar_surf = Covar.values[kc_surf]

    @inbounds for k in real_center_indicies(grid)
        D_env = 0.0

        @inbounds for i in xrange(self.n_updrafts)
            if self.UpdVar.Area.values[i, k] > self.minimum_area
                if Covar.name == "tke"
                    turb_entr =
                        interp2pt(self.frac_turb_entr[i, k - 1], self.frac_turb_entr[i, k]) / self.prandtl_nvec[k]
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

        a[k] = (-rho_ae_K_m[k - 1] * dzi * dzi)
        b[k] = (
            Ref.rho0_half[k] * ae[k] * dti - Ref.rho0_half[k] * ae[k] * whalf[k] * dzi +
            rho_ae_K_m[k] * dzi * dzi +
            rho_ae_K_m[k - 1] * dzi * dzi +
            D_env +
            Ref.rho0_half[k] * ae[k] * self.tke_diss_coeff * sqrt(max(self.EnvVar.TKE.values[k], 0)) /
            max(self.mixing_length[k], 1.0)
        )
        c[k] = (Ref.rho0_half[k + 1] * ae[k + 1] * whalf[k + 1] * dzi - rho_ae_K_m[k] * dzi * dzi)
        x[k] = (
            Ref.rho0_half[k] * ae_old[k] * Covar.values[k] * dti +
            Covar.press[k] +
            Covar.buoy[k] +
            Covar.shear[k] +
            Covar.entr_gain[k] +
            Covar.rain_src[k]
        ) #

        if is_surface_center(grid, k)
            a[k] = 0.0
            b[k] = 1.0
            c[k] = 0.0
            x[k] = Covar_surf
        end

        if is_toa_center(grid, k)
            b[k] += c[k]
            c[k] = 0.0
        end
    end
    # x .= tridiag_solve(x, a, b, c)
    cinterior = real_center_indicies(grid)
    x[cinterior] .= tridiag_solve(x[cinterior], a[cinterior], b[cinterior], c[cinterior])

    @inbounds for k in real_center_indicies(grid)
        if Covar.name == "thetal_qt_covar"
            Covar.values[k] = max(x[k], -sqrt(self.EnvVar.Hvar.values[k] * self.EnvVar.QTvar.values[k]))
            Covar.values[k] = min(x[k], sqrt(self.EnvVar.Hvar.values[k] * self.EnvVar.QTvar.values[k]))
        else
            Covar.values[k] = max(x[k], 0.0)
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
    ae = 1 .- self.UpdVar.Area.bulkvalues
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
                -self.horiz_K_eddy[i_last, k] * (self.EnvVar.W.values[k + 1] - self.EnvVar.W.values[k - 1]) /
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

        if is_surface_center(grid, k)
            Gmv_third_m.values[k] = 0.0 # this is here as first value is biased with BC area fraction
        else
            Gmv_third_m.values[k] =
                Upd_cubed + ae[k] * (env_mean.values[k]^3 + 3.0 * env_mean.values[k] * Envcov_) - GMVv_^3.0 -
                3.0 * GMVcov_ * GMVv_
        end
    end
    return
end


# Update the diagnosis of the inversion height, using the maximum temperature gradient method
function update_inversion(self::EDMF_PrognosticTKE, GMV::GridMeanVariables, option)
    theta_rho = center_field(self.Gr)
    ∇θ_liq_max = 0.0
    kc_surf = kc_surface(self.Gr)
    param_set = parameter_set(GMV)

    @inbounds for k in real_center_indicies(self.Gr)
        pp = TD.PhasePartition(GMV.QT.values[k], GMV.QL.values[k], 0.0)
        qv = TD.vapor_specific_humidity(pp)
        theta_rho[k] = theta_rho_c(self.Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], qv)
    end


    if option == "theta_rho"
        @inbounds for k in real_center_indicies(self.Gr)
            if theta_rho[k] > theta_rho[kc_surf]
                self.zi = self.Gr.z_half[k]
                break
            end
        end
    elseif option == "thetal_maxgrad"

        @inbounds for k in real_center_indicies(self.Gr)
            ∇θ_liq = ∇_upwind(GMV.H.values, self.Gr, k)
            if ∇θ_liq > ∇θ_liq_max
                ∇θ_liq_max = ∇θ_liq
                self.zi = self.Gr.z[k]
            end
        end
    elseif option == "critical_Ri"
        self.zi = get_inversion(param_set, theta_rho, GMV.U.values, GMV.V.values, self.Gr, Ri_bulk_crit(self))

    else
        error("INVERSION HEIGHT OPTION NOT RECOGNIZED")
    end

    return
end
