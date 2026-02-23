# == get_qc_stats_from_moments diagnostics (use downgradient diffusive fluxes) ==
    # aux_shoc = aux_tc.shoc
    use_prognostic_qc = edmf.moisture_model isa NonEquilibriumMoisture
    w_qt_cov_c = Ic.(aux_tc_f.diffusive_flux_qt)
    w_h_cov_c = Ic.(aux_tc_f.diffusive_flux_h)
    @inbounds for k in real_center_indices(grid)
        ts_k = aux_en.ts[k]
        qt_var_k = aux_en.QTvar[k]
        h_var_k = aux_en.Hvar[k]
        h_qt_cov_k = aux_en.HQTcov[k]

        w_qt_cov_k = w_qt_cov_c[k]
        w_h_cov_k = w_h_cov_c[k]

        ql_true_k = use_prognostic_qc ? aux_en.q_liq[k] : FT(NaN)
        qi_true_k = use_prognostic_qc ? aux_en.q_ice[k] : FT(NaN)

        ql_k, qi_k, cov_ql_qt_k, cov_ql_h_k, w_ql_k, w_qi_k, ql_qt_sd_k, qi_qt_sd_k, qc_mean_pdf_k, w_qc_cov_pdf_k, ql_pdf_k, qi_pdf_k, Φ_k, α_k, f_liq_k, f_ice_k, scale_l_k, scale_i_k, w_ql_pot_k, w_qi_pot_k =
            get_qc_stats_from_moments(
                param_set,
                ts_k,
                qt_var_k,
                h_var_k,
                h_qt_cov_k,
                w_qt_cov_k,
                w_h_cov_k,
                ql_true_k,
                qi_true_k;
            )

        # aux_shoc.ql_from_moments[k] = ql_k
        # aux_shoc.qi_from_moments[k] = qi_k
        # aux_shoc.cov_ql_qt_from_moments[k] = cov_ql_qt_k
        # aux_shoc.cov_ql_h_from_moments[k] = cov_ql_h_k
        # aux_shoc.w_ql_from_moments[k] = w_ql_k
        # aux_shoc.w_qi_from_moments[k] = w_qi_k
        # aux_shoc.ql_qt_sd[k] = ql_qt_sd_k
        # aux_shoc.qi_qt_sd[k] = qi_qt_sd_k
        # aux_shoc.qc_mean_pdf[k] = qc_mean_pdf_k
        # aux_shoc.w_qc_cov_pdf[k] = w_qc_cov_pdf_k
        # aux_shoc.ql_pdf[k] = ql_pdf_k
        # aux_shoc.qi_pdf[k] = qi_pdf_k
        # aux_shoc.pdf_cloud_frac[k] = Φ_k
        # aux_shoc.pdf_alpha[k] = α_k
        # aux_shoc.phase_frac_liq[k] = f_liq_k
        # aux_shoc.phase_frac_ice[k] = f_ice_k
        # aux_shoc.flux_scale_liq[k] = scale_l_k
        # aux_shoc.flux_scale_ice[k] = scale_i_k
        # aux_shoc.w_ql_pot[k] = w_ql_pot_k
        # aux_shoc.w_qi_pot[k] = w_qi_pot_k
    end
    # === end get_qc_stats_from_moments == #

    # == cSigma correlation/covariance matrix diagnostics ==
    # Interpolate face fluxes to centers for use in cSigma
    diffusive_flux_ql_c = Ic.(aux_tc_f.diffusive_flux_ql)
    diffusive_flux_qr_c = Ic.(aux_tc_f.diffusive_flux_qr)
    diffusive_flux_qi_c = Ic.(aux_tc_f.diffusive_flux_qi)
    diffusive_flux_qs_c = Ic.(aux_tc_f.diffusive_flux_qs)

    @inbounds for k in real_center_indices(grid)
        # Extract hydrometeor means from environment
        mean_ql = aux_en.q_liq[k]
        mean_qr = center_prog_precipitation(state).q_rai[k]
        mean_qi = aux_en.q_ice[k]
        mean_qs = center_prog_precipitation(state).q_sno[k]
        mean_qt = aux_en.q_tot[k]  # Total water in environment
        mean_h = aux_en.θ_liq_ice[k]  # Liquid-ice potential temperature in environment

        # Extract known variances
        var_qt = aux_en.QTvar[k]
        var_h = aux_en.Hvar[k]

        # Get TKE from environment
        tke_c = aux_en.tke[k]

        # Get turbulent fluxes (already interpolated to centers above)
        flux_ql = diffusive_flux_ql_c[k]
        flux_qr = diffusive_flux_qr_c[k]
        flux_qi = diffusive_flux_qi_c[k]
        flux_qs = diffusive_flux_qs_c[k]
        flux_qt = w_qt_cov_c[k]
        flux_h = w_h_cov_c[k]

        # Only compute if we have meaningful data
        if !any(isnan.([mean_ql, mean_qr, mean_qi, mean_qs, mean_qt, mean_h, tke_c]))
            Cov = hydrometeor_covariances(
                FT(mean_ql), FT(mean_qr), FT(mean_qi), FT(mean_qs), FT(mean_qt), FT(mean_h), FT(tke_c);
                var_qt=FT(var_qt), var_h=FT(var_h),
                flux_ql=FT(flux_ql), flux_qr=FT(flux_qr), flux_qi=FT(flux_qi), flux_qs=FT(flux_qs), flux_qt=FT(flux_qt), flux_h=FT(flux_h)
            )

            # Debug: Check why Cov might be NaN
            if any(isnan.(Cov))
                @warn "NaN in Cov at k=$k: mean_qt=$(mean_qt), var_qt=$(var_qt), mean_h=$(mean_h), var_h=$(var_h)"
            end

            # Only extract correlations if covariance matrix is valid
            if !any(isnan.(Cov))
                # Extract correlations and covariances using cSigma utility
                corr_nt = hydrometeor_correlations_from_covariance(Cov)

                # # Store correlations
                # aux_shoc.csigma_corr_w_ql[k] = corr_nt.corr_w_ql
                # aux_shoc.csigma_corr_w_qr[k] = corr_nt.corr_w_qr
                # aux_shoc.csigma_corr_w_qi[k] = corr_nt.corr_w_qi
                # aux_shoc.csigma_corr_w_qs[k] = corr_nt.corr_w_qs
                # aux_shoc.csigma_corr_w_qt[k] = corr_nt.corr_w_qt
                # aux_shoc.csigma_corr_w_h[k] = corr_nt.corr_w_h
                # aux_shoc.csigma_corr_ql_qr[k] = corr_nt.corr_ql_qr
                # aux_shoc.csigma_corr_ql_qi[k] = corr_nt.corr_ql_qi
                # aux_shoc.csigma_corr_ql_qs[k] = corr_nt.corr_ql_qs
                # aux_shoc.csigma_corr_ql_qt[k] = corr_nt.corr_ql_qt
                # aux_shoc.csigma_corr_ql_h[k] = corr_nt.corr_ql_h
                # aux_shoc.csigma_corr_qr_qi[k] = corr_nt.corr_qr_qi
                # aux_shoc.csigma_corr_qr_qs[k] = corr_nt.corr_qr_qs
                # aux_shoc.csigma_corr_qr_qt[k] = corr_nt.corr_qr_qt
                # aux_shoc.csigma_corr_qr_h[k] = corr_nt.corr_qr_h
                # aux_shoc.csigma_corr_qi_qs[k] = corr_nt.corr_qi_qs
                # aux_shoc.csigma_corr_qi_qt[k] = corr_nt.corr_qi_qt
                # aux_shoc.csigma_corr_qi_h[k] = corr_nt.corr_qi_h
                # aux_shoc.csigma_corr_qs_qt[k] = corr_nt.corr_qs_qt
                # aux_shoc.csigma_corr_qs_h[k] = corr_nt.corr_qs_h
                # aux_shoc.csigma_corr_qt_h[k] = corr_nt.corr_qt_h

                # # Store covariances (diagonal and unique off-diagonal elements)
                # aux_shoc.csigma_cov_w[k] = Cov[1, 1]
                # aux_shoc.csigma_cov_ql[k] = Cov[2, 2]
                # aux_shoc.csigma_cov_qr[k] = Cov[3, 3]
                # aux_shoc.csigma_cov_qi[k] = Cov[4, 4]
                # aux_shoc.csigma_cov_qs[k] = Cov[5, 5]
                # aux_shoc.csigma_cov_qt[k] = Cov[6, 6]
                # aux_shoc.csigma_cov_h[k] = Cov[7, 7]
                # aux_shoc.csigma_cov_w_ql[k] = Cov[1, 2]
                # aux_shoc.csigma_cov_w_qr[k] = Cov[1, 3]
                # aux_shoc.csigma_cov_w_qi[k] = Cov[1, 4]
                # aux_shoc.csigma_cov_w_qs[k] = Cov[1, 5]
                # aux_shoc.csigma_cov_w_qt[k] = Cov[1, 6]
                # aux_shoc.csigma_cov_w_h[k] = Cov[1, 7]
                # aux_shoc.csigma_cov_ql_qr[k] = Cov[2, 3]
                # aux_shoc.csigma_cov_ql_qi[k] = Cov[2, 4]
                # aux_shoc.csigma_cov_ql_qs[k] = Cov[2, 5]
                # aux_shoc.csigma_cov_ql_qt[k] = Cov[2, 6]
                # aux_shoc.csigma_cov_ql_h[k] = Cov[2, 7]
                # aux_shoc.csigma_cov_qr_qi[k] = Cov[3, 4]
                # aux_shoc.csigma_cov_qr_qs[k] = Cov[3, 5]
                # aux_shoc.csigma_cov_qr_qt[k] = Cov[3, 6]
                # aux_shoc.csigma_cov_qr_h[k] = Cov[3, 7]
                # aux_shoc.csigma_cov_qi_qs[k] = Cov[4, 5]
                # aux_shoc.csigma_cov_qi_qt[k] = Cov[4, 6]
                # aux_shoc.csigma_cov_qi_h[k] = Cov[4, 7]
                # aux_shoc.csigma_cov_qs_qt[k] = Cov[5, 6]
                # aux_shoc.csigma_cov_qs_h[k] = Cov[5, 7]
                # aux_shoc.csigma_cov_qt_h[k] = Cov[6, 7]
            end
        end
    end
    # === end cSigma correlation/covariance output ==