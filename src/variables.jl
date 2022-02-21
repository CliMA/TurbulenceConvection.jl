#####
##### Fields
#####

##### Auxiliary fields

# Center only
cent_aux_vars_en_2m(FT) = (;
    dissipation = FT(0),
    shear = FT(0),
    entr_gain = FT(0),
    detr_loss = FT(0),
    press = FT(0),
    buoy = FT(0),
    interdomain = FT(0),
    rain_src = FT(0),
)
cent_aux_vars_up(FT) = (;
    q_liq = FT(0),
    q_ice = FT(0),
    T = FT(0),
    RH = FT(0),
    s = FT(0),
    buoy = FT(0),
    area = FT(0),
    q_tot = FT(0),
    θ_liq_ice = FT(0),
    θ_liq_ice_tendency_precip_formation = FT(0),
    qt_tendency_precip_formation = FT(0),
    entr_sc = FT(0),
    detr_sc = FT(0),
    ε_nondim = FT(0),  # nondimensional entrainment
    δ_nondim = FT(0),  # nondimensional detrainment
    frac_turb_entr = FT(0),
    entr_turb_dyn = FT(0),
    detr_turb_dyn = FT(0),
    asp_ratio = FT(0),
    Π₁ = FT(0),
    Π₂ = FT(0),
    Π₃ = FT(0),
    Π₄ = FT(0),
)
cent_aux_vars_edmf(FT, edmf) = (;
    turbconv = (;
        ϕ_temporary = FT(0),
        ψ_temporary = FT(0),
        bulk = (;
            area = FT(0),
            θ_liq_ice = FT(0),
            RH = FT(0),
            buoy = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            T = FT(0),
            cloud_fraction = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            qt_tendency_precip_formation = FT(0),
        ),
        up = ntuple(i -> cent_aux_vars_up(FT), n_updrafts(edmf)),
        en = (;
            w = FT(0),
            area = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            θ_liq_ice = FT(0),
            θ_virt = FT(0),
            θ_dry = FT(0),
            RH = FT(0),
            s = FT(0),
            T = FT(0),
            buoy = FT(0),
            cloud_fraction = FT(0),
            tke = FT(0),
            Hvar = FT(0),
            QTvar = FT(0),
            HQTcov = FT(0),
            qt_tendency_precip_formation = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            unsat = (; q_tot = FT(0), θ_dry = FT(0), θ_virt = FT(0)),
            sat = (; T = FT(0), q_vap = FT(0), q_tot = FT(0), θ_dry = FT(0), θ_liq_ice = FT(0)),
            Hvar_rain_dt = FT(0),
            QTvar_rain_dt = FT(0),
            HQTcov_rain_dt = FT(0),
        ),
        θ_liq_ice_tendency_precip_sinks = FT(0),
        qt_tendency_precip_sinks = FT(0),
        qr_tendency_evap = FT(0),
        qs_tendency_melt = FT(0),
        qs_tendency_dep_sub = FT(0),
        qr_tendency_advection = FT(0),
        qs_tendency_advection = FT(0),
        en_2m = (;
            tke = cent_aux_vars_en_2m(FT),
            Hvar = cent_aux_vars_en_2m(FT),
            QTvar = cent_aux_vars_en_2m(FT),
            HQTcov = cent_aux_vars_en_2m(FT),
        ),
        term_vel_rain = FT(0),
        term_vel_snow = FT(0),
        KM = FT(0),
        KH = FT(0),
        mixing_length = FT(0),
        massflux_tendency_h = FT(0),
        massflux_tendency_qt = FT(0),
        diffusive_tendency_h = FT(0),
        diffusive_tendency_qt = FT(0),
        prandtl_nvec = FT(0),
        # Added by Ignacio : Length scheme in use (mls), and smooth min effect (ml_ratio)
        # Variable Prandtl number initialized as neutral value.
        mls = FT(0),
        b_exch = FT(0),
        ml_ratio = FT(0),
        w_up_c = FT(0),
        w_en_c = FT(0),
        Shear² = FT(0),
        ∂θv∂z = FT(0),
        ∂qt∂z = FT(0),
        ∂θl∂z = FT(0),
        ∂θv∂z_unsat = FT(0),
        ∂qt∂z_sat = FT(0),
        ∂θl∂z_sat = FT(0),
        l_entdet = FT(0),
        ϕ_gm = FT(0), # temporary for grid-mean variables
        ϕ_gm_cov = FT(0), # temporary for grid-mean covariance variables
        ϕ_en_cov = FT(0), # temporary for environmental covariance variables
        ϕ_up_cubed = FT(0), # temporary for cubed updraft variables in grid mean 3rd moment functions
    ),
)

# Face only
face_aux_vars_up(FT) = (;
    w = FT(0),
    nh_pressure = FT(0),
    nh_pressure_b = FT(0),
    nh_pressure_adv = FT(0),
    nh_pressure_drag = FT(0),
    massflux = FT(0),
)

face_aux_vars_edmf(FT, edmf) = (;
    turbconv = (;
        bulk = (; w = FT(0)),
        ρ_ae_KM = FT(0),
        ρ_ae_KH = FT(0),
        ρ_ae_K = FT(0),
        en = (; w = FT(0)),
        up = ntuple(i -> face_aux_vars_up(FT), n_updrafts(edmf)),
        massflux_h = FT(0),
        massflux_qt = FT(0),
        ϕ_temporary = FT(0),
        ψ_temporary = FT(0),
        diffusive_flux_h = FT(0),
        diffusive_flux_qt = FT(0),
        diffusive_flux_u = FT(0),
        diffusive_flux_v = FT(0),
    ),
)

##### Diagnostic fields

# Center only
cent_diagnostic_vars_edmf(FT, edmf) = (;
    turbconv = (;
        asp_ratio = FT(0),
        entr_sc = FT(0),
        ε_nondim = FT(0),
        detr_sc = FT(0),
        δ_nondim = FT(0),
        massflux = FT(0),
        frac_turb_entr = FT(0),
    ),
)

# Face only
face_diagnostic_vars_edmf(FT, edmf) =
    (; turbconv = (; nh_pressure = FT(0), nh_pressure_adv = FT(0), nh_pressure_drag = FT(0), nh_pressure_b = FT(0)))

##### Prognostic fields

# Center only
cent_prognostic_vars_up(FT, _) = (; ρarea = FT(0), ρaθ_liq_ice = FT(0), ρaq_tot = FT(0))
cent_prognostic_vars_up(FT, ::PrognosticNoisyRelaxationProcess) =
    (; ρarea = FT(0), ρaθ_liq_ice = FT(0), ρaq_tot = FT(0), ε_nondim = FT(0), δ_nondim = FT(0))
cent_prognostic_vars_en(FT) = (; ρatke = FT(0), ρaHvar = FT(0), ρaQTvar = FT(0), ρaHQTcov = FT(0))
cent_prognostic_vars_edmf(FT, edmf) = (;
    turbconv = (;
        en = cent_prognostic_vars_en(FT),
        up = ntuple(i -> cent_prognostic_vars_up(FT, edmf.entr_closure), n_updrafts(edmf)),
        pr = (; q_rai = FT(0), q_sno = FT(0)),
    ),
)
# cent_prognostic_vars_edmf(FT, edmf) = (;) # could also use this for empty model

# Face only
face_prognostic_vars_up(FT) = (; ρaw = FT(0))
face_prognostic_vars_edmf(FT, edmf) = (; turbconv = (; up = ntuple(i -> face_prognostic_vars_up(FT), n_updrafts(edmf))))
