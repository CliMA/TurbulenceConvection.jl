#####
##### Fields
#####
import ClimaCore.Geometry: ⊗

# Helpers for adding empty thermodynamic state fields:
thermo_state(FT, ::EquilibriumMoisture) = TD.PhaseEquil{FT}(0, 0, 0, 0, 0)
thermo_state(FT, ::NonEquilibriumMoisture) = TD.PhaseNonEquil{FT}(0, 0, TD.PhasePartition(FT(0), FT(0), FT(0)))

##### Auxiliary fields

# note, supersat variables are currently not stored in the grid mean though we do store cloud_sedimentation variables there.
supersat_variables(FT, ::EquilibriumMoisture) = NamedTuple()
supersat_variables(FT, ::NonEquilibriumMoisture) = (;
    q_vap_sat_liq = FT(0),
    q_vap_sat_ice = FT(0),
    τ_ice = FT(0), # time scale for ice supersaturation
    τ_liq = FT(0), # time scale for liquid supersaturation
    # N_i = FT(0), # number of ice crystals  [ moved to cloud sedimentation variables so that they exist even in EquilibriumMoisture ]
    # N_l = FT(0), # number of liquid droplets [ moved to cloud sedimentation variables so that they exist even in EquilibriumMoisture ]
)

cloud_sedimentation_variables(FT, ::CloudSedimentationModel) = (; # these are different from the rain snow ones bc they are not on the grid mean, but also in env/up separately
    term_vel_liq = FT(0), # rn we're not using this
    term_vel_ice = FT(0),
    N_i = FT(0), # should these be in supersat variables? no bc we want them stored no matter what.
    N_l = FT(0), # idk if we want this yet..
    r_i_mean = FT(0), # mean ice radius
    r_l_mean = FT(0), # mean liquid radius
    #
    dN_i_dz = FT(0), # vertical gradient of N_i, used to adjust the timescale for acnv
    # dN_l_dz = FT(0), # vertical gradient of N_l, used to adjust the timescale for acnv
    # dr_i_mean_dz = FT(0), # vertical gradient of r_i_mean, used to adjust the timescale for acnv
    # dr_l_mean_dz = FT(0), # vertical gradient of r_l_mean, used to adjust the timescale for acnv
    dqidz = FT(0), # vertical gradient of qi, used to adjust the timescale for acnv
    #
    N_i_no_boost = FT(0), # N_i without the massflux boost factor
)
cloud_sedimentation_variables(FT, ::CloudNoSedimentationModel) =  (;
    N_i = FT(0), # should these be in supersat variables? no bc we want them stored no matter what.
    N_l = FT(0), # idk if we want this yet..
    r_i_mean = FT(0), # mean ice radius
    r_l_mean = FT(0), # mean liquid radius
    #
    dN_i_dz = FT(0), # vertical gradient of N_i, used to adjust the timescale for acnv
    # dN_l_dz = FT(0), # vertical gradient of N_l, used to adjust the timescale for acnv
    # dr_i_mean_dz = FT(0), # vertical gradient of r_i_mean, used to adjust the timescale for acnv
    # dr_l_mean_dz = FT(0), # vertical gradient of r_l_mean, used to adjust the timescale for acnv
    dqidz = FT(0), # vertical gradient of qi, used to adjust the timescale for acnv
    #
    N_i_no_boost = FT(0), # N_i without the massflux boost factor
)
cloud_sedimentation_variables(FT, ::Nothing) =  NamedTuple()

# consolidate my additions so it's easier to edit. These will exist everywhere... (though they may not be used or defined in eq case...)
my_microphysics_additions(FT, moisture_model::AbstractMoistureModel, cloud_sedimentation_model::Union{AbstractCloudSedimentationModel, Nothing}) = (;
    #
    ql_tendency_cond_evap = FT(0), # tendency due to condensation/evaporation
    qi_tendency_sub_dep = FT(0), # tendency due to sublimation/deposition
    #
    ql_tendency_sedimentation = FT(0), # tendency due to sedimentation
    qi_tendency_sedimentation = FT(0), # tendency due to sedimentation
    qt_tendency_sedimentation = FT(0), # tendency due to sedimentation
    θ_liq_ice_tendency_sedimentation = FT(0),
    #
    ql_tendency_sedimentation_other = FT(0), # this is the sedimentation contribution to precip formation from the other updrafts
    qi_tendency_sedimentation_other = FT(0), # this is the sedimentation contribution to precip formation from the other updrafts
    qt_tendency_sedimentation_other = FT(0), # this is the sedimentation contribution to precip formation from the other updrafts
    θ_liq_ice_tendency_sedimentation_other = FT(0), # this is the sedimentation contribution to precip formation from the other updrafts
    #
    # # # diagnostic stuff (not used in calculations)
    #
    ql_tendency_acnv = FT(0), # tendency due to autoconversion
    qi_tendency_acnv = FT(0), # tendency due to autoconversion
    qi_tendency_acnv_dep = FT(0), # tendency due to sublimation/deposition
    qi_tendency_acnv_dep_is = FT(0), # tendency due to sublimation/deposition, only from crossing r_is.
    qi_tendency_acnv_dep_above = FT(0), # tendency due to sublimation/deposition, only from particles that were already above r_is
    qi_tendency_acnv_agg = FT(0), # tendency due to autoconversion aggregation
    qi_tendency_acnv_agg_other = FT(0), # tendency due to autoconversion aggregation from other updrafts
    qi_tendency_acnv_agg_mix = FT(0), # tendency due to autoconversion aggregation from other updrafts that goes into snow (i.e. not the one that goes into rain)
    qi_tendency_acnv_thresh = FT(0), # tendency due to autoconversion aggregation from other updrafts that goes into snow (i.e. not the one that goes into rain) but is thresholded by the microphysics
    #
    ql_tendency_accr_liq_rai = FT(0), # tendency due to accretion of liquid by rain
    ql_tendency_accr_liq_ice = FT(0), # tendency due to accretion of liquid by ice
    ql_tendency_accr_liq_sno = FT(0), # tendency due to accretion of liquid by snow
    #
    qi_tendency_accr_ice_liq = FT(0), # tendency due to accretion of ice by liquid
    qi_tendency_accr_ice_rai = FT(0), # tendency due to accretion of ice by rain
    qi_tendency_accr_ice_sno = FT(0), # tendency due to accretion of ice by snow
    # (liq -> ice)
    qi_tendency_hom_frz = FT(0), # tendency due to homogeneous freezing
    qi_tendency_het_frz = FT(0), # tendency due to heterogeneous freezing
    qi_tendency_het_nuc = FT(0), # tendency due to heterogeneous nucleation
    qi_tendency_melt = FT(0), # tendency due to melting
    #
    dqvdt = FT(0), # store the tendency of water vapor for diagnostic purposes
    dTdt = FT(0), # store the tendency of temperature for diagnostic purposes
    #
    dTdz = FT(0), # store the vertical gradient of temperature for diagnostic purposes. Store here instead of aux_tc because we maybe want separately for env and up, though we don't do boost in the updraft so idk.
    #
    # supersat
    supersat_variables(FT, moisture_model)...,
    #
    # cloud sedimentation
    cloud_sedimentation_variables(FT, cloud_sedimentation_model)...,
)

my_microphysics_additions(FT, mm::NonEquilibriumMoisture) = my_microphysics_additions(FT, mm, nothing)
my_microphysics_additions(FT, mm::EquilibriumMoisture) = my_microphysics_additions(FT, mm, nothing)

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
cent_aux_vars_up_moisture(FT, moisture_model::NonEquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel) = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model)...,
)
cent_aux_vars_up_moisture(FT, moisture_model::EquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel) = (;
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model)...,
)
cent_aux_vars_up(FT, local_geometry, edmf) = (;
    ts = thermo_state(FT, edmf.moisture_model),
    q_liq = FT(0),
    q_ice = FT(0),
    T = FT(0),
    p = FT(0), # store this instead of constantly recalculating it [ both phasenonequil and phaseequil store ρ but phaseequil doesn't store p ] -- note typically in our case it should be the same as p_c but better safe than sorry.
    RH = FT(0),
    RH_ice = FT(0), # RH for ice, i.e. the ratio of q_vap to q_vap_sat_ice
    s = FT(0),
    buoy = FT(0),
    area = FT(0),
    q_tot = FT(0),
    θ_liq_ice = FT(0),
    θ_liq_ice_tendency_precip_formation = FT(0),
    qt_tendency_precip_formation = FT(0),
    cent_aux_vars_up_moisture(FT, edmf.moisture_model, edmf.cloud_sedimentation_model)...,
    entr_sc = FT(0),
    detr_sc = FT(0),
    entr_ml = FT(0),
    detr_ml = FT(0),
    ε_nondim = FT(0),  # nondimensional entrainment
    δ_nondim = FT(0),  # nondimensional detrainment
    ε_ml_nondim = FT(0),  # nondimensional entrainment
    δ_ml_nondim = FT(0),  # nondimensional detrainment
    frac_turb_entr = FT(0),
    entr_turb_dyn = FT(0),
    detr_turb_dyn = FT(0),
    entr_rate_inv_s = FT(0),
    detr_rate_inv_s = FT(0),
    asp_ratio = FT(0),
    Π_groups = ntuple(i -> FT(0), n_Π_groups(edmf)),
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::NonEquilibriumMoisture) = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model)..., # don't need sedimentation in bulk
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::NonEquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel) = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model)..., # added sedimentation to bulk for storing N_i etc...
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::EquilibriumMoisture) = (;
    my_microphysics_additions(FT, moisture_model)..., # don't need sedimentation in bulk
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::EquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel) = (;
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model)..., # added sedimentation to bulk for storing N_i etc...
)
cent_aux_vars_edmf_en_moisture(FT, moisture_model::NonEquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel) = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model)...,
)
cent_aux_vars_edmf_en_moisture(FT, moisture_model::EquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel) = (;
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model)...,
)
cent_aux_vars_edmf_moisture(FT, ::NonEquilibriumMoisture,) = (;
    # massflux_tendency_ql = FT(0), # moved so is always there even for eq so can put in diagnostics output
    # massflux_tendency_qi = FT(0), # moved so is always there even for eq so can put in diagnostics output
    # diffusive_tendency_ql = FT(0), # moved so is always there even for eq so can put in diagnostics output
    # diffusive_tendency_qi = FT(0), # moved so is always there even for eq so can put in diagnostics output
)
cent_aux_vars_edmf_moisture(FT, ::EquilibriumMoisture) = NamedTuple()
cent_aux_vars_edmf(::Type{FT}, local_geometry, edmf) where {FT} = (;
    turbconv = (;
        ϕ_temporary = FT(0),
        ψ_temporary = FT(0),
        φ_temporary = FT(0),
        k̂ = CCG.Contravariant3Vector(CCG.WVector(FT(1)), local_geometry),
        bulk = (;
            area = FT(0),
            θ_liq_ice = FT(0),
            RH = FT(0),
            RH_ice = FT(0),
            buoy = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            T = FT(0),
            p = FT(0), # store this instead of constantly recalculating it [ both phasenonequil and phaseequil store ρ but phaseequil doesn't store p ] -- note typically in our case it should be the same as p_c but better safe than sorry.
            cloud_fraction = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            qt_tendency_precip_formation = FT(0),
            # cent_aux_vars_edmf_bulk_moisture(FT, edmf.moisture_model)...,
            cent_aux_vars_edmf_bulk_moisture(FT, edmf.moisture_model, edmf.cloud_sedimentation_model)..., # add sedimentation to bulk for storing N_i etc...
            # tendencies = cent_prognostic_vars_up(FT, edmf), # storage for tendencies limiter so we can reuse same memory
            # tendencies_adjustments = cent_prognostic_vars_up(FT, edmf), # storage for tendencies limiter so we can reuse same memory
            # prognostic = cent_prognostic_vars_up(FT, edmf), # storage for prognostic bulk w/ no clippings so we can be sure we have the true bulk for limiting
        ),
        up = ntuple(i -> cent_aux_vars_up(FT, local_geometry, edmf), Val(n_updrafts(edmf))),
        en = (;
            ts = thermo_state(FT, edmf.moisture_model),
            w = FT(0),
            area = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            θ_liq_ice = FT(0),
            θ_virt = FT(0),
            θ_dry = FT(0),
            RH = FT(0),
            RH_ice = FT(0),
            s = FT(0),
            T = FT(0),
            p = FT(0), # store this instead of constantly recalculating it [ both phasenonequil and phaseequil store ρ but phaseequil doesn't store p ] -- note typically in our case it should be the same as p_c but better safe than sorry.
            buoy = FT(0),
            cloud_fraction = FT(0),
            tke = FT(0),
            Hvar = FT(0),
            QTvar = FT(0),
            HQTcov = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            qt_tendency_precip_formation = FT(0),
            cent_aux_vars_edmf_en_moisture(FT, edmf.moisture_model, edmf.cloud_sedimentation_model)...,
            unsat = (; q_tot = FT(0), θ_dry = FT(0), θ_virt = FT(0)),
            sat = (; T = FT(0), q_vap = FT(0), q_tot = FT(0), θ_dry = FT(0), θ_liq_ice = FT(0)),
            Hvar_rain_dt = FT(0),
            QTvar_rain_dt = FT(0),
            HQTcov_rain_dt = FT(0),
            # prognostic = cent_prognostic_vars_up(FT, edmf), # storage for tracking what en is w/o any clippings
        ),
        θ_liq_ice_tendency_precip_sinks = FT(0),
        qt_tendency_precip_sinks = FT(0),
        qr_tendency_evap = FT(0),
        qs_tendency_melt = FT(0),
        qs_tendency_dep_sub = FT(0),
        qr_tendency_advection = FT(0),
        qs_tendency_advection = FT(0),
        qs_tendency_accr_rai_sno = FT(0), # tendency due to accretion QR by QS [QR -> QS]  (we store it this way, it's always to snow below freezing and [QS -> QR] above freezing.)
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
        KQ = FT(0),
        mixing_length = FT(0),
        massflux = FT(0),
        massflux_tendency_h = FT(0),
        massflux_tendency_qt = FT(0),
        massflux_tendency_ql = FT(0), # moving to here so they're always there even for eq so can put in diagnostics output
        massflux_tendency_qi = FT(0), #  moving to here so they're always there even for eq so can put in diagnostics output
        diffusive_tendency_h = FT(0),
        diffusive_tendency_qt = FT(0),
        diffusive_tendency_ql = FT(0), # my addition [ this is the diffusive tendency of ql, not including  massflux tendency ]
        diffusive_tendency_qi = FT(0), # my addition [ this is the
        diffusive_tendency_qr = FT(0), # my addition [ this is the diffusive tendency of qr, not including  massflux tendency ]
        diffusive_tendency_qs = FT(0), # my addition [ this is the diffusive tendency of qs, not including  massflux tendency ]
        cent_aux_vars_edmf_moisture(FT, edmf.moisture_model)...,
        prandtl_nvec = FT(0),
        # Added by Ignacio : Length scheme in use (mls), and smooth min effect (ml_ratio)
        # Variable Prandtl number initialized as neutral value.
        mls = FT(0),
        b_exch = FT(0),
        ml_ratio = FT(0),
        w_up_c = FT(0),
        w_en_c = FT(0),
        Shear² = FT(0),
        ∂M∂z = FT(0),
        ∂w∂z = FT(0),
        ∂M∂z_ρa = FT(0),
        ∂lnM∂z = FT(0),
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
    )
)

# Face only
face_aux_vars_up(FT, local_geometry) = (;
    w = FT(0),
    nh_pressure = FT(0),
    nh_pressure_b = FT(0),
    nh_pressure_adv = FT(0),
    nh_pressure_drag = FT(0),
    massflux = FT(0),
    area = FT(0),
)
face_aux_vars_edmf_moisture(FT, ::NonEquilibriumMoisture) = (;
    massflux_en = FT(0), # TODO: is this the right place for this?
    massflux_ql = FT(0),
    massflux_qi = FT(0),
    diffusive_flux_ql = FT(0),
    diffusive_flux_qi = FT(0),
)
face_aux_vars_edmf_moisture(FT, ::EquilibriumMoisture) = NamedTuple()
face_aux_vars_edmf(::Type{FT}, local_geometry, edmf) where {FT} = (;
    turbconv = (;
        bulk = (; w = FT(0), a_up = FT(0), 
        # tendencies = face_prognostic_vars_up(FT, local_geometry), # storage for tendencies limiter so we can reuse same memory
        # tendencies_adjustments = face_prognostic_vars_up(FT, local_geometry), # storage for tendencies limiter so we can reuse same memory
        # prognostic = face_prognostic_vars_up(FT, local_geometry), # storage for prognostic bulk w/ no clippings so we can be sure we have the true bulk for limiting. This can differ from aux_bulk bc of clippings a la minimum_area etc...
        ),
        ρ_ae_KM = FT(0),
        ρ_ae_KH = FT(0),
        ρ_ae_KQ = FT(0),
        ρ_ae_K = FT(0),
        en = (;
            w = FT(0),
            # prognostic = face_prognostic_vars_up(FT, local_geometry), # storage for tracking what en is w/o any clippings (prog_gm - prog_bulk) [deprecated]
        ),
        up = ntuple(i -> face_aux_vars_up(FT, local_geometry), Val(n_updrafts(edmf))),
        massflux = FT(0),
        massflux_h = FT(0),
        massflux_qt = FT(0),
        ϕ_temporary = FT(0),
        diffusive_flux_h = FT(0),
        diffusive_flux_qt = FT(0),
        diffusive_flux_qr = FT(0), # diffusive flux of rain [ my addition ]
        diffusive_flux_qs = FT(0), # diffusive flux of snow [ my addition ]
        face_aux_vars_edmf_moisture(FT, edmf.moisture_model)...,
        diffusive_flux_uₕ = CCG.Covariant3Vector(FT(0)) ⊗ CCG.Covariant12Vector(FT(0), FT(0)),
        uvw = CCG.Covariant123Vector(CCG.WVector(FT(0)), local_geometry),
    )
)

##### Diagnostic fields

# Center only
cent_diagnostic_vars_edmf(FT, local_geometry, edmf) = (;
    turbconv = (;
        asp_ratio = FT(0),
        entr_sc = FT(0),
        ε_nondim = FT(0),
        detr_sc = FT(0),
        δ_nondim = FT(0),
        entr_ml = FT(0),
        ε_ml_nondim = FT(0),
        detr_ml = FT(0),
        δ_ml_nondim = FT(0),
        massflux = FT(0),
        frac_turb_entr = FT(0),
        entr_rate_inv_s = FT(0),
        detr_rate_inv_s = FT(0),
        Π_groups = ntuple(i -> FT(0), n_Π_groups(edmf)),
        Π₁ = FT(0),
        Π₂ = FT(0),
        Π₃ = FT(0),
        Π₄ = FT(0),
        Π₅ = FT(0),
        Π₆ = FT(0),
    )
)

# Face only
face_diagnostic_vars_edmf(FT, local_geometry, edmf) = (;
    turbconv = (; nh_pressure = FT(0), nh_pressure_adv = FT(0), nh_pressure_drag = FT(0), nh_pressure_b = FT(0)),
    precip = (; rain_flux = FT(0), snow_flux = FT(0)),
)

# Single value per column diagnostic variables
single_value_per_col_diagnostic_vars_edmf(FT, edmf) = (;
    turbconv = (;
        env_cloud_base = FT(0),
        env_cloud_top = FT(0),
        env_cloud_cover = FT(0),
        env_lwp = FT(0),
        env_iwp = FT(0),
        updraft_cloud_cover = FT(0),
        updraft_cloud_base = FT(0),
        updraft_cloud_top = FT(0),
        updraft_lwp = FT(0),
        updraft_iwp = FT(0),
        Hd = FT(0),
    )
)

##### Prognostic fields

# Center only
cent_prognostic_vars_up_noisy_relaxation(::Type{FT}, ::PrognosticNoisyRelaxationProcess) where {FT} =
    (; ε_nondim = FT(0), δ_nondim = FT(0))
cent_prognostic_vars_up_noisy_relaxation(::Type{FT}, _) where {FT} = NamedTuple()
cent_prognostic_vars_up_moisture(::Type{FT}, ::EquilibriumMoisture) where {FT} = NamedTuple()
cent_prognostic_vars_up_moisture(::Type{FT}, ::NonEquilibriumMoisture) where {FT} = (; ρaq_liq = FT(0), ρaq_ice = FT(0))
cent_prognostic_vars_up(::Type{FT}, edmf) where {FT} = (;
    ρarea = FT(0),
    ρaθ_liq_ice = FT(0),
    ρaq_tot = FT(0),
    cent_prognostic_vars_up_noisy_relaxation(FT, edmf.entr_closure)...,
    cent_prognostic_vars_up_moisture(FT, edmf.moisture_model)...,
)

cent_prognostic_vars_en(::Type{FT}, edmf) where {FT} =
    (; ρatke = FT(0), cent_prognostic_vars_en_thermo(FT, edmf.thermo_covariance_model)...)
cent_prognostic_vars_en_thermo(::Type{FT}, ::DiagnosticThermoCovariances) where {FT} = NamedTuple()
cent_prognostic_vars_en_thermo(::Type{FT}, ::PrognosticThermoCovariances) where {FT} =
    (; ρaHvar = FT(0), ρaQTvar = FT(0), ρaHQTcov = FT(0))
cent_prognostic_vars_edmf(::Type{FT}, edmf) where {FT} = (;
    turbconv = (;
        en = cent_prognostic_vars_en(FT, edmf),
        up = ntuple(i -> cent_prognostic_vars_up(FT, edmf), Val(n_updrafts(edmf))),
        pr = (; q_rai = FT(0), q_sno = FT(0)),
    )
)
# cent_prognostic_vars_edmf(FT, edmf) = (;) # could also use this for empty model

# Face only
face_prognostic_vars_up(::Type{FT}, local_geometry) where {FT} = (; ρaw = FT(0))
face_prognostic_vars_edmf(::Type{FT}, local_geometry, edmf) where {FT} =
    (; turbconv = (; up = ntuple(i -> face_prognostic_vars_up(FT, local_geometry), Val(n_updrafts(edmf)))))
