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
)



cloud_sedimentation_variables(FT, ::CloudNoSedimentationModel) =  (;
    N_i = FT(0), # [[ These are here instead of in supersat variables so that they exist even in EquilibriumMoisture ]]
    N_l = FT(0), # [[ These are here instead of in supersat variables so that they exist even in EquilibriumMoisture ]]
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
    #
)
cloud_sedimentation_variables(FT, cloud_sedimentation_model::CloudSedimentationModel) = (; # these are different from the rain snow ones bc they are not on the grid mean, but also in env/up separately
    term_vel_liq = FT(0), # rn we're not using this
    term_vel_ice = FT(0),
    cloud_sedimentation_variables(FT, CloudNoSedimentationModel(cloud_sedimentation_model))...,
)
cloud_sedimentation_variables(FT, ::Nothing) =  NamedTuple()

cloak_variables(FT, ::Nothing, ::Val{calibrate_io}) where {calibrate_io} =  NamedTuple()
cloak_variables(FT, ::Nothing, ::AbstractMoistureModel, ::Val{calibrate_io}) where {calibrate_io} =  NamedTuple()
cloak_variables(FT, ::StandardAreaPartitionModel, ::Val{calibrate_io}) where {calibrate_io} = NamedTuple()
cloak_variables(FT, ::StandardAreaPartitionModel, ::AbstractMoistureModel, ::Val{calibrate_io}) where {calibrate_io} = NamedTuple()
cloak_variables(FT, ::CoreCloakAreaPartitionModel, moisture_model::AbstractMoistureModel, ::Val{calibrate_io}) where {calibrate_io} = (;
    a_cloak_up = FT(0), # env area in updraft cloak
    a_cloak_dn = FT(0), # env area in downdraft cloak
    a_en_remaining = FT(0), # env area not in cloak
    #
    #
    q_tot_cloak_up = FT(0), # env total water in updraft cloak
    q_tot_cloak_dn = FT(0), # env total water in downdraft cloak
    # q_tot_en_remaining = FT(0), # env total water not in cloak [[ we have the up and down cloak balance to not change this ]]
    #
    θ_liq_ice_cloak_up = FT(0), # env liquid-ice potential temperature in updraft cloak
    θ_liq_ice_cloak_dn = FT(0), # env liquid-ice potential temperature in downdraft cloak
    # θ_liq_ice_en_remaining = FT(0), # env liquid-ice potential temperature not in cloak [[ we have the up and down cloak balance to not change this ]]
    #
    # Liq and Ice are stored even in EquilibriumMoisture for other microphysics calculations
    #
    q_liq_cloak_up = FT(0), # env liquid water in updraft cloak
    q_liq_cloak_dn = FT(0), # env liquid water in downdraft cloak
    #
    q_ice_cloak_up = FT(0), # env ice water in updraft cloak
    q_ice_cloak_dn = FT(0), # env ice water in downdraft cloak
    #
    # TBD if we need ts, T, or other cloak state variables stored for efficiency.
    T_cloak_up = FT(0), # env temperature in updraft cloak
    T_cloak_dn = FT(0), # env temperature in downdraft cloak
    # thermodynamic state variables
    ts_cloak_up = thermo_state(FT, moisture_model),
    ts_cloak_dn = thermo_state(FT, moisture_model),

    # temporary
    (calibrate_io ? (; ) : (;
        RH_liq_cloak_up = FT(0), # env liquid relative humidity in updraft cloak
        RH_liq_cloak_dn = FT(0), # env liquid relative humidity in downdraft cloak
        RH_ice_cloak_up = FT(0), # env ice relative humidity in updraft cloak
        RH_ice_cloak_dn = FT(0), # env ice relative humidity in downdraft cloak
        RH_liq_en_remaining = FT(0), # env liquid relative humidity not in cloak
        RH_ice_en_remaining = FT(0), # env ice relative humidity not in cloak
    ))...,
    #
    qi_tendency_sub_dep_cloak_up = FT(0), # tendency due to sublimation/deposition in updraft cloak
    qi_tendency_sub_dep_cloak_dn = FT(0), # tendency due to sublimation/deposition in downdraft cloak
    qi_tendency_sub_dep_en_remaining = FT(0), # tendency due to sublimation/deposition not in cloak [[ we have the up and down cloak balance to not change this ]

    #
    q_liq_en_remaining = FT(0), # env liquid water not in cloak [[ we have the up and down cloak balance to not change this -- update i think it's better if we relax and include this for super sat reasons ]]
    q_ice_en_remaining = FT(0), # env ice water not in cloak [[ we have the up and down cloak balance to not change this -- update i think it's better if we relax and include this for super sat reasons ]]
    T_en_remaining = FT(0), # env temperature not in cloak [[ we have the up and down cloak balance to not change this -- update i think it's better if we relax and include this for super sat reasons ]]
    ts_en_remaining = thermo_state(FT, moisture_model), # env thermodynamic state not in cloak [[ we have the up and down cloak balance to not change this -- update i think it's better if we relax and include this for super sat reasons ]]

)

face_cloak_variables(FT, ::Nothing) =  NamedTuple()
face_cloak_variables(FT, ::StandardAreaPartitionModel) =  NamedTuple()
face_cloak_variables(FT, ::CoreCloakAreaPartitionModel) = (;
    w_cloak_up = FT(0), # env vertical velocity in updraft cloak
    w_cloak_dn = FT(0), # env vertical velocity in downdraft cloak
    # w_en_remaining = FT(0), # env vertical velocity not in cloak [[ this is either unchanged or 0 based on if `confine_all_downdraft_to_cloak` is true or false...]]
)

# consolidate my additions so it's easier to edit. These will exist everywhere... (though they may not be used or defined in eq case...)
my_microphysics_additions(FT, moisture_model::AbstractMoistureModel, cloud_sedimentation_model::Union{AbstractCloudSedimentationModel, Nothing}, area_partition_model::Union{AbstractAreaPartitionModel, Nothing}, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
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
    dTdz = FT(0), # store the vertical gradient of temperature for diagnostic purposes. Store here instead of aux_tc because we maybe want separately for env and up, though we don't do boost in the updraft so idk.
    #
    (calibrate_io ? (; ) : (;
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
        ))...,
    #
    # supersat
    supersat_variables(FT, moisture_model)...,
    #
    # cloud sedimentation
    cloud_sedimentation_variables(FT, cloud_sedimentation_model)...,
    # cloak variables
    cloak_variables(FT, area_partition_model, moisture_model, calibrate_io_val)..., # moisture_model is only for ts, so we'll see if we need it...
)

my_microphysics_additions(FT, mm::AbstractMoistureModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = my_microphysics_additions(FT, mm, nothing, nothing, calibrate_io_val)
my_microphysics_additions(FT, mm::AbstractMoistureModel, csm::AbstractCloudSedimentationModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = my_microphysics_additions(FT, mm, csm, nothing, calibrate_io_val)

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
cent_aux_vars_up_moisture(FT, moisture_model::NonEquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model, calibrate_io_val)...,
)
cent_aux_vars_up_moisture(FT, moisture_model::EquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model, calibrate_io_val)...,
)
cent_aux_vars_up(FT, local_geometry, edmf, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    ts = thermo_state(FT, edmf.moisture_model),
    q_liq = FT(0),
    q_ice = FT(0),
    T = FT(0),
    p = FT(0), # store this instead of constantly recalculating it [ both phasenonequil and phaseequil store ρ but phaseequil doesn't store p ] -- note typically in our case it should be the same as p_c but better safe than sorry.
    RH = FT(0),
    (calibrate_io ? (; ) : (; RH_liq = FT(0), RH_ice = FT(0)))..., # RH for liquid and ice, i.e. the ratio of q_vap to q_vap_sat_liq/ice
    s = FT(0),
    buoy = FT(0),
    area = FT(0),
    q_tot = FT(0),
    θ_liq_ice = FT(0),
    θ_liq_ice_tendency_precip_formation = FT(0),
    qt_tendency_precip_formation = FT(0),
    cent_aux_vars_up_moisture(FT, edmf.moisture_model, edmf.cloud_sedimentation_model, calibrate_io_val)...,
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
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::NonEquilibriumMoisture, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, calibrate_io_val)..., # don't need sedimentation in bulk
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::NonEquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model, calibrate_io_val)..., # added sedimentation to bulk for storing N_i etc...
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::EquilibriumMoisture, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    my_microphysics_additions(FT, moisture_model, calibrate_io_val)..., # don't need sedimentation in bulk
)
cent_aux_vars_edmf_bulk_moisture(FT, moisture_model::EquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model, calibrate_io_val)..., # added sedimentation to bulk for storing N_i etc...
)
cent_aux_vars_edmf_en_moisture(FT, moisture_model::NonEquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel, area_partition_model::AbstractAreaPartitionModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    ql_tendency_precip_formation = FT(0),
    qi_tendency_precip_formation = FT(0),
    ql_tendency_noneq = FT(0),
    qi_tendency_noneq = FT(0),
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model, area_partition_model, calibrate_io_val)..., # area partition for environment
)
cent_aux_vars_edmf_en_moisture(FT, moisture_model::EquilibriumMoisture, cloud_sedimentation_model::AbstractCloudSedimentationModel, area_partition_model::AbstractAreaPartitionModel, calibrate_io_val::Val{calibrate_io}) where {calibrate_io} = (;
    my_microphysics_additions(FT, moisture_model, cloud_sedimentation_model, area_partition_model, calibrate_io_val)..., # area partition for environment
)
cent_aux_vars_edmf_moisture(FT, ::NonEquilibriumMoisture, ::Val{calibrate_io}) where {calibrate_io} = (;
    # moved mass & diffusive ql,qi fluxes so is always there even for eq so can put in diagnostics output
)
cent_aux_vars_edmf_moisture(FT, ::EquilibriumMoisture, ::Val{calibrate_io}) where {calibrate_io} = NamedTuple()
cent_aux_vars_edmf(::Type{FT}, local_geometry, edmf, calibrate_io_val::Val{calibrate_io}) where {FT, calibrate_io} = (;
    turbconv = (;
        temporary_1 = FT(0), # temporary center variable for intermediate computations
        temporary_2 = FT(0), # temporary center variable for intermediate computations
        temporary_3 = FT(0), # temporary center variable for intermediate computations
        temporary_4 = FT(0), # temporary center variable for intermediate computations
        temporary_5 = FT(0), # temporary center variable for intermediate computations
        temporary_6 = FT(0), # temporary center variable for intermediate computations
        # temporary_locked = ntuple(i -> FT(1), 6), # whether the temporary variable is locked (i.e. in use) No way to use Bool
        k̂ = CCG.Contravariant3Vector(CCG.WVector(FT(1)), local_geometry),
        bulk = (;
            area = FT(0),
            θ_liq_ice = FT(0),
            RH = FT(0),
            (calibrate_io ? (; ) : (; RH_liq = FT(0), RH_ice = FT(0)))...,
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
            cent_aux_vars_edmf_bulk_moisture(FT, edmf.moisture_model, edmf.cloud_sedimentation_model, calibrate_io_val)..., # add sedimentation to bulk for storing N_i etc...
            # tendencies = cent_prognostic_vars_up(FT, edmf), # storage for tendencies limiter so we can reuse same memory
            # tendencies_adjustments = cent_prognostic_vars_up(FT, edmf), # storage for tendencies limiter so we can reuse same memory
            # prognostic = cent_prognostic_vars_up(FT, edmf), # storage for prognostic bulk w/ no clippings so we can be sure we have the true bulk for limiting
        ),
        up = ntuple(i -> cent_aux_vars_up(FT, local_geometry, edmf, calibrate_io_val), Val(n_updrafts(edmf))),
        en = (;
            ts = thermo_state(FT, edmf.moisture_model),
            area = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            θ_liq_ice = FT(0),
            θ_virt = FT(0),
            θ_dry = FT(0),
            RH = FT(0),
            (calibrate_io ? (; ) : (; RH_liq = FT(0), RH_ice = FT(0)))...,
            MSE = FT(0), # my addition
            ∂MSE∂z = FT(0), # my addition
            s = FT(0),
            T = FT(0),
            p = FT(0), # store this instead of constantly recalculating it [ both phasenonequil and phaseequil store ρ but phaseequil doesn't store p ] -- note typically in our case it should be the same as p_c but better safe than sorry.
            buoy = FT(0),
            cloud_fraction = FT(0),
            tke = FT(0),
            # Convert to convective TKE model soon
            tke_convective = FT(0),
            tke_convective_production = FT(0),
            tke_convective_advection = FT(0),
            tke_convective_dissipation = FT(0),
            latent_heating_pos = FT(0),
            latent_heating_neg = FT(0),
            instability = FT(0),
            stability = FT(0),
            frac_supersat = FT(0),
            frac_supersat_liq = FT(0),
            #
            Hvar = FT(0),
            QTvar = FT(0),
            HQTcov = FT(0),
            θ_liq_ice_tendency_precip_formation = FT(0),
            qt_tendency_precip_formation = FT(0),
            cent_aux_vars_edmf_en_moisture(FT, edmf.moisture_model, edmf.cloud_sedimentation_model, edmf.area_partition_model, calibrate_io_val)...,
            unsat = (; q_tot = FT(0), θ_dry = FT(0), θ_virt = FT(0)),
            sat = (; T = FT(0), q_vap = FT(0), q_tot = FT(0), θ_dry = FT(0), θ_liq_ice = FT(0), q_liq = FT(0), q_ice = FT(0)), # added liq and ice
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
        en_2m = (;
            tke = cent_aux_vars_en_2m(FT),
            Hvar = cent_aux_vars_en_2m(FT),
            QTvar = cent_aux_vars_en_2m(FT),
            HQTcov = cent_aux_vars_en_2m(FT),
            #
            # tke_convective = cent_aux_vars_en_2m(FT), # tke generated by convective processes that is non local so we do not dissipate it until it reaches strong stability... [[ this may be overkill because we don't track all the terms... It is only buoyancy generated, slowly dissipated and stability lost...]]
        ),
        term_vel_rain = FT(0),
        term_vel_snow = FT(0),
        # sedimentation = NoneqMoistureSource{FT}(0), # storage for sed [[ just gonna use a scratch field]]
        # sedimentation_other = NoneqMoistureSource{FT}(0) [[ just gonna use a scratch field ]]
        KM = FT(0),
        KH = FT(0),
        KQ = FT(0),
        mixing_length = FT(0),
        dqvdt = FT(0), # store the tendency of water vapor for diagnostic purposes [[ grid mean only ]]
        dTdt = FT(0), # store the tendency of temperature for diagnostic purposes [[ grid mean only ]]
        ∂b∂z = FT(0), # store this so we can try to use it to generate tke even if there's not existing tke
        (calibrate_io ? (; ) : (;
            qs_tendency_accr_rai_sno = FT(0), # tendency due to accretion QR by QS [QR -> QS]  (we store it this way, it's always to snow below freezing and [QS -> QR] above freezing.) Just for storage...
            #
            diffusive_tendency_h = FT(0),
            diffusive_tendency_qt = FT(0),
            diffusive_tendency_ql = FT(0), # my addition [ this is the diffusive tendency of ql, not including  massflux tendency ]
            diffusive_tendency_qi = FT(0), # my addition [ this is the
            diffusive_tendency_qr = FT(0), # my addition [ this is the diffusive tendency of qr, not including  massflux tendency ]
            diffusive_tendency_qs = FT(0), # my addition [ this is the diffusive tendency of qs, not including  massflux tendency ]
        ))...,
        massflux = FT(0),
        massflux_tendency_h = FT(0),
        massflux_tendency_qt = FT(0),
        massflux_tendency_ql = FT(0), # moving to here so they're always there even for eq so can put in diagnostics output
        massflux_tendency_qi = FT(0), # moving to here so they're always there even for eq so can put in diagnostics output

        #
        #
        cent_aux_vars_edmf_moisture(FT, edmf.moisture_model, calibrate_io_val)...,
        prandtl_nvec = FT(0),
        # Added by Ignacio : Length scheme in use (mls), and smooth min effect (ml_ratio)
        # Variable Prandtl number initialized as neutral value.
        mls = FT(0),
        b_exch = FT(0),
        ml_ratio = FT(0),
        w_up_c = FT(0), # gets reused for all updrafts...
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
        # l_entdet = FT(0), # not used
        ϕ_gm = FT(0), # temporary for grid-mean variables
        ϕ_gm_cov = FT(0), # temporary for grid-mean covariance variables
        ϕ_en_cov = FT(0), # temporary for environmental covariance variables
        (calibrate_io ? (;) : (; 
            ϕ_up_cubed = FT(0), # temporary for cubed updraft variables in grid mean 3rd moment functions
        ))...,
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
    # massflux_en = FT(0), # TODO: is this the right place for this?
    massflux_ql = FT(0),
    massflux_qi = FT(0),
    diffusive_flux_ql = FT(0),
    diffusive_flux_qi = FT(0),
)
face_aux_vars_edmf_moisture(FT, ::EquilibriumMoisture) = NamedTuple()
face_aux_vars_edmf(::Type{FT}, local_geometry, edmf, ::Val{calibrate_io}) where {FT, calibrate_io} = (;
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
            face_cloak_variables(FT, edmf.area_partition_model)...,
        ),
        up = ntuple(i -> face_aux_vars_up(FT, local_geometry), Val(n_updrafts(edmf))),
        massflux = FT(0),
        massflux_en = FT(0), # TODO: is this the right place for this?
        massflux_h = FT(0),
        massflux_qt = FT(0),
        temporary_f1 = FT(0), # temporary face variable for intermediate computations
        temporary_f2 = FT(0), # temporary face variable for intermediate computations
        temporary_f3 = FT(0), # temporary face variable for intermediate computations
        temporary_f4 = FT(0), # temporary face variable for intermediate computations
        # temporary_f_locked = ntuple(i -> FT(1), Val(4)), # temporary variable for locking temporary variables (No way to use bool)
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
cent_diagnostic_vars_edmf(FT, local_geometry, edmf, ::Val{calibrate_io}) where {calibrate_io} = (;
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
face_diagnostic_vars_edmf(FT, local_geometry, edmf, ::Val{calibrate_io}) where {calibrate_io} = (;
    turbconv = (; nh_pressure = FT(0), nh_pressure_adv = FT(0), nh_pressure_drag = FT(0), nh_pressure_b = FT(0)),
    precip = calibrate_io ? (;) : (; rain_flux = FT(0), snow_flux = FT(0)),
)

# Single value per column diagnostic variables
single_value_per_col_diagnostic_vars_edmf(FT, edmf, ::Val{calibrate_io}) where {calibrate_io} = (;
    turbconv = calibrate_io ? (;) : (;
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

##### Prognostic fields [[ These should not be calibrate_io dependent ]]

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
    (; ρatke = FT(0), cent_prognostic_vars_en_thermo(FT, edmf.thermo_covariance_model)..., ρatke_convective = FT(0))
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
