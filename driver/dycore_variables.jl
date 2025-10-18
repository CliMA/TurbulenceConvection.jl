#####
##### Fields
#####

import TurbulenceConvection as TC
import ClimaCore as CC
import ClimaCore.Geometry as CCG
import ClimaCore.Geometry: ⊗

##### Auxiliary fields

# Center only
cent_aux_vars_gm_moisture(FT, ::TC.NonEquilibriumMoisture) = (;
    ∇q_liq_gm = FT(0),
    ∇q_ice_gm = FT(0),
    dqldt_rad = FT(0),
    dqidt_rad = FT(0),
    dqldt = FT(0),
    dqidt = FT(0),
    dqldt_hadv = FT(0),
    dqidt_hadv = FT(0),
    ql_nudge = FT(0),
    qi_nudge = FT(0),
    dqldt_fluc = FT(0),
    dqidt_fluc = FT(0),
)
cent_aux_vars_gm_moisture(FT, ::TC.EquilibriumMoisture) = NamedTuple()
cent_aux_vars_gm(FT, local_geometry, edmf) = (;
    ts = TC.thermo_state(FT, edmf.moisture_model),
    tke = FT(0),
    Hvar = FT(0),
    QTvar = FT(0),
    HQTcov = FT(0),
    q_liq = FT(0),
    q_ice = FT(0),
    RH = FT(0),
    RH_liq = FT(0),
    RH_ice = FT(0), # RH for ice, i.e. the ratio of q_vap to q_vap_sat_ice
    s = FT(0),
    T = FT(0),
    buoy = FT(0),
    cloud_fraction = FT(0),
    H_third_m = FT(0),
    W_third_m = FT(0),
    QT_third_m = FT(0),
    # From RadiationBase
    dTdt_rad = FT(0), # horizontal advection temperature tendency
    dqtdt_rad = FT(0), # horizontal advection moisture tendency
    # From ForcingBase
    subsidence = FT(0), #Large-scale subsidence
    dTdt_hadv = FT(0), #Horizontal advection of temperature
    dqtdt_hadv = FT(0), #Horizontal advection of moisture
    H_nudge = FT(0), #Reference thetali profile for relaxation tendency
    qt_nudge = FT(0), #Reference qt profile for relaxation tendency
    dTdt_fluc = FT(0), #Vertical turbulent advection of temperature
    dqtdt_fluc = FT(0), #Vertical turbulent advection of moisture
    u_nudge = FT(0), #Reference u profile for relaxation tendency
    v_nudge = FT(0), #Reference v profile for relaxation tendency
    uₕ_g = CCG.Covariant12Vector(CCG.UVVector(FT(0), FT(0)), local_geometry), #Geostrophic u velocity
    ∇θ_liq_ice_gm = FT(0),
    ∇q_tot_gm = FT(0),
    cent_aux_vars_gm_moisture(FT, edmf.moisture_model)...,
    θ_virt = FT(0),
    Ri = FT(0),
    θ_liq_ice = FT(0),
    q_tot = FT(0),
    p = FT(0),
    qt_tendency_ls_vert_adv = FT(0), # my addition
    ql_tendency_ls_vert_adv = FT(0), # my addition
    qi_tendency_ls_vert_adv = FT(0), # my addition
    qr_tendency_ls_vert_adv = FT(0), # my addition [[ not sure how relevant this one is but ]]
    qs_tendency_ls_vert_adv = FT(0), # my addition [[ not sure how relevant this one is but ]]
    # qt_tendency_vert_adv = error("not implemented, see massflux_tendency_qt")
    ql_tendency_vert_adv = FT(0), # my addition [ For full flux, massflux is only differential between up/env and gm so sgs ]
    qi_tendency_vert_adv = FT(0), # my addition [ For full flux, massflux is only differential between up/env and gm so sgs ]
    qr_tendency_vert_adv = FT(0), # my addition 
    qs_tendency_vert_adv = FT(0), # my addition
    qr_tendency_sedimentation = FT(0), # storage
    qs_tendency_sedimentation = FT(0), # storage
    sgs_tendency_q_liq = FT(0), # testing here instead of below in face_aux_vars_gm_moisture() to see if it's the right spaces type... cause after deriv we're back on center no longer on faces...
    sgs_tendency_q_ice = FT(0), # testing here instead of below in face_aux_vars_gm_moisture() to see if it's the right spaces type... cause after deriv we're back on center no longer on faces...
    sgs_tendency_q_rai = FT(0), # my addition
    sgs_tendency_q_sno = FT(0), # my addition
    # term_vel_liq = FT(0), # for grid mean sed.. .would be nice to have that conditionally but that breaks eltype inference, could make a function w/ sedimentatinon_vars that is based on edmf.cloud_sedimentation_model
    # term_vel_ice = FT(0), # for grid mean sed.. .would be nice to have that conditionally but that breaks eltype inference, could make a function w/ sedimentatinon_vars that is based on edmf.cloud_sedimentation_model
    TC.cloud_sedimentation_variables(FT, edmf.cloud_sedimentation_model)...,
    f_ice_mult = FT(1),
)
cent_aux_vars(FT, local_geometry, edmf) =
    (; cent_aux_vars_gm(FT, local_geometry, edmf)..., TC.cent_aux_vars_edmf(FT, local_geometry, edmf)...)

# Face only
face_aux_vars_gm_moisture(FT, ::TC.NonEquilibriumMoisture) = (;
    sgs_flux_q_liq = FT(0),
    sgs_flux_q_ice = FT(0),
    ql_flux_vert_adv = FT(0), # my addition [ For full flux, massflux is only differential between up/env and gm so sgs , could also go in aux_tc_f [face_aux_vars_edmf_moisture(FT, ::NonEquilibriumMoisture) ]
    qi_flux_vert_adv = FT(0), # my addition [ For full flux, massflux is only differential between up/env and gm so sgs , could also go in aux_tc_f [face_aux_vars_edmf_moisture(FT, ::NonEquilibriumMoisture) ]
)
face_aux_vars_gm_moisture(FT, ::TC.EquilibriumMoisture) = NamedTuple()
face_aux_vars_gm(FT, local_geometry, edmf) = (;
    massflux_s = FT(0),
    diffusive_flux_s = FT(0),
    total_flux_s = FT(0),
    f_rad = FT(0),
    sgs_flux_θ_liq_ice = FT(0),
    sgs_flux_q_tot = FT(0),
    sgs_flux_q_rai = FT(0), # my addition
    sgs_flux_q_sno = FT(0), # my addition
    face_aux_vars_gm_moisture(FT, edmf.moisture_model)...,
    sgs_flux_uₕ = CCG.Covariant3Vector(FT(0)) ⊗ CCG.Covariant12Vector(FT(0), FT(0)),
    p = FT(0),
    ρ = FT(0),
)
face_aux_vars(FT, local_geometry, edmf) =
    (; face_aux_vars_gm(FT, local_geometry, edmf)..., TC.face_aux_vars_edmf(FT, local_geometry, edmf)...)

##### Diagnostic fields

# Center only
cent_diagnostic_vars_gm(FT, local_geometry) = NamedTuple()
cent_diagnostic_vars(FT, local_geometry, edmf) =
    (; cent_diagnostic_vars_gm(FT, local_geometry)..., TC.cent_diagnostic_vars_edmf(FT, local_geometry, edmf)...)

# Face only
face_diagnostic_vars_gm(FT, local_geometry) = NamedTuple()
face_diagnostic_vars(FT, local_geometry, edmf) =
    (; face_diagnostic_vars_gm(FT, local_geometry)..., TC.face_diagnostic_vars_edmf(FT, local_geometry, edmf)...)

# Single value per column diagnostic variables
single_value_per_col_diagnostic_vars_gm(FT) = (;
    Tsurface = FT(0),
    shf = FT(0),
    lhf = FT(0),
    ustar = FT(0),
    wstar = FT(0),
    lwp_mean = FT(0),
    iwp_mean = FT(0),
    rwp_mean = FT(0),
    swp_mean = FT(0),
    cutoff_precipitation_rate = FT(0),
    cloud_base_mean = FT(0),
    cloud_top_mean = FT(0),
    cloud_cover_mean = FT(0),
    integ_total_flux_qt = FT(0),
    integ_total_flux_s = FT(0),
)
single_value_per_col_diagnostic_vars(FT, edmf) =
    (; single_value_per_col_diagnostic_vars_gm(FT)..., TC.single_value_per_col_diagnostic_vars_edmf(FT, edmf)...)

##### Prognostic fields

# Center only
cent_prognostic_vars(::Type{FT}, local_geometry, edmf) where {FT} =
    (; cent_prognostic_vars_gm(FT, local_geometry, edmf)..., TC.cent_prognostic_vars_edmf(FT, edmf)...,)
cent_prognostic_vars_gm_moisture(::Type{FT}, ::TC.NonEquilibriumMoisture) where {FT} = (; q_liq = FT(0), q_ice = FT(0))
cent_prognostic_vars_gm_moisture(::Type{FT}, ::TC.EquilibriumMoisture) where {FT} = NamedTuple()
cent_prognostic_vars_gm(::Type{FT}, local_geometry, edmf) where {FT} = (;
    ρ = FT(0),
    uₕ = CCG.Covariant12Vector(CCG.UVVector(FT(0), FT(0)), local_geometry),
    ρθ_liq_ice = FT(0),
    ρq_tot = FT(0),
    cent_prognostic_vars_gm_moisture(FT, edmf.moisture_model)...,
)

# Face only
face_prognostic_vars(::Type{FT}, local_geometry, edmf) where {FT} =
    (; w = CCG.Covariant3Vector(FT(0)), TC.face_prognostic_vars_edmf(FT, local_geometry, edmf)...)

# TC.face_prognostic_vars_edmf(FT, edmf) = (;) # could also use this for empty model
