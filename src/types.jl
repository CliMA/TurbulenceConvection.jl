Base.@kwdef mutable struct pressure_drag_struct
    nh_pressure_adv::Float64 = 0
    nh_pressure_drag::Float64 = 0
end

Base.@kwdef mutable struct pressure_buoy_struct
    b_coeff::Float64 = 0
    nh_pressure_b::Float64 = 0
end

Base.@kwdef mutable struct pressure_in_struct
    updraft_top::Float64 = 0
    asp_label::Float64 = 0
    a_med::Float64 = 0
    a_kfull::Float64 = 0
    b_kfull::Float64 = 0
    rho0_kfull::Float64 = 0
    alpha1::Float64 = 0
    alpha2::Float64 = 0
    beta1::Float64 = 0
    beta2::Float64 = 0
    w_kfull::Float64 = 0
    w_kmfull::Float64 = 0
    w_kenv::Float64 = 0
    dzi::Float64 = 0
    drag_sign::Float64 = 0
    asp_ratio::Float64 = 0
end

Base.@kwdef mutable struct entr_struct
    entr_sc::Float64 = 0
    detr_sc::Float64 = 0
    sorting_function::Float64 = 0
    b_mix::Float64 = 0
end

Base.@kwdef mutable struct eos_struct
    T::Float64 = 0
    ql::Float64 = 0
end

Base.@kwdef mutable struct rain_struct
    qr::Float64 = 0
    ar::Float64 = 0
end

Base.@kwdef mutable struct mph_struct
    T::Float64 = 0
    thl::Float64 = 0
    th::Float64 = 0
    rho::Float64 = 0
    qt::Float64 = 0
    qv::Float64 = 0
    ql::Float64 = 0
    thl_rain_src::Float64 = 0
    qr_src::Float64 = 0
end

Base.@kwdef mutable struct entr_in_struct
    c_div::Float64 = 0
    M::Float64 = 0
    dMdz::Float64 = 0
    zi::Float64 = 0
    wstar::Float64 = 0
    z::Float64 = 0
    sort_pow::Float64 = 0
    c_det::Float64 = 0
    chi_upd::Float64 = 0
    tke_coef::Float64 = 0
    c_mu::Float64 = 0
    c_ed_mf::Float64 = 0
    c_mu0::Float64 = 0
    dz::Float64 = 0
    w_upd::Float64 = 0
    b_upd::Float64 = 0
    c_ent::Float64 = 0
    dt::Float64 = 0
    b_env::Float64 = 0
    a_upd::Float64 = 0
    a_env::Float64 = 0
    tke::Float64 = 0
    RH_upd::Float64 = 0
    RH_env::Float64 = 0
    qt_up::Float64 = 0
    ql_up::Float64 = 0
    T_env::Float64 = 0
    qt_env::Float64 = 0
    ql_env::Float64 = 0
    w_env::Float64 = 0
    nh_pressure::Float64 = 0
    L::Float64 = 0
    zbl::Float64 = 0
    poisson::Float64 = 0
    buoy_ed_flux::Float64 = 0
end

struct RainVariable{T}
    loc::String
    kind::String
    name::String
    units::String
    values::T
    new::T
    flux::T
    function RainVariable(grid, name, units)
        loc = "half"
        kind = "scalar"

        values = center_field(grid)
        new = center_field(grid)
        flux = center_field(grid)
        return new{typeof(values)}(loc, kind, name, units, values, new, flux)
    end
end

Base.@kwdef mutable struct RainVariables
    rain_model::String = "default_rain_model"
    max_supersaturation::Float64
    C_drag::Float64
    MP_n_0::Float64
    tau_cond_evap::Float64
    q_liq_threshold::Float64
    tau_acnv::Float64
    E_col::Float64
    a_vent::Float64
    b_vent::Float64
    mean_rwp::Float64 = 0
    env_rwp::Float64 = 0
    upd_rwp::Float64 = 0
    cutoff_rain_rate::Float64 = 0
    Gr::Grid
    QR::RainVariable
    RainArea::RainVariable
    Upd_QR::RainVariable
    Upd_RainArea::RainVariable
    Env_QR::RainVariable
    Env_RainArea::RainVariable
end
function RainVariables(namelist, Gr::Grid)

    QR = RainVariable(Gr, "qr_mean", "kg/kg")
    # temporary variables for diagnostics to know where the rain is coming from
    Upd_QR = RainVariable(Gr, "upd_qr", "kg/kg")
    Env_QR = RainVariable(Gr, "env_qr", "kg/kg")
    # in the future we could test prognostic equations for stratiform and updraft rain
    RainArea = RainVariable(Gr, "rain_area", "rain_area_fraction [-]")
    Upd_RainArea = RainVariable(Gr, "upd_rain_area", "updraft_rain_area_fraction [-]")
    Env_RainArea = RainVariable(Gr, "env_rain_area", "environment_rain_area_fraction [-]")

    mean_rwp = 0.0
    upd_rwp = 0.0
    env_rwp = 0.0

    rain_model = try
        string(namelist["microphysics"]["rain_model"])
    catch
        println("EDMF_Rain: defaulting to no rain")
        "None"
    end

    if !(rain_model in ["None", "cutoff", "clima_1m"])
        error("rain model not recognized")
    end

    max_supersaturation = try
        namelist["microphysics"]["max_supersaturation"]
    catch
        0.02
    end

    C_drag = try
        namelist["microphysics"]["C_drag"]
    catch
        0.55
    end

    MP_n_0 = try
        namelist["microphysics"]["MP_n_0"]
    catch
        16 * 1e6
    end

    tau_cond_evap = try
        namelist["microphysics"]["tau_cond_evap"]
    catch
        10.0
    end

    q_liq_threshold = try
        namelist["microphysics"]["q_liq_threshold"]
    catch
        5e-4
    end

    tau_acnv = try
        namelist["microphysics"]["tau_acnv"]
    catch
        1e3
    end

    E_col = try
        namelist["microphysics"]["E_col"]
    catch
        0.8
    end

    a_vent = try
        namelist["microphysics"]["a_vent"]
    catch
        1.5
    end

    b_vent = try
        namelist["microphysics"]["b_vent"]
    catch
        0.53
    end

    return RainVariables(;
        rain_model,
        max_supersaturation,
        C_drag,
        MP_n_0,
        tau_cond_evap,
        q_liq_threshold,
        tau_acnv,
        E_col,
        a_vent,
        b_vent,
        mean_rwp,
        env_rwp,
        upd_rwp,
        Gr,
        QR,
        RainArea,
        Upd_QR,
        Upd_RainArea,
        Env_QR,
        Env_RainArea,
    )
end

struct VariablePrognostic{T}
    values::T
    new::T
    mf_update::T
    tendencies::T
    loc::String
    bc::String
    kind::String
    name::String
    units::String
    function VariablePrognostic(grid, loc, kind, bc, name, units)
        # Value at the current timestep
        values = field(grid, loc)
        # Value at the next timestep, used for calculating turbulence tendencies
        new = field(grid, loc)
        mf_update = field(grid, loc)
        tendencies = field(grid, loc)
        if kind != "scalar" && kind != "velocity"
            print("Invalid kind setting for variable! Must be scalar or velocity")
        end
        return new{typeof(values)}(values, new, mf_update, tendencies, loc, bc, kind, name, units)
    end
end

struct VariableDiagnostic{T}
    values::T
    loc::String
    bc::String
    kind::String
    name::String
    units::String
    function VariableDiagnostic(grid, loc, kind, bc, name, units)
        # Value at the current timestep
        values = field(grid, loc)
        # Placement on staggered grid
        if kind != "scalar" && kind != "velocity"
            print("Invalid kind setting for variable! Must be scalar or velocity")
        end
        return new{typeof(values)}(values, loc, bc, kind, name, units)
    end
end

struct UpdraftVariable{A1, A2}
    values::A2
    old::A2
    new::A2
    tendencies::A2
    flux::A2
    bulkvalues::A1
    loc::String
    kind::String
    name::String
    units::String
    function UpdraftVariable(grid, nu, loc, kind, name, units)
        values = field(grid, loc, nu)
        old = field(grid, loc, nu)  # needed for prognostic updrafts
        new = field(grid, loc, nu) # needed for prognostic updrafts
        tendencies = field(grid, loc, nu)
        flux = field(grid, loc, nu)
        bulkvalues = field(grid, loc)
        if kind != "scalar" && kind != "velocity"
            print("Invalid kind setting for variable! Must be scalar or velocity")
        end
        A1 = typeof(bulkvalues)
        A2 = typeof(values)
        return new{A1, A2}(values, old, new, tendencies, flux, bulkvalues, loc, kind, name, units)
    end
end

mutable struct UpdraftVariables{A1}
    Gr::Grid
    n_updrafts::Int
    W::UpdraftVariable
    Area::UpdraftVariable
    QT::UpdraftVariable
    QL::UpdraftVariable
    RH::UpdraftVariable
    H::UpdraftVariable
    THL::UpdraftVariable
    T::UpdraftVariable
    B::UpdraftVariable
    prognostic::Bool
    updraft_fraction::Float64
    cloud_fraction::A1
    cloud_base::A1
    cloud_top::A1
    cloud_cover::A1
    updraft_top::A1
    lwp::Float64
    function UpdraftVariables(nu, namelist, grid::Grid)
        n_updrafts = nu

        W = UpdraftVariable(grid, nu, "full", "velocity", "w", "m/s")

        Area = UpdraftVariable(grid, nu, "half", "scalar", "area_fraction", "[-]")
        QT = UpdraftVariable(grid, nu, "half", "scalar", "qt", "kg/kg")
        QL = UpdraftVariable(grid, nu, "half", "scalar", "ql", "kg/kg")
        RH = UpdraftVariable(grid, nu, "half", "scalar", "RH", "%")
        H = UpdraftVariable(grid, nu, "half", "scalar", "thetal", "K")
        THL = UpdraftVariable(grid, nu, "half", "scalar", "thetal", "K")
        T = UpdraftVariable(grid, nu, "half", "scalar", "temperature", "K")
        B = UpdraftVariable(grid, nu, "half", "scalar", "buoyancy", "m^2/s^3")

        prognostic = true
        updraft_fraction = namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"]

        # cloud and rain diagnostics for output
        cloud_fraction = center_field(grid)

        cloud_base = pyzeros(nu)
        cloud_top = pyzeros(nu)
        cloud_cover = pyzeros(nu)
        updraft_top = pyzeros(nu)

        lwp = 0.0
        A1 = typeof(cloud_fraction)
        return new{A1}(
            grid,
            n_updrafts,
            W,
            Area,
            QT,
            QL,
            RH,
            H,
            THL,
            T,
            B,
            prognostic,
            updraft_fraction,
            cloud_fraction,
            cloud_base,
            cloud_top,
            cloud_cover,
            updraft_top,
            lwp,
        )
    end
end

Base.@kwdef mutable struct GridMeanVariables{PS}
    param_set::PS
    Gr::Grid
    Ref::ReferenceState
    lwp::Float64
    cloud_base::Float64
    cloud_top::Float64
    cloud_cover::Float64
    U::VariablePrognostic
    V::VariablePrognostic
    W::VariablePrognostic
    QT::VariablePrognostic
    RH::VariablePrognostic
    H::VariablePrognostic
    t_to_prog_fp::Function
    prog_to_t_fp::Function
    QL::VariableDiagnostic
    T::VariableDiagnostic
    B::VariableDiagnostic
    THL::VariableDiagnostic
    cloud_fraction::VariableDiagnostic
    calc_tke::Bool
    calc_scalar_var::Bool
    EnvThermo_scheme::String
    TKE::VariableDiagnostic
    W_third_m::VariableDiagnostic
    QTvar::VariableDiagnostic
    QT_third_m::VariableDiagnostic
    Hvar::VariableDiagnostic
    H_third_m::VariableDiagnostic
    HQTcov::VariableDiagnostic
end
function GridMeanVariables(namelist, Gr::Grid, Ref::ReferenceState, param_set::PS) where {PS}
    lwp = 0.0
    cloud_base = 0.0
    cloud_top = 0.0
    cloud_cover = 0.0

    U = VariablePrognostic(Gr, "half", "velocity", "sym", "u", "m/s")
    V = VariablePrognostic(Gr, "half", "velocity", "sym", "v", "m/s")
    # Just leave this zero for now!
    W = VariablePrognostic(Gr, "full", "velocity", "asym", "v", "m/s")

    # Create thermodynamic variables
    QT = VariablePrognostic(Gr, "half", "scalar", "sym", "qt", "kg/kg")
    RH = VariablePrognostic(Gr, "half", "scalar", "sym", "RH", "%")

    H = VariablePrognostic(Gr, "half", "scalar", "sym", "thetal", "K")
    t_to_prog_fp = t_to_thetali_c
    prog_to_t_fp = eos_first_guess_thetal

    # Diagnostic Variables--same class as the prognostic variables, but we append to diagnostics list
    QL = VariableDiagnostic(Gr, "half", "scalar", "sym", "ql", "kg/kg")
    T = VariableDiagnostic(Gr, "half", "scalar", "sym", "temperature", "K")
    B = VariableDiagnostic(Gr, "half", "scalar", "sym", "buoyancy", "m^2/s^3")
    THL = VariableDiagnostic(Gr, "half", "scalar", "sym", "thetal", "K")

    cloud_fraction = VariableDiagnostic(Gr, "half", "scalar", "sym", "cloud fraction", "-")

    # TKE   TODO   repeated from EDMF_Environment.pyx logic
    calc_tke = true
    try
        calc_tke = namelist["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"]
    catch
    end

    calc_scalar_var = try
        namelist["turbulence"]["EDMF_PrognosticTKE"]["calc_scalar_var"]
    catch
        false
    end

    EnvThermo_scheme = try
        string(namelist["thermodynamics"]["sgs"])
    catch
        "mean"
    end

    #Now add the 2nd moment variables
    if calc_tke
        TKE = VariableDiagnostic(Gr, "half", "scalar", "sym", "tke", "m^2/s^2")
        W_third_m = VariableDiagnostic(Gr, "half", "scalar", "sym", "W_third_m", "m^3/s^3")
    end

    if calc_scalar_var
        QTvar = VariableDiagnostic(Gr, "half", "scalar", "sym", "qt_var", "kg^2/kg^2")
        QT_third_m = VariableDiagnostic(Gr, "half", "scalar", "sym", "qt_third_m", "kg^3/kg^3")
        Hvar = VariableDiagnostic(Gr, "half", "scalar", "sym", "thetal_var", "K^2")
        H_third_m = VariableDiagnostic(Gr, "half", "scalar", "sym", "thetal_third_m", "-")
        HQTcov = VariableDiagnostic(Gr, "half", "scalar", "sym", "thetal_qt_covar", "K(kg/kg)")
    end

    return GridMeanVariables(;
        param_set,
        Gr,
        Ref,
        lwp,
        cloud_base,
        cloud_top,
        cloud_cover,
        U,
        V,
        W,
        QT,
        RH,
        H,
        t_to_prog_fp,
        prog_to_t_fp,
        QL,
        T,
        B,
        THL,
        cloud_fraction,
        calc_tke,
        calc_scalar_var,
        EnvThermo_scheme,
        TKE,
        W_third_m,
        QTvar,
        QT_third_m,
        Hvar,
        H_third_m,
        HQTcov,
    )
end


struct UpdraftThermodynamics{A1, A2}
    Gr::Grid
    Ref::ReferenceState
    n_updraft::Int
    t_to_prog_fp::Function
    prog_to_t_fp::Function
    prec_source_h::A2
    prec_source_qt::A2
    prec_source_h_tot::A1
    prec_source_qt_tot::A1
    function UpdraftThermodynamics(
        n_updraft::Int,
        Gr::Grid,
        Ref::ReferenceState,
        UpdVar::UpdraftVariables,
        Rain::RainVariables,
    )
        t_to_prog_fp = t_to_thetali_c
        prog_to_t_fp = eos_first_guess_thetal

        # rain source from each updraft from all sub-timesteps
        prec_source_h = center_field(Gr, n_updraft)
        prec_source_qt = center_field(Gr, n_updraft)

        # rain source from all updrafts from all sub-timesteps
        prec_source_h_tot = center_field(Gr)
        prec_source_qt_tot = center_field(Gr)
        A1 = typeof(prec_source_h_tot)
        A2 = typeof(prec_source_h)
        return new{A1, A2}(
            Gr,
            Ref,
            n_updraft,
            t_to_prog_fp,
            prog_to_t_fp,
            prec_source_h,
            prec_source_qt,
            prec_source_h_tot,
            prec_source_qt_tot,
        )
    end

end

struct EnvironmentVariable{T}
    values::T
    flux::T
    loc::String
    kind::String
    name::String
    units::String
    function EnvironmentVariable(grid, loc, kind, name, units)
        values = field(grid, loc)
        flux = field(grid, loc)
        if kind != "scalar" && kind != "velocity"
            println("Invalid kind setting for variable! Must be scalar or velocity")
        end
        return new{typeof(values)}(values, flux, loc, kind, name, units)
    end
end

struct EnvironmentVariable_2m{A1}
    values::A1
    dissipation::A1
    shear::A1
    entr_gain::A1
    detr_loss::A1
    press::A1
    buoy::A1
    interdomain::A1
    rain_src::A1
    loc::String
    kind::String
    name::String
    units::String
    function EnvironmentVariable_2m(grid, loc, kind, name, units)
        values = center_field(grid)
        dissipation = center_field(grid)
        entr_gain = center_field(grid)
        detr_loss = center_field(grid)
        buoy = center_field(grid)
        press = center_field(grid)
        shear = center_field(grid)
        interdomain = center_field(grid)
        rain_src = center_field(grid)
        if kind != "scalar" && kind != "velocity"
            println("Invalid kind setting for variable! Must be scalar or velocity")
        end
        return new{typeof(values)}(
            values,
            dissipation,
            shear,
            entr_gain,
            detr_loss,
            press,
            buoy,
            interdomain,
            rain_src,
            loc,
            kind,
            name,
            units,
        )
    end
end

Base.@kwdef mutable struct EnvironmentVariables{PS}
    param_set::PS
    Gr::Grid
    W::EnvironmentVariable
    Area::EnvironmentVariable
    QT::EnvironmentVariable
    QL::EnvironmentVariable
    H::EnvironmentVariable
    THL::EnvironmentVariable
    RH::EnvironmentVariable
    T::EnvironmentVariable
    B::EnvironmentVariable
    cloud_fraction::EnvironmentVariable
    TKE::EnvironmentVariable_2m
    Hvar::EnvironmentVariable_2m
    QTvar::EnvironmentVariable_2m
    HQTcov::EnvironmentVariable_2m
    calc_tke::Bool = false
    calc_scalar_var::Bool = false
    cloud_base::Float64 = 0
    cloud_top::Float64 = 0
    cloud_cover::Float64 = 0
    lwp::Float64 = 0
    EnvThermo_scheme::String = "default_EnvThermo_scheme"
end
function EnvironmentVariables(namelist, Gr::Grid, param_set::PS) where {PS}
    W = EnvironmentVariable(Gr, "full", "velocity", "w", "m/s")
    QT = EnvironmentVariable(Gr, "half", "scalar", "qt", "kg/kg")
    QL = EnvironmentVariable(Gr, "half", "scalar", "ql", "kg/kg")
    RH = EnvironmentVariable(Gr, "half", "scalar", "RH", "%")
    H = EnvironmentVariable(Gr, "half", "scalar", "thetal", "K")
    THL = EnvironmentVariable(Gr, "half", "scalar", "thetal", "K")
    T = EnvironmentVariable(Gr, "half", "scalar", "temperature", "K")
    B = EnvironmentVariable(Gr, "half", "scalar", "buoyancy", "m^2/s^3")
    Area = EnvironmentVariable(Gr, "half", "scalar", "env_area", "-")
    cloud_fraction = EnvironmentVariable(Gr, "half", "scalar", "env_cloud_fraction", "-")

    # TODO - the flag setting is repeated from Variables.pyx logic
    calc_tke = true
    calc_tke = try
        namelist["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"]
    catch
        nothing
    end

    calc_scalar_var = try
        namelist["turbulence"]["EDMF_PrognosticTKE"]["calc_scalar_var"]
    catch
        println("Defaulting to non-calculation of scalar variances")
        false
    end

    EnvThermo_scheme = try
        string(namelist["thermodynamics"]["sgs"])
    catch
        println("Defaulting to saturation adjustment and microphysics with respect to environmental means")
        "mean"
    end

    if calc_tke
        TKE = EnvironmentVariable_2m(Gr, "half", "scalar", "tke", "m^2/s^2")
    end

    if calc_scalar_var
        QTvar = EnvironmentVariable_2m(Gr, "half", "scalar", "qt_var", "kg^2/kg^2")
        Hvar = EnvironmentVariable_2m(Gr, "half", "scalar", "thetal_var", "K^2")
        HQTcov = EnvironmentVariable_2m(Gr, "half", "scalar", "thetal_qt_covar", "K(kg/kg)")
    end

    if EnvThermo_scheme == "quadrature"
        if (calc_scalar_var == false)
            error("EDMF_Environment.pyx: scalar variance has to be calculated for quadrature saturation and microphysics")
        end
    end

    return EnvironmentVariables{PS}(;
        param_set,
        Gr,
        W,
        Area,
        QT,
        QL,
        H,
        THL,
        RH,
        T,
        B,
        cloud_fraction,
        TKE,
        Hvar,
        QTvar,
        HQTcov,
        calc_tke,
        calc_scalar_var,
        EnvThermo_scheme,
    )
end

struct EnvironmentThermodynamics{A1}
    Gr::Grid
    Ref::ReferenceState
    quadrature_order::Int
    quadrature_type::String
    t_to_prog_fp::Function
    prog_to_t_fp::Function
    qt_dry::A1
    th_dry::A1
    t_cloudy::A1
    qv_cloudy::A1
    qt_cloudy::A1
    th_cloudy::A1
    Hvar_rain_dt::A1
    QTvar_rain_dt::A1
    HQTcov_rain_dt::A1
    prec_source_qt::A1
    prec_source_h::A1
    function EnvironmentThermodynamics(
        namelist,
        Gr::Grid,
        Ref::ReferenceState,
        EnvVar::EnvironmentVariables,
        Rain::RainVariables,
    )
        quadrature_order = try
            namelist["thermodynamics"]["quadrature_order"]
        catch
            3
        end
        quadrature_type = try
            namelist["thermodynamics"]["quadrature_type"]
        catch
            "gaussian"
        end
        t_to_prog_fp = t_to_thetali_c
        prog_to_t_fp = eos_first_guess_thetal

        qt_dry = center_field(Gr)
        th_dry = center_field(Gr)

        t_cloudy = center_field(Gr)
        qv_cloudy = center_field(Gr)
        qt_cloudy = center_field(Gr)
        th_cloudy = center_field(Gr)

        Hvar_rain_dt = center_field(Gr)
        QTvar_rain_dt = center_field(Gr)
        HQTcov_rain_dt = center_field(Gr)

        prec_source_qt = center_field(Gr)
        prec_source_h = center_field(Gr)
        A1 = typeof(qt_dry)
        return new{A1}(
            Gr,
            Ref,
            quadrature_order,
            quadrature_type,
            t_to_prog_fp,
            prog_to_t_fp,
            qt_dry,
            th_dry,
            t_cloudy,
            qv_cloudy,
            qt_cloudy,
            th_cloudy,
            Hvar_rain_dt,
            QTvar_rain_dt,
            HQTcov_rain_dt,
            prec_source_qt,
            prec_source_h,
        )
    end
end

mutable struct ParameterizationBase{T}
    turbulence_tendency::T
    Gr::Grid
    Ref::ReferenceState
    KM::VariableDiagnostic
    KH::VariableDiagnostic
    prandtl_number::Float64
    Ri_bulk_crit::Float64
    zi::Float64
    # A base class common to all turbulence parameterizations
    function ParameterizationBase(namelist, Gr::Grid, Ref::ReferenceState)
        turbulence_tendency = center_field(Gr)
        KM = VariableDiagnostic(Gr, "half", "scalar", "sym", "diffusivity", "m^2/s") # eddy viscosity
        KH = VariableDiagnostic(Gr, "half", "scalar", "sym", "viscosity", "m^2/s") # eddy diffusivity
        # get values from namelist
        prandtl_number = namelist["turbulence"]["prandtl_number_0"]
        Ri_bulk_crit = namelist["turbulence"]["Ri_bulk_crit"]

        return new{typeof(turbulence_tendency)}(turbulence_tendency, Gr, Ref, KM, KH, prandtl_number, Ri_bulk_crit, 0)
    end
end

# SurfaceMoninObukhovDry:
#     Needed for dry cases (qt=0). They have to be initialized with nonzero qtg for the
#     reference profiles. This surface subroutine sets the latent heat flux to zero
#     to prevent errors due to nonzero qtg in vars such as the obukhov_length.
# SurfaceSullivanPatton
#     Not fully implemented yet. Maybe not needed - Ignacio
struct SurfaceNone end
struct SurfaceFixedFlux end
struct SurfaceFixedCoeffs end
struct SurfaceMoninObukhov end
struct SurfaceMoninObukhovDry end
struct SurfaceSullivanPatton end

Base.@kwdef mutable struct SurfaceBase{T}
    zrough::Float64 = 0
    interactive_zrough::Bool = false
    Tsurface::Float64 = 0
    qsurface::Float64 = 0
    shf::Float64 = 0
    lhf::Float64 = 0
    cm::Float64 = 0
    ch::Float64 = 0
    cq::Float64 = 0
    bflux::Float64 = 0
    windspeed::Float64 = 0
    ustar::Float64 = 0
    rho_qtflux::Float64 = 0
    rho_hflux::Float64 = 0
    rho_uflux::Float64 = 0
    rho_vflux::Float64 = 0
    obukhov_length::Float64 = 0
    Ri_bulk_crit::Float64 = 0
    s_surface::Float64 = 0
    ustar_fixed::Bool = false
    Gr::Grid
    Ref::ReferenceState
end

function SurfaceBase(::Type{T}; Gr::Grid, Ref::ReferenceState, namelist::Dict) where {T}
    Ri_bulk_crit = namelist["turbulence"]["Ri_bulk_crit"]
    return SurfaceBase{T}(; Gr, Ref, Ri_bulk_crit)
end


struct RainPhysics{T}
    Gr::Grid
    Ref::ReferenceState
    rain_evap_source_h::T
    rain_evap_source_qt::T
    function RainPhysics(Gr::Grid, Ref::ReferenceState)
        rain_evap_source_h = center_field(Gr)
        rain_evap_source_qt = center_field(Gr)
        return new{typeof(rain_evap_source_h)}(Gr, Ref, rain_evap_source_h, rain_evap_source_qt)
    end
end

struct ForcingBaseType end
struct ForcingNone end
struct ForcingStandard end
struct ForcingDYCOMS_RF01 end

struct RadiationBaseType end
struct RadiationNone end
struct RadiationStandard end
struct RadiationDYCOMS_RF01 end

Base.@kwdef mutable struct ForcingBase{T}
    subsidence::AbstractArray{Float64, 1} = zeros(1)
    dTdt::AbstractArray{Float64, 1} = zeros(1) # horizontal advection temperature tendency
    dqtdt::AbstractArray{Float64, 1} = zeros(1) # horizontal advection moisture tendency
    apply_coriolis::Bool = false
    apply_subsidence::Bool = false
    coriolis_param::Float64 = 0
    ug::AbstractArray{Float64, 1} = zeros(1)
    vg::AbstractArray{Float64, 1} = zeros(1)
    # (*convert_forcing_prog_fp)(p0, qt, qv, T,::Float64
    #                                   qt_tendency, T_tendency) ::Float64
    convert_forcing_prog_fp::Function = x -> x
    Gr::Grid
    Ref::ReferenceState
end

Base.@kwdef mutable struct RadiationBase{T}
    dTdt::AbstractArray{Float64, 1} = zeros(1) # horizontal advection temperature tendency
    dqtdt::AbstractArray{Float64, 1} = zeros(1) # horizontal advection moisture tendency
    convert_forcing_prog_fp::Function = x -> x
    Gr::Grid
    Ref::ReferenceState
    divergence::Float64 = 0
    alpha_z::Float64 = 0
    kappa::Float64 = 0
    F0::Float64 = 0
    F1::Float64 = 0
    f_rad::AbstractArray{Float64, 1} = zeros(1)
end

Base.@kwdef mutable struct CasesBase{T}
    casename::String = "default_casename"
    inversion_option::String = "default_inversion_option"
    Sur::SurfaceBase
    Fo::ForcingBase
    Rad::RadiationBase
    rad_time::StepRangeLen = linspace(10, 360; num = 36) .* 60
    rad::AbstractMatrix{Float64} = zeros(1, 1)
    lhf0::Float64 = 0
    shf0::Float64 = 0
end

struct SimilarityED{PS}
    param_set::PS
    base::ParameterizationBase
    extrapolate_buoyancy::Bool
end

mutable struct EDMF_PrognosticTKE{PS, A1, A2}
    param_set::PS
    base::ParameterizationBase
    n_updrafts::Int
    calc_tke::Bool
    use_const_plume_spacing::Bool
    calc_scalar_var::Bool
    entr_detr_fp::Function
    pressure_func_buoy::Function
    drag_sign::Int
    pressure_func_drag::Function
    asp_label
    similarity_diffusivity
    extrapolate_buoyancy::Bool
    mixing_scheme::String
    surface_area::Float64
    max_area::Float64
    entrainment_factor::Float64
    entrainment_Mdiv_factor::Float64
    updraft_mixing_frac::Float64
    entrainment_sigma::Float64
    entrainment_smin_tke_coeff::Float64
    entrainment_ed_mf_sigma::Float64
    entrainment_scale::Float64
    constant_plume_spacing::Float64
    detrainment_factor::Float64
    sorting_power::Float64
    turbulent_entrainment_factor::Float64
    pressure_buoy_coeff::Float64
    aspect_ratio::Float64
    pressure_normalmode_buoy_coeff1::Float64
    pressure_normalmode_buoy_coeff2::Float64
    pressure_normalmode_adv_coeff::Float64
    pressure_normalmode_drag_coeff::Float64
    vel_buoy_coeff::Float64
    tke_ed_coeff::Float64
    tke_diss_coeff::Float64
    static_stab_coeff::Float64
    lambda_stab::Float64
    minimum_area::Float64
    Rain::RainVariables
    rainphysics::RainPhysics
    UpdVar::UpdraftVariables
    UpdThermo::UpdraftThermodynamics
    EnvVar::EnvironmentVariables
    EnvThermo::EnvironmentThermodynamics
    entr_sc::A2
    press::A2
    detr_sc::A2
    sorting_function::A2
    b_mix::A2
    frac_turb_entr::A2
    frac_turb_entr_full::A2
    turb_entr_W::A2
    turb_entr_H::A2
    turb_entr_QT::A2
    nh_pressure::A2
    nh_pressure_b::A2
    nh_pressure_adv::A2
    nh_pressure_drag::A2
    asp_ratio::A2
    b_coeff::A2
    m::A2
    mixing_length::A1
    horizontal_KM::A2
    horizontal_KH::A2
    tke_transport::A1
    tke_advection::A1
    area_surface_bc::A1
    w_surface_bc::A1
    h_surface_bc::A1
    qt_surface_bc::A1
    pressure_plume_spacing::A1
    massflux_tendency_h::A1
    massflux_tendency_qt::A1
    diffusive_tendency_h::A1
    diffusive_tendency_qt::A1
    massflux_h::A1
    massflux_qt::A1
    diffusive_flux_h::A1
    diffusive_flux_qt::A1
    diffusive_flux_u::A1
    diffusive_flux_v::A1
    massflux_tke::A1
    prandtl_nvec::A1
    prandtl_number::Float64
    mls::A1
    ml_ratio::A1
    l_entdet::A1
    b::A1
    wstar::Float64
    entr_surface_bc::Float64
    detr_surface_bc::Float64
    dt_upd::Float64
    function EDMF_PrognosticTKE(namelist, Gr::Grid, Ref::ReferenceState, param_set::PS) where {PS}
        # Initialize the base parameterization class
        base = ParameterizationBase(namelist, Gr, Ref)

        # Set the number of updrafts (1)
        n_updrafts = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"]
        catch
            println("Turbulence--EDMF_PrognosticTKE: defaulting to single updraft")
            1
        end

        calc_tke = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"]
        catch
            true
        end

        use_const_plume_spacing = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["use_constant_plume_spacing"]
        catch
            false
        end

        calc_scalar_var = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["calc_scalar_var"]
        catch
            false
        end
        if (calc_scalar_var == true && calc_tke == false)
            error("Turbulence--EDMF_PrognosticTKE: >>calculate_tke<< must be set to true when >>calc_scalar_var<< is true (to calculate the mixing length for the variance and covariance calculations")
        end

        entr_detr_fp = try
            if string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "inverse_z"
                entr_detr_inverse_z
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "dry"
                entr_detr_dry
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "inverse_w"
                entr_detr_inverse_w
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "b_w2"
                entr_detr_b_w2
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "entr_detr_tke"
                entr_detr_tke
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "suselj"
                entr_detr_suselj
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "buoyancy_sorting"
                entr_detr_buoyancy_sorting
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "moisture_deficit"
                entr_detr_env_moisture_deficit
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "moisture_deficit_div"
                entr_detr_env_moisture_deficit_div
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]) == "none"
                entr_detr_none
            else
                error("Turbulence--EDMF_PrognosticTKE: Entrainment rate namelist option is not recognized")
            end
        catch
            println("Turbulence--EDMF_PrognosticTKE: defaulting to cloudy entrainment formulation")
            entr_detr_b_w2
        end
        if (calc_tke == false && "tke" in string(namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]))
            error("Turbulence--EDMF_PrognosticTKE: >>calc_tke<< must be set to true when entrainment is using tke")
        end


        pressure_func_buoy = try
            if string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"]) == "tan18"
                pressure_tan18_buoy
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"]) == "normalmode"
                pressure_normalmode_buoy
            else
                error("Turbulence--EDMF_PrognosticTKE: pressure closure in namelist option is not recognized")
            end
        catch
            println("Turbulence--EDMF_PrognosticTKE: defaulting to pressure closure Tan2018")
            pressure_tan18_buoy
        end

        drag_sign = false
        pressure_func_drag = try
            if string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"]) == "tan18"
                pressure_tan18_drag
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"]) == "normalmode"
                pressure_normalmode_drag
            elseif string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"]) == "normalmode_signdf"
                drag_sign = true
                pressure_normalmode_drag
            else
                error("Turbulence--EDMF_PrognosticTKE: pressure closure in namelist option is not recognized")
            end
        catch
            println("Turbulence--EDMF_PrognosticTKE: defaulting to pressure closure Tan2018")
            pressure_tan18_drag
        end

        asp_label = try
            string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_asp_label"])
        catch
            println("Turbulence--EDMF_PrognosticTKE: H/2R defaulting to constant")
            "const"
        end

        similarity_diffusivity = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["use_similarity_diffusivity"]
        catch
            println("Turbulence--EDMF_PrognosticTKE: defaulting to TKE-based eddy diffusivity")
            false
        end
        if (similarity_diffusivity == false && calc_tke == false)
            error("Turbulence--EDMF_PrognosticTKE: either >>use_similarity_diffusivity<< or >>calc_tke<< flag is needed to get the eddy diffusivities")
        end

        if (similarity_diffusivity == true && calc_tke == true)
            println("TKE will be calculated but not used for eddy diffusivity calculation")
        end

        extrapolate_buoyancy = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["extrapolate_buoyancy"]
        catch
            println("Turbulence--EDMF_PrognosticTKE: defaulting to extrapolation of updraft buoyancy along a pseudoadiabat")
            true
        end

        mixing_scheme = try
            string(namelist["turbulence"]["EDMF_PrognosticTKE"]["mixing_length"])
        catch
            println("Using (Tan et al, 2018) default")
            "default"
        end

        # Get values from namelist
        # set defaults at some point?
        surface_area = namelist["turbulence"]["EDMF_PrognosticTKE"]["surface_area"]
        max_area = namelist["turbulence"]["EDMF_PrognosticTKE"]["max_area"]
        entrainment_factor = namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"]
        entrainment_Mdiv_factor = try
            namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"]
        catch
            0.0
        end
        updraft_mixing_frac = namelist["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"]
        entrainment_sigma = namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_sigma"]
        entrainment_smin_tke_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"]
        entrainment_ed_mf_sigma = namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"]
        entrainment_scale = namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"]
        constant_plume_spacing = namelist["turbulence"]["EDMF_PrognosticTKE"]["constant_plume_spacing"]
        detrainment_factor = namelist["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"]
        sorting_power = namelist["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"]
        turbulent_entrainment_factor = namelist["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"]
        pressure_buoy_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_buoy_coeff"]
        aspect_ratio = namelist["turbulence"]["EDMF_PrognosticTKE"]["aspect_ratio"]

        if string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"]) == "normalmode"
            pressure_normalmode_buoy_coeff1, pressure_normalmode_buoy_coeff2 = try
                (
                    namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"],
                    namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"],
                )
            catch
                println("Using (Tan et al, 2018) parameters as default for Normal Mode pressure formula buoyancy term")
                (pressure_buoy_coeff, 0.0)
            end
        end

        if string(namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"]) == "normalmode"
            pressure_normalmode_adv_coeff, pressure_normalmode_drag_coeff = try
                (
                    namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"],
                    namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"],
                )
            catch
                println("Using (Tan et al, 2018) parameters as default for Normal Mode pressure formula drag term")
                (0.0, 1.0)
            end
        end

        # "Legacy" coefficients used by the steady updraft routine
        vel_buoy_coeff = 1.0 - pressure_buoy_coeff
        if calc_tke == true
            tke_ed_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"]
            tke_diss_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"]
            static_stab_coeff = namelist["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"]
            # Latent heat stability effect
            lambda_stab = namelist["turbulence"]["EDMF_PrognosticTKE"]["lambda_stab"]
        end
        # Need to code up as namelist option?
        minimum_area = 1e-5

        # Create the class for rain
        Rain = RainVariables(namelist, Gr)

        # Create the updraft variable class (major diagnostic and prognostic variables)
        UpdVar = UpdraftVariables(n_updrafts, namelist, Gr)
        # Create the class for updraft thermodynamics
        UpdThermo = UpdraftThermodynamics(n_updrafts, Gr, Ref, UpdVar, Rain)

        # Create the environment variable class (major diagnostic and prognostic variables)
        EnvVar = EnvironmentVariables(namelist, Gr, param_set)
        # Create the class for environment thermodynamics
        EnvThermo = EnvironmentThermodynamics(namelist, Gr, Ref, EnvVar, Rain)

        # Entrainment rates
        entr_sc = center_field(Gr, n_updrafts)
        press = center_field(Gr, n_updrafts)

        # Detrainment rates
        detr_sc = center_field(Gr, n_updrafts)

        sorting_function = center_field(Gr, n_updrafts)
        b_mix = center_field(Gr, n_updrafts)

        # turbulent entrainment
        frac_turb_entr = center_field(Gr, n_updrafts)
        frac_turb_entr_full = face_field(Gr, n_updrafts)
        turb_entr_W = center_field(Gr, n_updrafts)
        turb_entr_H = center_field(Gr, n_updrafts)
        turb_entr_QT = center_field(Gr, n_updrafts)

        # Pressure term in updraft vertical momentum equation
        nh_pressure = center_field(Gr, n_updrafts)
        nh_pressure_b = center_field(Gr, n_updrafts)
        nh_pressure_adv = center_field(Gr, n_updrafts)
        nh_pressure_drag = center_field(Gr, n_updrafts)
        asp_ratio = center_field(Gr, n_updrafts)
        b_coeff = center_field(Gr, n_updrafts)

        # Mass flux
        m = center_field(Gr, n_updrafts)

        # mixing length
        mixing_length = center_field(Gr)
        horizontal_KM = center_field(Gr, n_updrafts)
        horizontal_KH = center_field(Gr, n_updrafts)

        # diagnosed tke budget terms
        tke_transport = center_field(Gr)
        tke_advection = center_field(Gr)

        # Near-surface BC of updraft area fraction
        area_surface_bc = pyzeros(n_updrafts)
        w_surface_bc = pyzeros(n_updrafts)
        h_surface_bc = pyzeros(n_updrafts)
        qt_surface_bc = pyzeros(n_updrafts)
        pressure_plume_spacing = pyzeros(n_updrafts)

        # Mass flux tendencies of mean scalars (for output)
        massflux_tendency_h = center_field(Gr)
        massflux_tendency_qt = center_field(Gr)
        rainphysics = RainPhysics(Gr, Ref)

        # (Eddy) diffusive tendencies of mean scalars (for output)
        diffusive_tendency_h = center_field(Gr)
        diffusive_tendency_qt = center_field(Gr)

        # Vertical fluxes for output
        massflux_h = center_field(Gr)
        massflux_qt = center_field(Gr)
        diffusive_flux_h = center_field(Gr)
        diffusive_flux_qt = center_field(Gr)
        diffusive_flux_u = center_field(Gr)
        diffusive_flux_v = center_field(Gr)
        if calc_tke
            massflux_tke = center_field(Gr)
        end

        # Added by Ignacio : Length scheme in use (mls), and smooth min effect (ml_ratio)
        # Variable Prandtl number initialized as neutral value.
        ones_vec = center_field(Gr)
        prandtl_nvec = base.prandtl_number .* ones_vec
        mls = center_field(Gr)
        ml_ratio = center_field(Gr)
        l_entdet = center_field(Gr)
        b = center_field(Gr)
        wstar = 0
        entr_surface_bc = 0
        detr_surface_bc = 0
        dt_upd = 0
        A1 = typeof(mixing_length)
        A2 = typeof(horizontal_KM)
        return new{PS, A1, A2}(
            param_set,
            base,
            n_updrafts,
            calc_tke,
            use_const_plume_spacing,
            calc_scalar_var,
            entr_detr_fp,
            pressure_func_buoy,
            drag_sign,
            pressure_func_drag,
            asp_label,
            similarity_diffusivity,
            extrapolate_buoyancy,
            mixing_scheme,
            surface_area,
            max_area,
            entrainment_factor,
            entrainment_Mdiv_factor,
            updraft_mixing_frac,
            entrainment_sigma,
            entrainment_smin_tke_coeff,
            entrainment_ed_mf_sigma,
            entrainment_scale,
            constant_plume_spacing,
            detrainment_factor,
            sorting_power,
            turbulent_entrainment_factor,
            pressure_buoy_coeff,
            aspect_ratio,
            pressure_normalmode_buoy_coeff1,
            pressure_normalmode_buoy_coeff2,
            pressure_normalmode_adv_coeff,
            pressure_normalmode_drag_coeff,
            vel_buoy_coeff,
            tke_ed_coeff,
            tke_diss_coeff,
            static_stab_coeff,
            lambda_stab,
            minimum_area,
            Rain,
            rainphysics,
            UpdVar,
            UpdThermo,
            EnvVar,
            EnvThermo,
            entr_sc,
            press,
            detr_sc,
            sorting_function,
            b_mix,
            frac_turb_entr,
            frac_turb_entr_full,
            turb_entr_W,
            turb_entr_H,
            turb_entr_QT,
            nh_pressure,
            nh_pressure_b,
            nh_pressure_adv,
            nh_pressure_drag,
            asp_ratio,
            b_coeff,
            m,
            mixing_length,
            horizontal_KM,
            horizontal_KH,
            tke_transport,
            tke_advection,
            area_surface_bc,
            w_surface_bc,
            h_surface_bc,
            qt_surface_bc,
            pressure_plume_spacing,
            massflux_tendency_h,
            massflux_tendency_qt,
            diffusive_tendency_h,
            diffusive_tendency_qt,
            massflux_h,
            massflux_qt,
            diffusive_flux_h,
            diffusive_flux_qt,
            diffusive_flux_u,
            diffusive_flux_v,
            massflux_tke,
            prandtl_nvec,
            base.prandtl_number,
            mls,
            ml_ratio,
            l_entdet,
            b,
            wstar,
            entr_surface_bc,
            detr_surface_bc,
            dt_upd,
        )
    end
end
get_grid(edmf::EDMF_PrognosticTKE) = edmf.base.Gr
reference_state(edmf::EDMF_PrognosticTKE) = edmf.base.Ref
prandtl_number(edmf::EDMF_PrognosticTKE) = edmf.base.prandtl_number
turbulence_tendency(edmf::EDMF_PrognosticTKE) = edmf.base.turbulence_tendency
diffusivity_m(edmf::EDMF_PrognosticTKE) = edmf.base.KM
diffusivity_h(edmf::EDMF_PrognosticTKE) = edmf.base.KH
Ri_bulk_crit(edmf::EDMF_PrognosticTKE) = edmf.base.Ri_bulk_crit
Ri_bulk_crit(base::ParameterizationBase) = base.Ri_bulk_crit
parameter_set(obj) = obj.param_set
