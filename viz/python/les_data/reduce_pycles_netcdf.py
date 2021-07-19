
import netCDF4 as nc
import pylab as plt
import argparse
import subprocess
import os
# command line:
# python reduce_pycles_netcdf.py input   output
def main():
    parser = argparse.ArgumentParser(prog="PyCLES")
    parser.add_argument("fullfilename")
    parser.add_argument("fname")
    args = parser.parse_args()
    fullfilename = args.fullfilename
    fname = args.fname

    data = nc.Dataset(fullfilename, "r")
    buoyancy_mean_  = np.array(data.groups["profiles"].variables["buoyancy_mean"])
    env_u_ = data.groups["profiles"].variables["env_u"]
    env_v_ = data.groups["profiles"].variables["env_v"]
    env_w_ = data.groups["profiles"].variables["env_w"]
    temperature_mean_ = data.groups["profiles"].variables["temperature_mean"]
    u_mean_ = data.groups["profiles"].variables["u_mean"]
    v_mean_ = data.groups["profiles"].variables["v_mean"]
    w_mean_ = data.groups["profiles"].variables["w_mean"]
    w_mean2_ = data.groups["profiles"].variables["w_mean2"]
    w_mean3_ = data.groups["profiles"].variables["w_mean3"]
    tke_mean_ = data.groups["profiles"].variables["tke_mean"]
    v_translational_mean_ = data.groups["profiles"].variables["v_translational_mean"]
    u_translational_mean_ = data.groups["profiles"].variables["u_translational_mean"]
    updraft_buoyancy_ = data.groups["profiles"].variables["updraft_b"]
    updraft_fraction_ = data.groups["profiles"].variables["updraft_fraction"]

    total_flux_u_ = data.groups["profiles"].variables["resolved_x_vel_flux"]
    total_flux_v_ = data.groups["profiles"].variables["resolved_y_vel_flux"]

    updraft_ddz_p_alpha_ = data.groups["profiles"].variables["updraft_ddz_p_alpha"]

    rho_ = data.groups["reference"].variables["rho0_half"]
    p0_ = data.groups["reference"].variables["p0_half"]

    t = np.array(data.groups["profiles"].variables["t"])
    z = np.array(data.groups["profiles"].variables["z_half"])

    # try the TKE diagnostics outputs
    tke_prod_A_ = data.groups["profiles"].variables["tke_prod_A"]
    tke_prod_B_ = data.groups["profiles"].variables["tke_prod_B"]
    tke_prod_D_ = data.groups["profiles"].variables["tke_prod_D"]
    tke_prod_P_ = data.groups["profiles"].variables["tke_prod_P"]
    tke_prod_T_ = data.groups["profiles"].variables["tke_prod_T"]
    tke_prod_S_ = data.groups["profiles"].variables["tke_prod_S"]
    tke_nd_mean_ = data.groups["profiles"].variables["tke_nd_mean"]
    try:
        cloud_fraction_ = data.groups["profiles"].variables["cloud_fraction"]
    except:
        cloud_fraction_ = np.zeros_like(env_w_)
    try:
        dqtdt_precip_mean_ = data.groups["profiles"].variables["dqtdt_precip_mean"]
    except:
        dqtdt_precip_mean_ = np.zeros_like(env_w_)

    try:
        resolved_z_flux_thetali_ = data.groups["profiles"].variables["resolved_z_flux_thetali"]
        resolved_z_flux_qt_ = data.groups["profiles"].variables["resolved_z_flux_qt"]
        sgs_z_flux_thetali_ = data.groups["profiles"].variables["sgs_z_flux_thetali"]
        sgs_z_flux_qt_ = data.groups["profiles"].variables["sgs_z_flux_qt"]
    except:
        resolved_z_flux_thetali_ = data.groups["profiles"].variables["resolved_z_flux_theta"]
        sgs_z_flux_thetali_ = data.groups["profiles"].variables["sgs_z_flux_theta"]
        sgs_z_flux_qt_ = np.zeros_like(env_w_)
        resolved_z_flux_qt_ = np.zeros_like(env_w_)
    try:
        qr_mean_ = data.groups["profiles"].variables["qr_mean"]
        env_qr_ =  data.groups["profiles"].variables["env_qr"]
        updraft_qr_ = data.groups["profiles"].variables["updraft_qr"]
    except:
        qr_mean_ = np.zeros_like(env_w_)
        env_qr_ = np.zeros_like(env_w_)
        updraft_qr_ = np.zeros_like(env_w_)
    try:
        qi_mean_ = data.groups["profiles"].variables["qi_mean"]
        env_qi_ = data.groups["profiles"].variables["env_qi"]
        updraft_qi_ = data.groups["profiles"].variables["updraft_qi"]
    except:
        qi_mean_ = np.zeros_like(env_w_)
        env_qi_ = np.zeros_like(env_w_)
        updraft_qi_ = np.zeros_like(env_w_)
    try:
        qt_mean_ = data.groups["profiles"].variables["qt_mean"]
        qt_mean2_ = data.groups["profiles"].variables["qt_mean2"]
        qt_mean3_ = data.groups["profiles"].variables["qt_mean3"]
        env_qt_ = data.groups["profiles"].variables["env_qt"]
        env_qt2_ = data.groups["profiles"].variables["env_qt2"]
        updraft_qt_ = data.groups["profiles"].variables["updraft_qt"]
    except:
        qt_mean_ = np.zeros_like(env_w_)
        qt_mean2_ = np.zeros_like(env_w_)
        qt_mean3_ = np.zeros_like(env_w_)
        env_qt_ = np.zeros_like(env_w_)
        env_qt2_ = np.zeros_like(env_w_)
        updraft_qt_ = np.zeros_like(env_w_)
    try:
        ql_mean_ = data.groups["profiles"].variables["ql_mean"]
        env_ql_ = data.groups["profiles"].variables["env_ql"]
        env_qv_ = data.groups["profiles"].variables["env_qv"]
        updraft_ql_ = data.groups["profiles"].variables["updraft_ql"]
        updraft_qv_ = data.groups["profiles"].variables["updraft_qv"]
    except:
        ql_mean_ = np.zeros_like(env_w_)
        env_ql_ = np.zeros_like(env_w_)
        env_qv_ = np.zeros_like(env_w_)
        updraft_ql_ = np.zeros_like(env_w_)
        updraft_qv_ = np.zeros_like(env_w_)
    try:
        env_qt_thetali_ = data.groups["profiles"].variables["env_qt_thetali"]
    except:
        env_qt_thetali_ = np.zeros_like(env_w_)

    try:
        thetali_mean_ = data.groups["profiles"].variables["thetali_mean"]
        thetali_mean2_ = data.groups["profiles"].variables["thetali_mean2"]
        thetali_mean3_ = data.groups["profiles"].variables["thetali_mean3"]
    except:
        thetali_mean_ = data.groups["profiles"].variables["theta_mean"]
        thetali_mean2_ = data.groups["profiles"].variables["theta_mean2"]
        thetali_mean3_ = data.groups["profiles"].variables["theta_mean3"]

    env_thetali_ = data.groups["profiles"].variables["env_thetali"]
    env_thetali2_ = data.groups["profiles"].variables["env_thetali2"]
    env_buoyancy_ = data.groups["profiles"].variables["env_b"]
    updraft_thetali_ = data.groups["profiles"].variables["updraft_thetali"]
    updraft_u_ = data.groups["profiles"].variables["updraft_u"]
    updraft_v_ = data.groups["profiles"].variables["updraft_v"]
    updraft_w_ = data.groups["profiles"].variables["updraft_w"]

    friction_velocity_mean_ = data.groups["timeseries"].variables["friction_velocity_mean"]
    shf_surface_mean_ = data.groups["timeseries"].variables["shf_surface_mean"]
    lhf_surface_mean_ = data.groups["timeseries"].variables["lhf_surface_mean"]
    try:
        cloud_fraction_mean_ = data.groups["timeseries"].variables["cloud_fraction"]
        cloud_base_ = data.groups["timeseries"].variables["cloud_base"]
        cloud_top_ = data.groups["timeseries"].variables["cloud_top"]
    except:
        cloud_fraction_mean_ = np.zeros_like(lhf_surface_mean_)
        cloud_base_ = np.zeros_like(lhf_surface_mean_)
        cloud_top_ = np.zeros_like(lhf_surface_mean_)
    try:
        lwp_ = data.groups["timeseries"].variables["lwp"]
    except:
    	lwp_ = np.zeros_like(lhf_surface_mean_)
    try:
        updraft_temperature_ = data.groups["profiles"].variables["updraft_temperature"]
        env_temperature_ = data.groups["profiles"].variables["env_temperature"]
    except:
        updraft_temperature_ = np.zeros_like(thetali_mean_)
        env_temperature_ = np.zeros_like(thetali_mean_)

    # thetali_srf_int_ = data.groups["timeseries"].variables["thetali_srf_int"] # this is here since

    z_half_ = data.groups["profiles"].variables["z_half"]
    t_ = data.groups["profiles"].variables["t"]


    # flux diagnosis (diffusive_flux = total_flux - massflux, diffusive_flux is weighted by env_fraction as in TurbulenceConvection TubrbProgTKE line 1723)
    rho_temp = np.tile(rho_,(np.shape(updraft_fraction_)[0],1))
    p0_temp = np.tile(p0_,(np.shape(updraft_fraction_)[0],1))
    a_1_a = np.multiply(rho_temp, np.multiply(updraft_fraction_, np.subtract(1.0,updraft_fraction_)))
    updraft_buoyancy_ -=buoyancy_mean_
    env_buoyancy_ -=buoyancy_mean_
    Wvar_mean_ = calc_covar(w_mean2_, w_mean_, w_mean_)
    Hvar_mean_ = calc_covar(thetali_mean2_, thetali_mean_, thetali_mean_)
    QTvar_mean_ = calc_covar(qt_mean2_,       qt_mean_,      qt_mean_)
    env_Hvar_ = calc_covar(env_thetali2_,   env_thetali_,  env_thetali_)
    env_QTvar_ = calc_covar(env_qt2_,        env_qt_,       env_qt_)
    env_HQTcov_ = calc_covar(env_qt_thetali_, env_qt_,       env_thetali_)
    massflux_          = np.multiply(updraft_fraction_,np.multiply(rho_,updraft_w_))
    massflux_h_        = np.multiply(a_1_a,np.multiply(np.subtract(updraft_w_, env_w_), np.subtract(updraft_thetali_, env_thetali_)))
    massflux_qt_       = np.multiply(a_1_a,np.multiply(np.subtract(updraft_w_, env_w_), np.subtract(updraft_qt_, env_qt_)))
    massflux_u_        = np.multiply(a_1_a,np.multiply(np.subtract(updraft_w_, env_w_), np.subtract(updraft_u_, env_u_)))
    massflux_v_        = np.multiply(a_1_a,np.multiply(np.subtract(updraft_w_, env_w_), np.subtract(updraft_v_, env_v_)))
    total_flux_h_      = np.add(resolved_z_flux_thetali_[:, :], sgs_z_flux_thetali_[:, :])
    total_flux_qt_     = np.add(resolved_z_flux_qt_[:, :], sgs_z_flux_qt_[:, :])
    massflux_          = np.multiply(rho_temp, massflux_)
    massflux_h_        = np.multiply(rho_temp, massflux_h_)
    massflux_qt_       = np.multiply(rho_temp, massflux_qt_) # already multiplied qt by 1000 above
    massflux_u_        = np.multiply(rho_temp, massflux_u_)
    massflux_v_        = np.multiply(rho_temp, massflux_v_)
    total_flux_h_      = np.multiply(rho_temp, total_flux_h_)
    total_flux_qt_     = np.multiply(rho_temp, total_flux_qt_)
    total_flux_u_      = np.multiply(rho_temp, total_flux_u_)
    total_flux_v_      = np.multiply(rho_temp, total_flux_v_)
    diffusive_flux_h_  = np.subtract(total_flux_h_,massflux_h_)
    diffusive_flux_qt_ = np.subtract(total_flux_qt_,massflux_qt_)
    diffusive_flux_u_  = np.subtract(total_flux_u_,massflux_u_)
    diffusive_flux_v_  = np.subtract(total_flux_v_,massflux_v_)

    H_third_m_ = calc_third_m(thetali_mean_, thetali_mean2_, thetali_mean3_, Hvar_mean_, updraft_fraction_)
    QT_third_m_ = calc_third_m(qt_mean_, qt_mean2_, qt_mean3_, QTvar_mean_, updraft_fraction_)
    W_third_m_ = calc_third_m(w_mean_, w_mean2_, w_mean3_, Wvar_mean_, updraft_fraction_)

    updraft_RH_ = np.zeros_like(env_qt_)
    env_RH_ = np.zeros_like(env_qt_)

    updraft_RH_ = relative_humidity(p0_temp, updraft_qt_, updraft_qv_, updraft_temperature_)
    env_RH_     = relative_humidity(p0_temp, env_qt_, env_qv_, env_temperature_)

    output = nc.Dataset(fname, "w", format="NETCDF4")
    output.createDimension("z", len(z_half_))
    output.createDimension("t", len(t_))
    output.createDimension("dim", None)
    output.createGroup("profiles")
    output.createGroup("timeseries")

    t = output.createVariable("t", "f4", "t")
    z_half = output.createVariable("z_half", "f4", "z")

    profiles_grp = output.groups["profiles"]

    rho = profiles_grp.createVariable("rho","f4",("z"))
    p0 = profiles_grp.createVariable("p0","f4",("z"))
    Hvar_mean = profiles_grp.createVariable("Hvar_mean","f4",("t","z"))
    QTvar_mean = profiles_grp.createVariable("QTvar_mean","f4",("t","z"))
    env_Hvar = profiles_grp.createVariable("env_Hvar","f4",("t","z"))
    env_QTvar = profiles_grp.createVariable("env_QTvar","f4",("t","z"))
    env_HQTcov = profiles_grp.createVariable("env_HQTcov","f4",("t","z"))
    W_third_m = profiles_grp.createVariable("W_third_m","f4",("t","z"))
    H_third_m = profiles_grp.createVariable("H_third_m","f4",("t","z"))
    QT_third_m = profiles_grp.createVariable("QT_third_m","f4",("t","z"))
    massflux = profiles_grp.createVariable("massflux","f4",("t","z"))
    massflux_h = profiles_grp.createVariable("massflux_h","f4",("t","z"))
    massflux_qt = profiles_grp.createVariable("massflux_qt","f4",("t","z"))
    total_flux_h = profiles_grp.createVariable("total_flux_h","f4",("t","z"))
    total_flux_qt = profiles_grp.createVariable("total_flux_qt","f4",("t","z"))
    diffusive_flux_h = profiles_grp.createVariable("diffusive_flux_h","f4",("t","z"))
    diffusive_flux_qt = profiles_grp.createVariable("diffusive_flux_qt","f4",("t","z"))
    total_flux_u = profiles_grp.createVariable("total_flux_u","f4",("t","z"))
    total_flux_v = profiles_grp.createVariable("total_flux_v","f4",("t","z"))
    diffusive_flux_u = profiles_grp.createVariable("diffusive_flux_u","f4",("t","z"))
    diffusive_flux_v = profiles_grp.createVariable("diffusive_flux_v","f4",("t","z"))
    massflux_u = profiles_grp.createVariable("massflux_u","f4",("t","z"))
    massflux_v = profiles_grp.createVariable("massflux_v","f4",("t","z"))
    buoyancy_mean = profiles_grp.createVariable("buoyancy_mean","f4",("t","z"))
    cloud_fraction = profiles_grp.createVariable("cloud_fraction","f4",("t","z"))
    dqtdt_precip_mean = profiles_grp.createVariable("dqtdt_precip_mean","f4",("t","z"))
    # resolved_z_flux_thetali = profiles_grp.createVariable("resolved_z_flux_thetali","f4",("t","z"))
    # resolved_z_flux_qt = profiles_grp.createVariable("resolved_z_flux_qt","f4",("t","z"))
    temperature_mean = profiles_grp.createVariable("temperature_mean","f4",("t","z"))
    updraft_ddz_p_alpha = profiles_grp.createVariable("updraft_ddz_p_alpha","f4",("t","z"))
    thetali_mean = profiles_grp.createVariable("thetali_mean","f4",("t","z"))
    qt_mean = profiles_grp.createVariable("qt_mean","f4",("t","z"))
    ql_mean = profiles_grp.createVariable("ql_mean","f4",("t","z"))
    u_mean = profiles_grp.createVariable("u_mean","f4",("t","z"))
    v_mean = profiles_grp.createVariable("v_mean","f4",("t","z"))
    tke_mean = profiles_grp.createVariable("tke_mean","f4",("t","z"))
    v_translational_mean = profiles_grp.createVariable("v_translational_mean","f4",("t","z"))
    u_translational_mean = profiles_grp.createVariable("u_translational_mean","f4",("t","z"))
    updraft_buoyancy = profiles_grp.createVariable("updraft_buoyancy","f4",("t","z"))
    updraft_fraction = profiles_grp.createVariable("updraft_fraction","f4",("t","z"))
    env_thetali = profiles_grp.createVariable("env_thetali","f4",("t","z"))
    updraft_thetali = profiles_grp.createVariable("updraft_thetali","f4",("t","z"))
    env_qt = profiles_grp.createVariable("env_qt","f4",("t","z"))
    updraft_qt = profiles_grp.createVariable("updraft_qt","f4",("t","z"))
    env_ql = profiles_grp.createVariable("env_ql","f4",("t","z"))
    env_buoyancy = profiles_grp.createVariable("env_buoyancy","f4",("t","z"))
    updraft_ql = profiles_grp.createVariable("updraft_ql","f4",("t","z"))
    qr_mean = profiles_grp.createVariable("qr_mean","f4",("t","z"))
    env_qr = profiles_grp.createVariable("env_qr","f4",("t","z"))
    updraft_qr = profiles_grp.createVariable("updraft_qr","f4",("t","z"))
    updraft_w = profiles_grp.createVariable("updraft_w","f4",("t","z"))
    updraft_RH = profiles_grp.createVariable("updraft_RH","f4",("t","z"))
    env_RH     = profiles_grp.createVariable("env_RH","f4",("t","z"))
    env_w = profiles_grp.createVariable("env_w","f4",("t","z"))
    thetali_mean2 = profiles_grp.createVariable("thetali_mean2","f4",("t","z"))
    qt_mean2 = profiles_grp.createVariable("qt_mean2","f4",("t","z"))
    env_thetali2 = profiles_grp.createVariable("env_thetali2","f4",("t","z"))
    env_qt2 = profiles_grp.createVariable("env_qt2","f4",("t","z"))
    env_qt_thetali = profiles_grp.createVariable("env_qt_thetali","f4",("t","z"))
    tke_prod_A = profiles_grp.createVariable("tke_prod_A","f4",("t","z"))
    tke_prod_B = profiles_grp.createVariable("tke_prod_B","f4",("t","z"))
    tke_prod_D = profiles_grp.createVariable("tke_prod_D","f4",("t","z"))
    tke_prod_P = profiles_grp.createVariable("tke_prod_P","f4",("t","z"))
    tke_prod_T = profiles_grp.createVariable("tke_prod_T","f4",("t","z"))
    tke_prod_S = profiles_grp.createVariable("tke_prod_S","f4",("t","z"))
    tke_nd_mean = profiles_grp.createVariable("tke_nd_mean","f4",("t","z"))

    timeseries_grp = output.groups["timeseries"]
    cloud_fraction_mean = timeseries_grp.createVariable("cloud_fraction_mean","f4","t")
    cloud_base_mean = timeseries_grp.createVariable("cloud_base_mean","f4","t")
    cloud_top_mean = timeseries_grp.createVariable("cloud_top_mean","f4","t")
    friction_velocity_mean = timeseries_grp.createVariable("friction_velocity_mean","f4","t")
    shf_surface_mean = timeseries_grp.createVariable("shf_surface_mean","f4","t")
    lhf_surface_mean = timeseries_grp.createVariable("lhf_surface_mean","f4","t")
    lwp_mean = timeseries_grp.createVariable("lwp_mean","f4","t")
    # thetali_srf_int = timeseries_grp.createVariable("thetali_srf_int","f4","t")

    rho[:] = rho_[:]
    p0[:] = p0_[:]
    Hvar_mean[:,:] = Hvar_mean_[:,:]
    QTvar_mean[:,:] = QTvar_mean_[:,:]
    W_third_m[:,:] = W_third_m_[:,:]
    H_third_m[:,:] = H_third_m_[:,:]
    QT_third_m[:,:] = QT_third_m_[:,:]
    env_Hvar[:,:] = env_Hvar_[:,:]
    env_QTvar[:,:] = env_QTvar_[:,:]
    env_HQTcov[:,:] = env_HQTcov_[:,:]
    massflux[:,:] = massflux_[:,:]
    massflux_h[:,:] = massflux_h_[:,:]
    massflux_qt[:,:] = massflux_qt_[:,:]
    massflux_u[:,:] = massflux_u_[:,:]
    massflux_v[:,:] = massflux_v_[:,:]
    total_flux_u[:,:] = total_flux_u_[:,:]
    total_flux_v[:,:] = total_flux_v_[:,:]
    total_flux_h[:,:] = total_flux_h_[:,:]
    total_flux_qt[:,:] = total_flux_qt_[:,:]
    diffusive_flux_h[:,:] = diffusive_flux_h_[:,:]
    diffusive_flux_qt[:,:] = diffusive_flux_qt_[:,:]
    diffusive_flux_u[:,:] = diffusive_flux_u_[:,:]
    diffusive_flux_v[:,:] = diffusive_flux_v_[:,:]
    buoyancy_mean[:,:] = buoyancy_mean_[:,:]
    # resolved_z_flux_thetali[:,:] = resolved_z_flux_thetali_[:,:]
    # resolved_z_flux_qt[:,:] = resolved_z_flux_qt_[:,:]
    temperature_mean[:,:] = temperature_mean_[:,:]
    updraft_ddz_p_alpha[:,:] = updraft_ddz_p_alpha_[:,:]
    thetali_mean[:,:] = thetali_mean_[:,:]
    qt_mean[:,:] = qt_mean_[:,:]
    ql_mean[:,:] = ql_mean_[:,:]
    u_mean[:,:] = u_mean_[:,:]
    v_mean[:,:] = v_mean_[:,:]
    tke_mean[:,:] = tke_mean_[:,:]
    v_translational_mean[:,:] = v_translational_mean_[:,:]
    u_translational_mean[:,:] = u_translational_mean_[:,:]
    updraft_buoyancy[:,:] = updraft_buoyancy_[:,:]
    updraft_RH[:,:] = updraft_RH_[:,:]
    env_RH[:,:] = env_RH_[:,:]
    updraft_fraction[:,:] = updraft_fraction_[:,:]
    env_thetali[:,:] = env_thetali_[:,:]
    updraft_thetali[:,:] = updraft_thetali_[:,:]
    env_qt[:,:] = env_qt_[:,:]
    updraft_qt[:,:] = updraft_qt_[:,:]
    env_ql[:,:] = env_ql_[:,:]
    updraft_ql[:,:] = updraft_ql_[:,:]
    qr_mean[:,:] = qr_mean_[:,:]
    env_qr[:,:] = env_qr_[:,:]
    updraft_qr[:,:] = updraft_qr_[:,:]
    updraft_w[:,:] = updraft_w_[:,:]
    env_w[:,:] = env_w_[:,:]
    thetali_mean2[:,:] = thetali_mean2_[:,:]
    qt_mean2[:,:] = qt_mean2_[:,:]
    env_thetali2[:,:] = env_thetali2_[:,:]
    env_qt2[:,:] = env_qt2_[:,:]
    env_qt_thetali[:,:] = env_qt_thetali_[:,:]
    env_buoyancy[:,:] = env_buoyancy_[:,:]
    tke_prod_A[:,:] = tke_prod_A_[:,:]
    tke_prod_B[:,:] = tke_prod_B_[:,:]
    tke_prod_D[:,:] = tke_prod_D_[:,:]
    tke_prod_P[:,:] = tke_prod_P_[:,:]
    tke_prod_T[:,:] = tke_prod_T_[:,:]
    tke_prod_S[:,:] = tke_prod_S_[:,:]
    tke_nd_mean[:,:] = tke_nd_mean_[:,:]
    dqtdt_precip_mean[:,:] = dqtdt_precip_mean_[:,:]
    cloud_fraction[:,:] = cloud_fraction_[:,:]

    cloud_fraction_mean[:] = cloud_fraction_mean_[:]
    cloud_base_mean[:] = cloud_base_[:]
    cloud_top_mean[:] = cloud_top_[:]
    friction_velocity_mean[:] = friction_velocity_mean_[:]
    shf_surface_mean[:] = shf_surface_mean_[:]
    lhf_surface_mean[:] = lhf_surface_mean_[:]
    lwp_mean[:] = lwp_[:]

    z_half[:] = z_half_[:]
    t[:] = t_[:]
    output.close()

def calc_covar(var_sq, var1, var2):

    covar = np.subtract(var_sq,np.multiply(var1,var2))
    return covar

def calc_third_m(var, var2, var3, covar, upd_frac):
    A = np.multiply(3.0,np.multiply(var,covar))
    B = np.power(var,3.0)
    C = np.add(A,B)
    third_m = np.subtract(var3,C)
    return third_m

def relative_humidity(p0, qt, qv, T):
    eps_vi   = 1.60745384883
    Tc       = np.subtract(T, 273.15)
    pv       = np.multiply(p0, np.divide(np.multiply(eps_vi , qv) ,(np.subtract(np.add(1.0,np.multiply(eps_vi , qv)),qt))))
    pv_star_ = np.multiply(100.0, 6.1094*np.exp((17.625*np.divide(Tc,np.add(Tc,243.04)))))
    RH       = np.multiply(100.0, np.divide(pv,pv_star_))
    return RH

if __name__ == "__main__":
    main()

