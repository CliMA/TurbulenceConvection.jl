# Entrainment Rates
function entr_detr_dry(entr_in::entr_struct)
    _ret = entr_struct()
    eps = 1.0 # to avoid division by zero when z = 0 or z_i
    # Following Soares 2004
    _ret.entr_sc = 0.5*(1.0/entr_in.z + 1.0/fmax(entr_in.zi - entr_in.z, 10.0)) #vkb/(z + 1.0e-3)
    _ret.detr_sc = 0.0

    return  _ret
end

function entr_detr_inverse_z(entr_in::entr_struct)
    _ret = entr_struct()
    _ret.entr_sc = vkb/entr_in.z
    _ret.detr_sc= 0.0
    return _ret
end

function entr_detr_inverse_w(entr_in::entr_struct)
    _ret = entr_struct()

    eps_w = 1.0/(fmax(fabs(entr_in.w_upd),1.0)* 1000)
    if entr_in.a_upd>0.0
        sorting_function  = buoyancy_sorting(entr_in)
        _ret.entr_sc = sorting_function*eps_w/2.0
        _ret.detr_sc = (1.0-sorting_function/2.0)*eps_w
    else
        _ret.entr_sc = 0.0
        _ret.detr_sc = 0.0
    end
    return _ret
end

function entr_detr_env_moisture_deficit_b_ED_MF(entr_in::entr_struct)
    _ret = entr_struct()

    moisture_deficit_d = (fmax((entr_in.RH_upd/100.0)^entr_in.sort_pow-(entr_in.RH_env/100.0)^entr_in.sort_pow,0.0))^(1.0/entr_in.sort_pow)
    moisture_deficit_e = (fmax((entr_in.RH_env/100.0)^entr_in.sort_pow-(entr_in.RH_upd/100.0)^entr_in.sort_pow,0.0))^(1.0/entr_in.sort_pow)
    _ret.sorting_function = moisture_deficit_e
    c_det = entr_in.c_det
    if (entr_in.ql_up+entr_in.ql_env) ≈ 0.0
        c_det = 0.0
    end

    dw   = entr_in.w_upd - entr_in.w_env
    if dw < 0.0
        dw -= 0.001
    else
        dw += 0.001
    end

    db = (entr_in.b_upd - entr_in.b_env)
    mu = entr_in.c_mu/entr_in.c_mu0

    inv_timescale = fabs(db/dw)
    logistic_e = 1.0/(1.0+exp(-mu*db/dw*(entr_in.chi_upd - entr_in.a_upd/(entr_in.a_upd+entr_in.a_env))))
    logistic_d = 1.0/(1.0+exp( mu*db/dw*(entr_in.chi_upd - entr_in.a_upd/(entr_in.a_upd+entr_in.a_env))))

    #Logistic of buoyancy fluxes
    inv_timescale = fabs(db/dw)
    ed_mf_ratio = fabs(entr_in.buoy_ed_flux)/(fabs(entr_in.a_upd*entr_in.a_env*(entr_in.w_upd-entr_in.w_env)*(entr_in.b_upd - entr_in.b_env))+1e-8)
    logistic_e *= (1.0/(1.0+exp(entr_in.c_ed_mf*(ed_mf_ratio-1.0))))
    _ret.entr_sc = inv_timescale/dw*(entr_in.c_ent*logistic_e + c_det*moisture_deficit_e)
    _ret.detr_sc = inv_timescale/dw*(entr_in.c_ent*logistic_d + c_det*moisture_deficit_d)

    return _ret
end

function entr_detr_env_moisture_deficit(entr_in::entr_in_struct)
    _ret = entr_struct()
    l = pyzeros(2)

    moisture_deficit_d = (fmax((entr_in.RH_upd/100.0)^entr_in.sort_pow-(entr_in.RH_env/100.0)^entr_in.sort_pow,0.0))^(1.0/entr_in.sort_pow)
    moisture_deficit_e = (fmax((entr_in.RH_env/100.0)^entr_in.sort_pow-(entr_in.RH_upd/100.0)^entr_in.sort_pow,0.0))^(1.0/entr_in.sort_pow)
    _ret.sorting_function = moisture_deficit_e
    c_det = entr_in.c_det
    if (entr_in.ql_up+entr_in.ql_env)==0.0
        c_det = 0.0
    end

    dw   = entr_in.w_upd - entr_in.w_env
    if dw < 0.0
        dw -= 0.001
    else
        dw += 0.001
    end

    db = (entr_in.b_upd - entr_in.b_env)
    mu = entr_in.c_mu/entr_in.c_mu0

    inv_timescale = fabs(db/dw)
    logistic_e = 1.0/(1.0+exp(-mu*db/dw*(entr_in.chi_upd - entr_in.a_upd/(entr_in.a_upd+entr_in.a_env))))
    logistic_d = 1.0/(1.0+exp( mu*db/dw*(entr_in.chi_upd - entr_in.a_upd/(entr_in.a_upd+entr_in.a_env))))

    #smooth min
    l[0] = entr_in.tke_coef*fabs(db/sqrt(entr_in.tke+1e-8))
    l[1] = fabs(db/dw)
    inv_timescale = lamb_smooth_minimum(l, 0.1, 0.0005)
    _ret.entr_sc = inv_timescale/dw*(entr_in.c_ent*logistic_e + c_det*moisture_deficit_e)
    _ret.detr_sc = inv_timescale/dw*(entr_in.c_ent*logistic_d + c_det*moisture_deficit_d)

    return _ret
end

function entr_detr_env_moisture_deficit_div(entr_in::entr_struct)
    _ret = entr_struct()
    l = pyzeros(2)

    moisture_deficit_d = (fmax((entr_in.RH_upd/100.0)^entr_in.sort_pow-(entr_in.RH_env/100.0)^entr_in.sort_pow,0.0))^(1.0/entr_in.sort_pow)
    moisture_deficit_e = (fmax((entr_in.RH_env/100.0)^entr_in.sort_pow-(entr_in.RH_upd/100.0)^entr_in.sort_pow,0.0))^(1.0/entr_in.sort_pow)
    _ret.sorting_function = moisture_deficit_e
    c_det = entr_in.c_det
    if (entr_in.ql_up+entr_in.ql_env)==0.0
        c_det = 0.0
    end

    dw   = entr_in.w_upd - entr_in.w_env
    if dw < 0.0
        dw -= 0.001
    else
        dw += 0.001
    end

    db = (entr_in.b_upd - entr_in.b_env)
    mu = entr_in.c_mu/entr_in.c_mu0

    inv_timescale = fabs(db/dw)
    logistic_e = 1.0/(1.0+exp(-mu*db/dw*(entr_in.chi_upd - entr_in.a_upd/(entr_in.a_upd+entr_in.a_env))))
    logistic_d = 1.0/(1.0+exp( mu*db/dw*(entr_in.chi_upd - entr_in.a_upd/(entr_in.a_upd+entr_in.a_env))))

    entr_MdMdz = fmax( entr_in.dMdz/fmax(entr_in.M,1e-12),0.0)
    detr_MdMdz = fmax(-entr_in.dMdz/fmax(entr_in.M,1e-12),0.0)


    #smooth min
    l[0] = entr_in.tke_coef*fabs(db/sqrt(entr_in.tke+1e-8))
    l[1] = fabs(db/dw)
    inv_timescale = lamb_smooth_minimum(l, 0.1, 0.0005)

    _ret.entr_sc = inv_timescale/dw*(entr_in.c_ent*logistic_e + c_det*moisture_deficit_e) + entr_MdMdz * entr_in.c_div
    _ret.detr_sc = inv_timescale/dw*(entr_in.c_ent*logistic_d + c_det*moisture_deficit_d) + detr_MdMdz * entr_in.c_div


    return _ret
end

function entr_detr_buoyancy_sorting(entr_in::entr_struct)
    _ret = entr_struct()

    ret_b = buoyancy_sorting_mean(entr_in)
    b_mix = ret_b.b_mix
    eps_bw2 = entr_in.c_ent*fmax(entr_in.b_upd,0.0) / fmax(entr_in.w_upd * entr_in.w_upd, 1e-2)
    del_bw2 = entr_in.c_ent*fabs(entr_in.b_upd) / fmax(entr_in.w_upd * entr_in.w_upd, 1e-2)
    _ret.b_mix = b_mix
    _ret.sorting_function = ret_b.sorting_function
    _ret.entr_sc = eps_bw2
    if entr_in.ql_up>0.0
        D_ = 0.5*(1.0+entr_in.sort_pow*(ret_b.sorting_function))
        _ret.detr_sc = del_bw2*(1.0+entr_in.c_det*D_)
    else
        _ret.detr_sc = 0.0
    end

    return _ret
end

function buoyancy_sorting_mean(entr_in::entr_struct)

    sorting_function = 0.0
    sa = eos_struct()
    ret_b = buoyant_stract()

    sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, entr_in.qt_env, entr_in.H_env)
    qv_ = entr_in.qt_env - sa.ql
    T_env = sa.T
    ql_env = sa.ql
    rho_env = rho_c(entr_in.p0, sa.T, entr_in.qt_env, qv_)
    b_env = buoyancy_c(entr_in.rho0, rho_env)

    sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, entr_in.qt_up, entr_in.H_up)
    qv_ = entr_in.qt_up - sa.ql
    T_up = sa.T
    ql_up = sa.ql
    rho_up = rho_c(entr_in.p0, sa.T, entr_in.qt_up, qv_)
    b_up = buoyancy_c(entr_in.rho0, rho_up)

    b_mean = entr_in.a_upd*b_up +  (1.0-entr_in.a_upd)*b_env

    # qt_mix = (0.25*entr_in.qt_up + 0.75*entr_in.qt_env)
    # H_mix =  (0.25*entr_in.H_up  + 0.75*entr_in.H_env)
    qt_mix = (0.5*entr_in.qt_up + 0.5*entr_in.qt_env)
    H_mix =  (0.5*entr_in.H_up  + 0.5*entr_in.H_env)
    sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, qt_mix, H_mix)
    qv_ = (entr_in.qt_up+entr_in.qt_env)/2.0 - sa.ql
    rho_mix = rho_c(entr_in.p0, sa.T, qt_mix, qv_)
    b_mix = buoyancy_c(entr_in.rho0, rho_mix)-b_mean
    sorting_function = -(b_mix)/fmax(fabs(b_up-b_env),0.0000001)
    ret_b.b_mix = b_mix
    ret_b.sorting_function = sorting_function

    return ret_b
end

function buoyancy_sorting(entr_in::entr_struct)

    sqpi_inv = 1.0/sqrt(pi)
    sqrt2 = sqrt(2.0)
    sorting_function = 0.0
    inner_sorting_function = 0.0
    sa = eos_struct()

    abscissas, weights = np.polynomial.hermite.hermgauss(entr_in.quadrature_order)

    sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, entr_in.qt_env, entr_in.H_env)
    qv_ = entr_in.qt_env - sa.ql
    T_env = sa.T
    ql_env = sa.ql
    rho_env = rho_c(entr_in.p0, sa.T, entr_in.qt_env, qv_)
    b_env = buoyancy_c(entr_in.rho0, rho_env)

    sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, entr_in.qt_up, entr_in.H_up)
    qv_ = entr_in.qt_up - sa.ql
    T_up = sa.T
    ql_up = sa.ql
    rho_up = rho_c(entr_in.p0, sa.T, entr_in.qt_up, qv_)
    b_up = buoyancy_c(entr_in.rho0, rho_up)

    b_mean = entr_in.a_upd*b_up +  (1.0-entr_in.a_upd)*b_env

    if entr_in.env_QTvar != 0.0 && entr_in.env_Hvar != 0.0
        sd_q = sqrt(entr_in.env_QTvar)
        sd_h = sqrt(entr_in.env_Hvar)
        corr = fmax(fmin(entr_in.env_HQTcov/fmax(sd_h*sd_q, 1e-13),1.0),-1.0)

        # limit sd_q to prevent negative qt_hat
        sd_q_lim = (1e-10 - entr_in.qt_env)/(sqrt2 * abscissas[0])
        sd_q = fmin(sd_q, sd_q_lim)
        qt_var = sd_q * sd_q
        sigma_h_star = sqrt(fmax(1.0-corr*corr,0.0)) * sd_h

        @inbounds for m_q in xrange(entr_in.quadrature_order)
            qt_hat    = (entr_in.qt_env + sqrt2 * sd_q * abscissas[m_q] + entr_in.qt_up)/2.0
            mu_h_star = entr_in.H_env + sqrt2 * corr * sd_h * abscissas[m_q]
            inner_sorting_function = 0.0
            @inbounds for m_h in xrange(entr_in.quadrature_order)
                h_hat = (sqrt2 * sigma_h_star * abscissas[m_h] + mu_h_star + entr_in.H_up)/2.0
                # condensation - evaporation
                sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, qt_hat, h_hat)
                # calcualte buoyancy
                qv_ = qt_hat - sa.ql
                L_ = latent_heat(sa.T)
                dT = L_*((entr_in.ql_up+entr_in.ql_env)/2.0- sa.ql)/1004.0
                rho_mix = rho_c(entr_in.p0, sa.T, qt_hat, qv_)
                bmix = buoyancy_c(entr_in.rho0, rho_mix) - b_mean #- entr_in.dw2dz

                if bmix >0.0
                    inner_sorting_function  += weights[m_h] * sqpi_inv
                end
            end

            sorting_function  += inner_sorting_function * weights[m_q] * sqpi_inv
        end
    else
        h_hat = ( entr_in.H_env + entr_in.H_up)/2.0
        qt_hat = ( entr_in.qt_env + entr_in.qt_up)/2.0

        # condensation
        sa  = eos(t_to_thetali_c, eos_first_guess_thetal, entr_in.p0, qt_hat, h_hat)
        # calcualte buoyancy
        rho_mix = rho_c(entr_in.p0, sa.T, qt_hat, qt_hat - sa.ql)
        bmix = buoyancy_c(entr_in.rho0, rho_mix) - entr_in.b_mean
        if bmix  - entr_in.dw2dz >0.0
            sorting_function  = 1.0
        else
            sorting_function  = 0.0
        end
    end

    return sorting_function
end

function entr_detr_tke(entr_in::entr_struct)
    _ret = entr_struct()
    _ret.detr_sc = fabs(entr_in.b_upd)/ fmax(entr_in.w_upd * entr_in.w_upd, 1e-3)
    _ret.entr_sc = sqrt(entr_in.tke) / fmax(entr_in.w_upd, 0.01) / fmax(sqrt(entr_in.a_upd), 0.001) / 50000.0
    return  _ret
end

function entr_detr_b_w2(entr_in::entr_struct)
    _ret = entr_struct()
    # in cloud portion from Soares 2004
    if entr_in.z >= entr_in.zi
        _ret.detr_sc= 4.0e-3 + 0.12 *fabs(fmin(entr_in.b_upd,0.0)) / fmax(entr_in.w_upd * entr_in.w_upd, 1e-2)
    else
        _ret.detr_sc = 0.0
    end

    _ret.entr_sc = 0.12 * fmax(entr_in.b_upd,0.0) / fmax(entr_in.w_upd * entr_in.w_upd, 1e-2)

    return  _ret
end

function entr_detr_suselj(entr_in::entr_struct)
    _ret = entr_struct()
    entr_dry = 2.5e-3

    l0 = (entr_in.zbl - entr_in.zi)/10.0
    if entr_in.z >= entr_in.zi
        _ret.detr_sc= 4.0e-3 +  0.12* fabs(fmin(entr_in.b_upd,0.0)) / fmax(entr_in.w_upd * entr_in.w_upd, 1e-2)
        _ret.entr_sc = 0.002 # 0.1 / entr_in.dz * entr_in.poisson

    else
        _ret.detr_sc = 0.0
        _ret.entr_sc = 0.0 #entr_dry # Very low entrainment rate needed for Dycoms to work
    end

    return  _ret
end


function entr_detr_none(entr_in::entr_struct)
    _ret = entr_struct()
    _ret.entr_sc = 0.0
    _ret.detr_sc = 0.0

    return  _ret
end

function pressure_tan18_buoy(press_in::pressure_in_struct)
    _ret = pressure_buoy_struct()

    _ret.b_coeff = press_in.bcoeff_tan18
    _ret.nh_pressure_b = -1.0 * press_in.rho0_kfull * press_in.a_kfull * press_in.b_kfull * _ret.b_coeff

    return _ret
end

function pressure_tan18_drag(press_in::pressure_drag_struct)
    _ret = pressure_drag_struct()

    _ret.nh_pressure_adv = 0.0
    _ret.nh_pressure_drag = -1.0 * press_in.rho0_kfull * sqrt(press_in.a_kfull)* sqrt(press_in.a_kfull) * (1.0/press_in.rd
                          * (press_in.w_kfull - press_in.w_kenv)*fabs(press_in.w_kfull - press_in.w_kenv))
    return _ret
end

function pressure_normalmode_buoy(press_in::pressure_in_struct)
    _ret = pressure_buoy_struct()
    _ret.b_coeff = press_in.alpha1 / ( 1+press_in.alpha2*press_in.asp_ratio^2 )
    _ret.nh_pressure_b = -1.0 * press_in.rho0_kfull * press_in.a_kfull * press_in.b_kfull * _ret.b_coeff

    return _ret
end

function pressure_normalmode_drag(press_in::pressure_in_struct)
    _ret = pressure_drag_struct()

    _ret.nh_pressure_adv = press_in.rho0_kfull * press_in.a_kfull * press_in.beta1*press_in.w_kfull*(press_in.w_kfull
                          -press_in.w_kmfull)*press_in.dzi

    # drag as w_dif and account for downdrafts
    _ret.nh_pressure_drag = -1.0 * press_in.rho0_kfull * press_in.a_kfull * press_in.beta2 * (press_in.w_kfull -
                            press_in.w_kenv)*fabs(press_in.w_kfull - press_in.w_kenv)/fmax(press_in.updraft_top, 500.0)

    return _ret
end

# convective velocity scale
function get_wstar(bflux, zi)
    return cbrt(fmax(bflux * zi, 0.0))
end

# BL height
function get_inversion(theta_rho, u, v, z_half,
                          kmin, kmax, Ri_bulk_crit)
    theta_rho_b = theta_rho[kmin]
    k = kmin
    Ri_bulk=0.0
    Ri_bulk_low = 0.0


    # test if we need to look at the free convective limit
    if (u[kmin] * u[kmin] + v[kmin] * v[kmin]) <= 0.01
        @inbounds for k in xrange(kmin,kmax)
            if theta_rho[k] > theta_rho_b
                break
            end
        end
        h = (z_half[k] - z_half[k-1])/(theta_rho[k] - theta_rho[k-1]) * (theta_rho_b - theta_rho[k-1]) + z_half[k-1]
    else
        @inbounds for k in xrange(kmin,kmax)
            Ri_bulk_low = Ri_bulk
            Ri_bulk = g * (theta_rho[k] - theta_rho_b) * z_half[k]/theta_rho_b / (u[k] * u[k] + v[k] * v[k])
            if Ri_bulk > Ri_bulk_crit
                break
            end
        end
        h = (z_half[k] - z_half[k-1])/(Ri_bulk - Ri_bulk_low) * (Ri_bulk_crit - Ri_bulk_low) + z_half[k-1]
    end

    return h
end
# Teixiera convective tau
function get_mixing_tau(zi, wstar)
    # return 0.5 * zi / wstar
    #return zi / (fmax(wstar, 1e-5))
    return zi / (wstar + 0.001)
end




# MO scaling of near surface tke and scalar variance

function get_surface_tke(ustar, wstar, zLL, oblength)
    if oblength < 0.0
        return ((3.75 + cbrt(zLL/oblength * zLL/oblength)) * ustar * ustar)
    else
        return (3.75 * ustar * ustar)
    end
end
function get_surface_variance(flux1, flux2, ustar, zLL, oblength)
    c_star1 = -flux1/ustar
    c_star2 = -flux2/ustar
    if oblength < 0.0
        return 4.0 * c_star1 * c_star2 * pow(1.0 - 8.3 * zLL/oblength, -2.0/3.0)
    else
        return 4.0 * c_star1 * c_star2
    end
end


# Math-y stuff
function construct_tridiag_diffusion(
    nzg,
    gw,
    dzi,
    dt,
    rho_ae_K_m,
    rho,
    ae,
    a,
    b,
    c)
    nz = nzg - 2* gw
    @inbounds for k in xrange(gw,nzg-gw)
        X = rho[k] * ae[k]/dt
        Y = rho_ae_K_m[k] * dzi * dzi
        Z = rho_ae_K_m[k-1] * dzi * dzi
        if k == gw
            Z = 0.0
        elseif k == nzg-gw-1
            Y = 0.0
        end
        a[k-gw] = - Z/X
        b[k-gw] = 1.0 + Y/X + Z/X
        c[k-gw] = -Y/X
    end
    return
end

function construct_tridiag_diffusion_implicitMF(
    nzg,
    gw,
    dzi,
    dt,
    rho_ae_K_m,
    massflux,
    rho,
    alpha,
    ae,
    a,
    b,
    c)
    nz = nzg - 2* gw
    @inbounds for k in xrange(gw,nzg-gw)
        X = rho[k] * ae[k]/dt
        Y = rho_ae_K_m[k] * dzi * dzi
        Z = rho_ae_K_m[k-1] * dzi * dzi
        if k == gw
            Z = 0.0
        elseif k == nzg-gw-1
            Y = 0.0
        end
        a[k-gw] = - Z/X + 0.5 * massflux[k-1] * dt * dzi/rho[k]
        b[k-gw] = 1.0 + Y/X + Z/X + 0.5 * dt * dzi * (massflux[k-1]-massflux[k])/rho[k]
        c[k-gw] = -Y/X - 0.5 * dt * dzi * massflux[k]/rho[k]
    end
    return
end

function construct_tridiag_diffusion_dirichlet(
    nzg,
    gw,
    dzi,
    dt,
    rho_ae_K_m,
    rho,
    ae,
    a,
    b,
    c)
    nz = nzg - 2* gw
    @inbounds for k in xrange(gw,nzg-gw)
        X = rho[k] * ae[k]/dt
        Y = rho_ae_K_m[k] * dzi * dzi
        Z = rho_ae_K_m[k-1] * dzi * dzi
        if k == gw
            Z = 0.0
            Y = 0.0
        elseif k == nzg-gw-1
            Y = 0.0
        end
        a[k-gw] = - Z/X
        b[k-gw] = 1.0 + Y/X + Z/X
        c[k-gw] = -Y/X
    end
    return
end



function tridiag_solve(nz, x, a, b, c)
    scratch = pyzeros(length(x))
    scratch[0] = c[0]/b[0]
    x[0] = x[0]/b[0]
    @inbounds for i in xrange(1,nz)
        m = 1.0/(b[i] - a[i] * scratch[i-1])
        scratch[i] = c[i] * m
        x[i] = (x[i] - a[i] * x[i-1])*m
    end
    # TODO: verify translation
    @inbounds for i in revxrange(nz-2,0,-1)
        x[i] = x[i] - scratch[i] * x[i+1]
    end
    return
end

# Dustbin

function set_cloudbase_flag(ql, current_flag)
    if ql > 1.0e-8
        new_flag = true
    else
        new_flag = current_flag
    end
    return  new_flag
end
