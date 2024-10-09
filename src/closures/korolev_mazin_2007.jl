"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function korolev_mazin_2007(param_set::APS, area::FT, ρ::FT, Δt::Real, ts::TD.ThermodynamicState, w::FT) where {FT}
    """
    Need to add w in EDMF_Env and EDMF_Up
    # still no partitioning like tan et al in space
    # no immersion frz etc
    # no activation etc...
    """
    thermo_params = TCP.thermodynamics_params(param_set) # currently repeated in both places, pare down later
    microphys_params = TCP.microphysics_params(param_set)
    S_ql = 0
    S_qi = 0
    if area > 0

        q = TD.PhasePartition(thermo_params, ts)
        T = TD.air_temperature(thermo_params, ts)
        q_vap = TD.vapor_specific_humidity(thermo_params, ts)
        p = TD.air_pressure(thermo_params, T, ρ, q)

        q_v = q_vap
        qsvap_l = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Liquid())
        qsvap_i = TD.q_vap_saturation_generic(thermo_params, T, ρ, TD.Ice())

        c = 1 # ice particle shape factor characterizing ‘capacitance’ 0 < c ≤ 1 (c = 1 for spheres)
        D = FT(0.0000226) # TCP.CM.Parameters.D_vapor(param_set) # CP.diffusivity_of_water_vapor(param_set) # coefficient of water vapour diffusion in the air [0.0000226] # https://github.com/CliMA/CLIMAParameters.jl/blob/2f298661d527c1b9a11188ee2abc92ba4ba0c5ec/src/parameters.toml#L137
        e = TD.partial_pressure_vapor(thermo_params, p, q) # water vapour pressure
        c_p = TD.cp_m(thermo_params, ts)# specific heat capacity of moist air at constant pressure
        k = FT(0.024) # coefficient of air heat conductivity # https://github.com/CliMA/CLIMAParameters.jl/blob/2f298661d527c1b9a11188ee2abc92ba4ba0c5ec/src/parameters.toml#L132


        E_i = qsvap_l # saturation vapour pressure above flat surface of ice
        E_l = qsvap_i # saturation vapour pressure above flat surface of water
        ξ = E_l / E_i

        g = TCP.grav(param_set) # acceleration of gravity

        L_i = TD.latent_heat_sublim(thermo_params, ts) # Latent heat for ice
        L_l = TD.latent_heat_vapor(thermo_params, ts)  # Latent heat for water
        R_a = TD.gas_constant_air(thermo_params, q) # specific gas constant of moist air
        R_v = TCP.R_v(param_set)# specific gas constant of water vapour
        ρ_l = TCP.ρ_cloud_liq(param_set) # density of liquid water
        ρ_i = FT(916.7) # TCP.ρ_cloud_ice(param_Set) # density of an ice particle # https://github.com/CliMA/CLIMAParameters.jl/blob/2f298661d527c1b9a11188ee2abc92ba4ba0c5ec/src/parameters.toml#L513


        a_0 = (g / (R_a * T)) * ((L_l * R_a) / (c_p * R_v * T) - 1)
        a_1 = (1 / q_v) + (L_l * L_l) / (c_p * R_v * T^2)
        a_2 = (1 / q_v) + (L_l * L_i) / (c_p * R_v * T^2)

        A_i = ((ρ_i * L_i^2) / (k * R_v * T^2) + (ρ_i * R_v * T) / (E_i * D))^-1
        A_l = ((ρ_l * L_l^2) / (k * R_v * T^2) + (ρ_l * R_v * T) / (E_l * D))^-1
        B_i0 = (4 * π * ρ_i * A_i) / ρ
        B_i = ξ * c * B_i0
        B_is = (ξ - 1) * c * B_i0
        B_l = (4 * π * ρ_l * A_l) / ρ


        b_l = a_1 * B_l
        b_i = a_2 * B_i
        b_i0 = a_2 * B_i0
        b_is = a_2 * B_is

        # currently these effectively set the τ 

        # N_l, r_l = FT(10^6), FT(10^-5) # τ ~ 1/(Nr) = 1/(10) # unstable (maybe more stable than just setting tau)
        # N_i, r_i = FT(10^3), FT(10^-2) # τ ~ 1/(Nr) = 1/(10)

        # N_l, r_l = FT(10^4), FT(10^-5) # τ ~ 1/(Nr) = 1/(1/10) = 10 # plot looked same...
        # N_i, r_i = FT(10^1), FT(10^-2) # τ ~ 1/(Nr) = 1/(1/10) = 10 

        N_l, r_l = FT(10^2), FT(10^-5) # τ ~ 1/(Nr) = 1/(1/1000) = 1000 # plot looked same...
        N_i, r_i = FT(10^-1), FT(10^-2) # τ ~ 1/(Nr) = 1/(1/1000) = 1000

        denom = b_l * N_l * r_l + b_i * N_i * r_i
        dql_dt = (a_0 * w - b_is * N_i * r_i) * (B_l * N_l * r_l) / denom
        dqi_dt = ((a_0 * w - (1 - ξ) / ξ * (b_l * N_l * r_l)) * (B_i * N_i * r_i)) / denom

        S_ql += dql_dt
        S_qi += dqi_dt

        ## === LIMITERS === ## (left most external to this)

        # Hopefully everything converges and no temp limiters at T > 0

        # korolev bounds might be useful for limiters also...
        w_0 = FT(NaN)
        w_p = FT(NaN)
        w_s = FT(NaN)
        if w > w_s  # liq -->  vap <-- ice
        elseif w > w_p  # liq ==>  vap --> ice
        elseif w > w_0  # liq -->  vap ==? ice
        else            # liq <--  vap --> ice
        end

        ## === LIMITERS === ##

    end
    return S_ql, S_qi
end
