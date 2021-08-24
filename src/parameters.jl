#Adapated from PyCLES: https://github.com/pressel/pycles
const Rd = 287.1
const Rv = 461.5
const eps_v = 0.62210184182   # Rd / Rv
const eps_vi = 1.60745384883  # Rv / Rd
const cpd = 1004.0
const cpv = 1859.0
const kappa = 0.285956175299
const T_tilde = 298.15
const p_tilde = 100000.0
const sd_tilde = 6864.8
const sv_tilde = 10513.6
const omega = 7.29211514671e-05
const vkb = 0.4
const Pr0 = 0.74
const beta_m = 4.8
const beta_h = 6.5 # beta_h = beta_m / Pr0
const gamma_m = 15.0
const gamma_h = 9.0
# CLIMA microphysics parameters
const rho_cloud_liq = 1e3
const nu_air = 1.6e-5
const K_therm = 2.4e-2
const D_vapor = 2.26e-5
