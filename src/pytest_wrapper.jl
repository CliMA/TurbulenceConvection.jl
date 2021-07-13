# wrapper from cdef (not accessible from Python) to def (accessible from Python)
# written to be able to use the thermodynamic functions from Pytest

cimport thermodynamic_functions as fun
include "parameters.pxi"

cdef class TurbulenceConvection_constants
    def __init__(self)
        self.g  = g
        self.Rd = Rd

def pd_c(p0, qt, qv)
    return fun.pd_c(p0, qt, qv)

def pv_c(p0, qt, qv)
    return fun.pv_c(p0, qt, qv)

def pv_star(T)
    return fun.pv_star(T)

# assuming qt = qv + ql + qi
def qv_star_c(p0, qt, pv)
    return fun.qv_star_c(p0, qt, pv)

# assuming qt = qv
def qv_star_t(p0, T)
    return fun.qv_star_t(p0, T)

def latent_heat(T)
    return fun.latent_heat(T)

# t_to_prog
def t_to_thetali_c(p0, T, qt, ql, qi)
    return fun.t_to_thetali_c(p0, T, qt, ql, qi)

# prog_to_t (first guess assuming no ql)
def eos_first_guess_thetal(H, pd, pv, qt)
    return fun.eos_first_guess_thetal(H, pd, pv, qt)

# prog_to_t (saturation adjustment)
def eos(p0, qt, prog) 
    return fun.eos(fun.t_to_thetali_c, fun.eos_first_guess_thetal, p0, qt, prog)

def theta_c(p0, T)
    return fun.theta_c(p0, T)

