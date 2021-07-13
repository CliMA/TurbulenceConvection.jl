import sys
sys.path.insert(0, "./../")

import os
import subprocess
import json
import warnings

from netCDF4 import Dataset

import pprint as pp

import math as mt

import pytest
import pytest_wrapper as wrp

# https://hypothesis.readthedocs.io/en/latest/
from hypothesis import given, strategies as st


@given(p0   = st.floats(min_value = 12000, max_value = 101300), # pressure between 1013hPa - 300hPa
       qt   = st.floats(min_value = 0,     max_value = 0.040),  # total water specific humidity between 0 - 40 g/kg
       frac = st.floats(min_value = 0,     max_value = 1))      # fraction of qt that is water vapor specific humidity
def test_pressure(p0, qt, frac):
    """
    Tests functions pd_c and pd_v from thermodynamic_functions.pyx
    by checking if the total pressure = dry air pressure + water vapor pressure.
    """
    qv = frac * qt

    pd = wrp.pd_c(p0, qt, qv)
    pv = wrp.pv_c(p0, qt, qv)

    assert(np.isclose(p0, pd + pv, rtol = 1e-12))


@given(p0 = st.floats(min_value = 12000,  max_value = 101300),      # pressure between 1013hPa - 300hPa
       qt = st.floats(min_value = 0,      max_value = 0.040),       # total water specific humidity between 0 - 40 g/kg
       T  = st.floats(min_value = 273.15, max_value = 273.15 + 40)) # temperature between 0 - 40 C
def test_eos_first_guess(p0, qt, T):
    """
    Tests functions t_to_thetali and eos_first_guess from thermodynamic_functions.pyx
    by checking if calculating thetali using T ad thet calculating T using thetali
    results in the same temperature

    checks subsaturated cases
    """
    qv = qt
    ql = 0.
    qi = 0.

    pd = wrp.pd_c(p0, qt, qv)
    pv = wrp.pv_c(p0, qt, qv)

    # calculate thetali(T) assuming no liquid water or ice are present
    thetali = wrp.t_to_thetali_c(p0, T, qt, ql, qi)
    # calculate T(thetali) assuming no liquid water or ice are present
    T_new = wrp.eos_first_guess_thetal(thetali, pd, pv, qt)

    assert(np.isclose(T, T_new, rtol = 1e-12))


@given(z  = st.floats(min_value = 0,     max_value = 18000), # height in the atmosphere around which we will test
       rnd= st.floats(min_value = -0.05, max_value = 0.05),  # random factor
       ql = st.floats(min_value = 0,     max_value = 0.010)) # liquid water specific humidity between 0 - 10 g/kg
def test_eos_saturated(z, rnd, ql):
    """
    Check if the saturation adjustment scheme finds the correct T and ql
    """
    # constructing initial condition such that
    # a) it"s saturated
    # b) there is some ql
    # c) p and T are close to conditions possible in the atmosphere

    # TurbulenceConvection constants
    const = wrp.TurbulenceConvection_constants()

    # "standard atmosphere" profile around which we will test
    # https://en.wikipedia.org/wiki/International_Standard_Atmosphere
    dTdz = -0.0065
    T0 = 273.15 + 15
    T1 = 273.15 - 56.5
    p0 = 101300.
    p1 = 22632
    z1 = 11000
    def T_standard(z):
        T = T1
        if z < z1:
            T = T0 + dTdz * z
        return T
    def p_standard(z):
        p = p0
        if z < z1:
            p = p0 * mt.exp(-const.g / const.Rd / dTdz * mt.log(1. + dTdz / T0 * z))
        else:
            p = p1 * mt.exp(-const.g / const.Rd / T1 * (z - z1))
        return p

    T    = T_standard(z) / (1. - rnd)
    p    = p_standard(T) / (1. - rnd)
    pv_s = wrp.pv_star(T)
    # qv_star assuming that qt = qv + ql
    qv_s = wrp.qv_star_t(p, T) * (1 - ql)
    qt   = qv_s + ql
    # initial prognostic variable thetali
    thetali = wrp.t_to_thetali_c(p, T, qt, ql, 0)

    # saturation adjustment
    res = wrp.eos(p, qt, thetali)

    assert(np.isclose(T,  res["T"],  rtol = 1e-3))
    if (ql > 1e-7):
        assert(np.isclose(ql, res["ql"], rtol = 1e-2))
