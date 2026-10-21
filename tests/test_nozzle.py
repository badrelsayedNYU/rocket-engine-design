"""Independent isentropic checks and operating limits; no CEA dependency."""
from math import sqrt
import warnings

import numpy as np
import pytest

from engine_design.atmosphere import ambient_pressure
from engine_design.nozzle import (NozzleWarning, area_ratio, mach_from_area_ratio,
                                  pressure_ratio, thrust_coefficient, isp,
                                  optimal_eps, separation_check)


def test_textbook_mach_two():
    # gamma=1.4, M=2: area=(1/2)*(1.5)^3=1.6875;
    # p/p0=(1.8)^(-3.5)=0.1278045255. Spec Stage 2 / NASA isentropic relations.
    assert area_ratio(2, 1.4) == pytest.approx(1.6875, rel=1e-12)
    assert mach_from_area_ratio(1.6875, 1.4) == pytest.approx(2, rel=1e-12)
    assert pressure_ratio(2, 1.4) == pytest.approx(0.12780452546295096)


@pytest.mark.parametrize("gamma", [1.01, 1.2, 1.4, 5/3])
@pytest.mark.parametrize("eps", [1, 1 + 1e-10, 1.01, 10, 1000])
def test_inverse_and_supersonic_branch(gamma, eps):
    mach = mach_from_area_ratio(eps, gamma)
    assert mach >= 1
    assert area_ratio(mach, gamma) == pytest.approx(eps, rel=1e-10)


def test_sonic_and_stagnation_limits():
    assert mach_from_area_ratio(1, 1.4) == 1
    assert pressure_ratio(0, 1.4) == 1
    assert pressure_ratio(1, 1.4) == pytest.approx(0.5282817877171742)
    pc = 1e7
    assert optimal_eps(1.4, pc, pc * pressure_ratio(1, 1.4)) == pytest.approx(1)


def test_sutton_chart_vacuum_asymptotes():
    # MIT OCW 16.50 Lecture 7, PDF p.3, plots after Sutton & Biblarz 2001.
    # Explicit chart labels: max CF=2.246 (gamma=1.2), 1.964 (gamma=1.3).
    # https://ocw.mit.edu/courses/16-50-introduction-to-propulsion-systems-spring-2012/e2035e49d0312366d0772766f594b671_MIT16_50S12_lec7.pdf
    # Very large finite eps approximates the chart's infinite-expansion limit.
    assert thrust_coefficient(1.2, 1e20, 1e7, 0) == pytest.approx(2.246, abs=0.001)
    assert thrust_coefficient(1.3, 1e20, 1e7, 0) == pytest.approx(1.964, abs=0.001)


def test_thrust_from_mass_flow_and_exit_velocity():
    # Independently reconstruct F=mdot*ve+(pe-pa)*Ae for an analyst-defined
    # perfect gas. Throat density and sonic speed determine mass flux.
    gamma, r, tc, pc, at, mach = 1.4, 287.0, 3000.0, 1e7, 0.01, 2.0
    eps, pa = 1.6875, 101325.0
    tt = tc / 1.2
    pt = pc * 0.5282817877171742
    mdot = pt / (r * tt) * sqrt(gamma * r * tt) * at
    te = tc / 1.8
    pe = pc * 0.12780452546295096
    ve = mach * sqrt(gamma * r * te)
    thrust = mdot * ve + (pe - pa) * eps * at
    cf = thrust_coefficient(gamma, eps, pc, pa)
    assert cf == pytest.approx(thrust / (pc * at), rel=1e-12)
    assert isp(pc * at / mdot, cf) == pytest.approx(thrust / (mdot * 9.80665))


def test_vacuum_pressure_penalty_and_optimum():
    pc, pa, gamma = 1e7, 101325.0, 1.2
    eps = optimal_eps(gamma, pc, pa)
    assert pressure_ratio(mach_from_area_ratio(eps, gamma), gamma) * pc == pytest.approx(pa)
    best = thrust_coefficient(gamma, eps, pc, pa)
    assert best > thrust_coefficient(gamma, eps * 0.99, pc, pa)
    assert best > thrust_coefficient(gamma, eps * 1.01, pc, pa)
    vacuum = thrust_coefficient(gamma, eps, pc, 0)
    assert vacuum > best
    assert vacuum - best == pytest.approx(eps * pa / pc)


def test_separation_threshold_and_extrapolation():
    with pytest.warns(NozzleWarning, match="separation"):
        assert separation_check(39999, 100000)
    assert not separation_check(40000, 100000)
    assert not separation_check(0, 0)
    assert not separation_check(100000, 100000)
    with pytest.warns(NozzleWarning):
        cf = thrust_coefficient(1.2, 100, 1e6, 1e5)
    assert cf < 0  # Never clip a failed attached-flow extrapolation.
    assert isp(1800, cf) < 0


@pytest.mark.parametrize("call,args", [
    (mach_from_area_ratio, (0.9, 1.2)), (mach_from_area_ratio, (float('inf'), 1.2)),
    (area_ratio, (0, 1.2)), (pressure_ratio, (-1, 1.2)),
    (pressure_ratio, (2, 1)), (pressure_ratio, (2, float('nan'))),
    (thrust_coefficient, (1.2, 10, 0, 0)), (thrust_coefficient, (1.2, 10, 1e7, -1)),
    (thrust_coefficient, (1.2, 10, 1e7, 1e7)), (isp, (0, 1.5)),
    (isp, (1800, float('nan'))), (optimal_eps, (1.2, 1e7, 0)),
    (optimal_eps, (1.4, 1e7, 6e6)), (separation_check, (1, 2, 1.1)),
    (separation_check, (-1, 2)),
])
def test_invalid_inputs(call, args):
    with pytest.raises(ValueError):
        call(*args)


@pytest.mark.parametrize("height,expected", [(0, 101325), (11000, 22632.06),
    (20000, 5474.89), (32000, 868.02), (47000, 110.91), (51000, 66.94)])
def test_standard_atmosphere_reference(height, expected):
    # NASA/TM-2005-213659 Table A1; tolerance covers tabulated rounding.
    assert ambient_pressure(height) == pytest.approx(expected, abs=0.01)


def test_atmosphere_continuity_and_altitude_isp():
    for height in [11000, 20000, 32000, 47000]:
        assert ambient_pressure(height - 1e-5) == pytest.approx(ambient_pressure(height + 1e-5), rel=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', NozzleWarning)
        performance = [isp(1800, thrust_coefficient(1.2, 40, 1e7, ambient_pressure(h)))
                       for h in np.linspace(0, 50000, 30)]
    assert np.all(np.diff(performance) > 0)
    assert performance[-1] < isp(1800, thrust_coefficient(1.2, 40, 1e7, 0))


@pytest.mark.parametrize("height", [-1, 51001, float('nan')])
def test_atmosphere_range(height):
    with pytest.raises(ValueError):
        ambient_pressure(height)
