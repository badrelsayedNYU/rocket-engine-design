"""Stage 4: divergence factors, efficiencies, frozen vs shifting."""
from math import atan, degrees, radians, tan
import warnings

import pytest

from engine_design.constants import (BELL_DIVERGENCE_AREA_RATIOS, BELL_DIVERGENCE_FRACTIONS,
                                     BELL_DIVERGENCE_PERCENT)
from engine_design.losses import (Performance, apply_losses, bell_divergence_factor,
                                  conical_divergence_factor, ideal_performance)


def test_conical_fifteen_degrees():
    # Theory section 4 / NASA SP-125 eq. 4-8: (1+cos 15°)/2 = 0.98296.
    assert conical_divergence_factor(15) == pytest.approx(0.983, abs=5e-4)
    assert conical_divergence_factor(0) == 1


@pytest.mark.parametrize("fraction,reading", [(0.6, 95.70), (0.7, 96.67), (0.8, 97.43),
                                               (0.9, 97.92), (1.0, 98.28)])
def test_digitization_calibrated_by_cone_curve(fraction, reading):
    # The same pixel extraction that produced the bell table read SP-125
    # Fig.4-12's conical curve at these points. A cone of equal length has
    # tan(alpha) = tan(15°)/Lf; eq. 4-8 then predicts lambda independently.
    alpha = degrees(atan(tan(radians(15)) / fraction))
    assert 100 * conical_divergence_factor(alpha) == pytest.approx(reading, abs=0.1)


def test_bell_table_nodes_and_trends():
    for i, eps in enumerate(BELL_DIVERGENCE_AREA_RATIOS):
        for j, fraction in enumerate(BELL_DIVERGENCE_FRACTIONS):
            assert bell_divergence_factor(eps, fraction) == pytest.approx(BELL_DIVERGENCE_PERCENT[i][j] / 100)
    # A bell beats a cone of the same length; longer bells lose less.
    assert bell_divergence_factor(10, 0.8) > conical_divergence_factor(15)
    assert bell_divergence_factor(20, 0.6) < bell_divergence_factor(20, 0.8) < bell_divergence_factor(20, 1.0)
    # Bilinear midpoint of eps 10–20 and Lf 0.7–0.8 is the mean of four nodes.
    assert bell_divergence_factor(15, 0.75) == pytest.approx(
        (97.58 + 98.55 + 98.00 + 98.71) / 4 / 100, rel=1e-12)


@pytest.mark.parametrize("eps,fraction", [(9.9, 0.8), (40.1, 0.8), (20, 0.59), (20, 1.01)])
def test_bell_no_extrapolation(eps, fraction):
    with pytest.raises(ValueError):
        bell_divergence_factor(eps, fraction)


def test_apply_losses_hand_calculation():
    # Analyst-defined numbers. CF_vac = 1.9, eps*pe/pc = 40*20e3/10e6 = 0.08,
    # momentum = 1.82. CF_del = 0.98*(0.99*1.82 + 0.08) = 0.98*1.8818 = 1.844164.
    ideal = ideal_performance(1.9, 2300.0, 20e3, 40, 10e6, 101325)
    assert ideal.CF_amb == pytest.approx(1.9 - 40 * 101325 / 10e6)
    d = apply_losses(ideal, 40, 10e6, 101325, 0.99, 0.97, 0.98)
    assert d.CF_vac == pytest.approx(1.844164, rel=1e-12)
    assert d.c_star_m_s == pytest.approx(0.97 * 2300)
    # Exact ambient pressure force is not scaled by efficiencies.
    assert d.CF_vac - d.CF_amb == pytest.approx(ideal.CF_vac - ideal.CF_amb)
    assert d.p_e_Pa == ideal.p_e_Pa
    # Vacuum: Isp_del = eta_c* * eta_CF * Isp(with lambda), theory section 4.
    with_lambda = 2300 * (0.99 * 1.82 + 0.08) / 9.80665
    assert d.Isp_vac_s == pytest.approx(0.97 * 0.98 * with_lambda)
    assert d.Isp_vac_s < ideal.Isp_vac_s and d.Isp_amb_s < ideal.Isp_amb_s


def test_unit_efficiencies_are_identity():
    ideal = ideal_performance(1.85, 1800.0, 50e3, 20, 7e6, 101325)
    assert apply_losses(ideal, 20, 7e6, 101325, 1, 1, 1) == ideal


@pytest.mark.parametrize("args", [(0, 0.98, 0.98), (0.99, 1.01, 0.98),
                                  (0.99, 0.98, float("nan")), (1.2, 0.98, 0.98)])
def test_invalid_efficiencies(args):
    ideal = Performance(1800.0, 1.85, 1.7, 50e3)
    with pytest.raises(ValueError):
        apply_losses(ideal, 20, 7e6, 101325, *args)


@pytest.mark.cea
@pytest.mark.parametrize("pair,of,eps", [("LOX/LH2", 6, 40), ("LOX/RP-1", 2.5, 16),
                                         ("LOX/CH4", 3.5, 30)])
def test_frozen_below_shifting(pair, of, eps):
    # Shifting equilibrium recovers recombination energy during expansion.
    pytest.importorskip("rocketcea")
    from engine_design.thermo import expansion_performance
    shifting = expansion_performance(pair, 10e6, of, eps)
    frozen = expansion_performance(pair, 10e6, of, eps, frozen=True)
    assert frozen.Isp_vac_s < shifting.Isp_vac_s
    assert frozen.p_e_Pa < shifting.p_e_Pa  # less energy left in the gas


@pytest.mark.cea
def test_cea_momentum_split_independent(tmp_path):
    # CF_vac = v_e/c* + eps*p_e/p_c. Rebuild v_e from CEA exit Mach and
    # sonic velocity, independently of the Isp used by expansion_performance.
    pytest.importorskip("rocketcea")
    from engine_design.thermo import _cea, expansion_performance
    perf = expansion_performance("LOX/CH4", 10e6, 3.4, 20)
    with _cea("LOX/CH4") as cea:
        mach = cea.get_MachNumber(Pc=10e6, MR=3.4, eps=20)
        sonic = cea.get_SonicVelocities(Pc=10e6, MR=3.4, eps=20)[2]
    momentum = perf.CF_vac - 20 * perf.p_e_Pa / 10e6
    assert momentum == pytest.approx(mach * sonic / perf.c_star_m_s, rel=1e-3)


@pytest.fixture(scope="module")
def designs():
    pytest.importorskip("rocketcea")
    from engine_design import EngineRequirements, design_engine
    base = dict(propellants="LOX/CH4", thrust_N=100e3, p_c_Pa=10e6, mixture_ratio=3.4)
    return (design_engine(EngineRequirements(**base)),
            design_engine(EngineRequirements(**base, divergence_factor=1, eta_c_star=1, eta_CF=1)))


@pytest.mark.cea
def test_design_sized_on_delivered_values(designs):
    real, ideal = designs
    assert real.thrust_vac_N == pytest.approx(100e3, rel=1e-12)
    assert ideal.delivered == ideal.ideal
    assert real.delivered.Isp_vac_s < real.ideal.Isp_vac_s
    assert real.ideal_frozen.Isp_vac_s < real.ideal.Isp_vac_s
    # Same F and p_c: At scales with 1/CF_del; mdot = pc*At/(eta_c* c*).
    assert real.geometry.A_t_m2 / ideal.geometry.A_t_m2 == pytest.approx(
        ideal.delivered.CF_vac / real.delivered.CF_vac, rel=1e-12)
    assert real.geometry.mdot_kg_s / ideal.geometry.mdot_kg_s == pytest.approx(
        real.geometry.A_t_m2 / ideal.geometry.A_t_m2 / 0.98, rel=1e-12)
    assert real.divergence_factor == pytest.approx(bell_divergence_factor(real.eps, 0.8))


@pytest.mark.cea
def test_user_angles_require_divergence_factor():
    pytest.importorskip("rocketcea")
    from engine_design import EngineRequirements, design_engine
    req = EngineRequirements(propellants="LOX/CH4", thrust_N=100e3, p_c_Pa=10e6,
                             mixture_ratio=3.4, nozzle=20.0, theta_n_deg=30, theta_e_deg=9)
    with pytest.raises(ValueError, match="divergence_factor"):
        design_engine(req)


@pytest.mark.parametrize("kwargs", [dict(eta_c_star=0), dict(eta_CF=1.1),
                                    dict(divergence_factor=float("nan"))])
def test_invalid_design_efficiencies(kwargs):
    from engine_design import EngineRequirements, design_engine
    with pytest.raises(ValueError):
        design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=1e5, **kwargs))
