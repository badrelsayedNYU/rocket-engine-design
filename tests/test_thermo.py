"""Analytical checks, CEA integration, SI conversions, and optimization."""
import warnings

import numpy as np
import pytest

from engine_design.constants import C_STAR_SANITY_M_S, DEFAULT_OF_BOUNDS
from engine_design.thermo import (
    ThermoWarning, chamber_props, ideal_c_star, optimal_OF, sweep_OF, _cea, _isp,
)


@pytest.fixture
def cea_backend():
    pytest.importorskip("rocketcea")


def test_perfect_gas_hand_calculation():
    # For gamma=1.4: ((gamma+1)/2)^((gamma+1)/(2(gamma-1))) = 1.2^3.
    # Analyst-defined gas: R=287 J/(kg K), T=300 K -> c*=428.5300642895 m/s.
    assert ideal_c_star(300, 8314.462618 / 287, 1.4) == pytest.approx(428.5300642895)


@pytest.mark.parametrize("args", [(0, 20, 1.2), (300, 0, 1.2),
                                  (300, 20, 1), (float('nan'), 20, 1.2)])
def test_invalid_ideal_inputs(args):
    with pytest.raises(ValueError):
        ideal_c_star(*args)


@pytest.mark.cea
@pytest.mark.parametrize("pair,of", [("LOX/LH2", 6), ("LOX/RP-1", 2.5), ("LOX/CH4", 3.5)])
def test_nominal_properties_and_crosscheck(cea_backend, pair, of):
    p = chamber_props(pair, 10e6, of)
    low, high = C_STAR_SANITY_M_S[pair]
    assert low <= p.c_star_m_s <= high
    assert 1.13 <= p.gamma <= 1.25
    assert 0 < p.T_t_K < p.T_c_K
    assert ideal_c_star(p.T_c_K, p.M_kg_kmol, p.gamma) == pytest.approx(p.c_star_m_s, rel=0.02)
    # Pr = mu*cp/k independently checks all three transport conversions.
    for tr in (p.chamber_transport, p.throat_transport):
        assert tr.mu_Pa_s * tr.cp_J_kg_K / tr.conductivity_W_m_K == pytest.approx(tr.Pr, rel=1e-8)
        assert 1e-5 < tr.mu_Pa_s < 1e-3
    # The computational exit station cannot change chamber or throat properties.
    q = chamber_props(pair, 10e6, of, eps=10)
    assert q.c_star_m_s == pytest.approx(p.c_star_m_s, rel=1e-8)
    assert q.T_t_K == pytest.approx(p.T_t_K, rel=1e-8)


@pytest.mark.cea
@pytest.mark.parametrize("pair,of", [("LOX/LH2", 6), ("LOX/RP-1", 2.5), ("LOX/CH4", 3.5)])
def test_frozen_state_consistency(cea_backend, pair, of):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ThermoWarning)
        p = chamber_props(pair, 10e6, of, frozen=True)
    assert p.M_t_kg_kmol == pytest.approx(p.M_kg_kmol, rel=1e-8)
    # Fixed-composition ideal gas: gamma=cp/(cp-R) at each station.
    for tr, mw, gamma in [(p.chamber_transport, p.M_kg_kmol, p.gamma),
                          (p.throat_transport, p.M_t_kg_kmol, p.gamma_t)]:
        r = 8314.462618 / mw
        assert tr.cp_J_kg_K / (tr.cp_J_kg_K - r) == pytest.approx(gamma, rel=0.001)
    # Regression for RocketCEA's transport flow-state issue: frozen throat
    # viscosity must differ from the equilibrium throat's viscosity.
    eq = chamber_props(pair, 10e6, of)
    assert p.throat_transport.mu_Pa_s < eq.throat_transport.mu_Pa_s
    assert ideal_c_star(p.T_c_K, p.M_kg_kmol, p.gamma) == pytest.approx(p.c_star_m_s, rel=0.02)


@pytest.mark.cea
def test_published_rocketcea_example(cea_backend):
    # RocketCEA QuickStart: default LOX/LH2, Pc=100 psia, MR=1, eps=40.
    # https://rocketcea.readthedocs.io/en/latest/quickstart.html
    # 1 psi = 6894.757293168 Pa (international pound/inch, standard gravity).
    with _cea("LOX/LH2") as cea:
        assert _isp(cea, 100 * 6894.757293168, 1, 40, 0) == pytest.approx(374.30361765576265, rel=1e-6)


@pytest.mark.cea
@pytest.mark.parametrize("pair", list(DEFAULT_OF_BOUNDS))
@pytest.mark.parametrize("condition,eps,ambient", [("vacuum", 40, 0), ("sea_level", 10, 101325)])
def test_optimizer_against_dense_grid(cea_backend, pair, condition, eps, ambient):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ThermoWarning)
        best = optimal_OF(pair, 10e6, eps, condition)
    low, high = DEFAULT_OF_BOUNDS[pair]
    assert low < best < high
    with _cea(pair) as cea:
        samples = [_isp(cea, 10e6, of, eps, ambient) for of in np.linspace(low, high, 101)]
        assert _isp(cea, 10e6, best, eps, ambient) >= max(samples) - 1e-5


@pytest.mark.cea
def test_sweep_and_boundary_warning(cea_backend):
    s = sweep_OF("LOX/LH2", 10e6, [6, 5.5])
    assert s.OF.tolist() == [6, 5.5]
    assert s.T_c_K.shape == s.gamma.shape == (2,)
    assert s.c_star_m_s[0] == pytest.approx(chamber_props("LOX/LH2", 10e6, 6).c_star_m_s)
    with pytest.warns(ThermoWarning, match="boundary"):
        assert optimal_OF("LOX/LH2", 10e6, 40, "vacuum", OF_bounds=(2, 3)) == 3
    with pytest.warns(ThermoWarning, match="nominal"):
        chamber_props("LOX/LH2", 10e6, 2)


@pytest.mark.parametrize("pressure,of,eps", [(0, 6, 40), (1e7, -1, 40),
                                              (1e7, 6, 0.9), (1e7, float('inf'), 40)])
def test_invalid_chamber_inputs(pressure, of, eps):
    with pytest.raises(ValueError):
        chamber_props("LOX/LH2", pressure, of, eps)


@pytest.mark.parametrize("samples", [[], [1, -1], [[1, 2]], [float('nan')]])
def test_invalid_sweep(samples):
    with pytest.raises(ValueError):
        sweep_OF("LOX/LH2", 1e7, samples)


def test_invalid_optimizer_and_pair():
    with pytest.raises(ValueError):
        optimal_OF("LOX/LH2", 1e7, 40, "unknown")
    with pytest.raises(ValueError):
        optimal_OF("LOX/LH2", 1e7, 40, "vacuum", OF_bounds=(4, 2))
    with pytest.raises(ValueError):
        chamber_props("unknown", 1e7, 3)


@pytest.mark.cea
def test_published_chamber_and_transport(cea_backend):
    # RocketCEA Standard Examples, LOX/LH2: Pc=1000 psia, MR=6, eps=40.
    # https://rocketcea.readthedocs.io/en/latest/std_examples.html
    # Tolerances cover the rounding in the published output.
    p = chamber_props("LOX/LH2", 1000 * 6894.757293168, 6, 40)
    assert p.T_c_K == pytest.approx(3483.35, abs=0.02)
    assert p.M_kg_kmol == pytest.approx(13.458, abs=0.001)
    assert p.gamma == pytest.approx(1.1401, abs=0.0001)
    assert p.c_star_m_s == pytest.approx(7560.0 * 0.3048, abs=0.05)
    assert p.chamber_transport.mu_Pa_s == pytest.approx(1.0588e-4, abs=1e-8)
