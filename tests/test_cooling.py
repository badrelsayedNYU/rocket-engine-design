"""Stage 5: Bartz heat flux, heat load and coolant capacity."""
from math import pi, sqrt
import re

import numpy as np
import pytest

from engine_design.constants import COOLANT_INLET_T_K
from engine_design.cooling import (adiabatic_wall_temperature, bartz_h, bartz_sigma,
                                   coolant_capacity, coolant_enthalpy_rise,
                                   heat_flux_profile, wall_heat_load)
from engine_design.thermo import TransportProperties


def test_bartz_equals_turbulent_nusselt_form():
    # Unit check for the SI coefficient: with G = p_c/c*, k = mu*cp/Pr and
    # Re = G*D/mu, Bartz's group is (k/D)*0.026*Re^0.8*Pr^0.4*(D/r_c)^0.1.
    # Analyst-defined gas properties, not CEA output.
    D, mu, cp, Pr, pc, cstar, rc = 0.08, 1.1e-4, 2400.0, 0.67, 10e6, 1850.0, 0.04
    G = pc / cstar
    k = mu * cp / Pr
    nusselt = k / D * 0.026 * (G * D / mu) ** 0.8 * Pr ** 0.4 * (D / rc) ** 0.1
    assert bartz_h(D, mu, cp, Pr, pc, cstar, rc, 1.0, 1.0) == pytest.approx(nusselt, rel=1e-12)
    # Local scaling (A_t/A)^0.9 and linear in sigma.
    assert bartz_h(D, mu, cp, Pr, pc, cstar, rc, 4.0, 1.3) == pytest.approx(nusselt * 4 ** -0.9 * 1.3)


def test_sigma_hand_calculation():
    assert bartz_sigma(1.0, 0.0, 1.2) == pytest.approx(1.0)
    # gamma=1.2, M=1: stagnation = 1.1; Tw/Tc = 0.25: film = 0.5*0.25*1.1+0.5 = 0.6375.
    expected = 0.6375 ** -0.68 * 1.1 ** -0.12
    assert bartz_sigma(0.25, 1.0, 1.2) == pytest.approx(expected, rel=1e-12)
    # A cool wall raises sigma (denser, more viscous film), and it falls with Mach.
    assert bartz_sigma(0.25, 0.2, 1.2) > bartz_sigma(0.25, 1.0, 1.2) > bartz_sigma(0.25, 3.0, 1.2)


def test_adiabatic_wall_temperature_limits():
    assert adiabatic_wall_temperature(3500.0, 0.0, 1.2, 0.7) == pytest.approx(3500.0)
    assert adiabatic_wall_temperature(3500.0, 3.0, 1.2, 1.0) == pytest.approx(3500.0)
    # gamma=1.2, M=2: rise=0.4; Pr=0.729 -> r=0.9; T_aw = 3500*(1.36/1.4).
    assert adiabatic_wall_temperature(3500.0, 2.0, 1.2, 0.729) == pytest.approx(3400.0)


def test_heat_load_is_exact_on_a_frustum():
    # Uniform q=1 over a cone frustum r: 1 -> 3 over x: 0 -> 2 gives the
    # lateral area pi*(r1+r2)*slant = 4*pi*sqrt(8).
    x = np.linspace(0, 2, 7)
    r = 1 + x
    assert wall_heat_load(x, r, np.ones_like(x)) == pytest.approx(4 * pi * sqrt(8), rel=1e-12)


def _profile(T_w=800.0):
    # Analyst-defined converging-diverging duct; throat at x=0.
    x = np.linspace(-0.1, 0.2, 301)
    r = 0.04 + np.where(x < 0, 2.0, 0.6) * x ** 2
    tr = TransportProperties(2400.0, 1.1e-4, 0.4, 0.66)
    return heat_flux_profile(x, r, 3500.0, 1.2, tr, 10e6, 1850.0, T_w)


def test_profile_branches_and_peak_location():
    p = _profile()
    # Exclude round-off samples at the throat (linspace leaves x ~ 1e-17).
    assert np.all(p.mach[p.x_m < -1e-12] < 1) and np.all(p.mach[p.x_m > 1e-12] > 1)
    assert p.mach[np.argmin(np.abs(p.x_m))] == pytest.approx(1.0)
    # sigma and T_aw are larger subsonic, so the peak is at or just upstream.
    assert -0.01 <= p.x_peak_m <= 0
    assert np.all(p.q_W_m2 > 0)


def test_hotter_wall_lowers_flux_and_invalid_wall():
    assert _profile(1000.0).peak_q_W_m2 < _profile(600.0).peak_q_W_m2
    with pytest.raises(ValueError):
        _profile(4000.0)


# ---- Coolant: requires CoolProp ----

@pytest.mark.parametrize("fluid,T_in,T_out,nist_kJ_kg", [
    # NIST Chemistry WebBook fluid data, 10 MPa isobar (accessed 2026-10-21):
    # parahydrogen h(20 K)=98.9190, h(800 K)=11834.7 kJ/kg;
    # methane h(100 K)=-24.4783, h(600 K)=1712.17 kJ/kg.
    ("ParaHydrogen", 20.0, 800.0, 11834.7 - 98.9190),
    ("Methane", 100.0, 600.0, 1712.17 + 24.4783),
])
def test_coolant_enthalpy_matches_nist(fluid, T_in, T_out, nist_kJ_kg):
    pytest.importorskip("CoolProp")
    assert coolant_enthalpy_rise(fluid, T_in, T_out, 10e6) / 1e3 == pytest.approx(nist_kJ_kg, abs=0.1)


def test_coolant_caps_and_rp1():
    pytest.importorskip("CoolProp")
    h2, _ = coolant_capacity("LOX/LH2", 5.0, 10e6, 800.0)
    assert h2.T_out_K == 800.0 and h2.T_in_K == 20.27
    assert h2.capacity_W == pytest.approx(5.0 * h2.delta_h_J_kg)
    ch4, _ = coolant_capacity("LOX/CH4", 7.0, 10e6, 800.0)
    assert ch4.T_out_K == 625.0 and "EOS" in ch4.T_out_basis
    assert coolant_capacity("LOX/CH4", 7.0, 10e6, 500.0)[0].T_out_K == 500.0
    with pytest.raises(ValueError):
        coolant_capacity("LOX/CH4", 7.0, 10e6, 800.0, T_out_K=700.0)
    # RP-1: no verified c_p is set, so no capacity unless supplied.
    none, note = coolant_capacity("LOX/RP-1", 9.0, 10e6, 800.0)
    assert none is None and "c_p" in note
    rp1, _ = coolant_capacity("LOX/RP-1", 9.0, 10e6, 800.0, cp_J_kg_K=2000.0)
    assert rp1.T_out_K == 600.0
    assert rp1.capacity_W == pytest.approx(9.0 * 2000.0 * (600.0 - 298.15))


def test_inlet_temperatures_match_cea_cards():
    pytest.importorskip("rocketcea")
    from rocketcea import cea_obj
    from engine_design.propellants import PROPELLANTS
    for name, T in COOLANT_INLET_T_K.items():
        card = " ".join(cea_obj.fuelCards[PROPELLANTS[name].cea_fuel])
        assert float(re.search(r"t\(k\)=\s*([\d.]+)", card).group(1)) == T


# ---- design_engine: requires CEA and CoolProp ----

@pytest.fixture(scope="module")
def spec_design():
    pytest.importorskip("rocketcea")
    pytest.importorskip("CoolProp")
    from engine_design import EngineRequirements, design_engine
    return design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                                            thrust_condition="vacuum", p_c_Pa=10e6))


@pytest.mark.cea
def test_design_peak_location_and_band(spec_design):
    p, rt = spec_design.cooling.profile, spec_design.geometry.D_t_m / 2
    # Stage 5 spec: peak at or slightly upstream of the throat.
    assert -0.25 * rt <= p.x_peak_m <= 0
    # Theory section 6: tens of MW/m² for high-p_c engines.
    assert 10e6 < p.peak_q_W_m2 < 200e6


@pytest.mark.cea
def test_design_uses_sized_mass_flux_and_frozen_properties(spec_design):
    from engine_design.thermo import chamber_props
    d, g = spec_design, spec_design.geometry
    # p_c/c*_delivered is the sized throat mass flux.
    assert d.requirements.p_c_Pa / d.delivered.c_star_m_s == pytest.approx(g.mdot_kg_s / g.A_t_m2)
    frozen = chamber_props("LOX/CH4", 10e6, d.OF, d.eps, frozen=True)
    p = d.cooling.profile
    i = int(np.argmin(np.abs(p.x_m)))
    tr = frozen.chamber_transport
    sigma = bartz_sigma(800.0 / frozen.T_c_K, 1.0, frozen.gamma)
    rt = g.D_t_m / 2
    expected = bartz_h(g.D_t_m, tr.mu_Pa_s, tr.cp_J_kg_K, tr.Pr, 10e6, d.delivered.c_star_m_s,
                       (1.5 + 0.382) / 2 * rt, 1.0, sigma)
    assert p.h_g_W_m2_K[i] == pytest.approx(expected, rel=1e-9)


@pytest.mark.cea
def test_design_capacity_ratio(spec_design):
    c = spec_design.cooling
    assert c.capacity_ratio == pytest.approx(c.coolant.capacity_W / c.heat_load_W)
    assert c.coolant.mdot_kg_s == spec_design.geometry.mdot_fuel_kg_s
    assert "Cooling" in spec_design.summary()


@pytest.mark.cea
def test_design_rp1_without_cp_reports_missing():
    pytest.importorskip("CoolProp")
    from engine_design import EngineRequirements, design_engine
    d = design_engine(EngineRequirements(propellants="LOX/RP-1", thrust_N=100e3, p_c_Pa=10e6))
    assert d.cooling.coolant is None and d.cooling.capacity_ratio is None
    assert d.cooling.heat_load_W > 0
