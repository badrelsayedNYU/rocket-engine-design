"""Stage 3: hand-calculated sizing, contour continuity and design API."""
from math import pi, radians, sqrt, tan
import warnings

import numpy as np
import pytest

from engine_design.constants import (RAO_AREA_RATIOS, RAO_ANGLE_TABLE,
    UPSTREAM_THROAT_RADIUS_RATIO, DOWNSTREAM_THROAT_RADIUS_RATIO)
from engine_design.sizing import rao_angles, size_engine

# Analyst-defined inputs, not engine data.
F, CF, PC, CSTAR, OF, EPS, LSTAR = 100e3, 1.85, 10e6, 1840.0, 3.4, 40.0, 1.1


@pytest.fixture(scope="module")
def geom():
    return size_engine(F, CF, PC, CSTAR, OF, EPS, LSTAR, points_per_segment=400)


def test_hand_calculation(geom):
    # At = F/(CF*pc) = 1e5/1.85e7; mdot = pc*At/c*; theory section 3.
    at = 1e5 / 1.85e7
    assert geom.A_t_m2 == pytest.approx(at, rel=1e-14)
    assert geom.D_t_m == pytest.approx(2 * sqrt(at / pi), rel=1e-14)
    assert geom.mdot_kg_s == pytest.approx(1e7 * at / 1840, rel=1e-14)
    assert geom.mdot_ox_kg_s / geom.mdot_fuel_kg_s == pytest.approx(OF)
    assert geom.mdot_ox_kg_s + geom.mdot_fuel_kg_s == pytest.approx(geom.mdot_kg_s)
    assert geom.V_c_m3 == pytest.approx(LSTAR * at)


def test_thrust_recomputed(geom):
    assert CF * PC * geom.A_t_m2 == pytest.approx(F, rel=1e-12)


def test_radii_match_area_ratios(geom):
    rt = geom.D_t_m / 2
    assert geom.r_m.min() == pytest.approx(rt, rel=1e-12)
    assert geom.x_m[np.argmin(geom.r_m)] == pytest.approx(0, abs=1e-12)
    assert (geom.r_m[-1] / rt) ** 2 == pytest.approx(EPS, rel=1e-12)
    assert (geom.r_m[0] / rt) ** 2 == pytest.approx(4.0, rel=1e-12)


def test_contour_continuous_with_smooth_slope(geom):
    segments = geom.segments
    for a, b in zip(segments, segments[1:]):
        assert b.x_m[0] == pytest.approx(a.x_m[-1], abs=1e-12)
        assert b.r_m[0] == pytest.approx(a.r_m[-1], abs=1e-12)
        assert b.slope[0] == pytest.approx(a.slope[-1], abs=1e-9), (a.name, b.name)
    assert np.all(np.diff(geom.x_m) > 0)


def test_wall_angles_and_bell_length(geom):
    bell = geom.segments[-1]
    assert bell.slope[0] == pytest.approx(tan(radians(geom.theta_n_deg)))
    assert bell.slope[-1] == pytest.approx(tan(radians(geom.theta_e_deg)))
    # Throat arc radii, NASA SP-125 Fig. 4-13: curvature of a sampled circle.
    rt = geom.D_t_m / 2
    for segment, ratio in [(geom.segments[3], UPSTREAM_THROAT_RADIUS_RATIO),
                           (geom.segments[4], DOWNSTREAM_THROAT_RADIUS_RATIO)]:
        x0 = 0.0
        r_center = rt + ratio * rt
        distance = np.hypot(segment.x_m - x0, segment.r_m - r_center)
        assert np.allclose(distance, ratio * rt, rtol=1e-12)
    # 15° cone length (Re-Rt)/tan(15°); bell = 80% of it.
    l_cone = (geom.D_e_m - geom.D_t_m) / 2 / tan(radians(15))
    assert geom.L_cone_m == pytest.approx(l_cone)
    assert geom.L_nozzle_m == pytest.approx(0.8 * l_cone)
    assert geom.x_m[-1] == pytest.approx(geom.L_nozzle_m)


def test_volume_independent_check(geom):
    # Exact cylinder and frustum formulas; trapezoid of samples for arcs.
    cyl, fillet, cone, upstream = geom.segments[:4]
    rc = geom.D_c_m / 2
    assert cyl.volume_m3 == pytest.approx(pi * rc**2 * geom.L_cylinder_m, rel=1e-10)
    h = cone.x_m[-1] - cone.x_m[0]
    r1, r2 = cone.r_m[0], cone.r_m[-1]
    assert cone.volume_m3 == pytest.approx(pi * h / 3 * (r1*r1 + r1*r2 + r2*r2), rel=1e-10)
    for arc in (fillet, upstream):
        assert arc.volume_m3 == pytest.approx(np.trapezoid(pi * arc.r_m**2, arc.x_m), rel=1e-5)
    assert sum(s.volume_m3 for s in geom.segments[:4]) == pytest.approx(geom.V_c_m3, rel=1e-10)


def test_rao_table_nodes_and_interpolation():
    for fraction, (theta_n, theta_e) in RAO_ANGLE_TABLE.items():
        for i, eps in enumerate(RAO_AREA_RATIOS):
            assert rao_angles(eps, fraction) == (theta_n[i], theta_e[i])
    # Midpoint in both directions is the mean of four nodes (bilinear).
    n = RAO_ANGLE_TABLE[0.8][0][:2] + RAO_ANGLE_TABLE[0.9][0][:2]
    assert rao_angles(7.5, 0.85)[0] == pytest.approx(sum(n) / 4)
    # Chart trends: theta_n rises with eps, both fall with bell length.
    assert rao_angles(10)[0] < rao_angles(40)[0]
    assert rao_angles(20, 1.0)[1] < rao_angles(20, 0.6)[1]


@pytest.mark.parametrize("eps,fraction", [(4.9, 0.8), (50.1, 0.8), (20, 0.55), (20, 1.01)])
def test_rao_no_extrapolation(eps, fraction):
    with pytest.raises(ValueError):
        rao_angles(eps, fraction)


def test_explicit_angles_outside_chart():
    g = size_engine(F, CF, PC, CSTAR, OF, 80, LSTAR, theta_n_deg=34, theta_e_deg=7)
    assert "User-specified" in g.angle_source
    assert (g.D_e_m / g.D_t_m) ** 2 == pytest.approx(80)
    with pytest.raises(ValueError, match="Rao chart"):
        size_engine(F, CF, PC, CSTAR, OF, 80, LSTAR)


@pytest.mark.parametrize("kwargs", [
    dict(L_star_m=0.05), dict(contraction_ratio=1.0), dict(eps=1.0),
    dict(theta_n_deg=30), dict(theta_n_deg=8, theta_e_deg=10),
    dict(convergence_half_angle_deg=90), dict(points_per_segment=2),
    dict(thrust_N=float("nan")), dict(CF=-1),
])
def test_invalid_sizing(kwargs):
    args = dict(thrust_N=F, CF=CF, p_c_Pa=PC, c_star_m_s=CSTAR, OF=OF, eps=EPS, L_star_m=LSTAR)
    args.update(kwargs)
    with pytest.raises(ValueError):
        size_engine(**args)


# ---- design_engine: requires CEA ----

@pytest.fixture
def cea_backend():
    pytest.importorskip("rocketcea")


@pytest.fixture(scope="module")
def spec_design():
    pytest.importorskip("rocketcea")
    from engine_design import EngineRequirements, design_engine
    return design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                                            thrust_condition="vacuum", p_c_Pa=10e6))


@pytest.mark.cea
def test_design_meets_requirement(spec_design):
    d = spec_design
    # The throat is sized so DELIVERED vacuum thrust meets the requirement.
    assert d.thrust_vac_N == pytest.approx(100e3, rel=1e-12)
    # Matched nozzle: exit pressure equals ambient, so no pressure term.
    assert d.ideal.p_e_Pa == pytest.approx(101325, rel=1e-3)
    assert d.ideal.CF_vac - d.ideal.CF_amb == pytest.approx(d.eps * 101325 / 10e6)
    assert not d.likely_separation_at_p_amb
    assert 1.7 <= d.ideal.CF_vac <= 2.0  # theory section 6 sanity band
    assert d.ideal.Isp_vac_s == pytest.approx(d.ideal.c_star_m_s * d.ideal.CF_vac / 9.80665)


@pytest.mark.cea
def test_design_optimum_is_self_consistent(spec_design):
    # At the coupled fixed point, O/F maximizes vacuum Isp at the final eps,
    # and that eps is the matched eps of the final O/F.
    from engine_design.thermo import ThermoWarning, matched_eps, optimal_OF
    d = spec_design
    assert matched_eps("LOX/CH4", 10e6, d.OF, 101325) == pytest.approx(d.eps, rel=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ThermoWarning)
        assert optimal_OF("LOX/CH4", 10e6, d.eps, "vacuum") == pytest.approx(d.OF, abs=1e-3)


@pytest.mark.cea
def test_cea_cf_against_constant_gamma_frozen(cea_backend):
    # Independent model cross-check: frozen CEA CF vs the Stage 2 formula
    # with the frozen chamber gamma. Variable cp gives ~0.5% here (eps=13).
    from engine_design.nozzle import thrust_coefficient
    from engine_design.thermo import ThermoWarning, chamber_props, expansion_performance
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ThermoWarning)
        gamma = chamber_props("LOX/CH4", 10e6, 3.4, frozen=True).gamma
    perf = expansion_performance("LOX/CH4", 10e6, 3.4, 13, frozen=True)
    assert perf.CF_vac == pytest.approx(thrust_coefficient(gamma, 13, 10e6, 0), rel=0.01)


@pytest.mark.cea
def test_fixed_eps_sea_level_separation_warns(cea_backend):
    from engine_design import EngineRequirements, design_engine
    from engine_design.nozzle import NozzleWarning
    req = EngineRequirements(propellants="LOX/RP-1", thrust_N=500e3, thrust_condition="sea_level",
                             p_c_Pa=5e6, mixture_ratio=2.5, nozzle=40.0)
    with pytest.warns(NozzleWarning, match="separation"):
        d = design_engine(req)
    assert d.likely_separation_at_p_amb
    assert d.thrust_amb_N == pytest.approx(500e3, rel=1e-12)
    assert d.OF == 2.5 and d.eps == 40.0


@pytest.mark.cea
def test_fixed_eps_vacuum_and_summary(cea_backend):
    from engine_design import EngineRequirements, design_engine
    req = EngineRequirements(propellants="LOX/LH2", thrust_N=110e3, p_c_Pa=4.4e6,
                             mixture_ratio=5.9, nozzle=80.0, theta_n_deg=34, theta_e_deg=7,
                             divergence_factor=0.99)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # vacuum design: no separation warning
        d = design_engine(req)
    assert d.thrust_vac_N == pytest.approx(110e3, rel=1e-12)
    text = d.summary()
    for label in ("O/F", "Throat / exit diameter", "Isp vacuum", "User-specified"):
        assert label in text


@pytest.mark.parametrize("kwargs", [
    dict(propellants="LOX/N2"), dict(thrust_condition="orbit"), dict(nozzle="optimal"),
    dict(nozzle=0.5), dict(nozzle=True), dict(mixture_ratio=-1), dict(p_amb_Pa=2e7),
    dict(thrust_N=0),
])
def test_invalid_requirements(kwargs):
    from engine_design import EngineRequirements, design_engine
    args = dict(propellants="LOX/CH4", thrust_N=100e3, p_c_Pa=10e6)
    args.update(kwargs)
    with pytest.raises(ValueError):
        design_engine(EngineRequirements(**args))
