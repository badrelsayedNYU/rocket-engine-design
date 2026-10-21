"""Stage 5: Bartz gas-side heat flux along the wall and a regen capacity check.

Gas-side convection only: no wall conduction, coolant-side film, radiation,
film cooling or channel design. Bartz typically OVERPREDICTS throat heat flux
somewhat; treat the peak as a conservative first estimate, not a prediction
to design channels against. References: theory section 5; Bartz (1957),
Jet Propulsion 27(1); Huzel & Huang Ch. 4; Sutton & Biblarz Ch. 8.
"""
from dataclasses import dataclass
from math import isfinite, pi
import warnings

import numpy as np

from .constants import (BARTZ_COEFFICIENT, COOLANT_INLET_T_K, COOLPROP_FLUID,
                        DOWNSTREAM_THROAT_RADIUS_RATIO, HEAT_FLUX_SANITY_W_M2,
                        METHANE_EOS_T_MAX_K, RP1_CP_J_KG_K, RP1_DECOMPOSITION_LIMIT_K,
                        UPSTREAM_THROAT_RADIUS_RATIO)
from .nozzle import mach_from_area_ratio


class CoolingWarning(UserWarning):
    """A heat-transfer result falls outside its expected band."""


@dataclass(frozen=True)
class HeatFluxProfile:
    """Gas-side quantities at each wall sample, SI; throat at x=0."""
    x_m: np.ndarray
    r_m: np.ndarray
    mach: np.ndarray
    sigma: np.ndarray
    h_g_W_m2_K: np.ndarray
    T_aw_K: np.ndarray
    q_W_m2: np.ndarray
    T_w_K: float

    @property
    def peak_index(self):
        return int(np.argmax(self.q_W_m2))

    @property
    def peak_q_W_m2(self):
        return float(self.q_W_m2[self.peak_index])

    @property
    def x_peak_m(self):
        return float(self.x_m[self.peak_index])


@dataclass(frozen=True)
class CoolantCapacity:
    """Heat the fuel can absorb between inlet and allowable outlet, SI."""
    coolant: str
    p_Pa: float
    T_in_K: float
    T_out_K: float
    T_out_basis: str
    delta_h_J_kg: float
    mdot_kg_s: float

    @property
    def capacity_W(self):
        """Q_cap = mdot_fuel·(h_out - h_in) (theory section 5)."""
        return self.mdot_kg_s * self.delta_h_J_kg


@dataclass(frozen=True)
class CoolingResult:
    """Heat-flux profile, total wall heat load and fuel coolant capacity."""
    profile: HeatFluxProfile
    heat_load_W: float
    coolant: CoolantCapacity | None
    coolant_note: str

    @property
    def capacity_ratio(self):
        """Q_cap/Q. Above 1: the fuel can absorb the load within its limit.
        None when the coolant capacity is not available."""
        return None if self.coolant is None else self.coolant.capacity_W / self.heat_load_W


def _positive(name, value):
    if isinstance(value, bool) or not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _gamma(gamma):
    if not isfinite(gamma) or gamma <= 1:
        raise ValueError("gamma must be finite and > 1")


def bartz_sigma(T_w_over_T_c, M, gamma, omega=0.6):
    """Bartz property-variation correction (theory section 5):

    sigma = 1 / ([0.5·(Tw/Tc)·(1+(γ-1)/2·M²) + 0.5]^(0.8-ω/5)
                 · [1+(γ-1)/2·M²]^(ω/5)),
    with viscosity-temperature exponent ω=0.6, giving 0.68 and 0.12.
    Evaluates film properties at the mean of wall and free-stream static
    temperature. sigma=1 at Tw=Tc and M=0. Bartz (1957).
    """
    _positive("T_w_over_T_c", T_w_over_T_c)
    _gamma(gamma)
    stagnation = 1 + (gamma - 1) / 2 * np.asarray(M, dtype=float) ** 2
    film = 0.5 * T_w_over_T_c * stagnation + 0.5
    return film ** -(0.8 - omega / 5) * stagnation ** -(omega / 5)


def adiabatic_wall_temperature(T_c_K, M, gamma, Pr):
    """T_aw = T_c·(1 + r·(γ-1)/2·M²)/(1 + (γ-1)/2·M²), r = Pr^(1/3).

    Turbulent recovery factor; equals T_c at M=0 or r=1. Theory section 5;
    Huzel & Huang Ch. 4.
    """
    _positive("T_c_K", T_c_K)
    _positive("Pr", Pr)
    _gamma(gamma)
    rise = (gamma - 1) / 2 * np.asarray(M, dtype=float) ** 2
    return T_c_K * (1 + Pr ** (1 / 3) * rise) / (1 + rise)


def bartz_h(D_t_m, mu_Pa_s, cp_J_kg_K, Pr, p_c_Pa, c_star_m_s, r_c_m, area_ratio, sigma):
    """Bartz gas-side coefficient in SI (theory section 5):

    h_g = 0.026/D_t^0.2 · (μ^0.2·c_p/Pr^0.6)_c · (p_c/c*)^0.8
          · (D_t/r_c)^0.1 · (A_t/A)^0.9 · σ,   area_ratio = A/A_t.

    Turbulent pipe-flow form Nu = 0.026·Re^0.8·Pr^0.4 based on throat mass
    flux p_c/c* and diameter, scaled by local area. Properties at chamber
    stagnation conditions. Known to overpredict throat heat flux somewhat.
    """
    for name, value in dict(D_t_m=D_t_m, mu_Pa_s=mu_Pa_s, cp_J_kg_K=cp_J_kg_K, Pr=Pr,
                            p_c_Pa=p_c_Pa, c_star_m_s=c_star_m_s, r_c_m=r_c_m).items():
        _positive(name, value)
    group = (BARTZ_COEFFICIENT / D_t_m ** 0.2 * mu_Pa_s ** 0.2 * cp_J_kg_K / Pr ** 0.6
             * (p_c_Pa / c_star_m_s) ** 0.8 * (D_t_m / r_c_m) ** 0.1)
    return group * np.asarray(area_ratio, dtype=float) ** -0.9 * sigma


def heat_flux_profile(x_m, r_m, T_c_K, gamma, transport, p_c_Pa, c_star_m_s, T_w_K):
    """Return HeatFluxProfile q = h_g·(T_aw - T_w) at every wall sample.

    A/A_t = (r/R_t)² with R_t = min(r); constant-gamma Mach on the subsonic
    branch for x<0 and supersonic for x>=0 (throat at x=0). r_c is the mean
    of the upstream and downstream throat radii, (1.5+0.382)/2·R_t, as the
    spec states. The 1D Mach from area is poor near the injector, but h_g
    there is small. T_w is uniform and assumed, not solved.
    """
    x, r = np.asarray(x_m, dtype=float), np.asarray(r_m, dtype=float)
    if x.shape != r.shape or x.ndim != 1 or np.any(np.diff(x) <= 0) or np.any(r <= 0):
        raise ValueError("x_m must increase and r_m must be positive, same 1D shape")
    _positive("T_w_K", T_w_K)
    rt = float(r.min())
    area = np.maximum((r / rt) ** 2, 1.0)  # Guard round-off below 1 at the throat.
    mach = np.array([mach_from_area_ratio(a, gamma, "subsonic" if xi < 0 else "supersonic")
                     for a, xi in zip(area, x)])
    T_aw = adiabatic_wall_temperature(T_c_K, mach, gamma, transport.Pr)
    if np.any(T_aw <= T_w_K):
        raise ValueError("Wall temperature must lie below the adiabatic wall temperature")
    sigma = bartz_sigma(T_w_K / T_c_K, mach, gamma)
    r_c = (UPSTREAM_THROAT_RADIUS_RATIO + DOWNSTREAM_THROAT_RADIUS_RATIO) / 2 * rt
    h = bartz_h(2 * rt, transport.mu_Pa_s, transport.cp_J_kg_K, transport.Pr,
                p_c_Pa, c_star_m_s, r_c, area, sigma)
    arrays = (x, r, mach, sigma, h, T_aw, h * (T_aw - T_w_K))
    for a in arrays:
        a.flags.writeable = False
    return HeatFluxProfile(*arrays, float(T_w_K))


def wall_heat_load(x_m, r_m, q_W_m2):
    """Q = ∫ q dA, dA = 2πr ds over the wetted side wall (theory section 5).

    Trapezoidal rule on straight chords between samples, so a frustum's
    lateral area is exact. The injector face is not included.
    """
    x, r, q = (np.asarray(a, dtype=float) for a in (x_m, r_m, q_W_m2))
    ds = np.hypot(np.diff(x), np.diff(r))
    integrand = 2 * pi * r * q
    return float(np.sum(ds * (integrand[1:] + integrand[:-1]) / 2))


def coolant_enthalpy_rise(fluid, T_in_K, T_out_K, p_Pa):
    """h(T_out, p) - h(T_in, p) from the CoolProp reference EOS, J/kg.

    Isobaric; the jacket pressure drop is neglected. Real-fluid enthalpy
    captures the large c_p peak near the critical point, which a constant
    c_p cannot.
    """
    try:
        from CoolProp.CoolProp import PropsSI
    except ImportError as exc:
        raise ImportError("Install CoolProp for coolant properties: pip install -e .") from exc
    if not T_in_K < T_out_K:
        raise ValueError("Coolant outlet temperature must exceed inlet temperature")
    return float(PropsSI("H", "T", T_out_K, "P", p_Pa, fluid)
                 - PropsSI("H", "T", T_in_K, "P", p_Pa, fluid))


def coolant_capacity(propellants, mdot_fuel_kg_s, p_Pa, T_w_K, *, T_out_K=None, cp_J_kg_K=None):
    """Return (CoolantCapacity or None, note) for the fuel as sole coolant.

    The bulk coolant cannot leave hotter than the wall it cools, so T_w is
    the default outlet cap for every fuel; the capacity is then an UPPER
    bound. Methane is also capped at its EOS limit (625 K) and RP-1 at its
    ~600 K decomposition limit (constants.py). No verified methane coking
    limit is encoded. LH2 and CH4 use CoolProp enthalpy; RP-1 uses a
    constant c_p, which must be supplied because none is set by default.
    p_Pa is the jacket pressure, by default p_c (a lower bound).
    """
    _positive("mdot_fuel_kg_s", mdot_fuel_kg_s)
    _positive("p_Pa", p_Pa)
    _positive("T_w_K", T_w_K)
    T_in = COOLANT_INLET_T_K[propellants]
    if propellants == "LOX/CH4":
        cap, basis = min(T_w_K, METHANE_EOS_T_MAX_K), (
            "wall temperature" if T_w_K <= METHANE_EOS_T_MAX_K else "methane EOS limit 625 K")
    elif propellants == "LOX/RP-1":
        cap, basis = min(T_w_K, RP1_DECOMPOSITION_LIMIT_K), (
            "wall temperature" if T_w_K <= RP1_DECOMPOSITION_LIMIT_K
            else "RP-1 decomposition ~600 K (NASA CR-165405)")
    else:
        cap, basis = T_w_K, "wall temperature"
    if T_out_K is None:
        T_out = cap
    else:
        _positive("T_out_K", T_out_K)
        if T_out_K > cap:
            raise ValueError(f"Coolant outlet {T_out_K:g} K exceeds the {basis} cap of {cap:g} K")
        T_out, basis = T_out_K, "user-specified"
    if propellants in COOLPROP_FLUID:
        fluid = COOLPROP_FLUID[propellants]
        dh = coolant_enthalpy_rise(fluid, T_in, T_out, p_Pa)
        note = f"CoolProp {fluid}, isobaric at {p_Pa/1e6:g} MPa"
    else:
        cp = RP1_CP_J_KG_K if cp_J_kg_K is None else cp_J_kg_K
        if cp is None:
            return None, ("RP-1 c_p not set: no verified source accessible "
                          "(constants.py); supply coolant_cp_J_kg_K")
        _positive("coolant c_p", cp)
        if not T_in < T_out:
            raise ValueError("Coolant outlet temperature must exceed inlet temperature")
        dh = cp * (T_out - T_in)
        fluid, note = "RP-1", f"constant c_p = {cp:g} J/(kg K), user-specified"
    return CoolantCapacity(fluid, p_Pa, T_in, T_out, basis, dh, mdot_fuel_kg_s), note


def cooling_analysis(propellants, geometry, chamber_frozen, p_c_Pa, c_star_m_s, T_w_K, *,
                     coolant_p_Pa=None, T_out_K=None, cp_J_kg_K=None):
    """Bartz heat flux along the sized contour, heat load and coolant ratio.

    Uses FROZEN chamber properties (T_c, γ, μ, c_p, Pr): equilibrium c_p
    includes reaction heat, about 2.7x the frozen value for LOX/CH4 at
    10 MPa, and would triple h_g. Frozen properties keep the correlation in
    the non-reacting gas form it was fitted in. c_star_m_s should be the
    delivered c*, so p_c/c* equals the sized throat mass flux mdot/A_t.
    Warns if the peak lies outside the spec's sanity band.
    """
    transport = chamber_frozen.chamber_transport
    profile = heat_flux_profile(geometry.x_m, geometry.r_m, chamber_frozen.T_c_K,
                                chamber_frozen.gamma, transport, p_c_Pa, c_star_m_s, T_w_K)
    low, high = HEAT_FLUX_SANITY_W_M2
    if not low <= profile.peak_q_W_m2 <= high:
        warnings.warn(f"Peak heat flux {profile.peak_q_W_m2/1e6:.1f} MW/m² is outside the "
                      f"spec band {low/1e6:g}–{high/1e6:g} MW/m²; investigate.",
                      CoolingWarning, stacklevel=2)
    load = wall_heat_load(profile.x_m, profile.r_m, profile.q_W_m2)
    coolant, note = coolant_capacity(propellants, geometry.mdot_fuel_kg_s,
                                     p_c_Pa if coolant_p_Pa is None else coolant_p_Pa, T_w_K,
                                     T_out_K=T_out_K, cp_J_kg_K=cp_J_kg_K)
    return CoolingResult(profile, load, coolant, note)
