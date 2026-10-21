"""Top-level preliminary design: requirements -> EngineDesign.

Ideal performance comes from CEA equilibrium (shifting) expansion, with
CF = Isp*g0/c*; the same objective selects O/F. Delivered performance applies
lambda, eta_c* and eta_CF (losses.py), and the throat is sized so the
DELIVERED thrust meets the requirement. Frozen expansion is reported for
comparison only. Cooling is a later stage.
"""
from dataclasses import dataclass
from math import isfinite
import warnings

from .constants import (CHARACTERISTIC_LENGTH_M, CHARACTERISTIC_LENGTH_SOURCE,
    DEFAULT_BELL_FRACTION, DEFAULT_CONTRACTION_RATIO, DEFAULT_CONVERGENCE_ANGLE_DEG,
    STANDARD_ATMOSPHERE_PA, CEA_DEFAULT_EPS, COMBUSTION_EFFICIENCY, NOZZLE_EFFICIENCY)
from .losses import Performance, apply_losses, bell_divergence_factor, ideal_performance
from .nozzle import NozzleWarning, separation_check
from .propellants import PROPELLANTS
from .sizing import EngineGeometry, size_engine
from .thermo import (ChamberProperties, ThermoWarning, chamber_props,
                     expansion_performance, matched_eps, optimal_OF)

SEA_LEVEL_OPTIMIZED = "sea_level_optimized"
# Numerical tolerance and iteration cap for the coupled O/F–eps search.
_OF_TOLERANCE = 1e-4
_MAX_COUPLED_ITERATIONS = 30


@dataclass(frozen=True)
class EngineRequirements:
    """Top-level requirements in SI. Fields mirror the project brief API.

    nozzle is 'sea_level_optimized' (p_e = p_amb_Pa) or a fixed area ratio.
    mixture_ratio=None maximizes Isp at thrust_condition. thrust_condition
    'sea_level' means ambient pressure p_amb_Pa. L_star_m=None uses the
    propellant default in constants.py. Wall angles are required outside
    the digitized Rao chart (5 <= eps <= 50). divergence_factor=None uses
    SP-125 Fig.4-12 for chart-angle bells (10 <= eps <= 40); otherwise supply
    it. Efficiencies default to the documented assumptions in constants.py;
    set all three to 1 for an ideal design.
    """
    propellants: str
    thrust_N: float
    thrust_condition: str = "vacuum"
    p_c_Pa: float = 10e6
    mixture_ratio: float | None = None
    nozzle: str | float = SEA_LEVEL_OPTIMIZED
    p_amb_Pa: float = STANDARD_ATMOSPHERE_PA
    L_star_m: float | None = None
    contraction_ratio: float = DEFAULT_CONTRACTION_RATIO
    bell_fraction: float = DEFAULT_BELL_FRACTION
    convergence_half_angle_deg: float = DEFAULT_CONVERGENCE_ANGLE_DEG
    theta_n_deg: float | None = None
    theta_e_deg: float | None = None
    divergence_factor: float | None = None
    eta_c_star: float = COMBUSTION_EFFICIENCY
    eta_CF: float = NOZZLE_EFFICIENCY


@dataclass(frozen=True)
class EngineDesign:
    """Design point with ideal, frozen and delivered performance, all SI.

    ideal: CEA shifting equilibrium, no losses. ideal_frozen: CEA with
    composition frozen at the chamber, no losses. delivered: ideal with
    lambda, eta_c* and eta_CF; the geometry is sized from delivered values.
    """
    requirements: EngineRequirements
    OF: float
    OF_source: str
    eps: float
    eps_source: str
    chamber: ChamberProperties
    ideal: Performance
    ideal_frozen: Performance
    delivered: Performance
    divergence_factor: float
    divergence_source: str
    likely_separation_at_p_amb: bool
    L_star_source: str
    geometry: EngineGeometry

    @property
    def thrust_vac_N(self):
        """Delivered F = CF_del·p_c·A_t in vacuum."""
        return self.delivered.CF_vac * self.requirements.p_c_Pa * self.geometry.A_t_m2

    @property
    def thrust_amb_N(self):
        """Delivered F at p_amb (attached-flow extrapolation if separated)."""
        return self.delivered.CF_amb * self.requirements.p_c_Pa * self.geometry.A_t_m2

    def summary(self):
        """Return a plain-text table of the design point."""
        r, g, c = self.requirements, self.geometry, self.chamber
        amb = f"{r.p_amb_Pa/1e3:g} kPa"
        perf = (self.ideal_frozen, self.ideal, self.delivered)
        rows = [
            ("Requirements", None),
            ("Propellants", r.propellants),
            ("Thrust requirement", f"{r.thrust_N/1e3:.3f} kN ({r.thrust_condition})"),
            ("Chamber pressure", f"{r.p_c_Pa/1e6:.3f} MPa"),
            ("Ambient pressure", amb),
            ("Combustion (CEA equilibrium)", None),
            ("O/F (mass)", f"{self.OF:.4f}  [{self.OF_source}]"),
            ("Chamber temperature", f"{c.T_c_K:.1f} K"),
            ("Molecular weight", f"{c.M_kg_kmol:.3f} kg/kmol"),
            ("Isentropic exponent", f"{c.gamma:.4f}"),
            ("c* ideal", f"{self.ideal.c_star_m_s:.1f} m/s"),
            ("Nozzle (attached flow)", None),
            ("Area ratio eps", f"{self.eps:.3f}  [{self.eps_source}]"),
            ("Exit pressure (shifting)", f"{self.ideal.p_e_Pa/1e3:.3f} kPa"),
            (f"Separation flag at {amb}", "LIKELY (Summerfield)" if self.likely_separation_at_p_amb else "no"),
            ("Performance", None),
            ("", f"{'frozen':>9} {'shifting':>9} {'delivered':>10}"),
            ("c* [m/s]", "  ".join(f"{p.c_star_m_s:8.1f}" for p in perf) + " "),
            ("CF vacuum", "  ".join(f"{p.CF_vac:8.4f}" for p in perf) + " "),
            ("CF ambient", "  ".join(f"{p.CF_amb:8.4f}" for p in perf) + " "),
            ("Isp vacuum [s]", "  ".join(f"{p.Isp_vac_s:8.1f}" for p in perf) + " "),
            ("Isp ambient [s]", "  ".join(f"{p.Isp_amb_s:8.1f}" for p in perf) + " "),
            ("Divergence factor lambda", f"{self.divergence_factor:.4f}  [{self.divergence_source}]"),
            ("eta_c* / eta_CF", f"{r.eta_c_star:.3f} / {r.eta_CF:.3f}"),
            ("Delivered thrust vac / amb", f"{self.thrust_vac_N/1e3:.2f} / {self.thrust_amb_N/1e3:.2f} kN"),
            ("Mass flow", None),
            ("Total", f"{g.mdot_kg_s:.3f} kg/s"),
            ("Oxidizer / fuel", f"{g.mdot_ox_kg_s:.3f} / {g.mdot_fuel_kg_s:.3f} kg/s"),
            ("Geometry", None),
            ("Throat / exit diameter", f"{g.D_t_m*1e3:.2f} / {g.D_e_m*1e3:.2f} mm"),
            ("Chamber diameter (CR)", f"{g.D_c_m*1e3:.2f} mm  (CR={g.A_c_m2/g.A_t_m2:g})"),
            ("L* / chamber volume", f"{g.L_star_m:.3f} m / {g.V_c_m3*1e6:.1f} cm³"),
            ("L* source", self.L_star_source),
            ("Cylinder / converging length", f"{g.L_cylinder_m*1e3:.1f} / {g.L_converging_m*1e3:.1f} mm"),
            ("Bell length (vs 15° cone)", f"{g.L_nozzle_m*1e3:.1f} mm ({g.L_nozzle_m/g.L_cone_m:.0%} of {g.L_cone_m*1e3:.1f} mm)"),
            ("Wall angles theta_n / theta_e", f"{g.theta_n_deg:.2f}° / {g.theta_e_deg:.2f}°"),
            ("Angle source", g.angle_source),
            ("Overall length (injector to exit)", f"{(g.L_chamber_m + g.L_nozzle_m)*1e3:.1f} mm"),
        ]
        width = max(len(k) for k, v in rows if v is not None) + 2
        lines = []
        for key, value in rows:
            lines.append(f"\n{key}" if value is None else f"  {key:<{width}}{value}")
        return "\n".join(lines).lstrip("\n")

    def plot_contour(self, ax=None):
        """Draw the wall contour to scale with key dimensions; return figure."""
        from .plots import plot_contour
        return plot_contour(self, ax=ax)


def _validate(req):
    if req.propellants not in PROPELLANTS:
        raise ValueError(f"Unknown propellants {req.propellants!r}; choose {tuple(PROPELLANTS)}")
    if req.thrust_condition not in ("vacuum", "sea_level"):
        raise ValueError("thrust_condition must be 'vacuum' or 'sea_level'")
    for name in ("thrust_N", "p_c_Pa", "p_amb_Pa"):
        value = getattr(req, name)
        if not isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if req.p_amb_Pa >= req.p_c_Pa:
        raise ValueError("p_amb_Pa must be below p_c_Pa")
    if isinstance(req.nozzle, str):
        if req.nozzle != SEA_LEVEL_OPTIMIZED:
            raise ValueError(f"nozzle must be {SEA_LEVEL_OPTIMIZED!r} or a numeric area ratio")
    elif isinstance(req.nozzle, bool) or not isfinite(req.nozzle) or req.nozzle <= 1:
        raise ValueError("A fixed nozzle area ratio must be finite and > 1")
    if req.mixture_ratio is not None and (not isfinite(req.mixture_ratio) or req.mixture_ratio <= 0):
        raise ValueError("mixture_ratio must be None or finite and positive")
    for name in ("eta_c_star", "eta_CF") + (("divergence_factor",) if req.divergence_factor is not None else ()):
        value = getattr(req, name)
        if isinstance(value, bool) or not isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"{name} must lie within (0, 1]")


def _optimize_OF(req, eps):
    """optimal_OF at the thrust condition; its sea-level warning is replaced
    by this module's own separation check at the final design point."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Sea-level optimization", category=ThermoWarning)
        return optimal_OF(req.propellants, req.p_c_Pa, eps, req.thrust_condition,
                          p_amb_Pa=req.p_amb_Pa)


def _select_OF_and_eps(req):
    """Return (OF, OF_source, eps, eps_source).

    With an optimized O/F and a matched nozzle, O/F and eps are coupled:
    the best O/F depends on eps and the matched eps depends on O/F. Alternate
    the two one-dimensional problems until O/F changes by < 1e-4. At the fixed
    point, O/F maximizes Isp at the matched eps it produces.
    """
    fixed_eps = not isinstance(req.nozzle, str)
    if req.mixture_ratio is not None:
        of = float(req.mixture_ratio)
        eps = float(req.nozzle) if fixed_eps else matched_eps(req.propellants, req.p_c_Pa, of, req.p_amb_Pa)
        return of, "user", eps, "user" if fixed_eps else f"matched p_e = {req.p_amb_Pa:g} Pa (CEA)"
    if fixed_eps:
        eps = float(req.nozzle)
        return _optimize_OF(req, eps), f"max {req.thrust_condition} Isp at fixed eps", eps, "user"
    of = _optimize_OF(req, CEA_DEFAULT_EPS)
    for _ in range(_MAX_COUPLED_ITERATIONS):
        eps = matched_eps(req.propellants, req.p_c_Pa, of, req.p_amb_Pa)
        new = _optimize_OF(req, eps)
        if abs(new - of) < _OF_TOLERANCE:
            of = new
            eps = matched_eps(req.propellants, req.p_c_Pa, of, req.p_amb_Pa)
            return (of, f"max {req.thrust_condition} Isp, coupled with matched eps",
                    eps, f"matched p_e = {req.p_amb_Pa:g} Pa (CEA)")
        of = new
    raise RuntimeError("Coupled O/F and matched-eps iteration did not converge")


def _divergence(req, eps):
    """Return (lambda, source). Fig.4-12 applies only to chart-angle bells."""
    if req.divergence_factor is not None:
        return req.divergence_factor, "user-specified"
    if req.theta_n_deg is not None:
        raise ValueError("User wall angles are not Rao chart contours; supply divergence_factor")
    return (bell_divergence_factor(eps, req.bell_fraction),
            "NASA SP-125 Fig.4-12 digitization, +/-0.0005")


def design_engine(req):
    """Design an engine from EngineRequirements; return EngineDesign.

    1. O/F and eps (user, optimized, or coupled with matched expansion);
       optimization uses ideal shifting Isp. Constant efficiencies do not
       move the optimum, but lambda's eps dependence is not included there.
    2. CEA expansion, shifting and frozen: CF_vac = Isp_vac*g0/c*,
       CF_amb = CF_vac - eps*p_amb/p_c (theory section 2).
    3. Delivered = apply_losses(ideal shifting) (theory section 4).
    4. A_t = F/(CF_del*p_c) at the thrust condition; mdot = p_c*A_t/(eta_c*·c*).
       Then sizing.size_engine for flows, chamber and contour (section 3).
    Separation is screened with the Summerfield criterion at p_amb; it warns
    only when the thrust requirement itself is at that ambient pressure.
    """
    _validate(req)
    of, of_source, eps, eps_source = _select_OF_and_eps(req)
    divergence, divergence_source = _divergence(req, eps)
    args = (req.propellants, req.p_c_Pa, of, eps)
    shifting = expansion_performance(*args)
    frozen = expansion_performance(*args, frozen=True)
    ideal, ideal_frozen = (ideal_performance(e.CF_vac, e.c_star_m_s, e.p_e_Pa, eps,
                                             req.p_c_Pa, req.p_amb_Pa)
                           for e in (shifting, frozen))
    delivered = apply_losses(ideal, eps, req.p_c_Pa, req.p_amb_Pa, divergence,
                             req.eta_c_star, req.eta_CF)
    cf_design = delivered.CF_vac if req.thrust_condition == "vacuum" else delivered.CF_amb
    if cf_design <= 0:
        raise ValueError("Nonpositive thrust coefficient at the thrust condition; reduce eps")
    with warnings.catch_warnings():
        if req.thrust_condition == "vacuum":
            warnings.simplefilter("ignore", NozzleWarning)
        separated = separation_check(ideal.p_e_Pa, req.p_amb_Pa)
    chamber = chamber_props(req.propellants, req.p_c_Pa, of, eps)
    if req.L_star_m is None:
        l_star = CHARACTERISTIC_LENGTH_M[req.propellants]
        l_star_source = CHARACTERISTIC_LENGTH_SOURCE[req.propellants]
    else:
        l_star, l_star_source = req.L_star_m, "User-specified"
    geometry = size_engine(req.thrust_N, cf_design, req.p_c_Pa, delivered.c_star_m_s, of, eps,
                           l_star, contraction_ratio=req.contraction_ratio,
                           bell_fraction=req.bell_fraction,
                           convergence_half_angle_deg=req.convergence_half_angle_deg,
                           theta_n_deg=req.theta_n_deg, theta_e_deg=req.theta_e_deg)
    return EngineDesign(req, of, of_source, eps, eps_source, chamber, ideal, ideal_frozen,
                        delivered, divergence, divergence_source, separated,
                        l_star_source, geometry)
