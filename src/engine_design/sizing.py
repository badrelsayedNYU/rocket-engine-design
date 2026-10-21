"""Axisymmetric preliminary dimensions and a Rao parabolic approximation.

References: theory section 3; NASA SP-125, section 4.3, Figs. 4-13/4-14.
Coordinates are inner gas-wall surfaces, not wall thickness or hardware.
"""
from dataclasses import dataclass
from math import cos, isfinite, pi, radians, sin, sqrt, tan

import numpy as np
from scipy.integrate import quad

from .constants import (RAO_AREA_RATIOS, RAO_ANGLE_TABLE,
    UPSTREAM_THROAT_RADIUS_RATIO, DOWNSTREAM_THROAT_RADIUS_RATIO,
    CHAMBER_FILLET_RADIUS_RATIO, DEFAULT_CONTRACTION_RATIO,
    DEFAULT_CONVERGENCE_ANGLE_DEG, DEFAULT_BELL_FRACTION, REFERENCE_CONE_ANGLE_DEG)


@dataclass(frozen=True)
class ContourSegment:
    """Sampled inner wall with analytic dr/dx and integrated volume of revolution."""
    name: str
    x_m: np.ndarray
    r_m: np.ndarray
    slope: np.ndarray
    volume_m3: float


@dataclass(frozen=True)
class EngineGeometry:
    """SI sizing result; chamber volume includes all geometry upstream of throat."""
    A_t_m2: float
    A_e_m2: float
    A_c_m2: float
    D_t_m: float
    D_e_m: float
    D_c_m: float
    mdot_kg_s: float
    mdot_ox_kg_s: float
    mdot_fuel_kg_s: float
    V_c_m3: float
    L_cylinder_m: float
    L_converging_m: float
    L_chamber_m: float
    L_nozzle_m: float
    L_cone_m: float
    L_star_m: float
    theta_n_deg: float
    theta_e_deg: float
    angle_source: str
    segments: tuple[ContourSegment, ...]

    @property
    def x_m(self):
        """Axial wall samples, duplicate joint endpoints removed; throat x=0."""
        return np.concatenate([s.x_m if i == 0 else s.x_m[1:]
                               for i, s in enumerate(self.segments)])

    @property
    def r_m(self):
        """Wall radius at each x_m sample."""
        return np.concatenate([s.r_m if i == 0 else s.r_m[1:]
                               for i, s in enumerate(self.segments)])


def _positive(name, value):
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def rao_angles(eps, bell_fraction=DEFAULT_BELL_FRACTION):
    """Return (theta_n, theta_e) in degrees by piecewise bilinear interpolation.

    Table is this project's manual digitization of NASA SP-125 Fig.4-14;
    precision ~0.5 degree, 5<=eps<=50 and 0.6<=bell_fraction<=1.0.
    Approximates Rao's near-optimum parabolic contours, not a new optimization.
    No extrapolation beyond the chart. See constants.py and docs/references/.
    """
    if not isfinite(eps) or not RAO_AREA_RATIOS[0] <= eps <= RAO_AREA_RATIOS[-1]:
        raise ValueError("Rao chart supports 5 <= eps <= 50; supply explicit wall angles outside it")
    fractions = sorted(RAO_ANGLE_TABLE)
    if not isfinite(bell_fraction) or not fractions[0] <= bell_fraction <= fractions[-1]:
        raise ValueError("Rao chart supports 0.6 <= bell_fraction <= 1.0")
    return tuple(float(np.interp(bell_fraction, fractions,
                 [np.interp(eps, RAO_AREA_RATIOS, RAO_ANGLE_TABLE[f][j]) for f in fractions]))
                 for j in (0, 1))


def _segment(name, curve, count):
    """Sample a parametric wall and integrate V=integral(pi*r²*dx).

    curve(t) returns x,r,dx/dt,dr/dt on [0,1]. Axisymmetric geometry;
    no flow assumptions. Volume is evaluated by adaptive quadrature,
    independently of plotting resolution (theory section 3, Vc=L*At).
    """
    values = np.array([curve(t) for t in np.linspace(0, 1, count)])
    x, r, dx, dr = values.T
    if np.any(dx <= 0) or np.any(r <= 0):
        raise ValueError(f"Invalid nonmonotonic or nonpositive {name} geometry")
    volume = quad(lambda t: pi * curve(t)[1]**2 * curve(t)[2], 0, 1,
                  epsabs=1e-13, epsrel=1e-10)[0]
    slope = dr / dx
    for a in (x, r, slope):
        a.flags.writeable = False
    return ContourSegment(name, x, r, slope, volume)


def _line(name, start, end, count):
    """Straight wall interpolation; dx and dr are constant geometric derivatives."""
    x0, r0 = start
    dx, dr = end[0] - x0, end[1] - r0
    return _segment(name, lambda t: (x0 + t*dx, r0 + t*dr, dx, dr), count)


def size_engine(thrust_N, CF, p_c_Pa, c_star_m_s, OF, eps, L_star_m, *,
                contraction_ratio=DEFAULT_CONTRACTION_RATIO,
                bell_fraction=DEFAULT_BELL_FRACTION,
                convergence_half_angle_deg=DEFAULT_CONVERGENCE_ANGLE_DEG,
                theta_n_deg=None, theta_e_deg=None, points_per_segment=100):
    """Size throat, flows, chamber and wall contour using theory section 3.

    At=F/(CF*pc), Ae=eps*At, mdot=pc*At/c*, Vc=L*At, Ac=CR*At.
    mdot_ox=mdot*OF/(1+OF), mdot_fuel=mdot/(1+OF). Axisymmetric geometry.
    Vc includes cylinder, chamber fillet, cone and upstream throat arc.
    Cylinder length=(Vc-Vconverging)/Ac; impossible volumes are rejected.
    Throat arcs have radii 1.5Rt and 0.382Rt (NASA SP-125 Fig.4-13).
    The bell is a quadratic Bezier (rotated parabola) whose control point is
    the intersection of endpoint tangents. Ln=bell_fraction*(Re-Rt)/tan(15°)
    follows the SPEC's straight-cone convention, omitting the small rounded
    cone correction in SP-125 eq.4-7. A 0.5Rt chamber fillet is an explicit
    geometric assumption to ensure slope continuity at the cylinder.
    """
    for name, value in dict(thrust_N=thrust_N, CF=CF, p_c_Pa=p_c_Pa,
        c_star_m_s=c_star_m_s, OF=OF, eps=eps, L_star_m=L_star_m,
        contraction_ratio=contraction_ratio, bell_fraction=bell_fraction).items():
        _positive(name, value)
    if eps <= 1 or contraction_ratio <= 1:
        raise ValueError("eps and contraction_ratio must exceed 1 for a chamber and bell")
    if not isfinite(convergence_half_angle_deg) or not 0 < convergence_half_angle_deg < 90:
        raise ValueError("convergence half angle must lie between 0 and 90 degrees")
    if not 0 < bell_fraction <= 1:
        raise ValueError("bell_fraction must be within (0,1]")
    if isinstance(points_per_segment, bool) or not isinstance(points_per_segment, int) or points_per_segment < 3:
        raise ValueError("points_per_segment must be an integer >= 3")
    if (theta_n_deg is None) != (theta_e_deg is None):
        raise ValueError("Supply both wall angles, or neither")
    if theta_n_deg is None:
        theta_n_deg, theta_e_deg = rao_angles(eps, bell_fraction)
        angle_source = "NASA SP-125 Fig.4-14; manual digitization +/-0.5 deg"
    else:
        angle_source = "User-specified wall angles; not validated by Rao chart"
    if not all(isfinite(a) for a in (theta_n_deg, theta_e_deg)) or not 0 < theta_e_deg < theta_n_deg < 90:
        raise ValueError("Wall angles must satisfy 0 < theta_e < theta_n < 90 degrees")
    at = thrust_N / (CF * p_c_Pa)
    ae, ac = eps * at, contraction_ratio * at
    rt, re, rc = sqrt(at/pi), sqrt(ae/pi), sqrt(ac/pi)
    mdot = p_c_Pa * at / c_star_m_s
    volume = L_star_m * at
    alpha, tn, te = map(radians, (convergence_half_angle_deg, theta_n_deg, theta_e_deg))
    ru = UPSTREAM_THROAT_RADIUS_RATIO * rt
    rd = DOWNSTREAM_THROAT_RADIUS_RATIO * rt
    rf = CHAMBER_FILLET_RADIUS_RATIO * rt
    # A joins cone to upstream throat arc; B joins chamber fillet to cone.
    xa, ra = -ru*sin(alpha), rt + ru*(1-cos(alpha))
    rb = rc - rf*(1-cos(alpha))
    if rb <= ra:
        raise ValueError("Contraction ratio too small for selected converging angle and fillets")
    xb = xa - (rb-ra)/tan(alpha)
    x0 = xb - rf*sin(alpha)
    n = points_per_segment
    fillet = _segment("chamber fillet", lambda t: (
        x0+rf*sin(alpha*t), rc-rf*(1-cos(alpha*t)),
        rf*alpha*cos(alpha*t), -rf*alpha*sin(alpha*t)), n)
    cone = _line("converging cone", (xb,rb), (xa,ra), n)
    upstream = _segment("upstream throat arc", lambda t: (
        ru*sin(alpha*(t-1)), rt+ru*(1-cos(alpha*(t-1))),
        ru*alpha*cos(alpha*(t-1)), ru*alpha*sin(alpha*(t-1))), n)
    v_converging = sum(s.volume_m3 for s in (fillet, cone, upstream))
    l_cylinder = (volume-v_converging)/ac
    if l_cylinder <= 0:
        raise ValueError("L* volume is too small for converging geometry; increase L* or change geometry")
    cylinder = _line("cylinder", (x0-l_cylinder,rc), (x0,rc), n)
    downstream = _segment("downstream throat arc", lambda t: (
        rd*sin(tn*t), rt+rd*(1-cos(tn*t)),
        rd*tn*cos(tn*t), rd*tn*sin(tn*t)), n)
    xn, rn = downstream.x_m[-1], downstream.r_m[-1]
    l_cone = (re-rt)/tan(radians(REFERENCE_CONE_ANGLE_DEG))
    length = bell_fraction*l_cone
    mn, me = tan(tn), tan(te)
    qx = (re-rn+mn*xn-me*length)/(mn-me)
    qr = rn+mn*(qx-xn)
    if not xn < qx < length or not rn < qr < re:
        raise ValueError("Wall angles and bell length cannot form a monotonic parabolic contour")
    def bezier(t):
        """B=(1-t)^2*N+2t(1-t)*Q+t^2*E; exact tangent derivatives."""
        u = 1-t
        return (u*u*xn+2*t*u*qx+t*t*length, u*u*rn+2*t*u*qr+t*t*re,
                2*(u*(qx-xn)+t*(length-qx)), 2*(u*(qr-rn)+t*(re-qr)))
    bell = _segment("parabolic bell", bezier, n)
    return EngineGeometry(at, ae, ac, 2*rt, 2*re, 2*rc, mdot,
        mdot*OF/(1+OF), mdot/(1+OF), volume, l_cylinder, -x0,
        l_cylinder-x0, length, l_cone, L_star_m, theta_n_deg, theta_e_deg,
        angle_source, (cylinder, fillet, cone, upstream, downstream, bell))
