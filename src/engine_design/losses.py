"""Stage 4: divergence, combustion and nozzle efficiencies; delivered values.

Delivered performance corrects ideal CEA results with three factors:
lambda (nonaxial exit momentum), eta_c* (incomplete combustion) and eta_CF
(other nozzle losses). References: theory section 4; NASA SP-125 eqs. 1-39
to 1-43 and 4-8, Fig. 4-12; Sutton & Biblarz Ch. 3.
"""
from dataclasses import dataclass
from math import cos, isfinite, radians

import numpy as np

from .constants import (BELL_DIVERGENCE_AREA_RATIOS, BELL_DIVERGENCE_FRACTIONS,
                        BELL_DIVERGENCE_PERCENT, G0_M_S2)


@dataclass(frozen=True)
class Performance:
    """c*, vacuum and ambient CF, and exit pressure for one model, SI units."""
    c_star_m_s: float
    CF_vac: float
    CF_amb: float
    p_e_Pa: float

    @property
    def Isp_vac_s(self):
        """Isp = c*·CF/g0 in vacuum."""
        return self.c_star_m_s * self.CF_vac / G0_M_S2

    @property
    def Isp_amb_s(self):
        """Isp = c*·CF/g0 at ambient pressure (attached flow, signed)."""
        return self.c_star_m_s * self.CF_amb / G0_M_S2


def _fraction(name, value):
    if not isfinite(value) or not 0 < value <= 1:
        raise ValueError(f"{name} must lie within (0, 1]")


def conical_divergence_factor(half_angle_deg):
    """lambda = (1 + cos(alpha))/2 for a conical nozzle of half angle alpha.

    Ratio of axial exit momentum from a point source spreading over a
    spherical exit cap to that of uniform axial flow. Theory section 4;
    NASA SP-125 eq. 4-8; Sutton & Biblarz Ch. 3. lambda(15°) = 0.983.
    """
    if not isfinite(half_angle_deg) or not 0 <= half_angle_deg < 90:
        raise ValueError("Cone half angle must lie within [0, 90) degrees")
    return (1 + cos(radians(half_angle_deg))) / 2


def bell_divergence_factor(eps, bell_fraction):
    """Return lambda for a Rao bell by bilinear interpolation of SP-125 Fig.4-12.

    The chart gives thrust efficiency of optimum bell contours as a percent
    of ideal CF, versus length relative to a 15° cone. Valid only for
    10 <= eps <= 40 and 0.6 <= bell_fraction <= 1.0; no extrapolation.
    It describes a near-optimum contour; other wall angles need their own value.
    """
    ratios, fractions = BELL_DIVERGENCE_AREA_RATIOS, BELL_DIVERGENCE_FRACTIONS
    if not isfinite(eps) or not ratios[0] <= eps <= ratios[-1]:
        raise ValueError("SP-125 Fig.4-12 covers 10 <= eps <= 40; supply divergence_factor")
    if not isfinite(bell_fraction) or not fractions[0] <= bell_fraction <= fractions[-1]:
        raise ValueError("SP-125 Fig.4-12 covers 0.6 <= bell_fraction <= 1.0")
    by_ratio = [np.interp(bell_fraction, fractions, row) for row in BELL_DIVERGENCE_PERCENT]
    return float(np.interp(eps, ratios, by_ratio)) / 100


def ideal_performance(CF_vac, c_star_m_s, p_e_Pa, eps, p_c_Pa, p_amb_Pa):
    """Ideal Performance; CF_amb = CF_vac - eps*p_a/p_c (theory section 2)."""
    return Performance(c_star_m_s, CF_vac, CF_vac - eps * p_amb_Pa / p_c_Pa, p_e_Pa)


def apply_losses(ideal, eps, p_c_Pa, p_amb_Pa, divergence_factor, eta_c_star, eta_CF):
    """Return delivered Performance from an ideal one.

    Split ideal vacuum CF into momentum and exit-pressure parts,
    CF_vac = CF_mom + eps*p_e/p_c. Apply lambda to the momentum part only:
    CF_vac,del = eta_CF*(lambda*CF_mom + eps*p_e/p_c) (theory section 4).
    The ambient term eps*p_a/p_c is an exact pressure force on the exit area,
    so it is subtracted after the efficiencies: CF_amb,del = CF_vac,del - eps*p_a/p_c.
    c*_del = eta_c*·c*. In vacuum this gives Isp_del = eta_c*·eta_CF·Isp,
    where Isp includes divergence, as in theory section 4.
    """
    for name, value in (("divergence_factor", divergence_factor),
                        ("eta_c_star", eta_c_star), ("eta_CF", eta_CF)):
        _fraction(name, value)
    pressure = eps * ideal.p_e_Pa / p_c_Pa
    momentum = ideal.CF_vac - pressure
    if momentum <= 0:
        raise ValueError("Ideal CF is inconsistent with its exit pressure")
    cf_vac = eta_CF * (divergence_factor * momentum + pressure)
    return Performance(eta_c_star * ideal.c_star_m_s, cf_vac,
                       cf_vac - eps * p_amb_Pa / p_c_Pa, ideal.p_e_Pa)
