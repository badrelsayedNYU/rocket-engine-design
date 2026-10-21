"""Steady 1D isentropic, calorically perfect, choked-nozzle relations.

No variable-composition expansion, shocks, losses, or separated-flow solution.
Reference: docs/spec/01_THEORY.md section 2; Sutton & Biblarz, Ch. 3;
https://www.grc.nasa.gov/WWW/k-12/airplane/rktthsum.html
"""
from math import exp, expm1, isfinite, log, log1p, sqrt
import warnings

from scipy.optimize import brentq

from .constants import G0_M_S2, SEPARATION_PRESSURE_RATIO


class NozzleWarning(UserWarning):
    """The attached-flow prediction may not describe the actual flow."""


def _gamma(gamma):
    if not isfinite(gamma) or gamma <= 1:
        raise ValueError("gamma must be finite and > 1")


def _pressure(name, value, *, allow_zero=False):
    if not isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")


def _log_area(M, gamma):
    # log1p preserves accuracy near the sonic throat and gamma=1.
    return -log(M) + (gamma + 1) / (2 * (gamma - 1)) * log1p(
        (gamma - 1) * (M * M - 1) / (gamma + 1))


def area_ratio(M, gamma):
    """A/A* = (1/M)[2/(gamma+1)*(1+(gamma-1)*M²/2)]^((gamma+1)/(2(gamma-1))).

    Steady 1D isentropic perfect gas, constant gamma, M>0. Works on either
    branch; A*=sonic area. Reference: theory section 2; Sutton Ch. 3.
    """
    _gamma(gamma)
    if not isfinite(M) or M <= 0:
        raise ValueError("M must be finite and positive")
    return exp(_log_area(M, gamma))


def mach_from_area_ratio(eps, gamma):
    """Invert the area–Mach relation on the supersonic branch using Brent's method.

    eps=A/A* >= 1; eps=1 returns M=1 exactly. Constant-gamma, choked,
    isentropic flow; theory section 2 and Sutton Ch. 3. The bracket expands
    until it encloses the root; the residual is logarithmic for stability.
    """
    _gamma(gamma)
    if not isfinite(eps) or eps < 1:
        raise ValueError("eps must be finite and >= 1")
    if eps == 1:
        return 1.0
    target = log(eps)
    upper = 2.0
    # Numerical safeguard, not a physical upper Mach limit.
    for _ in range(128):
        if _log_area(upper, gamma) >= target:
            return float(brentq(lambda m: _log_area(m, gamma) - target,
                                1.0, upper, xtol=1e-12, rtol=1e-12))
        upper *= 2
    raise ValueError("Could not bracket Mach number for these extreme inputs")


def pressure_ratio(M, gamma):
    """p/p_c = [1+(gamma-1)*M²/2]^(-gamma/(gamma-1)).

    Chamber pressure is stagnation pressure; steady isentropic perfect gas
    with constant gamma. M>=0. Reference: theory section 2; Sutton Ch. 3.
    """
    _gamma(gamma)
    if not isfinite(M) or M < 0:
        raise ValueError("M must be finite and nonnegative")
    return exp(-gamma / (gamma - 1) * log1p((gamma - 1) * M * M / 2))


def separation_check(p_e, p_a, k_sep=SEPARATION_PRESSURE_RATIO):
    """Flag likely separation when p_e < k_sep*p_a and emit NozzleWarning.

    Summerfield screening rule, default k_sep=0.4. Exit pressure approximates
    minimum wall pressure; contour and boundary-layer effects are omitted.
    False is not proof of attachment. Vacuum always returns False.
    Reference: theory section 2; Schmucker, NASA TM-77396, eq. (7), p. 37
    (PDF p. 43). This is a heuristic, not a separated-flow calculation.
    """
    _pressure("p_e", p_e, allow_zero=True)
    _pressure("p_a", p_a, allow_zero=True)
    if not isfinite(k_sep) or not 0 < k_sep < 1:
        raise ValueError("k_sep must lie strictly between 0 and 1")
    separated = p_e < k_sep * p_a
    if separated:
        warnings.warn("Likely nozzle separation: p_e < k_sep*p_a. Attached-flow "
                      "performance is an extrapolation; separated flow is not modeled.",
                      NozzleWarning, stacklevel=2)
    return separated


def thrust_coefficient(gamma, eps, p_c, p_a):
    """CF = sqrt[2γ²/(γ−1)*(2/(γ+1))^((γ+1)/(γ−1))
    * (1−(p_e/p_c)^((γ−1)/γ))] + eps*(p_e−p_a)/p_c.

    Constant-gamma, steady 1D isentropic choked expansion on the supersonic
    branch. Reference: theory section 2; Sutton Ch. 3. Warns for likely
    separation and returns the unmodified attached-flow formula, even if
    negative. Does not predict shocks or actual separated performance.
    Requires p_a < p_c; this alone does not guarantee the assumed flow.
    """
    _pressure("p_c", p_c)
    _pressure("p_a", p_a, allow_zero=True)
    if p_a >= p_c:
        raise ValueError("p_a must be below chamber stagnation pressure")
    mach = mach_from_area_ratio(eps, gamma)
    pe_pc = pressure_ratio(mach, gamma)
    # Equivalent to 1-(pe/pc)^((gamma-1)/gamma), without subtraction loss.
    velocity_factor = -expm1(-log1p((gamma - 1) * mach * mach / 2))
    choke_factor = exp(-(gamma + 1) / (gamma - 1) * log1p((gamma - 1) / 2))
    momentum = sqrt(2 * gamma * gamma / (gamma - 1) * choke_factor * velocity_factor)
    separation_check(pe_pc * p_c, p_a)
    return momentum + eps * (pe_pc - p_a / p_c)


def isp(c_star, CF):
    """Isp [s] = c* [m/s] * CF / g0; theory section 2, Sutton Ch. 3.

    Uses the same chamber/nozzle model as the supplied c* and CF. Positive
    c* required; finite signed CF is preserved to expose invalid extrapolation
    rather than silently clipping an attached-flow prediction to zero.
    """
    _pressure("c_star", c_star)
    if not isfinite(CF):
        raise ValueError("CF must be finite")
    return c_star * CF / G0_M_S2


def optimal_eps(gamma, p_c, p_a):
    """Return the choked supersonic area ratio for matched expansion, p_e=p_a.

    M_e² = 2/(gamma−1)*[(p_c/p_a)^((gamma−1)/gamma)−1], followed by
    the area–Mach relation. Constant-gamma isentropic flow; theory sections
    2–3 and Sutton Ch. 3. Vacuum has no finite optimum: raises ValueError.
    Ambient pressure above the sonic pressure cannot be matched on this
    branch and is rejected. Exactly sonic expansion returns eps=1.
    """
    _gamma(gamma)
    _pressure("p_c", p_c)
    _pressure("p_a", p_a, allow_zero=True)
    if p_a == 0:
        raise ValueError("Vacuum has no finite optimum area ratio; supply a fixed eps")
    critical = pressure_ratio(1.0, gamma)
    if p_a / p_c > critical:
        raise ValueError("Ambient pressure exceeds the sonic pressure; no matched supersonic solution")
    mach_squared = 2 / (gamma - 1) * expm1((gamma - 1) / gamma * log(p_c / p_a))
    return area_ratio(sqrt(max(1.0, mach_squared)), gamma)
