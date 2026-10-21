"""Limited standard-atmosphere pressure model for the Stage 2 altitude plot."""
from math import exp, isfinite

from .constants import (G0_M_S2, STANDARD_ATMOSPHERE_PA, ATMOSPHERE_R_J_KG_K,
                        ATMOSPHERE_T0_K, ATMOSPHERE_LAYERS)


def ambient_pressure(geopotential_height_m):
    """Return ambient pressure [Pa] at geopotential height 0–51 km.

    Hydrostatic ideal gas: dp/p=-g0*dH/(R*T). For constant lapse rate L,
    p/p_b=(T/T_b)^(-g0/(R*L)); for L=0, exp(-g0*(H-H_b)/(R*T_b)).
    Uses U.S. Standard Atmosphere 1976 lower layers, dry air and prescribed
    temperatures; not weather or a geometric-altitude model. Reference:
    NASA/TM-2005-213659, Appendix A, eqs. A.2–A.5 and Table A1.
    Layer base pressures are propagated to preserve exact continuity.
    """
    if not isfinite(geopotential_height_m) or not 0 <= geopotential_height_m <= ATMOSPHERE_LAYERS[-1][0]:
        raise ValueError("geopotential height must be finite and within 0–51000 m")
    pressure, temperature, base = STANDARD_ATMOSPHERE_PA, ATMOSPHERE_T0_K, 0.0
    for ceiling, lapse in ATMOSPHERE_LAYERS:
        step = min(geopotential_height_m, ceiling) - base
        if lapse == 0:
            pressure *= exp(-G0_M_S2 * step / (ATMOSPHERE_R_J_KG_K * temperature))
        else:
            next_temperature = temperature + lapse * step
            pressure *= (next_temperature / temperature) ** (-G0_M_S2 / (ATMOSPHERE_R_J_KG_K * lapse))
            temperature = next_temperature
        if geopotential_height_m <= ceiling:
            return pressure
        base = ceiling
