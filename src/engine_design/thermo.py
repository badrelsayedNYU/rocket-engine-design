"""CEA thermochemistry in SI; no nozzle geometry or loss models.

Reference: NASA RP-1311 and RocketCEA 1.2.3 functions documentation:
https://rocketcea.readthedocs.io/en/latest/functions.html
CEA's shared Fortran state is serialized here. Do not mix direct RocketCEA
calls from other threads with this interface.
"""

from dataclasses import dataclass
from contextlib import contextmanager
from functools import lru_cache
from math import exp, isfinite, log1p, sqrt
from tempfile import TemporaryDirectory
from threading import RLock
import warnings

import numpy as np
from scipy.optimize import minimize_scalar

from .constants import (
    G0_M_S2, R_U_J_KMOL_K, STANDARD_ATMOSPHERE_PA, C_STAR_SANITY_M_S,
    GAMMA_SANITY, DEFAULT_OF_BOUNDS, CEA_DEFAULT_EPS,
    MILLIPOISE_TO_PA_S, MCAL_CM_S_K_TO_W_M_K, CAL_G_K_TO_J_KG_K,
)
from .propellants import PROPELLANTS

_LOCK = RLock()


class ThermoWarning(UserWarning):
    """A model result needs interpretation, not automatic rejection."""


@dataclass(frozen=True)
class TransportProperties:
    """Gas transport at a station, with SI units and dimensionless Prandtl."""
    cp_J_kg_K: float
    mu_Pa_s: float
    conductivity_W_m_K: float
    Pr: float


@dataclass(frozen=True)
class ChamberProperties:
    """Infinite-area chamber and sonic throat; molecular weights in kg/kmol.

    gamma is CEA's isentropic exponent (equilibrium by default), not necessarily
    cp/cv for a reacting mixture. frozen=True freezes composition at chamber.
    """
    T_c_K: float
    M_kg_kmol: float
    gamma: float
    c_star_m_s: float
    T_t_K: float
    M_t_kg_kmol: float
    gamma_t: float
    chamber_transport: TransportProperties
    throat_transport: TransportProperties
    frozen: bool


@dataclass(frozen=True)
class OFSweep:
    """One-dimensional arrays in the same order as the supplied O/F values."""
    OF: np.ndarray
    T_c_K: np.ndarray
    M_kg_kmol: np.ndarray
    gamma: np.ndarray
    c_star_m_s: np.ndarray


@lru_cache(maxsize=1)
def _backend():
    try:
        from rocketcea import cea_obj
        from rocketcea.cea_obj_w_units import CEA_Obj
    except ImportError as exc:
        raise ImportError("Install the CEA backend: pip install -e '.[cea]'") from exc
    # Retain the object so scratch files survive until interpreter shutdown.
    scratch = TemporaryDirectory(prefix="engine-cea-")
    return cea_obj, CEA_Obj, scratch


@contextmanager
def _cea(propellants):
    if propellants not in PROPELLANTS:
        raise ValueError(f"Unknown propellants {propellants!r}; choose {tuple(PROPELLANTS)}")
    with _LOCK:
        module, constructor, scratch = _backend()
        previous = module.ROCKETCEA_DATA_DIR
        module.ROCKETCEA_DATA_DIR = scratch.name
        try:
            pair = PROPELLANTS[propellants]
            yield constructor(
                oxName=pair.cea_oxidizer, fuelName=pair.cea_fuel,
                pressure_units="Pa", temperature_units="K", cstar_units="m/s",
                specific_heat_units="cal/g-K", sonic_velocity_units="m/s",
                viscosity_units="millipoise", thermal_cond_units="mcal/cm-K-s",
            )
        finally:
            module.ROCKETCEA_DATA_DIR = previous


def _positive(name, value):
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _inputs(p_c_Pa, OF, eps):
    _positive("p_c_Pa", p_c_Pa)
    _positive("OF", OF)
    if not isfinite(eps) or eps < 1:
        raise ValueError("eps must be finite and >= 1")


def ideal_c_star(T_c_K, M_kg_kmol, gamma):
    """Return perfect-gas c* in m/s using theory section 1.

    c* = sqrt(R*T/gamma)*((gamma+1)/2)**((gamma+1)/(2*(gamma-1))),
    R = R_u/M. Assumes steady choked isentropic flow, constant composition
    and constant gamma. Reference: Sutton & Biblarz, 9th ed., Ch. 3.
    This is an approximation when evaluated using equilibrium CEA properties.
    """
    _positive("T_c_K", T_c_K)
    _positive("M_kg_kmol", M_kg_kmol)
    if not isfinite(gamma) or gamma <= 1:
        raise ValueError("gamma must be finite and > 1")
    exponent = (gamma + 1) / (2 * (gamma - 1)) * log1p((gamma - 1) / 2)
    return sqrt(R_U_J_KMOL_K / M_kg_kmol * T_c_K / gamma) * exp(exponent)


def _transport(values):
    cp, mu, conductivity, pr = map(float, values)
    result = TransportProperties(cp * CAL_G_K_TO_J_KG_K, mu * MILLIPOISE_TO_PA_S,
                                 conductivity * MCAL_CM_S_K_TO_W_M_K, pr)
    for name, value in vars(result).items():
        if not isfinite(value) or value <= 0:
            raise RuntimeError(f"CEA returned invalid transport {name}={value}")
    return result


def _frozen_transport(cea, args):
    """Read transport at the chamber-frozen throat, not equilibrium throat.

    RocketCEA 1.2.3 get_Throat_Transport(frozen=1) selects frozen transport
    coefficients but calls setupCards without frozen=1. Explicitly set the
    flow mode, then immediately copy the Fortran outputs under our lock.
    Reference: rocketcea/cea_obj.py, setupCards and get_Throat_Transport.
    """
    from rocketcea import cea_obj
    raw = cea.cea_obj
    raw.setupCards(Pc=cea.Pc_U.uval_to_dval(args["Pc"]), MR=args["MR"],
                   eps=args["eps"], frozen=1, frozenAtThroat=0, show_transport=1)
    tr = cea_obj.py_cea.trpts
    stations = []
    for index in (raw.i_chm, raw.i_thrt):
        stations.append(_transport((cea.specific_heat_U.dval_to_uval(tr.cpfro[index]),
                                   tr.vis[index], tr.confro[index], tr.prfro[index])))
    gamma = float(cea_obj.py_cea.prtout.gammas[raw.i_chm])
    return stations[0], stations[1], gamma


def chamber_props(propellants, p_c_Pa, OF, eps=None, frozen=False):
    """Return chamber/throat thermochemistry and transport, entirely in SI.

    CEA solves adiabatic constant-pressure equilibrium combustion by Gibbs
    energy minimization with elemental conservation, then isentropic choked
    expansion (NASA RP-1311; theory section 1). c* = p_c*A_t/mdot.
    frozen=True holds chamber composition fixed during expansion; chamber
    combustion itself remains equilibrium. eps=None uses a computational
    exit ratio of 40; chamber/throat quantities do not depend on exit ratio.
    """
    eps = CEA_DEFAULT_EPS if eps is None else eps
    _inputs(p_c_Pa, OF, eps)
    if not isinstance(frozen, bool):
        raise ValueError("frozen must be a bool")
    args = dict(Pc=p_c_Pa, MR=OF, eps=eps)
    with _cea(propellants) as cea:
        _, c_star, temperature = cea.get_IvacCstrTc(**args, frozen=int(frozen))
        mw, gamma = cea.get_Chamber_MolWt_gamma(**args)
        mw_t, gamma_t = cea.get_Throat_MolWt_gamma(**args, frozen=int(frozen))
        t_t = cea.get_Temperatures(**args, frozen=int(frozen))[1]
        if frozen:
            chamber_tr, throat_tr, gamma = _frozen_transport(cea, args)
        else:
            chamber_tr = _transport(cea.get_Chamber_Transport(**args))
            throat_tr = _transport(cea.get_Throat_Transport(**args))
    scalars = tuple(map(float, (temperature, mw, gamma, c_star, t_t, mw_t, gamma_t)))
    if not all(isfinite(v) and v > 0 for v in scalars) or min(gamma, gamma_t) <= 1:
        raise RuntimeError("CEA returned nonphysical chamber/throat properties")
    result = ChamberProperties(*scalars, chamber_tr, throat_tr, frozen)
    low, high = C_STAR_SANITY_M_S[propellants]
    if not low <= c_star <= high or not GAMMA_SANITY[0] <= gamma <= GAMMA_SANITY[1]:
        warnings.warn(
            f"{propellants} at O/F={OF:g}: c*={c_star:.1f} m/s, gamma={gamma:.4f} "
            "outside the spec's nominal equilibrium bands; these are not universal "
            "limits, especially off-design or with frozen chemistry.",
            ThermoWarning, stacklevel=2,
        )
    return result


@dataclass(frozen=True)
class ExpansionPerformance:
    """CEA expansion to a fixed area ratio; ideal, no losses, SI units."""
    Isp_vac_s: float
    c_star_m_s: float
    p_e_Pa: float
    CF_vac: float


def expansion_performance(propellants, p_c_Pa, OF, eps, frozen=False):
    """Return CEA vacuum Isp, c*, exit pressure and vacuum CF at fixed eps.

    Isp_vac = (v_e + p_e*A_e/mdot)/g0 for equilibrium (or chamber-frozen)
    isentropic expansion to eps; p_e = p_c/(p_c/p_e). CF_vac = Isp_vac*g0/c*,
    from Isp = c*·CF/g0 (theory section 2). Ambient CF is then
    CF_vac - eps*p_a/p_c. Attached flow, no losses. NASA RP-1311.
    """
    _inputs(p_c_Pa, OF, eps)
    if eps == 1:
        raise ValueError("eps must exceed 1 for a supersonic exit")
    args = dict(Pc=p_c_Pa, MR=OF, eps=eps, frozen=int(bool(frozen)))
    with _cea(propellants) as cea:
        isp_vac, c_star, _ = cea.get_IvacCstrTc(**args)
        pc_over_pe = cea.get_PcOvPe(**args)
    isp_vac, c_star, pc_over_pe = map(float, (isp_vac, c_star, pc_over_pe))
    if not all(isfinite(v) and v > 0 for v in (isp_vac, c_star)) or not pc_over_pe > 1:
        raise RuntimeError("CEA returned nonphysical expansion properties")
    return ExpansionPerformance(isp_vac, c_star, p_c_Pa / pc_over_pe,
                                isp_vac * G0_M_S2 / c_star)


def matched_eps(propellants, p_c_Pa, OF, p_amb_Pa, frozen=False):
    """Return the area ratio at which CEA expansion gives p_e = p_amb.

    Matched (optimum) expansion, theory section 3: the pressure-thrust term
    vanishes and attached-flow CF is maximal at this ambient pressure.
    Uses CEA's variable-property expansion, not the constant-gamma formula
    in nozzle.optimal_eps, which remains an independent cross-check.
    """
    _inputs(p_c_Pa, OF, 1.0)
    _positive("p_amb_Pa", p_amb_Pa)
    if p_amb_Pa >= p_c_Pa:
        raise ValueError("p_amb_Pa must be below chamber pressure")
    with _cea(propellants) as cea:
        eps = float(cea.get_eps_at_PcOvPe(Pc=p_c_Pa, MR=OF, PcOvPe=p_c_Pa / p_amb_Pa,
                                          frozen=int(bool(frozen))))
    if not isfinite(eps) or eps <= 1:
        raise ValueError("No supersonic matched expansion; ambient pressure is too high")
    return eps


def sweep_OF(propellants, p_c_Pa, OF_range):
    """Evaluate equilibrium chamber properties at explicit O/F samples.

    Applies chamber_props independently at each point; same adiabatic CEA
    assumptions and NASA RP-1311 reference. OF_range is a 1D iterable of
    samples, not a (start, stop) specification. Arrays retain input order.
    """
    values = np.asarray(list(OF_range), dtype=float)
    if values.ndim != 1 or not values.size:
        raise ValueError("OF_range must be a nonempty one-dimensional iterable")
    if not np.all(np.isfinite(values) & (values > 0)):
        raise ValueError("Every O/F sample must be finite and positive")
    rows = [chamber_props(propellants, p_c_Pa, float(of)) for of in values]
    return OFSweep(values, *(np.array([getattr(row, field) for row in rows])
                           for field in ("T_c_K", "M_kg_kmol", "gamma", "c_star_m_s")))


def _isp(cea, p_c_Pa, OF, eps, p_amb_Pa):
    """Attached-flow Isp = Isp_vac - p_a*eps*c*/(p_c*g0).

    Derived from F=mdot*v_e+(p_e-p_a)*A_e and c*=p_c*A_t/mdot;
    theory section 2, Sutton & Biblarz Ch. 3. CEA equilibrium expansion,
    fixed exit ratio, no separation correction or losses.
    """
    isp_vac, c_star, _ = cea.get_IvacCstrTc(Pc=p_c_Pa, MR=OF, eps=eps)
    result = isp_vac - p_amb_Pa * eps * c_star / (p_c_Pa * G0_M_S2)
    if not isfinite(result) or result <= 0:
        raise ValueError("Nonpositive attached-flow Isp; choose a smaller eps or ambient pressure")
    return float(result)


def optimal_OF(propellants, p_c_Pa, eps, condition, *, OF_bounds=None,
               p_amb_Pa=STANDARD_ATMOSPHERE_PA):
    """Maximize equilibrium attached-flow Isp at fixed eps; return O/F.

    Objective is Isp_vac - p_a*eps*c*/(p_c*g0), as in theory section 2.
    condition is 'vacuum' or 'sea_level'; 'sea_level' uses p_amb_Pa
    (default: standard atmosphere). Search a
    documented fuel-specific interval, or supply OF_bounds=(low, high).
    A coarse grid locates local peaks, each refined with a bounded solver;
    endpoints are also compared. This is a numerical search within bounds,
    not a guarantee outside them. No nozzle sizing or separation model.
    References: NASA RP-1311; Sutton & Biblarz Ch. 3 and 5.
    """
    _inputs(p_c_Pa, 1.0, eps)
    if condition not in ("vacuum", "sea_level"):
        raise ValueError("condition must be 'vacuum' or 'sea_level'")
    if propellants not in PROPELLANTS:
        raise ValueError(f"Unknown propellants {propellants!r}")
    bounds = DEFAULT_OF_BOUNDS[propellants] if OF_bounds is None else OF_bounds
    if len(bounds) != 2:
        raise ValueError("OF_bounds must contain two values")
    low, high = map(float, bounds)
    _positive("lower O/F bound", low)
    _positive("upper O/F bound", high)
    if low >= high:
        raise ValueError("O/F bounds must be strictly increasing")
    _positive("p_amb_Pa", p_amb_Pa)
    p_amb = 0.0 if condition == "vacuum" else p_amb_Pa
    with _cea(propellants) as cea:
        objective = lambda of: -_isp(cea, p_c_Pa, float(of), eps, p_amb)
        # Numerical resolution and tolerance, not empirical physics constants.
        grid = np.linspace(low, high, 65)
        scores = [objective(of) for of in grid]
        candidates = [(scores[i], float(grid[i])) for i in range(len(grid))]
        for i in range(1, len(grid) - 1):
            if scores[i] <= scores[i - 1] and scores[i] <= scores[i + 1]:
                fit = minimize_scalar(objective, bounds=(grid[i - 1], grid[i + 1]),
                                      method="bounded", options={"xatol": 1e-5})
                if not fit.success:
                    raise RuntimeError(f"O/F optimization failed: {fit.message}")
                candidates.append((float(fit.fun), float(fit.x)))
        _, best = min(candidates)
    if best in (low, high):
        warnings.warn("Best O/F lies at a search boundary; consider wider bounds.",
                      ThermoWarning, stacklevel=2)
    if p_amb:
        warnings.warn("Sea-level optimization assumes attached flow. Separation must "
                      "be checked in Stage 2 before using this result for design.",
                      ThermoWarning, stacklevel=2)
    return best
