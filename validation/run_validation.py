"""Stage 6: predicted vs published Isp for real engines.

Run: python validation/run_validation.py
Writes validation/results.md and figures/stage6_validation.{png,svg}.

Each engine runs at its published p_c, O/F and eps through
design.evaluate_performance (CEA frozen and shifting, then lambda, eta_c*
and eta_CF; theory sections 2 and 4). No nozzle contour is assumed. lambda
comes from the SP-125 Fig.4-12 digitization only where the chart applies
(Merlin, eps=16); RS-25 and RL10B-2 lie outside it, so lambda is a labeled
assumption, swept over 0.985-0.995. Charts are never extrapolated.
Engines marked illustrative (Merlin: unsourced O/F) are reported as a band
over their O/F range and are not judged pass/fail.
"""
from dataclasses import dataclass, field, replace
from pathlib import Path
import warnings

import yaml

from engine_design.constants import (COMBUSTION_EFFICIENCY, LBF_TO_N, NOZZLE_EFFICIENCY,
                                     PSI_TO_PA, STANDARD_ATMOSPHERE_PA,
                                     VALIDATION_DIVERGENCE_ASSUMED, VALIDATION_DIVERGENCE_SWEEP)
from engine_design.design import evaluate_performance
from engine_design.losses import bell_divergence_factor
from engine_design.thermo import ThermoWarning

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).with_name("engines.yaml")
# Target error bands for DELIVERED vacuum Isp, upper ends of the approximate
# bands in docs/spec/03_VALIDATION.md ("within ~2%", "~2-3%"). The spec's
# Merlin band (~3-5%) is not applied: its inputs make it illustrative only.
TARGET_PERCENT = {"RS-25": 2.0, "RL10B-2": 3.0}
# SI factors for the units engines.yaml stores exactly as published.
_UNITS = {"psia": PSI_TO_PA, "MPa": 1e6, "lbf": LBF_TO_N, "kN": 1e3, "s": 1.0, "-": 1.0}


@dataclass(frozen=True)
class Engine:
    """Published engine inputs in SI, with each value's quality label."""
    name: str
    propellants: str
    cycle: str
    p_c_Pa: float
    OF: float
    eps: float
    isp_vac_s: float
    isp_sl_s: float | None
    thrust_vac_N: float | None
    thrust_sl_N: float | None
    bell_fraction: float | None
    quality: dict = field(compare=False)
    illustrative: bool = False
    OF_band: tuple | None = None


def _si(entry):
    if entry is None or entry.get("value") is None:
        return None
    return float(entry["value"]) * _UNITS[entry["unit"]]


def load_engines(path=DATA):
    """Read engines.yaml and convert every published value to SI."""
    raw = yaml.safe_load(Path(path).read_text())
    engines = []
    for e in raw["engines"]:
        keys = ("p_c", "mixture_ratio", "area_ratio", "isp_vac", "isp_sl", "thrust_vac", "thrust_sl")
        isp_sl = _si(e["isp_sl"])
        if isp_sl is None and e["isp_sl"].get("quality") == "derived":
            # Same mass flow at both conditions: Isp_sl = Isp_vac * F_sl / F_vac.
            isp_sl = _si(e["isp_vac"]) * _si(e["thrust_sl"]) / _si(e["thrust_vac"])
        engines.append(Engine(
            e["name"], e["propellants"], e["cycle"], _si(e["p_c"]), _si(e["mixture_ratio"]),
            _si(e["area_ratio"]), _si(e["isp_vac"]), isp_sl, _si(e["thrust_vac"]),
            _si(e["thrust_sl"]), e["divergence"].get("bell_fraction"),
            {k: e[k].get("quality") for k in keys}, bool(e.get("illustrative", False)),
            tuple(e["mixture_ratio_band"]) if "mixture_ratio_band" in e else None))
    return engines


@dataclass(frozen=True)
class Prediction:
    """Ideal frozen, ideal shifting and delivered Performance for one case."""
    engine: Engine
    divergence_factor: float
    divergence_source: str
    eta_c_star: float
    eta_CF: float
    frozen: object
    shifting: object
    delivered: object
    band: dict | None = None

    def error_percent(self, isp_s, published_s):
        return 100 * (isp_s / published_s - 1)

    @property
    def vac_error_percent(self):
        return self.error_percent(self.delivered.Isp_vac_s, self.engine.isp_vac_s)

    def isp_range(self, model, condition="vac"):
        """(min, max) Isp over the O/F band, or (value, value) without one."""
        if self.band is not None:
            return self.band[(model, condition)]
        value = getattr(getattr(self, model), f"Isp_{condition if condition == 'vac' else 'amb'}_s")
        return value, value


def divergence(engine):
    """Return (lambda, source); the chart only where it applies, else assumption."""
    if engine.bell_fraction is not None:
        return (bell_divergence_factor(engine.eps, engine.bell_fraction),
                f"SP-125 Fig.4-12, eps={engine.eps:g}, assumed {engine.bell_fraction:.0%} bell")
    return VALIDATION_DIVERGENCE_ASSUMED, "assumption (outside Fig.4-12; swept 0.985-0.995)"


def predict(engine, *, divergence_factor=None, eta_c_star=COMBUSTION_EFFICIENCY,
            eta_CF=NOZZLE_EFFICIENCY):
    """Run one engine at its published operating point; return Prediction."""
    if divergence_factor is None:
        lam, source = divergence(engine)
    else:
        lam, source = divergence_factor, "sensitivity case"
    def run(OF):
        with warnings.catch_warnings():
            # Band warnings describe nominal equilibrium; real engines may sit outside.
            warnings.simplefilter("ignore", ThermoWarning)
            ideal, frozen, delivered = evaluate_performance(
                engine.propellants, engine.p_c_Pa, OF, engine.eps,
                STANDARD_ATMOSPHERE_PA, lam, eta_c_star, eta_CF)
        return {"frozen": frozen, "shifting": ideal, "delivered": delivered}

    nominal = run(engine.OF)
    band = None
    if engine.OF_band is not None:
        cases = [nominal] + [run(of) for of in engine.OF_band]
        band = {(m, c): (min(getattr(r[m], f"Isp_{a}_s") for r in cases),
                         max(getattr(r[m], f"Isp_{a}_s") for r in cases))
                for m in nominal for c, a in (("vac", "vac"), ("sl", "amb"))}
    return Prediction(engine, lam, source, eta_c_star, eta_CF, nominal["frozen"],
                      nominal["shifting"], nominal["delivered"], band)


def throat_diameter_m(prediction):
    """D_t from A_t = F_vac/(CF_vac,del*p_c) (theory section 3); None without thrust."""
    from math import pi, sqrt
    e = prediction.engine
    if e.thrust_vac_N is None:
        return None
    return 2 * sqrt(e.thrust_vac_N / (prediction.delivered.CF_vac * e.p_c_Pa) / pi)


def sensitivity(engine):
    """Return [(case, delivered Isp_vac, change from baseline in s)].

    Spec study for every engine: +/-5% p_c, +/-0.2 O/F, +/-1% eta_c*.
    Plus lambda 0.985-0.995 where it is an assumption, RS-25 eps=77.5, and
    RL10B-2 eta_CF 0.97-0.99 (large-eps boundary-layer and kinetic losses).
    """
    base = predict(engine)
    b = base.delivered.Isp_vac_s
    cases = [("baseline", engine, {})]
    for f in (0.95, 1.05):
        cases.append((f"p_c x{f:g}", replace(engine, p_c_Pa=engine.p_c_Pa * f), {}))
    for d in (-0.2, 0.2):
        cases.append((f"O/F {d:+g}", replace(engine, OF=engine.OF + d), {}))
    for f in (0.99, 1.01):
        cases.append((f"eta_c* x{f:g}", engine, {"eta_c_star": COMBUSTION_EFFICIENCY * f}))
    if engine.bell_fraction is None:
        for lam in VALIDATION_DIVERGENCE_SWEEP:
            if lam != VALIDATION_DIVERGENCE_ASSUMED:
                cases.append((f"lambda {lam:g}", engine, {"divergence_factor": lam}))
    if engine.name == "RS-25":
        cases.append(("eps 77.5 (alt. source)", replace(engine, eps=77.5), {}))
    if engine.name == "RL10B-2":
        for eta in (0.97, 0.99):
            cases.append((f"eta_CF {eta:g}", engine, {"eta_CF": eta}))
        cases.append(("p_c 644 psia (NRC table)", replace(engine, p_c_Pa=644 * PSI_TO_PA), {}))
    rows = []
    for label, eng, kwargs in cases:
        isp = predict(eng, **kwargs).delivered.Isp_vac_s
        rows.append((label, isp, isp - b))
    return rows


def _fmt(x, spec="6.1f"):
    return "—" if x is None else format(x, spec)


def _span(low, high, spec):
    """'x' for a point, 'low–high' for a band."""
    joint = " to " if spec.startswith("+") else "–"  # "+2.47 to +3.01", not "+2.47–+3.01"
    return format(low, spec) if f"{low:{spec}}" == f"{high:{spec}}" else \
        f"{low:{spec}}{joint}{high:{spec}}"


def _verdict(p, condition="vac"):
    """Pass/fail text against TARGET_PERCENT; illustrative engines get none."""
    if p.engine.illustrative:
        lo, hi = p.engine.OF_band
        return f"illustrative (O/F {lo:g}–{hi:g}, unsourced)"
    target = TARGET_PERCENT[p.engine.name]
    return f"±{target:g} {'✅' if abs(p.vac_error_percent) <= target else '❌'}"


def markdown(predictions, sensitivities):
    """Results tables in Markdown; numbers only, discussion lives in README."""
    lines = ["# Stage 6 validation results", "",
             "Generated by `python validation/run_validation.py`. Vacuum Isp in seconds; "
             "error is delivered vs published. η_c* = η_CF = 0.98 (assumptions). "
             "Ranges are over an illustrative engine's O/F band.", "",
             "| Engine | Cycle | Published | Predicted (frozen) | Predicted (shifting) "
             "| Delivered (with η) | % error | Target | Data quality |",
             "|---|---|---|---|---|---|---|---|---|"]
    for p in predictions:
        e = p.engine
        quality = ", ".join(f"{k} {v}" for k, v in e.quality.items()
                            if v not in (None, "primary") and k != "isp_sl")
        d_lo, d_hi = p.isp_range("delivered")
        err = _span(p.error_percent(d_lo, e.isp_vac_s), p.error_percent(d_hi, e.isp_vac_s), "+.2f")
        lines.append(f"| {e.name} | {e.cycle} | {e.isp_vac_s:.1f} | "
                     f"{_span(*p.isp_range('frozen'), '.1f')} | {_span(*p.isp_range('shifting'), '.1f')} | "
                     f"{_span(d_lo, d_hi, '.1f')} | {err} | {_verdict(p)} | {quality or 'all primary'} |")
    lines += ["", "## Sea level (attached flow, 101.325 kPa)", "",
              "| Engine | Published | Delivered | % error | Note |", "|---|---|---|---|---|"]
    for p in predictions:
        e = p.engine
        if e.isp_sl_s is None:
            continue
        note = ("derived: Isp_vac·F_sl/F_vac (engines.yaml)" if e.quality["isp_sl"] == "derived"
                else e.quality["isp_sl"])
        lo, hi = p.isp_range("delivered", "sl")
        err = _span(p.error_percent(lo, e.isp_sl_s), p.error_percent(hi, e.isp_sl_s), "+.2f")
        lines.append(f"| {e.name} | {e.isp_sl_s:.1f} | {_span(lo, hi, '.1f')} | {err} | {note} |")
    lines += ["", "## Model inputs", "",
              "| Engine | p_c [MPa] | O/F | ε | λ | λ source | Implied Isp_pub/Isp_shifting | Predicted D_t [mm] |",
              "|---|---|---|---|---|---|---|---|"]
    for p in predictions:
        e, dt = p.engine, throat_diameter_m(p)
        of = e.OF if e.OF_band is None else f"{e.OF:g} ({e.OF_band[0]:g}–{e.OF_band[1]:g})"
        lines.append(f"| {e.name} | {e.p_c_Pa/1e6:.3f} | {of} | {e.eps:g} | "
                     f"{p.divergence_factor:.4f} | {p.divergence_source} | "
                     f"{e.isp_vac_s / p.shifting.Isp_vac_s:.4f} | {_fmt(dt and dt*1e3)} |")
    lines += ["", "Predicted D_t uses published vacuum thrust; no primary source publishes "
              "a throat diameter to compare against. Merlin thrust is from a different engine "
              "version than its Isp, so its D_t is context only.", ""]
    for name, rows in sensitivities.items():
        lines += [f"## Sensitivity: {name}", "", "| Case | Delivered Isp_vac [s] | Change [s] |",
                  "|---|---|---|"]
        lines += [f"| {label} | {isp:.2f} | {d:+.2f} |" for label, isp, d in rows]
        lines.append("")
    return "\n".join(lines)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from engine_design.plots import plot_validation

    engines = load_engines()
    predictions = [predict(e) for e in engines]
    sensitivities = {e.name: sensitivity(e) for e in engines}
    for p in predictions:
        if p.isp_range("shifting")[0] <= p.engine.isp_vac_s:
            # 03_VALIDATION.md: ideal shifting below published means a model error.
            warnings.warn(f"{p.engine.name}: ideal shifting Isp is below the published value")
    text = markdown(predictions, sensitivities)
    (Path(__file__).with_name("results.md")).write_text(text + "\n")
    print(text)
    fig = plot_validation(predictions, TARGET_PERCENT)
    output = ROOT / "figures"
    output.mkdir(exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(output / f"stage6_validation.{extension}", dpi=180,
                    metadata={"Date": None} if extension == "svg" else None)
    plt.close(fig)
    print(f"\nWrote validation/results.md and figures/stage6_validation.png")


if __name__ == "__main__":
    main()
