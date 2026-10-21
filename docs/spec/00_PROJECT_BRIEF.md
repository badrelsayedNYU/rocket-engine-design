# Liquid Rocket Engine Design Tool: Project Brief

> **To the AI assistant building this:** read all five files (00–04) before writing any code. Build one stage at a time (see `02_STAGES.md`). Stop after each stage so I can review it, run it, and make sure I understand the physics. I need to be able to defend every number in an interview, so explain your choices as you go, keep the code readable, and never hide the physics inside a black box.

## 1. What this is

A Python tool that designs a liquid bipropellant rocket engine from top-level requirements and then validates the design against real engines.

**Input (example):** "100 kN vacuum thrust, LOX/methane, 10 MPa chamber pressure, sea-level-optimized nozzle."

**Output:** a full preliminary design:
- combustion properties (chamber temperature, molecular weight, γ, c*)
- nozzle performance (exit Mach, thrust coefficient, Isp at sea level and in vacuum, separation check)
- geometry (throat/exit diameters, chamber volume and length, bell nozzle contour)
- loss-corrected performance
- throat heat flux and a regenerative cooling feasibility check
- plots, plus a validation table comparing predictions to the RS-25, RL10, and Merlin

## 2. Who I am and why I'm building it

I'm Badr Elsayed, a mechanical engineering senior (aerospace minor) at NYU Tandon, applying to aerospace M.S. programs for Fall 2027. I've taken compressible flow, propulsion, heat transfer, and thermodynamics. I've already built and published three validated Python solvers: a vortex panel method, a 1D Euler shock tube solver, and a reentry heat-shield model.

This project should show **propulsion system design and validation**. It is **not** meant to show CFD. I want to work as a propulsion/aerodynamics engineer on launch vehicles, not as a CFD specialist. Keep everything 0D/1D, built on analytical and thermochemical models.

## 3. Principles (non-negotiable)

1. **Physics first, transparent.** Every equation in the code must match `01_THEORY.md`, with a docstring citing the equation and its assumptions.
2. **Validate everything.** Each module needs tests against a known answer: hand calculation, textbook example, or published engine data.
3. **SI units internally.** Convert only at input/output. Name variables with units where it helps (`p_c_Pa`, `T_c_K`).
4. **No magic numbers.** Every empirical constant (L*, Bartz coefficient, separation criterion, efficiencies) lives in one `constants.py` or config file, with a source comment.
5. **Small, readable modules.** I'll read all of it.
6. **Reproducible.** One command reproduces every figure and the validation table.

## 4. Tech stack

- Python 3.11+
- **Thermochemistry:** `rocketcea` (NASA CEA wrapper) as the primary engine. If it won't install, fall back to `cantera` with an equilibrium calculation. Hide this behind one interface (`thermo.py`) so the backend can swap.
- `numpy`, `scipy` (root finding), `matplotlib`
- `pytest` for tests
- **Web app (final stage):** Streamlit, matching the style of my existing airfoil tool
- Packaging: `pyproject.toml`, installable with `pip install -e .`

## 5. Repository layout

```
rocket-engine-design/
├── README.md                 # results-first: validation table + key plots at the top
├── pyproject.toml
├── src/engine_design/
│   ├── __init__.py
│   ├── constants.py          # g0, R_u, empirical constants with sources
│   ├── propellants.py        # propellant pair definitions, CEA names
│   ├── thermo.py             # Stage 1: CEA interface (Tc, M, gamma, c*, transport)
│   ├── nozzle.py             # Stage 2: isentropic nozzle, CF, Isp, separation
│   ├── sizing.py             # Stage 3: throat/exit, mdot, chamber, contour
│   ├── losses.py             # Stage 4: efficiencies, frozen vs shifting
│   ├── cooling.py            # Stage 5: Bartz heat flux along the nozzle
│   ├── design.py             # top-level: requirements -> EngineDesign dataclass
│   └── plots.py
├── validation/
│   ├── engines.yaml          # published engine data (see 03_VALIDATION.md)
│   └── run_validation.py     # produces the validation table + plots
├── tests/
├── examples/                 # example designs as scripts
├── figures/
└── app.py                    # Stage 7: Streamlit app
```

## 6. Top-level API (target)

```python
from engine_design import EngineRequirements, design_engine

req = EngineRequirements(
    propellants="LOX/CH4",
    thrust_N=100e3,
    thrust_condition="vacuum",      # or "sea_level"
    p_c_Pa=10e6,
    mixture_ratio=None,             # None -> optimize for max Isp
    nozzle="sea_level_optimized",   # or a fixed area ratio, e.g. eps=40
    p_amb_Pa=101325,
)
eng = design_engine(req)
print(eng.summary())                # human-readable table
eng.plot_contour(); eng.plot_heat_flux()
```

## 7. Definition of done

- All seven stages in `02_STAGES.md` are complete, with passing tests.
- `python validation/run_validation.py` reproduces the validation table. Predicted vacuum Isp is within the target bands in `03_VALIDATION.md` for the RS-25, RL10, and Merlin, and every discrepancy is explained in the README.
- README written results-first, with figures.
- Streamlit app runs locally with `streamlit run app.py`.
