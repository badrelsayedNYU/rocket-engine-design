# Build Plan: Stage by Stage

Build in order. **After each stage, stop and give me:**
1. what you built (files and functions)
2. how to run it
3. test results
4. 2–3 sentences on the physics I should be able to explain

Don't start the next stage until I say so.

---

## Stage 0: Scaffold
- Create the repo structure from `00_PROJECT_BRIEF.md`, with `pyproject.toml`, an empty README, and a pytest config.
- `constants.py`: g0 = 9.80665 m/s², R_u = 8314.462618 J/(kmol·K), and placeholders for empirical constants with source comments.
- `propellants.py`: definitions for **LOX/LH2, LOX/RP-1, LOX/CH4**, including CEA propellant names. **Verify the exact names RocketCEA accepts** (e.g., `LOX`, `LH2`, `RP-1` or `RP1`, `CH4`) by testing them; don't guess.
- **Done when:** `pip install -e .` works and `pytest` runs (even with zero tests).

## Stage 1: Combustion (`thermo.py`)
- Function: `chamber_props(propellants, p_c_Pa, OF, eps=None, frozen=False)`. Returns T_c, M, γ, c*, plus transport (μ, c_p, Pr) at chamber and throat.
- Function: `sweep_OF(propellants, p_c_Pa, OF_range)`. Returns arrays for plotting.
- Function: `optimal_OF(propellants, p_c_Pa, eps, condition)`. The mixture ratio that maximizes Isp. Note that this is usually fuel-rich of stoichiometric; explain why (lower molecular weight).
- Cross-check: ideal c* from the formula in `01_THEORY.md` vs. CEA's c*. They should agree within about 1–2%.
- **Plots:** T_c, M, γ, and c* vs. O/F for all three pairs at p_c = 10 MPa.
- **Tests:** values inside the sanity bands; ideal-vs-CEA c* agreement.

## Stage 2: Nozzle (`nozzle.py`)
- `mach_from_area_ratio(eps, gamma)`: supersonic branch, bracketed root finder.
- `pressure_ratio(M, gamma)`, `thrust_coefficient(gamma, eps, p_c, p_a)`, `isp(c_star, CF)`.
- `optimal_eps(gamma, p_c, p_a)`: ε for p_e = p_a.
- `separation_check(p_e, p_a, k_sep=0.4)`.
- **Plots:** C_F vs. ε at sea level and in vacuum (Sutton-style chart); Isp vs. altitude for a fixed nozzle.
- **Tests:** reproduce a textbook isentropic table value (e.g., γ = 1.4, M = 2 gives A/A* = 1.6875); reproduce a Sutton C_F chart point; vacuum C_F > sea-level C_F.

## Stage 3: Sizing (`sizing.py`, `design.py`)
- From requirements, compute A_t, D_t, A_e, D_e, ṁ, ṁ_ox, ṁ_f, chamber volume (L*), chamber diameter and length (contraction ratio).
- Generate the contour: converging arc, throat arcs, Rao parabolic bell (80% default). Digitized θ_n and θ_e tables, with source.
- `EngineDesign` dataclass with a `.summary()` method that prints a clean table.
- **Plots:** engine contour to scale, with key dimensions annotated.
- **Tests:** F = C_F·p_c·A_t holds when recomputed; contour is continuous with a smooth slope at the joints; exit radius matches ε.

## Stage 4: Losses (`losses.py`)
- Divergence factor (conical and bell), η_c*, η_CF, delivered Isp.
- Frozen vs. shifting equilibrium Isp from CEA. Report both.
- **Tests:** λ(15°) ≈ 0.983; delivered Isp < ideal Isp.

## Stage 5: Cooling (`cooling.py`)
- Bartz h_g, σ, T_aw, and q along the contour. Peak q and its location.
- Total heat load. Coolant-capacity ratio using fuel as coolant.
- **Plots:** q(x) along the engine, overlaid on the contour.
- **Tests:** peak q is at or slightly upstream of the throat; magnitude falls inside the sanity band.

## Stage 6: Validation (`validation/`)
- Load `engines.yaml` (data and sources in `03_VALIDATION.md`).
- For each engine: run the tool with that engine's p_c, O/F, ε, and propellants. Compare predicted vs. published vacuum Isp (plus sea-level Isp and thrust where available), and throat diameter if published.
- Output a markdown table: **engine | published | predicted (frozen) | predicted (shifting) | delivered (with η) | % error**.
- **Plot:** predicted vs. published Isp with the frozen–shifting band.
- In the README, explain *why* each discrepancy exists (e.g., the RS-25 runs a staged combustion cycle; the Merlin is a gas-generator cycle whose GG exhaust lowers effective Isp; the RL10 has a huge ε and expander-cycle details). **This explanation is the most important part of the project.**

## Stage 7: App and README
- Streamlit `app.py`: sidebar inputs (propellants, thrust, p_c, O/F or auto, nozzle option); outputs the summary table, contour plot, heat-flux plot, and C_F/Isp plots. Include a "compare to real engines" tab.
- README written results-first: validation table and contour image at the top, then method, assumptions, limitations, and references.
- Publish the completed project to GitHub and deploy the Streamlit app from that repository; add the live demo link to the README.
