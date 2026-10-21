# Theory Reference

Implement these equations as written. Assumptions are stated for each. Main references:
- Sutton & Biblarz, *Rocket Propulsion Elements* (9th ed.), Ch. 3 (nozzle theory), Ch. 5 (thermochemistry), Ch. 8 (thrust chambers)
- Huzel & Huang, *Modern Engineering for Design of Liquid-Propellant Rocket Engines* (AIAA)
- Bartz, D. R. (1957), "A Simple Equation for Rapid Estimation of Rocket Nozzle Convective Heat Transfer Coefficients," *Jet Propulsion* 27(1)
- Rao, G. V. R. (1958/1960), thrust-optimized bell nozzle contours

**Notation:** p_c chamber (stagnation) pressure, T_c chamber temperature, γ ratio of specific heats, M molecular weight, R = R_u/M, ε = A_e/A_t area ratio, g0 = 9.80665 m/s².

---

## 1. Combustion (Stage 1)

Use NASA CEA (through RocketCEA) for T_c, M, γ, c*, and transport properties (μ, c_p, Pr) at the chamber and throat. CEA assumes adiabatic, equilibrium combustion at constant pressure.

**Ideal c\*** (use as a cross-check against CEA's c*):

$$c^* = \frac{\sqrt{\gamma R T_c}}{\gamma\,\sqrt{\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}}}$$

---

## 2. Isentropic nozzle (Stage 2)

Assumptions: steady, 1D, isentropic, calorically perfect gas with a single (frozen) γ, choked throat. Stage 4 adds the equilibrium (shifting) comparison from CEA.

**Area–Mach relation** (solve for M_e given ε, supersonic branch, using a bracketed root finder):

$$\frac{A}{A^*} = \frac{1}{M}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}M^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}$$

**Pressure ratio:**

$$\frac{p_e}{p_c} = \left(1+\frac{\gamma-1}{2}M_e^2\right)^{-\frac{\gamma}{\gamma-1}}$$

**Thrust coefficient:**

$$C_F = \sqrt{\frac{2\gamma^2}{\gamma-1}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}\left[1-\left(\frac{p_e}{p_c}\right)^{\frac{\gamma-1}{\gamma}}\right]} + \frac{p_e-p_a}{p_c}\,\varepsilon$$

- Vacuum: p_a = 0. Sea level: p_a = 101325 Pa.
- Optimum expansion: p_e = p_a, so the pressure term is zero.

**Performance:**

$$F = C_F\,p_c\,A_t \qquad \dot m = \frac{p_c A_t}{c^*} \qquad I_{sp} = \frac{c^* C_F}{g_0} \qquad c_{eff} = c^* C_F$$

**Flow separation check (overexpanded nozzle):** flag likely separation when

$$p_e < k_{sep}\, p_a, \qquad k_{sep} \approx 0.25\text{–}0.4 \;\;(\text{Summerfield criterion; default } 0.4)$$

Report a warning; don't try to model the separated flow.

---

## 3. Sizing (Stage 3)

**Throat area from required thrust:**

$$A_t = \frac{F}{C_F\,p_c}$$

(Use the C_F at the specified thrust condition.)

**Exit area:** A_e = ε·A_t. For a "sea-level-optimized" nozzle, choose ε so that p_e = p_a.

**Propellant flows:** total mass flow ṁ from above, then

$$\dot m_{ox} = \dot m\,\frac{O/F}{1+O/F}, \qquad \dot m_{f} = \dot m\,\frac{1}{1+O/F}$$

**Chamber volume, using characteristic length L\*:**

$$V_c = L^*\,A_t$$

Typical L* (Huzel & Huang / Sutton; put these in `constants.py` with sources):

| Propellants | L* (m) |
|---|---|
| LOX/LH2 | 0.76–1.02 |
| LOX/RP-1 | 1.02–1.27 |
| LOX/CH4 | ~0.9–1.3 (less data; treat as an assumption) |

**Contraction ratio** ε_c = A_c/A_t: larger for small engines, smaller for big ones. Default to an empirical correlation or let the user set it (typical 2–8). Chamber length then follows from V_c, A_c, and the converging-section geometry.

**Nozzle contour:**
- Converging section: arc of radius 1.5·R_t upstream of the throat (default convergence half-angle 30–45°).
- Diverging section: arc of radius 0.382·R_t downstream of the throat, then a **Rao thrust-optimized parabola** (default 80% bell length relative to a 15° cone).
- Initial and exit wall angles (θ_n, θ_e) come from Rao's charts as a function of ε and bell %. Tabulate digitized values in `constants.py` and interpolate. Cite the source of the digitized data.
- 15° cone length for comparison: L_cone = (R_e − R_t)/tan(15°).

---

## 4. Losses (Stage 4)

**Divergence (2D flow) loss.** Conical nozzle with half-angle α:

$$\lambda = \frac{1+\cos\alpha}{2}$$

For a bell nozzle use λ ≈ 0.985–0.99 (or estimate from the exit angle θ_e). Apply to the momentum term of C_F.

**Combustion efficiency:** η_c* = c*_actual / c*_ideal, typically 0.96–0.99. User input; default 0.98.

**Nozzle efficiency:** η_CF (beyond divergence: boundary layer, kinetics), typically 0.97–0.99.

**Delivered values:**

$$I_{sp,del} = \eta_{c^*}\,\eta_{C_F}\,I_{sp,ideal}$$

**Frozen vs. shifting equilibrium:** report both CEA values (RocketCEA supports frozen and equilibrium expansion). Real engines usually fall between the two. Show this band on the validation plot.

---

## 5. Cooling: Bartz correlation (Stage 5)

Gas-side convective heat transfer coefficient (SI form, all inputs at chamber stagnation conditions unless noted):

$$h_g = \frac{0.026}{D_t^{0.2}}\left(\frac{\mu^{0.2}c_p}{Pr^{0.6}}\right)_{c}\left(\frac{p_c}{c^*}\right)^{0.8}\left(\frac{D_t}{r_c}\right)^{0.1}\left(\frac{A_t}{A}\right)^{0.9}\sigma$$

where r_c is the throat radius of curvature (use the mean of the upstream and downstream radii), and σ is the property-variation correction:

$$\sigma = \frac{1}{\left[\frac{1}{2}\frac{T_w}{T_c}\left(1+\frac{\gamma-1}{2}M^2\right)+\frac{1}{2}\right]^{0.68}\left[1+\frac{\gamma-1}{2}M^2\right]^{0.12}}$$

**Adiabatic wall temperature,** with recovery factor r = Pr^{1/3} (turbulent):

$$T_{aw} = T_c\,\frac{1+r\frac{\gamma-1}{2}M^2}{1+\frac{\gamma-1}{2}M^2}$$

**Heat flux:**

$$q = h_g\,(T_{aw}-T_w)$$

Compute q along the contour (M from the local area ratio). Peak is near the throat. Assume a wall temperature T_w (user input; default ~800 K for copper-alloy liners, cited as an assumption).

**Regen feasibility (first-order):** total heat load Q = ∫ q dA over the wetted wall. Compare it with the coolant's capacity to absorb heat, ṁ_coolant·c_p·ΔT_allowable (fuel as coolant). Report the ratio. This is a sanity check, not a channel design.

Bartz typically **overpredicts** throat heat flux somewhat. State that in the docs.

---

## 6. Sanity bands (use in tests)

- c* for LOX/LH2 is about 2300–2400 m/s; LOX/RP-1 about 1750–1850 m/s; LOX/CH4 about 1800–1900 m/s (equilibrium, high p_c).
- γ in the chamber is about 1.13–1.25.
- C_F,vac is about 1.7–2.0 for typical upper-stage area ratios.
- Throat heat flux is tens of MW/m² for high-p_c engines (up to roughly 100+ MW/m² for RS-25-class).

If an output falls outside these bands, raise a warning and investigate before moving on.
