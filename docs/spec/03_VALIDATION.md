# Validation Data

> **To the AI assistant:** the values below are approximate starting points. **Before using any number, verify it against a primary or authoritative source** (NASA technical reports, manufacturer data sheets such as Aerojet Rocketdyne/L3Harris, SpaceX, or Sutton's tables). Record the exact source and URL for every value in `engines.yaml`. Where sources disagree, note the range. Don't invent missing values; mark them `null`.

## Engines

### RS-25 (Space Shuttle Main Engine / SLS core stage)
- Propellants: LOX/LH2
- Cycle: fuel-rich **staged combustion**
- Chamber pressure: ~20.6 MPa (~2,990 psia) at 109% rated power level (verify the power level the figure refers to)
- Mixture ratio (O/F): ~6.0
- Area ratio ε: ~69
- Vacuum Isp: ~452 s; sea-level Isp: ~366 s
- Vacuum thrust: ~2,280 kN (109%)

### RL10 (Centaur / upper stages)
- Propellants: LOX/LH2
- Cycle: **expander**
- Pick one variant and stay consistent. Suggested: **RL10B-2** (extendible nozzle, very large ε)
- Chamber pressure: ~4.4 MPa (~630 psia)
- Mixture ratio: ~5.9
- Area ratio ε: ~280 (with the nozzle extension deployed)
- Vacuum Isp: ~465 s
- Vacuum thrust: ~110 kN

### Merlin 1D (Falcon 9 first stage, sea-level version)
- Propellants: LOX/RP-1
- Cycle: **gas generator**
- Chamber pressure: ~9.7 MPa (published figures vary; verify)
- Mixture ratio: ~2.3–2.4
- Area ratio ε: ~16
- Sea-level Isp: ~282 s; vacuum Isp: ~311 s
- Sea-level thrust: ~845 kN (varies by block; verify)

### Optional stretch: Raptor (LOX/CH4, full-flow staged combustion)
Add only if reliable public data exists. Public figures are less authoritative. Label them clearly as such.

## Target accuracy

Ideal 1D CEA analysis overpredicts real engines. These are reasonable targets for **delivered** vacuum Isp (after η_c*, η_CF, and divergence losses):

| Engine | Target error band |
|---|---|
| RS-25 | within ~2% |
| RL10B-2 | within ~2–3% |
| Merlin 1D | within ~3–5% (the gas-generator exhaust isn't modeled; explain this) |

Ideal (no-loss) shifting-equilibrium Isp should come out **above** published values. If it comes out below, something is wrong.

## Expected discussion points (for the README)

- **Frozen vs. shifting:** real performance falls between them; explain why (recombination during expansion releases energy).
- **Cycle effects:** a gas-generator cycle dumps turbine exhaust at low Isp, so the engine's overall Isp is lower than the main chamber's. Staged combustion and expander cycles burn all the propellant in the main chamber.
- **Large-ε effects (RL10):** boundary-layer and kinetic losses grow; check how sensitive the result is to η_CF.
- **Sensitivity study:** for one engine, show how ±5% in p_c, ±0.2 in O/F, and ±1% in η_c* move the predicted Isp. That shows which inputs matter most.
