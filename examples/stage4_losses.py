"""Reproduce Stage 4 losses and frozen/shifting comparison: python examples/stage4_losses.py."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from engine_design import EngineRequirements, design_engine
from engine_design.losses import conical_divergence_factor
from engine_design.plots import plot_losses


def main():
    output = Path(__file__).resolve().parents[1] / "figures"
    output.mkdir(exist_ok=True)
    print(f"Conical lambda(15°) = {conical_divergence_factor(15):.5f}")
    # The project brief's example requirement; not a real engine.
    design = design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                                              thrust_condition="vacuum", p_c_Pa=10e6))
    print(design.summary())
    ideal = design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                                             thrust_condition="vacuum", p_c_Pa=10e6,
                                             divergence_factor=1, eta_c_star=1, eta_CF=1))
    g, gi = design.geometry, ideal.geometry
    print(f"\nLosses enlarge the throat by {g.A_t_m2 / gi.A_t_m2 - 1:+.2%} "
          f"(= CF_ideal/CF_delivered - 1) and mass flow by {g.mdot_kg_s / gi.mdot_kg_s - 1:+.2%}.")
    fig = plot_losses(design)
    for extension in ("png", "svg"):
        fig.savefig(output / f"stage4_losses.{extension}", dpi=180)
    plt.close(fig)
    print(f"Saved figures to {output}")


if __name__ == "__main__":
    main()
