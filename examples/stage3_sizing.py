"""Reproduce Stage 3 sizing and contour: python examples/stage3_sizing.py."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from engine_design import EngineRequirements, design_engine


def main():
    output = Path(__file__).resolve().parents[1] / "figures"
    output.mkdir(exist_ok=True)
    # The project brief's example requirement; not a real engine.
    req = EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                             thrust_condition="vacuum", p_c_Pa=10e6,
                             mixture_ratio=None, nozzle="sea_level_optimized",
                             p_amb_Pa=101325)
    design = design_engine(req)
    print(design.summary())
    g = design.geometry
    recomputed = design.delivered.CF_vac * req.p_c_Pa * g.A_t_m2
    print(f"\nCheck: delivered CF*pc*At = {recomputed:.6f} N (requirement {req.thrust_N:g} N)")
    fig = design.plot_contour()
    for extension in ("png", "svg"):
        fig.savefig(output / f"stage3_contour.{extension}", dpi=180)
    plt.close(fig)
    print(f"Saved figures to {output}")


if __name__ == "__main__":
    main()
