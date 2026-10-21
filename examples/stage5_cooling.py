"""Reproduce Stage 5 Bartz heat flux and coolant check: python examples/stage5_cooling.py."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from engine_design import EngineRequirements, design_engine


def main():
    output = Path(__file__).resolve().parents[1] / "figures"
    output.mkdir(exist_ok=True)
    # The project brief's example requirement; not a real engine.
    design = design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                                              thrust_condition="vacuum", p_c_Pa=10e6))
    print(design.summary())
    p, rt = design.cooling.profile, design.geometry.D_t_m / 2
    print(f"\nPeak q at x/Rt = {p.x_peak_m / rt:+.3f} (upstream of the throat is negative).")
    for T_w in (600.0, 800.0, 1000.0):
        q = design_engine(EngineRequirements(propellants="LOX/CH4", thrust_N=100e3,
                                             p_c_Pa=10e6, wall_temperature_K=T_w)).cooling
        print(f"T_w = {T_w:6.0f} K: peak q = {q.profile.peak_q_W_m2/1e6:6.1f} MW/m², "
              f"Q = {q.heat_load_W/1e6:.2f} MW, capacity/load = {q.capacity_ratio:.2f}")
    fig = design.plot_heat_flux()
    for extension in ("png", "svg"):
        fig.savefig(output / f"stage5_heat_flux.{extension}", dpi=180)
    plt.close(fig)
    print(f"Saved figures to {output}")


if __name__ == "__main__":
    main()
