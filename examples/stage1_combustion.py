"""Reproduce Stage 1 plots and cross-checks: python examples/stage1_combustion.py."""
from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from engine_design.constants import DEFAULT_OF_BOUNDS
from engine_design.plots import plot_OF_sweeps
from engine_design.thermo import ThermoWarning, chamber_props, ideal_c_star, optimal_OF, sweep_OF


def main():
    output = Path(__file__).resolve().parents[1] / "figures"
    output.mkdir(exist_ok=True)
    pressure = 10e6
    sweeps = {}
    for pair, bounds in DEFAULT_OF_BOUNDS.items():
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ThermoWarning)
            sweeps[pair] = sweep_OF(pair, pressure, np.linspace(*bounds, 81))
        print(f"{pair}: {len(caught)}/81 samples outside nominal spec bands; "
              "full sweep intentionally includes off-design mixtures.")
    fig = plot_OF_sweeps(sweeps, pressure)
    fig.savefig(output / "stage1_combustion.png", dpi=180)
    fig.savefig(output / "stage1_combustion.svg")
    plt.close(fig)
    for pair, of in [("LOX/LH2", 6), ("LOX/RP-1", 2.5), ("LOX/CH4", 3.5)]:
        p = chamber_props(pair, pressure, of)
        approx = ideal_c_star(p.T_c_K, p.M_kg_kmol, p.gamma)
        optimum = optimal_OF(pair, pressure, 40, "vacuum")
        print(f"{pair}: nominal O/F={of:g}, Tc={p.T_c_K:.2f} K, "
              f"c*={p.c_star_m_s:.2f} m/s, perfect-gas difference="
              f"{100 * (approx / p.c_star_m_s - 1):+.4f}%, "
              f"vacuum optimum O/F={optimum:.4f} (eps=40)")
    print(f"Saved figures to {output}")


if __name__ == "__main__":
    main()
