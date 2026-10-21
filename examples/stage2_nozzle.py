"""Reproduce Stage 2: python examples/stage2_nozzle.py."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from engine_design.constants import STANDARD_ATMOSPHERE_PA
from engine_design.nozzle import (mach_from_area_ratio, pressure_ratio,
                                  thrust_coefficient, optimal_eps, isp)
from engine_design.plots import plot_nozzle_performance


def main():
    # Analyst-defined demonstration gas and operating point, not engine data.
    gamma, c_star, pc, eps = 1.2, 1800.0, 10e6, 40.0
    output = Path(__file__).resolve().parents[1] / "figures"
    output.mkdir(exist_ok=True)
    fig = plot_nozzle_performance(gamma, c_star, pc, eps)
    for extension in ("png", "svg"):
        fig.savefig(output / f"stage2_nozzle.{extension}", dpi=180)
    plt.close(fig)
    mach = mach_from_area_ratio(eps, gamma)
    pe = pc * pressure_ratio(mach, gamma)
    matched = optimal_eps(gamma, pc, STANDARD_ATMOSPHERE_PA)
    print(f"Example perfect gas: gamma={gamma}, c*={c_star:g} m/s, pc={pc/1e6:g} MPa")
    print(f"eps={eps:g}: exit Mach={mach:.5f}, exit pressure={pe:.2f} Pa")
    print(f"Sea-level matched expansion: eps={matched:.5f}")
    for name, ambient in [("vacuum", 0.0), ("sea level", STANDARD_ATMOSPHERE_PA)]:
        cf = thrust_coefficient(gamma, eps, pc, ambient)
        print(f"{name}: attached-flow CF={cf:.6f}, Isp={isp(c_star, cf):.3f} s")
    print("Dashed segments are flagged extrapolations, not separated-flow predictions.")
    print(f"Saved figures to {output}")


if __name__ == "__main__":
    main()
