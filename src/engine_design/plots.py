"""Plots of computed properties; physics remains in thermo.py."""


def plot_OF_sweeps(sweeps, p_c_Pa):
    """Return a Matplotlib figure with four chamber-property panels."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    fields = [("T_c_K", "Chamber temperature [K]"),
              ("M_kg_kmol", "Molecular weight [kg/kmol]"),
              ("gamma", "CEA equilibrium isentropic exponent [–]"),
              ("c_star_m_s", "Characteristic velocity c* [m/s]")]
    for ax, (field, label) in zip(axes.flat, fields):
        for name, sweep in sweeps.items():
            ax.plot(sweep.OF, getattr(sweep, field), label=name, linewidth=2)
        ax.set(xlabel="Oxidizer/fuel mass ratio [–]", ylabel=label)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(frameon=False)
    fig.suptitle(f"Equilibrium combustion at {p_c_Pa / 1e6:g} MPa\n"
                 "RocketCEA · default liquid reactant cards · ideal chamber", fontsize=15)
    return fig


def plot_nozzle_performance(gamma, c_star_m_s, p_c_Pa, fixed_eps):
    """Return CF-vs-area and Isp-vs-altitude plots of the Stage 2 model.

    Fixed user-supplied perfect-gas properties; no chemistry or loss fitting.
    Dashed orange segments identify the Summerfield-flagged extrapolation.
    Physics and assumptions are documented in nozzle.py and atmosphere.py.
    """
    import warnings
    import numpy as np
    import matplotlib.pyplot as plt
    from .constants import STANDARD_ATMOSPHERE_PA, SEPARATION_PRESSURE_RATIO
    from .nozzle import (NozzleWarning, mach_from_area_ratio, pressure_ratio,
                         thrust_coefficient, optimal_eps, isp)
    from .atmosphere import ambient_pressure

    ratios = np.geomspace(1, 100, 300)
    heights = np.linspace(0, 50000, 301)
    pressures = np.array([ambient_pressure(h) for h in heights])
    exit_pressures = np.array([p_c_Pa * pressure_ratio(mach_from_area_ratio(e, gamma), gamma)
                              for e in ratios])
    fixed_pe = p_c_Pa * pressure_ratio(mach_from_area_ratio(fixed_eps, gamma), gamma)
    area_flag = exit_pressures < SEPARATION_PRESSURE_RATIO * STANDARD_ATMOSPHERE_PA
    height_flag = fixed_pe < SEPARATION_PRESSURE_RATIO * pressures
    # The flagged segments are explicitly displayed; avoid hundreds of warnings.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", NozzleWarning)
        cf_vac = [thrust_coefficient(gamma, e, p_c_Pa, 0) for e in ratios]
        cf_sl = np.array([thrust_coefficient(gamma, e, p_c_Pa, STANDARD_ATMOSPHERE_PA)
                          for e in ratios])
        isp_alt = np.array([isp(c_star_m_s, thrust_coefficient(gamma, fixed_eps, p_c_Pa, p))
                            for p in pressures])
        vacuum_isp = isp(c_star_m_s, thrust_coefficient(gamma, fixed_eps, p_c_Pa, 0))
        matched = optimal_eps(gamma, p_c_Pa, STANDARD_ATMOSPHERE_PA)
        matched_cf = thrust_coefficient(gamma, matched, p_c_Pa, STANDARD_ATMOSPHERE_PA)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    ax = axes[0]
    ax.plot(ratios, cf_vac, color="#2367a2", label="Vacuum")
    ax.plot(ratios, np.where(area_flag, np.nan, cf_sl), color="#d66820", label="Sea level")
    ax.plot(ratios, np.where(area_flag, cf_sl, np.nan), color="#d66820", linestyle="--",
            label="Likely separation: extrapolation")
    ax.scatter([matched], [matched_cf], color="black", s=28, zorder=3)
    ax.annotate(f"Matched expansion\nε = {matched:.2f}", (matched, matched_cf),
                xytext=(10, -45), textcoords="offset points", fontsize=9,
                arrowprops={"arrowstyle": "->", "color": "0.4"})
    ax.set(xscale="log", xlabel="Exit/throat area ratio ε [–]",
           ylabel="Thrust coefficient CF [–]", title="Area ratio and ambient pressure")
    ax.legend(fontsize=8, loc="lower left", frameon=False)
    ax = axes[1]
    ax.plot(heights / 1000, np.where(height_flag, np.nan, isp_alt), color="#2367a2",
            label="Attached-flow model")
    ax.plot(heights / 1000, np.where(height_flag, isp_alt, np.nan), color="#d66820", linestyle="--",
            label="Likely separation: extrapolation")
    ax.axhline(vacuum_isp, color="0.4", linestyle=":", label=f"Vacuum limit: {vacuum_isp:.1f} s")
    ax.set(xlabel="Geopotential height [km]", ylabel="Ideal specific impulse [s]",
           title=f"Fixed nozzle ε = {fixed_eps:g} · U.S. Standard Atmosphere")
    ax.legend(fontsize=8, loc="lower right", frameon=False)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle(f"Isentropic nozzle · γ = {gamma:g} · pc = {p_c_Pa/1e6:g} MPa · c* = {c_star_m_s:g} m/s\n"
                 "Constant properties; no losses or separated-flow correction", fontsize=13)
    return fig


def plot_contour(design, ax=None):
    """Return a to-scale wall contour with key dimensions annotated.

    Mirrors the inner gas-side wall about the axis, marks segment joints,
    and overlays the 15° reference cone used to define bell length.
    Geometry only; computed in sizing.py.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from math import radians, tan
    from .constants import REFERENCE_CONE_ANGLE_DEG

    g, mm = design.geometry, 1e3
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5.5), layout="constrained")
    else:
        fig = ax.figure
    colors = ["#2367a2", "#5a8fc0", "#7a5195", "#d66820", "#c0392b", "#2a9d8f"]
    for segment, color in zip(g.segments, colors):
        for sign in (1, -1):
            ax.plot(segment.x_m * mm, sign * segment.r_m * mm, color=color, linewidth=2,
                    label=segment.name if sign == 1 else None)
    joints = [s.x_m[0] for s in g.segments[1:]]
    ax.scatter(np.array(joints) * mm,
               np.interp(joints, g.x_m, g.r_m) * mm, color="black", s=12, zorder=3)
    rt, re, rc = g.D_t_m / 2, g.D_e_m / 2, g.D_c_m / 2
    cone_x = np.array([0, g.L_cone_m])
    cone_r = rt + cone_x * tan(radians(REFERENCE_CONE_ANGLE_DEG))
    ax.plot(cone_x * mm, cone_r * mm, color="0.6", linestyle=":",
            label=f"{REFERENCE_CONE_ANGLE_DEG:g}° reference cone")
    ax.axhline(0, color="0.5", linewidth=0.8, linestyle="-.")

    def vertical(x, r, text):
        ax.annotate("", (x * mm, -r * mm), (x * mm, r * mm),
                    arrowprops={"arrowstyle": "<->", "color": "0.25"})
        ax.text(x * mm, 0, f" {text}", va="center", ha="left", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1})

    def horizontal(x0, x1, y, text):
        ax.annotate("", (x0 * mm, y * mm), (x1 * mm, y * mm),
                    arrowprops={"arrowstyle": "<->", "color": "0.25"})
        ax.text((x0 + x1) / 2 * mm, y * mm, text, va="bottom", ha="center", fontsize=9)

    x_start = g.x_m[0]
    vertical(x_start + 0.25 * g.L_cylinder_m, rc, f"Dc = {g.D_c_m*mm:.1f} mm")
    vertical(0, rt, f"Dt = {g.D_t_m*mm:.1f} mm")
    vertical(g.L_nozzle_m, re, f"De = {g.D_e_m*mm:.1f} mm")
    low = -(max(re, rc) * 1.18)
    horizontal(x_start, 0, low, f"Chamber to throat {g.L_chamber_m*mm:.1f} mm")
    horizontal(0, g.L_nozzle_m, low, f"Bell {g.L_nozzle_m*mm:.1f} mm")
    ax.set_aspect("equal")
    ax.set(xlabel="Axial position from throat [mm]", ylabel="Radius [mm]",
           ylim=(low * mm * 1.12, max(re, rc) * mm * 1.45))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left", frameon=False, ncol=2)
    r = design.requirements
    ax.set_title(f"{r.propellants} · {r.thrust_N/1e3:g} kN {r.thrust_condition} · "
                 f"pc = {r.p_c_Pa/1e6:g} MPa · O/F = {design.OF:.3f} · ε = {design.eps:.2f}\n"
                 f"L* = {g.L_star_m:.3f} m · {g.L_nozzle_m/g.L_cone_m:.0%} Rao parabolic bell "
                 f"(θn = {g.theta_n_deg:.1f}°, θe = {g.theta_e_deg:.1f}°) · inner wall, to scale",
                 fontsize=11)
    return fig


def plot_losses(design):
    """Return divergence-chart and vacuum-Isp loss-breakdown plots.

    Left: digitized SP-125 Fig.4-12 bell lambda and the conical relation
    lambda=(1+cos a)/2 at equal length, Lf=tan15°/tan a. Right: ideal
    shifting vacuum Isp reduced by lambda, eta_CF and eta_c* in turn, with
    ideal frozen Isp for comparison. Physics is in losses.py.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from math import atan, degrees, radians, tan
    from .constants import (BELL_DIVERGENCE_AREA_RATIOS, BELL_DIVERGENCE_FRACTIONS,
                            BELL_DIVERGENCE_PERCENT, G0_M_S2, REFERENCE_CONE_ANGLE_DEG)
    from .losses import conical_divergence_factor

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    ax = axes[0]
    colors = ["#9cc3e6", "#5a8fc0", "#2367a2", "#123a63"]
    for eps, row, color in zip(BELL_DIVERGENCE_AREA_RATIOS, BELL_DIVERGENCE_PERCENT, colors):
        ax.plot(np.array(BELL_DIVERGENCE_FRACTIONS) * 100, row, "o-", color=color,
                markersize=4, label=f"Bell, ε = {eps:g} (digitized)")
    fractions = np.linspace(0.55, 1.0, 100)
    cone = [100 * conical_divergence_factor(degrees(atan(tan(radians(REFERENCE_CONE_ANGLE_DEG)) / f)))
            for f in fractions]
    ax.plot(fractions * 100, cone, color="#d66820", label="Cone, (1+cos α)/2 at equal length")
    r = design.requirements
    ax.scatter([r.bell_fraction * 100], [design.divergence_factor * 100], color="black",
               zorder=3, label=f"This design: λ = {design.divergence_factor:.4f}")
    ax.set(xlabel="Length relative to 15° cone [%]", ylabel="Thrust efficiency λ [% of ideal CF]",
           title="Divergence loss · NASA SP-125 Fig. 4-12")
    ax.legend(fontsize=8, frameon=False, loc="lower right")

    ax = axes[1]
    ideal, delivered = design.ideal, design.delivered
    eps, pc = design.eps, r.p_c_Pa
    pressure = eps * ideal.p_e_Pa / pc
    momentum = ideal.CF_vac - pressure
    after_lambda = ideal.c_star_m_s * (design.divergence_factor * momentum + pressure) / G0_M_S2
    steps = [("Ideal\nshifting", ideal.Isp_vac_s),
             (f"× λ\n{design.divergence_factor:.4f}", after_lambda),
             (f"× η_CF\n{r.eta_CF:.3f}", after_lambda * r.eta_CF),
             (f"× η_c*\n{r.eta_c_star:.3f}", delivered.Isp_vac_s)]
    values = [v for _, v in steps]
    ax.bar(range(4), values, color=["#2367a2", "#5a8fc0", "#5a8fc0", "#2a9d8f"], width=0.6)
    for i, v in enumerate(values):
        ax.text(i, v + 0.4, f"{v:.1f} s", ha="center", fontsize=9)
    ax.axhline(design.ideal_frozen.Isp_vac_s, color="#d66820", linestyle="--",
               label=f"Ideal frozen: {design.ideal_frozen.Isp_vac_s:.1f} s")
    ax.set_xticks(range(4), [k for k, _ in steps])
    low = min(values + [design.ideal_frozen.Isp_vac_s])
    ax.set(ylim=(low - 12, max(values) + 4), ylabel="Vacuum specific impulse [s]",
           title=f"{r.propellants} · ε = {eps:.2f} · O/F = {design.OF:.3f}")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle("Stage 4 losses: delivered = η_c* · η_CF · (λ applied to momentum CF)\n"
                 "Efficiencies are documented assumptions, not measured values", fontsize=12)
    return fig
