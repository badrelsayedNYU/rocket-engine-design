"""Web interface for the engine design tool (Stage 7). UI only: all physics
lives in engine_design (design, thermo, cooling, losses) and validation/."""

import json
from pathlib import Path
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from engine_design import EngineRequirements, design_engine  # noqa: E402
from engine_design.constants import (COMBUSTION_EFFICIENCY, DEFAULT_BELL_FRACTION,  # noqa: E402
                                     DEFAULT_CONTRACTION_RATIO, NOZZLE_EFFICIENCY,
                                     STANDARD_ATMOSPHERE_PA, WALL_TEMPERATURE_K)
from engine_design.propellants import PROPELLANTS  # noqa: E402
from engine_design import plots  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent / "validation"))
import run_validation as validation  # noqa: E402

st.set_page_config(page_title="Liquid Rocket Engine Design", layout="wide")
st.title("Liquid rocket engine design")
st.caption("Preliminary design of a liquid bipropellant engine from top-level requirements: "
           "NASA CEA thermochemistry (RocketCEA), isentropic nozzle theory, Rao bell contour, "
           "loss factors and Bartz heat flux, validated against the RS-25, RL10B-2 and Merlin 1D. "
           "Everything is computed by this project's own model.")
st.info("**Physics note.** Ideal performance is CEA shifting equilibrium; delivered performance "
        "applies a divergence factor λ and the assumed efficiencies η_c* and η_CF, and the throat "
        "is sized so *delivered* thrust meets the requirement. Real engines fall between frozen and "
        "shifting expansion. Bartz typically **overpredicts** throat heat flux. This is a 0D/1D "
        "preliminary-design model: no separated-flow, combustion-stability or channel design.")


@st.cache_data(show_spinner="Designing engine (CEA)...")
def run(fields):
    """Design once per input set; return the design and captured warnings."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        design = design_engine(EngineRequirements(**fields))
    messages = sorted({str(w.message) for w in caught if not issubclass(w.category, ResourceWarning)})
    return design, messages


@st.cache_data(show_spinner="Frozen chamber properties...")
def frozen_gamma(propellants, p_c_Pa, OF, eps):
    from engine_design.thermo import chamber_props
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return chamber_props(propellants, p_c_Pa, OF, eps, frozen=True).gamma


@st.cache_data(show_spinner="Running validation engines...")
def run_validation():
    engines = validation.load_engines()
    predictions = [validation.predict(e) for e in engines]
    return predictions, validation.markdown(predictions, {})


with st.sidebar:
    st.header("Requirements")
    propellants = st.selectbox("Propellants", list(PROPELLANTS), index=2)
    thrust_kN = st.number_input("Thrust (kN)", value=100.0, min_value=0.1, step=10.0)
    condition = st.radio("Thrust condition", ["vacuum", "sea_level"], horizontal=True)
    p_c_MPa = st.number_input("Chamber pressure (MPa)", value=10.0, min_value=0.2, max_value=40.0,
                              step=0.5)
    p_amb_kPa = st.number_input("Ambient pressure for sea-level / matched nozzle (kPa)",
                                value=STANDARD_ATMOSPHERE_PA / 1e3, min_value=0.1, step=1.0)

    st.header("Mixture ratio")
    auto_of = st.toggle("Optimize O/F for max Isp", value=True)
    of = None if auto_of else st.number_input("O/F (mass)", value=3.5, min_value=0.5, max_value=10.0,
                                              step=0.1)

    st.header("Nozzle")
    nozzle_mode = st.radio("Area ratio", ["Matched to ambient", "Fixed ε"], horizontal=True)
    nozzle = "sea_level_optimized" if nozzle_mode == "Matched to ambient" else \
        st.number_input("ε = A_e/A_t", value=40.0, min_value=1.5, step=5.0)
    bell = st.slider("Bell length (fraction of 15° cone)", 0.6, 1.0, DEFAULT_BELL_FRACTION, 0.05)
    st.caption("Rao angles are digitized for 5 ≤ ε ≤ 50 and λ for 10 ≤ ε ≤ 40. "
               "Outside these charts, supply wall angles and λ below; charts are never extrapolated.")
    with st.expander("Wall angles and λ (outside the charts)"):
        manual = st.checkbox("Supply θn, θe and λ", value=False)
        theta_n = st.number_input("θn (deg)", value=30.0, min_value=1.0, max_value=89.0, disabled=not manual)
        theta_e = st.number_input("θe (deg)", value=8.0, min_value=0.5, max_value=89.0, disabled=not manual)
        lam = st.number_input("λ divergence factor", value=0.99, min_value=0.9, max_value=1.0,
                              step=0.001, format="%.3f", disabled=not manual)

    st.header("Losses and cooling")
    eta_c = st.number_input("η_c* (combustion efficiency)", value=COMBUSTION_EFFICIENCY,
                            min_value=0.8, max_value=1.0, step=0.005, format="%.3f")
    eta_cf = st.number_input("η_CF (nozzle efficiency)", value=NOZZLE_EFFICIENCY,
                             min_value=0.8, max_value=1.0, step=0.005, format="%.3f")
    wall_T = st.number_input("Gas-side wall temperature (K, assumed)", value=WALL_TEMPERATURE_K,
                             min_value=300.0, max_value=1500.0, step=50.0)
    rp1_cp = None
    if propellants == "LOX/RP-1":
        if st.checkbox("Supply RP-1 c_p (no verified default)", value=False):
            rp1_cp = st.number_input("RP-1 c_p (J/(kg K))", value=2000.0, min_value=500.0, step=50.0)

    with st.expander("Chamber geometry"):
        cr = st.number_input("Contraction ratio A_c/A_t", value=DEFAULT_CONTRACTION_RATIO,
                             min_value=1.5, max_value=10.0, step=0.5)
        custom_lstar = st.checkbox("Override L*", value=False)
        l_star = st.number_input("L* (m)", value=1.0, min_value=0.2, max_value=3.0, step=0.05,
                                 disabled=not custom_lstar)

fields = dict(propellants=propellants, thrust_N=thrust_kN * 1e3, thrust_condition=condition,
              p_c_Pa=p_c_MPa * 1e6, mixture_ratio=of, nozzle=nozzle, p_amb_Pa=p_amb_kPa * 1e3,
              bell_fraction=bell, contraction_ratio=cr, eta_c_star=eta_c, eta_CF=eta_cf,
              wall_temperature_K=wall_T, coolant_cp_J_kg_K=rp1_cp,
              L_star_m=l_star if custom_lstar else None)
if manual:
    fields.update(theta_n_deg=theta_n, theta_e_deg=theta_e, divergence_factor=lam)

try:
    design, messages = run(fields)
except (ValueError, RuntimeError) as err:
    st.error(f"**Cannot design this engine:** {err}")
    st.stop()

for message in messages:
    st.warning(message)

r, g, c, cool = design.requirements, design.geometry, design.chamber, design.cooling
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Delivered Isp (vac)", f"{design.delivered.Isp_vac_s:.1f} s",
          f"{design.delivered.Isp_vac_s - design.ideal.Isp_vac_s:.1f} s vs ideal", delta_color="off")
m2.metric("O/F · ε", f"{design.OF:.3f} · {design.eps:.2f}")
m3.metric("Throat / exit D", f"{g.D_t_m*1e3:.1f} / {g.D_e_m*1e3:.0f} mm")
m4.metric("Peak heat flux", f"{cool.profile.peak_q_W_m2/1e6:.1f} MW/m²")
m5.metric("Coolant capacity / load", "n/a" if cool.capacity_ratio is None else f"{cool.capacity_ratio:.2f}")

tabs = st.tabs(["Summary", "Contour", "Heat flux", "CF & Isp", "Losses",
                "Compare to real engines", "Data"])
with tabs[0]:
    st.code(design.summary(), language=None)
with tabs[1]:
    st.pyplot(design.plot_contour())
with tabs[2]:
    st.pyplot(design.plot_heat_flux())
    st.caption("Bartz with frozen chamber properties and a uniform assumed wall temperature. The "
               "capacity ratio is an upper bound: the coolant outlet is capped at the wall "
               "temperature (methane also at its 625 K EOS limit, RP-1 at ~600 K decomposition).")
with tabs[3]:
    gamma_f = frozen_gamma(r.propellants, r.p_c_Pa, design.OF, design.eps)
    st.pyplot(plots.plot_nozzle_performance(round(gamma_f, 4), round(design.ideal.c_star_m_s, 1),
                                            r.p_c_Pa, round(design.eps, 2)))
    st.caption(f"Stage 2 constant-γ model with this design's frozen chamber γ = {gamma_f:.4f} and "
               "ideal c*. It is an independent check; the Summary tab uses CEA shifting equilibrium.")
with tabs[4]:
    if manual:
        st.write(f"User-supplied λ = {design.divergence_factor:.4f}; the SP-125 Fig. 4-12 chart "
                 "does not apply to user-supplied wall angles, so it is not shown.")
    else:
        st.pyplot(plots.plot_losses(design))
with tabs[5]:
    predictions, table = run_validation()
    st.pyplot(plots.plot_validation(predictions, validation.TARGET_PERCENT))
    st.markdown(table)
    st.caption("Published data and sources: validation/engines.yaml. The Merlin 1D is illustrative "
               "only: secondary-source data and an unsourced O/F, shown over O/F 2.3–2.4.")
with tabs[6]:
    profile = cool.profile
    df = pd.DataFrame({"x_m": profile.x_m, "r_m": profile.r_m, "mach": profile.mach,
                       "sigma": profile.sigma, "h_g_W_m2_K": profile.h_g_W_m2_K,
                       "T_aw_K": profile.T_aw_K, "q_W_m2": profile.q_W_m2})
    st.dataframe(df, width="stretch")
    st.download_button("Download contour and heat flux (CSV)", df.to_csv(index=False),
                       "contour_heat_flux.csv", "text/csv")
    summary = {
        "requirements": {k: v for k, v in fields.items()},
        "OF": design.OF, "eps": design.eps,
        "Isp_vac_s": {"frozen": design.ideal_frozen.Isp_vac_s, "shifting": design.ideal.Isp_vac_s,
                      "delivered": design.delivered.Isp_vac_s},
        "Isp_amb_s_delivered": design.delivered.Isp_amb_s,
        "c_star_m_s_delivered": design.delivered.c_star_m_s,
        "D_t_m": g.D_t_m, "D_e_m": g.D_e_m, "D_c_m": g.D_c_m,
        "mdot_kg_s": g.mdot_kg_s, "mdot_ox_kg_s": g.mdot_ox_kg_s, "mdot_fuel_kg_s": g.mdot_fuel_kg_s,
        "peak_q_W_m2": profile.peak_q_W_m2, "heat_load_W": cool.heat_load_W,
        "capacity_ratio": cool.capacity_ratio,
    }
    st.download_button("Download summary (JSON)", json.dumps(summary, indent=2, default=float),
                       "engine_summary.json", "application/json")
