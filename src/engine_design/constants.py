"""SI constants and explicitly unimplemented empirical inputs."""

# Conventional standard gravity, not local gravitational acceleration.
# NIST: https://physics.nist.gov/cgi-bin/cuu/Value?gn
G0_M_S2 = 9.80665

# Molar gas constant, rounded as specified in docs/spec/02_STAGES.md.
# NIST: https://physics.nist.gov/cgi-bin/cuu/Value?r
# Use molecular weight in kg/kmol with this value, not kg/mol.
R_U_J_KMOL_K = 8314.462618

# Empirical values remain unset until their stage and source are verified.
# Summerfield screening default; theory section 2. NASA TM-77396, eq. (7).
# https://ntrs.nasa.gov/api/citations/19840011402/downloads/19840011402.pdf
# The original 0.4 guideline is limited; this is not a precise separation model.
SEPARATION_PRESSURE_RATIO = 0.4
# NASA SP-125, 2nd ed., Table 4-1, printed p.87 (PDF p.96):
# LOX/LH2 with LIQUID injection 30–40 in; LOX/RP-1 40–50 in.
# Defaults are the interval midpoints, not measured optimum values.
# https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
CHARACTERISTIC_LENGTH_M = {"LOX/LH2": 35 * 0.0254,
                           "LOX/RP-1": 45 * 0.0254,
                           "LOX/CH4": 1.1}
# Methane 1.1 m is the midpoint of the SPEC ASSUMPTION 0.9–1.3 m;
# it is not a verified literature recommendation or stability guarantee.
CHARACTERISTIC_LENGTH_SOURCE = {
    "LOX/LH2": "NASA SP-125 Table 4-1, liquid hydrogen injection; midpoint",
    "LOX/RP-1": "NASA SP-125 Table 4-1; midpoint",
    "LOX/CH4": "Assumption: spec range 0.9–1.3 m; midpoint (not validated)",
}

# Original manual digitization for this project from NASA SP-125 Fig.4-14,
# printed p.91 (PDF p.100). Source image is retained under docs/references/.
# Read to nearest 0.5 degree, estimated reading uncertainty +/-0.5 degree;
# this is NOT an uncertainty bound on the underlying design approximation.
# Rows are (theta_n degrees, theta_e degrees) at RAO_AREA_RATIOS.
# Interpolate linearly in area ratio and bell fraction; do not extrapolate.
RAO_AREA_RATIOS = (5., 10., 15., 20., 30., 40., 50.)
RAO_ANGLE_TABLE = {
    0.6: ((28.5,32.,34.,35.5,37.5,38.5,39.), (20.,17.,15.5,14.5,13.5,13.,13.)),
    0.7: ((25.5,28.5,30.5,31.5,33.,34.,35.), (16.5,13.5,12.,11.5,10.5,10.,10.)),
    0.8: ((23.5,26.,28.,29.,30.5,31.5,32.5), (13.5,11.,9.5,9.,8.5,8.,8.)),
    0.9: ((21.5,24.,26.,27.,28.5,30.,31.), (10.5,8.5,7.5,7.,6.5,6.5,6.5)),
    1.0: ((20.5,22.5,24.5,26.,27.5,29.,29.5), (8.5,7.,5.5,5.,5.,5.,5.)),
}
# NASA SP-125 Fig.4-13: upstream and downstream throat curvature radii / Rt.
UPSTREAM_THROAT_RADIUS_RATIO = 1.5
DOWNSTREAM_THROAT_RADIUS_RATIO = 0.382
# Analyst-selected chamber fillet, contraction and convergence angle.
# Fillet closes the cylinder-to-cone joint with continuous slope.
CHAMBER_FILLET_RADIUS_RATIO = 0.5
DEFAULT_CONTRACTION_RATIO = 4.0
DEFAULT_CONVERGENCE_ANGLE_DEG = 30.0
DEFAULT_BELL_FRACTION = 0.8
REFERENCE_CONE_ANGLE_DEG = 15.0
# Stage 4 efficiencies: ASSUMPTIONS, not measured values; override per design.
# eta_c* = c*_actual/c*_ideal. Spec default 0.98 (theory section 4, typical
# 0.96–0.99); NASA SP-125 sample calc. 4-1 (PDF p.94) uses 0.975 for LOX/RP-1.
COMBUSTION_EFFICIENCY = 0.98
# eta_CF beyond divergence (boundary layer, kinetics). Midpoint of the
# spec's 0.97–0.99. Divergence is applied separately through lambda.
NOZZLE_EFFICIENCY = 0.98
# Original digitization for this project of NASA SP-125 Fig.4-12, printed
# p.91 (PDF p.100): bell-nozzle thrust efficiency lambda, percent of ideal CF.
# Curve centres located from pixels of a 600-dpi render, calibrated to the
# printed gridlines; reading uncertainty about +/-0.05 percentage points.
# Calibration check: the chart's 15° conical point reads 98.28 %, versus
# (1+cos 15°)/2 = 98.30 % (SP-125 eq. 4-8). Source image retained in
# docs/references/sp125_figure4-12.png. Curves merge above ~85% length.
# Interpolate linearly; do not extrapolate outside 10<=eps<=40, 0.6<=Lf<=1.
BELL_DIVERGENCE_AREA_RATIOS = (10., 20., 30., 40.)
BELL_DIVERGENCE_FRACTIONS = (0.6, 0.7, 0.8, 0.9, 1.0)
BELL_DIVERGENCE_PERCENT = (  # rows: area ratio; columns: bell fraction
    (96.26, 97.58, 98.55, 99.05, 99.16),
    (96.73, 98.00, 98.71, 99.05, 99.16),
    (97.17, 98.18, 98.81, 99.05, 99.16),
    (97.27, 98.24, 98.81, 99.05, 99.16),
)
# Stage 5: Bartz (1957), verify unit convention before setting the coefficient.
BARTZ_COEFFICIENT: float | None = None
WALL_TEMPERATURE_K: float | None = None

# Standard atmosphere by definition (NIST SP 811, Appendix B.8).
# https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
STANDARD_ATMOSPHERE_PA = 101325.0

# Exact dimensional conversions at the RocketCEA boundary:
# 1 P = 0.1 Pa s; 1 millipoise = 1e-4 Pa s.
# 1 thermochemical cal = 4.184 J; 1 mcal/(cm s K) = 0.4184 W/(m K).
# Reference: NIST SP 811 Appendix B.8 (poise, calorie).
MILLIPOISE_TO_PA_S = 1e-4
MCAL_CM_S_K_TO_W_M_K = 0.4184

# Nominal equilibrium sanity bands from docs/spec/01_THEORY.md section 6.
# These are heuristic warning bands, not universal physical bounds.
C_STAR_SANITY_M_S = {
    "LOX/LH2": (2300.0, 2400.0),
    "LOX/RP-1": (1750.0, 1850.0),
    "LOX/CH4": (1800.0, 1900.0),
}
GAMMA_SANITY = (1.13, 1.25)

# Analyst-selected numerical search domains, not measured engine values.
# Override with optimal_OF(..., OF_bounds=(low, high)) when needed.
DEFAULT_OF_BOUNDS = {"LOX/LH2": (2.0, 8.0), "LOX/RP-1": (1.5, 3.5),
                     "LOX/CH4": (2.0, 4.5)}
# Computational exit station; chamber/throat properties are independent of it.
CEA_DEFAULT_EPS = 40.0

# CEA transport Cp uses thermochemical cal/(g K), not IT Btu/(lbm R).
# RocketCEA std_examples.html, Transport Properties; NIST SP 811.
CAL_G_K_TO_J_KG_K = 4184.0

# U.S. Standard Atmosphere 1976, NASA/TM-2005-213659 Appendix A, Table A1.
# https://ntrs.nasa.gov/api/citations/20050207438/downloads/20050207438.pdf
# Retain the atmosphere model's historical gas constant for reproducibility.
ATMOSPHERE_R_J_KG_K = 8314.32 / 28.9644
ATMOSPHERE_T0_K = 288.15
# (upper geopotential height [m], temperature lapse rate [K/m]); 0–51 km.
ATMOSPHERE_LAYERS = ((11000.0, -0.0065), (20000.0, 0.0),
                     (32000.0, 0.001), (47000.0, 0.0028), (51000.0, 0.0))
