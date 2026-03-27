import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# fixed simulation parameters
k0 = 0.001
tau = 10.0
alpha = 100.0
n_genes = 10
dt = 0.01
nstep = 3000000
nequilrun = 1000000

# measurement time after equilibration
T = (nstep - nequilrun) * dt

# baseline transcription rate at Jbar = 0
kt0 = k0 / (1.0 + k0 * tau)

# load one data file
# col0 = Jbar/D, col1 = total number of events
def load_rate_file(filename):
    data = np.genfromtxt(filename, comments="!") # ignore comments
    data = data[~np.isnan(data).any(axis=1)]

    JD = data[:, 0]
    events = data[:, 1]

    # convert total events into transcription rate per gene
    kt = events / (T * n_genes)

    # normalise by baseline rate
    kt_scaled = kt / kt0

    return JD, kt_scaled

# load all 4 data sets
# JD1 is the seed = 62 case used for fitting/plotting
files = [
    "transcription_activity_vs_JD1.dat",
    "transcription_activity_vs_JD2.dat",
    "transcription_activity_vs_JD3.dat",
    "transcription_activity_vs_JD4.dat",
]

JD_list = []
kt_list = []

for f in files:
    JD, kt_scaled = load_rate_file(f)
    JD_list.append(JD)
    kt_list.append(kt_scaled)

JD_array = np.array(JD_list)
kt_array = np.array(kt_list)

JD = JD_array[0]

# seed = 62 data (file 1)
kt_seed62 = kt_array[0]

# mean and standard deviation across the 4 simulations
kt_mean4 = np.mean(kt_array, axis=0)
sigma_y = np.std(kt_array, axis=0, ddof=1)

# avoid division by zero in chi-squared if any point has zero spread
sigma_floor = 1e-8
sigma_eff = np.maximum(sigma_y, sigma_floor)

# 1) General analytical mean-field model
# Phi = k_on * tau
# h = k0*tau*(1 + alpha*(Jbar/D)/2) - 1
# kt = k_on / (1 + k_on*tau) = (Phi/tau)/(1+Phi)
def analytical_meanfield_rate_scaled(jd_vals):
    jd_vals = np.asarray(jd_vals)
    h = k0 * tau * (1.0 + alpha * jd_vals / 2.0) - 1.0
    Phi = (h + np.sqrt(h**2 + 4.0 * k0 * tau)) / 2.0
    kt = (Phi / tau) / (1.0 + Phi)
    return kt / kt0

# 2) Semi-analytical model
# sigma_p = (Jbar/D) * A * kon / (B*kon + 1)
# kon = k0 * (1 + alpha * sigma_p)
# kt  = kon / (1 + kon*tau)
def solve_kon(jd, A, B, k0=0.001, alpha=100.0, tau=10.0, n_iter=5000):
    kon = k0
    for _ in range(n_iter):
        sigma_p = jd * (A * kon) / (B * kon + 1.0)
        kon_new = k0 * (1.0 + alpha * sigma_p)
        if abs(kon_new - kon) < 1e-14:
            break
        kon = kon_new
    return kon

def semianalytical_rate_scaled(jd_vals, A, B):
    out = []
    for jd in jd_vals:
        kon = solve_kon(jd, A, B, k0=k0, alpha=alpha, tau=tau)
        kt = kon / (1.0 + kon * tau)
        out.append(kt / kt0)
    return np.array(out)

# only fit the seed = 62 data for semi-analytical model
def fit_residuals(params, jd_vals, ydata):
    logA, logB = params
    A = np.exp(logA)
    B = np.exp(logB)
    ymodel = semianalytical_rate_scaled(jd_vals, A, B)
    return ymodel - ydata

p0 = np.log([4.0, 10.0])

fit = least_squares(fit_residuals, p0, args=(JD, kt_seed62))
A_fit, B_fit = np.exp(fit.x)

print(f"Best-fit A = {A_fit:.6f}")
print(f"Best-fit B = {B_fit:.6f}")

# evaluate both models
JD_fine = np.linspace(0, JD.max(), 400)

kt_analytical_fine = analytical_meanfield_rate_scaled(JD_fine)
kt_semi_fine = semianalytical_rate_scaled(JD_fine, A_fit, B_fit)

kt_analytical_at_data = analytical_meanfield_rate_scaled(JD)
kt_semi_at_data = semianalytical_rate_scaled(JD, A_fit, B_fit)

# chi-squared and reduced chi-squared
# using sigma from the 4-run standard deviation
N = len(JD)

# analytical mean-field has no fitted parameters here
chi2_analytical = np.sum(((kt_seed62 - kt_analytical_at_data) / sigma_eff) ** 2)
red_chi2_analytical = chi2_analytical / N

# semi-analytical has 2 fitted parameters: A, B
chi2_semi = np.sum(((kt_seed62 - kt_semi_at_data) / sigma_eff) ** 2)
red_chi2_semi = chi2_semi / (N - 2)

print(f"Analytical mean-field chi-squared = {chi2_analytical:.6f}")
print(f"Analytical mean-field reduced chi-squared = {red_chi2_analytical:.6f}")

print(f"Semianalytical chi-squared = {chi2_semi:.6f}")
print(f"Semianalytical reduced chi-squared = {red_chi2_semi:.6f}")

# residuals
# model residuals: seed=62 data minus model prediction
resid_analytical = kt_seed62 - kt_analytical_at_data
resid_semi = kt_seed62 - kt_semi_at_data

# empirical residual: seed=62 data minus mean over 4 runs
resid_seed_vs_mean4 = kt_seed62 - kt_mean4

# standardised residuals
std_resid_analytical = resid_analytical / sigma_eff
std_resid_semi = resid_semi / sigma_eff
std_resid_seed_vs_mean4 = resid_seed_vs_mean4 / sigma_eff

# plotting settings
sim_color = "#A020F0"   # vivid purple
analytical_color = "#1f77b4"
semi_color = "#2ca02c"

# Plot 1: analytical mean-field vs simulation
plt.figure(figsize=(8, 5))

plt.plot(
    JD_fine,
    kt_analytical_fine,
    lw=2,
    color=analytical_color,
    label="Analytical mean-field model",
    zorder=1
)

plt.errorbar(
    JD,
    kt_seed62,
    yerr=sigma_y,
    fmt='o',
    color=sim_color,
    ecolor=sim_color,
    elinewidth=1.8,
    capsize=4,
    capthick=1.4,
    markersize=6,
    markerfacecolor=sim_color,
    markeredgecolor=sim_color,
    label="Simulation (seed = 62)",
    zorder=3
)

plt.title("Analytical Mean-Field Fit")
plt.xlabel(r"$\bar{J}/D$")
plt.ylabel("Normalised transcription rate per gene")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Plot 2: semi-analytical mean-field fit vs simulation
plt.figure(figsize=(8, 5))

plt.plot(
    JD_fine,
    kt_semi_fine,
    lw=2,
    color=semi_color,
    label="Fitted semianalytical model",
    zorder=1
)

plt.errorbar(
    JD,
    kt_seed62,
    yerr=sigma_y,
    fmt='o',
    color=sim_color,
    ecolor=sim_color,
    elinewidth=1.8,
    capsize=4,
    capthick=1.4,
    markersize=6,
    markerfacecolor=sim_color,
    markeredgecolor=sim_color,
    label="Simulation (seed = 62)",
    zorder=3
)

plt.title("Semianalytical Fit")
plt.xlabel(r"$\bar{J}/D$")
plt.ylabel("Normalised transcription rate per gene")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Plot 3: raw residuals
plt.figure(figsize=(8, 5))

plt.axhline(0.0, color='black', lw=1, linestyle='--')
plt.plot(JD, resid_analytical, 's-', label="Residual: seed62 - analytical mean field")
plt.plot(JD, resid_semi, 'o-', label="Residual: seed62 - semianalytical")
plt.plot(JD, resid_seed_vs_mean4, 'd-', label="Residual: seed62 - mean(4 runs)")

plt.title("Raw Residuals")
plt.xlabel(r"$\bar{J}/D$")
plt.ylabel("Residual")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Plot 4: standardised residuals
plt.figure(figsize=(8, 5))

plt.axhline(0.0, color='black', lw=1, linestyle='--')
plt.axhline(1.0, color='grey', lw=1, linestyle=':')
plt.axhline(-1.0, color='grey', lw=1, linestyle=':')
plt.axhline(2.0, color='grey', lw=1, linestyle=':')
plt.axhline(-2.0, color='grey', lw=1, linestyle=':')

plt.plot(JD, std_resid_analytical, 's-', label="Std. residual: analytical mean field")
plt.plot(JD, std_resid_semi, 'o-', label="Std. residual: semianalytical")
plt.plot(JD, std_resid_seed_vs_mean4, 'd-', label="Std. residual: seed62 - mean(4 runs)")

plt.title("Standardised Residuals")
plt.xlabel(r"$\bar{J}/D$")
plt.ylabel("Standardised residual")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
