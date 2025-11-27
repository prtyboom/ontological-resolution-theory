#!/usr/bin/env python3
"""
08_mond_comparison.py - Compare NFW (CDM/ORT) vs MOND for NGC 3198
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os

print("="*60)
print("NFW vs MOND COMPARISON - NGC 3198")
print("="*60)

# NGC 3198 data
r = np.array([2.72, 5.44, 8.16, 10.88, 13.60, 19.04, 24.48, 29.92])
v_obs = np.array([100, 146, 150, 149, 148, 147, 149, 150])
v_err = np.array([6, 4, 3, 3, 3, 4, 5, 6])
v_gas = np.array([28, 38, 45, 51, 56, 62, 68, 72])
v_disk = np.array([95, 132, 136, 131, 124, 108, 92, 79])
v_bar = np.sqrt(v_gas**2 + v_disk**2)

# MOND parameters
a0 = 1.2e-10  # m/s^2, Milgrom's constant
G = 4.302e-6  # kpc (km/s)^2 / Msun

def v_mond(r, v_bar, a0_scale=1.0):
    """MOND interpolation (simple formula)."""
    a0_eff = a0 * a0_scale * 3.086e16  # convert to kpc/s^2
    g_bar = (v_bar * 1e3)**2 / (r * 3.086e19)  # m/s^2
    # Simple interpolation function
    nu = 1 / (1 - np.exp(-np.sqrt(g_bar / (a0 * a0_scale))))
    v_mond = v_bar * nu**0.25
    return v_mond

def v_nfw(r, rho0, rs):
    """NFW halo velocity."""
    x = r / rs
    M_enc = 4 * np.pi * rho0 * rs**3 * (np.log(1 + x) - x/(1 + x))
    return np.sqrt(G * M_enc / r)

# Fit NFW
def chi2_nfw(params):
    log_rho0, rs = params
    if rs < 1 or rs > 100 or log_rho0 < 4 or log_rho0 > 9:
        return 1e10
    v_h = v_nfw(r, 10**log_rho0, rs)
    v_tot = np.sqrt(v_bar**2 + v_h**2)
    return np.sum(((v_obs - v_tot) / v_err)**2)

res_nfw = minimize(chi2_nfw, [6.5, 15], method='Nelder-Mead')
rho0_fit = 10**res_nfw.x[0]
rs_fit = res_nfw.x[1]
chi2_nfw_val = res_nfw.fun / (len(r) - 2)

# Fit MOND (only a0 scale factor)
def chi2_mond(params):
    a0_scale = params[0]
    if a0_scale < 0.1 or a0_scale > 10:
        return 1e10
    v_m = v_mond(r, v_bar, a0_scale)
    return np.sum(((v_obs - v_m) / v_err)**2)

res_mond = minimize(chi2_mond, [1.0], method='Nelder-Mead')
a0_scale_fit = res_mond.x[0]
chi2_mond_val = res_mond.fun / (len(r) - 1)

print(f"\nNFW (CDM/ORT) Results:")
print(f"  rho0 = {rho0_fit:.2e} Msun/kpc^3")
print(f"  rs = {rs_fit:.1f} kpc")
print(f"  chi2/dof = {chi2_nfw_val:.3f}")

print(f"\nMOND Results:")
print(f"  a0 scale = {a0_scale_fit:.2f}")
print(f"  chi2/dof = {chi2_mond_val:.3f}")

print(f"\n{'='*60}")
if chi2_nfw_val < chi2_mond_val:
    print(f"WINNER: NFW (CDM/ORT) - better fit by {chi2_mond_val/chi2_nfw_val:.1f}x")
else:
    print(f"WINNER: MOND - better fit by {chi2_nfw_val/chi2_mond_val:.1f}x")
print(f"{'='*60}")

# Plot
r_fine = np.linspace(2, 32, 100)
v_bar_fine = np.interp(r_fine, r, v_bar)
v_nfw_fine = v_nfw(r_fine, rho0_fit, rs_fit)
v_tot_nfw = np.sqrt(v_bar_fine**2 + v_nfw_fine**2)
v_mond_fine = v_mond(r_fine, v_bar_fine, a0_scale_fit)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: rotation curves
ax = axes[0]
ax.errorbar(r, v_obs, yerr=v_err, fmt='ko', capsize=3, ms=8, label='Observed')
ax.plot(r_fine, v_bar_fine, 'g--', lw=2, label='Baryonic')
ax.plot(r_fine, v_tot_nfw, 'b-', lw=2.5, label=f'NFW (χ²={chi2_nfw_val:.2f})')
ax.plot(r_fine, v_mond_fine, 'r-.', lw=2.5, label=f'MOND (χ²={chi2_mond_val:.2f})')
ax.set_xlabel('Radius (kpc)', fontsize=12)
ax.set_ylabel('Velocity (km/s)', fontsize=12)
ax.set_title('NGC 3198: NFW vs MOND', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 35)
ax.set_ylim(0, 180)

# Right: residuals
ax = axes[1]
v_nfw_at_r = np.sqrt(v_bar**2 + v_nfw(r, rho0_fit, rs_fit)**2)
v_mond_at_r = v_mond(r, v_bar, a0_scale_fit)
res_nfw_pts = (v_obs - v_nfw_at_r) / v_err
res_mond_pts = (v_obs - v_mond_at_r) / v_err

ax.axhline(0, color='k', lw=1)
ax.plot(r, res_nfw_pts, 'bs-', ms=8, lw=2, label='NFW residuals')
ax.plot(r, res_mond_pts, 'r^--', ms=8, lw=2, label='MOND residuals')
ax.axhline(1, color='gray', ls=':', lw=1)
ax.axhline(-1, color='gray', ls=':', lw=1)
ax.set_xlabel('Radius (kpc)', fontsize=12)
ax.set_ylabel('Residual (σ)', fontsize=12)
ax.set_title('Fit Residuals', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 35)
ax.set_ylim(-3, 3)

plt.tight_layout()
os.makedirs('../results/figures', exist_ok=True)
plt.savefig('../results/figures/nfw_vs_mond.png', dpi=150)
print(f"\nSaved: ../results/figures/nfw_vs_mond.png")
plt.show()