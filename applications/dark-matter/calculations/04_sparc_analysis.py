#!/usr/bin/env python3
"""
04_sparc_analysis.py - Rotation curve analysis with NFW halo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import json
import os

print("="*60)
print("ROTATION CURVE ANALYSIS")
print("="*60)

# Sample galaxy data (NGC 2403-like)
r_kpc = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 20.0])
v_obs = np.array([45, 75, 100, 115, 122, 128, 131, 133, 134, 135, 136, 136, 137])
v_err = np.array([8, 6, 5, 5, 4, 4, 4, 5, 6, 7, 9, 11, 13])
v_bar = np.array([40, 70, 90, 88, 78, 68, 58, 45, 38, 32, 25, 20, 18])

print(f"\nGalaxy: Sample (NGC 2403-like)")
print(f"Data points: {len(r_kpc)}")
print(f"Radial range: {r_kpc.min():.1f} - {r_kpc.max():.1f} kpc")

def v_nfw(r, v200, c):
    """NFW halo circular velocity."""
    r200 = v200 * 1.4  # rough approximation: r200 ~ v200 in kpc for typical H0
    x = r / r200 * c
    f_c = np.log(1 + c) - c / (1 + c)
    f_x = np.log(1 + x) - x / (1 + x)
    v2 = v200**2 * (c / x) * f_x / f_c
    return np.sqrt(np.maximum(v2, 0))

def v_total(r, v200, c, v_bar_interp):
    """Total velocity: sqrt(v_bar^2 + v_halo^2)"""
    v_h = v_nfw(r, v200, c)
    return np.sqrt(v_bar_interp**2 + v_h**2)

def chi2(params):
    """Chi-squared for fitting."""
    v200, c = params
    if v200 < 50 or v200 > 300 or c < 2 or c > 40:
        return 1e10
    v_model = v_total(r_kpc, v200, c, v_bar)
    return np.sum(((v_obs - v_model) / v_err)**2)

# Fit
result = minimize(chi2, [130, 10], method='Nelder-Mead')
v200_fit, c_fit = result.x
chi2_min = result.fun
dof = len(r_kpc) - 2
chi2_red = chi2_min / dof

print(f"\nNFW Fit Results:")
print(f"  V200 = {v200_fit:.1f} km/s")
print(f"  c = {c_fit:.1f}")
print(f"  chi2/dof = {chi2_red:.2f}")

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

r_fine = np.linspace(0.3, 25, 200)
v_bar_fine = np.interp(r_fine, r_kpc, v_bar)
v_halo_fine = v_nfw(r_fine, v200_fit, c_fit)
v_tot_fine = np.sqrt(v_bar_fine**2 + v_halo_fine**2)

ax.errorbar(r_kpc, v_obs, yerr=v_err, fmt='ko', capsize=3, label='Observed')
ax.plot(r_fine, v_bar_fine, 'g--', linewidth=2, label='Baryonic')
ax.plot(r_fine, v_halo_fine, 'r:', linewidth=2, label='NFW Halo')
ax.plot(r_fine, v_tot_fine, 'b-', linewidth=2, label='Total')

ax.set_xlabel('Radius (kpc)', fontsize=12)
ax.set_ylabel('Circular Velocity (km/s)', fontsize=12)
ax.set_title(f'Rotation Curve Fit (chi2/dof = {chi2_red:.2f})', fontsize=14)
ax.legend(loc='lower right')
ax.set_xlim(0, 25)
ax.set_ylim(0, 180)
ax.grid(True, alpha=0.3)

os.makedirs('../results/figures', exist_ok=True)
plt.savefig('../results/figures/rotation_curve.png', dpi=150)
plt.savefig('../results/figures/rotation_curve.pdf', dpi=300)
print(f"\nSaved: ../results/figures/rotation_curve.png")

plt.show()

# DM fraction at r=10 kpc
v_bar_10 = np.interp(10, r_kpc, v_bar)
v_halo_10 = v_nfw(10, v200_fit, c_fit)
f_DM = v_halo_10**2 / (v_bar_10**2 + v_halo_10**2)
print(f"\nDM fraction at 10 kpc: {f_DM:.1%}")

# Save
results = {
    'v200_km_s': float(v200_fit),
    'concentration': float(c_fit),
    'chi2_dof': float(chi2_red),
    'DM_fraction_10kpc': float(f_DM),
}
with open('../results/04_sparc_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Saved: ../results/04_sparc_analysis.json")