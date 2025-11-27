#!/usr/bin/env python3
"""
05_sparc_real_fit.py - NFW fit to real NGC 2403 data
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os

# Load data
data = np.loadtxt('../data/SPARC_NGC2403.csv', delimiter=',', comments='#')
r = data[:, 0]       # kpc
v_obs = data[:, 1]   # km/s
v_err = data[:, 2]
v_gas = data[:, 3]
v_disk = data[:, 4]
v_bul = data[:, 5]

# Baryonic velocity
v_bar = np.sqrt(v_gas**2 + v_disk**2 + v_bul**2)

print("="*60)
print("NGC 2403 - NFW FIT")
print("="*60)
print(f"Data points: {len(r)}")

def v_nfw(r, rho0, rs):
    """NFW velocity curve. rho0 in Msun/kpc^3, rs in kpc."""
    G = 4.302e-6  # kpc (km/s)^2 / Msun
    x = r / rs
    M_enc = 4 * np.pi * rho0 * rs**3 * (np.log(1 + x) - x/(1 + x))
    return np.sqrt(G * M_enc / r)

def chi2(params):
    log_rho0, rs = params
    if rs < 0.5 or rs > 50 or log_rho0 < 4 or log_rho0 > 10:
        return 1e10
    rho0 = 10**log_rho0
    v_halo = v_nfw(r, rho0, rs)
    v_tot = np.sqrt(v_bar**2 + v_halo**2)
    return np.sum(((v_obs - v_tot) / v_err)**2)

# Fit
res = minimize(chi2, [7.0, 10.0], method='Nelder-Mead')
log_rho0_fit, rs_fit = res.x
rho0_fit = 10**log_rho0_fit
chi2_min = res.fun
dof = len(r) - 2
chi2_red = chi2_min / dof

print(f"\nResults:")
print(f"  rho0 = {rho0_fit:.2e} Msun/kpc^3")
print(f"  rs = {rs_fit:.2f} kpc")
print(f"  chi2/dof = {chi2_red:.3f}")
print(f"  Fit quality: {'GOOD' if chi2_red < 2 else 'POOR'}")

# Plot
r_fine = np.linspace(0.3, 22, 200)
v_bar_fine = np.interp(r_fine, r, v_bar)
v_halo_fine = v_nfw(r_fine, rho0_fit, rs_fit)
v_tot_fine = np.sqrt(v_bar_fine**2 + v_halo_fine**2)

fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(r, v_obs, yerr=v_err, fmt='ko', capsize=3, label='Observed', ms=6)
ax.plot(r_fine, v_bar_fine, 'g--', lw=2, label='Baryonic')
ax.plot(r_fine, v_halo_fine, 'r:', lw=2, label='NFW halo')
ax.plot(r_fine, v_tot_fine, 'b-', lw=2, label='Total')

ax.set_xlabel('Radius (kpc)', fontsize=12)
ax.set_ylabel('Velocity (km/s)', fontsize=12)
ax.set_title(f'NGC 2403: NFW fit (χ²/dof = {chi2_red:.2f})', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 22)
ax.set_ylim(0, 160)

os.makedirs('../results/figures', exist_ok=True)
plt.savefig('../results/figures/ngc2403_nfw_fit.png', dpi=150)
print(f"\nSaved: ../results/figures/ngc2403_nfw_fit.png")
plt.show()