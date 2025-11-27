import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- ORT/NFW Physics ---
G = 4.302e-6  # kpc km^2 / M_sun s^2

def nfw_velocity_sq(r, rho0, rs):
    """Квадрат скорости от гало NFW."""
    x = r / rs
    mass_profile = 4 * np.pi * rho0 * rs**3 * (np.log(1 + x) - x / (1 + x))
    return G * mass_profile / r

def total_velocity(r, rho0, rs, M_L_disk, v_gas, v_disk, v_bul):
    """Модель полной скорости: V^2 = V_halo^2 + M/L * V_disk^2 + V_gas^2"""
    v_halo_2 = nfw_velocity_sq(r, rho0, rs)
    # V_baryon^2 = (M/L_disk * V_disk^2) + V_gas^2 + V_bul^2
    v_baryon_2 = (M_L_disk * v_disk**2) + v_gas**2 + v_bul**2
    return np.sqrt(np.abs(v_halo_2 + v_baryon_2))

# --- Загрузка данных ---
print("Loading SPARC_multi.csv...")
df = pd.read_csv('../data/SPARC_multi.csv', comment='#')

# ВАЖНОЕ ИСПРАВЛЕНИЕ: Удаляем пробелы из названий колонок
df.columns = df.columns.str.strip()

galaxies = df['Galaxy'].unique()
print(f"Found {len(galaxies)} galaxies: {galaxies}")

# --- Настройка графиков ---
if len(galaxies) < 2:
    fig, axs = plt.subplots(1, 1, figsize=(8, 6))
    axs = [axs]
else:
    fig, axs = plt.subplots(1, len(galaxies), figsize=(6 * len(galaxies), 5))

for i, gal_name in enumerate(galaxies):
    print(f"\nAnalyzing {gal_name}...")
    data = df[df['Galaxy'] == gal_name]
    
    r = data['Radius'].values
    v_obs = data['Vobs'].values
    v_err = data['Err'].values
    v_gas = data['Vgas'].values
    v_disk = data['Vdisk'].values
    v_bul = data['Vbul'].values

    # Обертка для фиттинга
    def fit_func(r, rho0, rs, ml):
        return total_velocity(r, rho0, rs, ml, v_gas, v_disk, v_bul)

    # Начальные приближения
    p0 = [1e6, 10.0, 0.5]
    bounds = ([0, 0, 0.1], [1e9, 100, 2.0]) 

    try:
        popt, pcov = curve_fit(fit_func, r, v_obs, p0=p0, bounds=bounds, sigma=v_err)
        rho0_fit, rs_fit, ml_fit = popt
        
        # Хи-квадрат
        residuals = v_obs - fit_func(r, *popt)
        chi2 = np.sum((residuals / v_err)**2)
        dof = len(r) - len(popt)
        print(f"  FIT RESULTS: rho0={rho0_fit:.2e}, rs={rs_fit:.2f}, M/L={ml_fit:.2f}")
        print(f"  Chi2/dof: {chi2/dof:.3f}")

        # График
        ax = axs[i]
        ax.errorbar(r, v_obs, yerr=v_err, fmt='ko', label='Observed')
        ax.plot(r, fit_func(r, *popt), 'r-', label=f'Total Fit (M/L={ml_fit:.2f})')
        
        v_halo = np.sqrt(nfw_velocity_sq(r, rho0_fit, rs_fit))
        v_baryons = np.sqrt(ml_fit * v_disk**2 + v_gas**2)
        
        ax.plot(r, v_halo, 'g--', label='Dark Matter Halo')
        ax.plot(r, v_baryons, 'b:', label='Baryons')
        
        ax.set_title(f"{gal_name}")
        ax.set_xlabel("Radius (kpc)")
        ax.set_ylabel("Velocity (km/s)")
        ax.legend()
        ax.grid(True, alpha=0.3)

    except Exception as e:
        print(f"  Fit failed: {e}")

plt.tight_layout()
plt.savefig('../results/figures/sparc_multi_fit.png')
print("\nSaved plot to ../results/figures/sparc_multi_fit.png")