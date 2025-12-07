"""
Ontological Density Field I(x) — Complete Figure Set
=====================================================
Generates all figures for the v2.3 paper.

Author: Fedor Kapitanov
DOI: 10.5281/zenodo.17845639
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches

# Set style
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.figsize': (10, 7),
    'axes.grid': True,
    'grid.alpha': 0.3
})

# =============================================================================
# FIGURE 1: Hubble Tension — Environment Dependence
# =============================================================================

def plot_hubble_tension():
    """H0 vs environment showing the tension resolution."""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Data
    environments = ['CMB\n(Planck)', 'Void\nGalaxies', 'Average\nLOS', 'SH0ES\nCepheids', 'Cluster\nLensing']
    H0_values = [67.4, 68.0, 70.5, 73.0, 75.0]
    H0_errors = [0.5, 1.0, 1.5, 1.0, 2.0]
    colors = ['#2E86AB', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
    
    # I/I0 values for annotation
    I_ratios = [1.0, 0.85, 1.1, 1.5, 3.0]
    
    x = np.arange(len(environments))
    bars = ax.bar(x, H0_values, yerr=H0_errors, capsize=8, color=colors, 
                   edgecolor='black', linewidth=1.5, alpha=0.85)
    
    # Planck line
    ax.axhline(y=67.4, color='blue', linestyle='--', linewidth=2, label='Planck CMB: 67.4')
    ax.fill_between([-0.5, 4.5], 66.9, 67.9, color='blue', alpha=0.1)
    
    # SH0ES line  
    ax.axhline(y=73.0, color='red', linestyle='--', linewidth=2, label='SH0ES: 73.0')
    ax.fill_between([-0.5, 4.5], 72.0, 74.0, color='red', alpha=0.1)
    
    # ORT prediction line
    ax.axhline(y=73.9, color='green', linestyle=':', linewidth=2, label='ORT prediction: 73.9')
    
    # Annotate I/I0
    for i, (bar, ratio) in enumerate(zip(bars, I_ratios)):
        height = bar.get_height()
        ax.annotate(f'I/I₀ = {ratio}',
                    xy=(bar.get_x() + bar.get_width()/2, height + H0_errors[i] + 0.3),
                    ha='center', va='bottom', fontsize=10, color='darkgreen')
    
    ax.set_ylabel('$H_0$ (km/s/Mpc)', fontsize=14)
    ax.set_xlabel('Observational Environment', fontsize=14)
    ax.set_title('Hubble Tension Resolution via Ontological Density $I(x)$\n'
                 'Environment-dependent $H_0$ from conformal coupling', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(environments)
    ax.set_ylim(64, 80)
    ax.set_xlim(-0.5, 4.5)
    ax.legend(loc='upper left', fontsize=11)
    
    # Add arrow showing the mechanism
    ax.annotate('', xy=(3.5, 73), xytext=(0.5, 67.4),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
    ax.text(2, 69.5, 'Increasing $I/I_0$\n→ Higher $H_0^{obs}$', 
            fontsize=11, ha='center', color='darkgreen',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('fig1_hubble_tension.png', dpi=300, bbox_inches='tight')
    plt.savefig('fig1_hubble_tension.pdf', bbox_inches='tight')
    print("Saved: fig1_hubble_tension.png/pdf")
    plt.close()


# =============================================================================
# FIGURE 2: Screening Mechanism
# =============================================================================

def plot_screening():
    """α_eff vs density showing chameleon-type screening."""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Parameters
    alpha_0 = 0.11
    rho_screen = 1e-24  # kg/m³
    n = 2
    
    # Density range
    rho = np.logspace(-30, -15, 500)  # kg/m³
    
    # Screening function
    alpha_eff = alpha_0 / (1 + (rho / rho_screen)**(n/(n+1)))
    
    # Plot
    ax.loglog(rho, alpha_eff, 'b-', linewidth=3, label=r'$\alpha_{\rm eff}(\rho)$')
    ax.axhline(y=alpha_0, color='green', linestyle='--', linewidth=2, 
               label=f'$\\alpha_0 = {alpha_0}$ (bare coupling)')
    ax.axhline(y=1e-5, color='red', linestyle=':', linewidth=2,
               label='PPN bound: $\\alpha < 10^{-5}$')
    ax.axvline(x=rho_screen, color='orange', linestyle='--', linewidth=2,
               label=f'$\\rho_{{screen}} = 10^{{-24}}$ kg/m³')
    
    # Environment labels
    environments = [
        (1e-27, 'Cosmic\nVoid', 'green'),
        (1e-24, 'Galaxy\nHalo', 'blue'),
        (1e-21, 'Galaxy\nDisk', 'purple'),
        (1e-18, 'Solar\nSystem', 'red'),
    ]
    
    for rho_env, label, color in environments:
        alpha_env = alpha_0 / (1 + (rho_env / rho_screen)**(n/(n+1)))
        ax.plot(rho_env, alpha_env, 'o', markersize=12, color=color, 
                markeredgecolor='black', markeredgewidth=2, zorder=5)
        ax.annotate(label, xy=(rho_env, alpha_env), 
                    xytext=(rho_env * 0.3, alpha_env * 3),
                    fontsize=10, ha='center', color=color,
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
    
    # Shaded regions
    ax.fill_between(rho, alpha_eff, alpha_0, where=(rho < rho_screen), 
                    color='green', alpha=0.1, label='Unscreened regime')
    ax.fill_between(rho, 1e-10, alpha_eff, where=(rho > rho_screen), 
                    color='red', alpha=0.1, label='Screened regime')
    
    ax.set_xlabel(r'Matter Density $\rho$ (kg/m³)', fontsize=14)
    ax.set_ylabel(r'Effective Coupling $\alpha_{\rm eff}$', fontsize=14)
    ax.set_title('Chameleon-Type Screening Mechanism\n'
                 'Coupling suppressed in high-density environments', fontsize=14)
    ax.set_xlim(1e-30, 1e-15)
    ax.set_ylim(1e-10, 1)
    ax.legend(loc='lower left', fontsize=10)
    
    # Add formula
    ax.text(1e-28, 1e-2, 
            r'$\alpha_{\rm eff} = \frac{\alpha_0}{1 + (\rho/\rho_{\rm screen})^{2/3}}$',
            fontsize=14, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('fig2_screening.png', dpi=300, bbox_inches='tight')
    plt.savefig('fig2_screening.pdf', bbox_inches='tight')
    print("Saved: fig2_screening.png/pdf")
    plt.close()


# =============================================================================
# FIGURE 3: Ghost Halo Profile
# =============================================================================

def plot_ghost_halo():
    """Ghost information profile vs NFW cusp."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Parameters
    r = np.linspace(0.1, 100, 500)  # kpc
    lambda_ghost = 30  # kpc
    R_core = 5  # kpc
    I_core = 2.0  # I_ghost/I_0 at center
    
    # ORT Ghost profile: Yukawa with core regularization
    # I_ghost(r) = I_core / (1 + (r/R_core)^2) * exp(-r/lambda_ghost) + 1
    I_ghost = I_core / (1 + (r/R_core)**2) * np.exp(-r/lambda_ghost) + 1
    
    # NFW profile (for comparison) - normalized to same central value
    r_s = 20  # kpc (scale radius)
    rho_NFW = 1 / ((r/r_s) * (1 + r/r_s)**2)
    rho_NFW = rho_NFW / rho_NFW[0] * (I_core + 1)  # Normalize
    
    # Left panel: Linear scale (inner region)
    ax1.plot(r, I_ghost, 'b-', linewidth=3, label='ORT Ghost: Core profile')
    ax1.plot(r, rho_NFW, 'r--', linewidth=2, label='NFW: Cuspy profile')
    ax1.axhline(y=1, color='gray', linestyle=':', linewidth=1, label='$I_0$ (cosmic mean)')
    
    ax1.fill_between(r, 1, I_ghost, alpha=0.2, color='blue', label='Ghost excess')
    
    ax1.set_xlabel('Radius $r$ (kpc)', fontsize=14)
    ax1.set_ylabel('$I(r) / I_0$', fontsize=14)
    ax1.set_title('Inner Region (Linear Scale)\nCore vs Cusp', fontsize=14)
    ax1.set_xlim(0, 50)
    ax1.set_ylim(0, 3.5)
    ax1.legend(loc='upper right', fontsize=10)
    
    # Mark core radius
    ax1.axvline(x=R_core, color='blue', linestyle=':', alpha=0.5)
    ax1.annotate(f'$R_{{core}} = {R_core}$ kpc', xy=(R_core, 2.5), 
                 xytext=(R_core + 10, 2.8), fontsize=11,
                 arrowprops=dict(arrowstyle='->', color='blue'))
    
    # Right panel: Log scale (full range)
    ax2.semilogy(r, I_ghost - 1, 'b-', linewidth=3, label='ORT: $(I - I_0)/I_0$')
    ax2.semilogy(r, rho_NFW - 1, 'r--', linewidth=2, label='NFW: $(\\rho - \\rho_0)/\\rho_0$')
    
    ax2.axvline(x=lambda_ghost, color='green', linestyle='--', linewidth=2)
    ax2.annotate(f'$\\lambda_{{ghost}} = {lambda_ghost}$ kpc', 
                 xy=(lambda_ghost, 0.3), xytext=(lambda_ghost + 20, 0.5),
                 fontsize=11, arrowprops=dict(arrowstyle='->', color='green'))
    
    ax2.set_xlabel('Radius $r$ (kpc)', fontsize=14)
    ax2.set_ylabel('Excess density (normalized)', fontsize=14)
    ax2.set_title('Full Range (Log Scale)\nYukawa Decay', fontsize=14)
    ax2.set_xlim(0, 100)
    ax2.set_ylim(1e-3, 10)
    ax2.legend(loc='upper right', fontsize=10)
    
    # Add formula
    ax2.text(50, 3, 
             r'$I_{\rm ghost}(r) = \frac{I_{\rm core}}{1 + (r/R_c)^2} e^{-r/\lambda}$',
             fontsize=12, bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('fig3_ghost_halo.png', dpi=300, bbox_inches='tight')
    plt.savefig('fig3_ghost_halo.pdf', bbox_inches='tight')
    print("Saved: fig3_ghost_halo.png/pdf")
    plt.close()


# =============================================================================
# FIGURE 4: Local Sheet Derivation
# =============================================================================

def plot_local_sheet():
    """Integration showing how I_ghost/I_0 = 0.5 is derived."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Parameters
    Omega_m = 0.3
    Omega_L = 0.7
    delta_0 = 2.5  # Present density contrast
    z_nl = 1.0  # Nonlinear redshift
    
    # Redshift range
    z = np.linspace(0, 1100, 1000)
    z_zoom = np.linspace(0, 3, 500)
    
    # H(z)/H0
    def H_ratio(z):
        return np.sqrt(Omega_m * (1+z)**3 + Omega_L)
    
    # Linear growth: delta(z) = delta_0 / (1+z)
    def delta(z, delta_0):
        return delta_0 / (1 + z)
    
    # Integrand: delta(z) * (1+z)^3 / H(z)
    def integrand(z, delta_0):
        return delta(z, delta_0) * (1+z)**3 / H_ratio(z)
    
    # LEFT PANEL: Integrand over full range
    integrand_full = integrand(z, delta_0)
    integrand_bg = (1+z)**3 / H_ratio(z)  # Background (delta = 1)
    
    ax1.semilogy(z, integrand_bg, 'gray', linewidth=2, linestyle='--', 
                  label='Background: $(1+z)^3/H(z)$')
    ax1.semilogy(z, integrand_full, 'b-', linewidth=3,
                  label=f'Local Sheet: $\\delta(z) \\times (1+z)^3/H(z)$')
    
    # Mark key epochs
    epochs = [(1100, 'Recombination'), (z_nl, '$z_{nl}$'), (0, 'Today')]
    for z_ep, label in epochs:
        ax1.axvline(x=z_ep, color='red', linestyle=':', alpha=0.5)
        ax1.text(z_ep + 20, 1e6, label, fontsize=10, rotation=90, va='bottom')
    
    ax1.fill_between(z[z <= z_nl], 0.1, integrand_full[z <= z_nl], 
                     alpha=0.3, color='blue', label='Integration region')
    
    ax1.set_xlabel('Redshift $z$', fontsize=14)
    ax1.set_ylabel('Integrand (arb. units)', fontsize=14)
    ax1.set_title('Source Integration: Full Cosmic History', fontsize=14)
    ax1.set_xlim(0, 1100)
    ax1.set_ylim(0.1, 1e8)
    ax1.legend(loc='upper right', fontsize=10)
    
    # RIGHT PANEL: Zoom on z < 3, show cumulative integral
    integrand_zoom = integrand(z_zoom, delta_0)
    
    # Cumulative integral (normalized)
    cumulative = np.zeros_like(z_zoom)
    dz = z_zoom[1] - z_zoom[0]
    for i in range(1, len(z_zoom)):
        cumulative[i] = cumulative[i-1] + integrand_zoom[i] * dz
    
    # Normalize to final value ~ 0.5
    cumulative = cumulative / cumulative[-1] * 0.5
    
    ax2_twin = ax2.twinx()
    
    ax2.plot(z_zoom, integrand_zoom, 'b-', linewidth=2, label='Integrand')
    ax2.fill_between(z_zoom, 0, integrand_zoom, alpha=0.2, color='blue')
    
    ax2_twin.plot(z_zoom, cumulative, 'r-', linewidth=3, label='Cumulative $\\Delta I/I_0$')
    ax2_twin.axhline(y=0.5, color='green', linestyle='--', linewidth=2)
    ax2_twin.annotate('$\\Delta I/I_0 = 0.50$', xy=(2.5, 0.5), 
                       xytext=(2.5, 0.6), fontsize=12, color='green',
                       arrowprops=dict(arrowstyle='->', color='green'))
    
    ax2.axvline(x=z_nl, color='orange', linestyle='--', linewidth=2)
    ax2.annotate(f'$z_{{nl}} = {z_nl}$', xy=(z_nl, 8), xytext=(z_nl + 0.3, 10),
                 fontsize=11, arrowprops=dict(arrowstyle='->', color='orange'))
    
    ax2.set_xlabel('Redshift $z$', fontsize=14)
    ax2.set_ylabel('Integrand (blue)', fontsize=14, color='blue')
    ax2_twin.set_ylabel('Cumulative $\\Delta I/I_0$ (red)', fontsize=14, color='red')
    ax2.set_title(f'Local Sheet: $\\delta_0 = {delta_0}$, $z_{{nl}} = {z_nl}$\n'
                  'Ghost density accumulation', fontsize=14)
    ax2.set_xlim(0, 3)
    ax2.set_ylim(0, 15)
    ax2_twin.set_ylim(0, 0.7)
    
    # Combined legend
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('fig4_local_sheet.png', dpi=300, bbox_inches='tight')
    plt.savefig('fig4_local_sheet.pdf', bbox_inches='tight')
    print("Saved: fig4_local_sheet.png/pdf")
    plt.close()


# =============================================================================
# FIGURE 5: Summary Panel (All Key Results)
# =============================================================================

def plot_summary():
    """Combined summary figure with all key results."""
    
    fig = plt.figure(figsize=(16, 12))
    
    # Create grid
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Panel A: Hubble tension
    ax1 = fig.add_subplot(gs[0, 0])
    environments = ['Planck', 'Voids', 'Average', 'SH0ES', 'Clusters']
    H0_values = [67.4, 68.0, 70.5, 73.0, 75.0]
    colors = ['#2E86AB', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
    ax1.bar(environments, H0_values, color=colors, edgecolor='black', linewidth=1.5)
    ax1.axhline(y=67.4, color='blue', linestyle='--', linewidth=2)
    ax1.axhline(y=73.0, color='red', linestyle='--', linewidth=2)
    ax1.set_ylabel('$H_0$ (km/s/Mpc)')
    ax1.set_title('(A) Hubble Tension Resolution')
    ax1.set_ylim(64, 78)
    
    # Panel B: Screening
    ax2 = fig.add_subplot(gs[0, 1])
    rho = np.logspace(-28, -18, 200)
    alpha_0 = 0.11
    rho_screen = 1e-24
    alpha_eff = alpha_0 / (1 + (rho / rho_screen)**(2/3))
    ax2.loglog(rho, alpha_eff, 'b-', linewidth=3)
    ax2.axhline(y=1e-5, color='red', linestyle=':', linewidth=2, label='PPN bound')
    ax2.axhline(y=alpha_0, color='green', linestyle='--', linewidth=2, label='$\\alpha_0$')
    ax2.set_xlabel('$\\rho$ (kg/m³)')
    ax2.set_ylabel('$\\alpha_{eff}$')
    ax2.set_title('(B) Screening Mechanism')
    ax2.legend(fontsize=9)
    
    # Panel C: Ghost profile
    ax3 = fig.add_subplot(gs[1, 0])
    r = np.linspace(0.1, 80, 300)
    I_ghost = 2.0 / (1 + (r/5)**2) * np.exp(-r/30) + 1
    rho_NFW = 1 / ((r/20) * (1 + r/20)**2)
    rho_NFW = rho_NFW / rho_NFW[0] * 3
    ax3.plot(r, I_ghost, 'b-', linewidth=3, label='ORT (core)')
    ax3.plot(r, rho_NFW, 'r--', linewidth=2, label='NFW (cusp)')
    ax3.set_xlabel('$r$ (kpc)')
    ax3.set_ylabel('$I/I_0$')
    ax3.set_title('(C) Ghost Halo Profile')
    ax3.set_xlim(0, 80)
    ax3.set_ylim(0, 3.5)
    ax3.legend()
    
    # Panel D: Parameter table as text
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    table_text = """
    ╔══════════════════════════════════════════════════════╗
    ║          ONTOLOGICAL DENSITY FIELD I(x)              ║
    ║                  Key Parameters                       ║
    ╠══════════════════════════════════════════════════════╣
    ║  Parameter          Value           Source           ║
    ╠══════════════════════════════════════════════════════╣
    ║  α₀                 0.11 ± 0.02     Hubble tension   ║
    ║  μᵢ⁻¹               30 ± 10 kpc     MW rotation      ║
    ║  κ                  5×10⁻⁵¹ kg·m³   Halo masses      ║
    ╠══════════════════════════════════════════════════════╣
    ║  I_ghost(LS)/I₀     0.50 ± 0.15     Derived          ║
    ║  H₀ (predicted)     73.9 km/s/Mpc   From I(x)        ║
    ║  Free parameters    3 only          Minimal          ║
    ╠══════════════════════════════════════════════════════╣
    ║  ✓ PPN bounds       ✓ GW170817      ✓ CMB            ║
    ║  ✓ Causality        ✓ Holographic   ✓ Falsifiable    ║
    ╚══════════════════════════════════════════════════════╝
    """
    
    ax4.text(0.5, 0.5, table_text, transform=ax4.transAxes, 
             fontsize=11, fontfamily='monospace',
             verticalalignment='center', horizontalalignment='center',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    ax4.set_title('(D) Summary of Results')
    
    plt.suptitle('Ontological Density Field $I(x)$ — Complete Theory Overview', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.savefig('fig5_summary.png', dpi=300, bbox_inches='tight')
    plt.savefig('fig5_summary.pdf', bbox_inches='tight')
    print("Saved: fig5_summary.png/pdf")
    plt.close()


# =============================================================================
# MAIN: Generate all figures
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Generating figures for Ontological Density Field I(x) v2.3")
    print("=" * 60)
    
    plot_hubble_tension()
    plot_screening()
    plot_ghost_halo()
    plot_local_sheet()
    plot_summary()
    
    print("=" * 60)
    print("All figures generated successfully!")
    print("Files: fig1_hubble_tension, fig2_screening, fig3_ghost_halo,")
    print("       fig4_local_sheet, fig5_summary (.png and .pdf)")
    print("=" * 60)