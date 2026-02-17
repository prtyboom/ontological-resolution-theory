"""
ORT Simulation — Step 2: Dynamics T5.1
=======================================
S(x, i+1) = S(x, i) + U(x,i) * J * Delta_lat S(x,i)

State field: S = (q, p) per node (n=1, so 2D phase space).
Gating: U = exp(-W/Z)
Load: W = alpha_W * S^2
"""

import numpy as np
import math
from fcc_lattice import build_fcc, build_neighbours

# Canon constants
PI = math.pi
Z_BARE = PI**4 + 4*PI**2 + (PI - 3) + 1/144

# Symplectic matrix J for n=1: [[0,1],[-1,0]]
def apply_J(S):
    """S shape: (M, 2). Returns J @ S for each node."""
    out = np.empty_like(S)
    out[:, 0] = S[:, 1]
    out[:, 1] = -S[:, 0]
    return out


def discrete_laplacian(S, neighbours):
    """
    Delta_lat S(x) = (1/k) * sum_{y~x} [S(y) - S(x)]
    Continuum-normalised by C_FCC = k (in Planck units).
    So: Delta_lat S(x) = sum_{y~x} [S(y) - S(x)]
    """
    # S[neighbours] shape: (M, 12, 2)
    S_neigh = S[neighbours]  # (M, 12, 2)
    S_x = S[:, np.newaxis, :]  # (M, 1, 2)
    diff = S_neigh - S_x  # (M, 12, 2)
    return np.sum(diff, axis=1)  # (M, 2)


def compute_load(S, alpha_W=1.0):
    """W(x) = alpha_W * |S(x)|^2"""
    return alpha_W * np.sum(S**2, axis=1)  # (M,)


def gating(W, Z):
    """U(x) = exp(-W(x)/Z)"""
    return np.exp(-W / Z)


def step(S, neighbours, Z, alpha_W=1.0):
    """One tick of T5.1."""
    # Load and gating
    W = compute_load(S, alpha_W)
    U = gating(W, Z)  # (M,)

    # Laplacian
    Delta_S = discrete_laplacian(S, neighbours)  # (M, 2)

    # J * Delta_S
    J_Delta = apply_J(Delta_S)  # (M, 2)

    # Update: S(i+1) = S(i) + U * J * Delta S
    S_new = S + U[:, np.newaxis] * J_Delta

    return S_new, W, U


def run_simulation(N_lattice=4, N_ticks=200, source_mass=5.0,
                   alpha_W=1.0, save_every=10):
    """
    Run T5.1 on FCC lattice with a point source.
    Source: one node initialised with |S| = sqrt(source_mass).
    """
    # Build lattice
    coords, idx_map = build_fcc(N_lattice)
    neighbours = build_neighbours(coords, idx_map, N_lattice)
    M = len(coords)

    # Initialise state: vacuum (zero)
    S = np.zeros((M, 2), dtype=float)

    # Place point source at centre
    centre = np.array([N_lattice/2, N_lattice/2, N_lattice/2])
    dists = np.linalg.norm(coords - centre, axis=1)
    source_idx = np.argmin(dists)

    # Source: S = (sqrt(mass), 0)
    S[source_idx, 0] = math.sqrt(source_mass)

    print(f"Lattice: {M} nodes")
    print(f"Source at node {source_idx}, pos {coords[source_idx]}")
    print(f"Z_bare = {Z_BARE:.4f}")
    print(f"Source mass = {source_mass}")
    print(f"alpha_W = {alpha_W}")
    print(f"Ticks: {N_ticks}")
    print("-" * 50)

    # Storage for radial profile
    history = []

    for tick in range(N_ticks):
        S, W, U = step(S, neighbours, Z_BARE, alpha_W)

        if tick % save_every == 0:
            # Compute |S|^2 vs distance from source
            phi = np.sum(S**2, axis=1)
            snapshot = {
                'tick': tick,
                'phi': phi.copy(),
                'W': W.copy(),
                'U': U.copy(),
            }
            history.append(snapshot)

            # Quick stats
            phi_max = np.max(phi)
            phi_mean = np.mean(phi)
            U_min = np.min(U)
            print(f"  tick {tick:4d}  |S|^2_max={phi_max:.4e}"
                  f"  |S|^2_mean={phi_mean:.4e}  U_min={U_min:.4f}")

    return coords, source_idx, history, dists


if __name__ == "__main__":
    coords, src, history, dists = run_simulation(
        N_lattice=6,
        N_ticks=500,
        source_mass=0.001,
        alpha_W=0.001,
        save_every=50
    )

    # Radial profile at last snapshot
    snap = history[-1]
    phi = snap['phi']
    U = snap['U']

    # Bin by distance from source
    r = dists
    r_src = r[src]

    # Exclude source node
    mask = np.arange(len(r)) != src
    r_masked = r[mask]
    phi_masked = phi[mask]
    U_masked = U[mask]

    # Sort by distance
    order = np.argsort(r_masked)
    r_sorted = r_masked[order]
    phi_sorted = phi_masked[order]
    U_sorted = U_masked[order]

    # Print radial profile
    print("\n  RADIAL PROFILE (last tick)")
    print(f"  {'r':>8s}  {'|S|^2':>12s}  {'U':>8s}  {'r*|S|^2':>12s}")
    print("  " + "-" * 46)

    # Average in radial shells
    r_unique = np.unique(np.round(r_sorted, 4))
    for ru in r_unique[:15]:
        shell = np.abs(r_sorted - ru) < 0.01
        if np.sum(shell) == 0:
            continue
        phi_avg = np.mean(phi_sorted[shell])
        U_avg = np.mean(U_sorted[shell])
        rphi = ru * phi_avg
        print(f"  {ru:8.4f}  {phi_avg:12.4e}  {U_avg:8.4f}  {rphi:12.4e}")