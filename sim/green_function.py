"""
ORT Simulation — Lattice Green Function (T5.4)
================================================
Sparse solver. Subtract mean to remove zero mode.
Test: Phi(r) - <Phi> ~ 1/(4*pi*r)
"""

import numpy as np
import math
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from fcc_lattice import build_fcc, build_neighbours


def solve_green(N_lattice=16):

    coords, idx_map = build_fcc(N_lattice)
    neighbours = build_neighbours(coords, idx_map, N_lattice)
    M = len(coords)

    print(f"FCC {N_lattice}^3: {M} nodes")

    # Source at centre
    centre = np.array([N_lattice/2, N_lattice/2, N_lattice/2])
    dists = np.linalg.norm(coords - centre, axis=1)
    src = np.argmin(dists)

    # Farthest node for pinning
    pin = np.argmax(dists)

    print(f"Source node {src}, pos {coords[src]}")

    # Sparse Laplacian
    L = lil_matrix((M, M), dtype=float)
    for i in range(M):
        nn = neighbours[i]
        valid = nn[nn >= 0]
        L[i, i] = -len(valid)
        for j in valid:
            L[i, j] += 1.0

    # RHS: delta source
    rhs = np.zeros(M)
    rhs[src] = -1.0

    # Pin one node to break zero mode
    L[pin, :] = 0
    L[pin, pin] = 1.0
    rhs[pin] = 0.0

    print("Solving...")
    phi_raw = spsolve(L.tocsr(), rhs)

    # Subtract mean — this removes the constant background
    phi = phi_raw - np.mean(phi_raw)

    # Radial analysis
    r = dists
    pi = math.pi
    target = 1.0 / (4 * pi)

    # Exclude source
    mask = np.arange(M) != src
    r_m = r[mask]
    phi_m = phi[mask]

    # Shell averages
    r_bins = np.unique(np.round(r_m, 3))

    print(f"\n  Target: r*Phi = 1/(4*pi) = {target:.6f}")
    print(f"\n  {'r':>8s}  {'Phi':>12s}  {'1/4pi*r':>12s}"
          f"  {'r*Phi':>10s}  {'r*Phi/tgt':>10s}")
    print("  " + "-" * 60)

    data = []
    for ru in r_bins:
        if ru < 0.5 or ru > N_lattice / 2.5:
            continue
        shell = np.abs(r_m - ru) < 0.01
        n_shell = int(np.sum(shell))
        if n_shell == 0:
            continue
        phi_avg = np.mean(phi_m[shell])
        theory = 1.0 / (4 * pi * ru)
        rphi = ru * phi_avg
        ratio = rphi / target
        data.append((ru, phi_avg, rphi, ratio))
        print(f"  {ru:8.4f}  {phi_avg:12.8f}  {theory:12.8f}"
              f"  {rphi:10.8f}  {ratio:10.4f}")

    # Constancy in mid-range
    mid = [(ru, rphi, ratio) for ru, _, rphi, ratio in data
           if 2.0 < ru < N_lattice / 3]
    if mid:
        rphi_vals = [x[1] for x in mid]
        ratio_vals = [x[2] for x in mid]
        print(f"\n  MID-RANGE (2 < r < {N_lattice/3:.1f}):")
        print(f"    r*Phi mean   = {np.mean(rphi_vals):.8f}")
        print(f"    r*Phi std    = {np.std(rphi_vals):.8f}")
        print(f"    CV           = {np.std(rphi_vals)/abs(np.mean(rphi_vals)):.4f}")
        print(f"    ratio mean   = {np.mean(ratio_vals):.4f}")
        print(f"    Target       = 1.0000")


if __name__ == "__main__":
    solve_green(N_lattice=20)