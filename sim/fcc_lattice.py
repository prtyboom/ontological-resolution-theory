"""
ORT Simulation — Step 1: FCC Lattice Generator
===============================================
Axiom A3: FCC with k=12 coordination.
Builds a finite FCC lattice and neighbour table.
"""

import numpy as np

# FCC basis vectors (conventional cell, lattice constant a=1)
# Four atoms per conventional cell at:
#   (0,0,0), (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5)
FCC_BASIS = np.array([
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
])

# 12 nearest neighbours of origin in FCC (distance = a/sqrt(2))
FCC_NEIGHBOURS = np.array([
    [0.5,  0.5,  0.0],
    [0.5, -0.5,  0.0],
    [-0.5, 0.5,  0.0],
    [-0.5,-0.5,  0.0],
    [0.5,  0.0,  0.5],
    [0.5,  0.0, -0.5],
    [-0.5, 0.0,  0.5],
    [-0.5, 0.0, -0.5],
    [0.0,  0.5,  0.5],
    [0.0,  0.5, -0.5],
    [0.0, -0.5,  0.5],
    [0.0, -0.5, -0.5],
])


def build_fcc(N):
    """
    Build FCC lattice of N x N x N conventional cells.
    Returns:
      coords: (M, 3) array of node positions
      idx_map: dict mapping (x,y,z) tuple -> index
    """
    coords = []
    idx_map = {}

    for ix in range(N):
        for iy in range(N):
            for iz in range(N):
                origin = np.array([ix, iy, iz], dtype=float)
                for b in FCC_BASIS:
                    pos = origin + b
                    key = tuple(np.round(pos * 2).astype(int))
                    if key not in idx_map:
                        idx_map[key] = len(coords)
                        coords.append(pos)

    return np.array(coords), idx_map


def build_neighbours(coords, idx_map, N):
    """
    Build neighbour table with periodic boundary conditions.
    Returns:
      neighbours: (M, 12) array of neighbour indices
    """
    M = len(coords)
    neighbours = np.full((M, 12), -1, dtype=int)

    for i, pos in enumerate(coords):
        count = 0
        for dv in FCC_NEIGHBOURS:
            npos = pos + dv
            # Periodic wrapping
            npos_wrapped = npos % N
            key = tuple(np.round(npos_wrapped * 2).astype(int))
            if key in idx_map:
                neighbours[i, count] = idx_map[key]
                count += 1

    return neighbours


if __name__ == "__main__":
    N = 4  # 4x4x4 conventional cells
    coords, idx_map = build_fcc(N)
    neighbours = build_neighbours(coords, idx_map, N)

    M = len(coords)
    k_avg = np.mean(np.sum(neighbours >= 0, axis=1))

    print(f"FCC lattice: {N}x{N}x{N} cells")
    print(f"Total nodes: {M}")
    print(f"Expected: {4 * N**3}")
    print(f"Average coordination: {k_avg:.1f}")
    print(f"Target (A3): k=12")

    # Verify all nodes have 12 neighbours
    full = np.all(np.sum(neighbours >= 0, axis=1) == 12)
    print(f"All nodes k=12: {full}")