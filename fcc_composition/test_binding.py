from fcc_lattice import FCCLattice
from patterns import Electron, Proton
from operators import UpdateRules
import numpy as np

lattice = FCCLattice(1000)

# Place electron and proton
e = Electron()
p = Proton()

e_ok = e.apply(lattice, 100)
p_ok = p.apply(lattice, 200)

print("Formation results:")
print(f"  Electron at 100: {e_ok}")
print(f"  Proton at 200: {p_ok}")

if p_ok:
    print(f"  Proton quark positions: {p.positions}")

# Check if electron and first quark are neighbors
neighbors_100 = lattice.get_neighbors(100)
quark_0 = p.positions[0] if p_ok else None

print(f"\nNode 100 neighbors: {neighbors_100}")
if quark_0:
    print(f"First quark at: {quark_0}")
    print(f"Are they neighbors? {quark_0 in neighbors_100}")

# Count occupied neighbors for electron
occupied = sum(1 for n in neighbors_100 if np.sum(lattice.states[n]) > 0)
print(f"\nElectron has {occupied} occupied neighbors")

# Check stability
stable_e = UpdateRules.stability_check(lattice, 100, "e")
print(f"Electron stable: {stable_e}")

if p_ok:
    for i, pos in enumerate(p.positions):
        stable_q = UpdateRules.stability_check(lattice, pos, "q")
        print(f"Quark {i} at {pos} stable: {stable_q}")