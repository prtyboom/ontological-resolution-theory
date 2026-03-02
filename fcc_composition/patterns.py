import numpy as np

class Pattern:
    """Base class for lattice patterns (particles)"""
    
    def __init__(self, name, mass_cost, charge):
        self.name = name
        self.mass_cost = mass_cost
        self.charge = charge
        self.signature = None
    
    def can_replicate(self, node_load, K_cell=133):
        return (node_load + self.mass_cost) <= K_cell
    
    def apply(self, lattice, node_idx):
        raise NotImplementedError


class Electron(Pattern):
    """Electron: nodal defect, charge -1"""
    
    def __init__(self):
        super().__init__("electron", mass_cost=4, charge=-1)
        self.signature = np.array([
            1, 0, 0, 0,  # charge=-1
            0, 1, 0, 0,  # nodal defect
            0, 0, 0, 0,
            0, 0, 0, 1
        ], dtype=np.int8)
    
    def apply(self, lattice, node_idx):
        if self.can_replicate(np.sum(lattice.states[node_idx])):
            lattice.states[node_idx] = self.signature.copy()
            return True
        return False


class UpQuark(Pattern):
    """Up quark: void defect, charge +2/3"""
    
    def __init__(self):
        super().__init__("up_quark", mass_cost=2, charge=2/3)
        self.signature = np.array([
            0, 1, 0, 0,  # charge=+2/3
            0, 0, 1, 0,  # void defect
            0, 0, 0, 0,
            0, 0, 1, 0
        ], dtype=np.int8)
    
    def apply(self, lattice, node_idx):
        current_load = np.sum(lattice.states[node_idx])
        if self.can_replicate(current_load):
            lattice.states[node_idx] = self.signature.copy()
            return True
        return False


class DownQuark(Pattern):
    """Down quark: void defect, charge -1/3"""
    
    def __init__(self):
        super().__init__("down_quark", mass_cost=3, charge=-1/3)
        self.signature = np.array([
            1, 0, 0, 1,  # charge=-1/3
            0, 0, 1, 0,  # void defect
            0, 0, 0, 0,
            0, 1, 0, 0
        ], dtype=np.int8)
    
    def apply(self, lattice, node_idx):
        current_load = np.sum(lattice.states[node_idx])
        if self.can_replicate(current_load):
            lattice.states[node_idx] = self.signature.copy()
            return True
        return False


class Proton:
    """Composite: 2 up + 1 down (uud)"""
    
    def __init__(self):
        self.quarks = [UpQuark(), UpQuark(), DownQuark()]
        self.charge = sum(q.charge for q in self.quarks)  # +1
        self.mass_cost = sum(q.mass_cost for q in self.quarks)  # 7 bits
        self.positions = []
    
    def can_form(self, lattice, center_idx):
        """Check if 3 nodes near center are available"""
        neighbors = lattice.get_neighbors(center_idx)
        
        # Need 3 empty neighbors
        available = []
        for n in neighbors:
            if np.sum(lattice.states[n]) == 0:
                available.append(n)
                if len(available) >= 3:
                    return True
        
        return False
    
    def apply(self, lattice, center_idx):
        """Form proton: place 3 quarks on neighboring nodes"""
        neighbors = lattice.get_neighbors(center_idx)
        
        # Find 3 empty neighbors
        available = []
        for n in neighbors:
            if np.sum(lattice.states[n]) == 0:
                available.append(n)
                if len(available) >= 3:
                    break
        
        if len(available) < 3:
            return False
        
        # Place quarks
        success_count = 0
        self.positions = []
        
        for i, quark in enumerate(self.quarks):
            if quark.apply(lattice, available[i]):
                self.positions.append(available[i])
                success_count += 1
        
        # Verify all 3 placed
        if success_count == 3:
            return True
        else:
            # Rollback if partial failure
            for pos in self.positions:
                lattice.states[pos] = np.zeros(16, dtype=np.int8)
            return False


# Test
if __name__ == "__main__":
    from fcc_lattice import FCCLattice
    
    print("Testing pattern definitions...")
    
    lattice = FCCLattice(1000)
    
    e = Electron()
    u = UpQuark()
    d = DownQuark()
    p = Proton()
    
    print(f"\nElectron:")
    print(f"  Mass: {e.mass_cost} bits, Charge: {e.charge}")
    e_ok = e.apply(lattice, 100)
    print(f"  Applied to node 100: {e_ok}")
    print(f"  Node 100 state sum: {np.sum(lattice.states[100])}")
    
    print(f"\nProton (uud):")
    print(f"  Mass: {p.mass_cost} bits, Charge: {p.charge}")
    p_ok = p.apply(lattice, 200)
    print(f"  Applied to node 200: {p_ok}")
    
    if p_ok:
        print(f"  Quark positions: {p.positions}")
        for pos in p.positions:
            print(f"    Node {pos} state sum: {np.sum(lattice.states[pos])}")