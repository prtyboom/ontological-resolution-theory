import numpy as np
from fcc_lattice import FCCLattice
from patterns import Electron, Proton
from operators import TickEngine
import time

class Hydrogen:
    """Hydrogen atom: proton + electron composite"""
    
    def __init__(self):
        self.proton = Proton()
        self.electron = Electron()
        self.total_charge = self.proton.charge + self.electron.charge  # 0
        self.total_mass = self.proton.mass_cost + self.electron.mass_cost  # 11 bits
        self.center = None
        self.e_position = None
        self.p_nodes = []
    
    def can_form(self, lattice, center_idx):
        """Check if hydrogen can form at center"""
        neighbors = lattice.get_neighbors(center_idx)
        void_available = sum(1 for n in neighbors 
                            if np.sum(lattice.states[n]) == 0)
        return void_available >= 4
    
    def apply(self, lattice, center_idx):
        """Form hydrogen at center"""
        if not self.can_form(lattice, center_idx):
            return False
        
        proton_formed = self.proton.apply(lattice, center_idx)
        if not proton_formed:
            return False
        
        # Track which nodes belong to proton
        neighbors = lattice.get_neighbors(center_idx)
        self.p_nodes = []
        for n in neighbors[:3]:
            if lattice.states[n][6] == 1:  # void defect flag
                self.p_nodes.append(n)
        
        # Place electron
        for n in neighbors[3:]:
            if np.sum(lattice.states[n]) == 0:
                electron_formed = self.electron.apply(lattice, n)
                if electron_formed:
                    self.center = center_idx
                    self.e_position = n
                    return True
        
        return False
    
    def is_stable(self, lattice):
        """Check if structure remains stable"""
        if self.center is None or self.e_position is None:
            return False
        
        # Check proton integrity
        quark_count = sum(1 for n in self.p_nodes 
                         if lattice.states[n][6] == 1)
        
        # Check electron present
        electron_present = (lattice.states[self.e_position][1] == 1)
        
        return (quark_count == 3) and electron_present


class CompositionExperiment:
    """Test pattern composition on FCC lattice"""
    
    def __init__(self, N=100000):
        print("="*60)
        print("COMPOSITION EXPERIMENT: Hydrogen Formation")
        print("="*60)
        
        self.lattice = FCCLattice(N)
        self.engine = TickEngine(self.lattice)
        self.atoms = []
        self.stats = {
            'attempted': 0,
            'successful': 0,
            'stable_after_100': 0,
            'stable_after_1000': 0
        }
    
    def seed_hydrogen(self, count=10):
        """Attempt to form multiple hydrogen atoms"""
        print(f"\nAttempting to form {count} hydrogen atoms...")
        
        # Select well-separated nodes
        spacing = int(self.lattice.N / (count * 2))
        indices = np.random.choice(
            range(1000, self.lattice.N - 1000), 
            count*5, 
            replace=False
        )
        
        for idx in indices:
            if self.stats['successful'] >= count:
                break
            
            h = Hydrogen()
            self.stats['attempted'] += 1
            
            if h.apply(self.lattice, idx):
                self.atoms.append(h)
                self.stats['successful'] += 1
                print(f"  Atom {self.stats['successful']}: "
                      f"formed at node {idx}")
        
        print(f"\nFormation success rate: "
              f"{self.stats['successful']}/{self.stats['attempted']}")
    
    def evolve(self, ticks):
        """Run lattice evolution and track stability"""
        print(f"\nEvolving for {ticks} Ticks...")
        
        checkpoints = {100: 'stable_after_100', 1000: 'stable_after_1000'}
        
        start = time.time()
        for t in range(ticks):
            self.engine.execute()
            
            if (t+1) in checkpoints:
                stable = sum(1 for atom in self.atoms 
                           if atom.is_stable(self.lattice))
                self.stats[checkpoints[t+1]] = stable
                print(f"  Tick {t+1}: {stable}/{len(self.atoms)} atoms stable")
        
        elapsed = time.time() - start
        print(f"\nEvolution completed in {elapsed:.2f} sec "
              f"({elapsed/ticks:.4f} sec/Tick)")
    
    def report(self):
        """Final statistics"""
        print("\n" + "="*60)
        print("RESULTS")
        print("="*60)
        print(f"Lattice size: {self.lattice.N} nodes")
        print(f"Atoms formed: {self.stats['successful']}/{self.stats['attempted']}")
        print(f"Stable after 100 Ticks: {self.stats['stable_after_100']}")
        print(f"Stable after 1000 Ticks: {self.stats['stable_after_1000']}")
        
        if self.stats['successful'] > 0:
            retention = self.stats['stable_after_1000'] / self.stats['successful']
            print(f"Retention rate: {retention*100:.1f}%")


# Run experiment
if __name__ == "__main__":
    exp = CompositionExperiment(N=100000)
    exp.seed_hydrogen(count=10)
    exp.evolve(ticks=1000)
    exp.report()