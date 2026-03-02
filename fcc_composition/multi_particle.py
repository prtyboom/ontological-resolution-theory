import numpy as np
from fcc_lattice import FCCLattice
import time

class GaussianCluster:
    def __init__(self, center, sigma=2.0, total_mass=1.0, charge=0, particle_type="unknown", lattice=None):
        self.center = np.array(center, dtype=float)
        self.sigma = sigma
        self.total_mass = total_mass
        self.charge = charge
        self.particle_type = particle_type
        self.lattice = lattice
        self.nodes = {}
        self.bound_to = None
        self.active = True
        self.assign_weights()

    def assign_weights(self):
        self.nodes = {}
        if not self.active:
            return
        
        radius_limit = 3 * self.sigma
        nearby_indices = self.lattice.get_nodes_in_region(self.center, radius_limit)
        
        raw_weights = []
        indices = []
        
        for idx in nearby_indices:
            pos = self.lattice.get_pos(idx)
            r_squared = np.sum((pos - self.center)**2)
            weight = np.exp(-r_squared / (2 * self.sigma**2))
            raw_weights.append(weight)
            indices.append(idx)
        
        if len(raw_weights) == 0:
            return
        
        norm_factor = self.total_mass / np.sum(raw_weights)
        for i, idx in enumerate(indices):
            self.nodes[idx] = raw_weights[i] * norm_factor

    def drift_step(self, field_W):
        if not self.active or len(self.nodes) == 0:
            return
        
        drift_vector = np.zeros(3)
        
        for idx, w in self.nodes.items():
            current_pos = self.lattice.get_pos(idx)
            neighbors = self.lattice.get_neighbors(idx)
            
            best_neighbor_dir = np.zeros(3)
            max_delta_W = 0
            
            for n_idx in neighbors:
                delta_W = field_W[n_idx] - field_W[idx]
                if delta_W > max_delta_W:
                    max_delta_W = delta_W
                    best_neighbor_dir = self.lattice.get_pos(n_idx) - current_pos
            
            drift_vector += best_neighbor_dir * w
        
        if np.linalg.norm(drift_vector) > 0:
            step = drift_vector / np.linalg.norm(drift_vector)
            self.center += step
            self.assign_weights()


class FieldComputer:
    def __init__(self, lattice, Z=137.036):
        self.lattice = lattice
        self.Z = Z
        self.W_field = np.zeros(lattice.N, dtype=np.float32)
    
    def compute_field(self, sources):
        self.W_field.fill(0)
        
        for source_pos, mass in sources:
            distances = np.linalg.norm(
                self.lattice.positions - source_pos, 
                axis=1
            )
            distances = np.maximum(distances, 0.1)
            contribution = (2 * self.Z * mass) / distances
            self.W_field += contribution
        
        return self.W_field


class Atom:
    def __init__(self, proton, electron, tick_formed):
        self.proton = proton
        self.electron = electron
        self.tick_formed = tick_formed
        self.center = proton.center.copy()
        self.total_mass = proton.total_mass + electron.total_mass
        self.charge = 0


class PrimordialSoup:
    def __init__(self, N_lattice=50000, N_protons=50, N_electrons=50):
        print("="*60)
        print("PRIMORDIAL SOUP: Multi-Particle Dynamics")
        print("="*60)
        
        self.lattice = FCCLattice(N_lattice)
        self.field_computer = FieldComputer(self.lattice)
        
        self.protons = []
        self.electrons = []
        self.atoms = []
        
        self.N_protons = N_protons
        self.N_electrons = N_electrons
        
        self.tick_count = 0
        self.R_crit = 4.0
        
        self.stats = {
            'mergers_per_tick': [],
            'active_electrons': [],
            'active_protons': []
        }
    
    def setup(self):
        print(f"\nInitializing particles:")
        print(f"  Lattice: {self.lattice.L}x{self.lattice.L}x{self.lattice.L}")
        print(f"  Protons: {self.N_protons}")
        print(f"  Electrons: {self.N_electrons}")
        print(f"  Merge radius: {self.R_crit}")
        
        margin = 2.0
        L = self.lattice.L - margin * 2
        
        for i in range(self.N_protons):
            pos = np.array([
                margin + np.random.random() * L,
                margin + np.random.random() * L,
                margin + np.random.random() * L
            ])
            p = GaussianCluster(
                center=pos,
                sigma=2.0,
                total_mass=1836.0,
                charge=+1,
                particle_type="proton",
                lattice=self.lattice
            )
            self.protons.append(p)
        
        for i in range(self.N_electrons):
            pos = np.array([
                margin + np.random.random() * L,
                margin + np.random.random() * L,
                margin + np.random.random() * L
            ])
            e = GaussianCluster(
                center=pos,
                sigma=2.0,
                total_mass=1.0,
                charge=-1,
                particle_type="electron",
                lattice=self.lattice
            )
            self.electrons.append(e)
        
        print(f"\n  Particles placed randomly in volume")
        print(f"  Total proton nodes: {sum(len(p.nodes) for p in self.protons)}")
        print(f"  Total electron nodes: {sum(len(e.nodes) for e in self.electrons)}")
    
    def compute_proton_field(self):
        sources = []
        for p in self.protons:
            if p.active:
                sources.append((p.center, p.total_mass))
        return self.field_computer.compute_field(sources)
    
    def check_mergers(self):
        mergers_this_tick = 0
        
        for e in self.electrons:
            if not e.active:
                continue
            
            for p in self.protons:
                if not p.active:
                    continue
                
                dist = np.linalg.norm(e.center - p.center)
                
                if dist < self.R_crit:
                    atom = Atom(p, e, self.tick_count)
                    self.atoms.append(atom)
                    
                    e.active = False
                    e.bound_to = p
                    p.active = False
                    p.bound_to = e
                    
                    mergers_this_tick += 1
                    break
        
        return mergers_this_tick
    
    def count_active(self):
        active_p = sum(1 for p in self.protons if p.active)
        active_e = sum(1 for e in self.electrons if e.active)
        return active_p, active_e
    
    def evolve(self, ticks=500):
        print(f"\nEvolving for {ticks} Ticks...")
        
        checkpoints = [10, 50, 100, 200, 500]
        start = time.time()
        
        for t in range(ticks):
            self.tick_count += 1
            
            active_p, active_e = self.count_active()
            if active_e == 0:
                print(f"\n  All electrons bound at Tick {self.tick_count}")
                break
            
            W_field = self.compute_proton_field()
            
            for e in self.electrons:
                if e.active:
                    e.drift_step(W_field)
            
            mergers = self.check_mergers()
            
            self.stats['mergers_per_tick'].append(mergers)
            self.stats['active_electrons'].append(active_e)
            self.stats['active_protons'].append(active_p)
            
            if (t+1) in checkpoints or mergers > 0:
                print(f"  Tick {t+1}: atoms={len(self.atoms)}, "
                      f"active_e={active_e}, active_p={active_p}, "
                      f"mergers={mergers}")
        
        elapsed = time.time() - start
        print(f"\nCompleted in {elapsed:.2f} sec ({elapsed/self.tick_count:.4f} sec/Tick)")
    
    def report(self):
        print("\n" + "="*60)
        print("RESULTS")
        print("="*60)
        
        active_p, active_e = self.count_active()
        
        print(f"Total Ticks: {self.tick_count}")
        print(f"Atoms formed: {len(self.atoms)}")
        print(f"Remaining free protons: {active_p}")
        print(f"Remaining free electrons: {active_e}")
        
        if len(self.atoms) > 0:
            formation_ticks = [a.tick_formed for a in self.atoms]
            print(f"\nFormation statistics:")
            print(f"  First atom: Tick {min(formation_ticks)}")
            print(f"  Last atom: Tick {max(formation_ticks)}")
            print(f"  Mean formation time: {np.mean(formation_ticks):.1f} Ticks")
        
        efficiency = len(self.atoms) / min(self.N_protons, self.N_electrons) * 100
        print(f"\nBinding efficiency: {efficiency:.1f}%")
        
        if active_p > 0 or active_e > 0:
            print(f"\nUnbound particles:")
            if active_p > 0:
                print(f"  {active_p} protons seeking electrons")
            if active_e > 0:
                print(f"  {active_e} electrons seeking protons")


if __name__ == "__main__":
    soup = PrimordialSoup(N_lattice=50000, N_protons=50, N_electrons=50)
    soup.setup()
    soup.evolve(ticks=500)
    soup.report()