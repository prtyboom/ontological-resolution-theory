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
        self.tick_created = 0
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
            neighbors = self.lattice.get_neighbors(idx)
            best_neighbor_dir = np.zeros(3)
            max_delta_W = 0
            
            for n_idx in neighbors:
                delta_W = field_W[n_idx] - field_W[idx]
                if delta_W > max_delta_W:
                    max_delta_W = delta_W
                    best_neighbor_dir = self.lattice.get_pos(n_idx) - self.lattice.get_pos(idx)
            
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
            distances = np.linalg.norm(self.lattice.positions - source_pos, axis=1)
            distances = np.maximum(distances, 0.1)
            contribution = (2 * self.Z * mass) / distances
            self.W_field += contribution
        return self.W_field


class Reaction:
    def __init__(self, reaction_type, reactants, products, tick_occurred):
        self.reaction_type = reaction_type
        self.reactants = reactants
        self.products = products
        self.tick_occurred = tick_occurred


class Nucleosynthesis:
    def __init__(self, N_lattice=50000, N_protons=100, N_electrons=100, N_neutrinos=50):
        print("="*60)
        print("PRIMORDIAL NUCLEOSYNTHESIS")
        print("="*60)
        
        self.lattice = FCCLattice(N_lattice)
        self.field_computer = FieldComputer(self.lattice)
        
        self.protons = []
        self.electrons = []
        self.neutrinos = []
        self.neutrons = []
        
        self.hydrogen_atoms = []
        self.deuterium = []
        self.helium = []
        
        self.reactions = []
        self.tick_count = 0
        
        self.R_merge = 4.0
        self.R_nucleus = 3.0
        
        self.N_protons = N_protons
        self.N_electrons = N_electrons
        self.N_neutrinos = N_neutrinos
    
    def setup(self):
        print(f"\nInitial conditions:")
        print(f"  Lattice: {self.lattice.L}³")
        print(f"  Protons: {self.N_protons}")
        print(f"  Electrons: {self.N_electrons}")
        print(f"  Neutrinos: {self.N_neutrinos}")
        
        margin = 2.0
        L = self.lattice.L - margin * 2
        
        for i in range(self.N_protons):
            pos = margin + np.random.random(3) * L
            p = GaussianCluster(pos, sigma=2.0, total_mass=1836.0, charge=+1, 
                               particle_type="proton", lattice=self.lattice)
            p.tick_created = 0
            self.protons.append(p)
        
        for i in range(self.N_electrons):
            pos = margin + np.random.random(3) * L
            e = GaussianCluster(pos, sigma=2.0, total_mass=1.0, charge=-1,
                               particle_type="electron", lattice=self.lattice)
            e.tick_created = 0
            self.electrons.append(e)
        
        for i in range(self.N_neutrinos):
            pos = margin + np.random.random(3) * L
            nu = GaussianCluster(pos, sigma=1.0, total_mass=0.01, charge=0,
                                particle_type="neutrino", lattice=self.lattice)
            nu.tick_created = 0
            self.neutrinos.append(nu)
        
        print(f"\n  All particles placed randomly")
    
    def compute_proton_field(self):
        sources = [(p.center, p.total_mass) for p in self.protons if p.active]
        sources += [(n.center, n.total_mass) for n in self.neutrons if n.active]
        return self.field_computer.compute_field(sources)
    
    def check_neutron_synthesis(self):
        """p + e + ν → n (inverse beta decay, neutrino commits)"""
        reactions = 0
        
        for nu in self.neutrinos:
            if not nu.active:
                continue
            
            for p in self.protons:
                if not p.active:
                    continue
                
                for e in self.electrons:
                    if not e.active:
                        continue
                    
                    d_pe = np.linalg.norm(p.center - e.center)
                    d_pnu = np.linalg.norm(p.center - nu.center)
                    d_enu = np.linalg.norm(e.center - nu.center)
                    
                    if d_pe < self.R_merge and d_pnu < self.R_merge and d_enu < self.R_merge:
                        neutron_pos = (p.center + e.center) / 2
                        n = GaussianCluster(
                            neutron_pos, sigma=2.0, total_mass=1839.0, charge=0,
                            particle_type="neutron", lattice=self.lattice
                        )
                        n.tick_created = self.tick_count
                        self.neutrons.append(n)
                        
                        p.active = False
                        e.active = False
                        nu.active = False
                        
                        self.reactions.append(Reaction(
                            "neutron_synthesis", [p, e, nu], [n], self.tick_count
                        ))
                        reactions += 1
                        break
                if not p.active:
                    break
        
        return reactions
    
    def check_hydrogen_formation(self):
        """p + e → H"""
        reactions = 0
        
        for e in self.electrons:
            if not e.active:
                continue
            
            for p in self.protons:
                if not p.active:
                    continue
                
                dist = np.linalg.norm(e.center - p.center)
                
                if dist < self.R_merge:
                    atom_pos = p.center.copy()
                    atom = GaussianCluster(
                        atom_pos, sigma=2.5, total_mass=1837.0, charge=0,
                        particle_type="hydrogen", lattice=self.lattice
                    )
                    atom.tick_created = self.tick_count
                    self.hydrogen_atoms.append(atom)
                    
                    e.active = False
                    p.active = False
                    
                    self.reactions.append(Reaction(
                        "hydrogen", [p, e], [atom], self.tick_count
                    ))
                    reactions += 1
                    break
        
        return reactions
    
    def check_deuterium_formation(self):
        """p + n → D"""
        reactions = 0
        
        for n in self.neutrons:
            if not n.active:
                continue
            
            for p in self.protons:
                if not p.active:
                    continue
                
                dist = np.linalg.norm(n.center - p.center)
                
                if dist < self.R_nucleus:
                    nucleus_pos = (p.center + n.center) / 2
                    d = GaussianCluster(
                        nucleus_pos, sigma=2.5, total_mass=3675.0, charge=+1,
                        particle_type="deuterium", lattice=self.lattice
                    )
                    d.tick_created = self.tick_count
                    self.deuterium.append(d)
                    
                    n.active = False
                    p.active = False
                    
                    self.reactions.append(Reaction(
                        "deuterium", [p, n], [d], self.tick_count
                    ))
                    reactions += 1
                    break
        
        return reactions
    
    def count_active(self):
        return {
            'protons': sum(1 for p in self.protons if p.active),
            'electrons': sum(1 for e in self.electrons if e.active),
            'neutrinos': sum(1 for nu in self.neutrinos if nu.active),
            'neutrons': sum(1 for n in self.neutrons if n.active),
            'hydrogen': len(self.hydrogen_atoms),
            'deuterium': len(self.deuterium),
            'helium': len(self.helium)
        }
    
    def evolve(self, ticks=500):
        print(f"\nEvolving for {ticks} Ticks...")
        
        checkpoints = [1, 10, 50, 100, 200, 500]
        start = time.time()
        
        for t in range(ticks):
            self.tick_count += 1
            
            W_field = self.compute_proton_field()
            
            for e in self.electrons:
                if e.active:
                    e.drift_step(W_field)
            
            for nu in self.neutrinos:
                if nu.active:
                    nu.drift_step(W_field)
            
            for n in self.neutrons:
                if n.active:
                    n.drift_step(W_field)
            
            r_neutron = self.check_neutron_synthesis()
            r_hydrogen = self.check_hydrogen_formation()
            r_deuterium = self.check_deuterium_formation()
            
            total_reactions = r_neutron + r_hydrogen + r_deuterium
            
            if (t+1) in checkpoints or total_reactions > 0:
                stats = self.count_active()
                print(f"  Tick {t+1}: p={stats['protons']} e={stats['electrons']} "
                      f"ν={stats['neutrinos']} n={stats['neutrons']} | "
                      f"H={stats['hydrogen']} D={stats['deuterium']} | "
                      f"reactions={total_reactions}")
        
        elapsed = time.time() - start
        print(f"\nCompleted in {elapsed:.2f} sec ({elapsed/self.tick_count:.4f} sec/Tick)")
    
    def report(self):
        print("\n" + "="*60)
        print("FINAL COMPOSITION")
        print("="*60)
        
        stats = self.count_active()
        
        print(f"\nFree particles:")
        print(f"  Protons: {stats['protons']}")
        print(f"  Electrons: {stats['electrons']}")
        print(f"  Neutrinos: {stats['neutrinos']}")
        print(f"  Neutrons: {stats['neutrons']}")
        
        print(f"\nBound states:")
        print(f"  Hydrogen (H): {stats['hydrogen']}")
        print(f"  Deuterium (D): {stats['deuterium']}")
        print(f"  Helium (He): {stats['helium']}")
        
        print(f"\nReaction summary:")
        reaction_types = {}
        for r in self.reactions:
            reaction_types[r.reaction_type] = reaction_types.get(r.reaction_type, 0) + 1
        
        for rtype, count in reaction_types.items():
            print(f"  {rtype}: {count}")


if __name__ == "__main__":
    sim = Nucleosynthesis(N_lattice=50000, N_protons=100, N_electrons=100, N_neutrinos=50)
    sim.setup()
    sim.evolve(ticks=500)
    sim.report()