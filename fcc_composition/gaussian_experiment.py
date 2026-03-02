import numpy as np
from fcc_lattice import FCCLattice
import time

class GaussianCluster:
    def __init__(self, center, sigma=10.0, total_mass=1.0, lattice=None):
        self.center = np.array(center, dtype=float)
        self.sigma = sigma
        self.total_mass = total_mass
        self.lattice = lattice
        self.nodes = {}
        self.assign_weights()

    def assign_weights(self):
        self.nodes = {}
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
            print(f"WARNING: No nodes found near center {self.center}")
            return
        
        norm_factor = self.total_mass / np.sum(raw_weights)
        for i, idx in enumerate(indices):
            self.nodes[idx] = raw_weights[i] * norm_factor

    def compute_center(self):
        if len(self.nodes) == 0:
            return self.center
        
        weighted_pos_sum = np.zeros(3)
        total_w = 0
        for idx, w in self.nodes.items():
            weighted_pos_sum += self.lattice.get_pos(idx) * w
            total_w += w
        
        if total_w > 0:
            self.center = weighted_pos_sum / total_w
        return self.center

    def drift_step(self, field_W):
        if len(self.nodes) == 0:
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

    def get_total_load_at_node(self, node_idx):
        return self.nodes.get(node_idx, 0.0)


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


class HydrogenFormationExperiment:
    def __init__(self, N=10000):
        print("="*60)
        print("HYDROGEN FORMATION: Gaussian Cluster Dynamics")
        print("="*60)
        
        self.lattice = FCCLattice(N)
        self.field_computer = FieldComputer(self.lattice)
        
        self.proton = None
        self.electron = None
        
        self.distances = []
        self.tick_count = 0
        self.merged = False
    
    def setup(self, separation=5.0):
        center_idx = self.lattice.N // 2
        center_pos = self.lattice.positions[center_idx]
        
        max_sep = self.lattice.L / 3.0
        actual_sep = min(separation, max_sep)
        
        print(f"\nInitializing particles:")
        print(f"  Lattice size: {self.lattice.L}x{self.lattice.L}x{self.lattice.L}")
        print(f"  Proton mass: 1836 (electron mass units)")
        print(f"  Electron mass: 1")
        print(f"  Sigma: 2.0 lattice units")
        print(f"  Requested separation: {separation:.1f}")
        print(f"  Actual separation: {actual_sep:.1f}")
        
        self.proton = GaussianCluster(
            center=center_pos,
            sigma=2.0,
            total_mass=1836.0,
            lattice=self.lattice
        )
        
        electron_pos = center_pos + np.array([actual_sep, 0, 0])
        
        self.electron = GaussianCluster(
            center=electron_pos,
            sigma=2.0,
            total_mass=1.0,
            lattice=self.lattice
        )
        
        print(f"\n  Proton center: {self.proton.center}")
        print(f"  Proton nodes: {len(self.proton.nodes)}")
        print(f"  Electron center: {self.electron.center}")
        print(f"  Electron nodes: {len(self.electron.nodes)}")
        
        initial_dist = np.linalg.norm(
            self.electron.center - self.proton.center
        )
        self.distances.append(initial_dist)
        print(f"  Measured separation: {initial_dist:.2f}")
        
        self.R_crit = 2.0 * (self.proton.sigma + self.electron.sigma)
        print(f"  Critical merge radius: {self.R_crit:.2f}")
    
    def check_merger(self):
        dist = np.linalg.norm(self.electron.center - self.proton.center)
        if dist < self.R_crit:
            return True
        return False
    
    def evolve(self, ticks=1000):
        print(f"\nEvolving for {ticks} Ticks...")
        
        checkpoints = [10, 50, 100, 500, 1000]
        start = time.time()
        
        for t in range(ticks):
            self.tick_count += 1
            
            if self.merged:
                break
            
            sources = [(self.proton.center, self.proton.total_mass)]
            W_field = self.field_computer.compute_field(sources)
            
            self.electron.drift_step(W_field)
            
            dist = np.linalg.norm(
                self.electron.center - self.proton.center
            )
            self.distances.append(dist)
            
            if self.check_merger():
                print(f"\n  *** MERGER at Tick {t+1}! Distance = {dist:.4f} ***")
                self.merged = True
                break
            
            if (t+1) in checkpoints:
                print(f"  Tick {t+1}: distance = {dist:.4f}, "
                      f"electron nodes = {len(self.electron.nodes)}")
        
        elapsed = time.time() - start
        print(f"\nCompleted in {elapsed:.2f} sec ({elapsed/self.tick_count:.4f} sec/Tick)")
    
    def report(self):
        print("\n" + "="*60)
        print("RESULTS")
        print("="*60)
        
        final_dist = self.distances[-1]
        initial_dist = self.distances[0]
        
        print(f"Initial separation: {initial_dist:.4f}")
        print(f"Final separation: {final_dist:.4f}")
        print(f"Change: {initial_dist - final_dist:.4f} units")
        
        if initial_dist > 0:
            pct = 100 * (initial_dist - final_dist) / initial_dist
            print(f"Percent reduction: {pct:.2f}%")
        
        if self.merged:
            print(f"\nStatus: MERGED INTO HYDROGEN ATOM at Tick {self.tick_count}")
            print(f"  Coherent state I_atom < M_p + m_e")
        elif len(self.distances) >= 100:
            last_100 = self.distances[-100:]
            variance = np.var(last_100)
            
            if variance < 0.0001:
                print(f"\nStatus: CONVERGED (variance = {variance:.6f})")
            elif final_dist < initial_dist:
                print(f"\nStatus: APPROACHING (variance = {variance:.4f})")
            else:
                print(f"\nStatus: NO MOVEMENT")


if __name__ == "__main__":
    exp = HydrogenFormationExperiment(N=10000)
    exp.setup(separation=5.0)
    exp.evolve(ticks=1000)
    exp.report()