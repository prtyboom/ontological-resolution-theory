import numpy as np
from fcc_lattice import FCCLattice
from patterns import Electron, Proton
import matplotlib.pyplot as plt
import time

# ORT Constants
Z = 137.036  # Geometric impedance
M_PROTON = 1836  # Proton mass in electron mass units
M_ELECTRON = 1  # Electron mass (reference)


class InformationField:
    """Calculate information weight field W_eff(r) from massive defects"""
    
    def __init__(self, lattice):
        self.lattice = lattice
        self.W = np.zeros(lattice.N, dtype=np.float32)
        self.sources = []  # List of (position, mass) tuples
    
    def add_source(self, node_idx, mass):
        """Add a massive defect as field source"""
        self.sources.append((node_idx, mass))
        self._update_field()
    
    def _update_field(self):
        """Compute W_eff(r) = sum over sources of (2*Z*M/r)"""
        self.W.fill(0)
        
        for source_idx, mass in self.sources:
            source_pos = self.lattice.positions[source_idx]
            
            # Calculate distance from source to all nodes
            distances = np.linalg.norm(
                self.lattice.positions - source_pos, 
                axis=1
            )
            
            # Avoid division by zero at source position
            distances[source_idx] = 1e-10
            
            # W_eff = 2*Z*M/r
            contribution = (2 * Z * mass) / distances
            
            self.W += contribution
    
    def get_weight(self, node_idx):
        """Get effective weight at node"""
        return self.W[node_idx]
    
    def get_update_probability(self, node_idx):
        """Calculate U = (1 - W/Z)^2"""
        w = self.W[node_idx]
        ratio = w / Z
        
        if ratio >= 1.0:
            return 0.0  # Saturated - no updates possible
        
        U = (1.0 - ratio) ** 2
        return U


class MigrationEngine:
    """Handle pattern migration based on field gradients"""
    
    def __init__(self, lattice, field):
        self.lattice = lattice
        self.field = field
        self.pattern_positions = {}  # {pattern_id: node_idx}
        self.tick_count = 0
    
    def register_pattern(self, pattern_id, node_idx, mass):
        """Register a pattern for migration tracking"""
        self.pattern_positions[pattern_id] = {
            'position': node_idx,
            'mass': mass,
            'history': [node_idx]
        }
    
    def migrate_step(self, pattern_id):
        """Execute one migration step for pattern"""
        if pattern_id not in self.pattern_positions:
            return False
        
        current_pos = self.pattern_positions[pattern_id]['position']
        
        # Get update probability at current position
        U = self.field.get_update_probability(current_pos)
        
        # Roll dice - can this pattern move this Tick?
        if np.random.random() > U:
            return False  # Too much latency - no move
        
        # Get 12 neighbors
        neighbors = self.lattice.get_neighbors(current_pos)
        
        # Find neighbor with minimum W (steepest descent)
        weights = [self.field.get_weight(n) for n in neighbors]
        min_idx = np.argmin(weights)
        target_node = neighbors[min_idx]
        
        # Check if target is better than current
        current_w = self.field.get_weight(current_pos)
        target_w = weights[min_idx]
        
        if target_w < current_w:
            # Migrate
            self.pattern_positions[pattern_id]['position'] = target_node
            self.pattern_positions[pattern_id]['history'].append(target_node)
            return True
        
        return False  # Already at local minimum
    
    def evolve(self, ticks):
        """Run migration for all patterns over multiple ticks"""
        migration_counts = {pid: 0 for pid in self.pattern_positions}
        
        for t in range(ticks):
            self.tick_count += 1
            
            for pattern_id in self.pattern_positions:
                moved = self.migrate_step(pattern_id)
                if moved:
                    migration_counts[pattern_id] += 1
        
        return migration_counts
    
    def get_distance_to_source(self, pattern_id, source_idx):
        """Calculate distance between pattern and field source"""
        if pattern_id not in self.pattern_positions:
            return None
        
        pattern_pos = self.lattice.positions[
            self.pattern_positions[pattern_id]['position']
        ]
        source_pos = self.lattice.positions[source_idx]
        
        distance = np.linalg.norm(pattern_pos - source_pos)
        return distance


# Experiment: Electron falling toward Proton
class HydrogenFormationExperiment:
    """Test electron migration toward proton via field gradient"""
    
    def __init__(self, N=10000):
        print("="*60)
        print("HYDROGEN FORMATION: Field Dynamics")
        print("="*60)
        
        self.lattice = FCCLattice(N)
        self.field = InformationField(self.lattice)
        self.engine = MigrationEngine(self.lattice, self.field)
        
        self.proton_pos = None
        self.electron_id = "e1"
    
    def setup(self):
        """Place proton at center, electron at distance"""
        # Place proton near center
        center_idx = self.lattice.N // 2
        self.proton_pos = center_idx
        
        print(f"\nPlacing proton at node {center_idx}")
        self.field.add_source(center_idx, M_PROTON)
        
        # Place electron at distance ~50 cells
        electron_pos = center_idx + 2000  # arbitrary offset
        if electron_pos >= self.lattice.N:
            electron_pos = center_idx - 2000
        
        print(f"Placing electron at node {electron_pos}")
        self.engine.register_pattern(self.electron_id, electron_pos, M_ELECTRON)
        
        # Initial distance
        initial_dist = self.engine.get_distance_to_source(
            self.electron_id, 
            self.proton_pos
        )
        print(f"Initial separation: {initial_dist:.2f} cells")
    
    def run(self, ticks=1000):
        """Run migration simulation"""
        print(f"\nEvolving for {ticks} Ticks...")
        
        checkpoints = [100, 500, 1000]
        start = time.time()
        
        for checkpoint in checkpoints:
            remaining = checkpoint - self.engine.tick_count
            if remaining > 0:
                self.engine.evolve(remaining)
            
            dist = self.engine.get_distance_to_source(
                self.electron_id, 
                self.proton_pos
            )
            migrations = len(
                self.engine.pattern_positions[self.electron_id]['history']
            ) - 1
            
            print(f"  Tick {checkpoint}: distance = {dist:.2f}, "
                  f"migrations = {migrations}")
        
        elapsed = time.time() - start
        print(f"\nCompleted in {elapsed:.2f} sec")
    
    def report(self):
        """Final analysis"""
        print("\n" + "="*60)
        print("RESULTS")
        print("="*60)
        
        final_dist = self.engine.get_distance_to_source(
            self.electron_id, 
            self.proton_pos
        )
        
        history = self.engine.pattern_positions[self.electron_id]['history']
        total_migrations = len(history) - 1
        
        print(f"Final separation: {final_dist:.2f} cells")
        print(f"Total migrations: {total_migrations}")
        print(f"Migration rate: {total_migrations / self.engine.tick_count * 100:.1f}%")
        
        # Check if reached equilibrium (orbital)
        if total_migrations > 0:
            last_10_positions = history[-10:]
            unique_last_10 = len(set(last_10_positions))
            
            if unique_last_10 < 5:
                print("Status: STABLE ORBIT (oscillating)")
            elif final_dist < 5.0:
                print("Status: BOUND STATE")
            else:
                print("Status: STILL APPROACHING")


# Run
if __name__ == "__main__":
    exp = HydrogenFormationExperiment(N=10000)
    exp.setup()
    exp.run(ticks=1000)
    exp.report()