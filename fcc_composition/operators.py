import numpy as np

# Fundamental constants (ORT Canon v22.0)
PHI = (1 + np.sqrt(5)) / 2  # 1.618033988749895
E = np.e                      # 2.718281828459045
PI = np.pi                    # 3.141592653589793
K = 12                        # FCC coordination number

# Derived operators
DISCRETE_FLOW = 4 * PHI**3           # 16.944 (tau -> mu)
CONTINUOUS_SWEEP = 12 * 2 * PI * E   # 204.95 (mu -> e)
APERTURE_GAP = 2 * PI * E - 4 * PHI**3  # 0.135 (boundary)
ALPHA_INV = PI**4 + 4*PI**2 + (PI-3) + 1/144  # 137.036
K_CELL = 133  # bits per node

print("ORT Operators loaded:")
print(f"  PHI = {PHI:.6f}")
print(f"  4*PHI^3 = {DISCRETE_FLOW:.3f}")
print(f"  12*2*PI*e = {CONTINUOUS_SWEEP:.2f}")
print(f"  Aperture Gap = {APERTURE_GAP:.6f}")
print(f"  Alpha^-1 = {ALPHA_INV:.6f}")
print(f"  K_cell = {K_CELL} bits")


class UpdateRules:
    """Lattice update rules based on ORT operators"""
    
    @staticmethod
    def binding_energy(pattern_a, pattern_b, distance):
        if distance == 0:
            return 0
        charge_product = pattern_a['charge'] * pattern_b['charge']
        energy = charge_product / distance
        return energy * APERTURE_GAP
    
    @staticmethod
    def can_replicate(mass_cost, available_bits):
        return mass_cost <= available_bits
    
    @staticmethod
    def replication_probability(mass_cost):
        load_fraction = mass_cost / K_CELL
        prob = 1.0 - load_fraction
        return max(0.0, min(1.0, prob))
    
    @staticmethod
    def stability_check(lattice, node_idx, pattern_type):
        state = lattice.states[node_idx]
        if np.sum(state) == 0:
            return True
        
        neighbors = lattice.get_neighbors(node_idx)
        bound_count = 0
        total_binding = 0.0
        
        for n in neighbors:
            n_state = lattice.states[n]
            if np.sum(n_state) > 0:
                my_charge = UpdateRules._decode_charge(state)
                n_charge = UpdateRules._decode_charge(n_state)
                binding = my_charge * n_charge * (-1)
                total_binding += binding
                if binding < 0:
                    bound_count += 1
        
        is_stable = (bound_count >= 1) or (total_binding < 0)
        return is_stable
    
    @staticmethod
    def _decode_charge(state):
        if state[0] == 1 and state[3] == 0:
            return -1.0
        elif state[0] == 0 and state[1] == 1:
            return 2/3
        elif state[0] == 1 and state[3] == 1:
            return -1/3
        else:
            return 0.0
    
    @staticmethod
    def decay_probability(mass_cost, bound_count):
        mass_factor = mass_cost / K_CELL
        binding_factor = 1.0 / (1 + bound_count)
        prob = mass_factor * binding_factor
        return max(0.0, min(1.0, prob))


class TickEngine:
    """Execute one Tick of lattice evolution"""
    
    def __init__(self, lattice):
        self.lattice = lattice
        self.rules = UpdateRules()
        self.tick_count = 0
        self.decay_events = 0
        self.stable_events = 0
    
    def execute(self):
        self.tick_count += 1
        new_states = self.lattice.states.copy()
        
        for idx in range(self.lattice.N):
            state = self.lattice.states[idx]
            if np.sum(state) == 0:
                continue
            
            is_stable = self.rules.stability_check(
                self.lattice, idx, "generic"
            )
            
            if is_stable:
                self.stable_events += 1
                new_states[idx] = state
            else:
                mass_cost = np.sum(state)
                neighbors = self.lattice.get_neighbors(idx)
                bound_count = sum(1 for n in neighbors 
                                 if np.sum(self.lattice.states[n]) > 0)
                
                decay_prob = self.rules.decay_probability(mass_cost, bound_count)
                
                if np.random.random() < decay_prob:
                    new_states[idx] = np.zeros(16, dtype=np.int8)
                    self.decay_events += 1
                else:
                    new_states[idx] = state
                    self.stable_events += 1
        
        self.lattice.states = new_states
        
        return {
            'tick': self.tick_count,
            'decays': self.decay_events,
            'stable': self.stable_events
        }
    
    def reset_counters(self):
        self.decay_events = 0
        self.stable_events = 0


# Test
if __name__ == "__main__":
    from fcc_lattice import FCCLattice
    from patterns import Electron, Proton
    
    print("\n" + "="*60)
    print("TESTING OPERATORS")
    print("="*60)
    
    lattice = FCCLattice(1000)
    engine = TickEngine(lattice)
    
    # Place electron and proton nearby
    e = Electron()
    p = Proton()
    
    e_formed = e.apply(lattice, 100)
    p_formed = p.apply(lattice, 101)
    
    print(f"\nElectron formed: {e_formed}")
    print(f"Proton formed: {p_formed}")
    
    initial_patterns = sum(1 for i in range(lattice.N) 
                          if np.sum(lattice.states[i]) > 0)
    print(f"Initial patterns: {initial_patterns}")
    
    print("\nRunning 100 Ticks...")
    for t in range(100):
        stats = engine.execute()
        if (t+1) % 20 == 0:
            patterns = sum(1 for i in range(lattice.N) 
                          if np.sum(lattice.states[i]) > 0)
            print(f"  Tick {t+1}: {patterns} patterns remain")
    
    final_patterns = sum(1 for i in range(lattice.N) 
                        if np.sum(lattice.states[i]) > 0)
    print(f"\nFinal patterns: {final_patterns}")
    print(f"Retention: {100*final_patterns/max(1,initial_patterns):.1f}%")