import numpy as np
from scipy.spatial import cKDTree
import time

class FCCLattice:
    """Face-Centered Cubic lattice with k=12 coordination"""
    
    def __init__(self, N):
        """
        N: approximate number of nodes (will adjust to fit FCC structure)
        """
        cells_per_side = int((N / 4) ** (1/3)) + 1
        self.L = cells_per_side
        self.N = 4 * cells_per_side ** 3
        
        print(f"Initializing FCC lattice:")
        print(f"  Target nodes: {N}")
        print(f"  Actual nodes: {self.N}")
        print(f"  Grid size: {self.L}x{self.L}x{self.L}")
        
        print(f"  Generating positions...")
        self.positions = self._generate_positions()
        
        print(f"  Building neighbor index...")
        self.tree = cKDTree(self.positions)
        
        print(f"  Computing 12 nearest neighbors...")
        self._precompute_neighbors()
        
        self.states = np.zeros((self.N, 16), dtype=np.int8)
        
        mem_kb = (self.N * 16 + self.N * 12 * 4) / 1024
        print(f"  Memory: ~{mem_kb:.1f} KB")
    
    def _generate_positions(self):
        """Generate FCC lattice positions"""
        positions = []
        
        basis = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5]
        ])
        
        for i in range(self.L):
            for j in range(self.L):
                for k in range(self.L):
                    origin = np.array([i, j, k])
                    for b in basis:
                        positions.append(origin + b)
        
        return np.array(positions[:self.N], dtype=np.float32)
    
    def _precompute_neighbors(self):
        """Pre-compute 12 nearest neighbors for each node"""
        distances, indices = self.tree.query(self.positions, k=13)
        self.neighbors = indices[:, 1:13].astype(np.int32)
        self.neighbor_distances = distances[:, 1:13].astype(np.float32)
    
    def get_neighbors(self, node_idx):
        """Get 12 nearest neighbors for FCC node"""
        return self.neighbors[node_idx]
    
    def get_neighbor_distances(self, node_idx):
        """Get distances to 12 nearest neighbors"""
        return self.neighbor_distances[node_idx]
    
    def get_pos(self, node_idx):
        """Get position of node"""
        return self.positions[node_idx]
    
    def get_nodes_in_region(self, center, radius):
        """Find all nodes within radius of center"""
        distances = np.linalg.norm(self.positions - center, axis=1)
        indices = np.where(distances <= radius)[0]
        return indices
    
    def count_occupied_neighbors(self, node_idx):
        """Count how many neighbors have non-zero state"""
        neighbors = self.neighbors[node_idx]
        count = 0
        for n in neighbors:
            if np.sum(self.states[n]) > 0:
                count += 1
        return count
    
    def tick(self):
        """One global lattice update (placeholder)"""
        return True


if __name__ == "__main__":
    start = time.time()
    lattice = FCCLattice(100000)
    init_time = time.time() - start
    
    print(f"\nInitialization time: {init_time:.2f} sec")
    
    print(f"\nVerifying FCC structure...")
    
    for idx in [0, 100, 1000]:
        neighbors = lattice.get_neighbors(idx)
        distances = lattice.get_neighbor_distances(idx)
        
        print(f"\n  Node {idx}:")
        print(f"    Position: {lattice.positions[idx]}")
        print(f"    Neighbors: {neighbors[:4]}... (showing 4 of 12)")
        print(f"    Distances: {distances[:4]}... (should be ~0.707)")
    
    expected_dist = np.sqrt(2) / 2
    tolerance = 0.01
    
    correct = 0
    for idx in range(min(1000, lattice.N)):
        distances = lattice.get_neighbor_distances(idx)
        if np.all(np.abs(distances - expected_dist) < tolerance):
            correct += 1
    
    print(f"\n  FCC verification: {correct}/1000 nodes have correct neighbor distances")
    
    print(f"\nTesting new methods...")
    center = lattice.positions[5000]
    region = lattice.get_nodes_in_region(center, 5.0)
    print(f"  Nodes within radius 5.0 of node 5000: {len(region)}")
    
    pos = lattice.get_pos(100)
    print(f"  Position of node 100: {pos}")
    
    print(f"\nRunning 10 Ticks...")
    start = time.time()
    for t in range(10):
        lattice.tick()
    tick_time = (time.time() - start) / 10
    print(f"Average time per Tick: {tick_time:.4f} sec")