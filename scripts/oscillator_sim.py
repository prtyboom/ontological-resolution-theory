"""
Oscillator Simulation v1.0
Один осциллятор: (phi, sigma) -> Sep -> Phase -> (phi', sigma')
Решётка на T^2, проверка (pi-3) дефекта
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# КОНСТАНТЫ
# ============================================================

PI = np.pi
PI_MINUS_3 = PI - 3  # 0.14159265...

# ============================================================
# ОСЦИЛЛЯТОР
# ============================================================

class Oscillator:
    """Один осциллятор на мембране."""
    
    def __init__(self, phi=None, sigma=0):
        """
        phi: фаза [0, 2*pi)
        sigma: дискретное состояние (0 или 1)
        """
        if phi is None:
            phi = np.random.uniform(0, 2 * PI)
        self.phi = phi % (2 * PI)
        self.sigma = sigma
        
    def sep(self, neighbors_phi):
        """
        Акт сепарации (бадаль).
        Смотрит на разность фаз с соседями.
        Возвращает бит различения.
        """
        if len(neighbors_phi) == 0:
            return 0
        
        # Средняя разность фаз с соседями
        delta = np.mean([self._phase_diff(self.phi, n) for n in neighbors_phi])
        
        # Порог различения: если delta > pi/2, различаем (бит=1)
        bit = 1 if abs(delta) > PI / 2 else 0
        
        # Обновляем sigma
        self.sigma = (self.sigma + bit) % 2
        
        return bit
    
    def phase_update(self, neighbors_phi, alpha=1/137):
        """
        Обновление фазы с учётом соседей.
        alpha — константа связи.
        """
        if len(neighbors_phi) == 0:
            return
        
        # Курамото-подобная динамика
        coupling = 0
        for n_phi in neighbors_phi:
            coupling += np.sin(n_phi - self.phi)
        coupling /= len(neighbors_phi)
        
        # Сдвиг фазы
        self.phi = (self.phi + alpha * coupling) % (2 * PI)
    
    def _phase_diff(self, phi1, phi2):
        """Разность фаз на окружности [-pi, pi]."""
        d = (phi2 - phi1) % (2 * PI)
        if d > PI:
            d -= 2 * PI
        return d


# ============================================================
# МЕМБРАНА (РЕШЁТКА НА T^2)
# ============================================================

class Membrane:
    """Решётка осцилляторов на торе T^2."""
    
    def __init__(self, N, alpha=1/137):
        """
        N: размер решётки (N x N)
        alpha: константа связи
        """
        self.N = N
        self.alpha = alpha
        self.grid = [[Oscillator() for _ in range(N)] for _ in range(N)]
        self.history = []  # история битов
        
    def get_neighbors_phi(self, i, j):
        """Фазы 4 соседей (периодические границы)."""
        N = self.N
        neighbors = [
            self.grid[(i-1) % N][j],
            self.grid[(i+1) % N][j],
            self.grid[i][(j-1) % N],
            self.grid[i][(j+1) % N],
        ]
        return [osc.phi for osc in neighbors]
    
    def step(self):
        """Один такт G: Sep + Phase для всех осцилляторов."""
        bits = []
        new_phis = [[0.0] * self.N for _ in range(self.N)]
        
        # Фаза 1: Sep (все одновременно)
        for i in range(self.N):
            for j in range(self.N):
                osc = self.grid[i][j]
                neighbors_phi = self.get_neighbors_phi(i, j)
                bit = osc.sep(neighbors_phi)
                bits.append(bit)
        
        # Фаза 2: Phase (все одновременно, используем старые фазы)
        for i in range(self.N):
            for j in range(self.N):
                osc = self.grid[i][j]
                neighbors_phi = self.get_neighbors_phi(i, j)
                osc.phase_update(neighbors_phi, self.alpha)
        
        self.history.append(sum(bits))
        return bits
    
    def get_phase_matrix(self):
        """Матрица фаз."""
        return np.array([[self.grid[i][j].phi for j in range(self.N)] for i in range(self.N)])
    
    def compute_spectral_defect(self):
        """
        Спектральный дефект: разность между
        континуальным (pi) и дискретным спектром.
        """
        phases = self.get_phase_matrix()
        
        # Дискретный лапласиан на торе
        laplacian = np.zeros_like(phases)
        N = self.N
        for i in range(N):
            for j in range(N):
                laplacian[i,j] = (
                    phases[(i+1)%N, j] + phases[(i-1)%N, j] +
                    phases[i, (j+1)%N] + phases[i, (j-1)%N] -
                    4 * phases[i, j]
                )
        
        # Спектр (собственные значения)
        eigenvalues = np.linalg.eigvalsh(laplacian)
        
        # Дефект: отклонение от идеального спектра
        # Идеальный спектр на континуальном торе: 4*pi^2*(n^2 + m^2)
        # Мы смотрим на gap между первым ненулевым и нулём
        sorted_eig = np.sort(np.abs(eigenvalues))
        gap = sorted_eig[1] if len(sorted_eig) > 1 else 0
        
        return gap


# ============================================================
# СИМУЛЯЦИЯ
# ============================================================

def run_simulation(N=8, alpha=1/137, steps=100):
    """Запуск симуляции."""
    print("=" * 60)
    print(f"OSCILLATOR SIMULATION: N={N}, alpha={alpha:.6f}, steps={steps}")
    print(f"pi - 3 = {PI_MINUS_3:.8f}")
    print("=" * 60)
    
    membrane = Membrane(N, alpha)
    
    defects = []
    for t in range(steps):
        membrane.step()
        if t % 10 == 0:
            defect = membrane.compute_spectral_defect()
            defects.append(defect)
            print(f"  t={t:4d}  defect={defect:.6f}")
    
    # Финальный анализ
    final_defect = defects[-1] if defects else 0
    
    # Ищем делитель, дающий (pi-3)
    best_div = 1
    best_match = 0
    for div in range(1, 1000):
        val = final_defect / div
        match = 1 - abs(val - PI_MINUS_3) / PI_MINUS_3
        if match > best_match:
            best_match = match
            best_div = div
    
    print()
    print(f"Final defect: {final_defect:.6f}")
    print(f"Best divisor: {best_div}")
    print(f"|Defect|/{best_div} = {final_defect/best_div:.6f}")
    print(f"Match to (pi-3): {best_match*100:.2f}%")
    
    # График
    plt.figure(figsize=(10, 4))
    
    plt.subplot(1, 2, 1)
    plt.plot(defects)
    plt.xlabel("Step / 10")
    plt.ylabel("Spectral Defect")
    plt.title("Defect Evolution")
    plt.axhline(y=PI_MINUS_3, color='r', linestyle='--', label=f'π-3={PI_MINUS_3:.4f}')
    plt.legend()
    
    plt.subplot(1, 2, 2)
    phases = membrane.get_phase_matrix()
    plt.imshow(phases, cmap='hsv', vmin=0, vmax=2*PI)
    plt.colorbar(label='Phase')
    plt.title(f"Phase Pattern (t={steps})")
    
    plt.tight_layout()
    plt.savefig(r"C:\Users\admin\Documents\ontological-resolution-theory\figures\oscillator_sim.png", dpi=150)
    print(f"\nSaved: oscillator_sim.png")
    
    return membrane, defects


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    # Базовый тест
    run_simulation(N=8, alpha=1/137, steps=100)
    
    print()
    print("=" * 60)
    print("ALPHA SCAN")
    print("=" * 60)
    
    # Сканируем alpha
    for alpha in [1/137, 1/100, 1/50, 1/20]:
        membrane = Membrane(N=8, alpha=alpha)
        for _ in range(100):
            membrane.step()
        defect = membrane.compute_spectral_defect()
        print(f"  alpha=1/{int(1/alpha):3d}  defect={defect:.6f}")
    
    print()
    print("=== SIMULATION COMPLETE ===")