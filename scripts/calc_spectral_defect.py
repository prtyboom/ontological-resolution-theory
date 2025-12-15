import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import os

os.makedirs(r'C:\Users\admin\Documents\ontological-resolution-theory\figures', exist_ok=True)

PI_MINUS_3 = np.pi - 3

def dirac_continuous(N_modes):
    k = np.arange(-N_modes, N_modes + 1)
    return np.sort((k + 0.5))

def dirac_discrete(N):
    a = 2 * np.pi / N
    shift = np.roll(np.eye(N), -1, axis=1)
    shift[-1, 0] = -1
    nabla = (shift - shift.T) / (2 * a)
    zero = np.zeros((N, N))
    D = np.block([[zero, nabla], [nabla.T, zero]])
    return D

def spectral_action(spectrum, Lambda):
    return np.sum(np.exp(-spectrum**2 / Lambda**2))

def compute(N, Lambda):
    D = dirac_discrete(N)
    eigs = np.sqrt(np.maximum(linalg.eigvalsh(D @ D), 0))
    spec_c = np.abs(dirac_continuous(N)[:2*N])
    S_c = spectral_action(spec_c, Lambda)
    S_d = spectral_action(eigs, Lambda)
    return S_c - S_d, S_c, S_d, eigs, spec_c

print("=" * 60)
print("OPERATION SPECTRAL GAP v3")
print("pi - 3 =", PI_MINUS_3)
print("=" * 60)

Lambda = 10.0
N_list = [3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]

print("\n  N       S_cont       S_disc       Defect    Def/N^2")
print("-" * 58)

results = []
for N in N_list:
    d, sc, sd, _, _ = compute(N, Lambda)
    results.append((N, d, sc, sd))
    print(f"{N:3d}  {sc:11.4f}  {sd:11.4f}  {d:11.4f}  {d/N**2:9.6f}")

print("\n" + "=" * 60)
print("FOCUS: N = 12")
print("=" * 60)

d12, sc12, sd12, spec_d, spec_c = compute(12, Lambda)

print(f"\nDefect at N=12:  {d12:.8f}")
print(f"pi - 3:          {PI_MINUS_3:.8f}")
print(f"Defect/(pi-3):   {d12/PI_MINUS_3:.4f}")
print(f"|Defect|/52:     {abs(d12)/52:.8f}")
print(f"Match to pi-3:   {100*(1 - abs(abs(d12)/52 - PI_MINUS_3)/PI_MINUS_3):.2f}%")

print("\n" + "=" * 60)
print("INTEGER SEARCH")
print("=" * 60)

ratio = abs(d12) / PI_MINUS_3
print(f"\n|Defect|/(pi-3) = {ratio:.4f} ~ {round(ratio)}")

n = round(ratio)
print(f"\nFactors of {n}:")
temp = n
for p in [2, 3, 5, 7, 11, 13]:
    while temp % p == 0:
        print(f"  {p}")
        temp //= p
if temp > 1:
    print(f"  {temp}")

fig, ax = plt.subplots(2, 2, figsize=(12, 10))

Ns = [r[0] for r in results]
defs = [r[1] for r in results]

ax[0,0].plot(Ns, defs, 'bo-', lw=2, ms=8)
ax[0,0].axhline(-52*PI_MINUS_3, color='r', ls='--', label='-52(pi-3)')
ax[0,0].axvline(12, color='gold', ls=':')
ax[0,0].set_xlabel('N')
ax[0,0].set_ylabel('Defect')
ax[0,0].set_title('Spectral Defect vs N')
ax[0,0].legend()
ax[0,0].grid(True, alpha=0.3)

ax[0,1].plot(Ns, [abs(d)/52 for d in defs], 'mo-', lw=2, ms=8)
ax[0,1].axhline(PI_MINUS_3, color='r', ls='--', label='pi-3')
ax[0,1].axvline(12, color='gold', ls=':')
ax[0,1].set_xlabel('N')
ax[0,1].set_ylabel('|Defect|/52')
ax[0,1].set_title('Normalized: |Defect|/52')
ax[0,1].legend()
ax[0,1].grid(True, alpha=0.3)

ax[1,0].hist(spec_d, bins=20, alpha=0.6, label='Discrete', density=True)
ax[1,0].hist(spec_c, bins=20, alpha=0.6, label='Continuous', density=True)
ax[1,0].set_xlabel('|lambda|')
ax[1,0].set_ylabel('Density')
ax[1,0].set_title('Spectra at N=12')
ax[1,0].legend()
ax[1,0].grid(True, alpha=0.3)

ax[1,1].plot(Ns, [abs(d)/PI_MINUS_3 for d in defs], 'co-', lw=2, ms=8)
ax[1,1].axhline(52, color='r', ls='--', label='52')
ax[1,1].axvline(12, color='gold', ls=':')
ax[1,1].set_xlabel('N')
ax[1,1].set_ylabel('|Defect|/(pi-3)')
ax[1,1].set_title('Ratio to (pi-3)')
ax[1,1].legend()
ax[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(r'C:\Users\admin\Documents\ontological-resolution-theory\figures\spectral_v3.png', dpi=150)
print("\nSaved: spectral_v3.png")
plt.show()

print("\n" + "=" * 60)
print("DONE")
print("=" * 60)