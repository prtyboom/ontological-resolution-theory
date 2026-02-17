import math

r =   [1.0, 1.414, 2.0, 2.828, 4.0, 5.0, 6.0, 7.0]
phi = [0.01714, 0.01222, 0.00784, 0.00502, 0.00295, 0.00196, 0.00128, 0.00082]

print("  r1 -> r2       n = -dlog(phi)/dlog(r)")
print("  " + "-" * 42)
for i in range(len(r) - 1):
    n = -(math.log(phi[i+1]) - math.log(phi[i])) / (math.log(r[i+1]) - math.log(r[i]))
    print(f"  {r[i]:.3f} -> {r[i+1]:.3f}    n = {n:.4f}")

print(f"\n  Expected for 1/r: n = 1.0000")