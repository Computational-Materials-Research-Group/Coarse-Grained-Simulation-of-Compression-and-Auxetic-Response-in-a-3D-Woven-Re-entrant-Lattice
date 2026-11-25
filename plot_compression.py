"""
plot_compression.py - stress-strain, lateral strain and Poisson's ratio of the
woven re-entrant lattice from stress_strain.txt (written by in.compression_lattice).

    python plot_compression.py                     # -> compression_results.png
    python plot_compression.py run1/stress_strain.txt run2/stress_strain.txt   # overlay runs
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

files = sys.argv[1:] or ["stress_strain.txt"]
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for f in files:
    d = np.loadtxt(f, comments="#")
    eps, sig, lat, nb = d[:, 1], d[:, 2], d[:, 4], d[:, 6]
    w = max(1, len(eps) // 40)                    # moving average: the dynamic stress is noisy
    sig = np.convolve(sig, np.ones(w) / w, mode="same")
    lab = f if len(files) > 1 else None
    ax[0, 0].plot(100 * eps, sig, lw=1.6, label=lab)
    ax[0, 1].plot(100 * eps, 100 * lat, lw=1.6, label=lab)
    ok = eps > 0.01                               # Poisson ratio is noisy at tiny strain
    ax[1, 0].plot(100 * eps[ok], -lat[ok] / -eps[ok], lw=1.6, label=lab)
    ax[1, 1].plot(100 * eps, 100 * nb / nb[0], lw=1.6, label=lab)
    # summary
    print(f"{f}: max strain {100*eps.max():.1f} %, peak stress {sig.max():.3g}, "
          f"lateral strain at end {100*lat[-1]:.2f} %, "
          f"Poisson ratio at end {(-lat[-1] / -eps[-1]):.3f}, bonds intact {100*nb[-1]/nb[0]:.1f} %")

ax[0, 0].set(xlabel="compressive strain (%)", ylabel="compressive stress (lj)", title="Stress-strain (moving average)")
ax[0, 1].set(xlabel="compressive strain (%)", ylabel="lateral strain $\\epsilon_x$ (%)",
             title="Lateral strain at mid-height (< 0: shrinks = auxetic)")
ax[0, 1].axhline(0, color="k", lw=0.6)
ax[1, 0].set(xlabel="compressive strain (%)", ylabel="$\\nu = -\\epsilon_x / \\epsilon_z$",
             title="Poisson's ratio")
ax[1, 0].axhline(0, color="k", lw=0.6)
ax[1, 1].set(xlabel="compressive strain (%)", ylabel="intact bonds (%)", title="Damage")
for a in ax.flat:
    a.grid(alpha=0.3)
    if len(files) > 1:
        a.legend(fontsize=8)
plt.tight_layout()
plt.savefig("compression_results.png", dpi=150)
print("wrote compression_results.png")
