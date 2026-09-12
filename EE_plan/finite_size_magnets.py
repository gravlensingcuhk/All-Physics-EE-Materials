"""
Part D - finite-size magnet check ("magnetic charge" / Gilbert model).
Two identical coaxial cylindrical magnets (diameter D, thickness t), like poles facing,
centre-to-centre separation d. Each magnet is replaced by two uniformly charged discs
(+ on one face, - on the other). The axial force is obtained by numerical integration
and gamma(d) = -dF/dd is compared with the point-dipole law gamma ∝ d^-5 over YOUR
range d = 21..35 mm.  Replace D and t with your magnet's dimensions.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
OUT = Path(__file__).parent / "figs"; OUT.mkdir(exist_ok=True)

def disc_disc_force(z, R, nr=24, nphi=36):
    r = (np.arange(nr) + 0.5) / nr * R
    phi = (np.arange(nphi) + 0.5) / nphi * 2 * np.pi
    dA = (R / nr) * (2 * np.pi / nphi) * r
    rr, pp = np.meshgrid(r, phi, indexing="ij")
    x1, y1 = rr * np.cos(pp), rr * np.sin(pp)
    dAA = np.broadcast_to(dA[:, None], rr.shape)
    F = 0.0
    for i in range(nr):                      # disc 2 is axisymmetric -> one azimuth suffices
        x2 = r[i]
        dist2 = (x1 - x2) ** 2 + y1 ** 2 + z ** 2
        F += np.sum(dAA * (dA[i] * nphi) * z / dist2 ** 1.5)
    return F

def magnet_force(d, D, t):
    R = D / 2
    return (disc_disc_force(d - t, R) + disc_disc_force(d + t, R) - 2 * disc_disc_force(d, R))

def gamma(d, D, t, h=0.05):
    return -(magnet_force(d + h, D, t) - magnet_force(d - h, D, t)) / (2 * h)

d = np.array([21, 23, 25, 28, 31, 35.0])
geoms = [(6, 3), (10, 3), (10, 5), (12, 6), (15, 5), (20, 3), (20, 5), (25, 5)]
print("apparent gradient of ln(gamma) vs ln(d) over d = 21..35 mm")
res = []
for D, t in geoms:
    g = np.array([gamma(x, D, t) for x in d])
    slope = np.polyfit(np.log(d), np.log(g), 1)[0]
    res.append((D, t, slope))
    print(f"  D = {D:2d} mm, t = {t:2d} mm  ->  {slope:.2f}")

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axs[0]
for D, t, c in [(6, 3, "tab:green"), (12, 6, "tab:orange"), (20, 5, "tab:red")]:
    g = np.array([gamma(x, D, t) for x in d]); g /= g[2]
    ax.plot(np.log(d), np.log(g), "o-", color=c, label=f"discs D = {D} mm, t = {t} mm  (gradient {np.polyfit(np.log(d), np.log(g), 1)[0]:.2f})")
ax.plot(np.log(d), -5 * (np.log(d) - np.log(25)), "k--", label="point dipole (gradient −5)")
ax.set_xlabel("ln(d / mm)"); ax.set_ylabel("ln(γ / γ(25 mm))"); ax.legend(fontsize=8)
ax.set_title("Finite-size magnets flatten the ln–ln gradient", fontsize=10)
ax = axs[1]
Ds = np.arange(4, 27, 1.5)
for t, c in [(3, "tab:blue"), (5, "tab:red"), (8, "tab:green")]:
    sl = []
    for D in Ds:
        g = np.array([gamma(x, D, t) for x in d]); sl.append(np.polyfit(np.log(d), np.log(g), 1)[0])
    ax.plot(Ds, sl, "o-", color=c, ms=4, label=f"thickness t = {t} mm")
ax.axhline(-5, color="k", ls="--", lw=1, label="point dipole")
ax.axhline(-3.75, color="tab:purple", ls=":", lw=1.5, label="your measured gradient (−3.75)")
ax.set_xlabel("magnet diameter D / mm"); ax.set_ylabel("apparent gradient over 21–35 mm"); ax.legend(fontsize=8)
ax.set_title("Apparent gradient vs magnet size (Gilbert disc model)", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / "D_finite_size_gradient.png", dpi=130)
print("saved", OUT / "D_finite_size_gradient.png")
