#!/usr/bin/env python3
"""
make_fig_release.py  --  Figure 13 of the theory section ("fig_release.pdf").

Linear two-mode model of the coupled oscillator, released from rest:

    x1(t) = (eta0/2) cos(w1 t) + (xi0/2) cos(w2 t),
    eta0 = x1(0) + x2(0),   xi0 = x1(0) - x2(0).

Row (a): ideal instantaneous release, x1(0) = A, x2(0) = 0
         -> eta0 = xi0 = A, equal FFT peaks, beat minima at zero.
Row (b): hold-and-release, spring 2 has already shifted to
         x2(0) = A*gamma/(k+gamma) while spring 1 is held
         -> peak ratio (A-x20)/(A+x20) = k/(k+2 gamma) = (f1/f2)^2.

Left column : x1(t) with its envelope.   Right column: |FFT| of x1(t).

Usage:  python3 make_fig_release.py          (writes fig_release.pdf and .png)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- parameters
f1, f2 = 0.802, 0.996        # Hz, the two mode frequencies measured at d = 23 mm
A = 8.0                      # mm, release displacement of spring 1
gk = ((f2 / f1) ** 2 - 1) / 2            # gamma/k  from (f2/f1)^2 = 1 + 2 gamma/k
x20_hold = A * gk / (1 + gk)             # shift of spring 2 while spring 1 is held
fs, T = 240.0, 120.0                     # sampling rate (fps) and record length for the FFT
t_show = 16.0                            # seconds of time series displayed

COPPER = "#b35900"   # same hue as the TikZ copper (orange!80!black)
ENV = "crimson"

t = np.arange(0.0, T, 1.0 / fs)
w1, w2 = 2 * np.pi * f1, 2 * np.pi * f2


def motion(x10, x20):
    """x1(t), its envelope, and the two mode amplitudes for a release from rest."""
    eta0, xi0 = x10 + x20, x10 - x20
    a, b = eta0 / 2, xi0 / 2
    x1 = a * np.cos(w1 * t) + b * np.cos(w2 * t)
    env = np.sqrt(a**2 + b**2 + 2 * a * b * np.cos((w2 - w1) * t))   # |a e^{iw1t} + b e^{iw2t}|
    return x1, env, eta0, xi0


def spectrum(x):
    """Amplitude spectrum with a Hann window, normalised so a cosine of amplitude a gives a peak of height a."""
    win = np.hanning(len(x))
    nfft = 8 * len(x)                       # zero-padding: removes the scalloping loss, so the
    X = np.fft.rfft(x * win, n=nfft)        # peak height equals the cosine amplitude to < 1 %
    f = np.fft.rfftfreq(nfft, 1.0 / fs)
    return f, 2.0 * np.abs(X) / win.sum()


cases = [
    ("(a) ideal release: $x_2(0)=0$", A, 0.0),
    (f"(b) hold-and-release: $x_2(0)={x20_hold:.1f}$ mm", A, x20_hold),
]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5})
fig, axes = plt.subplots(2, 2, figsize=(11.0, 5.4), gridspec_kw={"width_ratios": [1.6, 1.0]})

for row, (title, x10, x20) in enumerate(cases):
    x1, env, eta0, xi0 = motion(x10, x20)
    f, amp = spectrum(x1)

    # ---- time series
    ax = axes[row, 0]
    m = t <= t_show
    ax.plot(t[m], x1[m], color=COPPER, lw=1.0, label="$x_1(t)$")
    ax.plot(t[m], env[m], "--", color=ENV, lw=1.0, label="envelope")
    ax.plot(t[m], -env[m], "--", color=ENV, lw=1.0)
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlim(0, t_show)
    ax.set_ylim(-1.3 * A, 1.3 * A)
    ax.set_xlabel("$t$ / s")
    ax.set_ylabel("$x_1$ / mm")
    ax.set_title(title, loc="left")
    ax.text(0.99, 0.95, f"max $=(\\eta_0+\\xi_0)/2={(eta0+xi0)/2:.1f}$ mm",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=ENV)
    ax.text(0.99, 0.05, f"min $=|\\eta_0-\\xi_0|/2={abs(eta0-xi0)/2:.1f}$ mm",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=ENV)
    if row == 0:
        ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)

    # ---- spectrum
    ax = axes[row, 1]
    band = (f > 0.6) & (f < 1.2)
    ax.plot(f[band], amp[band], color="k", lw=1.0)
    for fk, hk, lab in [(f1, eta0 / 2, r"$f_1$: height $\eta_0/2$"), (f2, xi0 / 2, r"$f_2$: height $\xi_0/2$")]:
        ax.plot([fk - 0.03, fk + 0.03], [hk, hk], color=COPPER, lw=1.2)   # predicted height
        ax.annotate(lab, (fk, hk), xytext=(0, 5), textcoords="offset points",
                    ha="center", fontsize=8)
    ax.set_xlim(0.6, 1.2)
    ax.set_ylim(0, 1.35 * A / 2)
    ax.set_xlabel("frequency / Hz")
    ax.set_ylabel("|FFT| / mm")
    ratio = (xi0 / 2) / (eta0 / 2)
    ax.text(0.97, 0.93, f"peak ratio $= {ratio:.2f}$" + ("" if row == 0 else "\n$=(f_1/f_2)^2$"),
            transform=ax.transAxes, ha="right", va="top", fontsize=8)

fig.suptitle(f"Same $f_1={f1}$ Hz, $f_2={f2}$ Hz in both rows; only the starting condition differs "
             f"(linear model, $A={A:.0f}$ mm)", fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig("fig_release.pdf")
fig.savefig("fig_release.png", dpi=180)
print(f"gamma/k = {gk:.3f}, x20 = {x20_hold:.2f} mm, predicted peak ratio = {(f1/f2)**2:.3f}")
