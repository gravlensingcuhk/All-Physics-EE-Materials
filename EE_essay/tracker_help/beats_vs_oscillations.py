"""Illustration for the Tracker workflow: what a 'beat' is, and why a 6 s record cannot resolve f1, f2.
Parameters: d = 21 mm (f1 = 0.823, f2 = 1.063 Hz), hold-and-release with A = 10 mm, tau = 170 s."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
f1, f2, A, tau, fs = 0.823, 1.063, 10.0, 170.0, 240.0
gk = ((f2/f1)**2 - 1)/2; x20 = A*gk/(1+gk)            # spring 2 shift while spring 1 is held
eta0, xi0 = A + x20, A - x20; a, b = eta0/2, xi0/2
t = np.arange(0, 30, 1/fs); w1, w2 = 2*np.pi*f1, 2*np.pi*f2
x = (a*np.cos(w1*t) + b*np.cos(w2*t))*np.exp(-t/tau)
env = np.sqrt(a*a + b*b + 2*a*b*np.cos((w2-w1)*t))*np.exp(-t/tau)
Tbeat, Tavg = 1/(f2-f1), 2/(f1+f2)

def spec(sig, fs):
    sig = sig - sig.mean(); win = np.hanning(len(sig)); n = 2**18
    X = np.abs(np.fft.rfft(sig*win, n)); f = np.fft.rfftfreq(n, 1/fs); return f, X/X.max()

C = "#b35900"
fig = plt.figure(figsize=(11, 7.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.45, wspace=0.25)
ax = fig.add_subplot(gs[0, :])
ax.plot(t, x, color=C, lw=1.0); ax.plot(t, env, "--", color="crimson", lw=1); ax.plot(t, -env, "--", color="crimson", lw=1)
ax.axvspan(0, 6, color="gray", alpha=0.15); ax.axhline(0, color="0.6", lw=0.5)
ax.set_xlim(0, 30); ax.set_ylim(-15.5, 16.5); ax.set_xlabel("t / s"); ax.set_ylabel("x$_1$ / mm")
ax.set_title("d = 21 mm:  x$_1$(t) after releasing spring 1   (f$_1$ = 0.823 Hz, f$_2$ = 1.063 Hz)", loc="left", fontsize=10)
# one oscillation
t0 = 0.0; t1 = Tavg
t0 = 2*Tavg; t1 = 3*Tavg
ax.annotate("", (t0, 11.2), (t1, 11.2), arrowprops=dict(arrowstyle="<->", color="k"))
ax.text((t0+t1)/2, 11.7, "1 oscillation (T$_{avg}$ = 1.06 s)", ha="center", va="bottom", fontsize=8)
# one beat: between two envelope maxima (t = 0 and Tbeat)
ax.annotate("", (Tbeat, -11.8), (2*Tbeat, -11.8), arrowprops=dict(arrowstyle="<->", color="crimson"))
ax.text(1.5*Tbeat, -12.4, "1 beat = 1 swelling (T$_{beat}$ = 1/(f$_2$ - f$_1$) = 4.2 s)", ha="center", va="top", fontsize=8, color="crimson")
ax.text(6.3, 15.9, "shaded = first 6 s: 5.7 oscillations (11 turning points) but only 1.4 beats", ha="left", va="top", fontsize=8.5,
        bbox=dict(fc="white", ec="gray", lw=0.5))
ax.text(29.7, -13.3, "30 s: 28 oscillations, 7.2 beats", ha="right", va="bottom", fontsize=8.5, bbox=dict(fc="white", ec="gray", lw=0.5))
for k, (T, lab) in enumerate([(6, "FFT of the first 6 s"), (30, "FFT of the full 30 s")]):
    ax = fig.add_subplot(gs[1, k]); m = t < T; f, X = spec(x[m], fs)
    ax.plot(f, X, color="k", lw=1); ax.set_xlim(0.5, 1.4); ax.set_ylim(0, 1.15)
    for fk in (f1, f2): ax.axvline(fk, color="crimson", ls=":", lw=1)
    ax.set_xlabel("frequency / Hz"); ax.set_ylabel("|FFT| (normalised)")
    ax.set_title(f"{lab}:  bin width 1/T = {1/T:.3f} Hz,  f$_2$ - f$_1$ = {f2-f1:.2f} Hz = {(f2-f1)*T:.1f} bins", fontsize=9, loc="left")
    ax.text(0.98, 0.95, "one merged blob:\nf$_1$, f$_2$ NOT resolved" if T == 6 else "two clean peaks\nat f$_1$ and f$_2$",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color="crimson" if T == 6 else "green")
fig.savefig("beats_vs_oscillations.png", dpi=170, bbox_inches="tight")
print(f"x20 = {x20:.2f} mm, Tbeat = {Tbeat:.2f} s, Tavg = {Tavg:.3f} s, osc in 6 s = {6/Tavg:.1f}, beats in 6 s = {6/Tbeat:.2f}")
