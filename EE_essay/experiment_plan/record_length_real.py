"""Why 6 s is enough for the beats but not for the peak positions, and why 50 s is needed at 42 mm.
Real time base: f1 = 6.49 Hz, C = 1.33e-7, tau = 21 s, hold-and-release A = 5 mm, fs = 240."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
f1, C, tau, A, fs = 6.49, 1.33e-7, 21.0, 5.0, 240.0
def signal(d, Tmax):
    f2 = np.sqrt(f1**2 + C/d**5); gk = (f2**2/f1**2 - 1)/2; x20 = A*gk/(1 + gk)
    a, b = (A + x20)/2, (A - x20)/2; t = np.arange(0, Tmax, 1/fs)
    return t, (a*np.cos(2*np.pi*f1*t) + b*np.cos(2*np.pi*f2*t))*np.exp(-t/tau), f2
def spec(x):
    x = x - x.mean(); X = np.abs(np.fft.rfft(x*np.hanning(len(x)), 2**18)); f = np.fft.rfftfreq(2**18, 1/fs); return f, X/X.max()
C1, C2 = "#b35900", "crimson"
fig = plt.figure(figsize=(11, 6.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], hspace=0.5, wspace=0.22)
ax = fig.add_subplot(gs[0, :]); t, x, f2 = signal(0.0245, 6.0); Tb = 1/(f2 - f1)
ax.plot(t, x, color=C1, lw=0.9); ax.axhline(0, color="0.6", lw=0.5)
for k in range(1, 7): ax.axvline(k*Tb, color=C2, ls=":", lw=0.8)
ax.annotate("", (Tb, 4.3), (2*Tb, 4.3), arrowprops=dict(arrowstyle="<->", color=C2)); ax.text(1.5*Tb, 4.5, f"$T_{{beat}}$ = {Tb:.2f} s", ha="center", va="bottom", fontsize=8.5, color=C2)
ax.set_xlim(0, 6); ax.set_ylim(-5.5, 6.2); ax.set_xlabel("t / s"); ax.set_ylabel("x$_1$ / mm")
ax.set_title(f"d = 24.5 mm, real time base: f$_1$ = {f1:.2f} Hz, f$_2$ = {f2:.2f} Hz  —  6 s contain {6/Tb:.1f} beats and {6*(f1+f2)/2:.0f} oscillations", loc="left", fontsize=10)
for k, (dmm, Ts) in enumerate([(24.5, [6, 30]), (42, [6, 20, 50])]):
    ax = fig.add_subplot(gs[1, k]); d = dmm/1000
    for T, col, ls in zip(Ts, ["0.55", "k", "tab:blue"], ["-", "-", "-"]):
        t, x, f2 = signal(d, T); f, X = spec(x); ax.plot(f, X, color=col, lw=1.1, ls=ls, label=f"T = {T} s  (bin 1/T = {1/T:.3f} Hz)")
    for fk in (f1, f2): ax.axvline(fk, color=C2, ls=":", lw=1)
    ax.set_xlim(5.8, 8.3 if dmm < 30 else 7.2); ax.set_ylim(0, 1.5); ax.set_xlabel("frequency / Hz"); ax.set_ylabel("|FFT| (normalised, Hann)")
    ax.set_title(f"({'bc'[k]}) d = {dmm} mm:  f$_2$ - f$_1$ = {f2-f1:.3f} Hz,  T$_{{beat}}$ = {1/(f2-f1):.1f} s", fontsize=9.5, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper left" if dmm > 30 else "upper right")
fig.savefig("record_length_real.png", dpi=170, bbox_inches="tight"); print("ok")
