"""Design plot for the experimental plan (REAL time base: 240 fps clip read at 240 fps).
Constants from run 1 after the x8 time-base correction: f1 = 6.49 Hz, C = 1.33e-7 m^5 s^-2, tau ~ 21 s."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
f1, C, tau = 6.49, 1.33e-7, 21.0
d = np.linspace(23.5, 46, 300)/1000
f2 = np.sqrt(f1**2 + C/d**5); Tb = 1/(f2 - f1)
dsel = np.array([24.5, 26.5, 28.5, 31, 33.5, 36, 39, 42])/1000
Trec = np.array([45, 45, 45, 45, 45, 45, 45, 45])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
ax[0].plot(d*1e3, f2, color="crimson", lw=1.5, label="$f_2$ (anti-phase)")
ax[0].axhline(f1, color="#b35900", lw=1.5, label="$f_1$ (in-phase)")
ax[0].plot(dsel*1e3, np.sqrt(f1**2 + C/dsel**5), "o", color="crimson", ms=5)
ax[0].plot(dsel*1e3, np.full(len(dsel), f1), "o", color="#b35900", ms=5)
ax[0].set_xlabel("centre-to-centre separation $d$ / mm"); ax[0].set_ylabel("frequency / Hz")
ax[0].set_title("(a) predicted mode frequencies; dots = chosen separations", fontsize=10, loc="left")
ax[0].legend(frameon=False); ax[0].set_ylim(6.0, 8.6); ax[0].grid(alpha=0.3)
ax[1].semilogy(d*1e3, Tb, color="k", lw=1.5, label="beat period $T_{beat}=1/(f_2-f_1)$")
ax[1].semilogy(d*1e3, 3*Tb, "--", color="k", lw=1, label="$3\\,T_{beat}$ (minimum record)")
ax[1].step(np.r_[dsel*1e3, 46], np.r_[Trec, Trec[-1]], where="post", color="tab:blue", lw=2, label="recommended record: 45 s for every $d$")
ax[1].plot(dsel*1e3, Trec, "s", color="tab:blue", ms=5)
for x, T in zip(dsel*1e3, Trec): ax[1].annotate(f"{T:.0f} s", (x, T), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=7.5, color="tab:blue")
ax[1].axhline(6, color="gray", lw=0.8, ls=":"); ax[1].text(23.8, 6.5, "your 6 s record", fontsize=8, color="gray")
ax[1].axhline(2.5*tau, color="crimson", lw=0.8, ls=":"); ax[1].text(23.8, 2.5*tau*1.1, "2.5$\\tau$: amplitude down to 8 % (no point recording longer)", fontsize=8, color="crimson")
ax[1].set_xlabel("centre-to-centre separation $d$ / mm"); ax[1].set_ylabel("time / s")
ax[1].set_title("(b) how long to record", fontsize=10, loc="left"); ax[1].legend(frameon=False, fontsize=8, loc="lower right")
ax[1].set_ylim(0.5, 120); ax[1].grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("design_plot.pdf"); fig.savefig("design_plot.png", dpi=160)
print("ok")
