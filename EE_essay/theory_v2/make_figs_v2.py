"""Figures for theory_v2.tex (real time base: f1 = 6.49 Hz, C = 1.33e-7 m^5 s^-2, tau_d = 21 s).
fig_prediction.pdf  - f1, f2 against d and the log-log line of gradient -5
fig_release.pdf     - hold-and-release at d = 24.5 mm, A = 5 mm, with Kelvin-Voigt damping: traces + spectra
fig_nonlinear.pdf   - amplitude dependence of f2 (numerical vs perturbation) and the harmonics of xi
fig_damping.pdf     - how the beat minima evolve when the two modes decay at equal / unequal rates
"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
C1, C2, CG = "#b35900", "crimson", "0.45"
f1, C, tau = 6.49, 1.33e-7, 21.0
w1 = 2*np.pi*f1
def modes(d):
    gam_m = 2*np.pi**2*C/d**5; w2 = np.sqrt(w1**2 + 2*gam_m); return gam_m, w2
def spec(x, fs, n=2**20):
    x = x - x.mean(); X = np.abs(np.fft.rfft(x*np.hanning(len(x)), n))*2/np.sum(np.hanning(len(x)))
    return np.fft.rfftfreq(n, 1/fs), X

# ------------------------------------------------------------------ 1. prediction
d = np.linspace(20, 50, 400)/1000; f2 = np.sqrt(f1**2 + C/d**5)
dsel = np.array([24.5, 26.5, 28.5, 31, 33.5, 36, 39, 42])/1000
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(d*1e3, f2, color=C2, lw=1.6, label="$f_2$ (anti-phase)"); ax[0].axhline(f1, color=C1, lw=1.6, label="$f_1$ (in-phase)")
ax[0].plot(dsel*1e3, np.sqrt(f1**2 + C/dsel**5), "o", color=C2, ms=4); ax[0].plot(dsel*1e3, np.full(8, f1), "o", color=C1, ms=4)
# full non-linear force law: frequency of the relative motion for a hold-and-release start of amplitude A
def f2_nonlinear(dv, A):
    gam_m, w2 = modes(dv); xi0 = A/(1 + gam_m/w1**2); Fd_m = gam_m*dv/4
    rhs = lambda tt, y: [y[1], -w1**2*y[0] + 2*Fd_m*((1 + y[0]/dv)**-4 - 1)]
    T = 25/(w2/2/np.pi); tt = np.linspace(0, T, int(T*3000))
    sol = solve_ivp(rhs, (0, T), [xi0, 0], t_eval=tt, rtol=1e-10, atol=1e-13)
    pk, _ = find_peaks(sol.y[0]); return 1/np.mean(np.diff(tt[pk]))
dnl = np.linspace(21, 50, 16)/1000
for A, ls in [(5e-3, "--"), (10e-3, ":")]:
    fnl = np.array([f2_nonlinear(dv, A) for dv in dnl])
    ax[0].plot(dnl*1e3, fnl, ls, color=C2, lw=1.2, label=f"$f_2$, full force law, $A={A*1e3:.0f}$ mm")
    ax[1].loglog(dnl*1e3, fnl**2 - f1**2, ls, color="k", lw=1.2)
ax[0].axvline(24.5, color="0.7", lw=0.8, ls=":"); ax[0].text(24.9, 8.75, "smallest $d$\n(clamps)", fontsize=8, color="0.4", va="top")
ax[0].set_xlabel("centre-to-centre separation $d$ / mm"); ax[0].set_ylabel("frequency / Hz"); ax[0].set_ylim(6.0, 8.9)
ax[0].legend(frameon=False, loc="upper right", fontsize=8); ax[0].grid(alpha=0.3); ax[0].set_title("(a) mode frequencies", loc="left", fontsize=10)
ax[1].loglog(d*1e3, f2**2 - f1**2, color="k", lw=1.6, label="$f_2^2-f_1^2=C\\,d^{-5}$")
ax[1].loglog(dsel*1e3, C/dsel**5, "o", color="k", ms=4)
ax[1].set_xlabel("$d$ / mm (log scale)"); ax[1].set_ylabel("$f_2^2-f_1^2$ / Hz$^2$ (log scale)"); ax[1].grid(alpha=0.3, which="both")
ax[1].set_xticks([20, 25, 30, 40, 50]); ax[1].set_xticklabels(["20", "25", "30", "40", "50"])
ax[1].text(33, 4, "gradient $-5$", fontsize=9, rotation=-38); ax[1].legend(frameon=False)
ax[1].text(20.5, 0.55, "dashed / dotted: full force law,\nrelease amplitude $A=5$ / $10$ mm", fontsize=8, color="0.3")
ax[1].set_title("(b) the straight-line test", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("fig_prediction.pdf")

# ------------------------------------------------------------------ 2. hold-and-release with KV damping
d0 = 0.0245; gam_m, w2 = modes(d0); gk = gam_m/w1**2; A = 5.0
x20 = A*gk/(1 + gk); eta0, xi0 = -(A + x20), A - x20
fs = 240.0; t = np.arange(0, 45, 1/fs); dec = np.exp(-t/tau)
eta = eta0*np.cos(w1*t)*dec; xi = xi0*np.cos(w2*t)*dec
x1, x2 = (eta - xi)/2, (eta + xi)/2
# |x1| = (1/2)|eta0 cos w1t - xi0 cos w2t|: envelope (1/2)sqrt(eta0^2 + xi0^2 - 2 eta0 xi0 cos(dw t)); eta0 < 0 so this starts at A
env1 = 0.5*np.sqrt(eta0**2 + xi0**2 - 2*eta0*xi0*np.cos((w2 - w1)*t))*dec
env2 = 0.5*np.sqrt(eta0**2 + xi0**2 + 2*eta0*xi0*np.cos((w2 - w1)*t))*dec
fig, ax = plt.subplots(2, 2, figsize=(10.5, 5.6), gridspec_kw=dict(width_ratios=[1.7, 1]))
Tb = 2*np.pi/(w2 - w1)
for i, (x, env, lab, col) in enumerate([(x1, env1, "$x_1$ (released spring)", C1), (x2, env2, "$x_2$", "tab:blue")]):
    a = ax[i, 0]; m = t < 3.0
    a.plot(t[m], x[m], color=col, lw=0.9); a.plot(t[m], env[m], "--", color=C2, lw=1); a.plot(t[m], -env[m], "--", color=C2, lw=1)
    a.axhline(0, color="0.6", lw=0.5); a.set_ylim(-5.6, 5.6); a.set_ylabel(lab + " / mm")
    if i == 0:
        a.annotate("", (Tb, 4.6), (2*Tb, 4.6), arrowprops=dict(arrowstyle="<->", color=C2))
        a.text(1.5*Tb, 4.75, "$T_{beat}=1/(f_2-f_1)$", ha="center", va="bottom", fontsize=8.5, color=C2)
        a.annotate("", (1.5*Tb, -0.9), (1.5*Tb, 0.9), arrowprops=dict(arrowstyle="<->", color="k", lw=0.8))
        a.text(1.5*Tb + 0.04, -1.9, "beat minimum $=x_{20}=A\\gamma/(k+\\gamma)$", fontsize=8, va="center")
    if i == 1: a.set_xlabel("$t$ / s")
    fr, X = spec(x, fs); b = ax[i, 1]; mm = (fr > 5.5) & (fr < 8.6); X = X/X[mm].max()
    b.plot(fr[mm], X[mm], color=col, lw=1.1)
    r = (w1/w2)**2
    b.plot([w2/2/np.pi - 0.15, w2/2/np.pi + 0.15], [r, r], color=C2, lw=1.2)
    b.text(w2/2/np.pi + 0.18, r, "predicted ratio\n$\\xi_0/|\\eta_0|=(f_1/f_2)^2=%.2f$" % r, fontsize=8, va="center", color=C2)
    b.set_xlim(5.5, 8.6); b.set_ylim(0, 1.15); b.set_ylabel("|FFT| (normalised)"); b.grid(alpha=0.3)
    if i == 1: b.set_xlabel("frequency / Hz")
ax[0, 0].set_title(f"(a) $d=24.5$ mm, $A=5$ mm, $\\tau_d=21$ s:  $f_1={f1:.2f}$, $f_2={w2/2/np.pi:.2f}$ Hz", loc="left", fontsize=10)
ax[0, 1].set_title("(b) spectra of the 45 s records", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("fig_release.pdf")

# ------------------------------------------------------------------ 3. non-linear relative coordinate
def xi_num(d, xi0, T, tau_d=None):
    gam_m, w2 = modes(d); Fd_m = gam_m*d/4; lam = 0 if tau_d is None else 1/tau_d
    rhs = lambda tt, y: [y[1], -w1**2*y[0] + 2*Fd_m*((1 + y[0]/d)**-4 - 1) - 2*lam*y[1]]
    tt = np.linspace(0, T, int(T*3000)); sol = solve_ivp(rhs, (0, T), [xi0, 0], t_eval=tt, rtol=1e-11, atol=1e-14)
    return tt, sol.y[0]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
for dmm, col, mk in [(24.5, C2, "o"), (31.0, "tab:blue", "s")]:
    dd = dmm/1000; gam_m, w2 = modes(dd); g = gam_m/w2**2; f2l = w2/2/np.pi
    xs = np.linspace(0.0, 9, 60)/1000
    ax[0].plot(xs*1e3, 100*(xs/dd)**2*(15*g/4 - 125*g**2/12), color=col, lw=1.4, label=f"$d={dmm}$ mm, perturbation theory")
    pts = []
    for x0 in [1, 2, 3, 4, 5, 6, 7, 8]:
        tt, xx = xi_num(dd, x0/1000, 30/f2l); pk, _ = find_peaks(xx); fnum = 1/np.mean(np.diff(tt[pk]))
        w = 2*np.pi*fnum; M = np.column_stack([np.ones_like(tt), np.cos(w*tt), np.sin(w*tt)]); c, *_ = np.linalg.lstsq(M, xx, rcond=None)
        pts.append((np.hypot(c[1], c[2])*1e3, 100*(fnum/f2l - 1)))
    pts = np.array(pts); ax[0].plot(pts[:, 0], pts[:, 1], mk, color=col, ms=5, mfc="white", label=f"$d={dmm}$ mm, numerical")
ax[0].axvline(4.24, color="0.6", lw=0.8, ls=":"); ax[0].text(4.35, 3.45, "$\\xi_0$ for $A=5$ mm at 24.5 mm", fontsize=8, color="0.4", va="top", rotation=90)
ax[0].set_xlabel("amplitude $\\xi_0$ of the separation swing / mm"); ax[0].set_ylabel("shift of $f_2$ / %"); ax[0].set_ylim(0, 3.6)
ax[0].legend(frameon=False, fontsize=8); ax[0].grid(alpha=0.3); ax[0].set_title("(a) amplitude dependence of $f_2$", loc="left", fontsize=10)
dd = 0.0245; gam_m, w2 = modes(dd); g = gam_m/w2**2; f2l = w2/2/np.pi
tt, xx = xi_num(dd, 4.24e-3, 45); fsn = 1/(tt[1] - tt[0]); fr, X = spec(xx*1e3, fsn); mm = (fr > 0.5) & (fr < 26)
ax[1].semilogy(fr[mm], X[mm]/X[mm].max(), color="k", lw=1)
a0 = 4.24e-3
for n, rel, lab in [(1, 1, "$f_2$"), (2, 5*g*a0/(6*dd), "$2f_2$: $5g\\xi_0/6d$"), (3, (a0/dd)**2*(5*g/16 + 25*g**2/48), "$3f_2$: $(\\xi_0/d)^2(5g/16+25g^2/48)$")]:
    ax[1].plot([n*f2l - 0.6, n*f2l + 0.6], [rel, rel], color=C2, lw=1.2); ax[1].text(n*f2l + 0.8, rel*1.3, lab, fontsize=8, color=C2)
ax[1].set_ylim(1e-5, 3); ax[1].set_xlabel("frequency / Hz"); ax[1].set_ylabel("|FFT of $\\xi$| (normalised, log)"); ax[1].grid(alpha=0.3, which="both")
ax[1].set_title("(b) harmonics of $\\xi=x_2-x_1$: $d=24.5$ mm, $\\xi_0=4.2$ mm, undamped", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig("fig_nonlinear.pdf")

# ------------------------------------------------------------------ 4. damping test: equal vs unequal decay of the two modes
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
for a, (t1, t2, lab) in zip(ax, [(21, 21, "(a) equal decay times (internal friction, air): minima/maxima stay $=x_{20}/A$"),
                                  (21, 10.5, "(b) anti-phase mode decaying twice as fast (e.g. eddy currents)")]):
    e1 = abs(eta0)*np.exp(-t/t1); e2 = xi0*np.exp(-t/t2)
    env = 0.5*np.sqrt(e1**2 + e2**2 - 2*e1*e2*np.cos((w2 - w1)*t)); env_min = 0.5*np.abs(e1 - e2); env_max = 0.5*(e1 + e2)   # e1 = |eta0|e^-t/tau1, e2 = xi0 e^-t/tau2
    m = t < 30
    a.plot(t[m], env[m], color=C1, lw=0.7); a.plot(t[m], env_max[m], "--", color="k", lw=0.8, label="$\\frac{1}{2}(\\eta_0e^{-t/\\tau_1}+\\xi_0e^{-t/\\tau_2})$")
    a.plot(t[m], env_min[m], ":", color="k", lw=1.0, label="$\\frac{1}{2}|\\eta_0e^{-t/\\tau_1}-\\xi_0e^{-t/\\tau_2}|$")
    a.set_xlabel("$t$ / s"); a.set_title(lab, loc="left", fontsize=9); a.grid(alpha=0.3); a.legend(frameon=False, fontsize=8)
ax[0].set_ylabel("envelope of $x_1$ / mm")
fig.tight_layout(); fig.savefig("fig_damping.pdf")
print("figures written; x20 = %.2f mm, eta0 = %.2f, xi0 = %.2f, ratio %.3f, (f1/f2)^2 = %.3f, beat min/max = %.3f" %
      (x20, eta0, xi0, xi0/abs(eta0), (w1/w2)**2, x20/A))
