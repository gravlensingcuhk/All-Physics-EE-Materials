"""
make_figures.py  --  generates every data-driven figure and number used in the EE.

Run:  python3 make_figures.py          (takes ~2-4 min)
Outputs: figs/*.pdf  and  numbers.tex  (LaTeX macros with the computed values)

Everything marked ASSUMED is a placeholder that must be replaced by the measured value.
"""
import json, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, curve_fit
from scipy.integrate import solve_ivp
from pathlib import Path

HERE = Path(__file__).parent
FIG = HERE / "figs"; FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "legend.fontsize": 8,
                     "figure.dpi": 110, "savefig.bbox": "tight"})
NUM = {}                                   # name -> string, exported to numbers.tex
def put(name, val, fmt="{:.3g}"):
    NUM[name] = fmt.format(val) if not isinstance(val, str) else val
def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf"); plt.close(fig); print("  saved", name)

mu0 = 4e-7 * np.pi
# ----------------------------------------------------------------- measured data (means of 3 trials)
d_mm = np.array([21, 23, 25, 28, 31, 35, 40.0])
f1 = np.array([0.823, 0.808, 0.812, 0.810, 0.807, 0.806, 0.802])
f2 = np.array([1.063, 1.005, 0.953, 0.907, 0.881, 0.845, 0.811])
d_unres = np.array([45, 50, 60.0])
DF = 0.005            # ASSUMED uncertainty of each frequency / Hz (replace by half-range of trials)
DD = 0.5              # ASSUMED uncertainty of d / mm
T_REC, FPS = 30.0, 240.0      # user's acquisition
TAU = 25.0            # ASSUMED amplitude decay time / s
A_REL = 10.0          # release displacement used in the original runs / mm
six = slice(0, 6)

# =============================================================== 1. analysis of the measured data
print("1. analysis of measured data")
Y = f2**2 - f1**2
dY = 2 * (f1 * DF + f2 * DF)
lnY, dlnY = np.log(Y), dY / Y
lnd, dlnd = np.log(d_mm), DD / d_mm
p6, cov6 = np.polyfit(lnd[six], lnY[six], 1, cov=True)
se6 = np.sqrt(np.diag(cov6))
res6 = lnY[six] - np.polyval(p6, lnd[six])
p7 = np.polyfit(lnd, lnY, 1)
w = 1 / dlnY[six]**2
pw, covw = np.polyfit(lnd[six], lnY[six], 1, w=np.sqrt(w), cov="unscaled")
# best / worst lines through the extreme error bars of the first and last resolved points
x1, x6 = lnd[0], lnd[5]; y1, y6 = lnY[0], lnY[5]
g_steep = ((y6 - dlnY[5]) - (y1 + dlnY[0])) / ((x6 - dlnd[5]) - (x1 + dlnd[0]))
g_shallow = ((y6 + dlnY[5]) - (y1 - dlnY[0])) / ((x6 + dlnd[5]) - (x1 - dlnd[0]))
# Graph 2: Y vs d^-5 (SI) through the origin
X5 = (d_mm[six] * 1e-3) ** -5
C_exp = np.sum(X5 * Y[six]) / np.sum(X5**2)
C_se = np.sqrt(np.sum((Y[six] - C_exp * X5)**2) / 5 / np.sum(X5**2))
pfree = np.polyfit(X5, Y[six], 1)
# free power law on the untransformed Y
pw_n, cw_n = curve_fit(lambda d, C, n: C * d**-n, d_mm[six], Y[six], p0=[4e4, 3.7], sigma=dY[six])
# f1 drift
pf1, cf1 = np.polyfit(d_mm[six], f1[six], 1, cov=True)
# implied mu_m^2/m from each point and from C_exp
mm_ratio = Y * np.pi**3 * (d_mm * 1e-3)**5 / (3 * mu0)
mm_from_C = C_exp * np.pi**3 / (3 * mu0)
# static deflection delta = pi^2 Y d / (2 w_a^2) (independent of m and mu_m)
wa2 = (2 * np.pi * f1[6])**2
delta = np.pi**2 * Y * (d_mm) / (2 * wa2)             # mm
p_set = np.polyfit(np.log(d_mm[six] + 2 * delta[six]), lnY[six], 1)
for k, v in dict(gradSix=p6[0], gradSixErr=se6[0], intSix=p6[1], intSixErr=se6[1],
                 gradSeven=p7[0], gradW=pw[0], gradWErr=np.sqrt(covw[0, 0]),
                 gradSteep=g_steep, gradShallow=g_shallow, Cexp=C_exp, CexpSE=C_se,
                 nFree=pw_n[1], nFreeErr=np.sqrt(cw_n[1, 1]), fOneSlope=pf1[0] * 1e3,
                 fOneSlopeErr=np.sqrt(cf1[0, 0]) * 1e3, mmFromC=mm_from_C,
                 gradSet=p_set[0], deltaTwentyOne=delta[0], deltaThirtyFive=delta[5]).items():
    put(k, v, "{:.4g}" if k in ("Cexp", "CexpSE") else "{:.3g}")
put("Cexp", C_exp, "{:.2e}"); put("CexpSE", C_se, "{:.1e}")
put("intSix", p6[1], "{:.2f}"); put("gradSix", p6[0], "{:.2f}"); put("gradSixErr", se6[0], "{:.2f}")
put("gradW", pw[0], "{:.2f}"); put("gradWErr", np.sqrt(covw[0, 0]), "{:.2f}")
put("gradSteep", g_steep, "{:.2f}"); put("gradShallow", g_shallow, "{:.2f}")
put("gradSeven", p7[0], "{:.2f}"); put("gradSet", p_set[0], "{:.2f}")
put("nFree", pw_n[1], "{:.2f}"); put("nFreeErr", np.sqrt(cw_n[1, 1]), "{:.2f}")
put("mmFromC", mm_from_C, "{:.3f}")
put("pctGrad", 100 * (5 - abs(p6[0])) / 5, "{:.0f}")
print(f"   6-pt gradient {p6[0]:.3f}+-{se6[0]:.3f}  int {p6[1]:.2f}; weighted {pw[0]:.3f}; 7-pt {p7[0]:.3f}")
print(f"   steep/shallow {g_steep:.2f}/{g_shallow:.2f}; C_exp {C_exp:.3e}+-{C_se:.1e}; free n {pw_n[1]:.2f}")
print(f"   f1 slope {pf1[0]*1e3:.2f} mHz/mm; mu^2/m from C {mm_from_C:.4f}; delta(mm) {delta.round(2)}; grad(d_set) {p_set[0]:.2f}")

# table rows for the essay (processed data)
rows = []
for i in range(7):
    rows.append(dict(d=d_mm[i], f1=f1[i], f2=f2[i], Y=Y[i], dY=dY[i], lnY=lnY[i], dlnY=dlnY[i],
                     lnd=lnd[i], dlnd=dlnd[i], X5=(d_mm[i]*1e-3)**-5, dX5=5*DD/d_mm[i]*(d_mm[i]*1e-3)**-5,
                     delta=delta[i], mm=mm_ratio[i]))
json.dump(rows, open(HERE / "processed_rows.json", "w"), indent=1)

# ---- Graph 1: ln Y vs ln d
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.errorbar(lnd[six], lnY[six], xerr=dlnd[six], yerr=dlnY[six], fmt="o", color="tab:red", ms=5,
            capsize=3, label="resolved points (21–35 mm)")
ax.errorbar(lnd[6], lnY[6], xerr=dlnd[6], yerr=dlnY[6], fmt="o", mfc="white", color="tab:red", ms=6,
            capsize=3, label="40 mm point (peaks not resolved, excluded)")
xx = np.linspace(2.98, 3.72, 50)
ax.plot(xx, np.polyval(p6, xx), "k-", lw=1.4, label=f"unweighted LINEST, gradient {p6[0]:.2f} ± {se6[0]:.2f}")
ax.plot(xx, np.polyval(pw, xx), "-", color="tab:blue", lw=1.1, label=f"weighted fit, gradient {pw[0]:.2f} ± {np.sqrt(covw[0,0]):.2f}")
ax.plot(xx, lnY[2] - 5 * (xx - lnd[2]), "--", color="tab:green", lw=1.3, label="dipole prediction: gradient −5 (through 25 mm point)")
ax.set_xlabel(r"$\ln(d\,/\,\mathrm{mm})$"); ax.set_ylabel(r"$\ln\left[(f_2^2-f_1^2)\,/\,\mathrm{Hz^2}\right]$")
ax.legend(loc="lower left"); ax.grid(alpha=0.3)
save(fig, "G1_lnln")

# ---- Graph 2: Y vs d^-5
fig, ax = plt.subplots(figsize=(6.2, 4.2))
X5all = (d_mm * 1e-3) ** -5
ax.errorbar(X5all[six]/1e8, Y[six], xerr=5*DD/d_mm[six]*X5all[six]/1e8, yerr=dY[six], fmt="o", color="tab:red", ms=5, capsize=3, label="resolved points")
ax.errorbar(X5all[6]/1e8, Y[6], yerr=dY[6], fmt="o", mfc="white", color="tab:red", ms=6, capsize=3, label="40 mm (excluded)")
xx = np.linspace(0, 2.6e8, 20)
ax.plot(xx/1e8, C_exp * xx, "k-", lw=1.4, label=f"through origin: $C_{{\\rm exp}}$ = {C_exp:.2e} m$^5$ s$^{{-2}}$")
ax.plot(xx/1e8, np.polyval(pfree, xx), ":", color="tab:blue", lw=1.2, label=f"free intercept: gradient {pfree[0]:.2e}, intercept {pfree[1]:.3f} Hz$^2$")
ax.set_xlabel(r"$d^{-5}\;/\;10^{8}\ \mathrm{m^{-5}}$"); ax.set_ylabel(r"$f_2^2-f_1^2\;/\;\mathrm{Hz^2}$")
ax.set_xlim(0, 2.6); ax.set_ylim(0, 0.55); ax.legend(loc="upper left"); ax.grid(alpha=0.3)
save(fig, "G2_Y_vs_d5")

# ---- Graph 3: f1, f2 vs d + residuals
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.8))
ax = axs[0]
ax.errorbar(d_mm, f1, xerr=DD, yerr=DF, fmt="s", color="tab:blue", ms=5, capsize=3, label="$f_1$ (in-phase)")
ax.errorbar(d_mm, f2, xerr=DD, yerr=DF, fmt="o", color="tab:red", ms=5, capsize=3, label="$f_2$ (anti-phase)")
ax.axhline(np.mean(f1[six]), color="tab:blue", ls="--", lw=1, label=f"mean $f_1$ = {np.mean(f1[six]):.3f} Hz")
dd = np.linspace(20, 62, 200)
ax.plot(dd, np.sqrt(np.mean(f1[six])**2 + C_exp * (dd*1e-3)**-5), "--", color="tab:red", lw=1, label=r"dipole: $\sqrt{f_1^2 + C_{\rm exp} d^{-5}}$")
for x in d_unres: ax.axvspan(x-0.4, x+0.4, color="0.85")
ax.text(52, 1.05, "single peak\nonly", ha="center", fontsize=8, color="0.3")
ax.set_xlabel("$d$ / mm"); ax.set_ylabel("frequency / Hz"); ax.legend(loc="upper right"); ax.grid(alpha=0.3)
ax.set_ylim(0.78, 1.10)
ax = axs[1]
ax.bar(range(6), res6, color=["tab:red" if r < 0 else "tab:blue" for r in res6])
ax.errorbar(range(6), res6, yerr=dlnY[six], fmt="none", ecolor="k", capsize=3)
ax.set_xticks(range(6)); ax.set_xticklabels([f"{int(x)}" for x in d_mm[six]])
ax.set_xlabel("$d$ / mm"); ax.set_ylabel(r"residual of $\ln Y$ (data − LINEST)"); ax.axhline(0, color="k", lw=0.8)
ax.set_title("residual pattern: − + − + + −", fontsize=9); ax.grid(alpha=0.3, axis="y")
save(fig, "G3_f_vs_d")

# ---- Graph 4: average frequency and beat frequency (the two derived quantities named in the research question)
f_avg_m, f_beat_m = (f1 + f2) / 2, f2 - f1
df_avg, df_beat = DF / np.sqrt(2), DF * np.sqrt(2)             # independent errors of f1 and f2
f1bar = np.mean(f1[six])
dd = np.linspace(20, 62, 300); f2_pred = np.sqrt(f1bar**2 + C_exp * (dd*1e-3)**-5)
p_beat, cov_beat = np.polyfit(np.log(d_mm[six]), np.log(f_beat_m[six]), 1, cov=True)
put("favgSmall", f_avg_m[0]); put("favgLarge", f_avg_m[6]); put("fbeatSmall", f_beat_m[0]); put("fbeatLarge", f_beat_m[6])
put("TbeatSmall", 1/f_beat_m[0], "{:.1f}"); put("TbeatLarge", 1/f_beat_m[6], "{:.0f}")
put("gradBeat", p_beat[0], "{:.2f}"); put("gradBeatErr", np.sqrt(cov_beat[0, 0]), "{:.2f}")
put("favgPredSmall", 0.5*(f1bar + np.sqrt(f1bar**2 + C_exp*0.021**-5))); put("fbeatPredSmall", np.sqrt(f1bar**2 + C_exp*0.021**-5) - f1bar)
put("fbeatPredForty", np.sqrt(f1bar**2 + C_exp*0.040**-5) - f1bar, "{:.3f}")
put("fbeatPredSixty", 1e3*(np.sqrt(f1bar**2 + C_exp*0.060**-5) - f1bar), "{:.1f}"); put("TbeatPredSixty", (np.sqrt(f1bar**2 + C_exp*0.060**-5) - f1bar)**-1/60, "{:.0f}")
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.8))
ax = axs[0]
ax.errorbar(d_mm, f_avg_m, xerr=DD, yerr=df_avg, fmt="D", color="tab:purple", ms=5, capsize=3, label=r"$f_{\rm avg}=(f_1+f_2)/2$")
ax.plot(dd, 0.5*(f1bar + f2_pred), "--", color="tab:purple", lw=1, label=r"dipole: $\frac{1}{2}\left(f_1+\sqrt{f_1^2+C_{\rm exp}d^{-5}}\right)$")
ax.axhline(f1bar, color="tab:blue", ls=":", lw=1, label=rf"mean $f_1$ = {f1bar:.3f} Hz")
for x in d_unres: ax.axvspan(x-0.4, x+0.4, color="0.85")
ax.text(52, 0.784, "single peak\nonly", ha="center", va="bottom", fontsize=8, color="0.3")
ax.set_xlabel("$d$ / mm"); ax.set_ylabel(r"$f_{\rm avg}$ / Hz"); ax.set_ylim(0.78, 0.98); ax.legend(loc="upper right"); ax.grid(alpha=0.3)
ax.set_title("average frequency", fontsize=9)
ax = axs[1]
ax.errorbar(d_mm, f_beat_m, xerr=DD, yerr=df_beat, fmt="o", color="tab:green", ms=5, capsize=3, label=r"$f_{\rm beat}=f_2-f_1$")
ax.plot(dd, f2_pred - f1bar, "--", color="tab:green", lw=1, label=r"dipole: $\sqrt{f_1^2+C_{\rm exp}d^{-5}}-f_1$")
ax.plot(dd, np.exp(np.polyval(p_beat, np.log(dd))), ":", color="k", lw=1, label=rf"power fit 21–35 mm: gradient ${p_beat[0]:.2f}\pm{np.sqrt(cov_beat[0,0]):.2f}$")
ax.errorbar(d_unres, [1/T_REC]*3, yerr=[[0.3/T_REC]*3, [0]*3], uplims=True, fmt="v", mfc="white", color="0.4", ms=5, label=rf"unresolved: $f_{{\rm beat}}<1/T$ = {1e3/T_REC:.0f} mHz")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([20, 25, 30, 40, 50, 60]); ax.set_xticklabels(["20", "25", "30", "40", "50", "60"])
ax.set_xlabel("$d$ / mm"); ax.set_ylabel(r"$f_{\rm beat}$ / Hz"); ax.set_ylim(4e-3, 0.5); ax.legend(loc="lower left", fontsize=7.5); ax.grid(alpha=0.3, which="both")
ax.set_title("beat frequency (log–log axes)", fontsize=9)
save(fig, "G4_favg_fbeat")

# =============================================================== 2. FFT background figures
print("2. FFT background")
def spectrum_db(win, nfft=2**14):
    X = np.abs(np.fft.rfft(win, nfft)); X /= X.max()
    return np.arange(len(X)) * len(win) / nfft, 20 * np.log10(X + 1e-12)   # x in bins
N = 256
wins = {"rectangular": np.ones(N), "Hann": np.hanning(N), "Blackman": np.blackman(N)}
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.4))
for (nm, wv), c in zip(wins.items(), ["tab:red", "tab:blue", "tab:green"]):
    axs[0].plot(np.arange(N) / N, wv, color=c, label=nm)
    b, db = spectrum_db(wv); axs[1].plot(b, db, color=c, label=nm)
axs[0].set_xlabel("$t/T$"); axs[0].set_ylabel("window value $w(t)$"); axs[0].legend(); axs[0].grid(alpha=0.3)
axs[1].set_xlim(0, 8); axs[1].set_ylim(-100, 3); axs[1].set_xlabel("frequency offset from the tone / bins ($1/T$)")
axs[1].set_ylabel("magnitude / dB"); axs[1].legend(); axs[1].grid(alpha=0.3)
axs[1].annotate("main-lobe half-width:\nrect 1 bin, Hann 2, Blackman 3", xy=(3, -60), fontsize=8)
save(fig, "F_windows")

# leakage demo: 0.85 Hz for 30 s (25.5 cycles) vs 26 cycles
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.4))
t = np.arange(0, T_REC, 1 / FPS)
for ax, fr, ttl in [(axs[0], 26 / T_REC, "integer number of cycles (26)"), (axs[1], 25.5 / T_REC, "non-integer number of cycles (25.5)")]:
    x = np.cos(2 * np.pi * fr * t)
    for nm, c in [("rectangular", "tab:red"), ("Hann", "tab:blue")]:
        wv = np.ones(len(t)) if nm == "rectangular" else np.hanning(len(t))
        X = np.abs(np.fft.rfft(x * wv, 2**18)); X /= X.max(); fq = np.fft.rfftfreq(2**18, 1 / FPS)
        ax.plot(fq, 20 * np.log10(X + 1e-12), color=c, lw=1, label=nm + " window")
        Xb = np.abs(np.fft.rfft(x * wv)); Xb /= Xb.max(); fb = np.fft.rfftfreq(len(t), 1 / FPS)
        ax.plot(fb, 20 * np.log10(Xb + 1e-12), "o", color=c, ms=3)
    ax.set_xlim(0.6, 1.1); ax.set_ylim(-70, 3); ax.axvline(fr, color="k", ls=":", lw=0.8)
    ax.set_title(ttl, fontsize=9); ax.set_xlabel("frequency / Hz"); ax.grid(alpha=0.3)
axs[0].set_ylabel("|FFT| / dB (dots: DFT bins, lines: zero-padded)"); axs[0].legend(loc="upper right")
save(fig, "F_leakage")

# =============================================================== 3. resolution / pulling at the real T = 30 s
print("3. two-tone resolution test")
def peaks_two(x, fs, window, nfft=2**19, fmin=0.4, fmax=2.0, thr=0.03):
    n = len(x); x = x - x.mean()
    wv = {"rect": np.ones(n), "hann": np.hanning(n), "blackman": np.blackman(n)}[window]
    X = np.abs(np.fft.rfft(x * wv, nfft)); fq = np.fft.rfftfreq(nfft, 1 / fs)
    m = (fq > fmin) & (fq < fmax); idx = np.where(m)[0]
    Xm = X[idx].max(); pk = []
    for i in idx[1:-1]:
        if X[i] > X[i-1] and X[i] >= X[i+1] and X[i] > thr * Xm:
            a, b, c = X[i-1], X[i], X[i+1]; p = 0.5 * (a - c) / (a - 2*b + c)
            pk.append((fq[i] + p * (fq[1] - fq[0]), b))
    pk.sort(key=lambda z: -z[1]); return sorted(q[0] for q in pk[:2]), (fq, X / Xm)

rng = np.random.default_rng(1)
pull = {}
fig, axs = plt.subplots(2, 3, figsize=(10.5, 6), sharex="col")
for col, dsel in enumerate([25, 35, 40]):
    i = list(d_mm).index(dsel)
    for row, T in enumerate([30, 90]):
        t = np.arange(0, T, 1 / FPS)
        x = 0.5 * np.exp(-t / TAU) * (np.cos(2*np.pi*f1[i]*t) + np.cos(2*np.pi*f2[i]*t)) * A_REL + rng.normal(0, 0.03, len(t))
        ax = axs[row, col]
        for win, c in [("rect", "tab:red"), ("hann", "tab:blue")]:
            pk, (fq, X) = peaks_two(x, FPS, win)
            ax.plot(fq, X, color=c, lw=1, label=win)
            ax.plot(pk, np.interp(pk, fq, X), "v", color=c, ms=6)
            pull[(dsel, T, win)] = pk
        for fr in (f1[i], f2[i]): ax.axvline(fr, color="k", ls=":", lw=0.8)
        ax.set_xlim(0.6, 1.25); ax.set_title(f"$d$ = {dsel} mm, $T$ = {T} s, bin = {1/T*1000:.1f} mHz", fontsize=9); ax.grid(alpha=0.3)
        if col == 0: ax.set_ylabel("|FFT| (normalised)")
        if row == 1: ax.set_xlabel("frequency / Hz")
axs[0, 0].legend(loc="upper right")
save(fig, "F_resolution")
# full pulling table for all d at T=30 and T=90, rect + hann
pull_rows = []
for i, dsel in enumerate(d_mm):
    r = dict(d=dsel, df=f2[i] - f1[i], bins30=(f2[i]-f1[i]) * 30, Tmin_rect=1/(f2[i]-f1[i]), Tmin_hann=2/(f2[i]-f1[i]))
    for T in (30, 90):
        t = np.arange(0, T, 1 / FPS)
        x = 0.5*np.exp(-t/TAU)*(np.cos(2*np.pi*f1[i]*t)+np.cos(2*np.pi*f2[i]*t))*A_REL + rng.normal(0, 0.03, len(t))
        for win in ("rect", "hann"):
            pk, _ = peaks_two(x, FPS, win)
            if len(pk) == 2:
                r[f"e1_{T}_{win}"] = (pk[0]-f1[i])*1e3; r[f"e2_{T}_{win}"] = (pk[1]-f2[i])*1e3
            else:
                r[f"e1_{T}_{win}"] = np.nan; r[f"e2_{T}_{win}"] = np.nan
    pull_rows.append(r)
json.dump(pull_rows, open(HERE / "pulling_rows.json", "w"), indent=1)
for r in pull_rows:
    print("   ", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})

# =============================================================== 4. non-linear (RK4-type) amplitude study
print("4. amplitude non-linearity (full dipole force, identical springs)")
wa = 2 * np.pi * 0.806
Kmm = Y[2] * 2 * np.pi**2 * d_mm[2]**5        # gamma/m = Kmm * d_mm^-5  (s^-2) anchored on the 25 mm point
K_SI = Kmm * 1e-15
def rhs_full(t, s, d, wa2_, wb2_, K, tau):
    x1, v1, x2, v2 = s
    r = d + x2 - x1
    dF = K / 4 * (r**-4 - d**-4)      # [F(r) - F(d)]/m
    a1 = -wa2_ * x1 - dF - (2/tau) * v1 if tau else -wa2_ * x1 - dF
    a2 = -wb2_ * x2 + dF - (2/tau) * v2 if tau else -wb2_ * x2 + dF
    return [v1, a1, v2, a2]
def freq_from_zero_crossings(t, x):
    x = x - x.mean(); idx = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = t[idx] - x[idx] * (t[idx+1] - t[idx]) / (x[idx+1] - x[idx])
    return (len(tc) - 1) / (tc[-1] - tc[0]) if len(tc) > 2 else np.nan
amps = [1, 2, 3, 5, 7, 10]
shift = {}; f2lin = {}
for dsel in d_mm:
    d = dsel * 1e-3; gam_m = K_SI * d**-5
    f2lin[dsel] = np.sqrt(wa**2 + 2 * gam_m) / (2*np.pi)
    for A in amps:
        sol = solve_ivp(rhs_full, (0, 60), [A*1e-3, 0, 0, 0], args=(d, wa**2, wa**2, K_SI, None),
                        rtol=1e-10, atol=1e-13, dense_output=True)
        t = np.arange(0, 60, 1/FPS); s = sol.sol(t)
        xi = s[0] - s[2]
        shift[(dsel, A)] = freq_from_zero_crossings(t, xi)
for dsel in d_mm:
    print(f"   d={dsel:4.0f}: f2(lin)={f2lin[dsel]:.4f}", " ".join(f"A={A}:{100*(shift[(dsel,A)]/f2lin[dsel]-1):+.2f}%" for A in amps))
# gradient with A = 10 and 3 mm (identical springs, f1 unchanged in this model)
def grad_from(fs2):
    yy = np.log(np.array([fs2[k]**2 - 0.806**2 for k in d_mm[six]])); return np.polyfit(lnd[six], yy, 1)[0]
g_lin = grad_from(f2lin); g10 = grad_from({k: shift[(k, 10)] for k in d_mm}); g3 = grad_from({k: shift[(k, 3)] for k in d_mm})
print(f"   gradient: linear {g_lin:.3f}, A=3 {g3:.3f}, A=10 {g10:.3f}")
put("gradAten", g10, "{:.2f}"); put("gradAthree", g3, "{:.2f}")
put("shiftTwentyOneAten", 100*(shift[(21, 10)]/f2lin[21]-1), "{:.1f}"); put("shiftThirtyFiveAten", 100*(shift[(35, 10)]/f2lin[35]-1), "{:.1f}")
put("shiftTwentyOneAthree", 100*(shift[(21, 3)]/f2lin[21]-1), "{:.2f}"); put("shiftThirtyFiveAthree", 100*(shift[(35, 3)]/f2lin[35]-1), "{:.2f}")
put("shiftTwentyOneAtwo", 100*(shift[(21, 2)]/f2lin[21]-1), "{:.2f}")
# Landau-Lifshitz prediction for comparison
def LL(dsel, A):
    g = 1 - (0.806/f2lin[dsel])**2; return (A/dsel)**2 * (15*g/8 - 125*g**2/48) * 100
put("LLTwentyOneAten", LL(21, 10), "{:.1f}"); put("LLThirtyFiveAten", LL(35, 10), "{:.1f}")
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
ax = axs[0]
for dsel, c in [(21, "tab:red"), (25, "tab:orange"), (31, "tab:green"), (35, "tab:blue")]:
    ax.plot(amps, [100*(shift[(dsel, A)]/f2lin[dsel]-1) for A in amps], "o-", color=c, label=f"simulation, $d$ = {dsel} mm")
    ax.plot(amps, [LL(dsel, A) for A in amps], "--", color=c, lw=0.9)
ax.set_xlabel("release amplitude $A$ / mm"); ax.set_ylabel(r"$(f_2(A)-f_2(0))/f_2(0)$ / %"); ax.grid(alpha=0.3)
ax.legend(); ax.set_title("solid: numerical integration; dashed: perturbation formula", fontsize=9)
ax = axs[1]
for A, c, lab in [(None, "k", f"linear theory (gradient {g_lin:.2f})"), (3, "tab:blue", f"A = 3 mm (gradient {g3:.2f})"), (10, "tab:red", f"A = 10 mm (gradient {g10:.2f})")]:
    fs2 = f2lin if A is None else {k: shift[(k, A)] for k in d_mm}
    ax.plot(lnd[six], np.log([fs2[k]**2 - 0.806**2 for k in d_mm[six]]), "o-" if A else "-", color=c, label=lab, ms=4)
ax.set_xlabel(r"$\ln(d/\mathrm{mm})$"); ax.set_ylabel(r"$\ln(f_2^2-f_1^2)$"); ax.legend(); ax.grid(alpha=0.3)
ax.set_title("amplitude makes the ln–ln line steeper, not flatter", fontsize=9)
save(fig, "F_amplitude")

# =============================================================== 5. finite-size (Gilbert) magnets
print("5. finite-size magnets")
def disc_disc_force(z, R, nr=24, nphi=36):
    r = (np.arange(nr) + 0.5) / nr * R; phi = (np.arange(nphi) + 0.5) / nphi * 2*np.pi
    dA = (R/nr) * (2*np.pi/nphi) * r
    rr, pp = np.meshgrid(r, phi, indexing="ij"); x1_, y1_ = rr*np.cos(pp), rr*np.sin(pp)
    dAA = np.broadcast_to(dA[:, None], rr.shape); F = 0.0
    for i in range(nr):
        dist2 = (x1_ - r[i])**2 + y1_**2 + z**2
        F += np.sum(dAA * (dA[i]*nphi) * z / dist2**1.5)
    return F
def magnet_force(d, D, t):
    R = D/2; return disc_disc_force(d-t, R) + disc_disc_force(d+t, R) - 2*disc_disc_force(d, R)
def gamma_fs(d, D, t, h=0.05):
    return -(magnet_force(d+h, D, t) - magnet_force(d-h, D, t)) / (2*h)
dgrid = np.linspace(15, 65, 26)
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
fs_grad = {}
for D, t, c in [(6, 3, "tab:green"), (10, 3, "tab:olive"), (15, 5, "tab:orange"), (20, 5, "tab:red")]:
    Fg = np.array([magnet_force(x, D, t) for x in dgrid])
    Fdip = Fg[-1] * (dgrid / dgrid[-1])**-4          # dipole law normalised at 65 mm
    axs[0].plot(dgrid, Fg / Fdip, "-", color=c, label=f"$D$ = {D} mm, $t$ = {t} mm")
    g = np.array([gamma_fs(x, D, t) for x in d_mm[six]])
    fs_grad[(D, t)] = np.polyfit(lnd[six], np.log(g), 1)[0]
axs[0].axvspan(21, 35, color="0.9", label="measured range 21–35 mm")
axs[0].set_xlabel("$d$ / mm"); axs[0].set_ylabel(r"$F_{\rm disc}(d)\,/\,F_{\rm dipole}(d)$"); axs[0].legend(fontsize=7); axs[0].grid(alpha=0.3)
axs[0].set_title("real magnets repel less than the point-dipole law at short range", fontsize=9)
Ds = np.arange(4, 27, 1.5); 
for t, c in [(3, "tab:blue"), (5, "tab:red"), (8, "tab:green")]:
    sl = [np.polyfit(lnd[six], np.log([gamma_fs(x, D, t) for x in d_mm[six]]), 1)[0] for D in Ds]
    axs[1].plot(Ds, sl, "o-", ms=3, color=c, label=f"thickness $t$ = {t} mm")
axs[1].axhline(-5, color="k", ls="--", lw=1, label="point dipole (−5)")
axs[1].axhline(p6[0], color="tab:purple", ls=":", lw=1.5, label=f"measured ({p6[0]:.2f})")
axs[1].set_xlabel("magnet diameter $D$ / mm"); axs[1].set_ylabel("apparent ln–ln gradient over 21–35 mm"); axs[1].legend(fontsize=7); axs[1].grid(alpha=0.3)
save(fig, "F_finite_size")
for k, v in fs_grad.items(): print(f"   D={k[0]}, t={k[1]}: apparent gradient {v:.2f}")
put("fsGradSix", fs_grad[(6, 3)], "{:.2f}"); put("fsGradTen", fs_grad[(10, 3)], "{:.2f}")
put("fsGradFifteen", fs_grad[(15, 5)], "{:.2f}"); put("fsGradTwenty", fs_grad[(20, 5)], "{:.2f}")

# =============================================================== 6. non-identical springs
print("6. detuned-spring model")
def modes(d, fa, fb, K, n=5):
    wa2_ = (2*np.pi*fa)**2; wb2_ = (2*np.pi*fb)**2; kap = K * d**-n
    wbar2 = (wa2_ + wb2_)/2; Dl = (wb2_ - wa2_)/2
    return np.sqrt(wbar2 + kap - np.sqrt(Dl**2 + kap**2))/(2*np.pi), np.sqrt(wbar2 + kap + np.sqrt(Dl**2 + kap**2))/(2*np.pi)
def resid(p, dd, F1, F2):
    a, b = modes(dd, *p); return np.concatenate([(a-F1)/DF, (b-F2)/DF])
fit = least_squares(resid, [0.80, 0.81, 1e7], args=(d_mm[six], f1[six], f2[six]))
fa_fit, fb_fit, K_fit = fit.x
J = fit.jac; covd = np.linalg.inv(J.T @ J) * (fit.fun @ fit.fun) / (len(fit.fun) - 3); sd = np.sqrt(np.diag(covd))
print(f"   fa={fa_fit:.4f}+-{sd[0]:.4f} fb={fb_fit:.4f}+-{sd[1]:.4f} K={K_fit:.3e}  chi2={fit.fun@fit.fun:.1f}")
put("faFit", fa_fit, "{:.3f}"); put("fbFit", fb_fit, "{:.3f}"); put("faFitErr", sd[0], "{:.3f}"); put("fbFitErr", sd[1], "{:.3f}")
put("KFit", K_fit, "{:.2e}"); put("chiDetuned", fit.fun @ fit.fun, "{:.1f}")
put("YsatDetuned", abs(fb_fit**2 - fa_fit**2), "{:.3f}")
m1_, m2_ = modes(d_mm, fa_fit, fb_fit, K_fit)
put("detFoneSmall", m1_[0], "{:.3f}"); put("detFoneLarge", m1_[6], "{:.3f}")
# also a pure-dipole (identical springs) chi2 for comparison: f1 = const, f2 from C
def resid_id(p, dd, F1, F2):
    f0, K = p; return np.concatenate([(f0 - F1)/DF, (np.sqrt(f0**2 + K*dd**-5) - F2)/DF])
fit_id = least_squares(resid_id, [0.81, 1e6], args=(d_mm[six], f1[six], f2[six]))
put("chiIdentical", fit_id.fun @ fit_id.fun, "{:.1f}")
print(f"   identical-spring dipole fit chi2 = {fit_id.fun@fit_id.fun:.1f}")
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
dd = np.linspace(20, 62, 300); a_, b_ = modes(dd, fa_fit, fb_fit, K_fit)
ax = axs[0]
ax.errorbar(d_mm, f1, yerr=DF, fmt="s", color="tab:blue", ms=5, capsize=3, label="$f_1$ data")
ax.errorbar(d_mm, f2, yerr=DF, fmt="o", color="tab:red", ms=5, capsize=3, label="$f_2$ data")
ax.plot(dd, a_, "-", color="tab:blue", lw=1.2, label=f"detuned model, $f_a$ = {fa_fit:.3f}, $f_b$ = {fb_fit:.3f} Hz")
ax.plot(dd, b_, "-", color="tab:red", lw=1.2)
ax.axhline(fa_fit, color="tab:blue", ls=":", lw=0.8); ax.axhline(fb_fit, color="tab:red", ls=":", lw=0.8)
ax.set_xlabel("$d$ / mm"); ax.set_ylabel("frequency / Hz"); ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_ylim(0.78, 1.10)
ax = axs[1]
ax.plot(np.log(dd), np.log(b_**2 - a_**2), "-", color="tab:purple", lw=1.4, label="detuned model")
ax.plot(np.log(dd), np.log(K_fit*dd**-5/(2*np.pi**2)), "k--", lw=1, label="identical springs (gradient −5)")
ax.errorbar(lnd, lnY, yerr=dlnY, fmt="o", color="tab:red", ms=4, capsize=2, label="data")
ax.axhline(np.log(abs(fb_fit**2-fa_fit**2)), color="0.4", ls=":", label=r"saturation $\ln|f_b^2-f_a^2|$")
ax.set_xlabel(r"$\ln(d/\mathrm{mm})$"); ax.set_ylabel(r"$\ln(f_2^2-f_1^2)$"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
save(fig, "F_detuned")

# =============================================================== 7. static deflection
print("7. static deflection / definition of d")
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.4))
axs[0].plot(d_mm, delta, "o-", color="tab:brown"); axs[0].set_xlabel("$d$ / mm"); axs[0].set_ylabel(r"static deflection $\delta$ of each spring / mm")
axs[0].grid(alpha=0.3); axs[0].set_title(r"$\delta=\pi^2 (f_2^2-f_1^2)\,d\,/\,(2\omega_a^2)$ from the data", fontsize=9)
axs[1].plot(lnd[six], lnY[six], "o", color="tab:red", label=f"against $\\ln d$ as recorded (gradient {p6[0]:.2f})")
axs[1].plot(np.log(d_mm[six]+2*delta[six]), lnY[six], "s", color="tab:green", label=f"against $\\ln(d+2\\delta)$ (gradient {p_set[0]:.2f})")
for i in range(6): axs[1].annotate("", xy=(np.log(d_mm[i]+2*delta[i]), lnY[i]), xytext=(lnd[i], lnY[i]), arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
axs[1].set_xlabel(r"$\ln(d/\mathrm{mm})$"); axs[1].set_ylabel(r"$\ln(f_2^2-f_1^2)$"); axs[1].legend(fontsize=7); axs[1].grid(alpha=0.3)
axs[1].set_title("if $d$ had been set before the magnets pushed apart", fontsize=9)
save(fig, "F_deflection")

# =============================================================== 8. full synthetic pipeline (SIMULATED further data)
print("8. synthetic re-measurement pipeline (SIMULATED)")
T_NEW = 90.0; SIG_TRACK = 0.1e-3   # ASSUMED tracking noise 0.1 mm (rms)
# simulation truth: identical springs (f_a = f_b = mean measured f_1), exact dipole force anchored on the 25 mm point
F0_S = float(np.mean(f1[six])); K_s = K_SI
wa2_s = wb2_s = (2*np.pi*F0_S)**2
d_all = np.array([21, 23, 25, 28, 31, 35, 40, 45, 50, 60.0])
def simulate(dsel, A_mm, T=T_NEW, tau=TAU, seed=0):
    d = dsel*1e-3
    sol = solve_ivp(rhs_full, (0, T), [A_mm*1e-3, 0, 0, 0], args=(d, wa2_s, wb2_s, K_s, tau), rtol=1e-9, atol=1e-12, dense_output=True)
    t = np.arange(0, T, 1/FPS); s = sol.sol(t); r = np.random.default_rng(seed)
    return t, s[0] + r.normal(0, SIG_TRACK, len(t)), s[2] + r.normal(0, SIG_TRACK, len(t))
def one_peak(x, fs, window="rect", nfft=2**20, fmin=0.4, fmax=2.0):
    pk, _ = peaks_two(x, fs, window, nfft, fmin, fmax, thr=0.5); return pk[0] if pk else np.nan
def two_tone_fit(t, x, g1, g2, tau0=TAU):
    def model(p, t):
        A1, A2, ph1, ph2, fa_, fb_, ta, tb, c = p
        return A1*np.exp(-t/ta)*np.cos(2*np.pi*fa_*t+ph1) + A2*np.exp(-t/tb)*np.cos(2*np.pi*fb_*t+ph2) + c
    p0 = [x.std(), x.std(), 0, 0, g1, g2, tau0, tau0, 0]
    r = least_squares(lambda p: model(p, t) - x, p0, x_scale="jac")
    J = r.jac; cov = np.linalg.inv(J.T@J) * (r.fun@r.fun)/(len(t)-9); se = np.sqrt(np.diag(cov))
    fl = sorted([(r.x[4], se[4]), (r.x[5], se[5])]); return fl[0][0], fl[1][0], fl[0][1], fl[1][1]
true_f = {dsel: (F0_S, np.sqrt(F0_S**2 + K_s*(dsel*1e-3)**-5/(2*np.pi**2))) for dsel in d_all}
sim_rows = []
t0 = time.time()
for dsel in d_all:
    # ---- old procedure: A = 10 mm, T = 30 s, FFT of x1 only
    t, x1, x2 = simulate(dsel, A_REL, T=T_REC, seed=int(dsel)+100)
    row = dict(d=dsel, f1_true=true_f[dsel][0], f2_true=true_f[dsel][1])
    for win in ("rect", "hann"):
        pk, _ = peaks_two(x1, FPS, win, nfft=2**20)
        row[f"f1_old_{win}"], row[f"f2_old_{win}"] = (pk[0], pk[1]) if len(pk) == 2 else (pk[0], np.nan)
    # ---- new procedure: A = 3 mm, T = 90 s, both coordinates
    t, x1, x2 = simulate(dsel, 3.0, seed=int(dsel))
    for win in ("rect", "hann", "blackman"):
        pk, _ = peaks_two(x1, FPS, win, nfft=2**20)
        row[f"f1_{win}"], row[f"f2_{win}"] = (pk[0], pk[1]) if len(pk) == 2 else (pk[0], np.nan)
    eta, xi = x1 + x2, x1 - x2
    row["f1_eta"] = one_peak(eta, FPS); row["f2_xi"] = one_peak(xi, FPS)
    g1 = row["f1_eta"]; g2 = row["f2_xi"] if np.isfinite(row["f2_xi"]) else g1 + 0.01
    try:
        row["f1_fit"], row["f2_fit"], row["e1_fit"], row["e2_fit"] = two_tone_fit(t, x1, g1, g2)
    except Exception as e:
        row["f1_fit"] = row["f2_fit"] = row["e1_fit"] = row["e2_fit"] = np.nan
    sim_rows.append(row)
    print(f"   d={dsel:3.0f} true {row['f1_true']:.4f}/{row['f2_true']:.4f} old-rect {row['f1_old_rect']:.4f}/{row['f2_old_rect']:.4f} "
          f"rect90 {row['f1_rect']:.4f}/{row['f2_rect']:.4f} eta/xi {row['f1_eta']:.4f}/{row['f2_xi']:.4f} "
          f"fit {row['f1_fit']:.4f}±{row['e1_fit']*1e3:.2f} / {row['f2_fit']:.4f}±{row['e2_fit']*1e3:.2f} mHz ({time.time()-t0:.0f}s)")
    if dsel == 45: keep45 = (t, x1, x2)
    if dsel == 25: keep25 = (t, x1, x2)
json.dump(sim_rows, open(HERE / "sim_rows.json", "w"), indent=1, default=float)
def grad_method(key1, key2, sl):
    yy = []; xx = []
    for r in np.array(sim_rows)[sl]:
        a, b = r[key1], r[key2]
        if np.isfinite(a) and np.isfinite(b) and b > a:
            yy.append(np.log(b**2 - a**2)); xx.append(np.log(r["d"]))
    if len(xx) < 3: return np.nan, np.nan, len(xx)
    p, c = np.polyfit(xx, yy, 1, cov=True); return p[0], np.sqrt(c[0, 0]), len(xx)
methods = [("oldrect", "f1_old_rect", "f2_old_rect"), ("oldhann", "f1_old_hann", "f2_old_hann"),
           ("rect", "f1_rect", "f2_rect"), ("hann", "f1_hann", "f2_hann"), ("blackman", "f1_blackman", "f2_blackman"),
           ("normal", "f1_eta", "f2_xi"), ("fit", "f1_fit", "f2_fit"), ("true", "f1_true", "f2_true")]
grad_table = {}
for nm, k1, k2 in methods:
    grad_table[nm] = dict(zip(["g6", "e6", "n6"], grad_method(k1, k2, slice(0, 6))))
    grad_table[nm].update(dict(zip(["gall", "eall", "nall"], grad_method(k1, k2, slice(0, 10)))))
    print(f"   method {nm:8s}: 21-35 gradient {grad_table[nm]['g6']:.3f}±{grad_table[nm]['e6']:.3f} (n={grad_table[nm]['n6']}); all {grad_table[nm]['gall']:.3f}±{grad_table[nm]['eall']:.3f} (n={grad_table[nm]['nall']})")
json.dump(grad_table, open(HERE / "grad_table.json", "w"), indent=1, default=float)
for nm in ("oldrect", "oldhann", "rect", "hann", "blackman", "normal", "fit", "true"):
    put(f"simGrad{nm}", grad_table[nm]["g6"], "{:.2f}"); put(f"simGradErr{nm}", grad_table[nm]["e6"], "{:.2f}")
    put(f"simGradAll{nm}", grad_table[nm]["gall"], "{:.2f}"); put(f"simGradAllErr{nm}", grad_table[nm]["eall"], "{:.2f}")
    put(f"simN{nm}", grad_table[nm]["nall"], "{:d}")
put("simFzero", F0_S, "{:.3f}")

# figure: illustrative x(t) with beats and its spectrum (d = 25 mm, T = 30 s window of the 90 s record)
t, x1, x2 = keep25
fig, axs = plt.subplots(1, 2, figsize=(9.8, 3.4))
m = t < 30
axs[0].plot(t[m], x1[m]*1e3, lw=0.7, color="tab:blue", label="$x_1(t)$ (released spring)")
axs[0].plot(t[m], x2[m]*1e3, lw=0.7, color="tab:orange", alpha=0.8, label="$x_2(t)$ (other spring)")
axs[0].set_xlabel("$t$ / s"); axs[0].set_ylabel("displacement / mm"); axs[0].legend(loc="upper right"); axs[0].grid(alpha=0.3)
axs[0].set_title("beats: energy moves back and forth between the springs", fontsize=9)
pk, (fq, X) = peaks_two(x1[m], FPS, "rect", nfft=2**19)
axs[1].plot(fq, X, color="tab:red", lw=1); axs[1].set_xlim(0.5, 1.3)
for p_, lab in zip(pk, ["$f_1$", "$f_2$"]): axs[1].annotate(f"{lab} = {p_:.3f} Hz", xy=(p_, np.interp(p_, fq, X)), xytext=(p_+0.03, 0.9 if lab=="$f_1$" else 0.7), arrowprops=dict(arrowstyle="->"))
axs[1].set_xlabel("frequency / Hz"); axs[1].set_ylabel("|FFT($x_1$)| (normalised)"); axs[1].grid(alpha=0.3)
axs[1].set_title("amplitude spectrum of the first 30 s (rectangular window)", fontsize=9)
save(fig, "F_beats_example")

# figure: normal coordinates at d = 45 mm
t, x1, x2 = keep45
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
for ax, sig, lab, c in [(axs[0], x1, "$x_1$ alone", "tab:red"), (axs[1], x1+x2, r"$\eta=x_1+x_2$", "tab:blue"), (axs[2], x1-x2, r"$\xi=x_1-x_2$", "tab:green")]:
    pk, (fq, X) = peaks_two(sig, FPS, "rect", nfft=2**20, thr=0.5)
    ax.plot(fq, X, color=c, lw=1); ax.set_xlim(0.76, 0.86); ax.set_title(f"FFT of {lab}", fontsize=9); ax.grid(alpha=0.3)
    for fr in true_f[45]: ax.axvline(fr, color="k", ls=":", lw=0.8)
    ax.set_xlabel("frequency / Hz")
    if pk: ax.text(0.762, 0.9, "peaks: " + ", ".join(f"{p_:.4f}" for p_ in pk) + " Hz", fontsize=7)
axs[0].set_ylabel("|FFT| (normalised)")
fig.suptitle(f"SIMULATED, $d$ = 45 mm, $T$ = 90 s, $A$ = 3 mm: true $f_1$ = {true_f[45][0]:.4f}, $f_2$ = {true_f[45][1]:.4f} Hz (dotted)", fontsize=9)
save(fig, "F_normal_coords")

# figure: two-tone fit example (d = 35 mm) + residuals
t, x1, x2 = simulate(35, 3.0, seed=35)
f1f, f2f, e1f, e2f = two_tone_fit(t, x1, true_f[35][0]+0.002, true_f[35][1]-0.002)
def model_tt(p, t):
    A1, A2, ph1, ph2, fa_, fb_, ta, tb, c = p
    return A1*np.exp(-t/ta)*np.cos(2*np.pi*fa_*t+ph1) + A2*np.exp(-t/tb)*np.cos(2*np.pi*fb_*t+ph2) + c
p0 = [x1.std(), x1.std(), 0, 0, true_f[35][0]+0.002, true_f[35][1]-0.002, TAU, TAU, 0]
rr = least_squares(lambda p: model_tt(p, t) - x1, p0, x_scale="jac")
fig, axs = plt.subplots(2, 1, figsize=(9.5, 4.2), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
axs[0].plot(t, x1*1e3, lw=0.6, color="0.3", label="synthetic track $x_1(t)$ (0.1 mm rms tracking noise added)")
axs[0].plot(t, model_tt(rr.x, t)*1e3, lw=0.8, color="tab:red", alpha=0.8, label=f"two-damped-sinusoid fit: $f_1$ = {f1f:.4f} Hz (±{e1f*1e3:.2f} mHz), $f_2$ = {f2f:.4f} Hz (±{e2f*1e3:.2f} mHz)")
axs[0].legend(loc="upper right"); axs[0].set_ylabel("$x_1$ / mm"); axs[0].grid(alpha=0.3)
axs[1].plot(t, rr.fun*1e3, lw=0.5, color="tab:blue"); axs[1].set_ylabel("residual / mm"); axs[1].set_xlabel("$t$ / s"); axs[1].grid(alpha=0.3)
fig.suptitle(f"SIMULATED, $d$ = 35 mm: true $f_1$ = {true_f[35][0]:.4f}, $f_2$ = {true_f[35][1]:.4f} Hz", fontsize=9)
save(fig, "F_twotone_fit")
put("ttFone", f1f, "{:.4f}"); put("ttFtwo", f2f, "{:.4f}"); put("ttEone", e1f*1e3, "{:.1f}"); put("ttEtwo", e2f*1e3, "{:.1f}")
put("ttTrueOne", true_f[35][0], "{:.4f}"); put("ttTrueTwo", true_f[35][1], "{:.4f}")

# figure: simulated further data ln-ln by method
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for nm, k1, k2, mk, c in [("single-coordinate FFT, rectangular", "f1_rect", "f2_rect", "s", "tab:red"),
                          ("normal-coordinate FFT", "f1_eta", "f2_xi", "^", "tab:blue"),
                          ("two-tone time-domain fit", "f1_fit", "f2_fit", "o", "tab:green")]:
    xx, yy = [], []
    for r in sim_rows:
        a, b = r[k1], r[k2]
        if np.isfinite(a) and np.isfinite(b) and b > a: xx.append(np.log(r["d"])); yy.append(np.log(b**2-a**2))
    ax.plot(xx, yy, mk, color=c, ms=5, label=nm + f" ({len(xx)} points)")
dd = np.linspace(20, 62, 200)
ax.plot(np.log(dd), np.log(Kmm*dd**-5/(2*np.pi**2)), "k-", lw=1, label="simulation truth: dipole law, gradient −5")
xx, yy = [], []
for r in sim_rows:
    a, b = r["f1_old_rect"], r["f2_old_rect"]
    if np.isfinite(a) and np.isfinite(b) and b > a: xx.append(np.log(r["d"])); yy.append(np.log(b**2-a**2))
ax.plot(xx, yy, "x", color="0.4", ms=7, mew=1.5, label=f"original procedure ($T$ = 30 s, $A$ = 10 mm, $x_1$ only) ({len(xx)} points)")
ax.set_xlabel(r"$\ln(d/\mathrm{mm})$"); ax.set_ylabel(r"$\ln(f_2^2-f_1^2)$"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
ax.set_title("SIMULATED records: identical springs, dipole force, 240 fps", fontsize=9)
save(fig, "F_sim_lnln")

# figure: amplitude extrapolation from segments of an A = 10 mm run at d = 21 mm (SIMULATED)
print("9. segment analysis for amplitude extrapolation")
t, x1, x2 = simulate(21, 10.0, T=120, seed=21)
seg = 15.0; segs = []
for k in range(int(120/seg)):
    m = (t >= k*seg) & (t < (k+1)*seg)
    xi = x1[m] - x2[m]; a = 0.5*(xi.max() - xi.min())
    # frequency of xi in the segment from a single damped sinusoid fit
    def m1f(p, tt): return p[0]*np.exp(-tt/p[3])*np.cos(2*np.pi*p[1]*tt+p[2]) + p[4]
    p0 = [a, true_f[21][1], 0, TAU, 0]
    r_ = least_squares(lambda p: m1f(p, t[m]-t[m][0]) - xi, p0, x_scale="jac")
    segs.append((a*1e3, r_.x[1]))
segs = np.array(segs)
pa = np.polyfit(segs[:, 0]**2, segs[:, 1], 1)
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.plot(segs[:, 0]**2, segs[:, 1], "o", color="tab:red", label="15 s segments of one $A$ = 10 mm record")
aa = np.linspace(0, segs[:, 0].max()**2*1.05, 10)
ax.plot(aa, np.polyval(pa, aa), "k-", lw=1, label=f"linear in $a^2$: intercept $f_2(0)$ = {pa[1]:.4f} Hz")
ax.axhline(true_f[21][1], color="tab:green", ls="--", lw=1, label=f"small-amplitude (linear) value {true_f[21][1]:.4f} Hz")
ax.set_xlabel(r"segment amplitude squared $a^2$ / mm$^2$"); ax.set_ylabel("$f_2$ of the segment / Hz"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title("SIMULATED, $d$ = 21 mm: extrapolating $f_2$ to zero amplitude", fontsize=9)
save(fig, "F_extrapolation")
put("extrapFtwo", pa[1], "{:.4f}"); put("extrapTrue", true_f[21][1], "{:.4f}"); put("extrapFirst", segs[0, 1], "{:.4f}")

# =============================================================== 10. damping linewidth illustration + numbers
put("lwAmp", np.sqrt(3)/(np.pi*TAU), "{:.3f}"); put("lwPow", 1/(np.pi*TAU), "{:.3f}")
put("binT", 1/T_REC*1e3, "{:.1f}"); put("binTnew", 1/T_NEW*1e3, "{:.1f}")
put("Nsamples", int(T_REC*FPS), "{:d}")
put("Kmm", Kmm, "{:.3e}")
# tilt numbers
L_ASSUMED = 150.0
put("tiltRad", 3*A_REL/(2*L_ASSUMED), "{:.2f}"); put("tiltDeg", np.degrees(3*A_REL/(2*L_ASSUMED)), "{:.1f}")
put("tiltPct", 100*(1-np.cos(3*A_REL/(2*L_ASSUMED))), "{:.2f}")

# =============================================================== 11. LaTeX table fragments
TAB = HERE / "tables"; TAB.mkdir(exist_ok=True)
def fmt(v, f="{:.3f}"):
    return "--" if (v is None or (isinstance(v, float) and not np.isfinite(v))) else f.format(v)
# processed data table
with open(TAB / "tab_processed.tex", "w") as fh:
    for r in rows:
        fh.write(f"{r['d']:.0f} & {r['lnd']:.3f} & {r['dlnd']:.3f} & {r['f1']:.3f} & {r['f2']:.3f} & {r['Y']:.4f} & {r['dY']:.4f} & "
                 f"{r['lnY']:.3f} & {r['dlnY']:.3f} & {r['X5']/1e8:.3f} & {r['dX5']/1e8:.3f} \\\\\n")
# resolution / pulling table
with open(TAB / "tab_resolution.tex", "w") as fh:
    for r in pull_rows:
        fh.write(f"{r['d']:.0f} & {r['df']:.3f} & {r['bins30']:.1f} & {r['Tmin_rect']:.0f} & {r['Tmin_hann']:.0f} & "
                 f"{fmt(r['e1_30_rect'],'{:+.1f}')} & {fmt(r['e2_30_rect'],'{:+.1f}')} & {fmt(r['e1_30_hann'],'{:+.1f}')} & {fmt(r['e2_30_hann'],'{:+.1f}')} & "
                 f"{fmt(r['e1_90_rect'],'{:+.1f}')} & {fmt(r['e2_90_rect'],'{:+.1f}')} \\\\\n")
# amplitude table
with open(TAB / "tab_amplitude.tex", "w") as fh:
    for dsel in d_mm:
        fh.write(f"{dsel:.0f} & {f2lin[dsel]:.4f} & " + " & ".join(f"{100*(shift[(dsel,A)]/f2lin[dsel]-1):+.2f}" for A in amps) + f" & {LL(dsel,10):+.1f} \\\\\n")
# finite-size table
with open(TAB / "tab_finite.tex", "w") as fh:
    for (D, t), g in fs_grad.items():
        fh.write(f"{D} & {t} & {g:.2f} \\\\\n")
# simulated methods table
names = {"oldrect": "original procedure: $T=30$\\,s, $A=10$\\,mm, FFT of $x_1$, rectangular",
         "oldhann": "original procedure, Hann window",
         "rect": "new records ($T=90$\\,s, $A=3$\\,mm): FFT of $x_1$, rectangular",
         "hann": "new records: FFT of $x_1$, Hann", "blackman": "new records: FFT of $x_1$, Blackman",
         "normal": "new records: FFT of normal coordinates $\\eta$, $\\xi$",
         "fit": "new records: two-damped-sinusoid fit of $x_1(t)$", "true": "exact frequencies used to generate the records"}
with open(TAB / "tab_simmethods.tex", "w") as fh:
    for nm in ["oldrect", "oldhann", "rect", "hann", "blackman", "normal", "fit", "true"]:
        g = grad_table[nm]
        fh.write(f"{names[nm]} & {fmt(g['g6'],'{:.2f}')} $\\pm$ {fmt(g['e6'],'{:.2f}')} ({g['n6']}) & {fmt(g['gall'],'{:.2f}')} $\\pm$ {fmt(g['eall'],'{:.2f}')} ({g['nall']}) \\\\\n")
# appendix: simulated per-d frequencies
with open(TAB / "tab_simdata.tex", "w") as fh:
    for r in sim_rows:
        fh.write(f"{r['d']:.0f} & {r['f1_true']:.4f} & {r['f2_true']:.4f} & {fmt(r['f1_old_rect'],'{:.4f}')} & {fmt(r['f2_old_rect'],'{:.4f}')} & "
                 f"{fmt(r['f1_rect'],'{:.4f}')} & {fmt(r['f2_rect'],'{:.4f}')} & {fmt(r['f1_eta'],'{:.4f}')} & {fmt(r['f2_xi'],'{:.4f}')} & "
                 f"{fmt(r['f1_fit'],'{:.4f}')} & {fmt(r['f2_fit'],'{:.4f}')} \\\\\n")
# static deflection table
with open(TAB / "tab_deflection.tex", "w") as fh:
    for r in rows:
        fh.write(f"{r['d']:.0f} & {r['Y']:.4f} & {r['delta']:.2f} & {r['d']+2*r['delta']:.1f} & {100*2*r['delta']/r['d']:.1f} \\\\\n")

with open(HERE / "numbers.tex", "w") as fh:
    fh.write("% auto-generated by make_figures.py -- do not edit by hand\n")
    for k, v in NUM.items():
        fh.write(f"\\newcommand{{\\n{k.replace('_', '')}}}{{{v}}}\n")
json.dump(NUM, open(HERE / "numbers.json", "w"), indent=1)
print("done:", len(NUM), "numbers exported")
