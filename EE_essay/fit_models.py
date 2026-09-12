"""
fit_models.py -- weighted non-linear least-squares fits of four physical models to the
measured (d, f1, f2) data (Section 12.3 of the essay) + a few extra simulation numbers.
Writes numbers2.tex and tables/tab_modelfits.tex.
"""
import numpy as np, json
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from pathlib import Path
HERE = Path(__file__).parent
exec(open(HERE / "make_figures.py").read().split("# =============================================================== 1.")[0])  # reuse data + helpers
NUM2 = {}
def put2(k, v, f="{:.3g}"): NUM2[k] = f.format(v)

d6, F1, F2 = d_mm[six], f1[six], f2[six]
# ---------- finite-size force (same disc model as make_figures)
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

def modes_general(kap, fa, fb):
    wa2_ = (2*np.pi*fa)**2; wb2_ = (2*np.pi*fb)**2
    wbar2 = (wa2_ + wb2_)/2; Dl = (wb2_ - wa2_)/2
    return np.sqrt(wbar2 + kap - np.sqrt(Dl**2 + kap**2))/(2*np.pi), np.sqrt(wbar2 + kap + np.sqrt(Dl**2 + kap**2))/(2*np.pi)

fits = {}
# (a) identical springs, dipole n = 5 : params f0, K
def ra(p): 
    f0, K = p; a, b = modes_general(K*d6**-5, f0, f0); return np.concatenate([(a-F1)/DF, (b-F2)/DF])
r = least_squares(ra, [0.81, 4e7]); fits["a"] = dict(chi2=r.fun@r.fun, dof=12-2, p=r.x, note=f"$f_0={r.x[0]:.3f}$ Hz")
# (b) identical springs, free exponent n : params f0, K, n
def rb(p):
    f0, K, n = p; a, b = modes_general(K*d6**-n, f0, f0); return np.concatenate([(a-F1)/DF, (b-F2)/DF])
r = least_squares(rb, [0.81, 1e6, 3.7]); J = r.jac; cov = np.linalg.inv(J.T@J)*(r.fun@r.fun)/(12-3)
fits["b"] = dict(chi2=r.fun@r.fun, dof=9, p=r.x, note=f"$n={r.x[2]:.2f}\\pm{np.sqrt(cov[2,2]):.2f}$")
put2("nJoint", r.x[2], "{:.2f}"); put2("nJointErr", np.sqrt(cov[2,2]), "{:.2f}")
# (c) detuned, n = 5 : params fa, fb, K
def rc(p):
    fa, fb, K = p; a, b = modes_general(K*d6**-5, fa, fb); return np.concatenate([(a-F1)/DF, (b-F2)/DF])
r = least_squares(rc, [0.80, 0.82, 4e7]); J = r.jac; cov = np.linalg.inv(J.T@J)*(r.fun@r.fun)/(12-3)
fits["c"] = dict(chi2=r.fun@r.fun, dof=9, p=r.x, note=f"$f_a={r.x[0]:.3f}$, $f_b={r.x[1]:.3f}$ Hz")
# (d) finite-size force, identical springs : params f0, scale, D (t fixed 5 mm)
T_MAG = 5.0
cache = {}
def gam_vec(D):
    key = round(D, 2)
    if key not in cache: cache[key] = np.array([gamma_fs(x, D, T_MAG) for x in d6])
    return cache[key]
def rd(p):
    f0, S, D = p; D = abs(D); a, b = modes_general(S*gam_vec(D), f0, f0); return np.concatenate([(a-F1)/DF, (b-F2)/DF])
best = None
for D0 in [8, 12, 16, 20, 24]:
    r = least_squares(rd, [0.81, 1.0, D0], diff_step=[1e-4, 1e-4, 0.02])
    if best is None or r.fun@r.fun < best.fun@best.fun: best = r
r = best; J = r.jac; cov = np.linalg.inv(J.T@J)*(r.fun@r.fun)/(12-3)
fits["d"] = dict(chi2=r.fun@r.fun, dof=9, p=r.x, note=f"$D={abs(r.x[2]):.1f}\\pm{np.sqrt(cov[2,2]):.1f}$ mm ($t=5$ mm), $f_0={r.x[0]:.3f}$ Hz")
put2("Dfit", abs(r.x[2]), "{:.1f}"); put2("DfitErr", np.sqrt(cov[2,2]), "{:.1f}")
# (e) finite size + detuned: params fa, fb, S, D
def re_(p):
    fa, fb, S, D = p; D = abs(D); a, b = modes_general(S*gam_vec(D), fa, fb); return np.concatenate([(a-F1)/DF, (b-F2)/DF])
best = None
for D0 in [8, 14, 20]:
    r = least_squares(re_, [0.80, 0.82, 1.0, D0], diff_step=[1e-4, 1e-4, 1e-4, 0.02])
    if best is None or r.fun@r.fun < best.fun@best.fun: best = r
r = best; J = r.jac; cov = np.linalg.inv(J.T@J)*(r.fun@r.fun)/(12-4)
fits["e"] = dict(chi2=r.fun@r.fun, dof=8, p=r.x, note=f"$f_a={r.x[0]:.3f}$, $f_b={r.x[1]:.3f}$ Hz, $D={abs(r.x[3]):.1f}\\pm{np.sqrt(cov[3,3]):.1f}$ mm")
put2("DfitE", abs(r.x[3]), "{:.1f}"); put2("faE", r.x[0], "{:.3f}"); put2("fbE", r.x[1], "{:.3f}")
names = {"a": "(a) identical springs, point dipole ($n=5$)", "b": "(b) identical springs, free exponent $n$",
         "c": "(c) detuned springs, point dipole", "d": "(d) identical springs, finite-size magnets",
         "e": "(e) detuned springs + finite-size magnets"}
with open(HERE / "tables" / "tab_modelfits.tex", "w") as fh:
    for k in "abcde":
        f = fits[k]; fh.write(f"{names[k]} & {len(f['p'])} & {f['chi2']:.1f} & {f['chi2']/f['dof']:.1f} & {f['note']} \\\\\n")
        put2(f"chi{k}", f["chi2"], "{:.1f}"); put2(f"chired{k}", f["chi2"]/f["dof"], "{:.1f}")
        print(k, f)

# ---------- extra simulation numbers: pulling-only gradient at T=30, A=3; old-procedure f1 drift
Y = f2**2 - f1**2
F0_S = float(np.mean(f1[six])); K_s = Y[2] * 2 * np.pi**2 * d_mm[2]**5 * 1e-15
def rhs_full(t, s, d, wa2_, wb2_, K, tau):
    x1, v1, x2, v2 = s; r = d + x2 - x1; dF = K/4*(r**-4 - d**-4)
    return [v1, -wa2_*x1 - dF - (2/tau)*v1, v2, -wb2_*x2 + dF - (2/tau)*v2]
def peaks_two(x, fs, window, nfft=2**19, fmin=0.4, fmax=2.0, thr=0.03):
    n = len(x); x = x - x.mean()
    wv = {"rect": np.ones(n), "hann": np.hanning(n), "blackman": np.blackman(n)}[window]
    X = np.abs(np.fft.rfft(x*wv, nfft)); fq = np.fft.rfftfreq(nfft, 1/fs)
    m = (fq > fmin) & (fq < fmax); idx = np.where(m)[0]; Xm = X[idx].max(); pk = []
    for i in idx[1:-1]:
        if X[i] > X[i-1] and X[i] >= X[i+1] and X[i] > thr*Xm:
            a, b, c = X[i-1], X[i], X[i+1]; p = 0.5*(a-c)/(a-2*b+c); pk.append((fq[i]+p*(fq[1]-fq[0]), b))
    pk.sort(key=lambda z: -z[1]); return sorted(q[0] for q in pk[:2])
w2 = (2*np.pi*F0_S)**2
res = {}
for A, T in [(3.0, 30.0), (10.0, 30.0), (10.0, 90.0)]:
    xs, ys, f1s = [], [], []
    for dsel in d_mm[six]:
        sol = solve_ivp(rhs_full, (0, T), [A*1e-3, 0, 0, 0], args=(dsel*1e-3, w2, w2, K_s, TAU), rtol=1e-9, atol=1e-12, dense_output=True)
        t = np.arange(0, T, 1/FPS); s = sol.sol(t); x1 = s[0] + np.random.default_rng(int(dsel)).normal(0, 1e-4, len(t))
        pk = peaks_two(x1, FPS, "rect", nfft=2**20)
        xs.append(np.log(dsel)); ys.append(np.log(pk[1]**2 - pk[0]**2)); f1s.append(pk[0])
    g = np.polyfit(xs, ys, 1)[0]; res[(A, T)] = (g, f1s)
    print(f"A={A} T={T}: gradient {g:.3f}; f1 recovered {np.round(f1s,4)}")
put2("gradPullOnly", res[(3.0, 30.0)][0], "{:.2f}")
put2("gradOldSim", res[(10.0, 30.0)][0], "{:.2f}")
put2("gradAtenTninety", res[(10.0, 90.0)][0], "{:.2f}")
put2("oldFoneSmall", res[(10.0, 30.0)][1][0], "{:.4f}"); put2("oldFoneLarge", res[(10.0, 30.0)][1][5], "{:.4f}")
put2("oldFoneDrift", (res[(10.0, 30.0)][1][0]-res[(10.0, 30.0)][1][5])*1e3, "{:.1f}")

with open(HERE / "numbers2.tex", "w") as fh:
    fh.write("% auto-generated by fit_models.py\n")
    for k, v in NUM2.items(): fh.write(f"\\newcommand{{\\n{k}}}{{{v}}}\n")
json.dump(NUM2, open(HERE / "numbers2.json", "w"), indent=1)
print(NUM2)
