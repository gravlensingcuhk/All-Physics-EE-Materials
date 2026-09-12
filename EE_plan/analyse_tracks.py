"""
analyse_tracks.py  -  re-analysis of your raw video-tracking time series (x(t) of ONE magnet)

For every file it does the four things the "Further methodology" section needs:
  1. FFT peak-picking with THREE windows (rectangular / Hann / Blackman), zero-padding x16
     and parabolic (3-bin) peak interpolation             -> f1, f2 per window
  2. Time-domain fit of TWO damped sinusoids (non-linear least squares, seeded by the FFT)
     -> f1, f2 with statistical uncertainties (beats the 1/T Fourier limit)
  3. Segment analysis: f2 in the early / middle / late third of the record vs amplitude
     -> amplitude dependence caused by the non-linear r^-4 force
  4. A summary CSV + one diagnostic figure per file, and the final ln(f2^2-f1^2) vs ln(d) fits.

USAGE
  python analyse_tracks.py "tracks/*.csv" --out results
  python analyse_tracks.py "tracks/*.txt" --tcol t --xcol x --fmin 0.3 --fmax 3 --out results

FILE NAMING: put the magnet separation in the filename, e.g.  d21_trial1.csv, d28.0_t2.txt
FILE FORMAT: anything Tracker exports (tab/comma separated, header lines are skipped).
             First numeric column = t, second = x unless --tcol/--xcol are given.
"""
import argparse, glob, re, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ----------------------------------------------------------------------------- loading
def load_track(path, tcol=None, xcol=None):
    rows, header = [], None
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = re.split(r"[,\t; ]+", line.strip().replace(",", "."))
            parts = [p for p in parts if p != ""]
            if not parts:
                continue
            try:
                vals = [float(p) for p in parts]
                rows.append(vals)
            except ValueError:
                if header is None and len(parts) >= 2:
                    header = parts
    if not rows:
        raise ValueError(f"no numeric data in {path}")
    n = min(len(r) for r in rows)
    arr = np.array([r[:n] for r in rows])
    ti, xi = 0, 1
    if header is not None and tcol in header and xcol in header:
        ti, xi = header.index(tcol), header.index(xcol)
    t, x = arr[:, ti], arr[:, xi]
    ok = np.isfinite(t) & np.isfinite(x)
    return t[ok], x[ok]

def resample_uniform(t, x):
    """Video frames are not always uniformly timed; FFT needs uniform sampling."""
    dt = np.diff(t)
    jitter = np.std(dt) / np.mean(dt)
    tu = np.arange(t[0], t[-1], np.median(dt))
    return tu, np.interp(tu, t, x), jitter

# ----------------------------------------------------------------------------- FFT
WINDOWS = {"rect": np.ones, "hann": np.hanning, "blackman": np.blackman}

def fft_peaks(t, x, window="hann", fmin=0.3, fmax=3.0, pad=16, npeaks=2, rel=0.05):
    fs = 1 / (t[1] - t[0])
    x = x - np.polyval(np.polyfit(t, x, 1), t)              # remove offset + drift
    w = WINDOWS[window](len(x))
    nfft = int(2 ** np.ceil(np.log2(len(x) * pad)))
    X = np.abs(np.fft.rfft(x * w, nfft))
    f = np.fft.rfftfreq(nfft, 1 / fs)
    m = (f >= fmin) & (f <= fmax)
    X = X / X[m].max()
    idx = np.where(m)[0]
    peaks = []
    for i in idx[1:-1]:
        if X[i] > X[i - 1] and X[i] >= X[i + 1] and X[i] > rel:
            a, b, c = X[i - 1], X[i], X[i + 1]
            p = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
            peaks.append((f[i] + p * (f[1] - f[0]), b))
    peaks.sort(key=lambda q: -q[1])
    fr = sorted(q[0] for q in peaks[:npeaks])
    return fr, f[m], X[m], 1 / (t[-1] - t[0])

# ----------------------------------------------------------------------------- time-domain fit
def two_tone(t, A1, tau1, f1, p1, A2, tau2, f2, p2, c):
    return (A1 * np.exp(-t / tau1) * np.cos(2 * np.pi * f1 * t + p1)
            + A2 * np.exp(-t / tau2) * np.cos(2 * np.pi * f2 * t + p2) + c)

def fit_two_tones(t, x, f_guess):
    t0 = t - t[0]
    A0 = (x.max() - x.min()) / 2
    T = t0[-1]
    best = None
    for tau0 in (T, T / 2, 5 * T):
        for p in (0.0, np.pi / 2):
            p0 = [A0 / 2, tau0, f_guess[0], p, A0 / 2, tau0, f_guess[1], p, x.mean()]
            lb = [0, 0.05 * T, 0.5 * f_guess[0], -2 * np.pi, 0, 0.05 * T, 0.5 * f_guess[1], -2 * np.pi, -np.inf]
            ub = [np.inf, 100 * T, 1.5 * f_guess[0], 2 * np.pi, np.inf, 100 * T, 1.5 * f_guess[1], 2 * np.pi, np.inf]
            try:
                popt, pcov = curve_fit(two_tone, t0, x, p0=p0, bounds=(lb, ub), maxfev=20000)
                rss = np.sum((two_tone(t0, *popt) - x) ** 2)
                if best is None or rss < best[0]:
                    best = (rss, popt, pcov)
            except Exception:
                continue
    if best is None:
        return None
    rss, popt, pcov = best
    perr = np.sqrt(np.diag(pcov))
    fa, fb = sorted([(popt[2], perr[2]), (popt[6], perr[6])])
    r2 = 1 - rss / np.sum((x - x.mean()) ** 2)
    return dict(f1=fa[0], df1=fa[1], f2=fb[0], df2=fb[1], tau=(popt[1], popt[5]), r2=r2, popt=popt)

# ----------------------------------------------------------------------------- segments
def segment_analysis(t, x, nseg=3, fmin=0.3, fmax=3.0):
    out = []
    edges = np.linspace(0, len(t), nseg + 1).astype(int)
    for a, b in zip(edges[:-1], edges[1:]):
        ts, xs = t[a:b], x[a:b]
        fr, *_ = fft_peaks(ts, xs, "hann", fmin, fmax)
        amp = np.sqrt(2) * np.std(xs - xs.mean())
        out.append((0.5 * (ts[0] + ts[-1]), amp, fr[-1] if fr else np.nan, fr[0] if fr else np.nan))
    return out

# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern")
    ap.add_argument("--tcol", default="t"); ap.add_argument("--xcol", default="x")
    ap.add_argument("--fmin", type=float, default=0.3); ap.add_argument("--fmax", type=float, default=3.0)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    files = sorted(glob.glob(a.pattern))
    if not files:
        sys.exit("no files matched " + a.pattern)

    rows = []
    for path in files:
        name = os.path.basename(path)
        mo = re.search(r"d(\d+(?:\.\d+)?)", name)
        d = float(mo.group(1)) if mo else np.nan
        t, x = load_track(path, a.tcol, a.xcol)
        t, x, jitter = resample_uniform(t, x)
        T = t[-1] - t[0]
        res = {"file": name, "d_mm": d, "T_s": T, "fs_Hz": 1 / (t[1] - t[0]), "frame_jitter": jitter}

        fig, axs = plt.subplots(1, 3, figsize=(15, 4))
        axs[0].plot(t, x, lw=0.6); axs[0].set_title(f"{name}: x(t), T = {T:.1f} s"); axs[0].set_xlabel("t / s")
        for win, c in [("rect", "tab:red"), ("hann", "tab:blue"), ("blackman", "tab:green")]:
            fr, f, X, binw = fft_peaks(t, x, win, a.fmin, a.fmax)
            res[f"f1_{win}"] = fr[0] if len(fr) == 2 else np.nan
            res[f"f2_{win}"] = fr[1] if len(fr) == 2 else (fr[0] if fr else np.nan)
            axs[1].plot(f, X, c, lw=1, label=f"{win}: " + (", ".join(f"{q:.4f}" for q in fr) if fr else "no peak"))
        res["bin_Hz"] = binw
        axs[1].set_xlim(a.fmin, a.fmax); axs[1].set_xlabel("f / Hz"); axs[1].set_title(f"FFT (bin = {binw:.4f} Hz, zero-padded x16)")
        axs[1].legend(fontsize=7)

        guess = [res["f1_hann"], res["f2_hann"]]
        if not np.isfinite(guess[0]):
            guess = [res["f2_hann"] * 0.9, res["f2_hann"]]
        fit = fit_two_tones(t, x, guess)
        if fit:
            res.update(f1_fit=fit["f1"], df1_fit=fit["df1"], f2_fit=fit["f2"], df2_fit=fit["df2"], fit_r2=fit["r2"])
            axs[0].plot(t, two_tone(t - t[0], *fit["popt"]), "r", lw=0.6, alpha=0.7,
                        label=f"2-tone fit: f1={fit['f1']:.4f}±{fit['df1']:.4f}, f2={fit['f2']:.4f}±{fit['df2']:.4f} Hz (R²={fit['r2']:.3f})")
            axs[0].legend(fontsize=7)

        seg = segment_analysis(t, x, 3, a.fmin, a.fmax)
        for i, (tc, amp, f2s, f1s) in enumerate(seg):
            res[f"seg{i+1}_amp"] = amp; res[f"seg{i+1}_f2"] = f2s; res[f"seg{i+1}_f1"] = f1s
        axs[2].plot([s[1] for s in seg], [s[2] for s in seg], "o-", color="tab:red", label="f2")
        axs[2].set_xlabel("amplitude in segment (tracker units)"); axs[2].set_ylabel("f2 / Hz")
        axs[2].set_title("f2 vs amplitude (early → late thirds)"); axs[2].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, name + ".png"), dpi=110); plt.close(fig)
        rows.append(res)
        print(f"{name:28s} T={T:6.1f}s bin={binw:.4f}Hz  hann: {res['f1_hann']:.4f}/{res['f2_hann']:.4f}  "
              f"rect: {res['f1_rect']:.4f}/{res['f2_rect']:.4f}  fit: "
              + (f"{fit['f1']:.4f}/{fit['f2']:.4f}" if fit else "failed"))

    # ------------------------------------------------------------- summary table
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "file", k != "d_mm", k))
    with open(os.path.join(a.out, "summary.csv"), "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r.get(k, "")) for k in keys) + "\n")

    # ------------------------------------------------------------- ln-ln fits per method
    ds = np.array([r["d_mm"] for r in rows])
    if np.all(np.isfinite(ds)) and len(set(ds)) >= 3:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        print("\nln(f2^2 - f1^2) vs ln(d)  (trials at the same d are averaged):")
        for meth, c in [("rect", "tab:red"), ("hann", "tab:blue"), ("blackman", "tab:green"), ("fit", "k")]:
            xs, ys = [], []
            for dv in sorted(set(ds)):
                sel = [r for r in rows if r["d_mm"] == dv and np.isfinite(r.get(f"f1_{meth}", np.nan))]
                if not sel:
                    continue
                f1 = np.mean([r[f"f1_{meth}"] for r in sel]); f2 = np.mean([r[f"f2_{meth}"] for r in sel])
                if f2 > f1:
                    xs.append(np.log(dv)); ys.append(np.log(f2 ** 2 - f1 ** 2))
            if len(xs) >= 3:
                p, cov = np.polyfit(xs, ys, 1, cov=True)
                ax.plot(xs, ys, "o", color=c)
                xx = np.linspace(min(xs), max(xs), 10)
                ax.plot(xx, np.polyval(p, xx), "-", color=c, lw=1, label=f"{meth}: gradient {p[0]:.2f} ± {np.sqrt(cov[0,0]):.2f}")
                print(f"  {meth:9s} gradient = {p[0]:.3f} ± {np.sqrt(cov[0,0]):.3f}, intercept = {p[1]:.3f} ± {np.sqrt(cov[1,1]):.3f}  (n = {len(xs)})")
        ax.set_xlabel("ln(d / mm)"); ax.set_ylabel("ln[(f2² − f1²) / Hz²]"); ax.legend(fontsize=8)
        ax.set_title("Power-law test with each analysis method (dipole theory: gradient −5)")
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "lnln_by_method.png"), dpi=120)
    print("\nwritten to", a.out)

if __name__ == "__main__":
    main()
