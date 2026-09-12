"""
plot_beats.py -- annotated beat figure (and FFT) from Tracker exports.

usage:
    python3 plot_beats.py springA.txt                      # one spring
    python3 plot_beats.py springA.txt springB.txt          # both springs (complementary envelopes)
    options: --tmin 54 --tmax 98   (time window to plot)   --out figs/beats_measured.pdf

Tracker export: File > Export > Data (or copy the t, x columns into a text file).
The loader accepts tab / comma / space separated files with any number of header lines;
the first two numeric columns are taken as t (s) and x (m or mm -- see --unit).
"""
import argparse, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert, find_peaks
from scipy.optimize import curve_fit

# ------------------------------------------------------------------ loading
def load_tracker(fn, col=1):
    """Read a Tracker export (t in column 0, x in column `col`; header lines and blank cells are skipped).
    A file name may carry the column: 'both_masses.txt:3' reads x from column 3 (t, xA, yA, xB, yB -> xB = 3)."""
    if ":" in fn and fn.rsplit(":", 1)[1].isdigit():
        fn, col = fn.rsplit(":", 1)[0], int(fn.rsplit(":", 1)[1])
    t, x = [], []
    for line in open(fn, encoding="utf-8", errors="ignore"):
        parts = [p for p in line.replace(",", " ").replace("\t", " ").replace(";", " ").split() if p]
        if len(parts) < col + 1:
            continue
        try:
            a, b = float(parts[0]), float(parts[col])
        except ValueError:
            continue                      # header line or empty cell
        t.append(a); x.append(b)
    t, x = np.array(t), np.array(x)
    # uniform resampling (Tracker frame times can jitter / drop frames)
    dt = np.median(np.diff(t))
    if np.max(np.diff(t)) < 1.5 * dt: dt = (t[-1] - t[0]) / (len(t) - 1)   # no dropped frames: exact mean spacing (Tracker rounds t to 1 ms)
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, x), 1.0 / dt

def release_index(t, x, fs):
    """Index of the release: end of the flat 'hold' segment at the start of the record (if there is one).
    The hold level is the median of the first 0.1 s; the motion is taken to start where |x - hold| first
    exceeds 25 % of the half range, backed up to the last sample within 5 % of the hold level."""
    n0 = max(3, int(0.1 * fs)); x_hold = np.median(x[:n0]); rng = 0.5 * (x.max() - x.min())
    if rng <= 0: return 0
    beyond = np.where(np.abs(x - x_hold) > 0.25 * rng)[0]
    if len(beyond) == 0: return 0
    j = beyond[0]
    while j > 0 and abs(x[j] - x_hold) > 0.05 * rng: j -= 1
    return j

# ------------------------------------------------------------------ spectrum
def make_window(n, name):
    name = (name or "rect").lower()
    if name in ("hann", "hanning"): return np.hanning(n)
    if name == "hamming": return np.hamming(n)
    if name == "blackman": return np.blackman(n)
    if name == "flattop":
        k = np.arange(n) / (n - 1); return (0.21557895 - 0.41663158 * np.cos(2 * np.pi * k) + 0.277263158 * np.cos(4 * np.pi * k)
                                            - 0.083578947 * np.cos(6 * np.pi * k) + 0.006947368 * np.cos(8 * np.pi * k))
    return np.ones(n)                                          # rect / boxcar

def spectrum(x, fs, nfft=2**20, window=None):
    w = make_window(len(x), window)
    X = np.abs(np.fft.rfft((x - x.mean()) * w, nfft)); f = np.fft.rfftfreq(nfft, 1 / fs)
    return f, X / X.max()

def peak_ratio(x, fs):
    """|X(f2)|/|X(f1)| from a Hann-windowed spectrum: the Hann side lobes are negligible 18 bins away, so the
    height of one peak is not biased by the other peak's leakage (with a rectangular window this bias is ~2-3 %
    and of opposite sign for the two springs, which would fake a small detuning)."""
    fh, Xh = spectrum(x, fs, window="hann"); p = two_peaks(fh, Xh); return p[1][1] / p[0][1], p

def two_peaks(f, X, fmin=3.0, fmax=15.0):
    m = (f > fmin) & (f < fmax); idx = np.where(m)[0]
    pk, _ = find_peaks(X[idx], prominence=0.05)
    pk = idx[pk]; pk = sorted(pk, key=lambda i: -X[i])[:2]
    out = []
    for i in pk:                                   # parabolic refinement
        a, b, c = X[i - 1], X[i], X[i + 1]; p = 0.5 * (a - c) / (a - 2 * b + c)
        out.append((f[i] + p * (f[1] - f[0]), b))
    out.sort(); return out

# ------------------------------------------------------------------ main
ap = argparse.ArgumentParser()
ap.add_argument("files", nargs="+"); ap.add_argument("--tmin", type=float); ap.add_argument("--tmax", type=float)
ap.add_argument("--unit", default="m", help="unit of the x column in the file: m or mm")
ap.add_argument("--fmin", type=float, default=3.0); ap.add_argument("--fmax", type=float, default=15.0)
ap.add_argument("--no-autotrim", action="store_true", help="keep the flat hold segment at the start (default: it is detected and removed)")
ap.add_argument("--window", default="rect", choices=["rect", "hann", "hamming", "blackman"], help="FFT window for the peak POSITIONS (peak heights always use Hann)")
ap.add_argument("--compare-windows", action="store_true", help="also print f1, f2 for rect/hann/hamming/blackman (window-choice uncertainty)")
ap.add_argument("--skip", type=float, default=0.0, help="extra seconds to discard after the detected release (e.g. 0.3 to drop the first two swings)")
ap.add_argument("--csv", default=None, help="append one summary row per run to this csv (for the results table)")
ap.add_argument("--d", default="", help="separation label written to the csv, e.g. 24.5"); ap.add_argument("--trial", default="", help="trial label written to the csv")
ap.add_argument("--out", default="figs/beats_measured.pdf"); ap.add_argument("--labels", default="spring 1,spring 2")
ap.add_argument("--panels", default="all", choices=["all", "time", "fft"], help="all = trace(+zoom)+spectrum; time = trace(+zoom) only (e.g. for the Introduction); fft = spectrum only")
ap.add_argument("--titles", default=None, help="three titles separated by | (trace | zoom | spectrum); {s1}/{s2} are replaced by the labels")
ap.add_argument("--tavg", action="store_true", help="also mark T_avg on the trace (only meaningful if the two FFT peaks have equal height)")
args = ap.parse_args()
labels = args.labels.split(",")
_default_titles = ("Displacement\u2013time graph for {s1} in one of the trials|Zoomed plot of the beating phenomenon|"
                   + ("FFT spectrum for both springs in one of the trials" if len(args.files) > 1 else "FFT spectrum for {s1} in one of the trials"))
titles = [tt.format(s1=labels[0], s2=labels[1] if len(labels) > 1 else "") for tt in (args.titles or _default_titles).split("|")]
titles += [""] * (3 - len(titles))
scale = 1e3 if args.unit == "m" else 1.0          # plot in mm

data = []; t_release = None
for fn in args.files:
    t, x, fs = load_tracker(fn); x = x * scale
    if args.tmin is not None: m = t >= args.tmin; t, x = t[m], x[m]
    if args.tmax is not None: m = t <= args.tmax; t, x = t[m], x[m]
    if not args.no_autotrim and args.tmin is None:
        if t_release is None:                       # detect on the first file (the released spring); reuse for the others
            j = release_index(t, x, fs); t_release = t[j] + args.skip
            print(f"auto-trim: hold segment {t[0]:.2f}-{t[j]:.2f} s removed, analysis starts at t = {t_release:.2f} s")
        m = t >= t_release; t, x = t[m], x[m]
    data.append((t, x, fs))
if len(data) > 1:                                    # autotracked files end at different frames: cut all to the common span
    t_end = min(dd[0][-1] for dd in data); t_start = max(dd[0][0] for dd in data)
    if any(abs(dd[0][-1] - t_end) > 0.5 / dd[2] or abs(dd[0][0] - t_start) > 0.5 / dd[2] for dd in data):
        print(f"common span: files end at {', '.join(f'{dd[0][-1]:.2f}' for dd in data)} s -> all cut to {t_start:.2f}-{t_end:.2f} s")
    data = [(tt[(tt >= t_start - 1e-9) & (tt <= t_end + 1e-9)], xx[(tt >= t_start - 1e-9) & (tt <= t_end + 1e-9)], ff) for tt, xx, ff in data]
data = [(tt, xx - np.polyval(np.polyfit(tt, xx, 1), tt), ff) for tt, xx, ff in data]   # remove the mean (equilibrium position) and any slow drift

t, x, fs = data[0]
f, X = spectrum(x, fs, window=args.window); pk = two_peaks(f, X, args.fmin, args.fmax)
f1, f2 = pk[0][0], pk[1][0]
f_avg, T_beat = (f1 + f2) / 2, 1 / (f2 - f1)
env = np.abs(hilbert(x))
k = max(1, int(0.2 * fs / f_avg)); env_s = np.convolve(env, np.ones(k) / k, mode="same")   # light smoothing
edge = int(fs / f_avg); env_s[:edge] = np.nan; env_s[-edge:] = np.nan                        # hide Hilbert edge artefacts
env_s_full = np.where(np.isnan(env_s), np.nanmax(env_s), env_s)
hi = np.nanmax(env_s)
imax, _ = find_peaks(np.nan_to_num(env_s, nan=0), distance=int(0.6 * T_beat * fs), height=0.5 * hi, prominence=0.15 * hi)
imin, _ = find_peaks(-np.nan_to_num(env_s, nan=hi), distance=int(0.6 * T_beat * fs), prominence=0.15 * hi)
imax = imax[(imax > 1.5 * edge) & (imax < len(t) - 1.5 * edge)]; imin = imin[(imin > 1.5 * edge) & (imin < len(t) - 1.5 * edge)]
ratio = np.median(env_s[imin]) / np.median(env_s[imax]) if len(imin) and len(imax) else np.nan
tau = np.nan
if len(imax) >= 3:
    try:
        p, _ = curve_fit(lambda tt, a, tau: a * np.exp(-tt / tau), t[imax] - t[imax][0], env_s[imax], p0=[env_s[imax][0], 100])
        tau = p[1]
    except Exception:
        pass

# ---- beat period from ALL envelope maxima: straight-line fit t_n = t_0 + n T_beat (much more precise than one interval)
Nb = len(imax) - 1
if Nb >= 2:
    n_idx = np.arange(len(imax)); (T_beat_meas, t0_fit), cov_b = np.polyfit(n_idx, t[imax], 1, cov=True)
    dT_beat_meas = np.sqrt(cov_b[0, 0])
elif Nb == 1:
    T_beat_meas, dT_beat_meas = t[imax[1]] - t[imax[0]], np.nan
else:
    T_beat_meas, dT_beat_meas = np.nan, np.nan
# ---- average period: peaks of the fast oscillation inside ONE beat lobe (away from the nodes, where the
#      fast oscillation slips by half a cycle); up to 5 periods, centred on the 2nd envelope maximum
n_avg, T_avg_meas, sel = 0, np.nan, []
if len(imax) >= 2:
    c = imax[1]; half = int(0.36 * T_beat * fs)                  # lobe spans +-0.5 T_beat; stay inside +-0.36
    seg = slice(max(0, c - half), min(len(t), c + half))
    ip, _ = find_peaks(x[seg], distance=int(0.6 / f_avg * fs), height=0.3 * env_s_full[c]); ip = ip + seg.start
    if len(ip) >= 3:
        j = np.argmin(abs(ip - c)); n_side = min(j, len(ip) - 1 - j, 2)   # symmetric about the maximum, at most 5 periods
        sel = ip[j - n_side: j + n_side + 1]; n_avg = len(sel) - 1
        T_avg_meas = (t[sel[-1]] - t[sel[0]]) / n_avg
T_rec = t[-1] - t[0]; Nsamp = len(t)
print(f"record: {t[0]:.1f}-{t[-1]:.1f} s, T = {T_rec:.1f} s, fs = {fs:.1f} Hz, N = {Nsamp} samples, bin width 1/T = {1e3/T_rec:.1f} mHz, "
      f"zero-padded to 2^20 (spacing {fs/2**20*1e3:.3f} mHz), {args.window} window, linear trend removed")
if args.compare_windows:
    print("window comparison (peak positions):")
    for wn in ["rect", "hann", "hamming", "blackman"]:
        fw, Xw = spectrum(x, fs, window=wn); pw = two_peaks(fw, Xw, args.fmin, args.fmax)
        print(f"   {wn:9s} f1 = {pw[0][0]:.4f} Hz  f2 = {pw[1][0]:.4f} Hz  f2-f1 = {pw[1][0]-pw[0][0]:.4f} Hz  f2^2-f1^2 = {pw[1][0]**2-pw[0][0]**2:.4f} Hz^2")
print(f"FFT:   f1 = {f1:.4f} Hz, f2 = {f2:.4f} Hz  ->  f_avg = {f_avg:.4f} Hz, f_beat = f2 - f1 = {f2-f1:.4f} Hz  ->  T_avg = {1/f_avg:.3f} s, T_beat = 1/f_beat = {T_beat:.2f} s")
print(f"time domain: straight-line fit to {len(imax)} envelope maxima ({t[imax[0]]:.1f}-{t[imax[-1]]:.1f} s, {Nb} beats) -> T_beat = {T_beat_meas:.3f} +- {dT_beat_meas:.3f} s, f_beat = 1/T_beat = {1/T_beat_meas:.4f} +- {dT_beat_meas/T_beat_meas**2:.4f} Hz  (= f2 - f1 from the FFT: {f2-f1:.4f} Hz)")
if len(imax) >= 6:                                   # drift check: a beat period that shortens/lengthens as the amplitude decays = amplitude-dependent f2 (non-linear coupling)
    h = len(imax) // 2; Ta = np.polyfit(np.arange(h), t[imax[:h]], 1)[0]; Tb = np.polyfit(np.arange(len(imax) - h), t[imax[h:]], 1)[0]
    print(f"             first half of record: T_beat = {Ta:.3f} s; second half: {Tb:.3f} s  (a difference beyond the fit error = amplitude-dependent frequency, i.e. non-linearity)")
print(f"             {n_avg} fast periods inside one beat lobe -> T_avg = {T_avg_meas:.3f} s, f_avg = {1/T_avg_meas:.4f} Hz")
r_pk, _ = peak_ratio(x, fs); r_hold = (f1 / f2) ** 2
print(f"peak-height ratio r1 = |X(f2)|/|X(f1)| = {r_pk:.3f} for {labels[0]} (Hann-windowed spectrum).  Linear theory: 1.00 for an instantaneous release (x2(0) = 0);"
      f" (f1/f2)^2 = {r_hold:.3f} for a slow hold-and-release (spring 2 already sitting at its displaced equilibrium x2(0) = A*gamma/(k+gamma))")
print(f"envelope min/max ratio = {ratio:.3f}; predicted from r1: (1 - r1)/(1 + r1) = {(1 - r_pk) / (1 + r_pk):.3f}")
if len(data) > 1:
    r2, pkb = peak_ratio(data[1][1], data[1][2])
    rho = np.sqrt(r2 / r_pk); dg = (rho - 1 / rho) / 2; dfsq = dg / np.sqrt(1 + dg ** 2) * (f2 ** 2 - f1 ** 2)
    print(f"peak-height ratio r2 = {r2:.3f} for {labels[1]} (peaks at {pkb[0][0]:.4f}, {pkb[1][0]:.4f} Hz).  Identical oscillators => r2 = r1 for ANY initial condition."
          f"\n   sqrt(r1 r2) = {np.sqrt(r_pk * r2):.3f} = mode-amplitude ratio |c2/c1| (compare with (f1/f2)^2 = {r_hold:.3f} for a hold-and-release);"
          f"\n   r2/r1 = (a1/b1)^2 = {r2 / r_pk:.3f} => in-phase mode shape a1/b1 = {rho:.3f}, Delta/gamma = {dg:+.4f},"
          f" fb^2 - fa^2 = {dfsq:+.4f} Hz^2, fb - fa ~ {dfsq / (2 * f1):+.4f} Hz  ({'spring 2' if dg > 0 else 'spring 1'} is the stiffer oscillator; rough: +-0.005 Hz)")
print("NOTE: the local period at an envelope maximum equals 1/f_avg only when the two modes have EQUAL amplitude;"
      f" with r = {r_pk:.2f} it is 1/[f1 + (f2 - f1) r/(1 + r)] = {1/(f1+(f2-f1)*r_pk/(1+r_pk)):.3f} s, so f_avg is read from the FFT, not the trace (use --tavg to mark it anyway).")
print(f"decay time of the envelope maxima tau = {tau:.0f} s  ->  amplitude linewidth sqrt(3)/(pi tau) = {np.sqrt(3)/(np.pi*tau)*1e3:.1f} mHz")

# ------------------------------------------------------------------ figure
show_time, show_fft = args.panels in ("all", "time"), args.panels in ("all", "fft")
nrow = (1 + (len(data) > 1)) * show_time + show_fft
hr = ([1.3] + [1] * (len(data) > 1)) * show_time + [1] * show_fft
fig, axs = plt.subplots(nrow, 1, figsize=(9.5, 2.9 * nrow), gridspec_kw=dict(height_ratios=hr), squeeze=False); axs = axs[:, 0]
if not show_time: axs = [None] * (1 + (len(data) > 1)) + list(axs)
ax = axs[0] if show_time else plt.figure().add_subplot()   # dummy axis when the trace is not drawn
ax.set_title(titles[0], fontsize=10)
ax.plot(t, x, lw=0.7, color="0.25", label=f"$x(t)$, {labels[0]} (Tracker)")
ax.plot(t, env_s, "--", color="tab:red", lw=1.2, label="envelope")
ax.plot(t, -env_s, "--", color="tab:red", lw=1.2)
# beat period: arrow across up to 5 beats (envelope maximum to envelope maximum)
if Nb >= 1:
    nb = min(5, Nb); i0, i1 = imax[0], imax[nb]; y = 1.12 * np.nanmax(env_s)
    ax.annotate("", xy=(t[i1], y), xytext=(t[i0], y), arrowprops=dict(arrowstyle="<->", color="tab:red", lw=1.3))
    for i in (i0, i1): ax.plot([t[i], t[i]], [env_s_full[i], y], ":", color="tab:red", lw=0.8)
    ax.text((t[i0] + t[i1]) / 2, y * 1.04, rf"${nb}\,T_{{\rm beat}}$", ha="center", va="bottom", color="tab:red", fontsize=9)
    if Nb >= 2:
        ax.text(0.01, 0.03, rf"$T_{{\rm beat}}={T_beat_meas:.3f}\pm{dT_beat_meas:.3f}$ s (fit to {len(imax)} envelope maxima) "
                rf"$\Rightarrow f_{{\rm beat}}=1/T_{{\rm beat}}={1/T_beat_meas:.3f}$ Hz", transform=ax.transAxes, ha="left", va="bottom", color="tab:red", fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
# average period: n consecutive periods inside one beat lobe (optional, see note printed below)
if args.tavg and n_avg >= 1:
    y2 = -1.12 * np.nanmax(env_s)
    ax.annotate("", xy=(t[sel[-1]], y2), xytext=(t[sel[0]], y2), arrowprops=dict(arrowstyle="<->", color="tab:blue", lw=1.3))
    for i in sel: ax.plot([t[i], t[i]], [x[i], y2], ":", color="tab:blue", lw=0.6)
    ax.text((t[sel[0]] + t[sel[-1]]) / 2, y2 * 1.04,
            rf"${n_avg}\,T_{{\rm avg}}={t[sel[-1]]-t[sel[0]]:.2f}$ s $\Rightarrow T_{{\rm avg}}={T_avg_meas:.2f}$ s",
            ha="center", va="top", color="tab:blue", fontsize=9)
ax.set_ylim(-1.45 * np.nanmax(env_s), 1.45 * np.nanmax(env_s)); ax.set_ylabel("displacement / mm"); ax.set_xlabel("$t$ / s")
ax.legend(loc="upper right", fontsize=8, ncol=2); ax.grid(alpha=0.3)

if len(data) > 1 and show_time:                     # second spring: complementary envelope
    tb, xb, fsb = data[1]; envb = np.abs(hilbert(xb))
    ax = axs[1]
    tz0 = t[imax[0]] - 0.5 * T_beat if len(imax) else t[0]; tz1 = tz0 + 3.2 * T_beat
    ax.plot(t, x, lw=0.6, color="0.6", label=labels[0])
    ax.plot(tb, xb, lw=0.7, color="tab:green", label=labels[1])
    ax.set_xlim(tz0, tz1)
    envb_s = np.convolve(envb, np.ones(k) / k, mode="same"); envb_s[:edge] = np.nan; envb_s[-edge:] = np.nan
    ax.plot(tb, envb_s, "--", color="tab:green", lw=1, label=f"envelope of {labels[1]}")
    ax.plot(t, env_s, "--", color="0.5", lw=1, label=f"envelope of {labels[0]}")
    ax.set_ylabel("displacement / mm"); ax.set_xlabel("$t$ / s"); ax.grid(alpha=0.3)
    ax.set_ylim(-1.55 * np.nanmax(env_s), 1.55 * np.nanmax(env_s))          # head-room for the legend
    ax.legend(loc="upper right", fontsize=8, ncol=4, framealpha=0.9)
    ax.set_title(titles[1], fontsize=10)

ax = axs[-1] if show_fft else plt.figure().add_subplot()
ax.set_title(titles[2], fontsize=10)
ax.plot(f, X, color="tab:red", lw=1, label=f"|FFT| of {labels[0]}")
if len(data) > 1:
    fb_, Xb_ = spectrum(data[1][1], data[1][2], window=args.window); ax.plot(fb_, Xb_, color="tab:green", lw=1, alpha=0.8, label=f"|FFT| of {labels[1]}")
ax.axvline(f_avg, color="tab:blue", ls="--", lw=1)
yA = 1.10
for fq in (f1, f2): ax.plot([fq, fq], [np.interp(fq, f, X), yA], ":", color="tab:red", lw=0.8)
ax.annotate("", xy=(f2, yA), xytext=(f1, yA), arrowprops=dict(arrowstyle="<->", color="tab:red", lw=1.3))
bb = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)
ax.text(f_avg, yA + 0.03, rf"$f_{{\rm beat}}=f_2-f_1={f2-f1:.3f}$ Hz $=1/T_{{\rm beat}}$", ha="center", va="bottom", color="tab:red", fontsize=9, bbox=bb)
ax.text(f_avg + 0.008, 0.80, rf"$f_{{\rm avg}}=(f_1+f_2)/2$" + "\n" + rf"$={f_avg:.3f}$ Hz", ha="left", va="center", color="tab:blue", fontsize=9, bbox=bb)
ax.text(f1 - 0.02, np.interp(f1, f, X) * 0.9, f"$f_1$ = {f1:.3f} Hz", ha="right", va="center", fontsize=9)
ax.text(f2 + 0.02, np.interp(f2, f, X) * 0.9, f"$f_2$ = {f2:.3f} Hz", ha="left", va="center", fontsize=9)
ax.set_xlim(max(0, f1 - 0.35), f2 + 0.35); ax.set_ylim(0, 1.32)
ax.set_xlabel("frequency / Hz"); ax.set_ylabel("amplitude spectrum\n(normalised)"); ax.legend(loc="upper right", fontsize=8); ax.grid(alpha=0.3)
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)   # create the output folder if needed
fig.tight_layout(); fig.savefig(args.out, bbox_inches="tight"); print("saved", args.out)

# ---- one-line summary for the results table (append mode; header written once)
if args.csv:
    import csv
    r2_val = r2 if len(data) > 1 else np.nan
    row = dict(d_mm=args.d, trial=args.trial, T_rec_s=round(T_rec, 2), fs_Hz=round(fs, 2), bin_mHz=round(1e3 / T_rec, 1), window=args.window,
               f1_Hz=round(f1, 4), f2_Hz=round(f2, 4), f_avg_Hz=round(f_avg, 4), f_beat_Hz=round(f2 - f1, 4),
               f2sq_minus_f1sq_Hz2=round(f2 ** 2 - f1 ** 2, 4), T_beat_fit_s=round(T_beat_meas, 4), dT_beat_fit_s=round(dT_beat_meas, 4),
               n_beats=Nb, r1=round(r_pk, 3), r2=round(r2_val, 3), env_min_max=round(ratio, 3), tau_s=round(tau, 1), file=args.files[0])
    new_file = not os.path.exists(args.csv)
    with open(args.csv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row)); (w.writeheader() if new_file else None); w.writerow(row)
    print("appended to", args.csv)
