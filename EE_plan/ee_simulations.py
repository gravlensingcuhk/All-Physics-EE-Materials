"""
Simulations supporting the EE plan for the coupled magnetic-mechanical oscillator.

Part A  - FFT error analysis (frequency resolution, leakage, peak pulling) using
          synthetic two-tone signals built from YOUR measured f1, f2.
Part B  - Non-linear (un-linearised) magnetic force: RK4 integration of the full
          equations of motion + FFT, to quantify the amplitude-dependent frequency
          shift and its effect on the ln(f2^2 - f1^2) vs ln(d) gradient.
Part C  - Effect of a distance offset (surface gap vs centre-to-centre, or static
          outward bending) on the apparent gradient.

All parameters that I had to ASSUME are marked with  # ASSUMED  - replace with
your measured values.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).parent / "figs"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------------- your data
d_mm = np.array([21, 23, 25, 28, 31, 35, 40.0])
f1 = np.array([0.823, 0.808, 0.812, 0.810, 0.807, 0.806, 0.802])
f2 = np.array([1.063, 1.005, 0.953, 0.907, 0.881, 0.845, 0.811])


# ---------------------------------------------------------------- helpers
def spectrum(x, fs, window="hann", nfft=2 ** 16):
    n = len(x)
    w = np.hanning(n) if window == "hann" else np.ones(n)
    X = np.abs(np.fft.rfft((x - x.mean()) * w, nfft))
    f = np.fft.rfftfreq(nfft, 1 / fs)
    return f, X / X.max()


def two_peaks(f, X, fmin=0.3, fmax=3.0, rel_thresh=0.05):
    """Return the two highest *separate* local maxima with parabolic interpolation."""
    m = (f > fmin) & (f < fmax)
    idx = np.where(m)[0]
    peaks = []
    for i in idx[1:-1]:
        if X[i] > X[i - 1] and X[i] >= X[i + 1] and X[i] > rel_thresh:
            # parabolic (quadratic) interpolation of the peak position
            a, b, c = X[i - 1], X[i], X[i + 1]
            p = 0.5 * (a - c) / (a - 2 * b + c)
            peaks.append((f[i] + p * (f[1] - f[0]), b))
    peaks.sort(key=lambda t: -t[1])
    return sorted(pk[0] for pk in peaks[:2])


# ================================================================= PART A
def part_A():
    fs = 40.0          # ASSUMED sample rate (Hz)  (Edwin used 40 Hz)
    tau = 30.0         # ASSUMED amplitude decay time of the oscillation (s)
    print("\n=== PART A : FFT resolution / two-peak resolvability ===")
    print(f"{'d/mm':>5} {'f2-f1':>7} {'bins@20s':>9} {'T_min rect':>11} {'T_min hann':>11}")
    rows = []
    for d, a, b in zip(d_mm, f1, f2):
        df = b - a
        # minimum record length needed: separation >= 1 bin (rect) or >= 2 bins (hann)
        rows.append((d, df, df / (1 / 20), 1 / df, 2 / df))
        print(f"{d:5.0f} {df:7.3f} {df/(1/20):9.2f} {1/df:11.1f} {2/df:11.1f}")

    # --- spectra for three separations at 20 s and 60 s
    fig, axs = plt.subplots(2, 3, figsize=(13, 6.5), sharex=True)
    for col, d_sel in enumerate([25, 35, 40]):
        i = list(d_mm).index(d_sel)
        for row, T in enumerate([20, 60]):
            t = np.arange(0, T, 1 / fs)
            x = np.exp(-t / tau) * (np.cos(2 * np.pi * f1[i] * t) + np.cos(2 * np.pi * f2[i] * t))
            ax = axs[row, col]
            for win, c in [("rect", "tab:red"), ("hann", "tab:blue")]:
                f, X = spectrum(x, fs, win)
                ax.plot(f, X, c, lw=1.2, label=f"{win}")
                pk = two_peaks(f, X)
                ax.plot(pk, np.interp(pk, f, X), "o", color=c, ms=5)
            ax.axvline(f1[i], color="k", ls=":", lw=0.8)
            ax.axvline(f2[i], color="k", ls=":", lw=0.8)
            ax.set_xlim(0.6, 1.3)
            ax.set_title(f"d = {d_sel} mm, T = {T} s  (bin = {1/T:.3f} Hz, true Δf = {f2[i]-f1[i]:.3f} Hz)", fontsize=9)
            if row == 1:
                ax.set_xlabel("f / Hz")
            if col == 0:
                ax.set_ylabel("normalised |FFT|")
            ax.legend(fontsize=8)
    fig.suptitle("Synthetic two-tone signal built from your f1, f2 (dotted = true values). "
                 "Peaks merge when Δf is below ~1–2 bins.", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "A_fft_resolution.png", dpi=130)

    # --- peak-pulling bias: recovered vs true, for T = 20 s, both windows
    print("\nRecovered peak frequencies from a T = 20 s record (fs = 40 Hz, tau = 30 s):")
    print(f"{'d':>4} {'win':>5} {'f1_rec':>7} {'f2_rec':>7} {'err1/mHz':>9} {'err2/mHz':>9}")
    T = 20
    t = np.arange(0, T, 1 / fs)
    for i, d in enumerate(d_mm):
        x = np.exp(-t / tau) * (np.cos(2 * np.pi * f1[i] * t) + np.cos(2 * np.pi * f2[i] * t))
        for win in ["rect", "hann"]:
            f, X = spectrum(x, fs, win)
            pk = two_peaks(f, X)
            if len(pk) == 2:
                print(f"{d:4.0f} {win:>5} {pk[0]:7.3f} {pk[1]:7.3f} {1e3*(pk[0]-f1[i]):9.1f} {1e3*(pk[1]-f2[i]):9.1f}")
            else:
                print(f"{d:4.0f} {win:>5}   UNRESOLVED (single peak at {pk[0]:.3f} Hz)")


# ================================================================= PART B
def rk4(fun, y0, t):
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    h = t[1] - t[0]
    for i in range(len(t) - 1):
        k1 = fun(y[i])
        k2 = fun(y[i] + 0.5 * h * k1)
        k3 = fun(y[i] + 0.5 * h * k2)
        k4 = fun(y[i] + h * k3)
        y[i + 1] = y[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return y


def part_B():
    print("\n=== PART B : non-linear magnetic force (RK4 + FFT) ===")
    w1sq = (2 * np.pi * 0.81) ** 2          # k/m from your f1 (≈0.81 Hz)
    # calibrate C/m (F = C r^-4) so that the LINEAR model reproduces f2 at d = 25 mm
    d_cal, f2_cal, f1_cal = 0.025, 0.953, 0.812
    wm_sq_cal = (2 * np.pi) ** 2 * (f2_cal ** 2 - f1_cal ** 2)   # = 2*gamma/m = 8C/(m d^5)
    C_over_m = wm_sq_cal * d_cal ** 5 / 8
    print(f"calibrated C/m = {C_over_m:.3e} m^5 s^-2  (F/m = (C/m) r^-4)")

    def make_rhs(d):
        def rhs(y):
            x1, v1, x2, v2 = y
            r = d + x2 - x1
            dF = C_over_m * (r ** -4 - d ** -4)          # (F(r) - F(d))/m
            return np.array([v1, -w1sq * x1 - dF, v2, -w1sq * x2 + dF])
        return rhs

    def measure(d, A, T=150.0, fs_int=1000.0, fs_fft=50.0):
        t = np.arange(0, T, 1 / fs_int)
        y = rk4(make_rhs(d), np.array([A, 0, 0, 0]), t)
        step = int(fs_int / fs_fft)
        x1 = y[::step, 0]
        f, X = spectrum(x1, fs_fft, "hann", nfft=2 ** 18)
        pk = two_peaks(f, X, fmin=0.3, fmax=2.5, rel_thresh=0.02)
        return pk

    # (i) amplitude dependence of f2 at d = 21 mm and 35 mm
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    amps = np.array([1, 2, 3, 5, 7, 10]) * 1e-3
    for d_sel, c in [(0.021, "tab:red"), (0.035, "tab:blue")]:
        f2_lin = np.sqrt(w1sq + 8 * C_over_m / d_sel ** 5) / (2 * np.pi)
        f2_nl = []
        for A in amps:
            pk = measure(d_sel, A)
            f2_nl.append(pk[-1])
        f2_nl = np.array(f2_nl)
        ax.plot(amps * 1e3, f2_nl, "o-", color=c, label=f"d = {d_sel*1e3:.0f} mm (RK4)")
        ax.axhline(f2_lin, color=c, ls="--", lw=1, label=f"linearised prediction, d = {d_sel*1e3:.0f} mm")
        print(f"d = {d_sel*1e3:.0f} mm : f2_lin = {f2_lin:.4f} Hz ; f2 at A = 10 mm = {f2_nl[-1]:.4f} Hz "
              f"({100*(f2_nl[-1]/f2_lin-1):+.1f} %) ; at A = 2 mm = {f2_nl[1]:.4f} Hz ({100*(f2_nl[1]/f2_lin-1):+.2f} %)")
    ax.set_xlabel("initial displacement A of spring 1 / mm")
    ax.set_ylabel("f2 (antiphase mode) / Hz")
    ax.set_title("Amplitude dependence of f2 caused by the r^-4 magnetic force", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "B_f2_vs_amplitude.png", dpi=130)

    # (ii) ln-ln gradient: linear theory vs non-linear simulation for A = 10 mm and 2 mm
    print("\nln(f2^2-f1^2) vs ln(d) gradients:")
    dd = d_mm[:-1] * 1e-3            # 21..35 mm (the 40 mm point is FFT-limited)
    y_lin = (8 * C_over_m / dd ** 5) / (2 * np.pi) ** 2
    print(f"  linearised theory            : {np.polyfit(np.log(dd), np.log(y_lin), 1)[0]:.3f}")
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(np.log(dd * 1e3), np.log(y_lin), "k--", label="linearised theory (gradient −5)")
    for A, c in [(2e-3, "tab:green"), (10e-3, "tab:red")]:
        ys = []
        for d in dd:
            pk = measure(d, A)
            ys.append(pk[-1] ** 2 - pk[0] ** 2)
        ys = np.array(ys)
        g, b = np.polyfit(np.log(dd), np.log(ys), 1)
        print(f"  RK4 non-linear, A = {A*1e3:.0f} mm  : {g:.3f}")
        ax.plot(np.log(dd * 1e3), np.log(ys), "o-", color=c, label=f"non-linear simulation, A = {A*1e3:.0f} mm (gradient {g:.2f})")
    ymeas = f2[:-1] ** 2 - f1[:-1] ** 2
    ax.plot(np.log(d_mm[:-1]), np.log(ymeas), "s", color="tab:purple", label="your data (gradient −3.75)")
    ax.set_xlabel("ln(d / mm)")
    ax.set_ylabel("ln[(f2² − f1²) / Hz²]")
    ax.legend(fontsize=8)
    ax.set_title("Effect of the un-linearised force on the ln–ln gradient", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "B_lnln_gradient.png", dpi=130)


# ================================================================= PART C
def part_C():
    print("\n=== PART C : apparent gradient if the true law is (d+δ)^-5 but ln d is plotted ===")
    dd = d_mm[:-1]
    for delta in [0, 2, 4, 6, 8, 10, 12]:
        g = np.polyfit(np.log(dd), -5 * np.log(dd + delta), 1)[0]
        print(f"  δ = {delta:2d} mm -> apparent gradient {g:.2f}")
    print("  (δ = magnet thickness if you measured the surface gap, and/or 2x static outward bending)")


if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
    print("\nfigures written to", OUT)
