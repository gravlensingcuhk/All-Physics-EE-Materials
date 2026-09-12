"""Which spring length?  -- REAL time base.
Run 1 was analysed with Tracker reading the 240 fps slow-motion file at its 30 fps playback rate, so every
frequency was 8x too low (0.811 Hz -> 6.49 Hz; the user's own FFT with the clip set to 240 fps shows 6.51 Hz).
Ratios (f2/f1, 2 gamma/k, delta, x20, non-linear shifts, the exponent) are unchanged; C scales by 64, tau by 1/8.

The design depends only on the single-spring frequency f_a (same magnets, C fixed).  Gravity is now a
few-per-cent correction:  k_eff = k - (1.2 M + 0.375 mu L) g / L,  i.e. (1.2..1.6) g/(L w^2) ~ 3-6 % at 6.5 Hz.
"""
import numpy as np
C, f_run1, tau = 1.33e-7, 6.49, 21.0
ds = np.array([24.5, 26.5, 28.5, 31, 33.5, 36, 39, 42])/1000
Trec = 40.0

print("f_a   L/L0  df(24.5) df(42)  Tb(42)  3Tb(42)  dmax(40 s)  2g/k(24.5)  delta(24.5)  x20(A=5)")
rows = []
for fa in [4.5, 5.0, 5.5, 6.5, 8.0, 10.0]:
    f2 = np.sqrt(fa**2 + C/ds**5); df = f2 - fa
    dfmin = 3/Trec                                   # 3 T_beat must fit into the record
    dmax = (C/(dfmin*(2*fa + dfmin)))**0.2*1e3
    gk = C/ds[0]**5/fa**2; delta = gk/2*ds[0]/4*1e3; x20 = 5*(gk/2)/(1 + gk/2)
    LL = (f_run1/fa)**(2/3)
    print(f"{fa:4.1f}  {LL:4.2f}  {df[0]:6.2f}   {df[-1]:6.3f}  {1/df[-1]:5.1f}   {3/df[-1]:5.0f}    {dmax:5.1f}       {gk:5.2f}       {delta:4.1f}        {x20:4.2f}")
    rows.append((fa, LL, df[0], df[-1]*1e3, 1/df[-1], dmax, gk, delta))

print("\nLaTeX rows:")
for fa, LL, d0, d1, tb, dmax, gk, delta in rows:
    lab = f"{fa:.1f}" + (" (run~1)" if fa == 6.5 else "")
    print(f"{lab} & {LL:.2f} & {d0:.2f} Hz & {d1:.0f} mHz & {tb:.1f} s & {dmax:.0f} mm & {gk:.2f} & {delta:.1f} mm\\\\")

print("\nladder for f_a = 6.49 Hz:")
f1 = f_run1; f2 = np.sqrt(f1**2 + C/ds**5); df = f2 - f1; Tb = 1/df; gk = C/ds**5/f1**2
print("d      ", " & ".join(f"{d*1e3:g}" for d in ds))
print("f2     ", " & ".join(f"{v:.2f}" for v in f2))
print("f2-f1  ", " & ".join(f"{v:.3f}" for v in df))
print("Tbeat  ", " & ".join(f"{v:.2f}" if v < 10 else f"{v:.1f}" for v in Tb))
print("3Tbeat ", " & ".join(f"{v:.1f}" for v in 3*Tb))
print("beats in 40 s", " & ".join(f"{40*v:.0f}" for v in df))
print("delta  ", " & ".join(f"{v:.2f}" for v in gk/2*ds/4*1e3))
print("x20    ", " & ".join(f"{v:.2f}" for v in 5*(gk/2)/(1 + gk/2)))
print("half-bin error on f2^2-f1^2 for T = 40 s (%):", " & ".join(f"{v:.1f}" for v in 2*(0.5/40)*np.sqrt(f1**2 + f2**2)/(f2**2 - f1**2)*100))
print("same for T = 6 s (%):", " & ".join(f"{v:.0f}" for v in 2*(0.5/6)*np.sqrt(f1**2 + f2**2)/(f2**2 - f1**2)*100))
w2 = (2*np.pi*f1)**2
print(f"\ngravity fraction of k, L = 0.15 / 0.25 m: {1.2*9.81/0.15/w2*100:.1f}-{1.6*9.81/0.15/w2*100:.1f} % / {1.2*9.81/0.25/w2*100:.1f}-{1.6*9.81/0.25/w2*100:.1f} %")
print(f"tip lean per degree of base tilt: {9.81*np.radians(1)/w2*1e3:.2f} mm;  k for m = 3, 5, 10 g: {[round(m*w2, 1) for m in (0.003, 0.005, 0.01)]} N/m")
print(f"amplitude left after 40 s (tau = 21 s): {np.exp(-40/tau)*100:.0f} %;  free decay 5 mm -> 0.1 mm takes {tau*np.log(50):.0f} s")
print(f"spring 2 ringing left after a 1 s pull at 24.5 mm: {0.76/(2*np.pi*f1*1.0):.3f} mm (vs 0.76 mm for a jerk)")
