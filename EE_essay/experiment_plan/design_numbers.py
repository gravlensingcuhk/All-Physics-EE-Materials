"""Design numbers for the redo: expected frequencies, beat periods, record lengths and the
amplitude (non-linearity) effect, using the constants from the first data set
(f1 = 0.811 Hz, C = f2^2-f1^2 * d^5 = 2.08e-9 m^5 s^-2).  Full r^-4 force, hold-and-release."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

f1, C = 0.811, 2.08e-9
w1 = 2*np.pi*f1; k_m = w1**2                        # k/m

def gamma_m(d):  return 2*np.pi**2*C/d**5           # gamma/m  (from f2^2-f1^2 = gamma/(2 pi^2 m))
def F_m(r, d):   return gamma_m(d)*d/4*(d/r)**4      # F(r)/m   (F = gamma d/4 at r = d)

def simulate(d, A, T, fs=60.0):
    g = gamma_m(d)
    # spring 1 held pulled OUTWARD by A (x1 = -A, separation grows); spring 2 relaxes to x20
    x20 = brentq(lambda x2: k_m*x2 - (F_m(d + A + x2, d) - F_m(d, d)), -A, A)
    def rhs(t, y):
        x1, v1, x2, v2 = y; r = d + x2 - x1; dF = F_m(r, d) - F_m(d, d)
        return [v1, -k_m*x1 - dF, v2, -k_m*x2 + dF]
    t = np.arange(0, T, 1/fs)
    sol = solve_ivp(rhs, (0, T), [-A, 0, x20, 0], t_eval=t, rtol=1e-10, atol=1e-12)
    return t, sol.y[0], sol.y[2], x20

def two_tone(t, x, f1g, f2g):
    def model(p): return p[0]*np.cos(2*np.pi*p[1]*t+p[2]) + p[3]*np.cos(2*np.pi*p[4]*t+p[5]) + p[6]
    res = least_squares(lambda p: model(p)-x, [x.std(), f1g, 0, x.std(), f2g, 0, 0])
    return sorted([abs(res.x[1]), abs(res.x[4])])

print(f"{'d/mm':>5} {'f2lin':>6} {'fbeat':>6} {'Tbeat':>6} {'Trec':>5} | A=5mm: {'x20':>5} {'f2 shift%':>9} | A=10mm: {'x20':>5} {'f2 shift%':>9} | {'delta':>5} {'ratio':>5}")
rows = []
for d_mm in [20, 22, 24.5, 27.5, 30.5, 34, 38, 42]:
    d = d_mm/1000; Y = C/d**5; f2 = np.sqrt(f1**2+Y); fb = f2-f1; Tb = 1/fb
    Trec = max(30, 10*np.ceil(3*Tb/10))
    delta = np.pi**2*Y*d/(2*w1**2)*1000
    out = []
    for A_mm in (5, 10):
        A = A_mm/1000; T = max(60, 4*Tb)
        t, x1, x2, x20 = simulate(d, A, T)
        fa, fbv = two_tone(t, x1, f1, f2)
        out += [x20*1000, (fbv/f2-1)*100]
    rows.append((d_mm, f2, fb, Tb, Trec, *out, delta, (f1/f2)**2))
    print(f"{d_mm:5} {f2:6.3f} {fb:6.3f} {Tb:6.1f} {Trec:5.0f} | {out[0]:12.2f} {out[1]:9.2f} | {out[2]:13.2f} {out[3]:9.2f} | {delta:5.2f} {(f1/f2)**2:5.2f}")
np.save("design_rows.npy", np.array(rows))
