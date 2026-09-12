"""Numerical verification of the non-linear (extension) results used in theory_v2.tex.
Real time base: f1 = 6.49 Hz, C = 1.33e-7 m^5 s^-2  (run 1 x8), tau_d = 21 s.
Relative coordinate xi = x2 - x1 obeys  m xi'' = -k xi + 2[F(d+xi) - F(d)],  F = F_d (d/r)^4.
Perturbation theory (Lindstedt-Poincare, to O(a^2)) with g = gamma/(k+2gamma):
  <xi> = 5 g a^2/(2d);  A2/a = 5 g a/(6d);  A3/a = (a/d)^2 (5g/16 + 25 g^2/48);
  dw/w2 = (a/d)^2 (15 g/4 - 125 g^2/12)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
f1, C = 6.49, 1.33e-7
w1 = 2*np.pi*f1
def params(d):
    gam_m = 2*np.pi**2*C/d**5            # gamma/m
    w2 = np.sqrt(w1**2 + 2*gam_m); g = gam_m/w2**2
    return gam_m, w2, g
def xi_solution(d, xi0, T, tau=None):
    gam_m, w2, g = params(d)
    Fd_m = gam_m*d/4                         # F(d)/m
    lam = 0 if tau is None else 1/tau
    def rhs(t, y):
        xi, v = y
        return [v, -w1**2*xi + 2*Fd_m*((1+xi/d)**-4 - 1) - 2*lam*v]
    t = np.linspace(0, T, int(T*2000))
    sol = solve_ivp(rhs, (0, T), [xi0, 0], t_eval=t, rtol=1e-11, atol=1e-14)
    return t, sol.y[0]
print("d/mm  xi0/mm  a_num  f2lin   f_num   shift_num%  shift_pert%(a)  <xi>num  <xi>pert  A2/a num  pert   A3/a num  pert")
for d_mm in [24.5, 31.0]:
    d = d_mm/1000; gam_m, w2, g = params(d)
    for xi0_mm in [2, 4.24, 6, 8.5]:
        xi0 = xi0_mm/1000; T = 40/(w2/2/np.pi)     # 40 cycles
        t, xi = xi_solution(d, xi0, T)
        pk, _ = find_peaks(xi); f_num = 1/np.mean(np.diff(t[pk]))
        w = 2*np.pi*f_num
        M = np.column_stack([np.ones_like(t), np.cos(w*t), np.cos(2*w*t), np.cos(3*w*t), np.sin(w*t), np.sin(2*w*t), np.sin(3*w*t)])
        c, *_ = np.linalg.lstsq(M, xi, rcond=None)
        a = np.hypot(c[1], c[4]); A2 = np.hypot(c[2], c[5]); A3 = np.hypot(c[3], c[6]); dc = c[0]
        f2lin = w2/2/np.pi
        sh_p = (a/d)**2*(15*g/4 - 125*g**2/12)
        print(f"{d_mm:5.1f} {xi0_mm:6.2f} {a*1e3:6.3f} {f2lin:6.3f} {f_num:7.4f} {100*(f_num/f2lin-1):9.2f} {100*sh_p:12.2f}  {dc*1e3:7.3f} {5*g*a**2/(2*d)*1e3:8.3f}  {A2/a:8.4f} {5*g*a/(6*d):7.4f}  {A3/a:8.5f} {(a/d)**2*(5*g/16+25*g**2/48):8.5f}")
print("\ng at 24.5 / 31 / 42 mm:", [round(params(d/1000)[2],4) for d in (24.5, 31, 42)])
print("2gamma/k = 3 (w2 = 2 w1) at d* =", (C/(3*f1**2))**0.2*1e3, "mm")
print("\nFFT-averaged shift with tau = 21 s over 45 s (Hann):")
for d_mm, A_mm in [(24.5, 5), (24.5, 10), (31, 5), (31, 10)]:
    d = d_mm/1000; gam_m, w2, g = params(d); gk = gam_m/w1**2
    xi0 = A_mm/1000/(1+gk)
    t, xi = xi_solution(d, xi0, 45, tau=21)
    fs = 1/(t[1]-t[0]); x = (xi-xi.mean())*np.hanning(len(xi)); n = 2**22
    X = np.abs(np.fft.rfft(x, n)); f = np.fft.rfftfreq(n, 1/fs); m = (f > 5) & (f < 9); i = np.argmax(X[m]); fpk = f[m][i]
    f2lin = w2/2/np.pi; sh_init = (xi0/d)**2*(15*g/4 - 125*g**2/12)
    print(f"d={d_mm} A={A_mm}: xi0={xi0*1e3:.2f} mm, initial-amplitude shift {100*sh_init:.2f} %, FFT peak shift {100*(fpk/f2lin-1):.2f} %  ratio {(fpk/f2lin-1)/sh_init:.2f}")
