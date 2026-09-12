"""Systematic frequency error from the non-linear (r^-4) force: pure anti-phase mode started with
x1 = -a, x2 = +a (separation swings by 2a), period measured from zero crossings over many cycles.
In a release experiment the anti-phase mode has amplitude a = xi0/2 ~ A/2."""
import numpy as np
from scipy.integrate import solve_ivp
f1, C = 6.49, 1.33e-7      # real time base (240 fps); percentages are independent of the time base
w1 = 2*np.pi*f1; k_m = w1**2
def gamma_m(d):  return 2*np.pi**2*C/d**5
def F_m(r, d):   return gamma_m(d)*d/4*(d/r)**4
def f2_measured(d, a, ncyc=40):
    g = gamma_m(d); f2lin = np.sqrt(f1**2 + C/d**5)
    def rhs(t, y):
        x1, v1, x2, v2 = y; r = d + x2 - x1; dF = F_m(r, d) - F_m(d, d)
        return [v1, -k_m*x1 - dF, v2, -k_m*x2 + dF]
    T = ncyc/f2lin
    sol = solve_ivp(rhs, (0, T), [-a, 0, a, 0], rtol=1e-11, atol=1e-13, dense_output=True)
    t = np.linspace(0, T, 200000); x1 = sol.sol(t)[0]
    zc = t[:-1][(x1[:-1] < 0) & (x1[1:] >= 0)]           # upward zero crossings
    return f2lin, 1/np.mean(np.diff(zc))
print(" d/mm   f2lin   a=2.5mm(A=5)  a=5mm(A=10)   shift%  shift%")
for d_mm in [24.5, 26.5, 28.5, 31, 33.5, 36, 39, 42]:
    d = d_mm/1000
    f2l, fa = f2_measured(d, 2.5e-3); _, fb = f2_measured(d, 5e-3)
    print(f"{d_mm:5}  {f2l:6.3f}   {fa:8.4f}     {fb:8.4f}     {100*(fa/f2l-1):6.2f}  {100*(fb/f2l-1):6.2f}")
