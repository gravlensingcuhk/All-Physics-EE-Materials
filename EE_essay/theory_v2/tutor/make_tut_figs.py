"""Figures for the theory tutorial (tutor.typ).
House style copied from EE_essay/theory_v2/make_figs_v2.py:
  copper #b35900 (spring 1 / leaf springs), crimson (spring 2 / f2), gray base.
Outputs PDF vectors into figs/.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, Circle
import os

plt.rcParams.update({
    "font.size": 10, "axes.linewidth": 0.9, "axes.grid": True,
    "grid.alpha": 0.3, "mathtext.fontset": "cm", "axes.prop_cycle":
    plt.cycler(color=["#b35900", "crimson", "0.4"]),
})
C1, C2, CG, CB = "#b35900", "crimson", "0.4", "tab:blue"
MG_S, MG_N = (1.0, 0.30, 0.30), (0.30, 0.30, 1.0)   # red!70, blue!70 in rgb
FIG = "figs"
os.makedirs(FIG, exist_ok=True)

f1, C, tau = 6.49, 1.33e-7, 21.0
w1 = 2*np.pi*f1
def modes(d):
    gam_m = 2*np.pi**2*C/d**5
    return gam_m, np.sqrt(w1**2 + 2*gam_m)

def phi(u):            # shape function with s = u*L, phi(L)=1
    return (3*u**2 - u**3)/2
def draw_spring(ax, xc, Lpx, tip_dx, color=C1, dashed=False):
    """vertical cantilever drawn from (xc,0) to (xc,Lpx), tip displaced by tip_dx (px units)"""
    u = np.linspace(0, 1, 60)
    ax.plot(xc + tip_dx*phi(u), Lpx*u, color=color, lw=2.4 if not dashed else 1.1,
            ls="--" if dashed else "-", alpha=1 if not dashed else 0.55, zorder=3)
def draw_magnet(ax, x, y, w=0.55, h=0.42, facing="N"):
    """magnet of height h centred at (x,y); facing pole on its right edge"""
    ax.add_patch(Rectangle((x-w/2, y-h/2), w/2, h, fc=MG_S, ec="k", lw=0.6, zorder=5))
    ax.add_patch(Rectangle((x,     y-h/2), w/2, h, fc=MG_N, ec="k", lw=0.6, zorder=5))
    if facing == "N":   # N to the right, S to the left
        ax.text(x+w/4, y, "N", color="w", fontsize=7, ha="center", va="center", zorder=6)
        ax.text(x-w/4, y, "S", color="w", fontsize=7, ha="center", va="center", zorder=6)
    else:
        ax.add_patch(Rectangle((x-w/2, y-h/2), w/2, h, fc=MG_N, ec="k", lw=0.6, zorder=5))
        ax.add_patch(Rectangle((x,     y-h/2), w/2, h, fc=MG_S, ec="k", lw=0.6, zorder=5))
        ax.text(x-w/4, y, "N", color="w", fontsize=7, ha="center", va="center", zorder=6)
        ax.text(x+w/4, y, "S", color="w", fontsize=7, ha="center", va="center", zorder=6)

# ------------------------------------------------------------------ 1. apparatus
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.set_xlim(-4.1, 4.1); ax.set_ylim(-0.85, 5.5); ax.axis("off"); ax.grid(False)
ax.add_patch(Rectangle((-3.8, -0.5), 7.6, 0.5, fc="0.72", ec="k", lw=1))
ax.text(0, -0.25, "non-magnetic base", ha="center", va="center", fontsize=9)
y0, H = 0.5, 3.2
tip_y = y0 + H                       # height of the tip / magnet centre
mc = {}
for sx, col, sgn in [(-1.35, C1, -1), (1.35, C2, +1)]:
    ax.add_patch(Rectangle((sx-0.18, 0), 0.36, y0, fc="0.62", ec="k", lw=0.8))       # clamp
    ax.plot([sx, sx], [y0, tip_y+0.45], ls="--", color="0.55", lw=0.9)               # rest line
    u = np.linspace(0, 1, 60)
    ax.plot(sx + sgn*0.42*phi(u), y0 + H*u, color=col, lw=2.4, zorder=3)             # bent spring
    mc[sgn] = sx + sgn*0.42 + sgn*0.62
    draw_magnet(ax, mc[sgn], tip_y, facing="N" if sgn < 0 else "S")
ax.annotate("", (-0.75, 5.15), (-0.15, 5.15), arrowprops=dict(arrowstyle="->", color="k"))
ax.text(-0.85, 5.15, "$x_1>0$", ha="right", va="center", fontsize=10)
ax.annotate("", (1.95, 5.15), (1.35, 5.15), arrowprops=dict(arrowstyle="->", color="k"))
ax.text(2.05, 5.15, "$x_2>0$ (both positive to the right)", ha="left", va="center", fontsize=9.5)
ax.plot([mc[-1], mc[-1]], [tip_y+0.21, 4.55], ls=":", color="0.4", lw=0.8)
ax.plot([mc[1], mc[1]], [tip_y+0.21, 4.55], ls=":", color="0.4", lw=0.8)
ax.annotate("", (mc[-1], 4.55), (mc[1], 4.55), arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
ax.text(0, 4.7, "$d$ = centre-to-centre separation at equilibrium", ha="center", fontsize=9.5)
ax.text(0, tip_y, "N–N\nfaces:\nrepulsion", ha="center", va="center", fontsize=8, color="0.3")
ax.annotate("", (-2.35, y0), (-2.35, y0+H), arrowprops=dict(arrowstyle="<->", color="k", lw=1))
ax.text(-2.5, y0+H/2, "$L$ free length", rotation=90, va="center", ha="center", fontsize=9.5)
ax.text(3.3, 2.0, "even at rest both\nsprings bow outwards\nby $\\delta=F(d)/k$", fontsize=8.5, ha="center", color="0.3")
ax.set_title("The apparatus: two leaf springs coupled by magnetic repulsion", fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_setup.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 2. shape function + integrals
fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.1))
u = np.linspace(0, 1, 200)
axs[0].plot(u, phi(u), color=C1, lw=2)
axs[0].plot(u, u**1.5, color="0.6", lw=1, ls=":")
axs[0].annotate("$\\varphi(0)=\\varphi'(0)=0$\n(clamped: no motion, no slope)", (0.02, 0.03), (0.34, 0.28),
                arrowprops=dict(arrowstyle="->", color="0.3"), fontsize=8.5)
axs[0].annotate("$\\varphi''(1)=0$ (free tip carries\nno bending moment); $\\varphi(1)=1$", (0.99, 1.0), (0.42, 0.75),
                arrowprops=dict(arrowstyle="->", color="0.3"), fontsize=8.5)
axs[0].set_xlabel("$s/L$"); axs[0].set_ylabel("$\\varphi(s)$"); axs[0].set_ylim(-0.03, 1.35)
axs[0].set_title("deflection shape $y(s,t)=x(t)\\,\\varphi(s)$", fontsize=10, loc="left")
axs[1].plot(u, phi(u)**2, color="0.25", lw=2)
axs[1].fill_between(u, 0, phi(u)**2, alpha=0.25, color=C1)
axs[1].text(0.55, 0.42, "$\\int_0^L \\varphi^2\\,\\mathrm{d}s=\\frac{33}{140}L$\n\nmass $\\mu\\,\\mathrm{d}s$ moves at speed\n$\\dot x\\varphi(s)$: only the shaded\nshare counts as oscillating mass", fontsize=8.5, ha="left")
axs[1].set_xlabel("$s/L$"); axs[1].set_ylim(0, 1.45)
axs[1].set_title(r"where the moving mass is ($\propto\varphi^2$)", fontsize=10, loc="left")
k2 = (6 - 6*u)
axs[2].plot(u, k2, color=C2, lw=2)
axs[2].fill_between(u, 0, k2, alpha=0.2, color=C2)
axs[2].text(0.45, 3.4, "$\\int_0^L(\\varphi'')^2\\,\\mathrm{d}s=\\frac{3}{L^3}$\n\n$\\varphi''=3(L-s)/L^3$: curvature is largest\nat the clamp — that is where the\nbending energy lives", fontsize=8.5, ha="left")
axs[2].set_xlabel("$s/L$"); axs[2].set_title(r"where the bending is ($\propto\varphi''^{\,2}$)", fontsize=10, loc="left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_shape.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 3. dipole: field + power-law ladder
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.3))
# (a) sketch of the two anti-parallel moments
axs[0].set_xlim(-3, 3); axs[0].set_ylim(-1.5, 1.7); axs[0].axis("off"); axs[0].grid(False)
draw_magnet(axs[0], -1.5, -0.9, w=1.0, h=0.55, facing="N")
draw_magnet(axs[0],  1.5, -0.9, w=1.0, h=0.55, facing="S")
axs[0].annotate("", (-1.05, -0.35), (-1.95, -0.35), arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
axs[0].text(-1.5, -0.2, "$\\vec{\\mu}_1$", ha="center", fontsize=10)
axs[0].annotate("", (1.05, -0.35), (1.95, -0.35), arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
axs[0].text(1.5, -0.2, "$\\vec{\\mu}_2$", ha="center", fontsize=10)
axs[0].text(0, -0.55, "moments point towards each other\n(N faces N)", ha="center", fontsize=8.5, color="0.3")
# field of magnet 1 at the position of magnet 2
r = np.linspace(0.7, 2.6, 100)
axs[0].plot(r, 0.30 + 1.0/r**3, color=C2, lw=1.6)
axs[0].set_title("$B_1\\propto r^{-3}$ along the axis; $U=-\\mu_2\\cdot B_1$", fontsize=9.5, loc="left")
axs[0].set_xlabel("distance $r$ from the centre of magnet 1")
axs[0].set_ylabel("$B_1$ / (arb. units)"); axs[0].grid(False)
axs[0].annotate("$B_1(r)$", (2.5, 0.46), fontsize=9, color=C2)
# (b) ladder of inverse powers on log-log
axs[1].set_title("the whole model on one log–log picture", fontsize=9.5, loc="left")
dd = np.linspace(2.0, 5.5, 60)/100
for n, (lab, col) in {3: (r"$U,\;B\;\propto\;d^{-3}$", "0.45"), 4: (r"$F\;\propto\;d^{-4}$", C2),
                      5: (r"$\gamma\;\propto\;d^{-5}$", C1)}.items():
    axs[1].loglog(dd*1e3, dd**(-n)/dd[0]**(-n), color=col, lw=1.8, label=lab)
axs[1].set_xlabel("separation $d$ / mm"); axs[1].set_ylabel("normalised")
axs[1].legend(frameon=False, fontsize=9)
axs[1].text(4.3, 1.9, "each derivative\nadds one power", fontsize=8.5, color="0.35", ha="center")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_dipole.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 4. linearisation of F
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.2))
d0 = 1.0
r = np.linspace(0.55, 1.8, 300)/1.0
Fd = lambda rr: (d0/rr)**4
axs[0].plot(r, Fd(r), color=C1, lw=2, label="$F(r)=F(d)(d/r)^4$")
axs[0].plot(r, Fd(d0) - 4*Fd(d0)/d0*(r-d0), color="0.35", lw=1.4, ls="--", label="tangent $F(d)-\\gamma(r-d)$")
axs[0].axvline(d0, color="0.6", lw=0.8, ls=":")
axs[0].fill_between([0.85, 1.15], 0, 3.0, color="gold", alpha=0.18)
axs[0].text(1.0, 2.75, "validity band $|\\xi|\\ll d$", fontsize=8.5, ha="center", color="0.45")
axs[0].annotate("slope $=-\\gamma$", (1.18, Fd(1.0)*0.999+ -4*(1.18-1)*Fd(1)*0.25+0.3), (1.32, 1.1),
                fontsize=9, arrowprops=dict(arrowstyle="->", color="0.3"))
axs[0].plot([0.85, 1.15], [Fd(0.85), Fd(0.85)+ -4*(0.85-1)*(Fd(0.85))], "k:", lw=0)
axs[0].set_xlabel("separation $r$"); axs[0].set_ylabel("repulsive force $F$")
axs[0].set_xlim(0.55, 1.8); axs[0].set_ylim(0, 3.05); axs[0].legend(frameon=False, fontsize=8.5, loc="upper right")
axs[0].set_title("$F(r)$ convex: tangent slope defines the magnetic stiffness", fontsize=9.5, loc="left")
# (b) convexity asymmetry
xi = np.linspace(-0.2, 0.2, 200)*1.0
axs[1].plot(xi, Fd(1+xi) - Fd(1), color=C2, lw=2, label="real $F(d+\\xi)-F(d)$")
axs[1].plot(xi, -4*xi*Fd(1), color="0.35", lw=1.3, ls="--", label="linear $-\\gamma\\xi$")
axs[1].fill_between(xi, -4*xi*Fd(1), Fd(1+xi)-Fd(1), color=C2, alpha=0.12)
axs[1].text(0.045, 0.15, "pushing together\nadds MORE than pulling\napart subtracts\n $\\Rightarrow$ net stiffening", fontsize=8.5, ha="center", color=C2)
axs[1].set_xlabel("change of separation $\\xi$"); axs[1].set_ylabel("$\\Delta F$")
axs[1].set_xlim(-0.2, 0.2); axs[1].legend(frameon=False, fontsize=8.5)
axs[1].set_title("why the non-linearity raises $f_2$ (extension, Ch. 9)", fontsize=9.5, loc="left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_linear.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 5. normal modes: mass-spring cartoons
fig, axs = plt.subplots(1, 2, figsize=(9.0, 2.9))
for ax, (lab, s1, s2, sub) in zip(axs, [("in-phase mode", +1, +1, "$\\xi=0$: gap fixed $\\Rightarrow$ $\\gamma$ does nothing\n$\\omega_1^2=k/m$  (independent of $d$)"),
                                         ("anti-phase mode", -1, +1, "$\\eta=0$: gap changes by $2a$ $\\Rightarrow$ $\\gamma$ is a 3rd spring\n$\\omega_2^2=(k+2\\gamma)/m$")]):
    ax.set_xlim(-3.3, 3.3); ax.set_ylim(-1.3, 1.3); ax.axis("off"); ax.grid(False)
    for wall in (-3.05, 3.05):
        ax.plot([wall, wall], [-1.0, 1.0], color="k", lw=2)
        for yy in np.linspace(-0.9, 0.9, 7):
            ax.plot([wall+0.16*(-1 if wall < 0 else 1), wall], [yy-0.14, yy], color="0.4", lw=0.8)
    x1, x2 = -1.25 + 0.55*s1, 1.25 + 0.55*s2
    def spring(xa, xb, col, lw=1.4, n=6, amp=0.16):
        xx = np.linspace(xa, xb, 120); yy = amp*np.sin(np.pi*n*(xx-xa)/(xb-xa)+np.pi/2*(n % 2))
        ax.plot(xx, (xx != xx)*0+yy*0 + np.sin(np.pi*n*(xx-xa)/(xb-xa))*amp, color=col, lw=lw)
    spring(-3.05, x1-0.32, "k"); spring(x2+0.32, 3.05, "k"); spring(x1+0.32, x2-0.32, C1 if s1*s2 > 0 else C2)
    for xx, s in ((x1, s1), (x2, s2)):
        ax.add_patch(Rectangle((xx-0.32, -0.32), 0.64, 0.64, fc="0.85", ec="k", lw=1.2))
        ax.annotate("", (xx+0.42*s, 0.75), (xx-0.12*s, 0.75), arrowprops=dict(arrowstyle="->", color="k"))
        ax.text(xx+0.1*s, 0.9, "$a$", fontsize=10, ha="center")
    ax.text(x1, 0, "$m$", ha="center", va="center", fontsize=10)
    ax.text(x2, 0, "$m$", ha="center", va="center", fontsize=10)
    ax.text(-2.15, 0.42, "$k$", fontsize=10, ha="center"); ax.text(2.15, 0.42, "$k$", fontsize=10, ha="center")
    ax.text((x1+x2)/2, -0.55, "$\\gamma$", fontsize=10, ha="center", color="k")
    ax.text(0, -1.08, lab, fontsize=11, ha="center", fontweight="bold")
    ax.text(0, -1.5, sub, fontsize=8.5, ha="center", color="0.25")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_modes.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 6. beats from two cosines (equal amplitudes)
fig, axs = plt.subplots(3, 1, figsize=(8.2, 4.4), sharex=True)
t = np.linspace(0, 2*40/1.07, 4000)
fa, fb = 6.49, 6.49+1.07
axs[0].plot(t, np.cos(2*np.pi*fa*t), color=C1, lw=1.1); axs[0].set_ylabel("$\\cos\\omega_1 t$")
axs[1].plot(t, np.cos(2*np.pi*fb*t), color=C2, lw=1.1); axs[1].set_ylabel("$\\cos\\omega_2 t$")
s = np.cos(2*np.pi*fa*t)+np.cos(2*np.pi*fb*t)
axs[2].plot(t, s, color="0.15", lw=1.1)
axs[2].plot(t, 2*np.cos(np.pi*(fb-fa)*t), "k", lw=1.0, ls="--")
axs[2].plot(t, -2*np.cos(np.pi*(fb-fa)*t), "k", lw=1.0, ls="--")
axs[2].set_ylabel("$x_1$\n(sum)")
axs[2].annotate("", (1/(fb-fa)*0.5, 2.45), (1.5/(fb-fa), 2.45), arrowprops=dict(arrowstyle="<->", color="0.3"))
axs[2].text(1.0/(fb-fa), 2.6, "$T_{\\rm beat}=1/(f_2-f_1)$", fontsize=9, ha="center", color="0.3")
axs[2].annotate("cosines in step\n→ maxima", (1/(fb-fa)/2, -0.3), (1/(fb-fa)/2-0.18, -1.6),
                fontsize=8.5, color="0.35", arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
axs[2].annotate("out of step\n→ cancellation", (1.0/(fb-fa), 0.25), (0.90/(fb-fa), -1.6), fontsize=8.5,
                color="0.35", ha="right", arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
for a_ in axs: a_.set_ylim(-2.7, 3.0); a_.axhline(0, color="0.7", lw=0.6)
axs[2].set_xlabel("time")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_beats.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 7. hold and release initial conditions
d0 = 0.0245; gam_m, w2 = modes(d0); gk = gam_m/w1**2
A = 5.0; x20 = A*gk/(1+gk)
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.4))
ax = axs[0]
ax.set_xlim(-4.4, 4.4); ax.set_ylim(-1.7, 4.6); ax.axis("off"); ax.grid(False)
ax.add_patch(Rectangle((-4.2, -1.55), 8.4, 0.45, fc="0.72", ec="k", lw=0.8))
ax.text(0, -1.33, "base (vertical = equilibrium line)", ha="center", va="center", fontsize=7.5, color="0.35")
# spring 1: equilibrium at x=-1.6 (dashed), held at -A
ax.plot([-1.6, -1.6], [-1.1, 3.7], ls="--", color="0.55", lw=0.9)
draw_spring(ax, -1.6, 3.6, -0.95, color=C1)
draw_magnet(ax, -2.55-0.28, 3.78, facing="N")
ax.annotate("", (-1.6, -0.55), (-2.55, -0.55), arrowprops=dict(arrowstyle="->", color=C1, lw=1.3))
ax.text(-2.08, -0.32, "$x_1=-A$ (held)", fontsize=9, ha="center", color=C1)
ax.text(-3.55, 1.8, "spring 1\npulled aside\n& held", fontsize=8.5, ha="center", color="0.3")
# spring 2: equilibrium at 1.6 (dashed), drifts left by x20 -> drawn at +0.95
ax.plot([1.6, 1.6], [-1.1, 3.7], ls="--", color="0.55", lw=0.9)
draw_spring(ax, 1.6, 3.6, -0.25, color=C2)
draw_magnet(ax, 1.35, 3.78, facing="S")
ax.annotate("", (1.6, -0.55), (1.35, -0.55), arrowprops=dict(arrowstyle="->", color=C2, lw=1.3))
ax.text(0.72, -0.32, "$x_2=-x_{20}$", fontsize=9, ha="center", color=C2)
ax.text(3.35, 1.8, "spring 2 relaxes\ntowards it: repulsion\nnow weaker than at\nrest ($r$ grew by $A$)", fontsize=8.5, ha="center", color="0.3")
ax.text(0, 4.35, "both springs at rest, then released", fontsize=9.5, ha="center", color="0.15")
ax.set_title("$kx_2=\\gamma(x_1-x_2)$ while held  $\\Rightarrow$  $x_{20}=A\\,\\gamma/(k+\\gamma)$", fontsize=9.5, loc="left")
ax.set_ylim(-1.7, 4.9)
# right: mode amplitudes bars + envelope sketch
ax = axs[1]
tb = np.linspace(0, 2/(1/ (w2/2/np.pi - f1)) , 900)
Tb = 1/(w2/2/np.pi - f1)
eta0, xi0 = -(A+x20), A-x20
env = 0.5*np.sqrt(eta0**2 + xi0**2 - 2*eta0*xi0*np.cos(2*np.pi*(w2/2/np.pi-f1)*tb))
ax.plot(tb, env, color=C1, lw=1.8, label="$|x_1|$ envelope")
ax.plot(tb, 0.5*np.sqrt(eta0**2+xi0**2 + 2*eta0*xi0*np.cos(2*np.pi*(w2/2/np.pi-f1)*tb)), color=C2, lw=1.4, ls="--", label="$|x_2|$ envelope")
ax.axhline(A, color="0.6", lw=0.7, ls=":"); ax.axhline(x20, color="0.6", lw=0.7, ls=":")
ax.text(2*Tb*0.98, A+0.15, f"max $=A=5$ mm", fontsize=8.5, ha="right")
ax.text(2*Tb*0.98, x20+0.15, f"min $=x_{{20}}={x20:.2f}$ mm", fontsize=8.5, ha="right")
ax.annotate("", (Tb/2, 0.6), (Tb/2, x20-0.05), arrowprops=dict(arrowstyle="->", color="0.35"))
ax.text(Tb/2-0.05, 0.28, "first minimum at $T_{beat}/2$", fontsize=8.5, ha="center", color="0.35")
ax.set_ylim(0, 5.9); ax.set_xlim(0, 2*Tb); ax.set_xlabel("time after release / s")
ax.set_title(f"envelopes:  starts at max, $\\xi_0/|\\eta_0|=(f_1/f_2)^2={ (f1/(w2/2/np.pi))**2:.3f}$", fontsize=9.5, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="center right")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_holdrelease.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 8. logic chain
fig, ax = plt.subplots(figsize=(8.6, 4.8)); ax.axis("off"); ax.grid(False)
ax.set_xlim(0, 10); ax.set_ylim(0, 5.6)
def box(x, y, w, h, txt, fc, fs=8.6):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec="0.35", lw=1.0, zorder=2,
                 joinstyle="round"))
    ax.text(x+w/2, y+h/2, txt, ha="center", va="center", fontsize=fs, zorder=3)
def arrow(x1, y1, x2, y2):
    ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="->", color="0.35", lw=1.2), zorder=1)
box(0.5, 4.6, 4.0, 0.7, "A. one leaf spring\nRayleigh shape $\\varphi$ → $k$, $m$", "0.93")
box(5.5, 4.6, 4.0, 0.7, "B. two dipoles\n$U(r)$ → $F(r)$, $\\gamma=-F'(d)$", "0.93")
box(2.6, 3.5, 4.8, 0.7, "C. Lagrangian  $m\\ddot x_i=-kx_i\\mp F$,  equilibrium,  linearise", "0.93")
box(2.6, 2.4, 4.8, 0.7, "D. normal modes  $\\eta,\\xi$ decouple", "0.93")
box(1.3, 1.3, 7.4, 0.75, "KEY RESULT  $f_2^2-f_1^2=\\dfrac{\\gamma}{2\\pi^2 m}=C\\,d^{-5}$,  $C=\\dfrac{3\\mu_0\\mu_m^2}{\\pi^3 m}$", "#d9ead9")
box(0.15, 2.55, 2.1, 0.55, "E. damping\nsame $\\lambda$ both modes", "0.97")
box(7.75, 2.55, 2.1, 0.55, "F. release\n$\\eta_0,\\xi_0$, beats", "0.97")
box(0.5, 0.15, 9.0, 0.75, "measurables:  two FFT peaks at $f_1,f_2$  ·  peak heights $\\xi_0/|\\eta_0|=(f_1/f_2)^2$  ·  beat period $1/(f_2-f_1)$  ·  equal decay", "#fdf1dc")
arrow(2.5, 4.6, 4.0, 4.2); arrow(7.5, 4.6, 6.0, 4.2); arrow(5.0, 3.5, 5.0, 3.1)
arrow(5.0, 2.4, 5.0, 2.05); arrow(1.2, 2.55, 2.4, 1.9); arrow(8.8, 2.55, 7.6, 1.9)
arrow(5.0, 1.3, 5.0, 0.9)
ax.text(9.55, 3.95, "extension:\nkeep $r^{-4}$\nexactly", fontsize=8, ha="right", color=C2)
ax.annotate("", (9.0, 3.85), (7.2, 3.85), arrowprops=dict(arrowstyle="->", color=C2, ls="--", lw=1.0))
fig.tight_layout(); fig.savefig(f"{FIG}/fig_chain.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 9. what the spectra look like (synthetic, like fig_release but pedagogical)
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.3))
fs = 240.0; t = np.arange(0, 45, 1/fs); dec = np.exp(-t/tau)
eta0, xi0 = -(A+x20), A-x20
w2l = w2
eta = eta0*np.cos(w1*t)*dec; xi = xi0*np.cos(w2l*t)*dec
x1 = 0.5*(eta - xi); x2 = 0.5*(eta + xi)
def spec(x):
    x = x - x.mean(); X = np.abs(np.fft.rfft(x*np.hanning(len(x)), 2**20))*2/np.sum(np.hanning(len(x)))
    return np.fft.rfftfreq(2**20, 1/fs), X
axs[0].plot(t[t<3.2], x1[t<3.2], color=C1, lw=0.9, label="$x_1$")
axs[0].plot(t[t<3.2], x2[t<3.2], color=C2, lw=0.9, label="$x_2$")
axs[0].set_xlabel("time / s"); axs[0].set_ylabel("tip displacement / mm")
axs[0].legend(frameon=False, fontsize=9); axs[0].grid(False)
axs[0].set_title(f"release at $d=24.5$ mm: raw traces ($f_1={f1}$, $f_2={w2l/2/np.pi:.2f}$ Hz)", fontsize=9.5, loc="left")
for x, col, lab in ((x1, C1, "$x_1$"), (x2, C2, "$x_2$")):
    fr, X = spec(x); mm = (fr > 5.5) & (fr < 8.8)
    axs[1].plot(fr[mm], X[mm]/X[mm].max(), color=col, lw=1.1, label=f"spectrum of {lab}")
axs[1].axvline(f1, color="0.55", lw=0.8, ls=":"); axs[1].axvline(w2l/2/np.pi, color="0.55", lw=0.8, ls=":")
axs[1].text(f1, 1.14, "$f_1$", ha="center"); axs[1].text(w2l/2/np.pi, 1.14, "$f_2$", ha="center")
axs[1].set_xlabel("frequency / Hz"); axs[1].set_ylabel("|FFT| normalised")
axs[1].set_ylim(0, 1.25); axs[1].legend(frameon=False, fontsize=9, loc="lower left")
axs[1].set_title("two peaks, same positions in both spectra", fontsize=9.5, loc="left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_spectra.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 10. damping: what λ changes / doesn't
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.3))
tt = np.linspace(0, 40, 800)
axs[0].plot(tt, np.cos(w1*tt), lw=1.0, color="0.7", label="undamped")
axs[0].plot(tt, np.exp(-tt/tau)*np.cos(w1*tt*np.sqrt(1-1/(tau**2*w1**2))), lw=1.4, color=C1, label="damped: $e^{-\\lambda t}\\cos\\omega' t$")
axs[0].plot(tt, np.exp(-tt/tau), "k:", lw=1); axs[0].plot(tt, -np.exp(-tt/tau), "k:", lw=1)
axs[0].set_title(f"damping shifts the frequency by $\\omega'-\\omega\\approx-\\lambda^2/2\\omega$: here $-5\\times10^{-6}\\,\\omega$ — invisible", fontsize=8.8, loc="left")
axs[0].legend(frameon=False, fontsize=8.5); axs[0].set_xlabel("time / s")
# lorentzian widths
f = np.linspace(6.2, 6.8, 500)
for wd, lab, col in [(0.005, "$\\tau_d=60$ s", C1), (0.0151, "$\\tau_d=21$ s", C2)]:
    axs[1].plot(f, 1/((f-6.49)**2+(wd/2)**2)/ (1/(wd/2)**2), color=col, lw=1.6, label=f"width FWHM $=1/\\pi\\tau_d$ = {wd*1000:.0f} mHz")
axs[1].set_title("the width (not the position) of each peak knows about damping", fontsize=9.5, loc="left")
axs[1].set_xlabel("frequency / Hz"); axs[1].set_ylabel("power spectrum (arb.)"); axs[1].legend(frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_damping_tut.pdf"); plt.close(fig)

# ------------------------------------------------------------------ 11. amplitude series & internal-resonance map
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.3))
d_arr = np.linspace(21, 42, 80)/1e3
gam_m_arr = 2*np.pi**2*C/d_arr**5
k_arr = w1**2
axs[0].semilogx(d_arr*1e3, 2*gam_m_arr/k_arr, color=C1, lw=2, label="$2\\gamma/k$ (coupling fraction)")
ratio = 1+2*gam_m_arr/k_arr
axs[0].semilogx(d_arr*1e3, ratio, color=C2, lw=1.4, ls="--", label="$f_2^2/f_1^2=1+2\\gamma/k$")
axs[0].axvline(24.5, color="0.5", lw=0.9, ls=":"); axs[0].text(24.7, 2.4, "closest $d$\nallowed", fontsize=8, color="0.4")
axs[0].axhline(3, color="0.4", lw=0.8, ls=":")
axs[0].text(38.5, 3.08, "$f_2=2f_1$ needs $2\\gamma/k=3$, i.e. $d^*\\approx16$ mm — unreachable", fontsize=8, color="0.4", ha="center")
axs[0].set_ylim(0, 4.2); axs[0].set_xlabel("separation $d$ / mm"); axs[0].legend(frameon=False, fontsize=8.5)
axs[0].set_title("the coupling dies fast: log-ladder $\\propto d^{-5}$", fontsize=9.5, loc="left")
# nonlin shift curves
for Av, ls in [(5, "-"), (10, "--")]:
    shift = (Av/(1+ (2*gam_m_arr/(k_arr+2*gam_m_arr))*0.5)/1e0)**2
    a = (Av*1e3) - (Av*1e3)*((gam_m_arr/k_arr)/(1+gam_m_arr/k_arr))
    g = (gam_m_arr/k_arr)/(1+2*gam_m_arr/k_arr)
    sh = (a*1e-3/d_arr)**2*(15*g/4-125*g**2/12)*100
    axs[1].semilogx(d_arr*1e3, sh, ls, color=C2, lw=1.7, label=f"release amplitude $A={Av}$ mm")
axs[1].axhline(0, color="0.5", lw=0.7)
axs[1].set_ylim(-0.4, 8); axs[1].set_xlabel("separation $d$ / mm"); axs[1].set_ylabel("$\\Delta f_2/f_2$ (%)")
axs[1].legend(frameon=False, fontsize=8.5)
axs[1].set_title("extension: non-linear shift $\\propto A^2 d^{-7}$ grows near the magnets", fontsize=9.5, loc="left")
axs[1].text(33, 5.6, "measure $f_2(A)$ at fixed $d$:\nslope must be $\\propto a^2$", fontsize=8, color="0.35")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_extension.pdf"); plt.close(fig)

print("figures written:", sorted(os.listdir(FIG)))
