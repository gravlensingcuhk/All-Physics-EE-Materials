// =====================================================================
//  THEORY, STEP BY STEP — teaching companion to
//  EE_essay/theory_v2/theory_v2.tex and theory_companion.tex
//  (coupled magnetic–mechanical oscillator).   Build: python3 compile.py
// =====================================================================
#set page(paper: "a4", margin: (top: 2.1cm, bottom: 2.3cm, x: 2.2cm))
#set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.68em)
#set heading(numbering: "1.")
#show heading.where(level: 1): set text(size: 13.5pt)
#show heading.where(level: 2): set text(size: 11.5pt)
#show figure.caption: set text(size: 9.2pt)
#show raw: set text(size: 9.5pt)
#set document(title: "Theory, step by step — coupled magnetic-mechanical oscillator")
#show figure: set block(above: 1.1em, below: 1.0em)

// ---------- boxes (call as: #keybox[Title][body]) ----------
#let tintbox(fill, bar, title, body) = block(
  width: 100%, fill: fill, inset: (x: 9pt, y: 7pt), radius: 2pt,
  stroke: (left: 2.8pt + bar, top: 0.6pt + gray.lighten(55%), right: 0.6pt + gray.lighten(55%), bottom: 0.6pt + gray.lighten(55%)),
  above: 1.0em, below: 1.0em,
)[
  #set text(size: 9.9pt)
  #set par(leading: 0.62em)
  #text(weight: "bold", title)\
  #v(0.4em)
  #body
]
#let boxed(body) = box(body, inset: 2pt, radius: 1pt, stroke: 0.7pt)
#let keybox(t, b)   = tintbox(rgb("#eef3fb"), rgb("#3b5f9e"), t, b)
#let warnbox(t, b)  = tintbox(rgb("#fbeeee"), rgb("#a83a3a"), t, b)
#let toolkit(t, b)  = tintbox(rgb("#f4eefb"), rgb("#7451a1"), t, b)
#let trybox(t, b)   = tintbox(rgb("#eef8ef"), rgb("#3d7a43"), t, b)
#let notebox(t, b)  = tintbox(rgb("#f5f5f5"), rgb("#808080"), t, b)
#let essaybox(t, b) = tintbox(rgb("#fff7ea"), rgb("#b07a2a"), t, b)

// =====================================================================
#align(center)[
  #text(size: 21pt, weight: "bold")[The Theory, Step by Step] \
  #v(0.9em)
  #text(size: 12.5pt)[A first-year's guided derivation of the coupled magnetic–mechanical oscillator] \
  #v(0.45em)
  #text(size: 9.6pt, fill: gray)[
    Teaching companion to #raw("theory_v2.tex") and #raw("theory_companion.tex") — every formula derived, every approximation named, every prediction tied to something measurable. \
    Pilot constants used throughout, so every claim is checkable: $f_1 = 6.49$ Hz · $C = 1.33 times 10^(-7)$ m⁵ s⁻² · $tau_d = 21$ s · $d = 24.5$–$42$ mm · release $A = 5$ mm.
  ]
]

#keybox[Read me first — how to use this booklet][
  The essay's theory is a *chain of six steps*; this booklet walks the chain in order, showing every line the essay compresses away:
  *A* a leaf spring is a mass $m$ on a spring $k$ (Ch. 4) →
  *B* two magnets are a force $F(d)$ whose *slope* $gamma$ is a second, tunable spring (Ch. 5) →
  *C* the Lagrangian turns the pair into two coupled equations, linearised by Taylor (Ch. 6) →
  *D* symmetry decouples them into two *normal modes*, and out pops the tested law $f_2^2 - f_1^2 = C d^(-5)$ (Ch. 7) →
  *E* damping (Kelvin–Voigt) says which spectral features to trust (Ch. 8) →
  *F* the hold-and-release start fixes the beat pattern, the FFT peak-height ratio and the beat minima (Ch. 9).
  Ch. 10 is the *extension* that drops the linearisation; Ch. 11 handles non-identical springs; Ch. 12 connects the model to real FFT spectra; Ch. 13 is the writing plan (what goes where, which pitfall lives where); Ch. 14 is an oral self-test; Ch. 15 a one-page cheat sheet.
  The purple boxes (Ch. 3) are the five maths tools used over and over — learn them once and the derivations read like sentences. Green boxes are *Try it yourself*: do them on paper; they are the exercises that convert reading into writing. Blue boxes are results worth memorising. Red boxes are traps, mostly from the companion's own pitfall list, with the reasons spelled out.
  *One sentence decides whether you understand this theory:* be able to say why the $-5$ exponent is a property of the *magnets* and not of the springs.
]

= The experiment, and the four things the theory must explain

Two identical copper leaf springs (thin elastic strips) are clamped vertically to a non-magnetic base, spaced so that their tips are $24$–$42$ mm apart. A disc neodymium magnet is glued to each tip, like poles facing (N–N), so the magnets *repel without touching*; the repulsion bows both springs outwards by $delta$ even at rest.

#figure(
  image("figs/fig_setup.pdf", width: 80%),
  caption: [The apparatus (house-style drawing; the essay's TikZ figure is equivalent). $L$: free length. $x_1, x_2$: horizontal tip displacements from the *bowed equilibrium*, both positive to the right. $d$: centre-to-centre magnet separation *at equilibrium* — the independent variable, set by moving one clamp; $d$ = gap between facing surfaces + one magnet thickness.],
) <fig-setup>

Pull one spring aside by $A$ (spring 1 to the left, so $x_1(0) = -A$), *hold* it until spring 2 hangs completely still, release (a withdrawn card). Video at 240 fps → Tracker → the two traces $x_1(t), x_2(t)$ → FFTs. Four facts must come out of the model:

+ each tip oscillates *fast* ($6.5$–$7.6$ Hz) while its amplitude *breathes* — beats — and everything decays with envelope time $tau_d approx 21$ s (Ch. 8–9);
+ the FFT of *either* spring shows *two* peaks, at $f_1$ and $f_2$ — the system's natural motions come in pairs (Ch. 7);
+ $f_1$ does not move as $d$ changes; $f_2$ climbs steeply at small $d$; and $f_2^2 - f_1^2$ falls as $d^(-5)$: a straight line of gradient $-5$ on log–log axes — *the* measurement the essay was built for (Ch. 7);
+ the breathing period is exactly $1/(f_2 - f_1)$; the beat minima sit at a computable height $x_20$ (not zero); the two springs' swellings interlace by half a beat; the two FFT peaks have height ratio $(f_1/f_2)^2$ — four free predictions the data can veto (Ch. 9).

Why this apparatus (the essay §1, sharpened): the coupling is a *field*, not a wire, so its stiffness $gamma$ is tuned continuously by one length $d$ — from "two lone oscillators" (large $d$) to "strongly locked pair" (small $d$); *two frequencies* are measured with $< 1$ % precision from a phone camera, and they encode the full force law via a *derivative*. The source problem is IYPT 2023, problem 6 ("Magnetic-mechanical oscillator"); see #raw("iypt.txt") and the slides for the research context. The theory below even predicts *which trials will fail* (beats slower than the record length) — a good theory bounds its own experiment.

#toolkit[Meta: the three-question habit][
  At every modelling step ask: (1) *what did we just assume?* (2) *what would look different if it were false?* (3) *which measurement is sensitive to it?* The essay's whole evaluation section is these questions answered for six assumptions (one-mode beam, point dipoles, small $xi$, identical springs, strip-only damping, ideal FFT), sorted by size. Keep that list running while you read Ch. 4–12; it is the skeleton of §9 of the essay.
]

= Symbols, and pilot numbers (keep this page open)

#figure(
  table(
    columns: (0.165fr, 0.50fr, 0.335fr),
    align: (left, left, left),
    table.header([*Symbol*], [*meaning*], [*pilot value / note*]),
    [$L, b, h$], [free length, width, thickness of one strip], [$L approx 150$ mm; $b, h$ calliper; $h$ enters *cubed*],
    [$E$; $I = b h^3/12$], [Young's modulus; second moment of area], [$E approx 110$ GPa (rolled strip: $+-10$ %)],
    [$mu$], [strip mass per unit length], [$mu L approx M/3$ here],
    [$M$], [mass of one magnet (+ glue)], [a few g],
    [$k$], [spring stiffness $k = 3E I/L^3 = E b h^3/(4L^3)$], [calibrated statically; model value $approx 33$ N/m],
    [$m$], [effective oscillating mass of spring + magnet], [$m = M + (33/140) mu L = k/(2 pi f_a)^2$],
    [$omega_0, f_a$], [natural frequency of one lone spring], [$f_a = 6.49$ Hz ("Infinity" trials)],
    [$mu_0$], [permeability of free space], [$4 pi times 10^(-7)$ T m/A],
    [$mu_m$], [magnetic moment of one magnet], [$approx 0.2$ A m²; datasheet: $mu_m = B_r V/mu_0$],
    [$d$], [centre-to-centre separation at equilibrium], [24.5–42 mm; $d = d_0 + 2 delta$],
    [$delta$], [static outward bow, $k delta = F(d)$], [$1.1$ mm at 24.5 mm — bigger than the ruler error],
    [$F(d)$], [dipole repulsion $3 mu_0 mu_m^2/(2 pi d^4)$], [$= k delta approx 0.04$ N at 24.5 mm],
    [$gamma$], [magnetic stiffness $-F'(d) = 4F(d)/d = 6 mu_0 mu_m^2/(pi d^5)$], [$0.18 k$ at 24.5 mm; $0.012 k$ at 42 mm],
    [$g$], [magnetic fraction of anti-phase stiffness $gamma/(k+2 gamma)$], [$= 1 - (f_1/f_2)^2$: 0.13 → 0.01 over the range],
    [$f_1, f_2$], [in-phase / anti-phase mode frequencies], [$f_2 = sqrt(f_1^2 + C d^(-5))$: 7.56 Hz at 24.5 mm],
    [$C$], [coupling constant $3 mu_0 mu_m^2/(pi^3 m)$], [$1.33 times 10^(-7)$ m⁵ s⁻² (pilot intercept)],
    [$eta, xi$], [normal coordinates: $x_1 + x_2$ (in-phase), $x_2 - x_1$ (gap change)], [$xi > 0$ = gap widened],
    [$eta_0, xi_0$], [mode amplitudes set by the release], [$eta_0 = -(A + x_20) < 0$; $xi_0 = A - x_20 > 0$],
    [$lambda$; $tau_d$], [damping rate $= c/(2m)$; envelope decay time $1/lambda$], [$tau_d = 21$ s; $lambda/omega = 1.2 times 10^(-3)$],
    [$tau_v$], [Kelvin–Voigt retardation time of the copper], [enters only as $k tau_v$ inside $c$],
    [$c_a$], [air-drag coefficient of a tip magnet], [the only mildly $d$-dependent damping piece],
    [$x_20$], [spring-2 drift while spring 1 is held: $A gamma/(k+gamma)$], [$0.76$ mm at 24.5 mm],
    [$f_"beat"$, $f_"avg"$], [$f_2 - f_1$; $(f_1+f_2)/2$], [$1.07$ Hz; $7.02$ Hz at 24.5 mm],
  ),
  kind: table,
  caption: [The complete symbol table. Numbers are the pilot constants of #raw("make_figs_v2.py"), so every worked example in this booklet can be re-checked; the *algebra* holds for any apparatus.],
) <tab-symbols>

#warnbox[Two conventions live in this repo — never mix them][
  This booklet and #raw("theory_v2.tex") use $xi = x_2 - x_1$, with both $x_i$ positive to the right. The older #raw("sec05_theory.tex") and its appendices use $xi = x_1 - x_2$. All final boxed laws are identical under either; intermediate signs are not. Pick one in your essay and keep it in every equation, caption and sentence (companion pitfall: "sign conventions").
]

= The toolkit — five maths tools, used over and over

#toolkit[T1 — the simple harmonic oscillator (the system's alphabet)][
  $m x'' = -k x$ i.e. $x'' + omega^2 x = 0$ with $omega = sqrt(k/m)$. Try $x = e^(r t)$: the characteristic equation $r^2 + omega^2 = 0$ has roots $r = +- i omega$; with Euler, $e^(+- i omega t) = cos omega t +- i sin omega t$, so the real general solution is
  $ x(t) = P cos(omega t) + Q sin(omega t) = A cos(omega t - phi). $
  The start fixes $P, Q$ — *never* $omega$. That is the whole reason normal modes matter: any motion is a sum of SHOs whose frequencies are properties of the system alone. Energy companion: $E = 1/2 m dot(x)^2 + 1/2 k x^2$ is conserved; equations from $E$ (or from $cal(L) = T - V$) are the same physics with sign bookkeeping built in.
]

#toolkit[T2 — Taylor expansion; the binomial series (linearising is this and nothing else)][
  $f(a + u) = f(a) + f'(a) u + 1/2 f''(a) u^2 + dots$ — the tangent line plus a named error. Two specialisations do all the work here:
  $ (1+u)^p = 1 + p u + p(p-1)/2 u^2 + dots $ — so $(1+u)^(-4) = 1 - 4u + 10u^2 - 20u^3 + dots$ (note: *negative* $p$ gives *positive* higher coefficients: the force curve is convex, and the corrections alternate sign); and the half-power $sqrt(1 + v) approx 1 + v/2 - v^2/8$ used for the weak-coupling beat formula (Ch. 7) and the damping shift (Ch. 8). We apply $(1+u)^(-4)$ twice — Ch. 6 drops everything after $u$, Ch. 10 keeps $u^2$ and $u^3$.
]

#toolkit[T3 — trig identities (the beat mechanism, literally)][
  $cos a + cos b = 2 cos((a+b)/2) cos((a-b)/2)$: the sum of two nearby cosines is one *fast* cosine at the average frequency whose amplitude is modulated by a *slow* cosine at half the difference — swellings recur at the full difference because the *modulus* of the slow factor doubles its frequency. Squaring and cubing fold power into new frequencies: $cos^2 theta = (1 + cos 2theta)/2$ (a DC part plus a second harmonic) and $cos^3 theta = (3 cos theta + cos 3theta)/4$ (a piece *at* the original frequency — which is exactly what shifts the frequency; Ch. 10).
]

#toolkit[T4 — complex amplitudes (two lines that replace a page of trig)][
  A real signal is the real part of a complex one: $P cos omega_1 t + Q cos omega_2 t = op(Re) ( e^(i omega_"avg" t) (P e^(-i Delta t/2) + Q e^(+i Delta t/2)) )$, $Delta = omega_2 - omega_1$, $omega_"avg" = (omega_1+omega_2)/2$. The first factor is the fast oscillation; the parenthesis is a *slowly turning vector* — its modulus is the *envelope*. $|P e^(-i theta) + Q e^(i theta)|^2 = P^2 + Q^2 + 2 P Q cos(2 theta)$ (one expansion). Keeping $P, Q$ *signed* through this line is the difference between a correct figure and a quarter-beat-wrong one (Ch. 9).
]

#toolkit[T5 — order-of-magnitude judgement (your evaluation superpower)][
  Correction size = coefficient $times$ small-parameter power. Example worked: the tip tilt angle is $theta = 3x/(2L) < 3 dot.op 5/(2 dot.op 150) = 0.05$ rad; tilt corrections to the force law carry factors $gamma/k approx 0.18$ and $(d/L)^2 approx 0.03$, relative to the mode stiffness: $approx 10^(-3)$; the FFT bin over 45 s is $0.022$ Hz $= 0.3$ % of $6.5$ Hz; so tilt is "below resolution — derived, estimated, dropped" (the honest phrase). Never write "negligible" without naming the yardstick.
]

#trybox[Ten minutes of pencil before the physics][
  (a) Check $omega_0 = sqrt(k/m)$ by dimensions: $(N/m) / "kg"$ gives what? (b) For $x = A cos omega t$ write $dot(x)$, form $E$ and verify it is constant. (c) First-order guesses: $(1.05)^(-4) approx 0.8$ vs exact? $(1.2)^(-4) approx 0.2$ vs exact? — at 5 % the tangent is within 3 %, at 20 % it is *off by more than half*. (d) With $f_1 = 6.49$, $f_2 = 7.56$ Hz: use T3 to name the fast frequency and the swelling frequency with no Fourier transform. *(Answers: (a) $1/s^2$; (c) $0.8227$ and $0.4823$ — keep (c) in mind when Ch. 6 talks about validity; (d) $7.02$ Hz and $1.07$ Hz.)*
]

= Step A — one leaf spring: $k$, $m$, and the 33/140 (theory_v2 §5.1)

#keybox[A — the results this chapter derives][
  $ k = 3 E I / L^3 = E b h^3 / (4 L^3), quad m = M + (33/140) dot mu L, quad omega_0^2 = k/m, $
  from the one assumed shape $phi(s) = (3 L s^2 - s^3)/(2 L^3)$, plus two by-products used later: the tip tilt $theta = 3x/(2L)$ and gravity's softening $Delta k = -(1.2 M + 0.375 dot mu L) g/L$.
]

== Why a beam may be collapsed to one coordinate

A copper strip is a continuum with infinitely many vibration shapes. The model's *first* approximation: only the lowest mode matters, and the lowest mode is (to excellent accuracy, with a heavy tip mass) the *static deflection curve* of a cantilever under a tip force (Rayleigh's method — Euler–Bernoulli beam theory is reviewed in the essay's App. A). So we separate variables with one shape and one amplitude:

$ y(s,t) = x(t) dot phi(s), quad phi(s) = (3 L s^2 - s^3)/(2 L^3), quad s in [0, L], $

and call $x(t)$ — the tip position — the whole state of the spring. Check $phi$ against the four conditions a clamped-free beam must satisfy: $phi(0) = 0$ and $phi'(0) = 0$ (clamp: no motion, no slope — both vanish like $s^2$ ✓); $phi(L) = 1$ ($(3-1)/2 = 1$ ✓, the normalisation that makes $x$ the tip displacement); $phi''(L) = 0$ (free end carries no bending moment: $phi'' = 3(L - s)/L^3$ ✓). That the same cubic also follows from $E I y'' = F(L - s)$ is §4.2.

#figure(
  image("figs/fig_shape.pdf", width: 100%),
  caption: [The assumed shape and where each energy lives. Kinetic energy weights with $phi^2$: the tip third dominates — that is why only $33/140$ of the strip's mass counts. Elastic energy weights with $(phi'')^2$: bending stress concentrates *at the clamp* — why leaf springs widen there, and why a nick at the clamp is fatal.],
)

== Statics: $k = 3E I/L^3$ from the bending equation

A tip load $F$ presses on the section at $s$ with bending moment $M(s) = F (L - s)$. Beam theory: curvature is moment over flexural rigidity, $E I y'' = M(s)$. The rigidity's geometry factor is the *second moment of area* of the rectangle about its mid-plane (the neutral plane, which neither stretches nor compresses when bent):

$ I = b integral_(-h/2)^(h/2) z^2 dif z = b h^3/12. $

Integrate $E I y'' = F (L - s)$ twice, clamped conditions $y(0) = y'(0) = 0$:

$ E I y' = F (L s - s^2/2), quad E I y = F (L s^2/2 - s^3/6), $

so the tip ($s = L$) stands off by $y(L) = F L^3 / (3 E I)$, and

$ boxed(k equiv F / y(L) = 3 E I / L^3 = E b h^3 / (4 L^3)). $

Two by-products of the same curve matter later. *Tilt:* the tip *slope* $y'(L) = F L^2/(2 E I) = (3/2) dot y(L)/L$, i.e. the magnet rotates by $theta = 3x/(2L)$ as it displaces (Ch. 5's tilt correction). *Gravity:* moving sideways, the tip also *descends*: $h_"drop" = 1/2 integral_0^L (y')^2 dif s = (3/5) x^2/L$ (with $y' = dot(x) dot 3 s(2L-s)/(2L^3)$, the integral evaluates to $6/(5L)$, half of which times $x^2$ is the drop). Every mass above the neutral line falls by (3/5)$x^2/L$: $Delta V = -(1.2 M + 0.375 mu L) g x^2/(2 L)$ — an inverted-pendulum "negative spring" of size a few per cent of $k$ at our frequencies, folded silently into the *statically calibrated* $k$ (theory_v2's notebox; say that in one sentence and move on).

== Elastic energy: the same $k$, one integral

Each fibre at height $z$ stretches by $epsilon = z y''$ (plane sections stay plane — the one beam axiom); $sigma = E epsilon$; energy density $1/2 sigma epsilon$ integrated over the section gives $1/2 E I (y'')^2$ per unit length:

$ V_"strip" = 1/2 E I x^2 integral_0^L (phi'')^2 dif s, quad integral_0^L (phi'')^2 dif s = 9/L^6 integral_0^L (L - s)^2 dif s = 3/L^3, $

so $V = 1/2 (3E I/L^3) x^2$ — *the same* $k$ as statics. This is not a coincidence: Rayleigh's method defines the effective stiffness as that energy ratio, so the two routes are one route; the *quality* of $k$ equals the quality of the shape ($< 0.2$ % error here; the exact comparison is the essay's App. A.5).

== Kinetic energy: the 33/140, done once, honestly

Element $dif m = mu dif s$ moves at speed $dot(x) phi(s)$, so

$ T_"strip" = 1/2 mu dot(x)^2 integral_0^L phi^2 dif s; quad
integral_0^L phi^2 dif s = 1/(4 L^6) integral_0^L (9 L^2 s^4 - 6 L s^5 + s^6) dif s = L/4 (9/5 - 1 + 1/7) = (33/140) L. $

(The polynomial square is the whole derivation — expand it yourself once; $(3 L s^2 - s^3)^2 = 9 L^2 s^4 - 6 L s^5 + s^6$.) Adding the magnet's $1/2 M dot(x)^2$:

$ T = 1/2 m dot(x)^2, quad boxed(m = M + (33/140) mu L). $

*Read the factor:* $33/140 = 0.236$: the strip moves so unevenly that under a quarter of its mass counts. (Companion check list: $33/280$ appears in #raw("sec05_theory") — it is the *same fact*, $T_"strip" = (33/280) mu L dot(x)^2$, half the $m$ coefficient; quote one, not both — pitfall.)

== Why the essay measures rather than computes

Pilot sanity: $k = 33$ N/m with $m = 21$ g gives $f = (1/2 pi) sqrt(k/m) = 6.3$ Hz — within 3 % of the measured $f_a = 6.49$ Hz; but since $k prop h^3$, a 5 % error in $h$ (a few vernier ticks on a 0.3 mm strip) is a 15 % error in $k$, and $E$ of a rolled strip is itself $±10$ %. So: *model for structure, calibration for numbers*: static $k$ by load–deflection; $m = k/(2 pi f_a)^2$ from a lone-spring pluck (one $f$-measurement — no $h$ anywhere, cf. theory_v2’s calibration note); and the exponent test $-5$ (Ch. 7) needs neither. That division of labour is a *design choice of the experiment* and a sentence your essay should own.

#trybox[Do it, then say it][
  (a) Reproduce $integral (phi'')^2 = 3/L^3$ and $integral phi^2 = 33L/140$ from $phi$ alone. (b) One paragraph in your own words: *why is $m$ not $M + mu L$?* (c) Nominal $h = 0.30 ± 0.01$ mm: what % range does $k$ inherit? (Answer: $±10$ %.) This is the arithmetic behind "we calibrated $k$ instead".
]

#essaybox[What the essay writes for A (companion writebox)][
  About half a page + the geometry figure: the boxed $k, m, omega_0$; one sentence on Rayleigh adequacy ("$< 0.2$ % with the tip mass; next mode $> 6 omega_0$, not excited"); one sentence absorbing gravity into the calibrated $k$; announce the tilt by-product (used in B). Full integrals go to the appendix (essay App. A already has them).
]

= Step B — two magnets: $U$, $F$, and the stiffness $gamma$ (theory_v2 §5.2)

#keybox[B — the results this chapter derives][
  $ U(r) = mu_0 mu_m^2 / (2 pi r^3), quad F(r) = 3 mu_0 mu_m^2 / (2 pi r^4), quad gamma equiv -dif F/dif r at r = d = 4 F(d)/d = 6 mu_0 mu_m^2 / (pi d^5). $
  $gamma$ — a *tunable, contactless* spring constant: the heart of the experiment.
]

== Dipole field, dipole energy, our geometry

A point dipole $mu$ creates $B(r) = mu_0/(4 pi r^3) dot [ 3 (mu dot hat r) hat r - mu ]$ (Griffiths §5.4 — the far field of any small loop or magnet; note the $1/r^3$ and that on the axis the field is $2 mu_0 mu/(4 pi r^3)$ *along* $mu$). A second dipole in a field has energy $U = -mu_2 dot(B_1)$ (a compass wants alignment). Combining:

$ U = mu_0 / (4 pi r^3) dot [ mu_1 dot(mu_2) - 3 (mu_1 dot hat r) dot (mu_2 dot hat r) ]. $

*Our geometry — like poles facing, moments along the joining axis pointing at each other* (N facing N): $mu_1 = mu_m hat r$ and $mu_2 = -mu_m hat r$, so $mu_1 dot(mu_2) = -mu_m^2$ and the axial product $= (mu_m)(-mu_m) = -mu_m^2$ too. Bracket: $-mu_m^2 + 3 mu_m^2 = +2 mu_m^2$:

$ boxed(U(r) = mu_0 mu_m^2 / (2 pi r^3) > 0 quad "repulsion") $

— positive energy that *falls* with separation; the magnets climb down the hill apart. (Had the moments been *parallel* — N facing S — the bracket would be $mu_m^2 - 3 mu_m^2 = -2$: attraction, $U < 0$: a compass aligning. Same formula, opposite social life.)

== Force and its slope: the ladder of inverse powers

$ F(r) = -dif U/dif r = 3 mu_0 mu_m^2 / (2 pi r^4); quad
gamma = -dif F/dif r = 12 mu_0 mu_m^2/(2 pi d^5) = 6 mu_0 mu_m^2/(pi d^5). $

And note the free identity $gamma = (4/r) F(r) = 4 F(d)/d$ at $r = d$ (power rule: derivative of $r^(-4)$ pulls down 4 and $F$ itself is $-dif U/dif r$; $gamma = -F' = 4F/d$). *What $gamma$ means:* displace the gap by $dif r$ near $d$: the repulsion *changes* by $F'(d) dif r = -gamma dif r$ — an extra restoring push proportional to $dif r$: exactly Hooke, with "spring constant" $gamma$. A repulsive force that weakens with distance therefore *stabilises* the anti-phase motion; an *attracting* configuration would give $gamma < 0$: the anti-phase mode *softens* (and below $k + 2 gamma < 0$ the springs buckle into contact — no small oscillations exist; the experiment would die. This is why "like poles facing" is in the apparatus list and not an option.)

#figure(
  image("figs/fig_dipole.pdf", width: 100%),
  caption: [Left: the head-on configuration and the sign bookkeeping (the bracket $mu_1 dot(mu_2) - 3(mu_1 dot hat r)(mu_2 dot hat r) = -mu_m^2 + 3 mu_m^2 = +2 mu_m^2 > 0$). Right: the ladder — field $ prop  r^(-3)$, energy $ prop  r^(-3)$, force $ prop  r^(-4)$ (one more power), stiffness $ prop  r^(-5)$ (one more): each differentiation takes one power; the experiment measures $gamma$, i.e. the *slope of the force*. This is why the essay's expected gradient is $-5$ and not $-4$.],
) <fig-dipole>

*Sanity arithmetic:* at $d = 24.5$ mm, pilot $gamma/k = 0.18$: the equilibrium bow $delta = F/k = gamma d/(4k) = 0.18 times 24.5/4 = 1.1$ mm ✓ (Ch. 2 table). And $gamma(42)/gamma(24.5) = (24.5/42)^5 = 0.066$: the coupling spans a factor 15 across the clamp range — that is the experiment's dynamic range, bought with two ruler marks.

*Datasheet moment:* a uniformly magnetised block has $mu_m = M_v V = B_r V / mu_0$. E.g. disc $phi 12 times 5$ mm³ at $B_r = 1.2$ T: $V = 565$ mm³ $=> mu_m = 1.2 times 5.65 times 10^(-7) / (4 pi times 10^(-7)) = 0.54$ A m². *One* number, $mu_m$, then sets $C$ (Ch. 7) — so the essay's *intercept* check is an absolute prediction from a datasheet line and a scale.

== The tilt refinement: derived, estimated, dropped

Because of Ch. 4, each moment *also rotates* with its tip, $theta_i = 3 x_i / (2 L)$. Redoing the bracket with tilted moments (theory_v2 eq. 5) gives

$ U(x_1, x_2) = mu_0 mu_m^2 / (4 pi r^3) dot [ 2 cos theta_1 cos theta_2 - sin theta_1 sin theta_2 ], quad r = d + x_2 - x_1, $

(the untilted bracket $=2$; theory_v2's "form used in most treatments"). Expanding to second order in $theta$ and re-running the whole mode analysis, theory_v2 §5.2 reports *relative* frequency shifts $-(9/32) dot (gamma/k) dot (d/L)^2$ and $-(3/32) dot g dot (d/L)^2$ — each below $10^(-3)$ for all $d <= 42$ mm, $L >= 150$ mm: below one FFT bin (Ch. 3 T5 did this arithmetic with you). The essay's App. B works the cycle-averaged version (the *mean* weakening $1 - theta^2/2$ of the effective moment, with $theta$ built from $delta + x$, is a few $times 10^(-4)$ in $F$). Write: "tilt corrections were derived, are negligible ($< 10^(-3)$ relative on both mode frequencies), and are dropped; they appear in the evaluation's list of neglected effects." That sentence is worth more than dropping silently.

#trybox[Do it, then say it][
  (a) Re-derive $U = +mu_0 mu_m^2/(2 pi r^3)$ from the general bracket without looking. (b) One line each: $F$ from $U$; $gamma$ from $F$; then $gamma = 4F/d$. (c) In two sentences: why is a *weakening* repulsion a *positive* stiffness? (d) If the essay's fit gave $gamma prop d^(-3.7)$ instead of $d^(-5)$, which assumption is the first suspect, and why is "the mode algebra" not on the list?
]

#essaybox[What the essay writes for B (½ page + equilibrium figure)][
  Boxed $U, F, gamma$; one "compressed spring" analogy sentence; the ladder sentence $r^(-3) → r^(-4) → r^(-5)$; the "measured, not computed, $d$" cross-reference to §1.3; one clause on point-dipole idealisation (finite-size model in the evaluation); one clause "tilt derived and negligible". Cite Griffiths.
]

= Step C — the Lagrangian, the equilibrium, and the linearisation (theory_v2 §5.3)

#keybox[C — the results this chapter derives][
  Exact (within A+B): $m x_1'' = -k x_1 - F(r)$, $m x_2'' = -k x_2 + F(r)$, $r = d_0 + x_2 - x_1$; equilibrium $k delta = F(d)$, $d = d_0 + 2 delta$.
  Linearised about equilibrium ($xi = x_2 - x_1$, $F(d+xi) - F(d) approx -gamma xi$):
  $ boxed(m x_1'' = -k x_1 + gamma (x_2 - x_1)) quad quad boxed(m x_2'' = -k x_2 - gamma (x_2 - x_1)). $
  *Two identical oscillators joined by a spring $gamma$.*
]

== The Lagrangian (why, and then in three lines)

Why bother with $cal(L)$ when Ch. 3's T1 was just $F = m a$? Because here the interaction force must appear on *both* masses with opposite signs (Newton 3) — a place where a careless $m a$ list goes wrong — and because the energy form is what extends (tilt, damping projection, nonlinear terms) without re-deriving anything. Lagrange for each coordinate $q$:

$ dif/dif t ( dif cal(L) / dif dot(q) ) = dif cal(L) / dif q. $

Inventory: $T = 1/2 m (dot(x_1)^2 + dot(x_2)^2)$ (Ch. 4, both tips; vertical motion is second order); $V = 1/2 k (x_1^2 + x_2^2) + U(d_0 + x_2 - x_1)$ (two springs from natural length + interaction). Hence, with $dif r / dif x_1 = -1$ and $dif r/dif x_2 = +1$:

$ m x_1'' = -k x_1 - dif U / dif r, quad m x_2'' = -k x_2 + dif U / dif r, quad dif U / dif r = -F(r), $

i.e. $m x_1'' = -k x_1 + F$... *watch the signs*: the chain rule gives $dif V/dif x_1 = (-dif U/dif r) dot.op (-1) = +F(r) dot.op ...$ — do it once with colours. $dif V/dif x_1 = (dif U/dif r)(dif r/dif x_1) = (-F)(-1) = +F$. So $m x_1'' = -k x_1 - F(r)$: magnet 1 pushed in $-x$; and $m x_2'' = -k x_2 + F(r)$: magnet 2 pushed in $+x$. *Equal and opposite, both outward* ✓ Newton 3 satisfied without being imposed. These equations are **exact within the two assumptions made so far** and are nonlinear because $F prop r^(-4)$.

== Equilibrium: why $d$ is what it is, and why it is measured

Set accelerations and velocities to zero: $k x_1 = -F$, $k x_2 = +F$: both springs bow *outward*, $|x_1| = |x_2| = delta$ with

$ k delta = F(d), quad d = d_0 + 2 delta. $

Solving $k (d - d_0)/2 = 3 mu_0 mu_m^2/(2 pi d^4)$ for $d$ given $d_0$ is a quintic; the apparatus definition *evades* it: $d$ is whatever the equilibrium gap is, and it is *measured* there (essay §1.3). Now *redefine coordinates*: let $x_i$ measure from the bowed equilibrium, so $r = d + (x_2 - x_1)$ and each spring's constant pre-tension force $k delta$ cancels $F(d)$: the equations of motion keep the form of §6.1 but with $F(d + xi) - F(d)$ (only the *change* of repulsion) as coupling. The system's dynamics never depends on the constant force — only its *slope* — which is the quiet reason everything below involves $gamma$ and not $F$.

#figure(
  image("figs/fig_holdrelease.pdf", width: 100%),
  caption: [Left: the *held* pre-history — with spring 1 displaced outward by $A$ (gap widened), the repulsion weakens and spring 2 drifts *toward* spring 1 by $x_20$ while at rest ($k x_20 = gamma (A - x_20)$... sign bookkeeping in Ch. 9). Right: $d$ vs $d_0$ with the $delta = F(d)/k$ bow that the definition absorbs — a 4 % effect at closest $d$ on $gamma$ itself if ignored ($gamma prop d^(-5)$: $(24.5/23.4)^5 = 1.26$, so a 1 mm slip in $d$ is 26 % in $gamma$: ruler care matters).],
) <fig-holdrelease>

== The linearisation, to the last decimal (theory_v2 eq. 13)

For $|xi| << d$, expand the force (T2):

$ F(d + xi) - F(d) = F'(d) xi + 1/2 F''(d) xi^2 + dots = -gamma xi + (5 gamma/2d) xi^2 + dots $

(using $F = F_0 (d/r)^4$ near $d$: $F' = -4F/d$, $F'' = +20 F/d^2 = 5 gamma/d$ — one line each by the power rule). Keeping the tangent only, the pair becomes the boxed system: the magnetic coupling *is* a spring of stiffness

$ gamma = -dif F/dif r at d = +dif^2 U/dif r^2 |_(d) $

(the second equality: curvature of the *energy* at the equilibrium is the stiffness — the general recipe for "small oscillations about a minimum"; Ch. 10 returns to the neglected $(5 gamma/2d) xi^2$ term). Equivalently the whole quadratic potential reads

$ V = 1/2 k (x_1^2 + x_2^2) + 1/2 gamma (x_2 - x_1)^2, $

three springs, two masses: "two copper springs to ground, one magnetic spring between" — and the linear eoms are $m x_i'' = -dif V/dif x_i$ in two lines. *Validity, honestly:* at $d = 24.5$ mm with $xi$ up to $4.2$ mm the neglected/retained *force* ratio is $[(5 gamma/2d) xi^2] / [gamma xi] = (5/2) xi/d approx 0.43$ at the swing extreme. Yet the *frequency* is disturbed only at order $(xi/d)^2$ with a small coefficient (Ch. 10 computes $+0.9$ %): the cycle-average of an even error shifts $omega^2$ via the mean stiffness — the odd part cancels over a period. Keep both numbers; the contrast (43 % force vs 1 % frequency) is the best one-sentence explanation of why linear models work: *errors that are antisymmetric across a cycle average away.*

#figure(
  image("figs/fig_linear.pdf", width: 100%),
  caption: [Left: the tangent replaces the curve — valid in the shaded band; at the closest trial the swing just touches the band's edge. Right: the tangent lies *below* the convex curve ($10u^2$ coefficient of §6.3, positive): squeezing the gap by $xi$ adds more restoring than stretching refunds, so the *average* stiffness of a big swing exceeds the tangent's — $f_2$ *grows with amplitude* (Ch. 10). The shaded area is a *cycle-averaged* number: $ prop  xi^2$, hence the $(a/d)^2$ law.],
) <fig-linear>

#warnbox[Pitfalls in this step (companion + ours)][
  (i) The *definition* chain: $gamma = -dif F/dif r$ carries a minus (so $gamma > 0$ for repulsion); writing $gamma = +dif F/dif r$ flips the coupling sign and "predicts" an *unstable* anti-phase mode — a great tell for a reviewer if you ever see $k - 2 gamma$ in a boxed formula. (ii) $d$ vs $d_0$: the $-5$ law applies at the equilibrium separation $d$ (which the experiment measures); relating to clamp spacing $d_0$ needs the quintic of §6.2 — the essay sidesteps it, and says so. (iii) Linearisation happens *after* shifting to equilibrium; linearising $U$ about the *unstressed* positions leaves a dangling $-F(d) xi$ and a wrong-looking "force at equilibrium". (iv) 33/280 vs 33/140: one statement, quoted once. (v) $k = 3E I/L^3$; the $48$ is for centre-loaded simply-supported beams — mis-citing Gere is a classic.
]

#trybox[Do it, then say it][
  (a) Write $cal(L)$, derive the exact pair and the linear pair without notes (energy route, not force route). (b) Show that *with identical springs* the constant $F(d)$ can never enter any *frequency*, only the equilibrium position — one sentence each for "why the experiment measures $d$ at equilibrium". (c) Reproduce the 43 %-force / 1 %-frequency pair and explain it in your own two sentences. (d) If $d$ were mis-measured by $+1$ mm at 24.5, what % error lands in $gamma$, hence in $f_2^2 - f_1^2$?
]

#essaybox[What the essay writes for C (¾ page)][
  The Lagrangian (one line); the exact equations (one line each); the equilibrium $k delta = F(d)$ and the "measure, don't quintic" remark; the Taylor step with the boxed linear pair; the quadratic-potential picture (three springs); the validity sentence with $|xi| << d$ and the promise "the extension keeps $F$ exact" (forward reference to §6). Full derivation lives in the companion's C1–C3; the essay main text *skips* the intermediate force-vs-chain-rule lines, and the appendix can hold them.
]

= Step D — normal modes, and the $d^(-5)$ law (theory_v2 §5.4)

#keybox[D — the results this chapter derives][
  With $eta = x_1 + x_2$ (in-phase) and $xi = x_2 - x_1$ (anti-phase):
  $ m eta'' = -k eta, quad m xi'' = -(k + 2 gamma) xi, quad
  boxed(omega_1^2 = k/m), quad boxed(omega_2^2 = (k + 2 gamma)/m). $
  Eliminating $k, m$ in favour of the measured frequencies and $gamma$:
  $ boxed(f_2^2 - f_1^2 = gamma / (2 pi^2 m) = C d^(-5)), quad C = 3 mu_0 mu_m^2 / (pi^3 m). $
  Log form: $ln(f_2^2 - f_1^2) = ln C - 5 ln d$ — *the essay's main graph*.
]

== Decoupling: add and subtract

Apply $+$ to the two boxed linear equations of Ch. 6: the coupling terms $+gamma xi$ and $-gamma xi$ annihilate (they are opposite by Newton 3 — *any* internal coupling dies in the sum):

$ m (x_1 + x_2)'' = -k (x_1 + x_2) quad=> quad m eta'' = -k eta. $

Apply $-$ (second from first): each mass now gets *both* coupling pushes in the same relative direction: $-gamma xi - gamma xi$:

$ m (x_1 - x_2)'' = -k (x_1 - x_2) - 2 gamma (x_2 - x_1) quad=> quad m xi'' = -(k + 2 gamma) xi. $

Two *independent* SHOs (T1). This worked because the system is symmetric under $1 ↔ 2$: the coordinates $(eta, xi)$ are the even/odd combinations, which every symmetric $2 times 2$ linear system admits; the matrix form is §7.4.

== The two modes, in words you can keep

*In-phase* ($xi = 0$, $x_1 = x_2$): the gap never changes, the magnetic force never changes, it does no work and adds no stiffness: $omega_1 = sqrt(k/m)$ — *identical to a lone spring, at any $d$*. This single sentence explains (essay §1.2, §5.4): (i) $f_1$ is $d$-independent — prediction (i) of the hypothesis; (ii) $f_1$ equals the single-spring calibration value — a *separate* check that both springs are really equal; (iii) later (Ch. 10): $f_1$ is even *amplitude*-independent, exactly — because the statement used only internal-force bookkeeping, never the linearisation.

*Anti-phase* ($eta = 0$, $x_1 = -x_2$): each tip at $a$ moves the gap by $2a$ (two tips, each contributing $a$); the coupling spring, compressed by $2a$, adds force $2 gamma a$ to *each* mass; effective stiffness per mass: $k + 2 gamma$. Extra stiffness ⇒ extra frequency: $f_2 > f_1$ always (for repulsion, $gamma > 0$), and *the only channel by which $d$ talks to the pair is $gamma(d)$*.

*Where the 2 lives, twice — derive, don't memorise:* (i) geometric: $xi = x_2 - x_1 = 2a$; (ii) energy: on the anti-phase line $V = 2 dot.op 1/2 k a^2 + 1/2 gamma (2a)^2 = (k + 2 gamma) a^2$ and $T = m dot(a)^2$ (two masses), so the $a$-oscillator has $omega^2 = (k + 2 gamma)/m$ ✓ (check the kinetic coefficient: $T = 2 dot.op 1/2 m dot(a)^2 = m dot(a)^2 = 1/2 (2m) dot(a)^2$ — mass $2m$; potential $1/2 (2k + 4 gamma) a^2$ — stiffness $2k + 4 gamma$; ratio cancels the $2$'s ✓). When both derivations agree, you own the factor of 2.

== The law, with every constant derived

$ omega_2^2 - omega_1^2 = 2 gamma/m $ — the *spring itself cancelled* — and with $omega = 2 pi f$:

$ f_2^2 - f_1^2 = gamma / (2 pi^2 m) = [ 6 mu_0 mu_m^2 / (pi d^5) ] / (2 pi^2 m) = underbrace( 3 mu_0 mu_m^2 / (pi^3 m) )_(= C) d^(-5). $

*Dimensional audit* (theory_v2 §5.5 does this; copy the habit): take $mu_0 = 4 pi times 10^(-7)$ N/A² exactly — the force-between-parallel-currents definition — then $[mu_0 mu_m^2] = (N A^(-2)) dot.op (A m^2)^2 = N m^4 = "kg" m^5 s^(-2)$, and dividing by $[m dot d^5] = "kg" dot.op m^5$ leaves $s^(-2)$: exactly the units of $omega^2$ ✓. $C$ therefore has m⁵ s⁻² ✓ (and $C d^(-5)$, s⁻² ✓). Do this audit whenever you meet a new constant — it costs 30 seconds and has caught more algebra slips than any other single check.

*Three consequences, which are literally the essay's hypotheses (§3):* (i) $f_1(d)$ flat; (ii) $ln(f_2^2 - f_1^2)$ vs $ln d$ straight with gradient exactly $-5$; (iii) intercept $ln C$ fixed by datasheet $mu_m$ and calibrated $m$ — an *absolute* check, no fitting constant. And a bonus reading: the *gap changes from $2a$ to $2 gamma a$* picture (essay §5.4 "gap changes... to $2 gamma a$") says the anti-phase mode's extra restoring is magnetic, so the *entire $d$-dependence of the spectrum* is carried by one power of the clamp position.

#figure(
  image("figs/fig_prediction.pdf", width: 100%),
  caption: [The predicted spectrum in both forms the essay draws (this is #raw("fig_prediction.pdf") from theory_v2, pilot constants): (a) $f_1$ flat (black), $f_2$ (copper) climbing below $approx 30$ mm; dotted marks $f_2 = 2 f_1$, the 2:1 internal resonance that would wreck the clean picture — *never reached* in $d >= 24.5$ mm (Ch. 10). (b) The ln–ln graph: gradient exactly $-5$ through the measurement separations (dots); the *nonlinear* corrections (Ch. 10) bend the true curve up at small $d$, up to $+0.9$ % (dashed: release $A = 5$ mm; dotted: $A = 10$ mm).],
) <fig-prediction>

*The pilot table* — compute every entry yourself from $f_2 = sqrt(f_1^2 + C d^(-5))$; these are the numbers used in the figures:

#figure(
  table(
    columns: (0.16fr, 0.14fr, 0.14fr, 0.14fr, 0.14fr, 0.12fr, 0.16fr),
    align: center,
    table.header([*$d$ / mm*], [*$f_2$ / Hz*], [*$Delta f = f_2 - f_1$ / Hz*], [*$T_"beat" = 1/Delta f$ / s*], [*beats in 10 s*], [*2 gamma/k*], [*$delta$ / mm*]),
    [24.5], [7.56], [1.07], [0.94], [10.6], [0.36], [1.10],
    [26.5], [7.23], [0.74], [1.35], [7.4], [0.24], [0.80],
    [28.5], [7.01], [0.52], [1.92], [5.2], [0.16], [0.55],
    [31.0], [6.84], [0.35], [2.87], [3.5], [0.11], [0.35],
    [33.5], [6.73], [0.24], [4.16], [2.4], [0.075], [0.21],
  ),
  kind: table,
  caption: [Predicted values from the boxed law with pilot constants ($f_1 = 6.49$ Hz, $C = 1.33 times 10^(-7)$ m⁵ s⁻², $k = 33$ N/m). The "beats in 10 s" column ($= 10 / T_"beat"$) shows why large-$d$ trials are resolution-limited (Ch. 12).],
) <tab-predictions>

== Optional but recommended: the matrix version (one general method)

Linear system $m dot(dot(bold(x))) = -sans(K) bold(x)$ with $sans(K) = mat(k + gamma, -gamma; -gamma, k + gamma)$ — diagonal $k + gamma$ ("my spring + the coupling's pull-back"), off-diagonal $-gamma$ ("the neighbour tugs"). Try $bold(x) = bold(X) e^(i omega t)$: eigenproblem $(sans(K) - omega^2 m I) bold(X) = 0$, nonzero $bold(X)$ when

$ det mat(k + gamma - omega^2 m, -gamma; -gamma, k + gamma - omega^2 m) = 0
quad=> quad (k + gamma - omega^2 m)^2 = gamma^2 quad=> quad omega^2 m = k + gamma -+ gamma, $

so $omega^2 in {k/m, (k + 2 gamma)/m}$ with eigenvectors $(1, 1)$ and $(1, -1)$ — the symmetry's ±. Keep the *recipe* (write the stiffness matrix, det $= 0$): Ch. 11 (non-identical springs) is the same two lines, and any future coupled-oscillator problem is $N$ lines of it.

= Step E — damping: Kelvin–Voigt, and what it does and does not move (theory_v2 §5.5)

#keybox[E — the results this chapter derives][
  Internal friction $sigma = E (epsilon + tau_v dot(epsilon))$ projects (with the *same* $integral (phi'')^2$ as Ch. 4) to a drag $-k tau_v dot(x)$ per spring; with air drag $-c_a dot(x)$ the normal coordinates obey
  $ q'' + 2 lambda q' + omega_i^2 q = 0, quad lambda = (k tau_v + c_a)/(2m) = 1/tau_d, $
  with solution $q = e^(-lambda t) cos(omega_i' t + phi)$, $omega_i' = sqrt(omega_i^2 - lambda^2)$. Consequences used by the analysis: frequencies safe to $approx lambda^2/omega^2 approx 10^(-6)$; both modes share one $tau_d$ (a *testable* statement); power-peak FWHM $= 1/(pi tau_d) approx 15$ mHz.
]

== What a Kelvin–Voigt solid is

A dashpot resists *rate* ($F = c v$); a spring resists *amount*. Put them *in parallel* — the Voigt cell — and both share the strain while stresses add: $sigma = E epsilon + c dot(epsilon) = E (epsilon + tau_v dot(epsilon))$ with the *retardation time* $tau_v = c/E$ ("the solid remembers how fast, for a while"). For bending, stress→moment and strain→curvature: $M_b = E I (y'' + tau_v dif/dif t y'')$, and the beam equation $mu dif^2 y/dif t^2 = -dif^2 M_b/dif s^2$ picks up a *viscous curvature* term.

Why does that become the *tip-coordinate* drag $-k tau_v dot(x)$? Because of the projection trick: with $y = x(t) phi(s)$, the viscous moment is the elastic moment with $x$ replaced by $tau_v dot(x)$ — every integral identical ($integral (phi'')^2 = 3/L^3$ again), so the energy-like route that produced $k x$ produces $k tau_v dot(x)$ for the drag. Add air drag on the magnet $-c_a dot(x)$ (Stokes-ish at these tiny Reynolds numbers; the essay keeps it as a lumped $c_a$). Per mass, per spring: total drag coefficient $c = k tau_v + c_a$, identical on both springs since the strips are identical. The coupling stays *conservative* (it has a potential) so adding/subtracting still decouples — each normal coordinate simply acquires the same $2 lambda dot(q)$, $lambda = c/2m$.

== Solve the damped SHO once, then never worry about frequency again

Characteristic equation: $r^2 + 2 lambda r + omega_i^2 = 0$, $r = -lambda +- sqrt(lambda^2 - omega_i^2)$; in our regime $lambda << omega_i$ the root is $-lambda +- i omega_i'$, $omega_i' = sqrt(omega_i^2 - lambda^2)$, so (real part) $q(t) = e^(-lambda t) cos(omega_i' t + phi)$: *exponential envelope, nearly the same period*. Now the three consequences:

+ *Peak positions are safe:* $omega_i' - omega_i = -lambda^2/(2 omega_i)$ (T2 on the sqrt): at $lambda/omega = 1.2 times 10^(-3)$, relative shift $-7 times 10^(-7)$ — invisible next to a 3 mHz FFT bin ($5 times 10^(-4)$ of $f$). Damping moves *widths*, not *positions*: the reason FFT frequencies can test the undamped theory at all (companion E2's point, "the peak *width* is set by damping, not the *position*").
+ *Width carries the decay:* a decaying cosine's power spectrum is a Lorentzian, FWHM in frequency $= 1/(pi tau_d)$ ($1/(pi tau_d) approx 15$ mHz; amplitude-spectrum half-max differs by $sqrt(3)$ — the essay's §6.1 convention: quote power FWHM and say so). Two peaks are resolvable iff $Delta f >= $ a few $times 1/(pi tau_d)$ — Ch. 12's resolution limit is *this* number dressed in window functions.
+ *Equal decay is a prediction:* both modes share $lambda$ *only because* strips are identical AND coupling is conservative. A mechanism damping the *relative* coordinate (eddy currents induced by one magnet in the other's field — the obvious suspect; squeeze-film air in the gap) would add to the $xi$ equation alone. Test the picture without new hardware: track the beat minima. Equal $lambda$ ⇒ minima/maxima ratio constant in time (the companion E2 line: "measure beat minima/maxima ratio vs time — constant ⇒ equal damping"). Gap-localised damping ⇒ minima rise with time. Theory_v2's damping figure draws both hypotheses.

#figure(
  image("figs/fig_damping.pdf", width: 100%),
  caption: [The envelope test (this is #raw("fig_damping.pdf") from theory_v2): both modes decaying at the same rate keeps the beat minima at a constant fraction of the maxima (left-style panel); anti-phase-selective damping lifts them with time. The pilot record matches the first (essay §8.1), which is a *real* check that the coupling is (near-)conservative — one more assumption tested rather than assumed.],
) <fig-damping>

*A measured value of $tau_d$ vs $d$ decomposes further:* $lambda = k tau_v/(2m) + c_a/(2m)$: the KV part is $d$-independent (a material constant of the copper); only $c_a$ varies with geometry — a (weak) monotonic trend with $d$. A $d$-independent measured $tau_d$ across trials *is therefore itself a result*: it says the damping is strip-internal (essay §8.1 measures $tau_d = 21.0 ± 0.4$ s for all five pilot trials). The 45 s records ($approx 2 tau_d$) are the sweet spot: enough decay for clean Lorentzians, not so much that late-time noise dominates the peak fit.

#trybox[Do it, then say it][
  (a) Solve $q'' + 2 lambda q' + omega^2 q = 0$ from scratch (ansatz $e^(r t)$) and read off the three consequences above. (b) Why can the FFT *frequency* not tell you $tau_d$ but the FFT *width* can — one sentence, then one sentence on what breaks that if the two tones are 2 bins apart? (c) Your record's beat minima grow linearly in time while maxima decay cleanly: name the mechanism and write the *one* modified equation that models it. *(Answers: (c) relative-coordinate damping: $xi'' + 2 lambda_"eff" xi' + omega_2^2 xi = 0$ with $lambda_"eff" > lambda$, eta-equation untouched.)*
]

= Step F — the release: beats, peak heights, and the sign lesson (theory_v2 §5.6)

#keybox[F — the results this chapter derives][
  Held-and-released start ($x_1 = -A$ held, spring 2 settled at $-x_20$, both at rest) gives
  $ eta_0 = -(A + x_20), quad xi_0 = A - x_20, quad x_20 = A gamma/(k + gamma), $
  hence, per tip (with the $e^(-lambda t)$ envelope on both modes),
  $ x_(1,2)(t) = 1/2 (eta_0 cos omega_1 t -+ xi_0 cos omega_2 t) e^(-lambda t); $
  FFT peaks: $|eta_0|/2$ at $f_1$, $xi_0/2$ at $f_2$, ratio $= (f_1/f_2)^2$; envelope $= 1/2 sqrt(eta_0^2 + xi_0^2 - 2 eta_0 xi_0 cos(Delta omega t))$: starts at its **maximum** $A$, first **minimum** $= x_20$ at $T_"beat"/2$; $T_"beat" = 1/(f_2 - f_1)$; and the always-valid identity $f_2^2 - f_1^2 = 2 f_"beat" f_"avg"$ (the $d$-law, measured from the *time domain*).
]

== The held state, and the initial conditions (theory_v2 Step F1)

While spring 1 is *held* at $-A$, spring 2 sits at the minimum of *its* energy in the displaced field: its spring pull $k x_2$ balances the (weakened) repulsion change $gamma(x_1 - x_2)$: $k x_2 = gamma(x_1 - x_2) = gamma(-A - x_2)$, i.e. $x_2 = -A gamma/(k + gamma) equiv -x_20$: *spring 2 drifts toward spring 1 by $x_20$* — sign: with magnet 1 farther away, less push on 2, so 2 relaxes inward toward its natural position... *draw it*. Pilot: $x_20 = 5 dot.op 0.18/1.18 = 0.76$ mm at 24.5 mm.

Release ($t = 0$): no damping has acted yet during the hold; both tips *start at rest*, so both modes start at their extremes with zero velocity — each mode is a *pure cosine* (T1), amplitude = initial value:

$ eta(0) = x_1 + x_2 = -(A + x_20) equiv eta_0; quad xi(0) = x_2 - x_1 = A - x_20 equiv xi_0; quad dot(eta) = dot(xi) = 0. $

*Note the signs* — $eta_0$ is negative (both tips displaced left), $xi_0$ positive (gap widened). They differ in size by $2 x_20$: the hold *splits* the excitation. If spring 2 were *also* moving when you released (a sloppy trial — the repo's data-checker flags exactly this), the pure-cosine start dies and every amplitude-based prediction below is void for that trial (frequencies survive).

== The two traces and their spectra (theory_v2 Step F2, eq. 20)

Invert $eta, xi$: $x_1 = (eta - xi)/2$, $x_2 = (eta + xi)/2$. Each trace: two cosines of known frequency, known amplitude, zero phase. *The FFT does the rest:* the amplitude spectrum of $A cos(2 pi f t) e^(-lambda t)$ is a Lorentzian of height $ prop  A$ at $f$ — so both records show *both* peaks, at $f_1$ with height $|eta_0|/2$ and at $f_2$ with height $xi_0/2$, and (because both modes share $lambda$) those heights *decay together*:

$ "peak ratio" = xi_0/|eta_0| = (A - x_20)/(A + x_20) = k/(k + 2 gamma) = (f_1/f_2)^2. $

Middle-to-right algebra: $A +- x_20 = A [1 +- gamma/(k+gamma)] = A (k + gamma +- gamma)/(k + gamma)$ — the ratio is $(k)/(k + 2 gamma)$ ✓ and Ch. 7's mode formulas give $(f_1/f_2)^2 = (k/m)/((k+2 gamma)/m)$ ✓✓. At $d = 24.5$ mm: $(6.49/7.56)^2 = 0.74$ — *the anti-phase peak should stand at 74 % of the in-phase peak*, a number that needs neither the release amplitude nor $k$ nor $m$ — parameter-free given $f_1, f_2$ alone. (The pilot table in the repo's #raw("results_DEMO") shows measured ratios 0.48–0.66 vs predicted 0.56–0.67 — the *variation* of the prediction across trials is itself $(f_1/f_2)^2$ evaluated at each trial's own $f_2$: never plug one trial's ratio into another's frequencies.) The one measured number here most sensitive to *how you released* — this is the essay §8's justification for the strict protocol.

== Beats and the envelope — with the signs carried (theory_v2 Step F3)

Write $x_1 = 1/2 op(Re) ( e^(i omega_"avg" t) ( eta_0 e^(-i Delta t/2) - xi_0 e^(i Delta t/2) ) )$, $Delta = omega_2 - omega_1$ (T4). The envelope is the modulus of the slow bracket; its square:

$ 1/4 [ eta_0^2 + xi_0^2 - 2 eta_0 xi_0 cos(Delta t) ] $

— a *cross-term carrying the product $eta_0 xi_0$ with a minus in front*. Since $eta_0 < 0 < xi_0$ the product term is *positive* at $cos = 1$: at $t = 0$ the square root collapses to $|eta_0| + xi_0$ and the envelope starts at $1/2 (|eta_0| + xi_0) = A$ — *maximum* at release ✓ (sanity: the spring was pulled to $-A$ and let go; it can never reach further out). The *first minimum* is where $cos(Delta t) = -1$, i.e. $t = pi/Delta = T_"beat"/2$, value $1/2 (|eta_0| - xi_0) = x_20$ — the swellings *never vanish*, because the start was unbalanced ($|eta_0| - xi_0 = 2 x_20 > 0$; perfect silence at the minima would need the push-pull start $|eta_0| = xi_0$ — theory_v2's second example). The squared envelope has period $2 pi/Delta$ and the modulus touches its extremes twice per period, so the *observable* swellings recur at $T_"beat" = 2 pi/Delta = 1/(f_2 - f_1)$ ✓ — no extra factor hiding anywhere.

*The sign lesson (companion §0, "the most common error"):* if you quote the identity in the memorised $2 cos((omega_1+omega_2)t/2) cos(Delta t/2)$ form with $|eta_0| = xi_0$ you get an envelope of $|eta_0| cos(Delta t/2)$, which vanishes at $t = pi/Delta$ *and at $3 pi/Delta$*, i.e. swellings every $2 pi/Delta$ — the same period, but *zero* minima instead of $x_20$ ones. If instead you drop the minus between the two mode terms (writing $+ xi_0$ where the geometry demands $-xi_0$), the *phase* shifts by a quarter-beat: minima appear where maxima belong — a figure error the companion reports as a real debugging story. The lesson generalises: **envelope algebra keeps the signs; only then does "minimum at release or maximum at release" come out of the derivation rather than of a wish.** The essay's Figure 3 and #raw("fig_release.pdf") have spring 1 at max and spring 2 at min at $t = 0$ — that is *this* derivation's result, and a figure you can regenerate (make_figs_v2.py panel A) and verify against §9.3's envelope formula.

*Energy transfer:* the *amplitude* of spring 2's anti-phase component is the same $xi_0$; spring 2's envelope is $1/2 sqrt(eta_0^2 + xi_0^2 + 2 eta_0 xi_0 cos Delta t)$: maxima when spring 1 has minima — a half-beat offset ✓ (theory_v2: "energy transfer at half-beat spacing"; the total energy $ prop  eta_0^2 + xi_0^2$ sloshes between springs without leaving the modes — *it is not dissipation*; the sum of squares of the two envelopes is constant). The beat minima depth gives a *third* handle on $gamma/k$ (via $x_20/A$) and the *time-domain* handle: $f_2^2 - f_1^2 = 2 f_"beat" f_"avg"$ *always* (difference of squares — no approximation; the weak-coupling form $f_"beat" approx C d^(-5)/(2 f_1)$ was Ch. 7's). Graph 4 of the essay measures $f_"beat", f_"avg"$ from the envelope and cross-validates the FFT peaks against $C d^(-5)$ *without a spectrum*.

#figure(
  image("figs/fig_beats.pdf", width: 100%),
  caption: [Equal-amplitude idealisation (T3): the fast oscillation $cos(omega_"avg" t)$ breathing with $|cos(Delta t/2)|$ — swellings every $2 pi/Delta$. The actual hold-and-release has $|eta_0| - xi_0 = 2 x_20 > 0$, so minima sit at $x_20$ (they never touch zero), and spring 2's envelope is the mirror image, shifted by half a beat. The half-factor question — call $f_"beat"$ the swelling rate or its half? — is convention; the physics is one sentence: *the observable swellings recur at $1/(2 pi/Delta) = f_2 - f_1$*. State the convention when you define $f_"beat"$ (the RQ uses the observable swellings: $f_2 - f_1$, never the half).],
) <fig-beats>

#warnbox[Pitfalls in this step][
  (1) $x_20 = A gamma/(k + gamma)$ — with $+ gamma$; the natural mis-derivation ("spring 2 also feels $gamma$ from the *grounded* equivalent"?) gives $k + 2 gamma$ — wrong: in the held static problem spring 2's *in-phase* combination does not run, no factor 2 exists. Balance forces on 2 directly: $k x_20$ vs $gamma(A - x_20)$ — that is the 3-line derivation of Ch. 6's held state; do it in the essay. (2) The envelope-sign story above. (3) Peak *ratio* $(f_1/f_2)^2$ requires the held-and-released start AND equal strips; for "push-pull" starts the *ideal* $|eta_0| = xi_0$ the in-phase peak would equal the anti-phase one — a different (also valid) test; know which start your data used. (4) $T_"beat"$ measured over the *first* beat only — after $2 tau_d$ the trace is at noise; the repo's demo data (10 s records, $tau_d approx 21$ s) is honest to a few beats, and the essay's longer 45–145 s records are what the improved method rides on.
]

= The extension: keep the full $r^(-4)$ law (theory_v2 §6; companion X1–X5)

#keybox[X — the results this chapter derives][
  Modal equations *exact*: $m eta'' = -k eta$ (unchanged, *any* force law) and
  $ m xi'' = -k xi + 2 [ F(d + xi) - F(d) ] $ — anharmonic. Lindstedt at $xi = a cos omega t + dots$:
  $ boxed( (omega - omega_2)/omega_2 = (a/d)^2 (15 g/4 - 125 g^2/12) ) $, $g = gamma/(k + 2 gamma)$ (positive: amplitude *stiffens* the anti-phase mode; $+0.9$ % at pilot release), a $2 f_2$ harmonic at $(5g/6)(a/d) approx 1.9$ %, a $3f_2$ line at $approx 0.15$ %, a mean *outward* separation shift $= (5g/2) a^2/d approx 0.24$ mm; $f_1$ exactly amplitude-blind; and the $f_2 = 2 f_1$ internal resonance needs $2 gamma/k = 3$ ⇒ $d^* approx 16$ mm, never reached.
]

== Do not linearise: the exact modal pair

The exact equations of Ch. 6 (equilibrium-shifted: each mass feels the *change* $-[F(d + xi) - F(d)]$ on 1 and $+[dots]$ on 2 — signs from §6.1's chain rule) split add/subtract *without any approximation*:

$ m eta'' = -k eta quad ("internal force: cancellation"), quad
m xi'' = -k xi + 2 [F(d + xi) - F(d)] quad ("they add"). $

Read the first box: *the in-phase mode is an exact SHO at $omega_1 = sqrt(k/m)$, for any amplitude and any force law* — "any $d$- or $A$-dependence of $f_1$ can only be measurement error, never physics" is theory_v2's strongest single line, and the essay uses it as a *systematic detector* (§6.3, 8.1). The second box is a Duffing-like oscillator with the full dipole nonlinearity; its linear part: $2 F'(d) xi = -2 gamma xi$ recovers Ch. 6 ✓ (sign check: $xi > 0$ (gap widened) ⇒ $F$ falls ⇒ the *excess* push-back on the inner-facing magnets is... $2[F(d+xi) - F(d)] < 0$: force on the $xi$ coordinate is negative — restoring ✓).

== Expand one order further: 5, 10 and the two small parameters

Binomial (T2) with $u = xi/d$: $F(d+xi) = F(d)(1 - 4u + 10u^2 - 20u^3 + dots)$, so

$ 2[F(d+xi) - F(d)] = -2 gamma xi + (5 gamma/d) xi^2 - (10 gamma/d^2) xi^3 + dots $

(the $-4$ term: $2 F(d)(-4 xi/d) = -8 F(d) xi/d = -2 gamma xi$ ✓; $10 F(d) dot.op 2/d^2 xi^2 = (5 gamma/d) xi^2$; $-20 dot.op 2 F(d) xi^3/d^3 = -(10 gamma/d^2) xi^3$). Dividing the $xi$-equation by $m$ and moving the linear part left ($omega_2^2 = (k + 2 gamma)/m$):

$ boxed(xi'' + omega_2^2 xi = omega_2^2 dot g dot ( 5 xi^2/d - 10 xi^3/d^2 + dots )), quad g = gamma/(k + 2 gamma). $

Two small parameters: the *swing* $a/d$ and the *magnetic fraction* $g$ (0.13 at 24.5 mm → 0.01 at 42 mm; the $g$ row of Ch. 2's table is exactly this $g$). The nonlinear drive is thus a $5 g / d$ (quadratic, *asymmetric*) and $10 g / d^2$ (cubic, symmetric) perturbation: the quadratic one knows about the convexity of fig @fig-linear right.

== Lindstedt–Poincaré in one page (companion X3)

*Method.* The nonlinearity shifts the period; instead of perturbing the *solution* at the linear frequency (which grows secularly — a term $ prop  t sin omega t$: unphysical), let the *frequency itself* adjust: $xi = a cos omega t + xi_2 + xi_3$, $omega = omega_2 + delta$, $a = O(1)$, $xi_2 = O(a^2/d)$, and demand the elimination of resonant forcing at every order ("no term at $omega$ may drive the $omega$-oscillator").

*Second order.* $5 g omega_2^2/d dot.op xi^2$ with $xi_1 = a cos omega t$: the $cos^2$ identity (T3) gives $(5 g omega_2^2 a^2/2d)(1 + cos 2 omega t)$. A constant drive shifts the equilibrium: mean shift $xi_"mean" = (5g/2)(a^2/d)$ — the *mean* separation widens (rectification: the magnets linger far apart where the force is weak). The $2 omega$ drive against the natural $omega$ responds with gain $1/(omega_2^2 - 4 omega^2) approx -1/(3 omega_2^2)$: $xi_2 = -(5g/6)(a^2/d) cos 2 omega t$. *Two results, the $1/6$ and the minus sign explained* (companion's exact wording) — the minus: a drive *above* resonance moves in opposition.

*Third order / the frequency shift.* Collect the resonant $cos omega t$ pieces produced by plugging the second-order solution (mean part $xi_"mean" = (5g/2)(a^2/d)$ and second harmonic $xi_"2h" = -(5g/6)(a^2/d)$, from the previous paragraph) back into the drives. Source (i), the cubic drive: $-10 g omega_2^2/d^2 dot.op (a cos omega t)^3$; with T3's $cos^3 = (3 cos + cos 3)/4$ its resonant share is $R_1 = -(15/2) g omega_2^2 a^3/d^2$. Source (ii), the quadratic drive $+5 g omega_2^2/d dot.op xi^2$ evaluated on the cross products of the fundamental with each correction: $2 xi_1 xi_"mean"$ contributes $2a dot.op (5g/2)(a^2/d) = 5 g a^3/d$, and $2 xi_1 (xi_"2h" cos 2 omega t)$ contributes (T3 again: $cos theta cos 2theta = (cos theta + cos 3theta)/2$) a resonant $a xi_"2h" = -(5g/6) a^3/d$. The projection of the mean-shift product keeps its factor 2, the harmonic product acquires $1/2$, so together: $R_2 = (5 g omega_2^2/d) dot ( 5 g a^3/d - (5g/12) a^3/d ) = +(125/6) g^2 omega_2^2 a^3/d^2$ (the $25 - 25/6 = 125/6$ split: $+25 g^2$ from $xi_"mean"$, $-25 g^2/6$ from $xi_"2h"$). The left-hand side, at the shifted frequency, contributes the mismatch $(omega_2^2 - omega^2) a approx -2 omega_2 delta omega dot.op a$, and a periodic solution whose amplitude stays $approx a$ must leave no resonant residue: $ -2 omega_2 delta omega a = R_1 + R_2 = omega_2^2 (a^3/d^2) ( -15 g/2 + 125 g^2/6 ), $ hence $ boxed(delta omega/omega_2 = (a/d)^2 (15 g/4 - 125 g^2/12)) $ — *the published law*. For $g <= 0.13$ the first term rules: positive ⇒ stiffer ⇒ frequency rises with amplitude (*convexity set the sign; Lindstedt set the coefficient*). The leftover $cos 3 omega t$ drives, divided by the response factor $omega_2^2 - 9 omega^2 approx -8 omega_2^2$, give the third-harmonic amplitude $+(5g/16 + 25 g^2/48) a^3/d^2$ of the key box. *Bookkeeping warning (this is where the $g^2$ coefficient is lost or mis-signed):* *both* the mean part and the second harmonic of the second-order solution cross the fundamental inside $xi^2$; keeping only one halves the $g^2$ term or flips it. The repo settles it numerically: #raw("nonlinear_check.py") integrates the exact $xi$-equation and gives $+0.88$ % where the boxed formula says $+0.94$ % (pilot release, $d = 24.5$ mm, $a = 4.24$ mm): sign and size confirmed; the small residual gap is higher order in $a/d$. *That is the habit this whole chapter teaches: every coefficient you cannot re-derive under exam conditions, you simulate once and then trust.*

*Third harmonic* comes from (i)'s leftover $(1/4) cos 3 omega t$ plus the $2 xi_1 xi_2$'s $cos 3 omega$ piece, with response $1/(omega_2^2 - 9 omega^2) approx -1/(8 omega_2^2)$: amplitude $(5g/16 + 25 g^2/48)(a^3/d^2)$ (collect the signs yourself as a check; the ratio to $a$ quoted in the box is what the simulated spectrum in fig @fig-nonlin's panel B marks with red bars).

== The measurable consequences — every one signed

#figure(
  image("figs/fig_nonlinear.pdf", width: 100%),
  caption: [The extension, verified (from #raw("make_figs_v2.py")). (a) $f_2$ shift vs swing amplitude at 24.5 and 31 mm: perturbation lines through direct-integration points (agreement a few % — the series is converging); dotted: the standard pilot release. (b) Simulated spectrum of $xi(t)$, log scale: second and third harmonics at $2f_2, 3f_2$ sit at the predicted $(5g/6)(a/d)$ and third-order levels (red = analytic).],
) <fig-nonlin>

1. *$f_2$ rises with amplitude; $f_1$ cannot.* At pilot release ($a = xi_0 = A - x_20 = 4.24$ mm; $a/d = 0.17$): $+0.94$ % at 24.5 mm $= +0.07$ Hz (about 1 FFT bin at 10 s, *inside* the current error bar), $+0.4$ % at 31 mm, $< 0.1$ % past 36 mm; doubling $A$ quadruples the shift ($a^2$ law). Consequences: the amplitude control of the essay's §4 is *theory-mandated*, not fussing; the measured ln–ln line should bend *upward* at small $d$ (residuals positive on the left) by $approx$ the shift size — a *signed* prediction; and an amplitude series $f_2(a)$ at fixed $d$ should be a straight line in $a^2$ with slope $(15g/4 - 125g^2/12) omega_2/d^2$ — the extension's most direct test. *Decay caveat:* over a 45 s record the swing relaxes $a: 4.2 → 1.2$ mm, so the FFT sits at a *time-averaged* shift, $approx 0.15$–$0.2$ of the initial-amplitude value; a two-tone fit of the *first beats* sees the full shift. (The companion's numerical note; both numbers appear in theory_v2 §6.5.) And: the amplitude test of $f_1$ (a "flat line with $A$") is a *systematics monitor* — it must be exactly flat; if it isn't, suspect peak pulling or leakage, not physics.
2. *Harmonics and DC.* In *both* $x_1$ and $x_2$ spectra: $2 f_2$ at $(5g/6)(xi_0/d) = 1.9$ % and $3 f_2$ at $approx 0.15$ % — *nothing at $2 f_1$ or $f_1 +- f_2$* (the in-phase motion is exactly harmonic; combination tones need the tilt coupling, $< 10^(-3)$ of that — "absent" is itself a prediction). The *mean shift* $(5g/2)(xi_0^2/d) = 0.24$ mm is DC — invisible in FFT (it shifts the pixel baseline, removed in detrending — one reason detrending is mandatory) but *visible in the equilibrium position*, in principle (essay: 2 px resolution vs $0.24$ mm $= 3$ px — "borderline"; the companion lists it as a possible future test).
3. *Internal resonance, by the way, is absent.* If $omega_2 = 2 omega_1$, the $2f$ harmonic of the anti-phase motion would drive the in-phase mode *resonantly* and the modal picture would mix at small amplitude. Condition: $omega_2^2/omega_1^2 = 1 + 2 gamma/k = 4$, $2 gamma/k = 3$, $d^* = (C/(3 f_1^2))^(1/5)$: pilot $d^* approx 16$ mm — inside the magnets' physical limit ($d >= 24.5$ mm, where $2 gamma/k = 0.36$). "The $f_2 = 2 f_1$ line of other IYPT groups' 'critical value' is at 16 mm here; this apparatus never reaches it, so the modal superposition stays valid at every measured $d$" — one paragraph that answers the *other* research question (the resonance hunt) by *avoiding* it. The condition also tells you the dipole model is itself questionable at $d^*$ (gap $≈ 10$ mm $< 2D$): the avoidance is doubly good.

#figure(
  image("figs/fig_extension.pdf", width: 100%),
  caption: [Ch. 10 in one view (from #raw("make_tut_figs.py")): left, the coupling ladder $2 gamma/k(d)$ and $f_2^2/f_1^2$ — the resonance condition $2 gamma/k = 3$ is off the top of the scale: the apparatus lives in the weak-coupling tail; right, the predicted non-linear $f_2$ shift vs $d$ for two release amplitudes: a measurement of the *slope in $A^2$* is the extension's direct test.],
) <fig-extension>

== The neglected-effects list (theory_v2 §6 closes with it; §9 inherits)

(i) *Finite magnets, not points:* the dipole law is the far-field expansion; the disc-vs-point correction is of order $(D/d)^2$ at our separations. The essay's charged-disc model (App. E) quantifies it: at surface gaps of 6–20 mm the true force lies a few per cent *below* the point-dipole curve and rises *less steeply* as the magnets approach (magnetisation spread over a volume cannot concentrate its flux at a point). A flatter force ⇒ flatter $gamma(d)$ ⇒ a log–log gradient *shallower than* $-5$: this is the leading suspect for the pilot's $-3.75$, with size corrections to $Y$ of 3–25 % at the closest separations. (ii) *Detuning:* $f_a$ and $f_b$ differ — Ch. 11: a $d$-independent floor that *flattens the ln–ln line at large $d$ only* (the pilot $f_a, f_b$ pin its size). (iii) *Tilt:* relative $10^(-3)$ on both mode frequencies — Ch. 5. (iv) *Gravity softening:* a few % of $k$, invisible because $k$ is calibrated with the strip vertical — Ch. 4. (v) *Higher beam modes:* above $6 omega_0$, not excited by a quasi-static release — Ch. 4. (vi) *Gap-localised damping:* would lift the beat minima — Ch. 8's envelope test says no. And one that *pushes the other way*: (vii) *nonlinearity* (this chapter) *steepens* the apparent gradient at small $d$ (up-bent residuals at closest approach) — the sign clash between suspects (i) and (vii) is why the evaluation can *adjudicate* rather than list: a measured gradient shallower than $-5$ with up-bending small-$d$ residuals still points at finite size. Suspects with signs, not a wish list — that is the evaluation chapter in embryo.

= When the two springs are not identical: detuning (essay App. C; companion D1–D2)

Pilot "Infinity" plucks: $f_a, f_b$ — suppose they differ by $Delta f / f approx 1$ %. The Ch. 6 linear system generalises with $m x_1'' = -k_1 x_1 + gamma xi$, $m x_2'' = -k_2 x_2 - gamma xi$; the ± combinations are no longer exact. Work in frequencies ($omega_a^2 = k_1/m$, $omega_b^2 = k_2/m$, $kappa = gamma/m$) and use Ch. 7's matrix recipe:

$ det mat(omega_a^2 + kappa - omega^2, -kappa; -kappa, omega_b^2 + kappa - omega^2) = 0. $

With mean $omega_m^2 = (omega_a^2 + omega_b^2)/2 + kappa$ and half-detuning $Delta = (omega_b^2 - omega_a^2)/2$: the determinant is $(omega_m^2 - omega^2)^2 - Delta^2 - kappa^2 = 0$, i.e.

$ boxed(omega_(1,2)^2 = omega_m^2 -+ sqrt(Delta^2 + kappa^2)), quad omega_2^2 - omega_1^2 = 2 sqrt(Delta^2 + kappa^2). $

In measured frequencies: $(f_2^2 - f_1^2)^2 = (f_b^2 - f_a^2)^2 + (C d^(-5))^2$ — *quadrature*: the magnetic term adds to the $d$-independent *floor* $f_b^2 - f_a^2$. Consequences: (a) large-$d$ flattening of the ln–ln line — the *same sign* as the resolution problem, but *predictable in size from the $f_a, f_b$ measurements alone* — measure the lone frequencies and you know how much of the flattening is detuning vs FFT (this is the whole point of the two calibration plucks: they *disentangle* — the essay's Graph-4/detuning figure); (b) the mode shapes mix: the lower mode tilts toward the *softer* spring, and the two records' peak-height ratios differ from the Ch. 9 equality (a measurable — the ratio *split* between the two spectra is a direct readout of $Delta$); (c) $f_1$ picks up a weak $d$-dependence (no longer perfectly blind). *Weak-coupling limit* ($kappa << |Delta|$, $d$ large): $omega_1^2 approx omega_a^2 + kappa (kappa/|Delta|)$-type shifts: each spring barely feels the other (the IYPT 2023 "weak coupling ⇒ $f_1, f_2 approx f_a, f_b$" statement); *strong coupling* recovers Ch. 7 ($Delta = 0$: $omega_m^2 = (omega_a^2 + omega_b^2)/2 + kappa$, $omega^2$ diff $= 2 kappa$ ✓). The essay's detuned-pilot data (numbers.tex: $f_a$ and $f_b$ 1 % apart, KFit $= 5.86$ vs $F("fitA") = 5.71$ — the detuning-aware fit shifts the coupling number) exercises exactly this formula; you only need to *own* the quadrature and the two limits.

= From theory to FFT: what the analysis can resolve (essay §6; Ch. 12 connects Ch. 7–9 to data)

The peaks that Ch. 7–9 predict are *what an FFT can or cannot deliver*; three numbers rule:

+ *Bin spacing $1/T$* (frequencies exist at $k/T$, $k$ integer — a 10 s record has 0.1 Hz bins; 45 s: 22 mHz). Two peaks are distinguishable iff $f_2 - f_1 >$ a few bins (rectangular window: "1/T" is the honest main-lobe half-width; Hann doubles the main lobe but kills sidelobes). Ch. 7's pilot table: at $d = 33.5$ mm, $f_2 - f_1 = 0.24$ Hz — just 2 bins at 10 s, 11 at 45 s: *the record length is a design parameter set by the theory*, not a habit. (The pilot's unresolved large-$d$ points: $T_"beat" = 4.2$ s vs 10 s record — the *time-domain* route $2 f_"beat" f_"avg"$ needs $approx 2$ beats: both methods fail at $d >= 36$ mm with 10 s records; the essay's §6 improvement (45–145 s + time-domain fit) extends the range.)
+ *Windowing & leakage:* truncation convolves the true spectrum with the window transform — rectangular's $sin x/x$ sidelobes pull a *nearby taller peak* on the smaller one: peak-pulling of up to half a bin when peaks are $approx 2$ bins apart (essay App. F simulates exactly the pilot's peak pair; the shift explains part of the large-$d$ scatter, and zero-padding does *not* fix it — interpolate the same wrong information more finely). A Hann window (×1.64 effective gain, $-31$ dB sidelobes) trades resolution for honesty; the essay §6.1 recommends peak *interpolation* (quadratic on log-magnitude) to beat the bin grid — $±0.1$ bin $= ±2$ mHz $= ±0.03$ %: fine for the $-5$ test, *not* fine enough for the 0.9 % non-linear shift: that needs the two-tone fit.
+ *Detrending:* a tracked-pixel offset (magnets not at same height, camera tilt) puts a DC pedestal under the peaks; the *spectral leakage of a step in the mean* is a $1/f$ skirt that biases close peaks asymmetrically — the essay subtracts mean/trend per trial; note the Ch. 10 *DC mean shift* (0.24 mm) would ride here — it is deliberately *removed* by detrending, since the mean-shift prediction is checked on the *held* displacement instead (Ch. 9's $x_20$) not the DC of the oscillating trace.

The falsifiable chain in data language: $f_1, f_2$ per trial (peak fit / two-tone fit) → $Y = f_2^2 - f_1^2$ → ln–ln gradient $-5 ± 1/T$-driven scatter; *plus*: the height ratio $→ (f_1/f_2)^2$; *plus*: beat minima $→ x_20/A$; *plus* (long records) $f_"beat" dot.op 2 f_"avg" → Y$ without spectra at all. Four measurements, three Ch. 7–9 formulas, *one* $gamma(d)$ — the summary table of theory_v2 §5.7 tabulates exactly these as predictions vs tests vs status; when you write your own §5 closing table, put the *equation number* in the middle column — that is what makes a review table a *contract*.

= Writing plan: main text, appendix, and the traps (companion §1–§3)

#essaybox[Main text of your theory section — five blocks, 3–4 pages][
  1 *One spring* (½ p + geometry fig): $T, V, k, m$ with *results only*, the by-product sentences (tilt, gravity→calibration). 2 *Two magnets* (½ p + equilibrium fig): $U, F, gamma = 4F/d$ boxed; the ladder sentence; $d$ defined at equilibrium (cross-ref §1.3); tilt "derived, negligible, dropped". 3 *Lagrangian & linearisation* (¾ p): exact eoms one line; equilibrium cancellation; Taylor step *explicitly* (one equation), validity $|xi| << d$ + forward reference to extension. 4 *Normal modes* (¾ p + modes fig): $eta, xi$; the two equations; boxed $omega_i^2$; *the key box* $f_2^2 - f_1^2 = C d^(-5)$ with the ln form; three predictions numbered as in §3 of your essay; coupling-size sentence ($2 gamma/k = 0.36$ at closest, $approx$ a fifth: one spring's stiffness). 5 *Damping + release* (¾ p + release fig + summary table): Ch. 8's three consequences in three sentences; Ch. 9's initial conditions with *signs shown*; the four predictions (ratio, minima, $T_"beat"$, $2 f_"beat" f_"avg"$); the predictions↔tests summary table (theory_v2 §5.7) closing the section.
]

#essaybox[Extension block (1–1.5 p) — the companion X-boxes][
  Exact modal pair (the $f_1$ amplitude-blindness *statement* earns its place); the $(5, 10, 20)$ binomial coefficients with a "convex ⇒ stiffening" sentence; *boxed* results only (shift, harmonics, mean shift — no Lindstedt machinery in main text; that goes to an appendix, exactly as theory_v2 does with eq. 22–27); the $f_2$ vs $A$ sentence ($+0.9$ % at closest — "comparable to but below present resolution, visible in the amplitude series"); the $d^* approx 16$ mm non-resonance paragraph; end with the neglected-effects list (Ch. 10's six suspects) so §9 inherits it *with numbers*.
]

#warnbox[The companion's pitfall list — each expanded to a habit][
  (1) $xi$-sign convention: one, always (pitfall of Ch. 2). (2) $f_"beat"$ convention: $f_2 - f_1$ (swellings), and $f_2^2 - f_1^2 = 2 f_"beat" f_"avg"$ *with* the factor 2, not 4 — the RQ's own algebraic footgun (numbers.tex was literally annotated about this). (3) The identity $(f_1/f_2)^2 = k/(k+2 gamma)$ for peak ratio holds *only* for the hold-release start on *identical* springs; state both conditions. (4) KV damping = strip material ⇒ equal mode decay; equal decay additionally needs conservative coupling — the "measure min/max ratio vs time" test exists for this. (5) If the fitted exponent is not $-5$, the suspect is the *dipole approximation for chunky magnets* (finite size) first — the mode algebra is exact and $k, m, E, h, I$ have cancelled from the gradient; do not defend the algebra, question the source. (6) $C$ uses *effective* mass $m$ (with 33/140), not $M$; $mu_m = B_r V/mu_0$ is per magnet, one. (7) $d$ is *defined* centre-to-centre at equilibrium; never quietly mix $d$ and clamp spacing $d_0$. (8) The $x_20$ factor 2 in $k + 2 gamma$ (modes) vs the $k + gamma$ in $x_20$ (held) — same coupling, different geometry; re-derive both per use. (9) Nonlinear frequency shift is *positive* (convexity); if your numbers say $f_2$ *drops* with $A$, look for a baseline/trend artefact in the trace. (10) $33/140$ or $33/280$, never both as separate facts. (11) $k = 3E I/L^3$, not the simply-supported $48$. (12) Envelope signs — derive (Ch. 9), the debugging story.
]

*Model sentences* (companion's "what to write", expanded; these are *skeletons* — rewrite in your voice, then check against the boxed equations you derived above):

- Beam: "Each spring is modelled as an Euler–Bernoulli cantilever of effective mass $m = M + (33/140) mu L$ and stiffness $k = 3E I/L^3$, its static deflection curve used as the fundamental-mode shape (Rayleigh); with the tip mass comparable to the strip this overestimates $omega$ by less than 0.2 %, and the next mode at $6 omega$ is not excited."
- Coupling: "Two coaxial dipoles facing like-pole-to-like-pole interact through $U = mu_0 mu_m^2/(2 pi r^3)$; for small changes of gap the force's slope acts as a spring of stiffness $gamma = 4F(d)/d = 6 mu_0 mu_m^2/(pi d^5)$ — tunable by the clamp position alone and *contactless*."
- Modes: "The equations decouple by the system's exchange symmetry: $m eta'' = -k eta$, $m xi'' = -(k+2 gamma) xi$. The in-phase mode never changes the gap, so its frequency is the single-spring value at *all* $d$; the anti-phase mode compresses the magnetic spring by twice the tip amplitude, so its frequency carries the coupling: $f_2^2 - f_1^2 = gamma/(2 pi^2 m) = C d^(-5)$ — the dipole gradient law measured through the *derivative* of a force, and free of every spring constant in the *gradient*."
- Release: "Holding one spring at $-A$ displaces the other by $x_20 = A gamma/(k+gamma)$ before release; both modes start at their extrema with $|eta_0| = A + x_20$, $xi_0 = A - x_20$. Each trace therefore shows *both* mode frequencies at heights $|eta_0|/2, xi_0/2$ — ratio $(f_1/f_2)^2$, parameter-free — while their beat has period $1/(f_2-f_1)$ with minima $x_20$ a half-beat later, and the two traces interlace by a quarter-beat: energy sloshing between the springs at exactly the beat rate."
- Damping: "Kelvin–Voigt internal friction adds $-k tau_v - c_a$ to each coordinate, so both modes decay at one rate $lambda = 1/tau_d$: peak *positions* are unaffected to $approx lambda^2/omega^2$, the *widths* are $1/(pi tau_d)$, and the constancy of the beat-minima depth in time tests the assumption that the coupling itself is conservative."
- Extension: "Keeping the $r^(-4)$ law exact leaves the in-phase mode *exactly* harmonic (so any $f_1$ amplitude dependence is systematic error), while the anti-phase frequency rises with the square of the swing, $(omega - omega_2)/omega_2 = (a/d)^2 (15g/4 - 125 g^2/12)$ ($+0.9$ % at the closest trial), generates $2f_2, 3f_2$ lines at 1.9 % and 0.15 % and a mean outward shift — and the 2:1 resonance ($2 gamma/k = 3$, $d^* approx 16$ mm) is never reached, so the modal picture holds over the whole clamp range."

= Self-test: the oral exam (answers at the end; be strict with yourself)

#trybox[Answer without notes; then check][
  + (1) Why may a copper strip be replaced by one coordinate, and *how good* is that? (2) Where does the 33/140 come from (one sentence + one integral)? (3) Why does $gamma$ scale as $d^(-5)$ when the force scales as $d^(-4)$ — which physical fact makes the experiment read the *gradient*? (4) Why must like poles face for the anti-phase mode to exist at all? (5) Do the Ch. 6 equilibrium cancellation in your head: what exactly removes the constant $F(d)$ from the *dynamics*? (6) Decouple the linear system two independent ways (forces; energy). (7) Prove $f_2^2 - f_1^2 = gamma/(2 pi^2 m)$ and then $= C d^(-5)$; audit units. (8) The release: write $eta_0, xi_0$ with signs, state the two FFT heights and the envelope extrema (with the sign argument!). (9) Why does equal decay of both modes *test* an assumption — which one, and against what alternative? (10) The extension: why is $f_1$'s amplitude-blindness *exact*, and what is the *sign* of $f_2$'s shift from convexity alone (no algebra)? (11) State the internal-resonance condition in terms of $2 gamma/k$ and evaluate with pilot $d = 24.5$ mm numbers. (12) Detuning: what does the ln–ln curve do at large $d$, and what *other* observable has the same sign?
]

*Answers (sketches — compare with the chapters).* (1) Rayleigh/static-curve ansatz; $< 0.2$ % with heavy tip, $1.5$ % bare; next mode $> 6 omega_0$. (2) $integral phi^2 dif s = 33L/140$ — kinetic-energy weighting by the shape square. (3) $gamma = -F'$: one more $1/r$ by the power rule; the modes respond to *force change*, i.e. the slope. (4) Attraction ⇒ $gamma < 0$ ⇒ $k + 2 gamma$ can reach 0: no stable anti-phase oscillation (soft-mode collapse into contact). (5) Spring pre-tension $k delta$ vs $F(d)$: only the *$xi$-dependence* of $F$ survives. (6) Add/subtract; or $V$ on the $x_2 = +- x_1$ lines. (7) $omega_2^2 - omega_1^2 = 2 gamma/m$; divide $(2 pi)^2$; insert $gamma(d)$; $C$: m⁵ s⁻² ✓. (8) $eta_0 = -(A+x_20)$, $xi_0 = A - x_20$; heights $|eta_0|/2, xi_0/2$, ratio $(f_1/f_2)^2$; envelope $1/2 sqrt(eta_0^2 + xi_0^2 - 2 eta_0 xi_0 cos Delta t)$: max $A$ at $t=0$, min $x_20$ at $T_"beat"/2$. (9) Equal $lambda$ rests on identical strips + conservative coupling; alternative: gap-localised damping ⇒ minima rise. (10) $m eta'' = -k eta$ used only internal-force symmetry; convexity: tangent below curve ⇒ mean stiffness rises ⇒ $delta omega > 0$. (11) $omega_2 = 2 omega_1 <=> 1 + 2 gamma/k = 4 <=> 2 gamma/k = 3$; pilot 0.36: never met. (12) Flattens toward the floor $sqrt(Delta^2 + kappa^2)$; peak-pulling / resolution has the same flattening sign (distinguish via $f_a, f_b$ measurements).

*The companion's checks table, reproduced for the impatient:* the integrals (33/140, $3/L^3$, $3/(2L)$ slope): exact polynomial; Rayleigh-vs-exact frequency: $<0.2$ %; $F, gamma$ ladder and $gamma = 4F/d$: differentiation; tilt: second-order series, $<10^(-3)$; binomial $-4$: series; *Lindstedt coefficients + full solution*: symbolic + numerical integration (agreement a few % at $a/d <= 0.17$, including the $g^2$ sign); peak-ratio identity: algebra (+ synthetic FFT: ratio $= (f_1/f_2)^2$ to 0.2 %); envelope formula: modulus of complex sum vs direct max of $|x_1|$; $x_20$ in held state: minimisation of full potential (repulsion $-$ spring at $x_1 = -A$, $x_2$ free).

= Cheat sheet (one page; when all else fails, rebuild from this)

#keybox[The chain][
  *Shape* $phi = (3 L s^2 - s^3)/2L^3$ → $integral phi^2 = 33L/140$, $integral (phi'')^2 = 3/L^3$ → \
  $k = 3E I/L^3 = E b h^3/4L^3$, $m = M + 33 mu L/140$, $omega_0^2 = k/m$ → \
  $U = mu_0 mu_m^2/2 pi r^3$, $F = 3 mu_0 mu_m^2/2 pi r^4$, $gamma = -F' = 4F/d = 6 mu_0 mu_m^2/pi d^5$ → \
  $k delta = F(d)$, $d$ at equilibrium → linear eoms: $m x_1'' = -k x_1 + gamma xi$, $m x_2'' = -k x_2 - gamma xi$ → \
  $eta'' + (k/m) eta = 0$; $xi'' + ((k+2 gamma)/m) xi = 0$ → $f_2^2 - f_1^2 = gamma/2 pi^2 m = C d^(-5)$, $C = 3 mu_0 mu_m^2/pi^3 m$ → \
  $ln(f_2^2 - f_1^2) = ln C - 5 ln d$; $f_1$ flat; $2 gamma/k = 0.36$ at 24.5 mm → \
  damping: $lambda = (k tau_v + c_a)/2m = 1/tau_d$ shared; FWHM $1/(pi tau_d)$; positions safe → \
  release: $x_20 = A gamma/(k+gamma)$; $eta_0 = -(A+x_20)$, $xi_0 = A - x_20$; $x_(1,2) = 1/2(eta_0 cos omega_1 t -+ xi_0 cos omega_2 t) e^(-lambda t)$; peak ratio $= (f_1/f_2)^2$; envelope $= 1/2 sqrt(eta_0^2 + xi_0^2 - 2 eta_0 xi_0 cos Delta omega t)$, min $x_20$ at $T_"beat"/2$; $T_"beat" = 1/(f_2-f_1)$; $f_2^2 - f_1^2 = 2 f_"beat" f_"avg"$ → \
  extension: $m eta'' = -k eta$ exact; $xi'' + omega_2^2 xi = omega_2^2 g(5 xi^2/d - 10 xi^3/d^2)$; $delta omega/omega_2 = (a/d)^2 (15g/4 - 125g^2/12)$; $2f_2$ at $(5g/6)(a/d)$; mean shift $(5g/2) a^2/d$; $f_1$: never; resonance needs $2 gamma/k = 3$ → $d^* approx 16$ mm (never); detuning: $(f_2^2 - f_1^2)^2 = (f_b^2 - f_a^2)^2 + (C d^(-5))^2$. \
] \

*Sources in the repo:* #raw("EE_essay/theory_v2/theory_v2.tex") (§5–6: the model) and #raw("theory_companion.tex") (the step chain, derivations A1–F2 + X1–X5, writeboxes, pitfall list, check table); the essay body #raw("sec05_theory.tex") and appendices A–D; #raw("sec06_methodology.tex") for FFT machinery; figures #raw("fig_setup/shape/dipole/linear/modes/holdrelease/beats/spectra/extension.pdf") are drawn for this booklet by #raw("make_tut_figs.py") (same copper/clamp grey–N-blue–S-red palette and the #raw("magnetSN")-style N–S striping as #raw("make_figs_v2.py")'s TikZ), while #raw("fig_prediction.pdf"), #raw("fig_release.pdf"), #raw("fig_damping.pdf"), #raw("fig_nonlinear.pdf") are theory_v2's own figures. The pilot constants are those of #raw("make_figs_v2.py"); measured values referenced from #raw("results_DEMO/summary.txt") and #raw("numbers.tex").

