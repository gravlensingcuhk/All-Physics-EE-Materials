# Physics EE – detailed structure plan
## "Investigation on the frequencies of oscillations of a coupled magnetic-mechanical oscillator"

This plan follows **Edwin's numbering and rhythm** (qualitative prediction → hypothesis → variables → methodology → data → theory → processed data → discussion → evaluation → appendices) and grafts on the **INTERNATIONAL essay's loop** (evaluation → *further* methodology → *further* data → *further* discussion → numerical simulation → conclusion). Your own idea is exactly right: main analysis with the **linearised** model (as in your S.T. Yau work), then the Evaluation exposes the linearisation as a systematic error and the "Further" sections switch to the **non-linear** treatment.

Word budget is for a 4000-word essay (equations, tables, captions, appendices, references do not count). Files referenced as `figs/…`, `ee_simulations.py`, `analyse_tracks.py`, `finite_size_magnets.py` are in this folder.

---

## 0. One-paragraph storyline (write this first, keep it on your wall)

> Two identical copper leaf springs with repelling magnets form two coupled oscillators. The in-phase mode (f₁) should not feel the magnets; the anti-phase mode (f₂) is stiffened by the magnetic coupling γ. Modelling the magnets as point dipoles gives **f₂² − f₁² = 3μ₀μₘ²/(π³ m d⁵)**, i.e. a −5 power law. I measured f₁, f₂ by video tracking + FFT for d = 21–60 mm. The ln–ln gradient was −3.75 ± 0.19, f₁ drifted by 2.6 %, and points beyond 35 mm could not be resolved. Evaluation identifies four systematic effects (FFT resolution/leakage, amplitude-dependent non-linearity of the r⁻⁴ force, finite magnet size, non-identical springs). A modified method (normal-coordinate tracking, longer records, small release, time-domain two-tone fitting, non-linear least squares) plus an RK4 simulation of the full non-linear equations quantifies each effect and recovers the dipole law within its range of validity.

---

## 1. Research question (sharpen the current title)

**RQ:** *How does the equilibrium centre-to-centre separation d of the magnets in a coupled magnetic–mechanical oscillator affect its two normal-mode frequencies f₁ and f₂, and does the coupling term f₂² − f₁² follow the inverse-fifth-power law predicted by the magnetic-dipole model?*

Why this wording: it names the IV (d), the DVs (f₁, f₂), the derived quantity that is actually tested (f₂² − f₁²) and the theoretical claim (d⁻⁵) — Criterion A wants exactly this. Keep the title as it is; put the RQ in §1.3.

---

## 2. Mapping: Edwin's structure → yours → INTERNATIONAL additions

| # | Edwin | Yours (title) | Words | Extra from INTERNATIONAL |
|---|---|---|---|---|
| 1 | Introduction 1.1–1.3 | 1.1 Coupled magnetic-mechanical oscillator (done) · 1.2 Modes of oscillation and beats · 1.3 Motivation & RQ | 350 | motivation paragraph (§1.4 style) |
| 2 | Qualitative prediction | 2 Qualitative prediction | 70 | |
| 3 | Hypothesis | 3 Hypothesis | 50 | |
| 4 | Variable table | 4 Variables (table) | 60 | list FFT settings as controlled variables |
| 9 (Edwin puts theory late) | Theoretical model | **5 Theoretical model** (put it *before* methodology, as INTERNATIONAL does) | 550 | 5.4 "Remarks on the model" incl. dimensional check and list of assumptions |
| 5–7 | Methodology, data collection, procedures | 6 Methodology: 6.1 FFT background · 6.2 Set-up & apparatus · 6.3 Method of data collection (why video + FFT) · 6.4 Procedures | 650 | 6.1 mirrors INTERNATIONAL §6.1 (DFT, resolution, leakage, windows) |
| 8 | Data | 7 Raw data | 200 | half-range uncertainty & justification |
| 10 | Processed data | 8 Analysis: 8.1 processed table · 8.2 ln–ln graph · 8.3 Y vs d⁻⁵ graph · 8.4 comparison with theory | 450 | |
| 11 | Discussion & conclusion | 9 Discussion (preliminary) | 250 | |
| 12 | Evaluation | 10 Evaluation: systematic (10.1–10.5) and random (10.6) | 550 | each error: sign, size, fix |
| — | — | **11 Further methodology** (modified procedure) | 200 | INTERNATIONAL §11 |
| — | — | **12 Further data & discussion** (windows / methods / models compared) | 350 | INTERNATIONAL §12–13 |
| — | — | **13 Numerical simulation** (RK4 of the non-linear equations) | 300 | INTERNATIONAL §14 |
| 11.2, 13 | Conclusion, extension | 14 Conclusion & extensions | 220 | |
| A–E | Appendices | A–H (see §10 below) | 0 | |

Total ≈ 4000. If you overrun, shrink 6.1 and 9 first; never shrink 10–13 (evaluation is the highest-weighted criterion in the new EE model — 8 of 30 marks for "Discussion and evaluation" [1](https://bespokelearning.io/blog/how-to-write-ib-extended-essay-2027); confirm with your supervisor which guide applies to your session).

---

## 3. Section-by-section content

### 1 Introduction
**1.2 Modes of oscillation and beats** (≈150 w). Same ideas as Edwin's 1.2 but for translational displacements x₁, x₂ of the magnets: in-phase mode (f₁, magnets move together, separation constant, magnets "invisible"), anti-phase mode (f₂, separation oscillates, magnetic force does work). Any release is a superposition ⇒ beats with f_avg = (f₁+f₂)/2 and f_beat = (f₂−f₁)/2. Include the identity **4 f_avg f_beat = f₂² − f₁²** — it ties your test quantity to Edwin's avg/beat language and explains why beats vanish at large d (f₂→f₁). Figure 2: your own sketch of the two modes (do not copy IYPT slide images). Figure 3: one of your own x(t) traces showing beats.

**1.3 Motivation & RQ** (≈150 w). Origin: IYPT 2023 problem 6 "Magnetic-mechanical oscillator" (cite the problem statement); analogy to magnetically coupled resonators / contactless coupling; personal angle (you can tune the coupling continuously by moving the springs, unlike a fixed spring). State the RQ.

### 2 Qualitative prediction · 3 Hypothesis
Repulsion ∝ d⁻⁴ ⇒ the stiffness added to the anti-phase mode, γ = −dF/dd ∝ d⁻⁵, collapses quickly with d, so f₂ → f₁ and the beat period → ∞; f₁ stays constant.
**Hypothesis:** f₁ independent of d; f₂² − f₁² ∝ d⁻⁵, i.e. ln(f₂² − f₁²) vs ln d is a straight line of gradient −5.

### 4 Variables (table, Edwin style)
IV: d (centre-to-centre, magnets mounted, measured at equilibrium). DV: f₁, f₂ (→ f₂² − f₁²). Controlled: same spring pair, clamped length L, clamp torque, same magnet pair & orientation, initial displacement (magnitude and which spring), release method, camera fps/resolution/distance, record length T, FFT settings (window, N, zero-padding) — INTERNATIONAL lists FFT settings as controlled variables; do the same. Uncontrolled: air currents, temperature (E of copper), residual clamp slip.

### 5 Theoretical model (≈550 w — the core of Criterion B/C)
5.1 *Leaf spring as an oscillator.* Cantilever tip stiffness k = 3EI/L³; effective mass m = M + (33/140)μL (Rayleigh) ⇒ ω₀² = k/m. (Derivation → Appendix A, or cite.) One sentence: the Euler–Bernoulli continuum model of the IYPT slides reduces to this for a heavy tip mass; a finite-element check shows < 1 % difference in f₂² − f₁² for M/(μL) ≥ 0.4 (Appendix A, ready-made numbers below in §9).
5.2 *Magnetic force.* Coaxial identical dipoles, like poles facing: F(r) = 3μ₀μₘ²/(2π r⁴), repulsive (from F = ∇(μ·B); Appendix B). Equilibrium: springs bend outward by δ with kδ = F(d). You measure d directly, so the IYPT quintic is not needed — say so.
5.3 *Linearised coupling and normal modes.* F(d + x₂ − x₁) ≈ F(d) − γ(x₂ − x₁), γ = −F′(d) = 6μ₀μₘ²/(π d⁵).
m ẍ₁ = −k x₁ + γ(x₂ − x₁), m ẍ₂ = −k x₂ − γ(x₂ − x₁). Normal coordinates η = x₁ + x₂, ξ = x₁ − x₂ (matrix form like Edwin's eq. 12–15):
ω₁² = k/m, ω₂² = (k + 2γ)/m ⇒ **f₂² − f₁² = γ/(2π² m) = 3μ₀μₘ² / (π³ m d⁵)**.
Linearised forms for graphs: (i) ln(f₂² − f₁²) = ln C − 5 ln d, C = 3μ₀μₘ²/(π³m); (ii) f₂² − f₁² = C·d⁻⁵ (straight line through the origin vs d⁻⁵ → compare its gradient with C from independently measured μₘ, m).
5.4 *Remarks* (INTERNATIONAL §5.3 style): dimensional check [μ₀μₘ²/(m d⁵)] = s⁻²; list the assumptions you will revisit: (a) point dipoles, (b) linearised force (small amplitude), (c) identical springs, (d) lumped mass, (e) no damping, (f) magnets stay coaxial, (g) the FFT returns the true frequencies.

### 6 Methodology (≈650 w)
**6.1 FFT background** (INTERNATIONAL §6.1 adapted — see §4 of this plan for the physics you need).
**6.2 Set-up & apparatus**: labelled photo (Edwin Fig. 3), apparatus list, magnet specs (grade, D × t, mass), spring dimensions (L, b, h), camera (model, fps, resolution), tracking software (Tracker / phyphox / Python).
**6.3 Method of data collection** — narrate like Edwin §6: Method 1 (read beat and average periods from x(t)) failed at large d because the envelope disappears when f₂ − f₁ is small; Method 2 (FFT of x(t)) adopted. Show one x(t) and its FFT with the two peaks labelled (Edwin Fig. 4–6).
**6.4 Procedures**: 6.4.1 spring stiffness k (static load–deflection, plot F vs x, LINEST — Edwin §7.1); 6.4.2 magnet moment μₘ (either datasheet B_r: μₘ = B_r V/μ₀, or better: put one magnet on an electronic balance, bring the other coaxially above it at measured r, force = Δm·g, fit F ∝ r⁻⁴ — this doubles as an independent check of the r⁻⁴ law); 6.4.3 masses M, μL; 6.4.4 f₁, f₂ vs d: set d, displace one spring by A = 10 mm, release, record T s at fs fps, track, FFT, read the two highest peaks, 3–5 trials per d.

### 7 Raw data (≈200 w)
Table 1 basic quantities (L, b, h, M, μ, k ± Δk, μₘ ± Δμₘ, m). Table 2 F vs x for k. Table 3/4 f₁ and f₂ per trial with mean and half-range (justify half-range vs SD as INTERNATIONAL §8.2 does: n = 3–5 too small for SD). **State T, fs, N and the bin width 1/T in the caption** and add a column "resolved? (f₂−f₁ ≥ 2/T)". Rows d = 40–60 mm: "peaks unresolved" — this is a result, not a gap.

### 8 Analysis (≈450 w)
8.1 Processed table: d, ln d ± Δd/d, Y = f₂² − f₁², ΔY = 2(f₁Δf₁ + f₂Δf₂), ln Y ± ΔY/Y, d⁻⁵ ± 5Δd/d·d⁻⁵. Sample calculation for one row (Edwin style).
8.2 Graph 1: ln Y vs ln d with error bars, best-fit and worst-fit lines (or LINEST if bars are too small). Report gradient −3.75 ± 0.19 (unweighted) and note the weighted value −3.48 ± 0.24 — the disagreement itself shows the ln-transform problem.
8.3 Graph 2: Y vs d⁻⁵ (theory: straight line through origin, gradient C). Graph 3: f₁ vs d (theory: horizontal) — you will see the 2.6 % drift.
8.4 Comparison table (Edwin Tables 13/14): theoretical gradient −5 (exact) vs experimental; theoretical C from μₘ, m vs experimental gradient of Graph 2; % differences.

### 9 Discussion — preliminary (≈250 w)
Trend confirmed (Y falls steeply, f₂ → f₁), but: exponent −3.75 ≠ −5; f₁ not constant; points beyond 35 mm unresolved; residuals systematic (curvature) rather than random. Frame the four candidate causes and hand over to §10.

### 10 Evaluation (≈550 w) — one sub-section per effect, each with *sign, size, remedy*
10.1 **FFT resolution & leakage** (methodological, procedural). Two tones separated by less than ~1/T (rect) or ~2/T (Hann) merge or "pull" each other; your f₂−f₁ = 0.039 Hz at 35 mm needs T ≥ 26–51 s, and 0.009 Hz at 40 mm needs T ≥ 110–220 s. Show the synthetic-spectrum figure (`figs/A_fft_resolution.png`) built from *your* f₁, f₂. Sign: biases the large-d end → flattens the gradient. Remedy → §11.
10.2 **Non-linear magnetic force** (model assumption b). Expand F(d − ξ): ξ̈ + ω₂²ξ + αξ² + βξ³ = 0 with α = 2.5 g ω₂²/d, β = 5 g ω₂²/d², g = (f₂² − f₁²)/f₂². Landau–Lifshitz: Δω/ω₂ = (a/d)²(15g/8 − 125g²/48), a = amplitude of ξ (Appendix D). With your A = 10 mm release: +7.6 % at 21 mm (RK4 exact: +9.8 %), +1.2 % at 35 mm (`figs/B_f2_vs_amplitude.png`). Sign: f₂ overestimated at small d → gradient *steeper* (−5.16 in simulation, `figs/B_lnln_gradient.png`), so it is a large error but not the cause of the flattening. Remedy: A ≤ 3 mm, analyse late (small-amplitude) segments, extrapolate f₂ to zero amplitude, RK4 correction (§13).
10.3 **Finite magnet size** (assumption a). Replace point dipoles by two uniformly "charged" discs per magnet (Gilbert model, Appendix E). For d only 2–4× the magnet diameter the near field is softer than d⁻⁵: apparent gradient over 21–35 mm is −4.9 (D = 6 mm), −4.5 (D = 12 mm), −3.75 (D = 20 mm) (`figs/D_finite_size_gradient.png`; run `finite_size_magnets.py` with *your* D, t). Sign: flattening — the leading candidate if your magnets are ≥ 12 mm. Remedy: exact cylindrical-magnet force in the model, or larger d (which collides with 10.1 — discuss the trade-off explicitly; examiners like that).
10.4 **Non-identical springs** (assumption c) — analogue of Edwin's 12.2/Appendix D. With uncoupled frequencies f_a ≠ f_b: ω₁,₂² = ω̄² + κ ∓ √(Δ² + κ²), κ = γ/m, Δ = (ω_b² − ω_a²)/2 ⇒ (f₂² − f₁²)² = (f_b² − f_a²)² + (C d⁻⁵)². Consequences: f₁ rises as d shrinks (observed: 0.802 → 0.823 Hz; the model with f_a = 0.793, f_b = 0.849 Hz gives 0.798 → 0.819), and Y saturates at |f_b² − f_a²| instead of → 0 → flattening at large d. Remedy: measure f_a, f_b with the other magnet removed; fit the detuned formula.
10.5 **Damping** (assumption e). Fit the envelope, get τ; frequency shift ≈ ½(1/(4πfτ))² — negligible for τ ~ 20–30 s — but the Lorentzian linewidth ≈ 1/(πτ) ≈ 0.01–0.02 Hz sets a floor on resolvability that no record length removes (feeds 10.1). Lumped-mass assumption d: < 1 % (Appendix A) — dismiss in one sentence, with the number.
10.6 **Random errors**: d (±0.5 mm; Δ(ln d) = 2 % at 21 mm), tracking noise (sub-pixel), frame timing jitter, release repeatability → half-range across trials; bin-centring (scalloping) ±1/(2T) unless interpolated.

### 11 Further methodology (≈200 w) — numbered modified procedure, INTERNATIONAL §11.2 style
1. Measure f_a, f_b of each spring alone (other magnet removed).
2. Track **both** magnets; form η = x₁ + x₂ and ξ = x₁ − x₂ before the FFT. Each normal coordinate contains essentially **one** frequency, so the two-peak resolution problem disappears — f₁ from η, f₂ from ξ, even at 45–60 mm. (Physically motivated, free, and the single best improvement.)
3. Record ≥ 90 s at 60 fps (or 30 fps), tripod, scale bar in frame; release A ≤ 3 mm; also one run at A = 10 mm per d for the amplitude study.
4. Re-analyse every record with rectangular / Hann / Blackman windows + 16× zero-padding + parabolic peak interpolation (INTERNATIONAL Table 10 style comparison), **and** with a time-domain two-damped-sinusoid least-squares fit (beats the 1/T limit; gives statistical Δf).
5. Fit models to (d, f₁, f₂) by weighted non-linear least squares: (a) dipole n = 5, (b) free exponent n, (c) finite-size disc model, (d) detuned model. Compare residuals/χ².
`analyse_tracks.py` does steps 4–5 for a folder of Tracker exports (see §6 of this plan).

### 12 Further data & discussion (≈350 w)
Table: gradient ± Δ for each method (rect/Hann/Blackman/2-tone fit/normal-coordinate FFT) — like INTERNATIONAL Fig. 26–27. Graph: f₂ vs amplitude (segments) with extrapolation to zero amplitude. Graph: ln Y vs ln d for old vs new data with both theory curves (point dipole and finite-size). Discuss which correction moved the gradient toward −5 and by how much; be explicit if a residual discrepancy remains and what it implies (e.g. near-field of real magnets ⇒ "dipole law valid only for d ≳ 4D").

### 13 Numerical simulation (≈300 w) — INTERNATIONAL §14 analogue
RK4 integration of the **full** equations m ẍ₁ = −k x₁ − [F(d + x₂ − x₁) − F(d)], m ẍ₂ = −k x₂ + [F(d + x₂ − x₁) − F(d)] with F = C r⁻⁴ (and optionally the finite-size F), using your measured k, m, C. Two uses: (i) predict f₂(A, d) → the amplitude correction used in §12; (ii) generate synthetic x(t) at your fps and T, push them through the *same* FFT pipeline → isolates the bias caused purely by the analysis (in my test at T = 40 s, window peak-picking gave gradients −3.3 to −3.5 while the two-tone fit recovered the true value). Code → Appendix G (`ee_simulations.py`).

### 14 Conclusion & extensions (≈220 w)
Answer the RQ in three sentences (f₁ constant within…; f₂² − f₁² falls as d⁻ⁿ with n = … ± …; dipole law holds for d ≥ …). Name the dominant systematic error and what fixed it. Extensions: vary μₘ (magnet grade) at fixed d to test the μₘ² dependence; attractive orientation (softening → instability at a critical d, a nice bifurcation); three springs; damping by eddy currents with a copper plate.

---

## 4. FFT error analysis — what to actually write (adapted from INTERNATIONAL to *frequency* measurement)

INTERNATIONAL needed **amplitudes** (SNR), so its corrections were energy-correction factors and baseline-noise subtraction. You need **frequencies**, so the relevant FFT physics is different — say this explicitly, it shows understanding:

1. **Resolution.** DFT X_k = Σ xₙ e^{−2πi kn/N}; bin spacing Δf_bin = fs/N = 1/T. Two tones are separable only if |f₂ − f₁| ≳ 1/T (rectangular) or ≳ 2/T (Hann, mainlobe twice as wide). Table from your data:

| d/mm | f₂−f₁ /Hz | bins at T = 20 s | T_min rect /s | T_min Hann /s |
|---|---|---|---|---|
| 21 | 0.240 | 4.8 | 4 | 8 |
| 25 | 0.141 | 2.8 | 7 | 14 |
| 28 | 0.097 | 1.9 | 10 | 21 |
| 31 | 0.074 | 1.5 | 14 | 27 |
| 35 | 0.039 | 0.8 | 26 | 51 |
| 40 | 0.009 | 0.2 | 111 | 222 |

2. **Leakage & windows.** Rectangular: narrowest mainlobe (best for two close peaks of similar height) but −13 dB sidelobes; Hann: 4-bin mainlobe, −31 dB; Blackman: 6-bin, −58 dB. Because both your peaks have similar amplitude and there is no strong interferer, the rectangular (or Hann) window is the right choice — the opposite of the usual "Hann by default" advice, and unlike INTERNATIONAL's case. Amplitude correction factors are **not** needed because you only read peak positions.
3. **Peak pulling.** When mainlobes overlap, the two maxima are pushed apart or merge; in my synthetic test (T = 20 s) the errors were ±2–7 mHz at 21–31 mm, ±10 mHz at 35 mm and total failure at 40 mm — quote your own numbers from `ee_simulations.py` Part A with your real T and fs.
4. **Scalloping / bin-centring.** True frequency anywhere within ±½ bin ⇒ ±1/(2T) unless you zero-pad and interpolate (parabolic on the three bins around the peak) — zero-padding interpolates the spectrum, it does *not* add resolution; make that distinction.
5. **Damping linewidth.** e^{−t/τ}cos(2πft) has FWHM ≈ 1/(πτ); with τ ≈ 25 s that is 0.013 Hz — comparable to f₂−f₁ at 40 mm, so beyond ~35 mm the modes cannot be separated by any spectral method with a single coordinate → justification for normal-coordinate tracking (§11 step 2) and time-domain fitting.
6. **Sampling.** fs = fps (30/60 Hz) ≫ 2f, no aliasing; but video frame timing jitters and frames drop → resample to a uniform grid before the FFT (the script reports the jitter).
7. **Propagation into your graph.** Δ(ln Y) = 2(f₁Δf₁ + f₂Δf₂)/Y grows from 0.04 to 1.1 across your range for Δf = 5 mHz ⇒ error bars must be drawn in ln-space and the fit weighted; otherwise the noisy large-d points dominate the gradient.

---

## 5. "Non-linearisation" improvements — three levels (use all three, in this order)

**Level 1 – analysis (cheap, do it now).** Replace "take logs then LINEST" with a weighted non-linear least-squares fit of Y = C d⁻ⁿ (and Y = C d⁻⁵ with n fixed) to the un-transformed data. Advantages to state: no point is discarded because ln(0) or ln(negative) fails; heteroscedastic errors are handled by the weights; C and n come with proper uncertainties. With your 6 points: free-exponent fit n = 3.5 ± 0.2 (vs ln–ln −3.75) — report both and explain the difference.

**Level 2 – physics (the systematic-error story).** The model linearised F(r) about d. Keep the full r⁻⁴ (or finite-size) force: perturbation formula (Appendix D) for the essay body, RK4 (§13) for the exact numbers, amplitude-segment analysis + extrapolation to A → 0 for the data. Report the size of the correction at each d in a table.

**Level 3 – model (finite size + detuning).** Fit (c) finite-size disc model and (d) detuned model from §11 step 5. Present the comparison of residual patterns as the decisive evidence for which assumption failed.

---

## 6. Lab / data checklist (you have lab access + raw tracks)

- [ ] Magnet dimensions D, t, mass M, grade (B_r) — needed for §10.3 and the theoretical C.
- [ ] Spring L, b, h, μ (weigh a known length), k by static loading; f_a, f_b of each spring alone.
- [ ] Balance measurement of F(r) for 4–6 values of r (checks r⁻⁴ and gives μₘ).
- [ ] Re-run d = 21…60 mm with both magnets tracked, T ≥ 90 s, A = 3 mm (and one A = 10 mm run each). Name files `d21_trial1.csv` etc.
- [ ] `python analyse_tracks.py "tracks/*.csv" --out results` → per-file diagnostic figure (x(t) + 2-tone fit, three-window FFT, f₂ vs amplitude), `summary.csv`, and the ln–ln gradient per method. Then form η, ξ from the two tracked magnets (one extra column operation) and run the same script on them.
- [ ] Sanity check of magnitudes: γ = 2π² m Y from your data vs γ = 6μ₀μₘ²/(π d⁵) from μₘ. If they differ by more than a factor ~2, revisit the coaxial-dipole assumption or the definition of d before writing §8.4.
- [ ] Re-derive f₂−f₁ requirement with your real T: if your original records were ~20–40 s, say in §7 that 35 mm was marginal and 40 mm unresolvable — it turns a weakness into evidence of understanding.

---

## 7. Figures & tables (final list)

Figures: 1 set-up sketch (done) · 2 mode sketch · 3 set-up photo · 4 x(t) with beats (small d) · 5 x(t) with no visible envelope (large d) · 6 FFT with two peaks · 7 F vs x (k) · 8 F vs r on balance (μₘ) · 9 ln Y vs ln d (Graph 1) · 10 Y vs d⁻⁵ (Graph 2) · 11 f₁ vs d (Graph 3) · 12 synthetic-spectrum resolution figure · 13 f₂ vs amplitude · 14 gradient by method (further data) · 15 old vs new data with dipole and finite-size theory · 16 RK4 vs linearised prediction.
Tables: 1 constants · 2 k data · 3/4 f₁, f₂ raw · 5 processed · 6 theory vs experiment (gradient, C) · 7 systematic-error summary (effect, sign, size, remedy) · 8 further raw data (appendix) · 9 gradient/C per method and per model · 10 simulated vs measured f₂ shift.

---

## 8. Appendices

A Cantilever stiffness 3EI/L³, effective mass 33/140 μL, and FEM check of lumped vs beam model (< 1 %).
B Dipole field → force between coaxial dipoles F = 3μ₀μₘ²/(2πr⁴) → γ = 6μ₀μₘ²/(πd⁵).
C Detuned normal modes: ω₁,₂² = ω̄² + κ ∓ √(Δ² + κ²) and its d⁻⁵ limit.
D Series expansion of F(d−ξ), α and β, Landau–Lifshitz amplitude shift (a/d)²(15g/8 − 125g²/48), comparison with RK4.
E Finite-size (charged-disc) magnet force and apparent gradient table.
F Uncertainty propagation (ΔY, Δ ln Y, Δ d⁻⁵, LINEST worst lines).
G Python: `analyse_tracks.py` (FFT windows, two-tone fit, segments), `ee_simulations.py` (RK4, synthetic FFT), `finite_size_magnets.py`.
H Further raw data tables (all windows/methods), as INTERNATIONAL Appendix A.

---

## 9. Ready-made numbers from my checks (replace with your own once you have real inputs)

- ln–ln fit, 6 points: gradient −3.75 ± 0.19, intercept 10.68 (matches your chart); weighted (Δf = 5 mHz): −3.48 ± 0.24. Free power law on raw Y: n = 3.5 ± 0.2.
- Non-linear force, A = 10 mm release: f₂ +9.8 % (21 mm), +1.1 % (35 mm); gradient −5.00 → −5.16. A = 2 mm: +0.3 % / +0.03 %.
- Finite-size discs, apparent gradient over 21–35 mm: D = 6 mm −4.87 · 10 mm −4.6 · 12 mm −4.5 · 15 mm −4.2 · 20 mm −3.75 · 25 mm −3.4 (t = 3–6 mm changes these by < 0.1).
- Lumped vs beam: Y differs by 0.03 % (M/μL = 3.3), 0.4 % (0.8), 1 % (0.4).
- Detuned model fitted with n = 5: f_a = 0.793 Hz, f_b = 0.849 Hz reproduces the f₁ drift 0.798 → 0.819 Hz.
- Synthetic tracks at T = 40 s: window peak-picking gradients −3.3 to −3.5, two-tone fit −4.9 (true −4.9) → analysis method alone can shift the gradient by > 1.

---

## 10. Pitfalls to avoid (examiner's eye)
- Do not present the 40 mm point as data; present the unresolved range as a measured limit of the method.
- Do not claim "f₁ constant" when it drifts 2.6 %; explain it (10.4).
- Do not copy IYPT slide figures or equations without derivation; cite the IYPT problem and derive the lumped model yourself.
- Keep all corrections quantitative (sign and size) — "small-angle-like" hand-waving loses Criterion D marks.
- Every graph: error bars, best/worst lines or LINEST uncertainty, units, gradient with uncertainty, theory line overlaid where possible.
- Word count: equations, tables, captions, appendices excluded — move derivations there ruthlessly.
