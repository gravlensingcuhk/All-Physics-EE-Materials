# Fact-check of `ee_plan_v2.pdf`

Verdict: the skeleton is sound (it is the previous plan plus eight edits). Three of the edits are
real improvements, but v2 introduces **five factual errors** — two of them are the *sign / direction* of a
systematic effect, which is exactly what costs marks in the evaluation criterion — and its
assessment-criteria map is built on the **old (pre-2027) EE criteria**.

## Must fix before using

| # | Where in v2 | v2 says | Reality | Fix |
|---|---|---|---|---|
| 1 | §8 Assessment-criteria map | C "Critical thinking (12 marks)", D "Presentation", E "Engagement" | Those are the legacy criteria (34 marks, sessions before May 2027). From May 2027: A Framework for the essay 6, B Knowledge & understanding 6, C Analysis & line of argument 6, D Discussion & evaluation 8, E Reflection 4 (total 30). Reflection is assessed from the 500-word RPF, not the essay body. | Replace the table (corrected version below). Confirm your session with your supervisor. |
| 2 | §10.2 | "the oscillation samples a **softer**-than-dipole force" | The r⁻⁴ force is convex, so the anti-phase swing samples a **stiffer** average force (hardening): effective stiffness/γ = 1.04 at a/d = 0.1, 1.42 at a/d = 0.3. That is why f₂ rises with amplitude (+9.8 % at 21 mm, A = 10 mm) and the ln–ln gradient gets *steeper* (−5.16), not flatter. "Softer than dipole" is the finite-size effect of §10.3 — v2 merged the two. | Rewrite: amplitude non-linearity is a *hardening* effect that biases f₂ upward at small d (steepens the slope); finite size is the softening effect (flattens). |
| 3 | §10.4 | "as the coupling grows, **both** modes are pushed towards ω̄" | Only ω₁ approaches ω̄ (from below). ω₂² = ω̄² + κ + √(Δ²+κ²) increases monotonically: with f_a = 0.793, f_b = 0.849 Hz, f₁ goes 0.793 → 0.821 Hz while f₂ goes 0.849 → 1.64 Hz as κ grows. | "ω₁ rises toward ω̄ while ω₂ increases without bound; the detuning is only visible in f₁ and in the saturation of Y at large d." |
| 4 | §1.3 | IYPT 2023 problem 6 "**Magnetic Pendulums**" | IYPT 2023 problem 6 is "**Magnetic-Mechanical Oscillator**" (official problem set, IOC 24 July 2022). | Correct the title; cite the official problem set. |
| 5 | §10.5 | "Δd = 0.5 mm gives **2.4 %** on Y at 21 mm" | 0.5/21 = 2.4 % is the error on d; on Y ∝ d⁻⁵ it is 5 × 2.4 % = **12 %**. (v2's own §4.7 states the ×5 rule correctly.) | "2.4 % on d, i.e. ≈ 12 % on Y." |

## Half-right — keep the idea, fix the reasoning

| # | Where | Issue | Fix |
|---|---|---|---|
| 6 | Change 3 / §5.1 / §6.4.3 — dynamic calibration | Plucking gives **k/m** for each spring, which is exactly what f₁ and the detuning analysis need. But the absolute test C = 3μ₀μₘ²/(π³m) still needs **m** separately, and plucking cannot give it; you still need static k (m = k/ω²) or M + (33/140)μL. So the static k is *not* "just a cross-check". Also, "pluck each **bare** spring" is wrong: pluck the spring **with its own magnet mounted** and the other magnet removed (M is part of m). | Keep dynamic calibration as primary for f_a, f_b; keep static k as a required measurement for m. |
| 7 | Change 4 / §4.6 — detrending | Detrending is correct practice (the analysis script already removes a linear trend). The justification is wrong: the lowest peak is f₁ ≈ 0.80–0.82 Hz at **every** d (≈ 16 bins from DC at T = 20 s); it is f₂ that moves, toward f₁, not toward DC. The real reason detrending matters: in pixel coordinates the DC offset can be ~50× the oscillation amplitude, and the rectangular-window leakage of that offset at bin 16 (≈ 1/(16π) ≈ −34 dB) is then comparable to the signal peak. Also "the offset term in eq. (4)" — eq. (4) has no offset term (F(d) was subtracted). | Keep the step; replace the reasoning with the offset/leakage argument. |
| 8 | Change 5 / §10.2 — tilt | θ = 3A/(2L) is right (cantilever tip slope), θ = 0.10 rad = 5.7°. "cos θ ≈ 0.5 %" should read "1 − cos θ ≈ 0.5 %". Dipole–dipole energy with tilts ∝ 3cosθ₁cosθ₂ − cos(θ₁−θ₂): mirror-symmetric (anti-phase) tilt weakens the interaction by 0.5 %, same-direction (in-phase) tilt by 1.5 %. Direction of the claim (force slightly overestimated) is right; the size is ≲ 1.5 % and oscillatory — one line is enough. | Fix the wording and give both numbers. |
| 9 | Change 8 — reduced χ² | Good, but χ²/dof is only meaningful with genuine standard uncertainties. Half-range of 3–5 trials is not a σ, so quote χ²/dof as a *relative* model-comparison figure, not as an absolute goodness-of-fit. | Add one sentence saying so. |
| 10 | §4.5 (and the previous plan — my error too) | "FWHM ≈ 1/(πτ)" is the **power**-spectrum linewidth. Most FFT displays (Tracker, Capstone) show **amplitude**, whose FWHM is √3/(πτ) ≈ 0.022 Hz at τ = 25 s. | Use 0.022 Hz — it makes the resolution argument stronger. |

## Correct and worth keeping (verified)

- Word budget: the previous plan really summed to 4250 (my mistake); v2's 3950 is right. (v2 then says "never shrink §10–13" while shrinking them by 130 words — harmless, just inconsistent.)
- d = g + t (face gap + one magnet thickness = centre-to-centre for identical magnets). Tip: read g from a video frame with a scale bar rather than with calipers between repelling magnets on flexible springs.
- Horizontal error bars Δ(ln d) = Δd/d (4.8 % at 21 mm for Δd = 1 mm), comparable to Δ(ln Y) ≈ 0.04 there; propagates as 5Δd/d into d⁻⁵.
- Fig. 3 (synthetic spectra: two peaks at 21 mm, one merged peak at 40 mm) matches my simulation. Caption says 0.245 Hz ≈ 5 bins; the table says 0.240 Hz = 4.8 bins — trivial.
- Fig. 4 caption "possible downward curvature at small d = finite-size signature": residuals of the 6-point fit are (−, +, −, +, +, −), i.e. concave-down, which is the finite-size direction (detuning would give concave-up). With 6 points it is only suggestive; "possible" is the right hedge.
- Writing order, raw-video archiving, FFT settings as controlled variables: fine.

## Corrected assessment-criteria map (first assessment May 2027, 30 marks)

| Criterion (marks) | Where the plan earns it |
|---|---|
| A Framework for the essay (6) | Focused RQ (§1.3); variables incl. FFT settings (§4); method choices justified, attempt 1 → attempt 2 narrative (§6.3); structure and presentation conventions (numbered figures/tables, appendices) — presentation is now part of A. |
| B Knowledge and understanding (6) | Full derivation chain (§5); FFT/DFT theory (§6.1); correct vocabulary (normal modes, coupling, spectral leakage, resolution). |
| C Analysis and line of argument (6) | Quantitative hypothesis; processed data with propagated uncertainties (§8); theory–experiment comparison (§8.4); the argument runs prediction → measurement → two error families → resolved conclusion (Fig. 1 of v2). |
| D Discussion and evaluation (8) | Each systematic effect with **correct sign**, size and remedy (§10); improvements quantified (§11–12); simulation as independent check (§13); honest statement of the residual discrepancy and the model's validity range (§14). |
| E Reflection (4) | Not the essay body: the 500-word reflective statement (RPF) — record the pivots as they happen (period method → FFT; FFT limit → normal coordinates; linear → non-linear model). |

*Check with your supervisor which model applies to your session; the legacy 34-mark model applies only to sessions before May 2027.*
