# Extended Essay draft – coupled magnetic-mechanical oscillator

## Files

| File | What it is |
|---|---|
| `EE_essay.pdf` | The compiled essay (55 pages: 43 body pages + bibliography + appendices). |
| `EE_full.tex` | **Single-file LaTeX source** (everything inlined, only needs `figs/`). Use this in Overleaf. |
| `main.tex` + `sec*.tex`, `appendices.tex`, `bibliography.tex` | The same essay split into one file per section (easier to edit). `main.tex` reads `numbers.tex`, `numbers2.tex` and `tables/*.tex`. |
| `figs/*.pdf` | All figures (vector PDF). |
| `make_figures.py` | Regenerates every data figure, `numbers.tex` and `tables/*.tex` from the data at the top of the file (≈1 min). |
| `fit_models.py` | The five model fits of §12.3 and a few simulation numbers (`numbers2.tex`, `tables/tab_modelfits.tex`). |
| `../EE_plan/analyse_tracks.py` | Analysis tool for real Tracker exports (§11). |

## Compile

```
tectonic -X compile EE_full.tex        # or: pdflatex EE_full (twice)
```
Packages: geometry, amsmath, amssymb, bm, graphicx, booktabs, array, multirow, caption,
subcaption, enumitem, xcolor, tikz, tcolorbox, listings, hyperref, setspace, lmodern.

## Things you must replace before submission

* **Red text in square brackets** = values to be measured or filled in (magnet size and mass,
  spring dimensions, k, f_a, f_b, μ_m, photo, Tracker screenshot, per-trial frequencies, references' access dates).
* **Everything labelled SIMULATED** (§12.1, §12.2, Fig. 15–18, Table 12, Table 14) was produced by passing
  synthetic records through the analysis pipeline. Replace with the real re-measurement, or – if you do not
  re-measure – keep it but say explicitly that it is a simulation of the improved method (the text is already
  written that way).
* Assumed constants in `make_figures.py` (marked `ASSUMED`): Δf = 5 mHz, Δd = 0.5 mm, τ = 25 s,
  tracking noise 0.1 mm, L = 150 mm (for the tilt estimate). Change them and re-run to update all numbers.
* Word count: the body prose is ≈10 000 words (IB counts main text only, not captions/tables/equations/appendices).
  Trim towards 4000: the FFT background (§6.1), the evaluation sub-sections and §12 are the natural places to cut;
  move derivations to the appendices rather than deleting them.


## Update (RQ revision)
* Research question now names f_avg=(f1+f2)/2 and f_beat=f2-f1 (title page + Section 1.3); the falsifiable test
  f2^2-f1^2 = 2 f_beat f_avg ∝ d^-5 is stated in the sentence after the RQ.
* Beat theory moved from the Introduction to Section 5.4 (eq:solution, eq:beats, eq:beats_identity);
  convention: f_beat = f2 - f1 = 1/T_beat (observable beat frequency), so f2^2-f1^2 = 2 f_beat f_avg (NOT 4).
* New Graph 4 (figs/G4_favg_fbeat.pdf, Section 8.5) + new macros in numbers.tex (favg*, fbeat*, Tbeat*, gradBeat*).
* Figure 3 (Introduction) is now the MEASURED beat trace: regenerate figs/beats_measured_time.pdf with
      python3 plot_beats.py springA.txt springB.txt --tmin 50 --tmax 142 --panels time --out figs/beats_measured_time.pdf
  (the file currently in figs/ is a synthetic PLACEHOLDER).
