# tutor/ — "The Theory, Step by Step" (teaching companion)

A 24-page first-year-level derivation booklet for the coupled magnetic–mechanical
oscillator, written as the study companion to `../theory_v2.tex` and
`../theory_companion.tex`.

* **`tutor.pdf`** — the compiled booklet (15 chapters, 13 figures, worked pilot
  numbers, self-test with answers, one-page cheat sheet).
* **`tutor.typ`** — its source. The sandbox has no LaTeX toolchain (no apt/GitHub
  access), so the document is authored in **Typst**, whose syntax mirrors the
  LaTeX source closely. If you want a LaTeX twin, `tutor.typ` translates
  section-for-section (boxes ↔ tcolorboxes, `$...$` ↔ `$...$`).
* **`compile.py`** — build: `python3 compile.py` (needs `pip install typst`;
  the wheel bundles the compiler, works offline).
* **`make_tut_figs.py`** — regenerates the teaching figures in `figs/`
  (matplotlib + scipy, house style: copper `orange!80!black`, clamp grey,
  N blue / S red). Four figures (`fig_prediction`, `fig_release`,
  `fig_damping`, `fig_nonlinear`) are copied from `../figs` so everything
  stays consistent with the essay's own numbers.

Pilot constants used throughout: f1 = 6.49 Hz, C = 1.33e-7 m^5/s^2,
tau_d = 21 s, d = 24.5–42 mm, A = 5 mm.
