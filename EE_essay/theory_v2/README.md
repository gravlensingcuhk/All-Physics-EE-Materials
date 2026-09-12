# theory_v2 — new Section 5 (theory) + Section 6 (extension)

* `theory_v2.tex` / `theory_v2.pdf` — 8 pages, self-contained (same house-style preamble as `theory_section/`).
  Build: `tectonic -X compile theory_v2.tex` (needs the four figure PDFs below).
* `make_figs_v2.py` — makes `fig_prediction.pdf`, `fig_release.pdf`, `fig_damping.pdf`, `fig_nonlinear.pdf`
  (real time base: f1 = 6.49 Hz, C = 1.33e-7 m^5 s^-2, tau_d = 21 s — pilot-run values; re-run with the new data).
* `nonlinear_check.py` — numerical check of every non-linear formula in Section 6 against direct integration.

To drop into the essay: copy the body between `\section{Theoretical model}` and `\end{thebibliography}`
into `sec05_theory.tex` (the extension goes to the evaluation/extension section); the macros
`\muz \mum \dd \un \Lag`, the colours and `\magnetSN/\magnetNS` must be in `main.tex`.
