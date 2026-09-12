# Section 5 — Theoretical model (standalone rewrite)

| file | what it is |
|---|---|
| `theory_section.tex` / `.pdf` | the complete Section 5 (14 pp, 15 figures, 38 equations, symbol table). Section counter is preset to 5. |
| `theory_figures.tex` / `.pdf` | **figures only** — the same 15 `\begin{figure}…\end{figure}` blocks (TikZ / pgfplots source), one after another, with the colour + magnet macros in the preamble. Copy-paste from here. |
| `make_fig_release.py` | Python (numpy + matplotlib) that generates Figure 13, `fig_release.pdf` / `.png` (ideal release vs hold-and-release, time series + FFT). Run `python3 make_fig_release.py`. |
| `fig_release.pdf` | output of the script; needed to compile both tex files. |

Compile: `tectonic -X compile theory_section.tex` (or `pdflatex` ×2). Same for `theory_figures.tex`.

## House style (identical to the Section 1 setup figure)
```latex
\colorlet{copper}{orange!80!black}   % leaf springs
\colorlet{clampgray}{gray!60}        % clamps
\colorlet{basegray}{gray!30}         % non-magnetic base
\colorlet{magnetS}{red!70}           % S pole
\colorlet{magnetN}{blue!70}          % N pole
\magnetSN{x}{y}{width}{height}       % magnet block, S left / N right, white pole letters
\magnetNS{x}{y}{width}{height}       % N left / S right
```
Magnets are always glued to the face of the strip at its tip (never sitting on top of it); the two facing poles are N–N.

## Merging into the main essay
Copy everything between `\section{Theoretical model}` and `\end{document}` into `sec05_theory.tex`; add to `main.tex`
`\usepackage{pgfplots}\pgfplotsset{compat=1.18}`, `\usetikzlibrary{patterns}`, the colour/magnet definitions above,
and the `notebox` tcolorbox; delete the `keybox` definition here (main.tex already has one); copy `fig_release.pdf` into `figs/`.
