# ZISAR workshop slides — GSSI, L'Aquila, 30 October 2026

- `zisar_workshop.tex` — Beamer source (16:9, 10pt). Compile with `tectonic zisar_workshop.tex`.
- `figures.R` — regenerates every figure in `figures/` (`Rscript figures.R`; needs ggplot2,
  deldir, patchwork). Monte Carlo numbers are copied from `../tex/montecarlo.tex`; the
  illustrative figures and the synthetic 300-municipality region are simulated from the model.
  Section 9 fits a didactic Laplace estimator in base R on that region (numbers of the
  "Estimation on the map" slide, `figures/est_results.txt`); the fit is cached in
  `figures/est_cache.rds` — delete it to re-run (~3 min).
- `guida_relatore.tex` → `guida_relatore.pdf` — speaker guide in Italian: what to say slide by slide, the numbers, likely challenges. Under each slide a box with its one-sentence message and the bridge sentence to the next; the "Schede per le prove" part (script, titles-only page for the slides-off test, rehearsal timing sheet) is generated from those boxes. Part II (appendices) holds the full explanations, cross-referenced from the text.
- `logos/` — drop `unicatt.png` and `gssi.png` here (PNG or PDF); the title page picks
  them up automatically and compiles without them.
- `zisar_workshop.pdf` — compiled deck: 33 main slides (one with 4 overlay steps) + 12 backup. Fonts: Helvetica Neue via fontspec (XeTeX/tectonic), Palatino math.

## Timing (50 min)

| Block | Minutes |
|---|---|
| Camagni + plan | 4 |
| 1. The data problem (4-step map build-up + signature) | 5 |
| 2. Toolkits A/B/C + the gap | 9 |
| 3. The ZISAR model (picture, formal, SAR-on-eta vs y, multiplier, nested cases) | 12 |
| 4. Estimation (integral, Laplace/TMB, asymptotics) | 8 |
| 5. Monte Carlo (design + 4 result slides) | 7 |
| 6. Effects + application / York | 4 |
| 7. Open problems + take-aways | 3, then discussion |

If running late, drop "Precision of rho" and "Overdispersion attribution, fit and speed":
the message survives on the two slides before them.

## Backup slides (for Q&A)

estimation algorithm · lattice with/without rho · zero-inflation vs hurdle · three channels of
overdispersion and identification · Laplace accuracy · SAR vs CAR · choice of W and diagnostics ·
computational cost · general model with spatial switch · space–time extension · MC tables · references
