# KPZ314 Camera-trap Practical

Hands-on lab for *Fauna of Tasmania* (KPZ314, Week 8).  
Students use an expert-validated MEWC camera-trap table to derive site- and community-level metrics in R.

## Repo structure

```
/_targets.R           # hidden pipeline (pre-built in _targets/)
R/                    # helper functions, one task each
practicum.qmd         # the workbook you knit
data/                 # master table + site + taxa
_targets/             # small (~1 MB) cache of pre-made targets

````

## Quick-start (students)

In RStudio, unzip the directory then:
```r
renv::restore() # installs all R packages via renv
quarto::quarto_render("practicum.qmd") # render the md and HTML
```

Optional (faster installs):
```r
source("dependencies.R")
pak::pkg_install(project_dependencies())
```

> Pipeline cache is pre-built; students do not run tar_make() unless they intentionally want to recompute.

### Independent events, effort, and RAI

The pipeline parameterises “independent events” via a time gap and derives sampling effort and RAIs (events per 100 trap‑nights):

- Change the default gap in `_targets.R` by editing `events_default = build_events(raw, gap_min = 5L)`.
- Effort is summarised from site metadata (`op_days` → `trap_nights`), and RAIs are computed as `100 * events / trap_nights`.
- A gap‑sensitivity diagnostic (`fig_gap`) explores how metrics vary for gaps {1, 5, 10, 30} minutes.

### Diversity with uncertainty

We report alpha diversity using Hill numbers (q = 0, 1, 2) with bootstrap CIs:
- q0 = richness; q1 = exp(Shannon); q2 = inverse Simpson. See `alpha` (point estimates) and `alpha_ci` (CIs).
- Habitat-level gamma with site-bootstrap CIs is available via `gamma_by_habitat`; overall gamma via `gamma_all_ci`.
- Rarefaction plots annotate sampling coverage (S_obs/Chao1) at the observed endpoint, and include a wet vs dry overlay.

### Community turnover (paired distances)

We quantify between-site turnover with both Bray–Curtis (abundance-sensitive) and Jaccard (presence–absence). Visuals: `fig_heat_bc`, `fig_heat_jac`, `fig_nmds` (Bray), and `fig_nmds_jac` (Jaccard), with NMDS fit stats in `nmds_stats`. Inference combines PERMANOVA (location; `R²`, `p_perm`) and a dispersion test (`disp_p`) in `beta_permanova`. If `disp_p < 0.05`, interpret the location effect cautiously.

### Habitat use models (effort‑offset GLMs)

Focal‑species GLMs model event counts with a log(trap_nights) offset. We select Poisson unless overdispersed (φ > 1.5), then prefer NB if it clearly improves AIC (≥ 2), otherwise use quasi‑Poisson for robust SEs. Numeric covariates are z‑scored so IRRs are “per +1 SD”. Partial‑dependence plots show expected rates per 100 trap‑nights, holding other covariates at z=0. Quasi CIs are Wald with dispersion‑inflated SEs. See `glm_diag_all`, `glm_irrs_all`, `fig_glm_coef`, and `fig_glm_pd_type`.

### Detectability (daily histories, naïve ψ, occupancy demo)

We surface imperfect detectability (p < 1) in three steps:

- Daily detection histories per species (capped at 21 days) show when sites first record a species, faceted by habitat.
- Naïve occupancy curves ψ̂(t) = proportion of sites detected by day t, by habitat, illustrate how ψ̂(t) rises with effort and can differ between wet vs dry.
- Auto-picks a focal species by windowed detection rate (closest to 50% of sites detected within the window) to keep the demo pedagogical.
- Demonstrator occupancy model (unmarked::occu, ψ ~ type, p ~ 1) reports ψ for wet/dry and the per‑day detection probability p.

Assumptions: days without detections are treated as non‑detections; cameras are assumed operating continuously for their `op_days`. The occupancy demo is pedagogical (no observation covariates). We cap detection histories at 21 days by default.

### Inter‑specific structure (diagnostics)

- Co‑occurrence: fixed‑margins independence residuals (hypergeometric), visualised as a z‑score heatmap. High |z| highlights hypotheses, not proof of interaction; check within habitats.
- Diel activity overlap: Δ̂ (Ridout & Linkie) with bootstrap 95% CIs, contrasted by habitat. Estimator rule: Δ4 for n ≥ 75 per species, else Δ1.
- Practicum prompts mechanism + falsification (e.g., within‑habitat tests, covariates, activity windows).
 - Note: Δ̂ bootstraps are over events, not sites; true site-level uncertainty would resample sites.

## Rebuilding the pipeline (optional)

```r
targets::tar_destroy(destroy = "objects")   # wipe old objects
targets::tar_make()                         # 7–10 s on a laptop
```

## FAQ

* **“tar_read can’t find objects”** – run `tar_make()` once or pull the latest commit with the `_targets/` folder.
* **“YAML parse error”** – ensure the first five lines of `practicum.qmd` match exactly the README sample.
* **“package not found”** – run `renv::restore()`.

Maintainers: Barry Brook (`barry.brook@utas.edu.au`).
License: CC-BY-SA 4.0.
