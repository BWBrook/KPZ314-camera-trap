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
