# KPZ314 Camera‑trap Practical — Student Guide

Welcome to Week‑8. You’ll use an expert‑validated camera‑trap dataset to learn how ecologists turn images into **inference** about communities, habitats, and species.

## What you will learn (outcomes)

- Why we use **independent events** and how the event gap affects inference.  
- How to compute and **interpret diversity** (α, β, γ) with **uncertainty**.  
- How to read **turnover** plots and tests (Bray–Curtis vs Jaccard, NMDS, PERMANOVA, dispersion).  
- How to fit **effort‑offset GLMs** and read **IRRs**.  
- Why **detectability** matters (p < 1): detection histories, naïve ψ̂(t), simple occupancy.  
- How to inspect **inter‑specific structure** (co‑occurrence residuals; activity overlap Δ̂) without claiming causation.

## Quick start

```r
renv::restore()                                # install packages for this project
quarto::quarto_render("practicum.qmd")         # open the HTML and answer prompts
```

Alternatively, you can try the below, but the restore() approach is recommended:

```r
source("dependencies.R")
pak::pkg_install(project_dependencies())
```

> The pipeline cache is pre‑built. Only run `targets::tar_make()` if the workbook tells you to recompute specific targets.

## Repo layout

```
/_targets.R         # pipeline definition (pre‑built cache in _targets/)
R/                  # helper functions (pure, unit‑tested)
practicum.qmd       # this workbook
data/               # observation table + site metadata + taxa
```

## Key concepts (glossary in one minute)

- **Event**: a detection separated by ≥ gap minutes from the previous detection of that species at that site.  
- **RAI**: events per 100 trap‑nights (controls for exposure).  
- **Hill numbers**: q=0 richness; q=1 exp(Shannon); q=2 inverse Simpson (dominance‑sensitive).  
- **γ diversity**: total diversity of a region (here: habitat), with **coverage** = S_obs/Chao1 indicating completeness.  
- **β turnover**: between‑site dissimilarity. **Jaccard** = composition; **Bray–Curtis** = composition + dominance.  
- **PERMANOVA**: tests **location** (centroid) differences; pair with a **dispersion** test to avoid confounding.  
- **IRR**: multiplicative change in event rate for a +1 unit (here **+1 SD**) change in a covariate; `typewet` compares wet to dry.  
- **Detectability (p)**: probability of detecting a species when present on a survey day. ψ̂(t) rises with effort if p < 1.  
- **Δ̂**: diel activity overlap (Ridout & Linkie), with bootstrap CIs.

## Troubleshooting (fast)

| Symptom | Fix |
| --- | --- |
| `tar_read()` can’t find an object | `targets::tar_make()` once; then re‑render |
| YAML error line 1–5 | Ensure the file starts with `---` and a valid header |
| Package not found | `renv::restore()` then **restart R** |
| Plotly not showing labels | Ignore; the static ggplot holds the same information |

## Frequently used targets

- **Event & effort**: `fig_gap`, `tab_effort`, `rai`  
- **Diversity**: `alpha`, `alpha_ci`, `gamma_all_ci`, `gamma_by_habitat`, `rare_raw`, `rare_region`  
- **Turnover**: `beta_permanova`, `fig_heat`, `fig_heat_jac`, `fig_nmds`, `fig_nmds_jac`, `nmds_stats`  
- **GLMs**: `glm_diag_all`, `glm_irrs_all`, `fig_glm_coef`, `fig_glm_pd_type`  
- **Detectability**: `fig_det_heat`, `psi_curves_all`, `fig_psi_curves`, `tab_occu`  
- **Inter‑specific**: `fig_coocc_heat`, `fig_coocc_heat_by_type`, `overlap_summary`, `fig_overlap_pairs`

## Rules-of-thumb for interpretation

- Quote **effect size + uncertainty** (CIs or p) every time.  
- Never claim causation from **co‑occurrence** or **overlap** alone—propose a falsification.  
- If **dispersion p < 0.05**, say so before interpreting PERMANOVA R².  
- If **p** is low, acknowledge that extra effort may change naïve ψ̂(t).

## Rebuild (optional)

```r
targets::tar_destroy(destroy = "objects")   # wipe old objects
targets::tar_make()                         # rebuild everything
```

Maintainers: Barry W. Brook (barry.brook@utas.edu.au)  
License: CC‑BY‑SA 4.0
