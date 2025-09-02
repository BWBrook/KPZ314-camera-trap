# Getting Started — KPZ314 Camera‑trap Practical (Revised)

**Goal:** get you from zero to a rendered workbook in **5 minutes**.

## 1) Open the project in RStudio
`File → Open Project… → KPZ314-camera-trap.Rproj`

## 2) Restore the project libraries
```r
renv::restore()        # installs the exact package set for this project
```

If installs are slow or failing, you can also try:
```r
source("dependencies.R")
pak::pkg_install(project_dependencies())
```

## 3) Render the workbook
Click **Render** in the `practicum.qmd` tab, or run:
```r
quarto::quarto_render("practicum.qmd")
```

The HTML opens automatically.  
Answer the **TASK** prompts inline or on separate notes (short, evidence‑first).

---

## What you’re looking at

- **Independent events**: we merge bursty triggers into events using a gap (default **5 min**).
- **Effort**: rates are **events per 100 trap‑nights (RAI)** so sites are comparable.
- **Diversity**: Hill numbers with **bootstrap CIs**; rarefaction with **coverage**.
- **Turnover**: paired distances (Bray vs Jaccard), NMDS, PERMANOVA + **dispersion**.
- **Habitat use**: GLMs with effort offset; coefficients reported as **IRRs**.
- **Detectability**: daily detection histories, ψ̂(t) curves, simple occupancy.
- **Inter‑specific**: co‑occurrence residuals; activity overlap Δ̂ with CIs.

---

## Changing the “independent event” gap

The default gap is **5 minutes**. To change it:

```r
# In _targets.R, edit:
events_default <- build_events(raw, gap_min = 5L)

# Then rebuild the affected targets only:
targets::tar_make(names = c("events_default","rai","alpha","comm","gap_sensitivity","fig_gap"))
```

**Tip:** Use `fig_gap` to justify your choice—pick the **smallest stable** gap.

---

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| `tar_read()` error | In Console: `targets::tar_make()` |
| YAMLException at top | Ensure the file begins with `---` and a valid header |
| “package XYZ not available” | `renv::restore()` then restart R |
| Quarto render hangs | Close all devices (`graphics.off()`), restart R, render again |

---

## After class

Want to explore further, or do something different for your assessment task?  
Create `explore.R` and call `targets::tar_read()` to pull any object (tables, matrices,  
ggplots) into your session. Everything is inspectable and reproducible.
