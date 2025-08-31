# Getting Started: KPZ314 Camera-trap Practical

*5-minute checklist for Week-8 Prac.*

1. **Open the project in RStudio**
   `File → Open Project… → KPZ314-camera-trap.Rproj`

2. **Restore the renv**

   ```r
   renv::restore()   # installs all R packages via renv
   ```

3. **Render the workbook**
   Click the **Render** button in the `practicum.qmd` tab.
   The HTML will open automatically; answer the TASK prompts inline.
   Alternatively, in the R console, type: `quarto::quarto_render("practicum.qmd")`

**Troubleshooting**

| Symptom                     | Fix                                              |
| --------------------------- | ------------------------------------------------ |
| `tar_read()` error          | In Console: `targets::tar_make()`                |
| YAMLException line 5        | Ensure the header starts with three dashes `---` |
| “package XYZ not available” | `renv::restore()` then restart R                 |

### Independent events: where to change the gap

- The default event gap is 5 minutes. To change it and rebuild affected outputs:

```r
# In _targets.R, edit this line
events_default <- build_events(raw, gap_min = 5L)

# Then rebuild the affected targets
targets::tar_make(names = c("events_default","rai","alpha","comm","gap_sensitivity","fig_gap"))
```

RAIs (rate of animal images) are per 100 trap‑nights and are available via:

```r
targets::tar_read(rai)
```
