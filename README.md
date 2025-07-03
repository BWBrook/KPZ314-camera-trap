# KPZ314 Camera-trap Practical

Hands-on lab for *Fauna of Tasmania* (KPZ314, Week 8).  
Students use an expert-validated MEWC camera-trap table to derive site- and community-level metrics in R.

## Repo structure

```
/bootstrap_env.R      # renv bootstrap
/_targets.R           # hidden pipeline (pre-built in _targets/)
R/                    # helper functions, one task each
practicum.qmd         # the workbook you knit
data/                 # master table + site + taxa
_targets/             # small (~1 MB) cache of pre-made targets

````

## Quick-start (students)

```bash
git clone https://github.com/UTasEco/KPZ314-camera-trap.git
cd KPZ314-camera-trap
Rscript bootstrap_env.R      # installs packages, links renv
quarto render practicum.qmd  # or click “Render” in RStudio
````

Or in RStudio, create a Project in the git directory then:
```r
source(bootstrap_env.R) # installs all R packages via renv
quarto::quarto_render("practicum.qmd") # render the md and HTML
```

> Pipeline cache is pre-built; students do not run tar_make() unless they intentionally want to recompute.

## Rebuilding the pipeline (optional)

```r
library(targets)
tar_destroy(destroy = "objects")   # wipe old objects
tar_make()                         # 7–10 s on a laptop
```

## FAQ

* **“tar_read can’t find objects”** – run `tar_make()` once or pull the latest commit with the `_targets/` folder.
* **“YAML parse error”** – ensure the first five lines of `practicum.qmd` match exactly the README sample.
* **“package not found”** – run `renv::restore()`.

Maintainers: Barry Brook (`barry.brook@utas.edu.au`).
License: CC-BY-SA 4.0.
