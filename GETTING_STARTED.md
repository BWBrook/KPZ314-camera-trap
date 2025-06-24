# Getting Started: KPZ314 Camera-trap Practical

*5-minute checklist for Week-8 lab.*

1. **Install Git (once)**  
   * Windows: <https://git-scm.com/download/win>  
   * macOS: `brew install git`

2. **Clone the repo**

   ```bash
   git clone https://github.com/UTasEco/KPZ314-camera-trap.git
   ````

3. **Open the project in RStudio**
   `File → Open Project… → KPZ314-camera-trap.Rproj`

4. **Run the bootstrap**

   ```r
   source("bootstrap_env.R")   # installs all R packages via renv
   ```

5. **Render the workbook**
   Click the **Render** button in the `practicum.qmd` tab.
   The HTML will open automatically; answer the TASK prompts inline.

**Troubleshooting**

| Symptom                     | Fix                                              |
| --------------------------- | ------------------------------------------------ |
| `tar_read()` error          | In Console: `targets::tar_make()`                |
| YAMLException line 5        | Ensure the header starts with three dashes `---` |
| “package XYZ not available” | `renv::restore()` then restart R                 |
