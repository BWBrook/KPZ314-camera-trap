# Clone the repo 
# Install Git for Windows or Git for Mac
system("git clone https://github.com/BWBrook/KPZ314-camera-trap.git")

# Run source("bootstrap_env.R") (or just open RStudio, and renv will activate)
renv::restore() # if prompted

# Knit the workbook
quarto::quarto_render("practicum.qmd")

